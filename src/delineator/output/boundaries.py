"""Watershed boundary generation from raster to vector."""

from pathlib import Path
from typing import Optional
import tempfile

import numpy as np
import rasterio
from rasterio.features import shapes
import geopandas as gpd
from shapely.geometry import shape
from whitebox import WhiteboxTools

from delineator.config import OutputFormat
from delineator.geometry.tin import StudyPoint
from delineator.exceptions import OutputError
from delineator.utils.logging import get_logger

# Conversion constants
SQFT_PER_ACRE = 43560
SQFT_PER_SQMILE = 27878400


class WatershedBoundaryGenerator:
    """Generates vector boundaries from watershed raster."""

    def __init__(self):
        """Initialize the boundary generator."""
        self.logger = get_logger(__name__)
        self.wbt = WhiteboxTools()
        self.wbt.set_verbose_mode(False)

    def generate(
        self,
        watershed_raster_path: Path,
        dem_path: Path,
        study_points: list[StudyPoint],
        output_path: Path,
        output_format: OutputFormat = OutputFormat.SHAPEFILE,
        flow_dir_path: Optional[Path] = None,
        flow_acc_path: Optional[Path] = None,
        pour_points_path: Optional[Path] = None,
        cumulative: bool = True,
    ) -> gpd.GeoDataFrame:
        """Generate watershed boundaries from raster.

        Args:
            watershed_raster_path: Path to watershed raster.
            dem_path: Path to DEM raster for elevation statistics.
            study_points: List of study points for labeling.
            output_path: Path for output vector file.
            output_format: Output format (shapefile, geojson, geopackage).
            flow_dir_path: Path to flow direction raster (for cumulative mode).
            flow_acc_path: Path to flow accumulation raster (for cumulative mode).
            pour_points_path: Path to snapped pour points (for cumulative mode).
            cumulative: Whether to compute cumulative (overlapping) watersheds.

        Returns:
            GeoDataFrame with watershed polygons and statistics.
        """
        self.logger.info("Generating watershed boundaries")

        # Read DEM for elevation statistics
        with rasterio.open(dem_path) as dem_src:
            dem_data = dem_src.read(1)
            dem_nodata = dem_src.nodata
            cell_size = dem_src.res[0]
            transform = dem_src.transform
            crs = dem_src.crs

        if cumulative and flow_dir_path and pour_points_path:
            # Cumulative mode: compute each watershed individually
            gdf = self._generate_cumulative(
                flow_dir_path=flow_dir_path,
                flow_acc_path=flow_acc_path,
                pour_points_path=pour_points_path,
                dem_data=dem_data,
                dem_nodata=dem_nodata,
                cell_size=cell_size,
                transform=transform,
                crs=crs,
                study_points=study_points,
            )
        else:
            # Exclusive mode: use the combined watershed raster
            gdf = self._generate_exclusive(
                watershed_raster_path=watershed_raster_path,
                dem_data=dem_data,
                dem_nodata=dem_nodata,
                cell_size=cell_size,
                study_points=study_points,
            )

        # Save to file
        self._save_output(gdf, output_path, output_format)

        return gdf

    def _generate_exclusive(
        self,
        watershed_raster_path: Path,
        dem_data: np.ndarray,
        dem_nodata: Optional[float],
        cell_size: float,
        study_points: list[StudyPoint],
    ) -> gpd.GeoDataFrame:
        """Generate boundaries from exclusive watershed raster."""
        with rasterio.open(watershed_raster_path) as src:
            watershed_data = src.read(1)
            transform = src.transform
            crs = src.crs
            nodata = src.nodata

        # Create mask for valid data
        if nodata is not None:
            mask = watershed_data != nodata
        else:
            mask = np.ones(watershed_data.shape, dtype=bool)

        # Vectorize watersheds
        self.logger.info("Vectorizing watershed raster")
        polygons = []
        watershed_ids = []

        for geom, value in shapes(watershed_data, mask=mask, transform=transform):
            if value == nodata or value == 0:
                continue

            polygon = shape(geom)
            if polygon.is_valid and not polygon.is_empty:
                polygons.append(polygon)
                watershed_ids.append(int(value))

        if not polygons:
            raise OutputError("No valid watershed polygons generated")

        self.logger.info(f"Generated {len(polygons)} raw watershed polygons")

        # Create GeoDataFrame
        gdf = gpd.GeoDataFrame(
            {"watershed_id": watershed_ids, "geometry": polygons},
            crs=crs,
        )

        # Dissolve polygons by watershed_id to merge disconnected parts
        gdf = gdf.dissolve(by="watershed_id", as_index=False)
        self.logger.info(f"Dissolved to {len(gdf)} unique watersheds")

        # Add study point names
        point_names = {sp.id: sp.name or f"Point_{sp.id}" for sp in study_points}
        gdf["point_name"] = gdf["watershed_id"].map(point_names)

        # Compute statistics (imperial units)
        self._compute_statistics(gdf, watershed_data, dem_data, dem_nodata, cell_size)

        return gdf

    def _generate_cumulative(
        self,
        flow_dir_path: Path,
        flow_acc_path: Optional[Path],
        pour_points_path: Path,
        dem_data: np.ndarray,
        dem_nodata: Optional[float],
        cell_size: float,
        transform,
        crs,
        study_points: list[StudyPoint],
    ) -> gpd.GeoDataFrame:
        """Generate cumulative watershed boundaries.

        Each watershed is computed independently, resulting in overlapping polygons.
        """
        self.logger.info("Computing cumulative watersheds")

        flow_dir_abs = str(Path(flow_dir_path).resolve())
        pour_points = gpd.read_file(pour_points_path)

        polygons = []
        watershed_ids = []
        areas_sqft = []
        areas_acres = []
        areas_sqmi = []
        min_elevs = []
        max_elevs = []
        mean_elevs = []
        point_names = []

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            for idx, row in pour_points.iterrows():
                pt_id = row['id']
                pt_name = row.get('name', f"Point_{pt_id}")

                self.logger.info(f"  Processing watershed for point {pt_id} ({pt_name})")

                # Create single-point shapefile
                single_pt_path = tmpdir / f"pt_{pt_id}.shp"
                single_gdf = gpd.GeoDataFrame(
                    {"id": [pt_id], "name": [pt_name], "geometry": [row.geometry]},
                    crs=pour_points.crs,
                )
                single_gdf.to_file(single_pt_path)

                # Delineate watershed for this point
                ws_output = tmpdir / f"ws_{pt_id}.tif"
                try:
                    result = self.wbt.watershed(
                        d8_pntr=flow_dir_abs,
                        pour_pts=str(single_pt_path.resolve()),
                        output=str(ws_output.resolve()),
                    )
                    if result != 0 or not ws_output.exists():
                        self.logger.warning(f"  Failed to delineate watershed for point {pt_id}")
                        continue
                except Exception as e:
                    self.logger.warning(f"  Error delineating watershed for point {pt_id}: {e}")
                    continue

                # Read and vectorize this watershed
                with rasterio.open(ws_output) as src:
                    ws_data = src.read(1)
                    ws_nodata = src.nodata
                    ws_transform = src.transform

                # Create mask
                if ws_nodata is not None:
                    ws_mask = ws_data != ws_nodata
                else:
                    ws_mask = ws_data > 0

                # Compute stats for this watershed
                cell_count = np.sum(ws_mask)
                area_sqft = cell_count * (cell_size ** 2)
                area_acres = area_sqft / SQFT_PER_ACRE
                area_sqmi = area_sqft / SQFT_PER_SQMILE

                # Elevation stats
                dem_values = dem_data[ws_mask]
                if dem_nodata is not None:
                    dem_values = dem_values[dem_values != dem_nodata]

                if len(dem_values) > 0:
                    min_elev = float(np.min(dem_values))
                    max_elev = float(np.max(dem_values))
                    mean_elev = float(np.mean(dem_values))
                else:
                    min_elev = max_elev = mean_elev = None

                # Vectorize
                for geom, value in shapes(ws_data.astype(np.int32), mask=ws_mask, transform=ws_transform):
                    polygon = shape(geom)
                    if polygon.is_valid and not polygon.is_empty:
                        polygons.append(polygon)
                        watershed_ids.append(pt_id)
                        point_names.append(pt_name)
                        areas_sqft.append(area_sqft)
                        areas_acres.append(area_acres)
                        areas_sqmi.append(area_sqmi)
                        min_elevs.append(min_elev)
                        max_elevs.append(max_elev)
                        mean_elevs.append(mean_elev)
                        break  # Only need one polygon per watershed

        if not polygons:
            raise OutputError("No valid watershed polygons generated")

        self.logger.info(f"Generated {len(polygons)} cumulative watershed polygons")

        # Create GeoDataFrame
        gdf = gpd.GeoDataFrame({
            "watershed_id": watershed_ids,
            "point_name": point_names,
            "area_sqft": areas_sqft,
            "area_acres": areas_acres,
            "area_sqmi": areas_sqmi,
            "elev_min": min_elevs,
            "elev_max": max_elevs,
            "elev_mean": mean_elevs,
            "geometry": polygons,
        }, crs=crs)

        return gdf

    def _compute_statistics(
        self,
        gdf: gpd.GeoDataFrame,
        watershed_data: np.ndarray,
        dem_data: np.ndarray,
        dem_nodata: Optional[float],
        cell_size: float,
    ) -> None:
        """Compute statistics for each watershed in imperial units."""
        self.logger.info("Computing watershed statistics")

        areas_sqft = []
        areas_acres = []
        areas_sqmi = []
        min_elevs = []
        max_elevs = []
        mean_elevs = []

        for idx, row in gdf.iterrows():
            watershed_id = row["watershed_id"]

            # Create mask for this watershed
            ws_mask = watershed_data == watershed_id

            # Get DEM values within watershed
            dem_values = dem_data[ws_mask]

            # Filter out nodata
            if dem_nodata is not None:
                dem_values = dem_values[dem_values != dem_nodata]

            # Compute area in imperial units (assuming cell_size is in feet)
            cell_count = np.sum(ws_mask)
            area_sqft = cell_count * (cell_size ** 2)
            area_acres = area_sqft / SQFT_PER_ACRE
            area_sqmi = area_sqft / SQFT_PER_SQMILE

            areas_sqft.append(area_sqft)
            areas_acres.append(area_acres)
            areas_sqmi.append(area_sqmi)

            # Compute elevation statistics (already in feet)
            if len(dem_values) > 0:
                min_elevs.append(float(np.min(dem_values)))
                max_elevs.append(float(np.max(dem_values)))
                mean_elevs.append(float(np.mean(dem_values)))
            else:
                min_elevs.append(None)
                max_elevs.append(None)
                mean_elevs.append(None)

        gdf["area_sqft"] = areas_sqft
        gdf["area_acres"] = areas_acres
        gdf["area_sqmi"] = areas_sqmi
        gdf["elev_min"] = min_elevs
        gdf["elev_max"] = max_elevs
        gdf["elev_mean"] = mean_elevs

    def _save_output(
        self,
        gdf: gpd.GeoDataFrame,
        output_path: Path,
        output_format: OutputFormat,
    ) -> None:
        """Save GeoDataFrame to specified format."""
        output_path.parent.mkdir(parents=True, exist_ok=True)

        if output_format == OutputFormat.SHAPEFILE:
            gdf.to_file(output_path, driver="ESRI Shapefile")
        elif output_format == OutputFormat.GEOJSON:
            gdf.to_file(output_path, driver="GeoJSON")
        elif output_format == OutputFormat.GEOPACKAGE:
            gdf.to_file(output_path, driver="GPKG")
        else:
            raise OutputError(f"Unsupported output format: {output_format}")

        self.logger.info(f"Boundaries saved to: {output_path}")
