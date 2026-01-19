"""Watershed boundary generation from raster to vector."""

from pathlib import Path
from typing import Optional
import tempfile

import numpy as np
import rasterio
from rasterio.features import shapes
import geopandas as gpd
from shapely.geometry import shape, Polygon, MultiPolygon, LineString
from whitebox import WhiteboxTools
import ezdxf
from rasterio.transform import rowcol

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
                dem_path=dem_path,
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
                dem_path=dem_path,
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
        dem_path: Path,
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

        # Compute longest flow path statistics
        with tempfile.TemporaryDirectory() as tmpdir:
            lfp_stats = self._compute_longest_flowpath_stats(
                dem_path=dem_path,
                basins_path=watershed_raster_path,
                work_dir=Path(tmpdir),
            )

        # Add flow path columns
        lfp_lengths = []
        lfp_slopes = []
        for _, row in gdf.iterrows():
            ws_id = int(row["watershed_id"])
            if ws_id in lfp_stats:
                length, slope = lfp_stats[ws_id]
                lfp_lengths.append(length)
                lfp_slopes.append(slope)
            else:
                lfp_lengths.append(None)
                lfp_slopes.append(None)

        gdf["lfp_length"] = lfp_lengths
        gdf["lfp_slope"] = lfp_slopes

        return gdf

    def _generate_cumulative(
        self,
        flow_dir_path: Path,
        flow_acc_path: Optional[Path],
        pour_points_path: Path,
        dem_path: Path,
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
        lfp_lengths = []
        lfp_slopes = []

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

                # Compute longest flow path stats for this watershed
                lfp_stats = self._compute_longest_flowpath_stats(
                    dem_path=dem_path,
                    basins_path=ws_output,
                    work_dir=tmpdir,
                )
                # Get stats for this watershed (ID will be 1 in single-watershed raster)
                if lfp_stats:
                    # The single watershed will have value from the raster
                    lfp_length, lfp_slope = next(iter(lfp_stats.values()), (None, None))
                else:
                    lfp_length, lfp_slope = None, None

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
                        lfp_lengths.append(lfp_length)
                        lfp_slopes.append(lfp_slope)
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
            "lfp_length": lfp_lengths,
            "lfp_slope": lfp_slopes,
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

    def _compute_longest_flowpath_stats(
        self,
        dem_path: Path,
        basins_path: Path,
        work_dir: Path,
    ) -> dict[int, tuple[float, float]]:
        """Compute longest flow path statistics for each watershed.

        Uses WhiteboxTools longest_flowpath to find the most hydrologically
        distant path to the pour point, then computes its length and average slope.

        Args:
            dem_path: Path to DEM raster.
            basins_path: Path to watershed/basins raster.
            work_dir: Working directory for temporary outputs.

        Returns:
            Dict mapping watershed_id to (path_length_ft, avg_slope_ft_per_ft).
        """
        self.logger.info("Computing longest flow path statistics")

        lfp_output = work_dir / "longest_flowpaths.shp"

        try:
            result = self.wbt.longest_flowpath(
                dem=str(dem_path.resolve()),
                basins=str(basins_path.resolve()),
                output=str(lfp_output.resolve()),
            )
            if result != 0 or not lfp_output.exists():
                self.logger.warning("Failed to compute longest flowpaths")
                return {}
        except Exception as e:
            self.logger.warning(f"Error computing longest flowpaths: {e}")
            return {}

        # Read the flowpath vectors
        try:
            lfp_gdf = gpd.read_file(lfp_output)
        except Exception as e:
            self.logger.warning(f"Error reading longest flowpath output: {e}")
            return {}

        if lfp_gdf.empty:
            return {}

        stats = {}

        with rasterio.open(dem_path) as dem_src:
            dem_data = dem_src.read(1)
            dem_nodata = dem_src.nodata
            transform = dem_src.transform

            for idx, row in lfp_gdf.iterrows():
                # Get watershed ID - WhiteboxTools uses 'BASIN' field
                ws_id = row.get("BASIN", row.get("basin", idx + 1))

                geom = row.geometry
                if geom is None or geom.is_empty:
                    continue

                # Handle both LineString and MultiLineString
                if hasattr(geom, "geoms"):
                    # MultiLineString - take the longest one
                    lines = list(geom.geoms)
                    geom = max(lines, key=lambda g: g.length)

                if not isinstance(geom, LineString):
                    continue

                # Calculate path length (in CRS units, assumed feet)
                path_length = geom.length

                if path_length <= 0:
                    continue

                # Sample elevations along the path at regular intervals
                # Use more sample points for longer paths
                num_samples = max(10, min(100, int(path_length / 50)))
                sample_distances = np.linspace(0, path_length, num_samples)

                elevations = []
                for dist in sample_distances:
                    point = geom.interpolate(dist)
                    x, y = point.x, point.y

                    # Convert to raster indices
                    try:
                        r, c = rowcol(transform, x, y)
                        if 0 <= r < dem_data.shape[0] and 0 <= c < dem_data.shape[1]:
                            elev = dem_data[r, c]
                            if dem_nodata is None or elev != dem_nodata:
                                elevations.append((dist, float(elev)))
                    except (IndexError, ValueError):
                        continue

                if len(elevations) < 2:
                    continue

                # Compute average slope along the path
                # Method: Average of segment slopes (weighted by segment length)
                total_weighted_slope = 0.0
                total_weight = 0.0

                for i in range(1, len(elevations)):
                    dist_prev, elev_prev = elevations[i - 1]
                    dist_curr, elev_curr = elevations[i]

                    segment_length = dist_curr - dist_prev
                    if segment_length > 0:
                        segment_slope = abs(elev_prev - elev_curr) / segment_length
                        total_weighted_slope += segment_slope * segment_length
                        total_weight += segment_length

                if total_weight > 0:
                    avg_slope = total_weighted_slope / total_weight
                else:
                    # Fallback: overall slope
                    elev_start = elevations[0][1]
                    elev_end = elevations[-1][1]
                    avg_slope = abs(elev_start - elev_end) / path_length

                stats[int(ws_id)] = (path_length, avg_slope)

        self.logger.info(f"Computed flow path stats for {len(stats)} watersheds")
        return stats

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
        elif output_format == OutputFormat.DXF:
            self._save_dxf(gdf, output_path)
        else:
            raise OutputError(f"Unsupported output format: {output_format}")

        self.logger.info(f"Boundaries saved to: {output_path}")

    def _save_dxf(self, gdf: gpd.GeoDataFrame, output_path: Path) -> None:
        """Save GeoDataFrame as DXF file with LWPOLYLINE entities."""
        doc = ezdxf.new("R2010")
        msp = doc.modelspace()

        # Register XDATA application ID
        doc.appids.add("DELINEATOR")

        # Create layers
        doc.layers.add("WATERSHEDS", color=7)  # White/default
        doc.layers.add("WATERSHEDS_HOLES", color=8)  # Gray

        # Color palette for watersheds (AutoCAD color indices)
        colors = [1, 2, 3, 4, 5, 6, 30, 40, 50, 60, 70, 80, 90, 100]

        for idx, row in gdf.iterrows():
            watershed_id = int(row["watershed_id"])
            point_name = row.get("point_name", f"Point_{watershed_id}")
            area_sqft = row.get("area_sqft", 0.0)
            area_acres = row.get("area_acres", 0.0)
            area_sqmi = row.get("area_sqmi", 0.0)
            elev_min = row.get("elev_min")
            elev_max = row.get("elev_max")
            elev_mean = row.get("elev_mean")
            lfp_length = row.get("lfp_length")
            lfp_slope = row.get("lfp_slope")

            # Select color based on watershed_id
            color = colors[watershed_id % len(colors)]

            geom = row.geometry
            if geom is None:
                continue

            # Handle both Polygon and MultiPolygon
            polygons = []
            if isinstance(geom, Polygon):
                polygons = [geom]
            elif isinstance(geom, MultiPolygon):
                polygons = list(geom.geoms)

            for polygon in polygons:
                # Add exterior ring as LWPOLYLINE
                exterior_coords = list(polygon.exterior.coords)
                if len(exterior_coords) >= 3:
                    # Remove closing point if present (ezdxf handles closing)
                    if exterior_coords[0] == exterior_coords[-1]:
                        exterior_coords = exterior_coords[:-1]

                    lwpoly = msp.add_lwpolyline(
                        exterior_coords,
                        dxfattribs={
                            "layer": "WATERSHEDS",
                            "color": color,
                        },
                        close=True,
                    )

                    # Attach XDATA
                    lwpoly.set_xdata(
                        "DELINEATOR",
                        [
                            (1000, f"watershed_id={watershed_id}"),
                            (1000, f"point_name={point_name}"),
                            (1040, area_sqft if area_sqft else 0.0),
                            (1040, area_acres if area_acres else 0.0),
                            (1040, area_sqmi if area_sqmi else 0.0),
                            (1040, elev_min if elev_min is not None else 0.0),
                            (1040, elev_max if elev_max is not None else 0.0),
                            (1040, elev_mean if elev_mean is not None else 0.0),
                            (1040, lfp_length if lfp_length is not None else 0.0),
                            (1040, lfp_slope if lfp_slope is not None else 0.0),
                        ],
                    )

                # Add interior rings (holes) as separate polylines
                for interior in polygon.interiors:
                    interior_coords = list(interior.coords)
                    if len(interior_coords) >= 3:
                        if interior_coords[0] == interior_coords[-1]:
                            interior_coords = interior_coords[:-1]

                        msp.add_lwpolyline(
                            interior_coords,
                            dxfattribs={
                                "layer": "WATERSHEDS_HOLES",
                                "color": 8,  # Gray
                            },
                            close=True,
                        )

        doc.saveas(output_path)
