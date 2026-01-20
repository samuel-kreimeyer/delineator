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
                flow_dir_path=None,  # Not available in exclusive mode
                pour_point_coords=None,
            )

        # Add flow path columns
        lfp_lengths = []
        lfp_slopes = []
        lfp_elev_starts = []
        lfp_elev_ends = []
        lfp_elev_drops = []
        for _, row in gdf.iterrows():
            ws_id = int(row["watershed_id"])
            if ws_id in lfp_stats:
                length, slope, elev_start, elev_end, elev_drop = lfp_stats[ws_id]
                lfp_lengths.append(length)
                lfp_slopes.append(slope)
                lfp_elev_starts.append(elev_start)
                lfp_elev_ends.append(elev_end)
                lfp_elev_drops.append(elev_drop)
            else:
                lfp_lengths.append(None)
                lfp_slopes.append(None)
                lfp_elev_starts.append(None)
                lfp_elev_ends.append(None)
                lfp_elev_drops.append(None)

        gdf["lfp_length"] = lfp_lengths
        gdf["lfp_slope"] = lfp_slopes
        gdf["lfp_elev_start"] = lfp_elev_starts
        gdf["lfp_elev_end"] = lfp_elev_ends
        gdf["lfp_elev_drop"] = lfp_elev_drops

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
        lfp_elev_starts = []
        lfp_elev_ends = []
        lfp_elev_drops = []

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            for idx, row in pour_points.iterrows():
                pt_id = row['id']
                pt_name = row.get('name', f"Point_{pt_id}")
                pt_geom = row.geometry

                self.logger.info(f"  Processing watershed for point {pt_id} ({pt_name})")

                # Create single-point shapefile
                single_pt_path = tmpdir / f"pt_{pt_id}.shp"
                single_gdf = gpd.GeoDataFrame(
                    {"id": [pt_id], "name": [pt_name], "geometry": [pt_geom]},
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
                # Pass the pour point coordinates for better accuracy
                pour_coords = (pt_geom.x, pt_geom.y) if pt_geom else None
                lfp_stats = self._compute_longest_flowpath_stats(
                    dem_path=dem_path,
                    basins_path=ws_output,
                    work_dir=tmpdir,
                    flow_dir_path=flow_dir_path,
                    pour_point_coords=pour_coords,
                )
                # Get stats for this watershed (ID will be 1 in single-watershed raster)
                if lfp_stats:
                    # The single watershed will have value from the raster
                    stats_tuple = next(iter(lfp_stats.values()), (None, None, None, None, None))
                    lfp_length, lfp_slope, lfp_elev_start, lfp_elev_end, lfp_elev_drop = stats_tuple
                else:
                    lfp_length, lfp_slope, lfp_elev_start, lfp_elev_end, lfp_elev_drop = None, None, None, None, None

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
                        lfp_elev_starts.append(lfp_elev_start)
                        lfp_elev_ends.append(lfp_elev_end)
                        lfp_elev_drops.append(lfp_elev_drop)
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
            "lfp_elev_start": lfp_elev_starts,
            "lfp_elev_end": lfp_elev_ends,
            "lfp_elev_drop": lfp_elev_drops,
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
        flow_dir_path: Optional[Path] = None,
        pour_point_coords: Optional[tuple[float, float]] = None,
    ) -> dict[int, tuple[float, float, float, float, float]]:
        """Compute longest flow path statistics for each watershed.

        Walks the D8 flow direction grid backwards from the pour point to find
        the furthest cell (longest hydraulic path). This represents the maximum
        distance a drop of rain would travel to reach the pour point.

        Args:
            dem_path: Path to DEM raster.
            basins_path: Path to watershed/basins raster.
            work_dir: Working directory for temporary outputs.
            flow_dir_path: Path to D8 flow direction raster (optional, for better accuracy).
            pour_point_coords: (x, y) coordinates of pour point in CRS units.

        Returns:
            Dict mapping watershed_id to (path_length_ft, avg_slope_ft_per_ft,
                                          elev_start_ft, elev_end_ft, elev_drop_ft).
        """
        self.logger.info("Computing longest flow path statistics via D8 backwards walk")

        stats = {}

        with rasterio.open(dem_path) as dem_src:
            dem_data = dem_src.read(1)
            dem_nodata = dem_src.nodata
            dem_transform = dem_src.transform
            cell_size = dem_src.res[0]  # Assumed square cells

        with rasterio.open(basins_path) as basin_src:
            basin_data = basin_src.read(1)
            basin_nodata = basin_src.nodata
            basin_transform = basin_src.transform

        # Load flow direction if provided
        flow_dir_data = None
        if flow_dir_path and flow_dir_path.exists():
            with rasterio.open(flow_dir_path) as flow_src:
                flow_dir_data = flow_src.read(1)
                flow_dir_nodata = flow_src.nodata

        # Get unique watershed IDs
        unique_basins = np.unique(basin_data)
        if basin_nodata is not None:
            unique_basins = unique_basins[unique_basins != basin_nodata]
        unique_basins = unique_basins[unique_basins > 0]

        for ws_id in unique_basins:
            ws_id = int(ws_id)

            # Create mask for this watershed
            ws_mask = basin_data == ws_id
            if not np.any(ws_mask):
                continue

            # Find pour point for this watershed
            # The pour point is typically the cell with the lowest elevation in the watershed
            # or the provided coordinates if available
            pour_r, pour_c = None, None

            if pour_point_coords and len(unique_basins) == 1:
                # Single watershed - use provided pour point coordinates
                pour_c, pour_r = ~basin_transform * pour_point_coords
                pour_r, pour_c = int(round(pour_r)), int(round(pour_c))
                if not (0 <= pour_r < basin_data.shape[0] and 0 <= pour_c < basin_data.shape[1]):
                    pour_r, pour_c = None, None
                elif basin_data[pour_r, pour_c] != ws_id:
                    pour_r, pour_c = None, None

            if pour_r is None or pour_c is None:
                # Find the lowest elevation cell in the watershed as pour point
                ws_elevs = np.where(ws_mask, dem_data, np.inf)
                min_idx = np.argmin(ws_elevs)
                pour_r, pour_c = np.unravel_index(min_idx, ws_elevs.shape)

            # Compute flow distances from all cells in watershed to pour point
            # using D8 backwards walk
            if flow_dir_data is not None:
                distances, elevations = self._compute_d8_distances(
                    flow_dir_data, dem_data, ws_mask, pour_r, pour_c,
                    cell_size, dem_nodata, flow_dir_nodata
                )
            else:
                # Fallback: use simple distance calculation
                self.logger.warning(
                    f"No flow direction provided for watershed {ws_id}, "
                    "using Euclidean distance approximation"
                )
                distances = self._compute_euclidean_distances(
                    ws_mask, pour_r, pour_c, cell_size
                )
                elevations = dem_data.copy()

            if distances is None or not np.any(distances > 0):
                continue

            # Find the cell with maximum flow distance (furthest from pour point)
            max_dist = np.max(distances)
            max_idx = np.argmax(distances)
            max_r, max_c = np.unravel_index(max_idx, distances.shape)

            # Get elevations
            elev_start = float(elevations[max_r, max_c])  # Elevation at furthest point
            elev_end = float(elevations[pour_r, pour_c])  # Elevation at pour point
            elev_drop = elev_start - elev_end  # Total elevation drop

            # Compute average slope
            if max_dist > 0:
                avg_slope = abs(elev_drop) / max_dist
            else:
                avg_slope = 0.0

            # Sanity check: minimum distance should be roughly sqrt(area)
            area_sqft = np.sum(ws_mask) * (cell_size ** 2)
            min_expected_dist = np.sqrt(area_sqft)

            if max_dist < min_expected_dist * 0.5:
                self.logger.warning(
                    f"Watershed {ws_id}: Computed hydraulic length {max_dist:.2f} ft "
                    f"is less than expected minimum {min_expected_dist:.2f} ft "
                    f"(sqrt of area {area_sqft:.2f} ft²)"
                )

            stats[ws_id] = (max_dist, avg_slope, elev_start, elev_end, elev_drop)
            self.logger.debug(
                f"Watershed {ws_id}: length={max_dist:.2f} ft, "
                f"slope={avg_slope:.6f}, drop={elev_drop:.2f} ft"
            )

        self.logger.info(f"Computed flow path stats for {len(stats)} watersheds")
        return stats

    def _compute_d8_distances(
        self,
        flow_dir: np.ndarray,
        dem: np.ndarray,
        ws_mask: np.ndarray,
        pour_r: int,
        pour_c: int,
        cell_size: float,
        dem_nodata: Optional[float],
        flow_nodata: Optional[float],
    ) -> tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        """Compute flow distances by walking D8 grid backwards from pour point.

        D8 encoding (WhiteboxTools):
        1=E, 2=SE, 4=S, 8=SW, 16=W, 32=NW, 64=N, 128=NE

        Args:
            flow_dir: D8 flow direction array.
            dem: DEM elevation array.
            ws_mask: Boolean mask for watershed.
            pour_r, pour_c: Row, col of pour point.
            cell_size: Cell size in feet.
            dem_nodata: DEM nodata value.
            flow_nodata: Flow direction nodata value.

        Returns:
            (distances, elevations) arrays, or (None, None) on error.
        """
        rows, cols = flow_dir.shape

        # D8 direction mappings: value -> (dr, dc)
        d8_map = {
            1: (0, 1),      # E
            2: (1, 1),      # SE
            4: (1, 0),      # S
            8: (1, -1),     # SW
            16: (0, -1),    # W
            32: (-1, -1),   # NW
            64: (-1, 0),    # N
            128: (-1, 1),   # NE
        }

        # Reverse mapping: to find cells that flow TO (r, c), we need cells
        # at (r + dr, c + dc) that have direction pointing to (r, c)
        # This means we need the opposite direction
        reverse_map = {
            1: (0, -1),     # Cell to the W flows E to current cell
            2: (-1, -1),    # Cell to the NW flows SE to current cell
            4: (-1, 0),     # Cell to the N flows S to current cell
            8: (-1, 1),     # Cell to the NE flows SW to current cell
            16: (0, 1),     # Cell to the E flows W to current cell
            32: (1, 1),     # Cell to the SE flows NW to current cell
            64: (1, 0),     # Cell to the S flows N to current cell
            128: (1, -1),   # Cell to the SW flows NE to current cell
        }

        # Initialize distance array
        distances = np.zeros((rows, cols), dtype=np.float32)
        visited = np.zeros((rows, cols), dtype=bool)

        # BFS from pour point backwards through the flow network
        from collections import deque
        queue = deque([(pour_r, pour_c, 0.0)])
        visited[pour_r, pour_c] = True

        while queue:
            r, c, dist = queue.popleft()
            distances[r, c] = dist

            # Find all cells that flow TO (r, c)
            for d8_val, (dr, dc) in reverse_map.items():
                nr, nc = r + dr, c + dc

                # Check bounds
                if not (0 <= nr < rows and 0 <= nc < cols):
                    continue

                # Check if already visited
                if visited[nr, nc]:
                    continue

                # Check if in watershed
                if not ws_mask[nr, nc]:
                    continue

                # Check if this cell flows to (r, c)
                cell_flow_dir = flow_dir[nr, nc]
                if flow_nodata is not None and cell_flow_dir == flow_nodata:
                    continue

                if cell_flow_dir != d8_val:
                    continue

                # Compute distance increment
                # Diagonal cells are sqrt(2) * cell_size apart
                if abs(dr) + abs(dc) == 2:  # Diagonal
                    step_dist = cell_size * np.sqrt(2)
                else:  # Cardinal
                    step_dist = cell_size

                new_dist = dist + step_dist

                # Add to queue
                queue.append((nr, nc, new_dist))
                visited[nr, nc] = True

        return distances, dem

    def _compute_euclidean_distances(
        self,
        ws_mask: np.ndarray,
        pour_r: int,
        pour_c: int,
        cell_size: float,
    ) -> np.ndarray:
        """Compute Euclidean distances from pour point (fallback method).

        Args:
            ws_mask: Boolean mask for watershed.
            pour_r, pour_c: Row, col of pour point.
            cell_size: Cell size in feet.

        Returns:
            Distance array.
        """
        rows, cols = ws_mask.shape

        # Create coordinate grids
        r_coords = np.arange(rows)[:, np.newaxis] - pour_r
        c_coords = np.arange(cols)[np.newaxis, :] - pour_c

        # Compute Euclidean distances
        distances = np.sqrt(r_coords**2 + c_coords**2) * cell_size

        # Mask to watershed only
        distances = np.where(ws_mask, distances, 0)

        return distances

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
            lfp_elev_start = row.get("lfp_elev_start")
            lfp_elev_end = row.get("lfp_elev_end")
            lfp_elev_drop = row.get("lfp_elev_drop")

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
                            (1040, lfp_elev_start if lfp_elev_start is not None else 0.0),
                            (1040, lfp_elev_end if lfp_elev_end is not None else 0.0),
                            (1040, lfp_elev_drop if lfp_elev_drop is not None else 0.0),
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
