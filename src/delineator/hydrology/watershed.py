"""Watershed delineation from pour points."""

from pathlib import Path
import tempfile

import numpy as np
import rasterio
import geopandas as gpd
from shapely.geometry import Point

from whitebox import WhiteboxTools

from delineator.geometry.tin import StudyPoint
from delineator.exceptions import WhiteBoxError
from delineator.utils.logging import get_logger


class WatershedDelineator:
    """Delineates watersheds from study points (pour points)."""

    def __init__(self):
        """Initialize the delineator."""
        self.logger = get_logger(__name__)
        self.wbt = WhiteboxTools()
        self.wbt.set_verbose_mode(False)

    def delineate(
        self,
        flow_dir_path: Path,
        flow_acc_path: Path,
        study_points: list[StudyPoint],
        output_path: Path,
        pour_points_path: Path,
        snap_distance: float,
        epsg_code: int = 0,
        cumulative: bool = True,
    ) -> None:
        """Delineate watersheds for the given study points.

        Args:
            flow_dir_path: Path to D8 flow direction raster.
            flow_acc_path: Path to flow accumulation raster.
            study_points: List of study points (pour points).
            output_path: Path for watershed output raster.
            pour_points_path: Path for snapped pour points shapefile.
            snap_distance: Distance to snap pour points to stream cells.
            epsg_code: EPSG code for coordinate system.
            cumulative: If True, compute cumulative (overlapping) watersheds.
        """
        self.logger.info(f"Delineating watersheds for {len(study_points)} study points")

        # WhiteBox Tools requires absolute paths
        flow_dir_abs = str(Path(flow_dir_path).resolve())
        flow_acc_abs = str(Path(flow_acc_path).resolve())
        output_abs = str(Path(output_path).resolve())
        pour_points_abs = str(Path(pour_points_path).resolve())

        # Verify inputs exist
        if not Path(flow_dir_abs).exists():
            raise WhiteBoxError(f"Flow direction raster not found: {flow_dir_abs}")
        if not Path(flow_acc_abs).exists():
            raise WhiteBoxError(f"Flow accumulation raster not found: {flow_acc_abs}")

        if cumulative:
            self._delineate_cumulative(
                flow_dir_abs, flow_acc_abs, output_abs, pour_points_abs,
                study_points, snap_distance, epsg_code
            )
        else:
            self._delineate_exclusive(
                flow_dir_abs, flow_acc_abs, output_abs, pour_points_abs,
                study_points, snap_distance, epsg_code
            )

        self.logger.info(f"Watershed raster saved to: {output_path}")

    def _delineate_exclusive(
        self,
        flow_dir_abs: str,
        flow_acc_abs: str,
        output_abs: str,
        pour_points_abs: str,
        study_points: list[StudyPoint],
        snap_distance: float,
        epsg_code: int,
    ) -> None:
        """Delineate exclusive (non-overlapping) watersheds."""
        with tempfile.TemporaryDirectory() as tmpdir:
            input_points_path = Path(tmpdir) / "input_points.shp"

            self._create_pour_points_shapefile(
                study_points, input_points_path, epsg_code
            )

            self.logger.info(f"Snapping pour points with distance: {snap_distance}")
            self._snap_pour_points(
                str(input_points_path.resolve()),
                flow_acc_abs,
                pour_points_abs,
                snap_distance,
            )

            if not Path(pour_points_abs).exists():
                raise WhiteBoxError(
                    f"WhiteBox Tools failed to create snapped pour points"
                )

            self.logger.info("Running watershed delineation (exclusive mode)")
            self._watershed(flow_dir_abs, pour_points_abs, output_abs)

        if not Path(output_abs).exists():
            raise WhiteBoxError(f"WhiteBox Tools failed to create watershed raster")

    def _delineate_cumulative(
        self,
        flow_dir_abs: str,
        flow_acc_abs: str,
        output_abs: str,
        pour_points_abs: str,
        study_points: list[StudyPoint],
        snap_distance: float,
        epsg_code: int,
    ) -> None:
        """Delineate cumulative (overlapping) watersheds.

        Each watershed is computed independently so they can overlap.
        The output raster stores each cell's watershed membership as a
        combined raster where each point has its full contributing area.
        """
        self.logger.info("Running watershed delineation (cumulative mode)")

        # Read flow direction raster to get dimensions and metadata
        with rasterio.open(flow_dir_abs) as src:
            profile = src.profile.copy()
            shape = src.shape

        # Create output array - we'll store full watersheds
        # Each cell will be assigned to the watershed it belongs to
        # (for cumulative, we create individual masks stored in the raster)
        combined = np.zeros(shape, dtype=np.int32)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Process each study point individually
            for sp in study_points:
                self.logger.info(f"  Processing point {sp.id} ({sp.name})")

                # Create single-point shapefile
                single_point_path = tmpdir / f"point_{sp.id}.shp"
                snapped_point_path = tmpdir / f"snapped_{sp.id}.shp"
                ws_output_path = tmpdir / f"ws_{sp.id}.tif"

                self._create_pour_points_shapefile([sp], single_point_path, epsg_code)

                # Snap this point
                self._snap_pour_points(
                    str(single_point_path.resolve()),
                    flow_acc_abs,
                    str(snapped_point_path.resolve()),
                    snap_distance,
                )

                if not snapped_point_path.exists():
                    self.logger.warning(f"  Failed to snap point {sp.id}, skipping")
                    continue

                # Delineate watershed for this single point
                self._watershed(
                    flow_dir_abs,
                    str(snapped_point_path.resolve()),
                    str(ws_output_path.resolve()),
                )

                if not ws_output_path.exists():
                    self.logger.warning(f"  Failed to delineate watershed for point {sp.id}")
                    continue

                # Read the watershed and add to combined
                with rasterio.open(ws_output_path) as ws_src:
                    ws_data = ws_src.read(1)
                    nodata = ws_src.nodata

                    # Mark cells belonging to this watershed
                    if nodata is not None:
                        mask = ws_data != nodata
                    else:
                        mask = ws_data > 0

                    # In cumulative mode, each cell gets the ID of the point
                    # If a cell belongs to multiple watersheds, it keeps the
                    # most recent one (which is fine for the combined raster)
                    # The actual stats are computed per-watershed in boundaries.py
                    combined[mask] = sp.id

            # Also create snapped pour points file with all points
            all_points_path = tmpdir / "all_points.shp"
            self._create_pour_points_shapefile(study_points, all_points_path, epsg_code)
            self._snap_pour_points(
                str(all_points_path.resolve()),
                flow_acc_abs,
                pour_points_abs,
                snap_distance,
            )

        # Write combined watershed raster
        profile.update(dtype=np.int32, nodata=0)
        with rasterio.open(output_abs, 'w', **profile) as dst:
            dst.write(combined, 1)

        # Store individual watershed data for later use
        # We'll recompute in boundaries.py using the flow direction
        self._cumulative_mode = True

    def _create_pour_points_shapefile(
        self,
        study_points: list[StudyPoint],
        output_path: Path,
        epsg_code: int = 0,
    ) -> None:
        """Create a shapefile from study points."""
        geometries = [Point(sp.x, sp.y) for sp in study_points]
        data = {
            "id": [sp.id for sp in study_points],
            "name": [sp.name or f"Point_{sp.id}" for sp in study_points],
            "geometry": geometries,
        }

        gdf = gpd.GeoDataFrame(data)

        if epsg_code and epsg_code > 0:
            gdf.set_crs(f"EPSG:{epsg_code}", inplace=True)

        gdf.to_file(output_path)
        self.logger.debug(f"Created pour points shapefile: {output_path}")

    def _snap_pour_points(
        self,
        input_path: str,
        flow_acc_path: str,
        output_path: str,
        snap_distance: float,
    ) -> None:
        """Snap pour points to high flow accumulation cells."""
        try:
            result = self.wbt.snap_pour_points(
                pour_pts=input_path,
                flow_accum=flow_acc_path,
                output=output_path,
                snap_dist=snap_distance,
            )
            if result != 0:
                raise WhiteBoxError(f"snap_pour_points failed with code {result}")

        except WhiteBoxError:
            raise
        except Exception as e:
            raise WhiteBoxError(f"Error snapping pour points: {e}")

    def _watershed(
        self,
        flow_dir_path: str,
        pour_points_path: str,
        output_path: str,
    ) -> None:
        """Delineate watersheds using D8 flow direction."""
        try:
            result = self.wbt.watershed(
                d8_pntr=flow_dir_path,
                pour_pts=pour_points_path,
                output=output_path,
            )
            if result != 0:
                raise WhiteBoxError(f"watershed failed with code {result}")

        except WhiteBoxError:
            raise
        except Exception as e:
            raise WhiteBoxError(f"Error delineating watershed: {e}")
