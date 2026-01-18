"""Flow direction and accumulation modeling."""

from pathlib import Path

from whitebox import WhiteboxTools

from delineator.exceptions import WhiteBoxError
from delineator.utils.logging import get_logger


class FlowModeler:
    """Computes flow direction and accumulation grids."""

    def __init__(self):
        """Initialize the flow modeler."""
        self.logger = get_logger(__name__)
        self.wbt = WhiteboxTools()
        self.wbt.set_verbose_mode(False)

    def compute_flow_direction(
        self,
        dem_path: Path,
        output_path: Path,
    ) -> None:
        """Compute D8 flow direction grid.

        Args:
            dem_path: Path to conditioned DEM raster.
            output_path: Path for flow direction output raster.
        """
        self.logger.info("Computing D8 flow direction")

        # WhiteBox Tools requires absolute paths
        dem_abs = str(Path(dem_path).resolve())
        output_abs = str(Path(output_path).resolve())

        # Verify input exists
        if not Path(dem_abs).exists():
            raise WhiteBoxError(f"Input DEM not found: {dem_abs}")

        try:
            result = self.wbt.d8_pointer(
                dem=dem_abs,
                output=output_abs,
            )
            if result != 0:
                raise WhiteBoxError(f"d8_pointer failed with code {result}")

        except WhiteBoxError:
            raise
        except Exception as e:
            raise WhiteBoxError(f"Error computing flow direction: {e}")

        # Verify output was created
        if not Path(output_abs).exists():
            raise WhiteBoxError(
                f"WhiteBox Tools failed to create flow direction: {output_path}"
            )

        self.logger.info(f"Flow direction saved to: {output_path}")

    def compute_flow_accumulation(
        self,
        dem_path: Path,
        output_path: Path,
        out_type: str = "cells",
    ) -> None:
        """Compute D8 flow accumulation grid.

        Args:
            dem_path: Path to conditioned DEM raster.
            output_path: Path for flow accumulation output raster.
            out_type: Output type - "cells" for cell count, "catchment area" for area.
        """
        self.logger.info("Computing D8 flow accumulation")

        # WhiteBox Tools requires absolute paths
        dem_abs = str(Path(dem_path).resolve())
        output_abs = str(Path(output_path).resolve())

        # Verify input exists
        if not Path(dem_abs).exists():
            raise WhiteBoxError(f"Input DEM not found: {dem_abs}")

        try:
            result = self.wbt.d8_flow_accumulation(
                i=dem_abs,
                output=output_abs,
                out_type=out_type,
            )
            if result != 0:
                raise WhiteBoxError(f"d8_flow_accumulation failed with code {result}")

        except WhiteBoxError:
            raise
        except Exception as e:
            raise WhiteBoxError(f"Error computing flow accumulation: {e}")

        # Verify output was created
        if not Path(output_abs).exists():
            raise WhiteBoxError(
                f"WhiteBox Tools failed to create flow accumulation: {output_path}"
            )

        self.logger.info(f"Flow accumulation saved to: {output_path}")

    def compute_all(
        self,
        dem_path: Path,
        flow_dir_path: Path,
        flow_acc_path: Path,
    ) -> None:
        """Compute both flow direction and accumulation.

        Args:
            dem_path: Path to conditioned DEM raster.
            flow_dir_path: Path for flow direction output raster.
            flow_acc_path: Path for flow accumulation output raster.
        """
        self.compute_flow_direction(dem_path, flow_dir_path)
        self.compute_flow_accumulation(dem_path, flow_acc_path)
