"""DEM preprocessing for hydrological analysis."""

from pathlib import Path

from whitebox import WhiteboxTools

from delineator.config import DepressionMethod
from delineator.exceptions import WhiteBoxError
from delineator.utils.logging import get_logger


class DEMPreprocessor:
    """Preprocessor for conditioning DEMs for hydrological analysis."""

    def __init__(self):
        """Initialize the preprocessor."""
        self.logger = get_logger(__name__)
        self.wbt = WhiteboxTools()
        self.wbt.set_verbose_mode(False)

    def condition_dem(
        self,
        input_path: Path,
        output_path: Path,
        method: DepressionMethod = DepressionMethod.BREACH,
    ) -> None:
        """Condition the DEM by handling depressions.

        Args:
            input_path: Path to input DEM raster.
            output_path: Path for conditioned output DEM.
            method: Depression handling method.
        """
        self.logger.info(f"Conditioning DEM using method: {method.value}")

        # WhiteBox Tools requires absolute paths
        input_abs = str(Path(input_path).resolve())
        output_abs = str(Path(output_path).resolve())

        # Verify input exists
        if not Path(input_abs).exists():
            raise WhiteBoxError(f"Input DEM not found: {input_abs}")

        try:
            if method == DepressionMethod.FILL:
                self._fill_depressions(input_abs, output_abs)
            elif method == DepressionMethod.BREACH:
                self._breach_depressions(input_abs, output_abs)
            elif method == DepressionMethod.BREACH_LEAST_COST:
                self._breach_depressions_least_cost(input_abs, output_abs)
            else:
                raise ValueError(f"Unknown depression method: {method}")

        except WhiteBoxError:
            raise
        except Exception as e:
            raise WhiteBoxError(f"Error conditioning DEM: {e}")

        # Verify output was created
        if not Path(output_abs).exists():
            raise WhiteBoxError(
                f"WhiteBox Tools failed to create output file: {output_path}"
            )

        self.logger.info(f"Conditioned DEM saved to: {output_path}")

    def _fill_depressions(self, input_path: str, output_path: str) -> None:
        """Fill depressions using standard filling algorithm.

        Faster but may alter terrain features significantly.
        """
        self.logger.debug("Using fill_depressions with fix_flats=True")
        result = self.wbt.fill_depressions(
            dem=input_path,
            output=output_path,
            fix_flats=True,
        )
        if result != 0:
            raise WhiteBoxError(f"fill_depressions failed with code {result}")

    def _breach_depressions(self, input_path: str, output_path: str) -> None:
        """Breach depressions using standard breaching algorithm.

        Preserves terrain features better than filling.
        This is the default method.
        """
        self.logger.debug("Using breach_depressions")
        result = self.wbt.breach_depressions(
            dem=input_path,
            output=output_path,
        )
        if result != 0:
            raise WhiteBoxError(f"breach_depressions failed with code {result}")

    def _breach_depressions_least_cost(
        self, input_path: str, output_path: str
    ) -> None:
        """Breach depressions using least-cost algorithm.

        Most accurate but slower. Best for high-quality analysis.
        """
        self.logger.debug("Using breach_depressions_least_cost")
        result = self.wbt.breach_depressions_least_cost(
            dem=input_path,
            output=output_path,
            dist=0,  # No limit on breach distance
            fill=True,  # Fill remaining depressions
        )
        if result != 0:
            raise WhiteBoxError(f"breach_depressions_least_cost failed with code {result}")
