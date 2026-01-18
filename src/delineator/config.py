"""Configuration dataclasses for the delineator package."""

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Optional


class DepressionMethod(str, Enum):
    """Method for handling depressions in the DEM."""

    FILL = "fill"
    BREACH = "breach"
    BREACH_LEAST_COST = "breach-least-cost"


class OutputFormat(str, Enum):
    """Output format for watershed boundaries."""

    SHAPEFILE = "shapefile"
    GEOJSON = "geojson"
    GEOPACKAGE = "geopackage"


class ReportFormat(str, Enum):
    """Output format for reports."""

    JSON = "json"
    CSV = "csv"
    TEXT = "text"
    ALL = "all"


@dataclass
class DelineatorConfig:
    """Configuration for watershed delineation."""

    input_file: Path
    output_dir: Path
    study_points_file: Optional[Path] = None
    use_cogo_points: bool = False
    cell_size: Optional[float] = None
    snap_distance: Optional[float] = None
    depression_method: DepressionMethod = DepressionMethod.BREACH
    output_format: OutputFormat = OutputFormat.SHAPEFILE
    report_formats: list[ReportFormat] = field(default_factory=lambda: [ReportFormat.JSON])
    keep_intermediates: bool = False
    epsg_code: Optional[int] = None
    verbose: bool = False
    cumulative: bool = True  # Compute cumulative (overlapping) watersheds by default

    def __post_init__(self) -> None:
        """Convert paths and validate configuration."""
        if isinstance(self.input_file, str):
            self.input_file = Path(self.input_file)
        if isinstance(self.output_dir, str):
            self.output_dir = Path(self.output_dir)
        if isinstance(self.study_points_file, str):
            self.study_points_file = Path(self.study_points_file)

        if not self.use_cogo_points and self.study_points_file is None:
            raise ValueError(
                "Either --use-cogo-points or --study-points must be specified"
            )

    @property
    def dem_path(self) -> Path:
        """Path to the raw DEM raster."""
        return self.output_dir / "dem.tif"

    @property
    def conditioned_dem_path(self) -> Path:
        """Path to the conditioned (depression-handled) DEM."""
        return self.output_dir / "conditioned_dem.tif"

    @property
    def flow_dir_path(self) -> Path:
        """Path to the flow direction raster."""
        return self.output_dir / "flow_dir.tif"

    @property
    def flow_acc_path(self) -> Path:
        """Path to the flow accumulation raster."""
        return self.output_dir / "flow_acc.tif"

    @property
    def watershed_path(self) -> Path:
        """Path to the watershed raster."""
        return self.output_dir / "watershed.tif"

    @property
    def pour_points_path(self) -> Path:
        """Path to the snapped pour points shapefile."""
        return self.output_dir / "pour_points.shp"

    @property
    def boundaries_path(self) -> Path:
        """Path to the watershed boundaries output."""
        ext_map = {
            OutputFormat.SHAPEFILE: ".shp",
            OutputFormat.GEOJSON: ".geojson",
            OutputFormat.GEOPACKAGE: ".gpkg",
        }
        return self.output_dir / f"watersheds{ext_map[self.output_format]}"

    def get_report_path(self, fmt: ReportFormat) -> Path:
        """Get the path for a report in the specified format."""
        ext_map = {
            ReportFormat.JSON: ".json",
            ReportFormat.CSV: ".csv",
            ReportFormat.TEXT: ".txt",
        }
        return self.output_dir / f"report{ext_map[fmt]}"

    def get_effective_snap_distance(self, cell_size: float) -> float:
        """Get the snap distance, defaulting to 2x cell size."""
        if self.snap_distance is not None:
            return self.snap_distance
        return cell_size * 2
