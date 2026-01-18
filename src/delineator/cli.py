"""Command-line interface for the delineator package."""

import argparse
import sys
from pathlib import Path

from delineator.config import (
    DelineatorConfig,
    DepressionMethod,
    OutputFormat,
    ReportFormat,
)
from delineator.exceptions import DelineatorError
from delineator.utils.logging import setup_logging, get_logger
from delineator.utils.validation import validate_input_file, validate_study_points
from delineator.parsers.landxml import LandXMLParser
from delineator.parsers.study_points import StudyPointLoader
from delineator.geometry.rasterize import TINRasterizer
from delineator.hydrology.preprocessing import DEMPreprocessor
from delineator.hydrology.flow import FlowModeler
from delineator.hydrology.watershed import WatershedDelineator
from delineator.output.boundaries import WatershedBoundaryGenerator
from delineator.output.reports import ReportGenerator


def parse_args(args: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        prog="delineator",
        description="Watershed delineation from LandXML TIN surfaces using WhiteBox Tools",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  delineator input.xml -o output/ --use-cogo-points
  delineator input.xml -o output/ --study-points points.csv
  delineator input.xml -o output/ --use-cogo-points --depression-method breach-least-cost
  delineator input.xml -o output/ --use-cogo-points --output-format geojson --report-format all
        """,
    )

    parser.add_argument(
        "input_file",
        type=Path,
        help="Path to LandXML file containing TIN surface",
    )

    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        required=True,
        help="Output directory for results",
    )

    # Study points options (mutually exclusive group)
    points_group = parser.add_mutually_exclusive_group()
    points_group.add_argument(
        "--study-points",
        type=Path,
        dest="study_points_file",
        help="External study points file (CSV, GeoJSON, or Shapefile)",
    )
    points_group.add_argument(
        "--use-cogo-points",
        action="store_true",
        help="Use COGO points from LandXML as study points",
    )

    # Processing options
    parser.add_argument(
        "--cell-size",
        type=float,
        help="Raster cell size in map units (default: auto-computed from TIN)",
    )
    parser.add_argument(
        "--snap-distance",
        type=float,
        help="Pour point snap distance in map units (default: 2x cell size)",
    )
    parser.add_argument(
        "--depression-method",
        type=str,
        choices=[m.value for m in DepressionMethod],
        default=DepressionMethod.BREACH.value,
        help="Depression handling method (default: breach)",
    )

    # Output options
    parser.add_argument(
        "--output-format",
        type=str,
        choices=[f.value for f in OutputFormat],
        default=OutputFormat.SHAPEFILE.value,
        help="Watershed boundary output format (default: shapefile)",
    )
    parser.add_argument(
        "--report-format",
        type=str,
        choices=[f.value for f in ReportFormat],
        default=ReportFormat.JSON.value,
        help="Report output format (default: json)",
    )
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Keep intermediate raster files",
    )
    parser.add_argument(
        "--exclusive",
        action="store_true",
        help="Compute exclusive (non-overlapping) watersheds instead of cumulative",
    )

    # Other options
    parser.add_argument(
        "--epsg",
        type=int,
        dest="epsg_code",
        help="Override EPSG code for coordinate system",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output",
    )

    return parser.parse_args(args)


def run_delineation(config: DelineatorConfig) -> None:
    """Run the complete watershed delineation workflow."""
    logger = get_logger(__name__)

    # Create output directory
    config.output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Parse LandXML
    logger.info(f"Parsing LandXML file: {config.input_file}")
    parser = LandXMLParser(config.input_file)
    tin_surface = parser.parse()

    # Get EPSG code
    epsg_code = config.epsg_code or parser.epsg_code
    if epsg_code is None:
        logger.warning("No EPSG code found, using 0 (undefined CRS)")
        epsg_code = 0

    # Step 2: Get study points
    if config.use_cogo_points:
        logger.info("Using COGO points from LandXML as study points")
        study_points = parser.get_cogo_points()
    else:
        logger.info(f"Loading study points from: {config.study_points_file}")
        loader = StudyPointLoader(config.study_points_file)
        study_points = loader.load()

    if not study_points:
        raise DelineatorError("No study points found")

    logger.info(f"Found {len(study_points)} study points")

    # Step 3: Rasterize TIN to DEM
    logger.info("Rasterizing TIN surface to DEM")
    rasterizer = TINRasterizer(tin_surface, epsg_code=epsg_code)
    cell_size = config.cell_size or rasterizer.compute_cell_size()
    logger.info(f"Using cell size: {cell_size:.4f}")
    rasterizer.rasterize(config.dem_path, cell_size=cell_size)

    # Step 4: Condition DEM (handle depressions)
    logger.info(f"Conditioning DEM using method: {config.depression_method.value}")
    preprocessor = DEMPreprocessor()
    preprocessor.condition_dem(
        config.dem_path,
        config.conditioned_dem_path,
        method=config.depression_method,
    )

    # Step 5: Compute flow direction and accumulation
    logger.info("Computing flow direction and accumulation")
    flow_modeler = FlowModeler()
    flow_modeler.compute_flow_direction(
        config.conditioned_dem_path,
        config.flow_dir_path,
    )
    flow_modeler.compute_flow_accumulation(
        config.conditioned_dem_path,
        config.flow_acc_path,
    )

    # Step 6: Delineate watersheds
    mode = "cumulative" if config.cumulative else "exclusive"
    logger.info(f"Delineating watersheds (mode: {mode})")
    snap_distance = config.get_effective_snap_distance(cell_size)
    logger.info(f"Using snap distance: {snap_distance:.4f}")

    delineator = WatershedDelineator()
    delineator.delineate(
        flow_dir_path=config.flow_dir_path,
        flow_acc_path=config.flow_acc_path,
        study_points=study_points,
        output_path=config.watershed_path,
        pour_points_path=config.pour_points_path,
        snap_distance=snap_distance,
        epsg_code=epsg_code,
        cumulative=config.cumulative,
    )

    # Step 7: Generate watershed boundaries
    logger.info("Generating watershed boundaries")
    boundary_generator = WatershedBoundaryGenerator()
    watersheds_gdf = boundary_generator.generate(
        watershed_raster_path=config.watershed_path,
        dem_path=config.dem_path,
        study_points=study_points,
        output_path=config.boundaries_path,
        output_format=config.output_format,
        flow_dir_path=config.flow_dir_path,
        flow_acc_path=config.flow_acc_path,
        pour_points_path=config.pour_points_path,
        cumulative=config.cumulative,
    )

    # Step 8: Generate reports
    logger.info("Generating reports")
    report_generator = ReportGenerator()

    if config.report_formats == [ReportFormat.ALL]:
        formats = [ReportFormat.JSON, ReportFormat.CSV, ReportFormat.TEXT]
    else:
        formats = config.report_formats

    for fmt in formats:
        report_path = config.get_report_path(fmt)
        report_generator.generate(
            watersheds_gdf=watersheds_gdf,
            output_path=report_path,
            report_format=fmt,
        )
        logger.info(f"Report saved: {report_path}")

    # Step 9: Clean up intermediate files if requested
    if not config.keep_intermediates:
        logger.info("Cleaning up intermediate files")
        intermediate_files = [
            config.conditioned_dem_path,
            config.flow_dir_path,
            config.flow_acc_path,
            config.watershed_path,
        ]
        for f in intermediate_files:
            if f.exists():
                f.unlink()

    logger.info("Watershed delineation complete!")
    logger.info(f"Boundaries saved to: {config.boundaries_path}")


def main(args: list[str] | None = None) -> int:
    """Main entry point for the CLI."""
    parsed_args = parse_args(args)

    # Set up logging
    setup_logging(verbose=parsed_args.verbose)
    logger = get_logger(__name__)

    try:
        # Validate inputs
        validate_input_file(parsed_args.input_file)
        if parsed_args.study_points_file:
            validate_study_points(parsed_args.study_points_file)

        # Parse report formats
        report_format = ReportFormat(parsed_args.report_format)
        if report_format == ReportFormat.ALL:
            report_formats = [ReportFormat.ALL]
        else:
            report_formats = [report_format]

        # Create configuration
        config = DelineatorConfig(
            input_file=parsed_args.input_file,
            output_dir=parsed_args.output_dir,
            study_points_file=parsed_args.study_points_file,
            use_cogo_points=parsed_args.use_cogo_points,
            cell_size=parsed_args.cell_size,
            snap_distance=parsed_args.snap_distance,
            depression_method=DepressionMethod(parsed_args.depression_method),
            output_format=OutputFormat(parsed_args.output_format),
            report_formats=report_formats,
            keep_intermediates=parsed_args.keep_intermediates,
            epsg_code=parsed_args.epsg_code,
            verbose=parsed_args.verbose,
            cumulative=not parsed_args.exclusive,
        )

        # Run delineation
        run_delineation(config)
        return 0

    except DelineatorError as e:
        logger.error(f"Delineation error: {e}")
        return 1
    except Exception as e:
        logger.exception(f"Unexpected error: {e}")
        return 2


if __name__ == "__main__":
    sys.exit(main())
