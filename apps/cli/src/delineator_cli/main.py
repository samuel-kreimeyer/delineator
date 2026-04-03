"""Command-line interface for the delineator package."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


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
        "-o",
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for results",
    )

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
        choices=["fill", "breach", "breach-least-cost"],
        default="breach",
        help="Depression handling method (default: breach)",
    )
    parser.add_argument(
        "--output-format",
        choices=["shapefile", "geojson", "geopackage", "dxf"],
        default="shapefile",
        help="Watershed boundary output format (default: shapefile)",
    )
    parser.add_argument(
        "--report-format",
        choices=["json", "csv", "text", "all"],
        default="json",
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
    parser.add_argument(
        "--epsg",
        type=int,
        dest="epsg_code",
        help="Override EPSG code for coordinate system",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose output",
    )

    return parser.parse_args(args)


def main(args: list[str] | None = None) -> int:
    """Main entry point for the CLI."""
    from delineator import DelineatorConfig, DepressionMethod, OutputFormat, ReportFormat
    from delineator.exceptions import DelineatorError
    from delineator.utils.logging import get_logger, setup_logging
    from delineator.utils.validation import validate_input_file, validate_study_points
    from delineator.workflow import run_delineation

    parsed_args = parse_args(args)
    setup_logging(verbose=parsed_args.verbose)
    logger = get_logger(__name__)

    try:
        validate_input_file(parsed_args.input_file)
        if parsed_args.study_points_file:
            validate_study_points(parsed_args.study_points_file)

        report_format = ReportFormat(parsed_args.report_format)
        if report_format == ReportFormat.ALL:
            report_formats = [ReportFormat.ALL]
        else:
            report_formats = [report_format]

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

        run_delineation(config)
        return 0
    except DelineatorError as exc:
        logger.error("Delineation error: %s", exc)
        return 1
    except Exception as exc:  # pragma: no cover
        logger.exception("Unexpected error: %s", exc)
        return 2


def entrypoint() -> None:
    """Console script entrypoint."""
    raise SystemExit(main())


if __name__ == "__main__":
    sys.exit(main())
