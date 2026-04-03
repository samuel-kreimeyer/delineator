"""Workflow orchestration for watershed delineation."""

from delineator.config import DelineatorConfig, ReportFormat
from delineator.exceptions import DelineatorError
from delineator.geometry.rasterize import TINRasterizer
from delineator.hydrology.flow import FlowModeler
from delineator.hydrology.preprocessing import DEMPreprocessor
from delineator.hydrology.watershed import WatershedDelineator
from delineator.output.boundaries import WatershedBoundaryGenerator
from delineator.output.reports import ReportGenerator
from delineator.parsers.landxml import LandXMLParser
from delineator.parsers.study_points import StudyPointLoader
from delineator.utils.logging import get_logger


def run_delineation(config: DelineatorConfig) -> None:
    """Run the complete watershed delineation workflow."""
    logger = get_logger(__name__)

    config.output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Parsing LandXML file: %s", config.input_file)
    parser = LandXMLParser(config.input_file)
    tin_surface = parser.parse()

    epsg_code = config.epsg_code or parser.epsg_code
    if epsg_code is None:
        logger.warning("No EPSG code found, using 0 (undefined CRS)")
        epsg_code = 0

    if config.use_cogo_points:
        logger.info("Using COGO points from LandXML as study points")
        study_points = parser.get_cogo_points()
    else:
        logger.info("Loading study points from: %s", config.study_points_file)
        loader = StudyPointLoader(config.study_points_file)
        study_points = loader.load()

    if not study_points:
        raise DelineatorError("No study points found")

    logger.info("Found %s study points", len(study_points))

    logger.info("Rasterizing TIN surface to DEM")
    rasterizer = TINRasterizer(tin_surface, epsg_code=epsg_code)
    cell_size = config.cell_size or rasterizer.compute_cell_size()
    logger.info("Using cell size: %.4f", cell_size)
    rasterizer.rasterize(config.dem_path, cell_size=cell_size)

    logger.info("Conditioning DEM using method: %s", config.depression_method.value)
    preprocessor = DEMPreprocessor()
    preprocessor.condition_dem(
        config.dem_path,
        config.conditioned_dem_path,
        method=config.depression_method,
    )

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

    mode = "cumulative" if config.cumulative else "exclusive"
    logger.info("Delineating watersheds (mode: %s)", mode)
    snap_distance = config.get_effective_snap_distance(cell_size)
    logger.info("Using snap distance: %.4f", snap_distance)

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
        logger.info("Report saved: %s", report_path)

    if not config.keep_intermediates:
        logger.info("Cleaning up intermediate files")
        intermediate_files = [
            config.conditioned_dem_path,
            config.flow_dir_path,
            config.flow_acc_path,
            config.watershed_path,
        ]
        for file_path in intermediate_files:
            if file_path.exists():
                file_path.unlink()

    logger.info("Watershed delineation complete")
    logger.info("Boundaries saved to: %s", config.boundaries_path)
