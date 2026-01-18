"""Watershed delineation automation from LandXML TIN surfaces using WhiteBox Tools."""

__version__ = "0.1.0"

from delineator.config import DelineatorConfig, DepressionMethod, OutputFormat, ReportFormat
from delineator.geometry.tin import Point3D, Triangle, TINSurface, StudyPoint
from delineator.parsers.landxml import LandXMLParser
from delineator.geometry.rasterize import TINRasterizer
from delineator.hydrology.preprocessing import DEMPreprocessor
from delineator.hydrology.flow import FlowModeler
from delineator.hydrology.watershed import WatershedDelineator
from delineator.output.boundaries import WatershedBoundaryGenerator
from delineator.output.reports import ReportGenerator

__all__ = [
    "__version__",
    "DelineatorConfig",
    "DepressionMethod",
    "OutputFormat",
    "ReportFormat",
    "Point3D",
    "Triangle",
    "TINSurface",
    "StudyPoint",
    "LandXMLParser",
    "TINRasterizer",
    "DEMPreprocessor",
    "FlowModeler",
    "WatershedDelineator",
    "WatershedBoundaryGenerator",
    "ReportGenerator",
]
