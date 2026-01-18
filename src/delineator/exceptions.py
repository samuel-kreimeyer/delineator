"""Custom exceptions for the delineator package."""


class DelineatorError(Exception):
    """Base exception for all delineator errors."""

    pass


class LandXMLParseError(DelineatorError):
    """Error parsing LandXML file."""

    pass


class InvalidTINError(DelineatorError):
    """Invalid TIN surface data."""

    pass


class StudyPointError(DelineatorError):
    """Error with study points."""

    pass


class RasterizationError(DelineatorError):
    """Error during TIN rasterization."""

    pass


class HydrologyError(DelineatorError):
    """Error during hydrological processing."""

    pass


class WhiteBoxError(HydrologyError):
    """Error from WhiteBox Tools."""

    pass


class OutputError(DelineatorError):
    """Error generating output files."""

    pass


class ValidationError(DelineatorError):
    """Input validation error."""

    pass


class CRSError(DelineatorError):
    """Coordinate reference system error."""

    pass
