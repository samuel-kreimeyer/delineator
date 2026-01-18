"""Input validation utilities."""

from pathlib import Path

from delineator.exceptions import ValidationError


def validate_input_file(file_path: Path) -> None:
    """Validate that the input file exists and is a LandXML file.

    Args:
        file_path: Path to the input file.

    Raises:
        ValidationError: If the file is invalid.
    """
    if not file_path.exists():
        raise ValidationError(f"Input file not found: {file_path}")

    if not file_path.is_file():
        raise ValidationError(f"Input path is not a file: {file_path}")

    # Check for common LandXML extensions
    valid_extensions = {".xml", ".landxml"}
    if file_path.suffix.lower() not in valid_extensions:
        raise ValidationError(
            f"Input file must be a LandXML file (.xml or .landxml), "
            f"got: {file_path.suffix}"
        )


def validate_study_points(file_path: Path) -> None:
    """Validate that the study points file exists and has a supported format.

    Args:
        file_path: Path to the study points file.

    Raises:
        ValidationError: If the file is invalid.
    """
    if not file_path.exists():
        raise ValidationError(f"Study points file not found: {file_path}")

    if not file_path.is_file():
        raise ValidationError(f"Study points path is not a file: {file_path}")

    # Check for supported extensions
    valid_extensions = {".csv", ".geojson", ".json", ".shp", ".gpkg"}
    if file_path.suffix.lower() not in valid_extensions:
        raise ValidationError(
            f"Study points file must be CSV, GeoJSON, Shapefile, or GeoPackage, "
            f"got: {file_path.suffix}"
        )


def validate_positive_float(value: float, name: str) -> None:
    """Validate that a value is a positive float.

    Args:
        value: The value to validate.
        name: Name of the parameter for error messages.

    Raises:
        ValidationError: If the value is not positive.
    """
    if value <= 0:
        raise ValidationError(f"{name} must be positive, got: {value}")


def validate_epsg_code(epsg_code: int) -> None:
    """Validate that an EPSG code is reasonable.

    Args:
        epsg_code: The EPSG code to validate.

    Raises:
        ValidationError: If the EPSG code is invalid.
    """
    if epsg_code < 0:
        raise ValidationError(f"EPSG code must be non-negative, got: {epsg_code}")

    # Common EPSG code ranges
    # 1-32767: Standard codes
    # 100000+: User-defined codes
    if epsg_code > 0 and epsg_code < 1000:
        raise ValidationError(
            f"EPSG code {epsg_code} seems unusually low. "
            "Common codes are typically 4-5 digits (e.g., 4326, 32618)."
        )
