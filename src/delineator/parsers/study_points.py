"""Loaders for external study point files."""

import csv
import json
from pathlib import Path
from typing import Optional

import geopandas as gpd

from delineator.geometry.tin import StudyPoint
from delineator.exceptions import StudyPointError
from delineator.utils.logging import get_logger


class StudyPointLoader:
    """Loader for study points from various file formats."""

    SUPPORTED_EXTENSIONS = {".csv", ".geojson", ".json", ".shp", ".gpkg"}

    def __init__(self, file_path: Path):
        """Initialize the loader.

        Args:
            file_path: Path to the study points file.
        """
        self.file_path = Path(file_path)
        self.logger = get_logger(__name__)

        if not self.file_path.exists():
            raise StudyPointError(f"Study points file not found: {self.file_path}")

        ext = self.file_path.suffix.lower()
        if ext not in self.SUPPORTED_EXTENSIONS:
            raise StudyPointError(
                f"Unsupported file format: {ext}. "
                f"Supported formats: {', '.join(self.SUPPORTED_EXTENSIONS)}"
            )

    def load(self) -> list[StudyPoint]:
        """Load study points from the file.

        Returns:
            List of study points.
        """
        ext = self.file_path.suffix.lower()

        if ext == ".csv":
            return self._load_csv()
        elif ext in {".geojson", ".json"}:
            return self._load_geojson()
        elif ext in {".shp", ".gpkg"}:
            return self._load_vector()
        else:
            raise StudyPointError(f"Unsupported file format: {ext}")

    def _load_csv(self) -> list[StudyPoint]:
        """Load study points from a CSV file.

        Expected columns: x/easting/lon, y/northing/lat, z/elevation (optional), name (optional)
        """
        study_points: list[StudyPoint] = []

        with open(self.file_path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            fieldnames = [fn.lower() for fn in reader.fieldnames or []]

            # Determine column names
            x_col = self._find_column(fieldnames, ["x", "easting", "lon", "longitude", "e"])
            y_col = self._find_column(fieldnames, ["y", "northing", "lat", "latitude", "n"])
            z_col = self._find_column(
                fieldnames, ["z", "elevation", "elev", "height", "alt", "altitude"]
            )
            name_col = self._find_column(fieldnames, ["name", "id", "label", "point_name"])

            if x_col is None or y_col is None:
                raise StudyPointError(
                    "CSV file must have x/easting and y/northing columns"
                )

            # Map lowercase back to actual column names
            original_fieldnames = reader.fieldnames or []
            col_map = {fn.lower(): fn for fn in original_fieldnames}

            # Reset reader
            f.seek(0)
            reader = csv.DictReader(f)

            for idx, row in enumerate(reader):
                try:
                    x = float(row[col_map[x_col]])
                    y = float(row[col_map[y_col]])
                    z: Optional[float] = None
                    name: Optional[str] = None

                    if z_col and col_map.get(z_col):
                        z_val = row.get(col_map[z_col])
                        if z_val:
                            z = float(z_val)

                    if name_col and col_map.get(name_col):
                        name = row.get(col_map[name_col])

                    study_points.append(
                        StudyPoint(x=x, y=y, z=z, name=name, id=idx + 1)
                    )

                except (ValueError, KeyError) as e:
                    self.logger.warning(f"Error parsing row {idx + 1}: {e}")

        self.logger.info(f"Loaded {len(study_points)} study points from CSV")
        return study_points

    def _load_geojson(self) -> list[StudyPoint]:
        """Load study points from a GeoJSON file."""
        with open(self.file_path, encoding="utf-8") as f:
            data = json.load(f)

        study_points: list[StudyPoint] = []

        if data.get("type") == "FeatureCollection":
            features = data.get("features", [])
        elif data.get("type") == "Feature":
            features = [data]
        else:
            raise StudyPointError("Invalid GeoJSON: must be Feature or FeatureCollection")

        for idx, feature in enumerate(features):
            geometry = feature.get("geometry", {})
            properties = feature.get("properties", {})

            if geometry.get("type") != "Point":
                continue

            coords = geometry.get("coordinates", [])
            if len(coords) < 2:
                continue

            x = float(coords[0])
            y = float(coords[1])
            z = float(coords[2]) if len(coords) > 2 else None

            # Try to get name from properties
            name = (
                properties.get("name")
                or properties.get("id")
                or properties.get("label")
            )
            if name is not None:
                name = str(name)

            study_points.append(StudyPoint(x=x, y=y, z=z, name=name, id=idx + 1))

        self.logger.info(f"Loaded {len(study_points)} study points from GeoJSON")
        return study_points

    def _load_vector(self) -> list[StudyPoint]:
        """Load study points from a vector file (Shapefile or GeoPackage)."""
        try:
            gdf = gpd.read_file(self.file_path)
        except Exception as e:
            raise StudyPointError(f"Error reading vector file: {e}")

        study_points: list[StudyPoint] = []

        for idx, row in gdf.iterrows():
            geom = row.geometry
            if geom is None or geom.geom_type != "Point":
                continue

            x = geom.x
            y = geom.y
            z: Optional[float] = None

            # Try to get z from geometry or attributes
            if geom.has_z:
                z = geom.z
            else:
                for col in ["z", "elevation", "elev", "height"]:
                    if col in gdf.columns:
                        val = row.get(col)
                        if val is not None:
                            z = float(val)
                            break

            # Try to get name
            name: Optional[str] = None
            for col in ["name", "id", "label", "point_name"]:
                if col in gdf.columns:
                    val = row.get(col)
                    if val is not None:
                        name = str(val)
                        break

            study_points.append(StudyPoint(x=x, y=y, z=z, name=name, id=idx + 1))

        self.logger.info(
            f"Loaded {len(study_points)} study points from {self.file_path.suffix}"
        )
        return study_points

    @staticmethod
    def _find_column(fieldnames: list[str], candidates: list[str]) -> Optional[str]:
        """Find the first matching column name from candidates."""
        for candidate in candidates:
            if candidate in fieldnames:
                return candidate
        return None
