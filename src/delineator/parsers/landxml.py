"""LandXML file parser for TIN surfaces and COGO points."""

from pathlib import Path
from typing import Optional

from lxml import etree

from delineator.geometry.tin import Point3D, Triangle, TINSurface, StudyPoint
from delineator.exceptions import LandXMLParseError, InvalidTINError
from delineator.utils.logging import get_logger


class LandXMLParser:
    """Parser for LandXML files containing TIN surfaces."""

    def __init__(self, file_path: Path):
        """Initialize the parser.

        Args:
            file_path: Path to the LandXML file.
        """
        self.file_path = Path(file_path)
        self.logger = get_logger(__name__)
        self.epsg_code: Optional[int] = None
        self._parsed = False
        self._tin_surface: Optional[TINSurface] = None
        self._cogo_points: list[StudyPoint] = []

    def _parse_file(self) -> None:
        """Parse the XML file and extract relevant data using iterparse."""
        if self._parsed:
            return

        self._parsed = True
        tin_surface: Optional[TINSurface] = None
        cogo_points: list[StudyPoint] = []
        cogo_idx = 0

        surface_found = False
        in_surface = False
        in_definition = False
        in_pnts = False
        in_faces = False
        in_cgpoints = False

        try:
            context = etree.iterparse(
                str(self.file_path), events=("start", "end"), recover=False
            )
        except etree.XMLSyntaxError as e:
            raise LandXMLParseError(f"Invalid XML syntax: {e}")
        except Exception as e:
            raise LandXMLParseError(f"Error parsing file: {e}")

        def _local_name(tag: str) -> str:
            if tag.startswith("{"):
                return tag.split("}", 1)[1]
            return tag

        for event, elem in context:
            tag = _local_name(elem.tag)

            if event == "start":
                if tag == "Surface" and not surface_found:
                    surface_found = True
                    in_surface = True
                    surface_name = elem.get("name", "Unnamed")
                    self.logger.info(f"Parsing surface: {surface_name}")
                    tin_surface = TINSurface(name=surface_name)
                elif tag == "Surface" and surface_found:
                    in_surface = False

                if in_surface and tag == "Definition":
                    in_definition = True
                if in_surface and in_definition and tag == "Pnts":
                    in_pnts = True
                if in_surface and in_definition and tag == "Faces":
                    in_faces = True
                if tag == "CgPoints":
                    in_cgpoints = True

            elif event == "end":
                if tag == "CoordinateSystem":
                    epsg_attr = elem.get("epsgCode")
                    if epsg_attr and self.epsg_code is None:
                        try:
                            self.epsg_code = int(epsg_attr)
                            self.logger.info(f"Found EPSG code: {self.epsg_code}")
                        except ValueError:
                            self.logger.warning(f"Invalid EPSG code: {epsg_attr}")

                if in_pnts and tag == "P" and tin_surface is not None:
                    point_id = elem.get("id")
                    if point_id is not None:
                        try:
                            point_id = int(point_id)
                        except ValueError:
                            self.logger.warning(f"Invalid point ID: {point_id}")
                            point_id = None
                        if point_id is not None and elem.text:
                            parts = elem.text.strip().split()
                            if len(parts) < 3:
                                self.logger.warning(
                                    f"Point {point_id} has insufficient coordinates"
                                )
                            else:
                                try:
                                    northing = float(parts[0])
                                    easting = float(parts[1])
                                    elevation = float(parts[2])
                                    point = Point3D(
                                        x=easting, y=northing, z=elevation, id=point_id
                                    )
                                    tin_surface.add_point(point)
                                except ValueError as e:
                                    self.logger.warning(
                                        f"Error parsing point {point_id}: {e}"
                                    )

                if in_faces and tag == "F" and tin_surface is not None:
                    if elem.text:
                        parts = elem.text.strip().split()
                        if len(parts) >= 3:
                            try:
                                p1_id = int(parts[0])
                                p2_id = int(parts[1])
                                p3_id = int(parts[2])
                                if all(
                                    pid in tin_surface.points for pid in [p1_id, p2_id, p3_id]
                                ):
                                    triangle = Triangle(
                                        p1_id=p1_id, p2_id=p2_id, p3_id=p3_id
                                    )
                                    tin_surface.add_triangle(triangle)
                                else:
                                    self.logger.warning(
                                        "Triangle references missing points: "
                                        f"{p1_id}, {p2_id}, {p3_id}"
                                    )
                            except ValueError as e:
                                self.logger.warning(f"Error parsing face: {e}")

                if in_cgpoints and tag == "CgPoint":
                    content = elem.text
                    if content:
                        parts = content.strip().split()
                        if len(parts) >= 2:
                            try:
                                northing = float(parts[0])
                                easting = float(parts[1])
                                elevation = float(parts[2]) if len(parts) > 2 else None
                                cogo_idx += 1
                                study_point = StudyPoint(
                                    x=easting,
                                    y=northing,
                                    z=elevation,
                                    name=elem.get("name"),
                                    id=cogo_idx,
                                )
                                cogo_points.append(study_point)
                            except ValueError as e:
                                self.logger.warning(f"Error parsing COGO point: {e}")

                if tag == "Pnts":
                    in_pnts = False
                if tag == "Faces":
                    in_faces = False
                if tag == "Definition":
                    in_definition = False
                if tag == "CgPoints":
                    in_cgpoints = False
                if tag == "Surface" and in_surface:
                    in_surface = False

                elem.clear()
                while elem.getprevious() is not None:
                    del elem.getparent()[0]

        self._tin_surface = tin_surface
        self._cogo_points = cogo_points

    def parse(self) -> TINSurface:
        """Parse the LandXML file and extract the TIN surface.

        Returns:
            The parsed TIN surface.
        """
        self._parse_file()

        if self._tin_surface is None:
            raise LandXMLParseError("No Surface elements found in LandXML file")

        self.logger.info(
            f"Parsed TIN surface: {self._tin_surface.num_points} points, "
            f"{self._tin_surface.num_triangles} triangles"
        )

        if self._tin_surface.num_points == 0:
            raise InvalidTINError("TIN surface has no points")
        if self._tin_surface.num_triangles == 0:
            raise InvalidTINError("TIN surface has no triangles")

        return self._tin_surface

    def get_cogo_points(self) -> list[StudyPoint]:
        """Extract COGO points from the LandXML file.

        COGO points are in format: "northing easting [elevation]"

        Returns:
            List of study points extracted from CgPoints.
        """
        self._parse_file()

        if not self._cogo_points:
            self.logger.warning("No CgPoints element found in LandXML file")
            return []

        self.logger.info(f"Found {len(self._cogo_points)} COGO points")
        return list(self._cogo_points)
