"""LandXML file parser for TIN surfaces and COGO points."""

import re
from pathlib import Path
from typing import Optional

from lxml import etree

from delineator.geometry.tin import Point3D, Triangle, TINSurface, StudyPoint
from delineator.exceptions import LandXMLParseError, InvalidTINError
from delineator.utils.logging import get_logger


class LandXMLParser:
    """Parser for LandXML files containing TIN surfaces."""

    # Namespace patterns for different LandXML versions
    NAMESPACE_PATTERN = re.compile(r"\{(http://www\.landxml\.org/schema/LandXML-\d+\.\d+)\}")

    def __init__(self, file_path: Path):
        """Initialize the parser.

        Args:
            file_path: Path to the LandXML file.
        """
        self.file_path = Path(file_path)
        self.logger = get_logger(__name__)
        self._tree: Optional[etree._ElementTree] = None
        self._root: Optional[etree._Element] = None
        self._namespace: Optional[str] = None
        self._nsmap: dict[str, str] = {}
        self.epsg_code: Optional[int] = None

    def _parse_file(self) -> None:
        """Parse the XML file and extract namespace."""
        if self._tree is not None:
            return

        try:
            self._tree = etree.parse(str(self.file_path))
            self._root = self._tree.getroot()
        except etree.XMLSyntaxError as e:
            raise LandXMLParseError(f"Invalid XML syntax: {e}")
        except Exception as e:
            raise LandXMLParseError(f"Error parsing file: {e}")

        # Extract namespace
        if self._root.tag.startswith("{"):
            match = self.NAMESPACE_PATTERN.match(self._root.tag)
            if match:
                self._namespace = match.group(1)
                self._nsmap["lx"] = self._namespace
                self.logger.debug(f"Detected namespace: {self._namespace}")

        # Extract EPSG code if present
        self._extract_epsg_code()

    def _extract_epsg_code(self) -> None:
        """Extract EPSG code from CoordinateSystem element."""
        if self._root is None:
            return

        # Try with namespace
        if self._nsmap:
            coord_sys = self._root.find(".//lx:CoordinateSystem", self._nsmap)
        else:
            coord_sys = self._root.find(".//CoordinateSystem")

        if coord_sys is not None:
            epsg_attr = coord_sys.get("epsgCode")
            if epsg_attr:
                try:
                    self.epsg_code = int(epsg_attr)
                    self.logger.info(f"Found EPSG code: {self.epsg_code}")
                except ValueError:
                    self.logger.warning(f"Invalid EPSG code: {epsg_attr}")

    def _find_elements(self, xpath: str) -> list[etree._Element]:
        """Find elements using XPath, handling namespaces."""
        if self._root is None:
            return []

        if self._nsmap:
            return self._root.findall(xpath, self._nsmap)
        else:
            # Try without namespace prefix
            xpath_no_ns = xpath.replace("lx:", "")
            return self._root.findall(xpath_no_ns)

    def parse(self) -> TINSurface:
        """Parse the LandXML file and extract the TIN surface.

        Returns:
            The parsed TIN surface.
        """
        self._parse_file()

        # Find surface elements
        surfaces = self._find_elements(".//lx:Surface")
        if not surfaces:
            raise LandXMLParseError("No Surface elements found in LandXML file")

        # Use the first surface
        surface_elem = surfaces[0]
        surface_name = surface_elem.get("name", "Unnamed")
        self.logger.info(f"Parsing surface: {surface_name}")

        tin_surface = TINSurface(name=surface_name)

        # Parse points
        self._parse_points(surface_elem, tin_surface)

        # Parse faces
        self._parse_faces(surface_elem, tin_surface)

        self.logger.info(
            f"Parsed TIN surface: {tin_surface.num_points} points, "
            f"{tin_surface.num_triangles} triangles"
        )

        if tin_surface.num_points == 0:
            raise InvalidTINError("TIN surface has no points")
        if tin_surface.num_triangles == 0:
            raise InvalidTINError("TIN surface has no triangles")

        return tin_surface

    def _parse_points(
        self, surface_elem: etree._Element, tin_surface: TINSurface
    ) -> None:
        """Parse point elements from the surface.

        Points are in format: "northing easting elevation"
        """
        # Find Pnts element
        if self._nsmap:
            pnts_elem = surface_elem.find(".//lx:Definition/lx:Pnts", self._nsmap)
        else:
            pnts_elem = surface_elem.find(".//Definition/Pnts")

        if pnts_elem is None:
            raise LandXMLParseError("No Pnts element found in surface")

        # Parse each P element
        if self._nsmap:
            p_elements = pnts_elem.findall("lx:P", self._nsmap)
        else:
            p_elements = pnts_elem.findall("P")

        for p_elem in p_elements:
            point_id = p_elem.get("id")
            if point_id is None:
                continue

            try:
                point_id = int(point_id)
            except ValueError:
                self.logger.warning(f"Invalid point ID: {point_id}")
                continue

            content = p_elem.text
            if content is None:
                continue

            try:
                parts = content.strip().split()
                if len(parts) < 3:
                    self.logger.warning(f"Point {point_id} has insufficient coordinates")
                    continue

                # LandXML format: northing easting elevation
                northing = float(parts[0])
                easting = float(parts[1])
                elevation = float(parts[2])

                point = Point3D(x=easting, y=northing, z=elevation, id=point_id)
                tin_surface.add_point(point)

            except ValueError as e:
                self.logger.warning(f"Error parsing point {point_id}: {e}")

    def _parse_faces(
        self, surface_elem: etree._Element, tin_surface: TINSurface
    ) -> None:
        """Parse face (triangle) elements from the surface.

        Faces are in format: "p1_id p2_id p3_id"
        """
        # Find Faces element
        if self._nsmap:
            faces_elem = surface_elem.find(".//lx:Definition/lx:Faces", self._nsmap)
        else:
            faces_elem = surface_elem.find(".//Definition/Faces")

        if faces_elem is None:
            raise LandXMLParseError("No Faces element found in surface")

        # Parse each F element
        if self._nsmap:
            f_elements = faces_elem.findall("lx:F", self._nsmap)
        else:
            f_elements = faces_elem.findall("F")

        for f_elem in f_elements:
            content = f_elem.text
            if content is None:
                continue

            try:
                parts = content.strip().split()
                if len(parts) < 3:
                    continue

                p1_id = int(parts[0])
                p2_id = int(parts[1])
                p3_id = int(parts[2])

                # Verify points exist
                if all(pid in tin_surface.points for pid in [p1_id, p2_id, p3_id]):
                    triangle = Triangle(p1_id=p1_id, p2_id=p2_id, p3_id=p3_id)
                    tin_surface.add_triangle(triangle)
                else:
                    self.logger.warning(
                        f"Triangle references missing points: {p1_id}, {p2_id}, {p3_id}"
                    )

            except ValueError as e:
                self.logger.warning(f"Error parsing face: {e}")

    def get_cogo_points(self) -> list[StudyPoint]:
        """Extract COGO points from the LandXML file.

        COGO points are in format: "northing easting [elevation]"

        Returns:
            List of study points extracted from CgPoints.
        """
        self._parse_file()

        study_points: list[StudyPoint] = []

        # Find CgPoints element
        cg_points_elems = self._find_elements(".//lx:CgPoints")
        if not cg_points_elems:
            self.logger.warning("No CgPoints element found in LandXML file")
            return study_points

        for cg_points_elem in cg_points_elems:
            # Parse each CgPoint element
            if self._nsmap:
                cg_point_elements = cg_points_elem.findall("lx:CgPoint", self._nsmap)
            else:
                cg_point_elements = cg_points_elem.findall("CgPoint")

            for idx, cg_elem in enumerate(cg_point_elements):
                name = cg_elem.get("name")
                content = cg_elem.text
                if content is None:
                    continue

                try:
                    parts = content.strip().split()
                    if len(parts) < 2:
                        continue

                    # LandXML format: northing easting [elevation]
                    northing = float(parts[0])
                    easting = float(parts[1])
                    elevation = float(parts[2]) if len(parts) > 2 else None

                    study_point = StudyPoint(
                        x=easting,
                        y=northing,
                        z=elevation,
                        name=name,
                        id=idx + 1,
                    )
                    study_points.append(study_point)

                except ValueError as e:
                    self.logger.warning(f"Error parsing COGO point: {e}")

        self.logger.info(f"Found {len(study_points)} COGO points")
        return study_points
