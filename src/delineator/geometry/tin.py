"""Data structures for TIN surfaces and study points."""

from dataclasses import dataclass, field
from typing import Optional
import numpy as np


@dataclass
class Point3D:
    """A 3D point with x (easting), y (northing), and z (elevation) coordinates."""

    x: float  # Easting
    y: float  # Northing
    z: float  # Elevation
    id: Optional[int] = None

    def __hash__(self) -> int:
        return hash((self.x, self.y, self.z, self.id))


@dataclass
class Triangle:
    """A triangle defined by three point indices."""

    p1_id: int
    p2_id: int
    p3_id: int

    def __hash__(self) -> int:
        return hash((self.p1_id, self.p2_id, self.p3_id))


@dataclass
class StudyPoint:
    """A study point (pour point) for watershed delineation."""

    x: float  # Easting
    y: float  # Northing
    z: Optional[float] = None  # Elevation (optional)
    name: Optional[str] = None
    id: Optional[int] = None

    def __hash__(self) -> int:
        return hash((self.x, self.y, self.z, self.name, self.id))


@dataclass
class TINSurface:
    """A Triangulated Irregular Network (TIN) surface."""

    name: str
    points: dict[int, Point3D] = field(default_factory=dict)
    triangles: list[Triangle] = field(default_factory=list)
    _point_arrays: Optional[tuple[np.ndarray, np.ndarray, np.ndarray]] = field(
        default=None, init=False, repr=False
    )
    _triangle_array: Optional[np.ndarray] = field(default=None, init=False, repr=False)

    def add_point(self, point: Point3D) -> None:
        """Add a point to the surface."""
        if point.id is None:
            point.id = len(self.points) + 1
        self.points[point.id] = point
        self._point_arrays = None
        self._triangle_array = None

    def add_triangle(self, triangle: Triangle) -> None:
        """Add a triangle to the surface."""
        self.triangles.append(triangle)
        self._triangle_array = None

    @property
    def num_points(self) -> int:
        """Number of points in the surface."""
        return len(self.points)

    @property
    def num_triangles(self) -> int:
        """Number of triangles in the surface."""
        return len(self.triangles)

    @property
    def bounds(self) -> tuple[float, float, float, float]:
        """Get the bounding box (min_x, min_y, max_x, max_y)."""
        if not self.points:
            raise ValueError("TIN surface has no points")

        xs = [p.x for p in self.points.values()]
        ys = [p.y for p in self.points.values()]

        return (min(xs), min(ys), max(xs), max(ys))

    @property
    def elevation_range(self) -> tuple[float, float]:
        """Get the elevation range (min_z, max_z)."""
        if not self.points:
            raise ValueError("TIN surface has no points")

        zs = [p.z for p in self.points.values()]
        return (min(zs), max(zs))

    def get_point_arrays(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Get numpy arrays of x, y, z coordinates."""
        if self._point_arrays is None:
            points_list = sorted(self.points.values(), key=lambda p: p.id)
            x = np.array([p.x for p in points_list])
            y = np.array([p.y for p in points_list])
            z = np.array([p.z for p in points_list])
            self._point_arrays = (x, y, z)
        return self._point_arrays

    def get_triangle_array(self) -> np.ndarray:
        """Get numpy array of triangle vertex indices (0-based)."""
        if self._triangle_array is None:
            # Create a mapping from point IDs to 0-based indices
            point_ids = sorted(self.points.keys())
            id_to_index = {pid: idx for idx, pid in enumerate(point_ids)}

            triangles = np.array([
                [
                    id_to_index[t.p1_id],
                    id_to_index[t.p2_id],
                    id_to_index[t.p3_id],
                ]
                for t in self.triangles
            ])
            self._triangle_array = triangles
        return self._triangle_array

    def compute_average_edge_length(self) -> float:
        """Compute the average edge length of all triangles."""
        if not self.triangles:
            raise ValueError("TIN surface has no triangles")

        edge_lengths = []
        for triangle in self.triangles:
            p1 = self.points[triangle.p1_id]
            p2 = self.points[triangle.p2_id]
            p3 = self.points[triangle.p3_id]

            # Compute edge lengths
            e1 = np.sqrt((p2.x - p1.x) ** 2 + (p2.y - p1.y) ** 2)
            e2 = np.sqrt((p3.x - p2.x) ** 2 + (p3.y - p2.y) ** 2)
            e3 = np.sqrt((p1.x - p3.x) ** 2 + (p1.y - p3.y) ** 2)

            edge_lengths.extend([e1, e2, e3])

        return np.mean(edge_lengths)
