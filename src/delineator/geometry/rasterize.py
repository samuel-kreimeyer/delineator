"""TIN to raster interpolation."""

from pathlib import Path
from typing import Optional

import numpy as np
import matplotlib.tri as mtri
import rasterio
from rasterio.transform import from_bounds

from delineator.geometry.tin import TINSurface
from delineator.exceptions import RasterizationError
from delineator.utils.logging import get_logger


class TINRasterizer:
    """Rasterizes a TIN surface to a DEM raster."""

    def __init__(self, tin_surface: TINSurface, epsg_code: int = 0):
        """Initialize the rasterizer.

        Args:
            tin_surface: The TIN surface to rasterize.
            epsg_code: EPSG code for the coordinate reference system.
        """
        self.tin_surface = tin_surface
        self.epsg_code = epsg_code
        self.logger = get_logger(__name__)

        # Create triangulation
        self._create_triangulation()

    def _create_triangulation(self) -> None:
        """Create matplotlib triangulation from TIN surface."""
        x, y, z = self.tin_surface.get_point_arrays()
        triangles = self.tin_surface.get_triangle_array()

        self.triangulation = mtri.Triangulation(x, y, triangles)
        self.z_values = z
        self.interpolator = mtri.LinearTriInterpolator(self.triangulation, z)

    def compute_cell_size(self) -> float:
        """Compute an appropriate cell size based on average triangle edge length.

        Returns:
            Cell size as half the average edge length.
        """
        avg_edge = self.tin_surface.compute_average_edge_length()
        cell_size = avg_edge / 2.0
        self.logger.info(
            f"Average edge length: {avg_edge:.4f}, computed cell size: {cell_size:.4f}"
        )
        return cell_size

    def rasterize(
        self,
        output_path: Path,
        cell_size: Optional[float] = None,
        nodata: float = -9999.0,
    ) -> None:
        """Rasterize the TIN surface to a GeoTIFF.

        Args:
            output_path: Path for the output GeoTIFF.
            cell_size: Raster cell size. If None, auto-computed.
            nodata: NoData value for the raster.
        """
        if cell_size is None:
            cell_size = self.compute_cell_size()

        # Get bounds
        min_x, min_y, max_x, max_y = self.tin_surface.bounds

        # Compute raster dimensions
        width = int(np.ceil((max_x - min_x) / cell_size))
        height = int(np.ceil((max_y - min_y) / cell_size))

        if width <= 0 or height <= 0:
            raise RasterizationError(
                f"Invalid raster dimensions: {width}x{height}. "
                "Check TIN bounds and cell size."
            )

        self.logger.info(f"Creating raster: {width}x{height} cells")

        # Create coordinate grids
        x_coords = np.linspace(min_x + cell_size / 2, max_x - cell_size / 2, width)
        y_coords = np.linspace(max_y - cell_size / 2, min_y + cell_size / 2, height)
        xx, yy = np.meshgrid(x_coords, y_coords)

        # Interpolate elevation values
        self.logger.info("Interpolating elevation values...")
        zz = self.interpolator(xx, yy)

        # Handle masked values (outside triangulation)
        if hasattr(zz, "mask"):
            zz = zz.filled(nodata)
        else:
            zz = np.where(np.isnan(zz), nodata, zz)

        # Create the transform
        transform = from_bounds(min_x, min_y, max_x, max_y, width, height)

        # Set up CRS
        if self.epsg_code and self.epsg_code > 0:
            crs = f"EPSG:{self.epsg_code}"
        else:
            crs = None

        # Write raster
        self.logger.info(f"Writing raster to: {output_path}")
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with rasterio.open(
            output_path,
            "w",
            driver="GTiff",
            height=height,
            width=width,
            count=1,
            dtype=zz.dtype,
            crs=crs,
            transform=transform,
            nodata=nodata,
        ) as dst:
            dst.write(zz, 1)

        self.logger.info(f"Raster created: {output_path}")
