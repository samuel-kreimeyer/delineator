# Getting Started with Delineator

Delineator is a command-line tool that automates watershed delineation from LandXML TIN surfaces using WhiteBox Tools.

## Installation

### Prerequisites

- Python 3.9 or higher
- pip package manager

### Install from Source

```bash
git clone https://github.com/example/delineator.git
cd delineator
pip install -e .
```

### Dependencies

The following packages will be installed automatically:

- `whitebox` - WhiteBox Tools for hydrological analysis
- `numpy` - Numerical operations
- `matplotlib` - TIN interpolation
- `rasterio` - Raster I/O (GeoTIFF)
- `geopandas` - Vector I/O and spatial operations
- `shapely` - Geometry operations
- `pyproj` - Coordinate transformations
- `lxml` - XML parsing

## Quick Start

### Basic Usage

The simplest way to run delineator is with a LandXML file containing both a TIN surface and COGO points:

```bash
delineator input.xml -o output/ --use-cogo-points
```

This will:
1. Parse the TIN surface from the LandXML file
2. Extract COGO points as study (pour) points
3. Rasterize the TIN to a DEM
4. Compute flow direction and accumulation
5. Delineate watersheds for each pour point
6. Generate watershed boundary shapefiles and reports

### Output Files

After running, you'll find these files in your output directory:

| File | Description |
|------|-------------|
| `dem.tif` | Raw DEM raster (if `--keep-intermediates`) |
| `watersheds.shp` | Watershed boundary polygons |
| `pour_points.shp` | Snapped pour point locations |
| `report.json` | Watershed statistics (area, elevation) |

### Example with Options

```bash
delineator input.xml -o output/ \
    --use-cogo-points \
    --cell-size 1 \
    --snap-distance 10 \
    --output-format geojson \
    --report-format all \
    --keep-intermediates
```

## Input Requirements

### LandXML File

Your LandXML file must contain:

1. **A TIN Surface** with:
   - `Surface/Definition/Pnts/P` elements containing point coordinates
   - `Surface/Definition/Faces/F` elements defining triangles
   - Format: `<P id="1">northing easting elevation</P>`

2. **COGO Points** (if using `--use-cogo-points`):
   - `CgPoints/CgPoint` elements
   - Format: `<CgPoint name="outlet">northing easting [elevation]</CgPoint>`

3. **Coordinate System** (optional but recommended):
   - `CoordinateSystem` element with `epsgCode` attribute

### External Study Points

Instead of COGO points, you can provide study points in:
- CSV (columns: x/easting, y/northing, optionally z/elevation, name)
- GeoJSON (Point features)
- Shapefile (.shp)
- GeoPackage (.gpkg)

```bash
delineator input.xml -o output/ --study-points points.csv
```

## Next Steps

- See [Configuration Guide](configuration.md) for all available options
- See [Understanding Results](understanding-results.md) for interpreting output
- See [Troubleshooting](troubleshooting.md) for common issues
