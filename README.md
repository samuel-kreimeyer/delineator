# Delineator

Automated watershed delineation from LandXML TIN surfaces using WhiteBox Tools.

## Overview

Delineator is a Python CLI tool that automates the process of delineating watersheds from civil engineering survey data. It takes a LandXML file containing a TIN (Triangulated Irregular Network) surface and study points, then produces watershed boundaries with area and elevation statistics.

## Features

- Parse LandXML TIN surfaces and COGO points
- Rasterize TIN to DEM with configurable resolution
- Hydrologically condition DEMs (breach/fill depressions)
- Compute D8 flow direction and accumulation
- Delineate cumulative or exclusive watersheds
- Generate vector boundaries (Shapefile, GeoJSON, GeoPackage, DXF)
- Compute longest flow path length and slope for each watershed
- Produce statistical reports (JSON, CSV, text)

## Installation

### Requirements

- Python 3.9+
- pip

### Install from Source

```bash
git clone https://github.com/example/delineator.git
cd delineator
pip install -e .
```

## Quick Start

```bash
# Basic usage with COGO points from LandXML
delineator input.xml -o output/ --use-cogo-points

# With custom resolution and snap distance
delineator input.xml -o output/ --use-cogo-points --cell-size 1 --snap-distance 10

# Using external study points
delineator input.xml -o output/ --study-points points.csv

# Full options
delineator input.xml -o output/ \
    --use-cogo-points \
    --cell-size 1 \
    --snap-distance 10 \
    --depression-method breach \
    --output-format geopackage \
    --report-format all \
    --keep-intermediates
```

## Command-Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `-o, --output-dir` | Output directory (required) | - |
| `--use-cogo-points` | Use COGO points from LandXML | - |
| `--study-points FILE` | External study points file | - |
| `--cell-size N` | Raster cell size in map units | Auto |
| `--snap-distance N` | Pour point snap distance | 2x cell size |
| `--depression-method` | `fill`, `breach`, or `breach-least-cost` | `breach` |
| `--output-format` | `shapefile`, `geojson`, `geopackage`, or `dxf` | `shapefile` |
| `--report-format` | `json`, `csv`, `text`, or `all` | `json` |
| `--exclusive` | Non-overlapping watersheds | Cumulative |
| `--keep-intermediates` | Keep intermediate rasters | No |
| `--epsg CODE` | Override EPSG code | From LandXML |
| `-v, --verbose` | Verbose output | No |

## Output Files

| File | Description |
|------|-------------|
| `watersheds.shp` | Watershed boundary polygons |
| `pour_points.shp` | Snapped pour point locations |
| `report.json` | Watershed statistics |
| `dem.tif` | DEM raster (with `--keep-intermediates`) |

## Watershed Statistics

Each watershed includes:
- Area in square feet, acres, and square miles
- Minimum, maximum, and mean elevation (feet)
- Longest flow path length (feet) and average slope (ft/ft)

## Processing Pipeline

```
LandXML → TIN Parser → Rasterizer → DEM
                                     ↓
                          Depression Handler
                                     ↓
                           Flow Direction
                                     ↓
                          Flow Accumulation
                                     ↓
               Study Points → Watershed Delineation
                                     ↓
                    Boundaries ← Vectorization
                                     ↓
                              Statistics Report
```

## Documentation

- [Getting Started](docs/user/getting-started.md)
- [Configuration Guide](docs/user/configuration.md)
- [Understanding Results](docs/user/understanding-results.md)
- [Troubleshooting](docs/user/troubleshooting.md)

## Tips for Best Results

### Grid Resolution

Use finer resolution (`--cell-size 1`) when:
- Delineating small watersheds
- Pour points are in narrow ditches or channels
- High accuracy is needed

Use coarser resolution for faster processing during testing.

### Snap Distance

Increase snap distance when:
- Watersheds have unexpectedly small areas
- Pour points are near but not exactly on channels

### Cumulative vs Exclusive

- **Cumulative** (default): Each watershed shows total contributing area. Watersheds can overlap if pour points are nested.
- **Exclusive** (`--exclusive`): Each cell belongs to one watershed only. Use when non-overlapping boundaries are required.

## Dependencies

- whitebox - WhiteBox Tools for hydrology
- numpy - Numerical operations
- matplotlib - TIN interpolation
- rasterio - Raster I/O
- geopandas - Vector I/O
- shapely - Geometry operations
- pyproj - Coordinate systems
- lxml - XML parsing
- ezdxf - DXF file output

## License

MIT License

## Contributing

Contributions are welcome! Please open an issue or submit a pull request.
