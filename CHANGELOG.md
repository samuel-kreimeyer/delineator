# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- DXF (AutoCAD Drawing Exchange Format) output support for watershed boundaries
  - New `--output-format dxf` option for CAD software compatibility
  - Watersheds exported as closed LWPOLYLINE entities on "WATERSHEDS" layer
  - Interior polygon holes exported on "WATERSHEDS_HOLES" layer
  - Attributes stored as XDATA (watershed_id, point_name, area, elevation stats)
  - Color-coded polygons based on watershed ID for visual distinction
  - Uses ezdxf library with DXF R2010 format

## [0.1.0] - 2026-01-18

### Added

- Initial release of Delineator
- LandXML TIN surface parsing with namespace support (1.0, 1.1, 1.2)
- COGO point extraction from LandXML files
- External study point loading (CSV, GeoJSON, Shapefile, GeoPackage)
- TIN to raster interpolation using matplotlib triangulation
- Auto-computed cell size based on average triangle edge length
- DEM conditioning with three depression handling methods:
  - `fill` - Fast depression filling
  - `breach` - Standard breach depressions (default)
  - `breach-least-cost` - Optimal least-cost breach paths
- D8 flow direction and accumulation computation
- Pour point snapping to high flow accumulation cells
- Watershed delineation with two modes:
  - Cumulative (default) - Overlapping watersheds showing total contributing area
  - Exclusive - Non-overlapping watersheds
- Watershed boundary vectorization with polygon dissolve
- Multiple output formats:
  - ESRI Shapefile
  - GeoJSON
  - OGC GeoPackage
- Statistical reports in JSON, CSV, and plain text formats
- Imperial units for all output (square feet, acres, square miles, feet)
- Coordinate system support via EPSG codes
- Comprehensive CLI with extensive options
- User documentation

### Fixed

- WhiteBox Tools path handling - now uses absolute paths to prevent silent failures
- Duplicate watershed entries - polygons are now dissolved by watershed ID
- Output file verification after WhiteBox operations

### Technical Details

- Python 3.9+ required
- Built on WhiteBox Tools for hydrological analysis
- Uses rasterio for GeoTIFF I/O
- Uses geopandas for vector operations
- Optimized statistics computation using vectorized numpy operations
