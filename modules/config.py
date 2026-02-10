"""
Configuration file for CHM Mapper.

This file contains all customizable parameters for your mapping workflow.
Copy this file and adjust the values to match your project requirements.
"""

# ============================================================================
# DATA PATHS
# ============================================================================

# Path to your CHM GeoTIFF raster file
CHM_RASTER_PATH = "data/chm_raster.tif"

# Path to your vector file
# Supported formats: .shp (Shapefile), .gpkg (GeoPackage), .geojson, .kml
# GeoPackage (.gpkg) is recommended: single file, faster, no size limits
VECTOR_SHAPEFILE_PATH = "data/forest_compartments.gpkg"  # or .shp, .kml, .geojson

# Output directory for generated maps
OUTPUT_DIRECTORY = "output"

# ============================================================================
# MAP SETTINGS
# ============================================================================

# Main title for the map
MAP_TITLE = "Baumhöhenkarte"

# Subtitle (can be set per area or use a field from the shapefile)
MAP_SUBTITLE = "Forstgebiet"

# Field name in shapefile to use for area-specific subtitles
# Set to None to use generic names
SUBTITLE_FIELD = "Name"  # e.g., "Name", "Compartment", "ID", etc.

# ============================================================================
# HEIGHT CLASSIFICATION
# ============================================================================

# Define height classes as list of tuples: (min_height, max_height, label)
# Heights are in meters
HEIGHT_CLASSES = [
    (0, 2, '0-2 m'),
    (2, 5, '2-5 m'),
    (5, 10, '5-10 m'),
    (10, 15, '10-15 m'),
    (15, 20, '15-20 m'),
    (20, 25, '20-25 m'),
    (25, 30, '25-30 m'),
    (30, 35, '30-35 m'),
    (35, 1000, '>35 m')
]

# Alternative: Forestry development stages
# HEIGHT_CLASSES = [
#     (0, 5, '0-5 m (Jungwuchs)'),
#     (5, 15, '5-15 m (Dickung)'),
#     (15, 25, '15-25 m (Stangenholz)'),
#     (25, 35, '25-35 m (Baumholz)'),
#     (35, 1000, '>35 m (Altholz)')
# ]

# ============================================================================
# COLOR SCHEME
# ============================================================================

# Custom colors for height classes (one per class)
# Colors in hex format. Must match the number of HEIGHT_CLASSES
# Default: gradient from light yellow-green to dark green
CUSTOM_COLORS = [
    '#FFF7BC',  # Very light yellow
    '#FEE391',  # Light yellow
    '#FEC44F',  # Yellow-orange
    '#D4EE9F',  # Light green
    '#A6D96A',  # Medium-light green
    '#66BD63',  # Medium green
    '#1A9850',  # Dark green
    '#006837',  # Very dark green
    '#00441B'   # Darkest green
]

# Alternative color schemes:
# 
# Viridis-like (purple to yellow-green):
# CUSTOM_COLORS = ['#440154', '#31688e', '#35b779', '#90d743', '#fde724']
#
# Brown-Green (soil to trees):
# CUSTOM_COLORS = ['#8C510A', '#D8B365', '#F6E8C3', '#C7EAE5', '#5AB4AC', '#01665E']

# ============================================================================
# OUTPUT SETTINGS
# ============================================================================

# Output resolution (DPI) for PDF
# 300 DPI is high quality for printing
# 150 DPI is good for digital viewing
OUTPUT_DPI = 300

# Figure size in inches (width, height)
# A4 portrait: (8.27, 11.69)
# A4 landscape: (11.69, 8.27)
# A3 portrait: (11.69, 16.53)
FIGURE_SIZE = (11.69, 16.53)  # A3 portrait

# ============================================================================
# VECTOR OVERLAY SETTINGS
# ============================================================================

# Boundary line color
BOUNDARY_COLOR = 'black'

# Boundary line width
BOUNDARY_LINEWIDTH = 1.5

# ============================================================================
# LEGEND SETTINGS
# ============================================================================

# Legend title
LEGEND_TITLE = 'Baumhöhenklassen'

# Number of columns in legend
LEGEND_NCOLS = 3

# Legend font size
LEGEND_FONTSIZE = 9

# ============================================================================
# PROCESSING OPTIONS
# ============================================================================

# Whether to clip CHM to each vector geometry
CLIP_TO_GEOMETRY = True

# Whether to process all areas in batch
BATCH_PROCESS = False

# If False and BATCH_PROCESS is False, which geometry index to process (0-based)
SINGLE_GEOMETRY_INDEX = 0

# ============================================================================
# ADVANCED SETTINGS
# ============================================================================

# NoData value for the CHM raster (will be set to NaN)
# Usually handled automatically by rasterio
NODATA_VALUE = None

# Whether to reproject vector to match CHM CRS automatically
AUTO_REPROJECT = True

# Map frame settings
FRAME_LINEWIDTH = 2
FRAME_COLOR = 'black'
