"""
CHM Mapper Modules

This package provides functionality for creating high-quality cartographic maps
from Canopy Height Model (CHM) GeoTIFF rasters with vector overlays.
"""

from .chm_mapper import CHMMapper, create_chm_map, load_vector_file
from .config import *
from .run_mapper import main as run_mapper_main

__version__ = "1.0.0"
__all__ = ['CHMMapper', 'create_chm_map', 'load_vector_file', 'run_mapper_main']
