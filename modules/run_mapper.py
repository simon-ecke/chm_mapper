"""
Main script to run the CHM mapping workflow using configuration file.

This script loads settings from config.py and generates maps accordingly.
"""

import os
import sys
from .chm_mapper import CHMMapper, create_chm_map
from . import config


def main():
    """Main execution function."""
    
    print("\n" + "=" * 70)
    print("CHM MAPPER - Automated Canopy Height Model Mapping")
    print("=" * 70 + "\n")
    
    # Validate input files
    if not os.path.exists(config.CHM_RASTER_PATH):
        print(f"ERROR: CHM raster file not found: {config.CHM_RASTER_PATH}")
        print("Please update the path in config.py")
        sys.exit(1)
    
    if config.VECTOR_SHAPEFILE_PATH and not os.path.exists(config.VECTOR_SHAPEFILE_PATH):
        print(f"WARNING: Vector shapefile not found: {config.VECTOR_SHAPEFILE_PATH}")
        print("Continuing without vector overlay...\n")
        vector_path = None
    else:
        vector_path = config.VECTOR_SHAPEFILE_PATH
    
    # Create output directory
    os.makedirs(config.OUTPUT_DIRECTORY, exist_ok=True)
    print(f"Output directory: {config.OUTPUT_DIRECTORY}\n")
    
    # Initialize mapper
    print("Loading data...")
    mapper = CHMMapper(config.CHM_RASTER_PATH, vector_path)
    mapper.load_data()
    print("✓ Data loaded successfully\n")
    
    # Determine processing mode
    if config.BATCH_PROCESS and vector_path:
        print(f"BATCH MODE: Processing all {len(mapper.vector_gdf)} areas...\n")
        mapper.process_all_areas(
            output_dir=config.OUTPUT_DIRECTORY,
            title_field=config.SUBTITLE_FIELD
        )
        
    elif vector_path and config.CLIP_TO_GEOMETRY:
        print(f"SINGLE AREA MODE: Processing geometry {config.SINGLE_GEOMETRY_INDEX}...\n")
        
        # Get subtitle from field if available
        if config.SUBTITLE_FIELD and config.SUBTITLE_FIELD in mapper.vector_gdf.columns:
            subtitle = str(mapper.vector_gdf.iloc[config.SINGLE_GEOMETRY_INDEX][config.SUBTITLE_FIELD])
        else:
            subtitle = config.MAP_SUBTITLE
        
        # Clip to geometry
        clipped_data, clipped_transform, clipped_bounds = mapper.clip_to_geometry(
            config.SINGLE_GEOMETRY_INDEX
        )
        
        # Generate output filename
        safe_subtitle = "".join(c if c.isalnum() or c in (' ', '_', '-') else '_' 
                               for c in subtitle)
        output_path = os.path.join(config.OUTPUT_DIRECTORY, 
                                   f"chm_map_{safe_subtitle}.pdf")
        
        # Create map
        mapper.create_map(
            data=clipped_data,
            transform=clipped_transform,
            bounds=clipped_bounds,
            title=config.MAP_TITLE,
            subtitle=subtitle,
            output_path=output_path,
            height_classes=config.HEIGHT_CLASSES,
            dpi=config.OUTPUT_DPI,
            figsize=config.FIGURE_SIZE
        )
        
    else:
        print("FULL EXTENT MODE: Processing entire CHM raster...\n")
        
        output_path = os.path.join(config.OUTPUT_DIRECTORY, "chm_map_full_extent.pdf")
        
        # Create map with full extent
        mapper.create_map(
            data=mapper.chm_data,
            transform=mapper.chm_transform,
            bounds=mapper.chm_bounds,
            title=config.MAP_TITLE,
            subtitle=config.MAP_SUBTITLE,
            output_path=output_path,
            height_classes=config.HEIGHT_CLASSES,
            dpi=config.OUTPUT_DPI,
            figsize=config.FIGURE_SIZE
        )
    
    print("\n" + "=" * 70)
    print("✓ COMPLETED SUCCESSFULLY")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nOperation cancelled by user.")
        sys.exit(1)
    except Exception as e:
        print(f"\n\nERROR: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
