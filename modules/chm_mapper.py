"""
CHM Mapper - Automated Canopy Height Model Cartographic Map Generator

This module provides functionality to create high-quality PDF maps from
Canopy Height Model (CHM) GeoTIFF rasters with vector overlays.
"""

import numpy as np
import rasterio
from rasterio.mask import mask
from rasterio.warp import calculate_default_transform, reproject, Resampling
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib import patheffects
import warnings
import os
from pathlib import Path


def load_vector_file(vector_path):
    """
    Load vector file supporting multiple formats (.shp, .gpkg, .geojson, .kml).
    
    Parameters
    ----------
    vector_path : str
        Path to the vector file
    
    Returns
    -------
    geopandas.GeoDataFrame
        Loaded vector data
    """
    file_ext = Path(vector_path).suffix.lower()
    
    # KML files require special handling
    if file_ext == '.kml':
        # Try to enable KML driver support
        try:
            import fiona
            # Modern fiona (1.9+) may have drvsupport
            if hasattr(fiona, 'drvsupport'):
                fiona.drvsupport.supported_drivers['KML'] = 'rw'  # type: ignore
                fiona.drvsupport.supported_drivers['LIBKML'] = 'rw'  # type: ignore
        except (ImportError, AttributeError):
            pass
        
        # Read the KML file
        try:
            # Try with explicit KML driver
            gdf = gpd.read_file(vector_path, driver='KML')
        except Exception as e1:
            try:
                # Try with LIBKML driver
                gdf = gpd.read_file(vector_path, driver='LIBKML')
            except Exception as e2:
                try:
                    # Try without specifying driver (let fiona auto-detect)
                    gdf = gpd.read_file(vector_path)
                except Exception as e3:
                    raise ValueError(
                        f"Failed to read KML file with multiple methods:\n"
                        f"  - KML driver: {e1}\n"
                        f"  - LIBKML driver: {e2}\n"
                        f"  - Auto-detect: {e3}\n"
                        f"Make sure GDAL is properly installed with KML support."
                    )
        
        return gdf
    
    # For other formats (.shp, .gpkg, .geojson, etc.), use standard read_file
    else:
        return gpd.read_file(vector_path)


class CHMMapper:
    """
    A class to create cartographic maps from Canopy Height Model rasters.
    """
    
    def __init__(self, chm_path, vector_path=None):
        """
        Initialize the CHM Mapper.
        
        Parameters
        ----------
        chm_path : str
            Path to the CHM GeoTIFF raster file
        vector_path : str, optional
            Path to the vector file for boundaries
            Supports: .shp, .gpkg, .geojson, .kml
        """
        self.chm_path = chm_path
        self.vector_path = vector_path
        self.chm_data = None
        self.chm_transform = None
        self.chm_crs = None
        self.chm_nodata = None
        self.vector_gdf = None
        self.vector_gdf_full = None  # Store full vector for overview map
        self.extent = None
        
    def load_data(self):
        """Load CHM raster and vector data."""
        # Load CHM raster with proper NoData handling
        with rasterio.open(self.chm_path) as src:
            # Read with masking to handle NoData
            chm_masked = src.read(1, masked=True)
            self.chm_nodata = src.nodata
            
            # Convert masked array to regular array with NaN for NoData
            if hasattr(chm_masked, 'filled'):
                self.chm_data = np.where(chm_masked.mask, np.nan, chm_masked.data)
            else:
                self.chm_data = chm_masked
                # If nodata value exists but not masked, replace it
                if self.chm_nodata is not None:
                    self.chm_data = np.where(self.chm_data == self.chm_nodata, np.nan, self.chm_data)
            
            self.chm_transform = src.transform
            self.chm_crs = src.crs
            self.chm_bounds = src.bounds
            
        # Load vector data if provided
        if self.vector_path:
            self.vector_gdf = load_vector_file(self.vector_path)
            print(f"  Vector CRS: {self.vector_gdf.crs}")
            print(f"  CHM CRS: {self.chm_crs}")
            
            # Reproject vector to match CHM CRS if needed
            if self.vector_gdf.crs != self.chm_crs:
                print(f"  → Reprojecting vector from {self.vector_gdf.crs} to {self.chm_crs}")
                self.vector_gdf = self.vector_gdf.to_crs(self.chm_crs)
            else:
                print("  → CRS match, no reprojection needed")
            
            # Store full vector for overview map (before any filtering)
            self.vector_gdf_full = self.vector_gdf.copy()
    
    def clip_to_geometry(self, geometry_index=0):
        """
        Clip CHM raster to a specific geometry from the vector layer.
        
        Parameters
        ----------
        geometry_index : int
            Index of the geometry to clip to (default: 0)
        
        Returns
        -------
        tuple
            (clipped_data, clipped_transform, clipped_bounds)
        """
        if self.vector_gdf is None:
            raise ValueError("No vector data loaded. Cannot clip.")
        
        # Get the geometry and its bounds
        geom = self.vector_gdf.geometry.iloc[geometry_index]
        geom_bounds = geom.bounds  # (minx, miny, maxx, maxy)
        
        print(f"\n  Clipping to geometry {geometry_index}:")
        print(f"    Geometry bounds: {geom_bounds}")
        print(f"    CHM bounds: {self.chm_bounds}")
        
        # Check if geometries overlap
        chm_minx, chm_miny, chm_maxx, chm_maxy = self.chm_bounds
        geom_minx, geom_miny, geom_maxx, geom_maxy = geom_bounds
        
        if (geom_maxx < chm_minx or geom_minx > chm_maxx or 
            geom_maxy < chm_miny or geom_miny > chm_maxy):
            raise ValueError(
                f"Geometry does not overlap with CHM raster!\n"
                f"  Geometry bounds: {geom_bounds}\n"
                f"  CHM bounds: {self.chm_bounds}\n"
                f"  Are they in the same CRS? Vector CRS: {self.vector_gdf.crs}, CHM CRS: {self.chm_crs}"
            )
        
        # Prepare geometry for clipping
        geometry = [geom.__geo_interface__]
        
        # Try clipping with error handling
        try:
            # Clip the raster with proper NoData handling
            with rasterio.open(self.chm_path) as src:
                # Double-check CRS match
                if self.vector_gdf.crs != src.crs:
                    print(f"  WARNING: CRS mismatch detected during clip!")
                    print(f"    Vector: {self.vector_gdf.crs}")
                    print(f"    Raster: {src.crs}")
                    # Re-reproject to exact CRS
                    temp_geom = gpd.GeoDataFrame([geom], columns=['geometry'], crs=self.vector_gdf.crs)
                    temp_geom = temp_geom.to_crs(src.crs)
                    geometry = [temp_geom.geometry.iloc[0].__geo_interface__]
                
                clipped_data, clipped_transform = mask(src, geometry, crop=True, nodata=np.nan, filled=False)
                clipped_data = clipped_data[0]
                
                # Convert masked array to regular array with NaN
                if hasattr(clipped_data, 'filled'):
                    clipped_data = np.where(clipped_data.mask, np.nan, clipped_data.data)
                elif src.nodata is not None:
                    # Replace nodata values with NaN if not already masked
                    clipped_data = np.where(clipped_data == src.nodata, np.nan, clipped_data)
                
        except ValueError as e:
            if "do not overlap" in str(e):
                # Provide helpful error message
                raise ValueError(
                    f"Clipping failed: Geometries do not overlap!\n"
                    f"Possible causes:\n"
                    f"  1. CRS mismatch (Vector: {self.vector_gdf.crs}, CHM: {self.chm_crs})\n"
                    f"  2. Coordinates in different units/systems\n"
                    f"  3. Geometry bounds: {geom_bounds}\n"
                    f"  4. CHM bounds: {self.chm_bounds}\n"
                    f"Suggestion: Check CRS definitions in both files using QGIS or gdalinfo"
                ) from e
            else:
                raise
            
        # Calculate bounds of clipped data
        height, width = clipped_data.shape
        clipped_bounds = rasterio.transform.array_bounds(height, width, clipped_transform)
        
        print(f"    ✓ Clipped successfully: {width} x {height} pixels")
        
        return clipped_data, clipped_transform, clipped_bounds
        
        return clipped_data, clipped_transform, clipped_bounds
    
    def classify_heights(self, data, height_classes=None):
        """
        Classify tree heights into discrete classes.
        
        Parameters
        ----------
        data : numpy.ndarray
            CHM data array
        height_classes : list of tuples, optional
            List of (min, max, label) tuples for height classes
            
        Returns
        -------
        tuple
            (classified_data, classes, colors, labels)
        """
        if height_classes is None:
            # Default height classes for forestry (in meters)
            height_classes = [
                (-100, 2, '0 - 2 m'),
                (2.001, 5, '2 - 5 m'),
                (5.001, 8, '5 - 8 m'),
                (8.001, 10, '8 - 10 m'),
                (10.001, 15, '10 - 15 m'),
                (15.001, 18, '15 - 18 m'),
                (18.001, 20, '18 - 20 m'),
                (20.001, 24, '20 - 24 m'),
                (24.001, 26, '24 - 26 m'),
                (26.001, 28, '26 - 28 m'),
                (28.001, 30, '28 - 30 m'),
                (30.001, 32, '30 - 32 m'),
                (32.001, 80, '> 32 m')
            ]
        
        # Custom color scheme matching the provided legend
        colors = [
            '#E8F5E9',  # Very light green (0-2m)
            '#AED581',  # Light green (2-5m)
            '#7CB342',  # Medium green (5-8m)
            '#558B2F',  # Dark green (8-10m)
            '#2196F3',  # Blue (10-15m)
            '#FFCDD2',  # Light pink/red (15-18m)
            '#EF5350',  # Medium red (18-20m)
            '#FF6F00',  # Orange (20-24m)
            '#D3D3D3',  # Light gray (24-26m)
            '#808080',  # Medium gray (26-28m)
            '#4A4A4A',  # Dark gray/black (28-30m)
            '#1A1A1A',  # Dark gray/black (30-32m)
            '#9C27B0'   # Purple/magenta (>32m)
        ]
        
        # Create classified array
        classified_data = np.zeros_like(data, dtype=np.float32)
        classified_data[:] = np.nan
        
        for i, (min_h, max_h, label) in enumerate(height_classes):
            mask_class = (data >= min_h) & (data < max_h)
            classified_data[mask_class] = i
        
        labels = [cls[2] for cls in height_classes]
        
        return classified_data, height_classes, colors[:len(height_classes)], labels
    
    def create_map(self, data, transform, bounds, title='Baumhöhenkarte',
                   subtitle='', output_path='chm_map.pdf', 
                   height_classes=None, dpi=300, figsize=(11.69, 16.53),
                   add_overview=True, config=None, export_geotiff=True,
                   export_geopdf=True):
        """
        Create a complete cartographic map with all elements.
        
        Parameters
        ----------
        data : numpy.ndarray
            CHM data to plot (can be raw or already clipped)
        transform : affine.Affine
            Raster transform
        bounds : tuple
            (minx, miny, maxx, maxy) bounds of the data
        title : str
            Map title
        subtitle : str
            Map subtitle
        output_path : str
            Path for output PDF file
        height_classes : list, optional
            Custom height classes
        dpi : int
            Output resolution (default: 300)
        figsize : tuple
            Figure size in inches (default: A4 portrait)
        add_overview : bool
            Whether to add overview/locator map (default: True)
        config : dict, optional
            Configuration dictionary for customizing map elements
        export_geotiff : bool
            Whether to export a georeferenced GeoTIFF alongside the PDF (default: True)
        export_geopdf : bool
            Whether to export a GeoPDF with embedded geospatial metadata (default: True)
        """
        # Use default config if none provided
        if config is None:
            config = {}
        
        # Classify the data
        classified_data, classes, colors, labels = self.classify_heights(data, height_classes)
        
        # Create figure and axis
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        fig.patch.set_facecolor('white')
        
        # Calculate extent for imshow
        minx, miny, maxx, maxy = bounds
        extent = [minx, maxx, miny, maxy]
        
        # Create colormap
        cmap = ListedColormap(colors)
        norm = BoundaryNorm(np.arange(len(colors) + 1) - 0.5, len(colors))
        
        # Plot the classified CHM
        im = ax.imshow(classified_data, cmap=cmap, norm=norm, 
                       extent=extent, interpolation='nearest', origin='upper')
        
        # Overlay vector boundaries if available
        if self.vector_gdf is not None:
            vector_lw = config.get('vector_linewidth', 1.5)
            self.vector_gdf.boundary.plot(ax=ax, color='black', linewidth=vector_lw)
        
        # Remove axes
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_frame_on(True)
        border_lw = config.get('border_linewidth', 2)
        for spine in ax.spines.values():
            spine.set_edgecolor('black')
            spine.set_linewidth(border_lw)
            spine.set_linestyle('-')  # Solid line
        
        # Set aspect ratio - use 'equal' to maintain correct spatial scale
        # This may adjust the axes within the figure to preserve real-world proportions
        ax.set_aspect('equal', adjustable='box')
        
        # Add title and subtitle (override with config if provided)
        display_title = config.get('title', title)
        display_subtitle = config.get('subtitle', subtitle)
        title_x = config.get('title_position_x', 0.5)
        title_y = config.get('title_position_y', 0.98)
        title_align = config.get('title_align', 'center')
        title_fs = config.get('title_fontsize', 18)
        subtitle_fs = config.get('subtitle_fontsize', 12)
        fig.text(title_x, title_y, display_title, ha=title_align, va='top', 
                fontsize=title_fs, fontweight='bold', transform=fig.transFigure)
        if display_subtitle:
            fig.text(title_x, title_y - 0.025, display_subtitle, ha=title_align, va='top',
                    fontsize=subtitle_fs, transform=fig.transFigure, style='italic')
        
        # Add legend
        self._add_legend(fig, colors, labels, config)
        
        # Add scale bar
        self._add_scalebar(ax, config)
        
        # Add scale text (Maßstab 1:xxxx)
        self._add_scale_text(fig, ax, bounds, figsize, config)
        
        # Add north arrow (with optional overview map) - now independent element
        self._add_north_arrow(fig, with_overview=add_overview, config=config)
        
        # Adjust layout - skip if we have inset axes (overview map)
        if not add_overview:
            plt.tight_layout(rect=[0, 0.05, 1, 0.95])
        else:
            # Manual adjustment for maps with overview inset
            map_left = config.get('map_left', 0.05)
            map_right = config.get('map_right', 0.95)
            map_top = config.get('map_top', 0.93)
            map_bottom = config.get('map_bottom', 0.07)
            plt.subplots_adjust(left=map_left, right=map_right, top=map_top, bottom=map_bottom)
            
            # Optional: Show box border for debugging layout
            if config.get('show_box_border', False):
                from matplotlib.patches import Rectangle
                box_rect = Rectangle((map_left, map_bottom), 
                                    map_right - map_left, 
                                    map_top - map_bottom,
                                    transform=fig.transFigure,
                                    fill=False, 
                                    edgecolor='black', 
                                    linewidth=1,
                                    linestyle='-',  # Solid line
                                    zorder=1000)
                fig.patches.append(box_rect)
        
        # Save to PDF (used as source for GeoPDF, then deleted)
        plt.savefig(output_path, format='pdf', 
                   dpi=dpi, facecolor='white')
        
        # Export GeoPDF (same visual as PDF, but with geo metadata)
        if export_geopdf:
            geopdf_path = output_path.replace('.pdf', '_geo.pdf')
            self._export_geopdf(fig, ax, bounds, dpi, geopdf_path, source_pdf_path=output_path)
            print(f"GeoPDF saved to: {geopdf_path}")
            # Remove the intermediate PDF (GeoPDF is the final output)
            import os
            if os.path.exists(output_path):
                os.remove(output_path)
        
        # # Export normal PDF if needed (uncomment to enable)
        # print(f"Map saved to: {output_path}")
        
        # # Export GeoTIFF if requested (uncomment to enable)
        # if export_geotiff:
        #     geotiff_path = output_path.replace('.pdf', '.tif')
        #     self._export_geotiff(classified_data, colors, transform, geotiff_path)
        #     print(f"GeoTIFF saved to: {geotiff_path}")
        
        plt.close()
    
    def _export_geopdf(self, fig, ax, bounds, dpi, output_path, source_pdf_path=None):
        """
        Export a GeoPDF by adding ISO 32000 geospatial metadata to the existing PDF.
        
        This preserves the exact visual appearance of the original matplotlib PDF
        while embedding coordinate reference information. The resulting GeoPDF can
        be opened in apps like Avenza Maps for GPS navigation.
        
        Uses pikepdf to inject a Viewport with Measure/GCS dictionaries into the
        PDF page, following the ISO 32000 geospatial PDF extension.
        
        Parameters
        ----------
        fig : matplotlib.figure.Figure
            The rendered figure (used to compute axes position)
        ax : matplotlib.axes.Axes
            The main map axes (used to compute geo-extent within the figure)
        bounds : tuple
            (minx, miny, maxx, maxy) geographic bounds of the map data
        dpi : int
            Resolution (unused, kept for API compatibility)
        output_path : str
            Path for output GeoPDF file
        source_pdf_path : str, optional
            Path to the source PDF to add geo metadata to
        """
        try:
            import pikepdf
            from pyproj import Transformer, CRS
            
            if source_pdf_path is None:
                print("       ⚠ GeoPDF export failed: no source PDF path provided")
                return
            
            # Open the already-saved matplotlib PDF
            pdf = pikepdf.open(source_pdf_path)
            page = pdf.pages[0]
            
            # Get page dimensions in PDF points (72 points/inch)
            mediabox = page.MediaBox
            page_width = float(mediabox[2]) - float(mediabox[0])
            page_height = float(mediabox[3]) - float(mediabox[1])
            
            # Get axes position in figure-fraction coordinates
            ax_bbox = ax.get_position()
            
            # Convert to PDF points (PDF origin is bottom-left, matching matplotlib)
            bbox_x1 = ax_bbox.x0 * page_width
            bbox_y1 = ax_bbox.y0 * page_height
            bbox_x2 = ax_bbox.x1 * page_width
            bbox_y2 = ax_bbox.y1 * page_height
            
            # Get geographic extent from axes (respects aspect='equal' adjustments)
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            geo_minx, geo_maxx = xlim
            geo_miny, geo_maxy = ylim
            
            # Determine source CRS
            if self.chm_crs is not None:
                src_crs = CRS(self.chm_crs)
            else:
                src_crs = CRS.from_epsg(4326)
            
            # Transform 4 corners to WGS84 (required by ISO 32000 GPTS)
            # Corner order: BL, TL, TR, BR (matching LPTS)
            corners_x = [geo_minx, geo_minx, geo_maxx, geo_maxx]
            corners_y = [geo_miny, geo_maxy, geo_maxy, geo_miny]
            
            wgs84 = CRS.from_epsg(4326)
            if src_crs != wgs84:
                transformer = Transformer.from_crs(src_crs, wgs84, always_xy=True)
                lons, lats = transformer.transform(corners_x, corners_y)
            else:
                lons, lats = corners_x, corners_y
            
            # Build GPTS: pairs of (latitude, longitude) for BL, TL, TR, BR
            gpts = []
            for lat, lon in zip(lats, lons):
                gpts.extend([float(lat), float(lon)])
            
            # LPTS: normalized viewport coordinates for BL, TL, TR, BR
            lpts = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0]
            
            # Bounds: clipping polygon in normalized viewport coordinates
            bounds_arr = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0]
            
            # Build GCS (Geographic/Projected Coordinate System) dictionary
            wkt_str = src_crs.to_wkt()
            gcs_type = "/PROJCS" if src_crs.is_projected else "/GEOGCS"
            
            gcs_dict = pikepdf.Dictionary({
                "/Type": pikepdf.Name(gcs_type),
                "/WKT": pikepdf.String(wkt_str)
            })
            
            # Build Measure dictionary (ISO 32000 geospatial extension)
            measure_dict = pikepdf.Dictionary({
                "/Type": pikepdf.Name("/Measure"),
                "/Subtype": pikepdf.Name("/GEO"),
                "/Bounds": pikepdf.Array(bounds_arr),
                "/GPTS": pikepdf.Array(gpts),
                "/LPTS": pikepdf.Array(lpts),
                "/GCS": gcs_dict
            })
            
            # Build Viewport dictionary
            viewport = pikepdf.Dictionary({
                "/Type": pikepdf.Name("/Viewport"),
                "/BBox": pikepdf.Array([
                    float(bbox_x1), float(bbox_y1),
                    float(bbox_x2), float(bbox_y2)
                ]),
                "/Measure": measure_dict
            })
            
            # Add VP (Viewport) array to the PDF page
            page[pikepdf.Name("/VP")] = pikepdf.Array([viewport])
            
            # Save the georeferenced PDF
            pdf.save(output_path)
            pdf.close()
            
        except ImportError:
            print("       ⚠ GeoPDF export requires pikepdf and pyproj.")
            print("         Install with: pip install pikepdf")
        except Exception as e:
            print(f"       ⚠ GeoPDF export failed: {e}")

    def _export_geotiff(self, classified_data, colors, transform, output_path):
        """
        Export classified CHM data as a georeferenced RGB GeoTIFF.
        
        Parameters
        ----------
        classified_data : numpy.ndarray
            Classified height data (indices into color array)
        colors : list
            List of color hex codes for each class
        transform : affine.Affine
            Raster geotransform
        output_path : str
            Path for output GeoTIFF file
        """
        from matplotlib.colors import to_rgb
        
        # Convert classified data to RGB
        height, width = classified_data.shape
        rgb_array = np.zeros((height, width, 3), dtype=np.uint8)
        
        # Fill RGB array based on classification
        for class_idx, color_hex in enumerate(colors):
            mask = (classified_data == class_idx)
            if np.any(mask):
                rgb = to_rgb(color_hex)
                rgb_array[mask] = [int(c * 255) for c in rgb]
        
        # Handle NoData (NaN values) - make them transparent/white
        nodata_mask = np.isnan(classified_data)
        rgb_array[nodata_mask] = [255, 255, 255]  # White for NoData
        
        # Write to GeoTIFF with georeferencing
        with rasterio.open(
            output_path,
            'w',
            driver='GTiff',
            height=height,
            width=width,
            count=3,  # RGB
            dtype=rasterio.uint8,
            crs=self.chm_crs,
            transform=transform,
            compress='lzw',  # Compress to save space
            photometric='RGB'
        ) as dst:
            # Write RGB bands
            dst.write(rgb_array[:, :, 0], 1)  # Red
            dst.write(rgb_array[:, :, 1], 2)  # Green
            dst.write(rgb_array[:, :, 2], 3)  # Blue
            
            # Add color interpretation
            dst.colorinterp = [rasterio.enums.ColorInterp.red,
                              rasterio.enums.ColorInterp.green,
                              rasterio.enums.ColorInterp.blue]
    
    def _add_legend(self, fig, colors, labels, config=None):
        """Add a discrete legend to the figure."""
        if config is None:
            config = {}
        
        # Get config values
        ncol = config.get('legend_ncol', 3)
        fontsize = config.get('legend_fontsize', 9)
        title_fontsize = config.get('legend_title_fontsize', 10)
        legend_x = config.get('legend_position_x', 0.5)
        legend_y = config.get('legend_position_y', 0.0)
        legend_loc = config.get('legend_loc', 'lower center')
        
        # Size control parameters
        labelspacing = config.get('legend_labelspacing', 0.5)  # Vertical space between entries
        handlelength = config.get('legend_handlelength', 2.0)  # Width of color boxes
        handletextpad = config.get('legend_handletextpad', 0.8)  # Space between box and text
        columnspacing = config.get('legend_columnspacing', 2.0)  # Space between columns
        legend_title = config.get('legend_title', 'Baumhöhenklassen')
        
        # Additional text lines below legend title
        legend_subtitle1 = config.get('legend_subtitle1', '')
        legend_subtitle2 = config.get('legend_subtitle2', '')
        legend_subtitle_fontsize = config.get('legend_subtitle_fontsize', 8)
        legend_border_linewidth = config.get('legend_border_linewidth', 1)  # Legend frame border thickness
        
        # Custom background box for legend (independent of matplotlib legend)
        legend_background_box_x = config.get('legend_background_box_x', None)  # Left edge in figure coords
        legend_background_box_y = config.get('legend_background_box_y', None)  # Bottom edge in figure coords
        legend_background_box_width = config.get('legend_background_box_width', None)  # Width in figure coords
        legend_background_box_height = config.get('legend_background_box_height', None)  # Height in figure coords
        legend_background_box_linewidth = config.get('legend_background_box_linewidth', 1)  # Border thickness
        
        # Create legend patches - add subtitle entries first (inside box)
        patches = []
        patch_labels = []
        
        # Draw custom background box BEFORE legend (if configured)
        if all([legend_background_box_x is not None, 
                legend_background_box_y is not None,
                legend_background_box_width is not None,
                legend_background_box_height is not None]):
            from matplotlib.patches import Rectangle
            background_box = Rectangle(
                (legend_background_box_x, legend_background_box_y),
                legend_background_box_width,
                legend_background_box_height,
                transform=fig.transFigure,
                facecolor='none',  # Transparent fill
                edgecolor='black',  # Black border
                linewidth=legend_background_box_linewidth,
                zorder=1  # Behind legend
            )
            fig.patches.append(background_box)
        
        # Add subtitles as text-only entries (no color patch, inside the legend box)
        if legend_subtitle1:
            patches.append(mpatches.Patch(color='none', linewidth=0))
            patch_labels.append(legend_subtitle1)
        if legend_subtitle2:
            patches.append(mpatches.Patch(color='none', linewidth=0))
            patch_labels.append(legend_subtitle2)
        
        # Add the actual color patches
        for i in range(len(colors)):
            patches.append(mpatches.Patch(color=colors[i], label=labels[i]))
            patch_labels.append(labels[i])
        
        # Add legend to figure with only main title (outside box)
        legend = fig.legend(handles=patches, labels=patch_labels, loc=legend_loc, 
                          ncol=min(ncol, len(patches)), 
                          title=legend_title,  # Only main title, outside box
                          frameon=True, fancybox=True,
                          bbox_to_anchor=(legend_x, legend_y),
                          fontsize=fontsize, title_fontsize=title_fontsize,
                          labelspacing=labelspacing,
                          handlelength=handlelength,
                          handletextpad=handletextpad,
                          columnspacing=columnspacing,
                          alignment='left')  # Align legend entries to the left
        legend.get_frame().set_facecolor('white')
        legend.get_frame().set_alpha(0.95)
        legend.get_frame().set_edgecolor('black')
        legend.get_frame().set_linewidth(legend_border_linewidth)
        
        # Make subtitle entries use smaller font and align left
        if legend_subtitle1 or legend_subtitle2:
            num_subtitles = (1 if legend_subtitle1 else 0) + (1 if legend_subtitle2 else 0)
            for i, text in enumerate(legend.get_texts()):
                if i < num_subtitles:
                    # This is a subtitle line - use smaller font, italic, and left-aligned
                    text.set_fontsize(legend_subtitle_fontsize)
                    text.set_style('italic')
                    text.set_ha('left')  # Set horizontal alignment
                    # Adjust x position to align left within legend box
                    text.set_x(-handlelength - handletextpad)
    
    def _add_scalebar(self, ax, config=None):
        """Add a metric scale bar to the map."""
        if config is None:
            config = {}
        
        # Get config values
        fontsize = config.get('scalebar_fontsize', 9)
        height_frac = config.get('scalebar_height', 0.015)
        length_frac = config.get('scalebar_length', 0.2)
        location = config.get('scalebar_location', 'lower right')
        pad_x = config.get('scalebar_pad_x', 0.5)
        pad_y = config.get('scalebar_pad_y', 0.5)
        
        # Create scale bar (matplotlib-scalebar automatically handles CRS)
        scalebar = ScaleBar(
            dx=1,  # 1 meter per pixel (will be adjusted automatically)
            units='m',
            location=location,
            length_fraction=length_frac,
            height_fraction=height_frac,
            pad=pad_x,
            sep=pad_y,
            box_alpha=0.8,
            box_color='white',
            color='black',
            font_properties={'size': fontsize}
        )
        ax.add_artist(scalebar)
    
    def _add_scale_text(self, fig, ax, bounds, figsize, config=None):
        """Add map scale text (Maßstab 1:xxxx) to the figure."""
        if config is None:
            config = {}
        
        # Get config values
        scale_text_x = config.get('scale_text_position_x', 0.88)
        scale_text_y = config.get('scale_text_position_y', 0.08)
        scale_text_fontsize = config.get('scale_text_fontsize', 10)
        scale_text_align = config.get('scale_text_align', 'left')
        
        # Calculate map scale
        # Get the actual displayed width of the axes in figure units
        bbox = ax.get_position()
        ax_width_fig = bbox.width  # Width in figure coordinates (0-1)
        ax_height_fig = bbox.height  # Height in figure coordinates (0-1)
        
        # Convert to inches
        fig_width_inches, fig_height_inches = figsize
        ax_width_inches = ax_width_fig * fig_width_inches
        ax_height_inches = ax_height_fig * fig_height_inches
        
        # Get map extent in real-world units (meters)
        minx, miny, maxx, maxy = bounds
        map_width_meters = maxx - minx
        map_height_meters = maxy - miny
        
        # Calculate scale based on width (meters per inch)
        meters_per_inch = map_width_meters / ax_width_inches
        
        # Convert to scale ratio (1:xxxx)
        # 1 inch = 0.0254 meters, so map_meters / 0.0254 gives map scale
        scale_ratio = meters_per_inch / 0.0254
        
        # Round to nice numbers
        if scale_ratio < 1000:
            scale_ratio = round(scale_ratio, -1)  # Round to nearest 10
        elif scale_ratio < 10000:
            scale_ratio = round(scale_ratio, -2)  # Round to nearest 100
        else:
            scale_ratio = round(scale_ratio, -3)  # Round to nearest 1000
        
        # Format the scale text
        scale_text = f"Maßstab 1:{int(scale_ratio)}"
        
        # Add text to figure
        fig.text(scale_text_x, scale_text_y, scale_text,
                ha=scale_text_align, va='bottom',
                fontsize=scale_text_fontsize,
                transform=fig.transFigure,
                fontweight='bold')
    
    def _add_north_arrow(self, fig, with_overview=False, config=None):
        """Add a north arrow to the figure as an independent element."""
        if config is None:
            config = {}
        
        # Get config values - now using figure coordinates (0-1)
        x_pos = config.get('north_arrow_position_x', 0.88)
        y_pos = config.get('north_arrow_position_y', 0.85)
        arrow_length = config.get('north_arrow_length', 0.04)  # As fraction of figure height
        arrow_width = config.get('north_arrow_width', 2)
        fontsize = config.get('north_arrow_fontsize', 14)
        pad = config.get('north_arrow_pad', 0.3)
        arrow_style = config.get('north_arrow_style', 'fancy')  # Arrow style: '->', 'fancy', 'wedge', 'simple'
        symbol = config.get('north_arrow_symbol', 'N')  # Symbol to display: 'N', '↑', '▲', '⬆', '★'
        
        # Calculate arrow positions: arrow points UP, with symbol below
        arrow_end_y = y_pos + arrow_length  # Arrow tip at top
        
        # Draw arrow using figure coordinates (from y_pos upward to arrow_end_y)
        arrow = mpatches.FancyArrowPatch(
            (x_pos, y_pos), (x_pos, arrow_end_y),
            arrowstyle=arrow_style, mutation_scale=20,
            linewidth=arrow_width, color='black',
            transform=fig.transFigure,
            zorder=1000
        )
        fig.patches.append(arrow)
        
        # Add symbol text with white box below the arrow
        fig.text(x_pos, y_pos - 0.01, symbol,
                ha='center', va='top',
                fontsize=fontsize, fontweight='bold',
                transform=fig.transFigure,
                bbox=dict(boxstyle=f'round,pad={pad}', 
                         facecolor='white', 
                         edgecolor='black', 
                         alpha=0.8),
                zorder=1001)
        
        # Add overview map if requested
        if with_overview and self.vector_gdf is not None:
            self._add_overview_map(fig, config)
    
    def _add_overview_map(self, fig, config=None):
        """Add a small overview/locator map showing CHM extent within full vector extent."""
        if config is None:
            config = {}
        
        # Get config values
        inset_width = config.get('overview_width', 0.12)
        inset_height = config.get('overview_height', 0.12)
        fontsize = config.get('overview_fontsize', 7)
        border_width = config.get('overview_border_width', 1)
        chm_box_width = config.get('overview_chm_box_width', 1.5)
        overview_label = config.get('overview_label', 'Übersicht')
        location_name = config.get('location_name', '')
        location_fontsize = config.get('location_fontsize', 10)
        location_y_offset = config.get('location_y_offset', 0.02)
        
        # Get position in figure coordinates (0-1)
        overview_x = config.get('overview_position_x', 0.88)  # X position in figure coords
        overview_y = config.get('overview_position_y', 0.35)  # Y position in figure coords (bottom left corner)
        
        # Create inset axes at specified position
        # [left, bottom, width, height] all in figure coordinates (0-1)
        inset_ax = fig.add_axes([overview_x - inset_width/2, overview_y, 
                                 inset_width, inset_height])
        
        # Use full vector extent (before filtering) if available, otherwise use current
        vector_for_overview = self.vector_gdf_full if self.vector_gdf_full is not None else self.vector_gdf
        vector_bounds = vector_for_overview.total_bounds  # (minx, miny, maxx, maxy)
        
        # Add buffer to extent so vectors don't touch the border
        overview_buffer = config.get('overview_buffer', 0.05)  # 5% buffer by default
        x_range = vector_bounds[2] - vector_bounds[0]
        y_range = vector_bounds[3] - vector_bounds[1]
        buffer_x = x_range * overview_buffer
        buffer_y = y_range * overview_buffer
        
        # Plot all vector polygons in light gray
        vector_for_overview.plot(ax=inset_ax, facecolor='lightgray', edgecolor='gray', 
                            linewidth=0.5, alpha=0.7)
        
        # Highlight CHM extent with red rectangle
        chm_minx, chm_miny, chm_maxx, chm_maxy = self.chm_bounds
        chm_rect = Rectangle((chm_minx, chm_miny), 
                             chm_maxx - chm_minx, 
                             chm_maxy - chm_miny,
                             linewidth=chm_box_width, edgecolor='red', facecolor='none')
        inset_ax.add_patch(chm_rect)
        
        # Set extent to full vector bounds with buffer
        inset_ax.set_xlim(vector_bounds[0] - buffer_x, vector_bounds[2] + buffer_x)
        inset_ax.set_ylim(vector_bounds[1] - buffer_y, vector_bounds[3] + buffer_y)
        
        # Remove ticks and labels
        inset_ax.set_xticks([])
        inset_ax.set_yticks([])
        
        # Add frame
        for spine in inset_ax.spines.values():
            spine.set_edgecolor('black')
            spine.set_linewidth(border_width)
        
        inset_ax.set_aspect('equal')
        
        # Add overview label below the map
        inset_ax.text(0.5, -0.15, overview_label, transform=inset_ax.transAxes,
                     ha='center', va='top', fontsize=fontsize, style='italic')
        
        # Add location name above the overview map if provided
        if location_name:
            # Position text above the overview box
            text_y = overview_y + inset_height + location_y_offset
            fig.text(overview_x, text_y, location_name,
                    ha='center', va='bottom',
                    fontsize=location_fontsize,
                    fontweight='bold',
                    transform=fig.transFigure)
    
    def process_all_areas(self, output_dir='output', title_field=None):
        """
        Process all geometries in the vector layer and create individual maps.
        
        Parameters
        ----------
        output_dir : str
            Directory for output PDF files
        title_field : str, optional
            Field name in vector attributes to use for map titles
        """
        import os
        
        if self.vector_gdf is None:
            raise ValueError("No vector data loaded.")
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Process each geometry
        for idx in range(len(self.vector_gdf)):
            print(f"\nProcessing area {idx + 1}/{len(self.vector_gdf)}...")
            
            # Get title from attribute if specified
            if title_field and title_field in self.vector_gdf.columns:
                area_name = str(self.vector_gdf.iloc[idx][title_field])
            else:
                area_name = f"Area_{idx + 1}"
            
            # Clip data
            clipped_data, clipped_transform, clipped_bounds = self.clip_to_geometry(idx)
            
            # Create output filename
            safe_name = "".join(c if c.isalnum() or c in (' ', '_', '-') else '_' 
                               for c in area_name)
            output_path = os.path.join(output_dir, f"chm_map_{safe_name}.pdf")
            
            # Create map
            self.create_map(
                data=clipped_data,
                transform=clipped_transform,
                bounds=clipped_bounds,
                title='Baumhöhenkarte',
                subtitle=area_name,
                output_path=output_path
            )
        
        print(f"\n✓ All maps generated successfully in '{output_dir}' directory")


def create_chm_map(chm_path, vector_path=None, output_path='chm_map.pdf',
                   title='Baumhöhenkarte', subtitle='', 
                   height_classes=None, clip_to_area=True, geometry_index=0):
    """
    Convenience function to create a CHM map in one call.
    
    Parameters
    ----------
    chm_path : str
        Path to CHM GeoTIFF file
    vector_path : str, optional
        Path to vector file (supports .shp, .gpkg, .geojson, .kml)
    output_path : str
        Output PDF path
    title : str
        Map title
    subtitle : str
        Map subtitle
    height_classes : list, optional
        Custom height classification
    clip_to_area : bool
        Whether to clip to a specific geometry in vector (if vector provided)
    geometry_index : int
        Index of geometry to clip to (default: 0). Only used if clip_to_area=True.
    """
    mapper = CHMMapper(chm_path, vector_path)
    mapper.load_data()
    
    if vector_path and clip_to_area:
        # Clip to specified geometry
        try:
            clipped_data, clipped_transform, clipped_bounds = mapper.clip_to_geometry(geometry_index)
            mapper.create_map(clipped_data, clipped_transform, clipped_bounds,
                             title, subtitle, output_path, height_classes)
        except ValueError as e:
            if "do not overlap" in str(e):
                print(f"\n⚠ Geometry {geometry_index} does not overlap with CHM!")
                print(f"  Use Cell 11 to find overlapping geometries, or set clip_to_area=False")
                raise
            else:
                raise
    else:
        # Use full extent
        mapper.create_map(mapper.chm_data, mapper.chm_transform, 
                         mapper.chm_bounds, title, subtitle, 
                         output_path, height_classes)
