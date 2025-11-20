import geopandas as gpd
import simplekml
import os
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from shapely.geometry import box
from rasterio import features
from rasterio.transform import from_bounds

# --- CONFIGURATION ---
DEFAULT_COLUMN = "EM" 
DEFAULT_PIXEL_SIZE = 10.0 
DEFAULT_MIN_CUTOFF = 64.0
DEFAULT_MAX_CUTOFF = 326.0

# --- YOUR SPECIFIC COLOR RAMP ---
# We map the colors to the normalized values (0.0 to 1.0) based on your Min/Max
def create_custom_cmap(min_val, max_val):
    # Your specific breaks
    breaks = [180, 197, 206, 228, 247, 260, 272]
    
    # Normalize these breaks to 0.0 - 1.0 range
    span = max_val - min_val
    norm_breaks = [(b - min_val) / span for b in breaks]
    
    # Ensure they stay within 0-1
    norm_breaks = [max(0, min(1, b)) for b in norm_breaks]

    # Define the transition points
    # (Value, Color)
    cdict = [
        (0.0, "red"),               # Bottom to 180
        (norm_breaks[0], "orangered"), # 180
        (norm_breaks[1], "darkorange"),# 197
        (norm_breaks[2], "orange"),    # 206
        (norm_breaks[3], "yellow"),    # 228
        (norm_breaks[4], "yellowgreen"),# 247
        (norm_breaks[5], "lime"),      # 260
        (norm_breaks[6], "green"),     # 272
        (1.0, "darkgreen")             # Top
    ]
    return mcolors.LinearSegmentedColormap.from_list("em_ramp", cdict)

def get_input(prompt, default_val):
    user_val = input(f"{prompt} [Default: {default_val}]: ").strip()
    if user_val == "":
        return default_val
    try:
        return type(default_val)(user_val)
    except ValueError:
        return default_val

def get_epsg(datum, zone):
    epsg_map = {
        ('GDA94', '54'): 28354, ('GDA94', '55'): 28355,
        ('GDA2020', '54'): 7854, ('GDA2020', '55'): 7855
    }
    return epsg_map.get((datum, zone))

def main():
    print("\n--- ArcView 3.1 Advanced Interpretation Generator ---\n")

    # 1. Inputs
    datum = input("Enter Datum (GDA94 / GDA2020) [Default: GDA94]: ").strip().upper() or "GDA94"
    zone = input("Enter Zone (54 / 55) [Default: 55]: ").strip() or "55"
    epsg = get_epsg(datum, zone)
    if not epsg: return

    # Data File
    shp_path = input("Full path to DATA Shapefile: ").strip().replace('"', '')
    base_path, ext = os.path.splitext(shp_path)
    if ext.lower() == '.dbf': shp_path = base_path + ".shp"
    
    if not os.path.exists(shp_path):
        print(f"File not found: {shp_path}")
        return

    # Boundary File (New Request)
    print("\n(Required for clipping)")
    bound_path = input("Full path to BOUNDARY Shapefile: ").strip().replace('"', '')

    # Settings
    col_name = input(f"Column Name [Default: {DEFAULT_COLUMN}]: ").strip() or DEFAULT_COLUMN
    resolution = get_input("Output Resolution (meters/pixel)", DEFAULT_PIXEL_SIZE)
    min_cut = get_input("Bottom Cutoff", DEFAULT_MIN_CUTOFF)
    max_cut = get_input("Top Cutoff", DEFAULT_MAX_CUTOFF)

    # 2. Load Data
    print("Reading Data...")
    try:
        gdf = gpd.read_file(shp_path)
        gdf.set_crs(epsg=epsg, allow_override=True, inplace=True)
        
        if col_name not in gdf.columns:
            print(f"Column '{col_name}' not found. Available: {list(gdf.columns)}")
            return

        gdf = gdf[(gdf[col_name] >= min_cut) & (gdf[col_name] <= max_cut)]
        print(f"Points to process: {len(gdf)}")

    except Exception as e:
        print(f"Error reading data: {e}")
        return
    
    # 3. Interpolate (Rasterize)
    print("Creating Interpolated Image...")
    
    # Get coordinates
    x = gdf.geometry.x.values
    y = gdf.geometry.y.values
    z = gdf[col_name].values

    # Determine bounds
    x_min, x_max = x.min(), x.max()
    y_min, y_max = y.min(), y.max()
    
    # Create grid dimensions
    width = int(np.ceil((x_max - x_min) / resolution))
    height = int(np.ceil((y_max - y_min) / resolution))
    
    # Create the grid coordinates
    # Note: We use 'complex' notation for mgrid to specify exact step counts if needed, 
    # but here simple steps work. We ensure alignment for rasterio.
    grid_x, grid_y = np.meshgrid(
        np.linspace(x_min, x_max, width),
        np.linspace(y_min, y_max, height)
    )
    
    # Interpolate data onto grid
    # method='linear' is smooth. method='nearest' is blocky (ArcView style).
    grid_z = griddata((x, y), z, (grid_x, grid_y), method='linear')
    
    # griddata returns (Height, Width). Flip standard image origin (Top-Left vs Bottom-Left)
    # We will align this carefully with imshow below.

    # 4. Clipping (Masking)
    if bound_path and os.path.exists(bound_path):
        print("Applying Boundary Clip...")
        try:
            bound_gdf = gpd.read_file(bound_path)
            if bound_gdf.crs is None: bound_gdf.set_crs(epsg=epsg, inplace=True)
            else: bound_gdf = bound_gdf.to_crs(epsg=epsg)
            
            # Create a Transform (Mapping pixels to coords) for Rasterio
            # from_bounds(west, south, east, north, width, height)
            transform = from_bounds(x_min, y_min, x_max, y_max, width, height)
            
            # Generate boolean mask from shapes
            shapes = [geom for geom in bound_gdf.geometry]
            # rasterio mask: True=Inside, False=Outside (or vice versa depending on invert)
            # We want pixels OUTSIDE the polygon to be True so we can mask them?
            # geometry_mask returns True where shapes overlap if invert=False?
            # Actually: invert=True creates a mask where pixels INSIDE shapes are False (unmasked).
            mask = features.geometry_mask(shapes, transform=transform, invert=True, out_shape=(height, width))
            
            # Apply mask (Set outside pixels to NaN)
            # Note: Rasterio image origin is usually Top-Left. Griddata with linspace/meshgrid is often Bottom-Left.
            # We need to flip the mask to match the grid_z orientation from meshgrid/griddata.
            mask = np.flipud(mask) 
            
            grid_z[~mask] = np.nan
            
        except Exception as e:
            print(f"Clipping warning: {e}")

    # 5. Generate PNG
    print("Generating PNG...")
    img_output = base_path + "_Overlay.png"
    
    # Create the correct Color Ramp
    cmap = create_custom_cmap(min_cut, max_cut)
    
    plt.figure(figsize=(10, 10), frameon=False)
    
    # origin='lower' aligns with x_min->x_max, y_min->y_max logic
    plt.imshow(grid_z, cmap=cmap, vmin=min_cut, vmax=max_cut, origin='lower')
    plt.axis('off')
    
    # Save transparent PNG
    plt.savefig(img_output, bbox_inches='tight', pad_inches=0, transparent=True)
    plt.close()

    # 6. Create KMZ
    print("Packaging into Google Earth KMZ...")
    kml = simplekml.Kml()
    
    # Calculate Lat/Long bounds for Google Earth
    bounds_poly = gpd.GeoSeries([box(x_min, y_min, x_max, y_max)], crs=epsg)
    bounds_wgs84 = bounds_poly.to_crs(epsg=4326)
    wgs_minx, wgs_miny, wgs_maxx, wgs_maxy = bounds_wgs84.total_bounds

    ground = kml.newgroundoverlay(name="EM Interpretation")
    ground.icon.href = img_output # Use relative path
    ground.latlonbox.north = wgs_maxy
    ground.latlonbox.south = wgs_miny
    ground.latlonbox.east = wgs_maxx
    ground.latlonbox.west = wgs_minx
    ground.color = "ccffffff" # 80% opacity

    kmz_output = base_path + "_Interpretation.kmz"
    kml.savekmz(kmz_output)
    
    print(f"Success! Saved to: {kmz_output}")
    print(f"(PNG image left at: {img_output} for inspection)")

if __name__ == "__main__":
    main()
