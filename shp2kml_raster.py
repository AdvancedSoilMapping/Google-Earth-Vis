import geopandas as gpd
import simplekml
import os
import numpy as np
import rasterio
from rasterio.transform import from_origin
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# --- CONFIGURATION ---
DEFAULT_COLUMN = "EM" 
DEFAULT_PIXEL_SIZE = 10.0  # Resolution of the resulting image (meters)
DEFAULT_MIN_CUTOFF = 64.0
DEFAULT_MAX_CUTOFF = 326.0

# Your Custom Color Ramp (Red -> Green)
# We define a custom colormap to match your specific ranges
COLORS = [
    (0.0, "red"),           # < 180
    (0.14, "blue"),         # 180-197 (Approximating Red-Orange/Blue mix from your previous code)
    (0.20, "darkorange"),   # 197-206
    (0.30, "orange"),       # 206-228
    (0.40, "yellow"),       # 228-247
    (0.50, "yellowgreen"),  # 247-260
    (0.60, "lime"),         # 260-272
    (1.0, "darkgreen")      # > 272
]
# Note: Matplotlib gradients are smooth. If you want hard steps, we use a ListedColormap.

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
    print("\n--- ArcView 3.1 Raster Overlay Generator ---\n")

    # 1. Inputs
    datum = input("Enter Datum (GDA94 / GDA2020) [Default: GDA94]: ").strip().upper() or "GDA94"
    zone = input("Enter Zone (54 / 55) [Default: 55]: ").strip() or "55"
    epsg = get_epsg(datum, zone)
    if not epsg: return

    shp_path = input("Full path to DATA Shapefile: ").strip().replace('"', '')
    base_path, ext = os.path.splitext(shp_path)
    if ext.lower() == '.dbf': shp_path = base_path + ".shp"
    
    if not os.path.exists(shp_path):
        print("File not found.")
        return

    col_name = input(f"Column Name [Default: {DEFAULT_COLUMN}]: ").strip() or DEFAULT_COLUMN
    resolution = get_input("Output Resolution (meters/pixel)", DEFAULT_PIXEL_SIZE)

    # 2. Load Data
    print("Reading Data...")
    gdf = gpd.read_file(shp_path)
    gdf.set_crs(epsg=epsg, allow_override=True, inplace=True)
    
    if col_name not in gdf.columns:
        print(f"Column {col_name} not found.")
        return

    # Filter Data
    gdf = gdf[(gdf[col_name] >= DEFAULT_MIN_CUTOFF) & (gdf[col_name] <= DEFAULT_MAX_CUTOFF)]
    
    # 3. Rasterize (Create Grid)
    print("Creating Interpolated Image (Grid)...")
    
    # Extract X, Y, Z
    x = gdf.geometry.x.values
    y = gdf.geometry.y.values
    z = gdf[col_name].values

    # Define Grid Bounds
    x_min, x_max = x.min(), x.max()
    y_min, y_max = y.min(), y.max()
    
    # Create Grid coordinates
    grid_x, grid_y = np.mgrid[x_min:x_max:resolution, y_min:y_max:resolution]
    
    # Interpolate (Nearest Neighbor is fastest and looks like blocks. 'linear' is smoother)
    grid_z = griddata((x, y), z, (grid_x, grid_y), method='nearest')

    # Flip Y axis for image generation
    grid_z = np.flipud(grid_z.T)

    # 4. Generate PNG Image
    print("Generating PNG...")
    img_output = base_path + "_Overlay.png"
    
    # Create Custom Colormap based on your ranges
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_em", COLORS)
    
    # Save the image without axes/borders
    plt.figure(figsize=(10, 10), frameon=False)
    plt.imshow(grid_z, cmap=cmap, vmin=DEFAULT_MIN_CUTOFF, vmax=DEFAULT_MAX_CUTOFF)
    plt.axis('off')
    plt.savefig(img_output, bbox_inches='tight', pad_inches=0, transparent=True)
    plt.close()

    # 5. Create KML GroundOverlay
    print("Creating Google Earth Overlay...")
    kml = simplekml.Kml()
    
    # We need the Lat/Long bounds for Google Earth
    # Create a bounding box polygon
    bounds_poly = gpd.GeoSeries([box(x_min, y_min, x_max, y_max)], crs=epsg)
    bounds_wgs84 = bounds_poly.to_crs(epsg=4326)
    wgs_minx, wgs_miny, wgs_maxx, wgs_maxy = bounds_wgs84.total_bounds

    # Add the image to KML
    ground = kml.newgroundoverlay(name="EM Survey Overlay")
    ground.icon.href = os.path.basename(img_output) # Reference the PNG we just made
    ground.latlonbox.north = wgs_maxy
    ground.latlonbox.south = wgs_miny
    ground.latlonbox.east = wgs_maxx
    ground.latlonbox.west = wgs_minx
    ground.color = "ccffffff" # Slight transparency (Alpha=cc)

    # 6. Save KMZ (Zips the KML and the PNG together)
    kmz_output = base_path + "_ImageOverlay.kmz"
    kml.savekmz(kmz_output)
    
    # Clean up the temporary PNG so it doesn't clutter
    if os.path.exists(img_output):
        os.remove(img_output)

    print(f"Success! Overlay saved to: {kmz_output}")

if __name__ == "__main__":
    main()
