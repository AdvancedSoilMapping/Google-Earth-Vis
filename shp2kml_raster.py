import geopandas as gpd
import simplekml
import os
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from shapely.geometry import box  # <--- This was missing!

# --- CONFIGURATION ---
DEFAULT_COLUMN = "EM" 
DEFAULT_PIXEL_SIZE = 10.0  # Resolution (meters)
DEFAULT_MIN_CUTOFF = 64.0
DEFAULT_MAX_CUTOFF = 326.0

# Define the Color Ramp (Value 0.0 to 1.0)
# You can tweak these positions to shift the color transition
COLORS = [
    (0.0, "red"),           # Low values
    (0.15, "blue"),         # Transition
    (0.30, "orange"),       
    (0.50, "yellow"),       
    (0.70, "lime"),         
    (1.0, "darkgreen")      # High values
]

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
    # Fix extension issues
    base_path, ext = os.path.splitext(shp_path)
    if ext.lower() == '.dbf': shp_path = base_path + ".shp"
    
    if not os.path.exists(shp_path):
        print(f"File not found: {shp_path}")
        return

    col_name = input(f"Column Name [Default: {DEFAULT_COLUMN}]: ").strip() or DEFAULT_COLUMN
    resolution = get_input("Output Resolution (meters/pixel)", DEFAULT_PIXEL_SIZE)

    # 2. Load Data
    print("Reading Data...")
    try:
        gdf = gpd.read_file(shp_path)
        gdf.set_crs(epsg=epsg, allow_override=True, inplace=True)
        
        if col_name not in gdf.columns:
            print(f"Column '{col_name}' not found. Available: {list(gdf.columns)}")
            return

        # Filter Data
        gdf = gdf[(gdf[col_name] >= DEFAULT_MIN_CUTOFF) & (gdf[col_name] <= DEFAULT_MAX_CUTOFF)]
        print(f"Points to process: {len(gdf)}")

    except Exception as e:
        print(f"Error reading shapefile: {e}")
        return
    
    # 3. Rasterize (Create Grid)
    print("Creating Interpolated Image (Grid)...")
    
    # Extract X, Y, Z
    x = gdf.geometry.x.values
    y = gdf.geometry.y.values
    z = gdf[col_name].values

    # Define Grid Bounds (add a small buffer so points aren't on the edge)
    x_min, x_max = x.min(), x.max()
    y_min, y_max = y.min(), y.max()
    
    # Create Grid coordinates
    # complex number 'j' in mgrid allows specifying the number of steps, 
    # but here we use float step size
    grid_x, grid_y = np.mgrid[x_min:x_max:resolution, y_min:y_max:resolution]
    
    # Interpolate
    # method='linear' looks smoother. method='nearest' looks like ArcView blocks.
    grid_z = griddata((x, y), z, (grid_x, grid_y), method='linear')

    # 4. Generate PNG Image
    print("Generating PNG...")
    img_output = base_path + "_Overlay.png"
    
    # Create Custom Colormap
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_em", COLORS)
    
    # Save the image
    # origin='lower' ensures the map isn't upside down
    # .T (Transpose) is needed because mgrid returns (x, y) but imshow expects (row, col)
    plt.figure(figsize=(10, 10), frameon=False)
    plt.imshow(grid_z.T, cmap=cmap, vmin=DEFAULT_MIN_CUTOFF, vmax=DEFAULT_MAX_CUTOFF, origin='lower')
    plt.axis('off')
    plt.savefig(img_output, bbox_inches='tight', pad_inches=0, transparent=True)
    plt.close()

    # 5. Create KML GroundOverlay
    print("Creating Google Earth Overlay...")
    kml = simplekml.Kml()
    
    # Calculate Lat/Long bounds for Google Earth placement
    # We create a box using the min/max of the original data
    bounds_poly = gpd.GeoSeries([box(x_min, y_min, x_max, y_max)], crs=epsg)
    bounds_wgs84 = bounds_poly.to_crs(epsg=4326)
    wgs_minx, wgs_miny, wgs_maxx, wgs_maxy = bounds_wgs84.total_bounds

    # Add the image to KML
    ground = kml.newgroundoverlay(name="EM Survey Heatmap")
    ground.icon.href = os.path.basename(img_output) 
    ground.latlonbox.north = wgs_maxy
    ground.latlonbox.south = wgs_miny
    ground.latlonbox.east = wgs_maxx
    ground.latlonbox.west = wgs_minx
    ground.color = "ccffffff" # 'cc' = approx 80% opacity

    # 6. Save KMZ
    kmz_output = base_path + "_ImageOverlay.kmz"
    kml.savekmz(kmz_output)
    
    # Clean up temp PNG
    if os.path.exists(img_output):
        os.remove(img_output)

    print(f"Success! Overlay saved to: {kmz_output}")

if __name__ == "__main__":
    main()
