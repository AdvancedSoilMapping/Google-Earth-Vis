import geopandas as gpd
import simplekml
import os
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.path as mpath
from shapely.geometry import Polygon, Point
from shapely.ops import unary_union

# --- CONFIGURATION ---
DEFAULT_COLUMN = "EM" 
DEFAULT_PIXEL_SIZE = 10.0 
DEFAULT_MIN_CUTOFF = 64.0
DEFAULT_MAX_CUTOFF = 326.0
DEFAULT_BOUND_WIDTH = 3.0

# --- HELPERS ---
def ring_coding(ob):
    n = len(ob.coords)
    codes = np.ones(n, dtype=mpath.Path.code_type) * mpath.Path.LINETO
    codes[0] = mpath.Path.MOVETO
    return codes

def pathify(polygon):
    path_data = []
    if polygon.geom_type == 'Polygon':
        for ring in [polygon.exterior] + list(polygon.interiors):
            path_data.append((np.array(ring.coords), ring_coding(ring)))
    elif polygon.geom_type == 'MultiPolygon':
        for poly in polygon.geoms:
            for ring in [poly.exterior] + list(poly.interiors):
                path_data.append((np.array(ring.coords), ring_coding(ring)))
    vertices = np.concatenate([p[0] for p in path_data])
    codes = np.concatenate([p[1] for p in path_data])
    return mpath.Path(vertices, codes)

def create_custom_cmap(min_val, max_val):
    breaks = [180, 197, 206, 228, 247, 260, 272]
    span = max_val - min_val
    norm_breaks = [(b - min_val) / span for b in breaks]
    norm_breaks = [max(0, min(1, b)) for b in norm_breaks]
    cdict = [
        (0.0, "red"), (norm_breaks[0], "orangered"), (norm_breaks[1], "darkorange"),
        (norm_breaks[2], "orange"), (norm_breaks[3], "yellow"), (norm_breaks[4], "yellowgreen"),
        (norm_breaks[5], "lime"), (norm_breaks[6], "green"), (1.0, "darkgreen")
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
    print("\n--- EM Survey Perfect Projection (Quad) Generator ---\n")

    # 1. Setup
    datum = input("Enter Datum (GDA94 / GDA2020) [Default: GDA94]: ").strip().upper() or "GDA94"
    zone = input("Enter Zone (54 / 55) [Default: 55]: ").strip() or "55"
    epsg = get_epsg(datum, zone)
    if not epsg: return

    shp_path = input("Full path to DATA Shapefile: ").strip().replace('"', '')
    base_path, ext = os.path.splitext(shp_path)
    if ext.lower() == '.dbf': shp_path = base_path + ".shp"
    if not os.path.exists(shp_path):
        print(f"File not found: {shp_path}")
        return

    bound_path = input("Full path to BOUNDARY Shapefile: ").strip().replace('"', '')
    if not os.path.exists(bound_path):
        print("Error: Boundary file required.")
        return

    # Settings
    col_name = input(f"Column Name [Default: {DEFAULT_COLUMN}]: ").strip() or DEFAULT_COLUMN
    resolution = get_input("Output Resolution (meters/pixel)", DEFAULT_PIXEL_SIZE)
    min_cut = get_input("Bottom Cutoff", DEFAULT_MIN_CUTOFF)
    max_cut = get_input("Top Cutoff", DEFAULT_MAX_CUTOFF)
    bound_width = get_input("Boundary Line Thickness", DEFAULT_BOUND_WIDTH)

    # 2. Load Data (Keep in MGA for processing)
    print("Reading Data (Processing in Grid Coordinates)...")
    try:
        gdf = gpd.read_file(shp_path)
        gdf.set_crs(epsg=epsg, allow_override=True, inplace=True)
        if col_name not in gdf.columns: return
        gdf = gdf[(gdf[col_name] >= min_cut) & (gdf[col_name] <= max_cut)]
        
        bound_gdf = gpd.read_file(bound_path)
        if bound_gdf.crs is None: bound_gdf.set_crs(epsg=epsg, inplace=True)
        else: bound_gdf = bound_gdf.to_crs(epsg=epsg)

    except Exception as e:
        print(f"Error: {e}")
        return
    
    # 3. Interpolate in MGA (This ensures colors and shapes are correct)
    print("Creating Interpolated Heatmap...")
    x, y, z = gdf.geometry.x.values, gdf.geometry.y.values, gdf[col_name].values
    x_min, x_max = x.min(), x.max()
    y_min, y_max = y.min(), y.max()
    
    # Add padding for the line width
    pad = resolution * 2
    x_min -= pad; x_max += pad
    y_min -= pad; y_max += pad
    
    width = int(np.ceil((x_max - x_min) / resolution))
    height = int(np.ceil((y_max - y_min) / resolution))
    
    grid_x, grid_y = np.meshgrid(np.linspace(x_min, x_max, width), np.linspace(y_min, y_max, height))
    grid_z = griddata((x, y), z, (grid_x, grid_y), method='linear')

    # 4. Vector Clip
    print("Applying Boundary Clip...")
    combined_poly = unary_union(bound_gdf.geometry)
    clip_path = pathify(combined_poly)
    clip_patch = mpatches.PathPatch(clip_path, facecolor='none', edgecolor='none')

    # 5. Render Image
    print("Rendering Composite Image...")
    img_output = base_path + "_Overlay.png"
    cmap = create_custom_cmap(min_cut, max_cut)
    
    fig, ax = plt.subplots(figsize=(10, 10), frameon=False)
    
    # Plot Heatmap
    im = ax.imshow(grid_z, cmap=cmap, vmin=min_cut, vmax=max_cut, 
                   origin='lower', interpolation='bilinear',
                   extent=[x_min, x_max, y_min, y_max])
    
    # Apply Clip
    ax.add_patch(clip_patch)
    im.set_clip_path(clip_patch)
    
    # Plot Black Boundary
    bound_gdf.plot(ax=ax, color='none', edgecolor='black', linewidth=bound_width)
    
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_axis_off()
    plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
    plt.savefig(img_output, dpi=150, transparent=True, pad_inches=0)
    plt.close()

    # 6. Calculate Quad Coordinates (The Anti-Skew Fix)
    print("Calculating Projection Corners...")
    
    # Create 4 corner points in MGA
    corners_mga = [
        Point(x_min, y_min), # Bottom-Left
        Point(x_max, y_min), # Bottom-Right
        Point(x_max, y_max), # Top-Right
        Point(x_min, y_max)  # Top-Left
    ]
    
    # Convert these 4 specific points to WGS84 Lat/Long
    corners_gdf = gpd.GeoSeries(corners_mga, crs=epsg)
    corners_wgs84 = corners_gdf.to_crs(epsg=4326)
    
    # Extract coords for SimpleKML (Longitude, Latitude)
    # Order: BL, BR, TR, TL
    quad_coords = [
        (corners_wgs84[0].x, corners_wgs84[0].y),
        (corners_wgs84[1].x, corners_wgs84[1].y),
        (corners_wgs84[2].x, corners_wgs84[2].y),
        (corners_wgs84[3].x, corners_wgs84[3].y)
    ]

    # 7. Create KML using gx:LatLonQuad
    print("Creating Google Earth KML (Quad-Warped)...")
    kml = simplekml.Kml()
    
    ground = kml.newgroundoverlay(name="EM Interpretation")
    ground.icon.href = img_output
    ground.color = "ffffffff"
    
    # Instead of LatLonBox, we use gxLatLonQuad to handle rotation/skew perfectly
    ground.gxlatlonquad.coords = quad_coords

    kmz_output = base_path + "_Quad.kmz"
    kml.savekmz(kmz_output)
    
    print(f"Success! Saved to: {kmz_output}")

if __name__ == "__main__":
    main()
