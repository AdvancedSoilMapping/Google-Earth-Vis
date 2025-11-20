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
DEFAULT_BOUND_WIDTH = 5.0  # Width of the black line in meters

# --- CUSTOM COLOR RAMP (Red -> Green) ---
def create_custom_cmap(min_val, max_val):
    breaks = [180, 197, 206, 228, 247, 260, 272]
    span = max_val - min_val
    norm_breaks = [(b - min_val) / span for b in breaks]
    norm_breaks = [max(0, min(1, b)) for b in norm_breaks]

    cdict = [
        (0.0, "red"),
        (norm_breaks[0], "orangered"),
        (norm_breaks[1], "darkorange"),
        (norm_breaks[2], "orange"),
        (norm_breaks[3], "yellow"),
        (norm_breaks[4], "yellowgreen"),
        (norm_breaks[5], "lime"),
        (norm_breaks[6], "green"),
        (1.0, "darkgreen")
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
    print("\n--- EM Survey Interpretation & Visualization System ---\n")

    # 1. Setup
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

    # Boundary File
    print("\n(Required for clipping and drawing the black line)")
    bound_path = input("Full path to BOUNDARY Shapefile: ").strip().replace('"', '')

    # Settings
    col_name = input(f"Column Name [Default: {DEFAULT_COLUMN}]: ").strip() or DEFAULT_COLUMN
    resolution = get_input("Output Resolution (meters/pixel)", DEFAULT_PIXEL_SIZE)
    min_cut = get_input("Bottom Cutoff", DEFAULT_MIN_CUTOFF)
    max_cut = get_input("Top Cutoff", DEFAULT_MAX_CUTOFF)
    bound_width = get_input("Boundary Line Width (meters)", DEFAULT_BOUND_WIDTH)

    # 2. Load Data
    print("Reading Data...")
    try:
        gdf = gpd.read_file(shp_path)
        gdf.set_crs(epsg=epsg, allow_override=True, inplace=True)
        if col_name not in gdf.columns:
            print(f"Column '{col_name}' not found.")
            return
        gdf = gdf[(gdf[col_name] >= min_cut) & (gdf[col_name] <= max_cut)]
    except Exception as e:
        print(f"Error reading data: {e}")
        return
    
    # 3. Generate Raster (Heatmap)
    print("Creating Interpolated Heatmap...")
    x, y, z = gdf.geometry.x.values, gdf.geometry.y.values, gdf[col_name].values
    x_min, x_max = x.min(), x.max()
    y_min, y_max = y.min(), y.max()
    
    width = int(np.ceil((x_max - x_min) / resolution))
    height = int(np.ceil((y_max - y_min) / resolution))
    
    grid_x, grid_y = np.meshgrid(np.linspace(x_min, x_max, width), np.linspace(y_min, y_max, height))
    grid_z = griddata((x, y), z, (grid_x, grid_y), method='linear')

    # 4. Clipping Mask
    if bound_path and os.path.exists(bound_path):
        print("Applying Boundary Clip...")
        try:
            bound_gdf = gpd.read_file(bound_path)
            if bound_gdf.crs is None: bound_gdf.set_crs(epsg=epsg, inplace=True)
            else: bound_gdf = bound_gdf.to_crs(epsg=epsg)
            
            transform = from_bounds(x_min, y_min, x_max, y_max, width, height)
            shapes = [geom for geom in bound_gdf.geometry]
            mask = features.geometry_mask(shapes, transform=transform, invert=True, out_shape=(height, width))
            mask = np.flipud(mask)
            grid_z[~mask] = np.nan
        except Exception as e:
            print(f"Clipping warning: {e}")

    # 5. Generate PNG Image
    print("Generating Overlay Image...")
    img_output = base_path + "_Overlay.png"
    cmap = create_custom_cmap(min_cut, max_cut)
    
    plt.figure(figsize=(10, 10), frameon=False)
    plt.imshow(grid_z, cmap=cmap, vmin=min_cut, vmax=max_cut, origin='lower')
    plt.axis('off')
    plt.savefig(img_output, bbox_inches='tight', pad_inches=0, transparent=True)
    plt.close()

    # 6. Create KML
    print("Constructing Google Earth KML...")
    kml = simplekml.Kml()
    
    # A. Add Ground Overlay (The Heatmap)
    bounds_poly = gpd.GeoSeries([box(x_min, y_min, x_max, y_max)], crs=epsg)
    bounds_wgs84 = bounds_poly.to_crs(epsg=4326)
    wgs_minx, wgs_miny, wgs_maxx, wgs_maxy = bounds_wgs84.total_bounds

    ground = kml.newgroundoverlay(name="EM Interpretation")
    ground.icon.href = img_output
    ground.latlonbox.north = wgs_maxy
    ground.latlonbox.south = wgs_miny
    ground.latlonbox.east = wgs_maxx
    ground.latlonbox.west = wgs_minx
    ground.color = "ccffffff" # Semi-transparent

    # B. Add Boundary Line (as a Thick Solid Polygon)
    if bound_path and os.path.exists(bound_path):
        print(f"Drawing {bound_width}m wide Boundary Line...")
        
        # Buffer the boundary to make it thick
        # buffer(width/2) makes it grow outwards by half the width in both directions
        thick_bound = bound_gdf.copy()
        thick_bound['geometry'] = thick_bound.geometry.buffer(bound_width / 2)
        
        # Reproject to Lat/Long for KML
        thick_bound_wgs84 = thick_bound.to_crs(epsg=4326)
        
        for geom in thick_bound_wgs84.geometry:
            # Handle MultiPolygons
            polys = [geom] if geom.geom_type == 'Polygon' else geom.geoms
            
            for poly in polys:
                kml_poly = kml.newpolygon(name="Boundary")
                kml_poly.outerboundaryis = list(poly.exterior.coords)
                
                # Solid Black (ff = Alpha, 000000 = Black)
                kml_poly.style.polystyle.color = "ff000000" 
                kml_poly.style.polystyle.outline = 0

    # 7. Save KMZ
    kmz_output = base_path + "_Map_V1.kmz"
    kml.savekmz(kmz_output)
    
    print(f"Success! Saved to: {kmz_output}")

if __name__ == "__main__":
    main()
