import geopandas as gpd
import simplekml
import os
from shapely.geometry import box

# --- CONFIGURATION ---
DEFAULT_COLUMN = "EM" 
DEFAULT_PIXEL_SIZE = 25.0  # INCREASED DEFAULT (Try matching your line spacing)
DEFAULT_MIN_CUTOFF = 64.0
DEFAULT_MAX_CUTOFF = 326.0
DEFAULT_BREAKS = [180, 197, 206, 228, 247, 260, 272, 326]

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
    print("\n--- ArcView 3.1 Final Surface Generator ---\n")

    # 1. Coordinate System
    datum = input("Enter Datum (GDA94 / GDA2020) [Default: GDA94]: ").strip().upper() or "GDA94"
    zone = input("Enter Zone (54 / 55) [Default: 55]: ").strip() or "55"
    epsg = get_epsg(datum, zone)
    if not epsg: return

    # 2. Input Files
    shp_path = input("Full path to DATA Shapefile: ").strip().replace('"', '')
    if not os.path.exists(shp_path):
        print("Error: Data file not found.")
        return
    
    print("\n(Required for clean edges)")
    bound_path = input("Full path to BOUNDARY Shapefile: ").strip().replace('"', '')
    
    # 3. Settings
    col_name = input(f"Column Name to plot [Default: {DEFAULT_COLUMN}]: ").strip() or DEFAULT_COLUMN
    pixel_size = get_input("Pixel Size (Should match line spacing)", DEFAULT_PIXEL_SIZE)
    min_cut = get_input("Bottom Cutoff", DEFAULT_MIN_CUTOFF)
    max_cut = get_input("Top Cutoff", DEFAULT_MAX_CUTOFF)

    # 4. Load and Process
    print("\nProcessing Geometry...")
    try:
        # Load Data
        gdf = gpd.read_file(shp_path)
        gdf.set_crs(epsg=epsg, allow_override=True, inplace=True)
        
        # Filter Values
        if col_name not in gdf.columns:
            print(f"Error: Column '{col_name}' not found.")
            return
        gdf = gdf[(gdf[col_name] >= min_cut) & (gdf[col_name] <= max_cut)].copy()

        # STEP A: EXPAND (Create Squares)
        # We assume square width = pixel_size. Buffer radius is half that.
        # cap_style=3 creates a square buffer.
        print(f"Creating {pixel_size}m squares...")
        gdf['geometry'] = gdf.geometry.buffer(pixel_size / 2, cap_style=3)

        # STEP B: CLIP (Cookie Cutter)
        if bound_path and os.path.exists(bound_path):
            print("Clipping squares to boundary...")
            bound_gdf = gpd.read_file(bound_path)
            # Force boundary to same projection
            bound_gdf.set_crs(epsg=epsg, allow_override=True, inplace=True)
            
            # Clip the SQUARES against the BOUNDARY
            gdf = gpd.clip(gdf, bound_gdf)
            print(f"Clip complete. Features remaining: {len(gdf)}")
        else:
            print("Warning: No boundary found. Edges will be jagged.")

        # Reproject to WGS84 for Google Earth
        gdf_wgs84 = gdf.to_crs(epsg=4326)

    except Exception as e:
        print(f"Error: {e}")
        return

    # 5. Generate KML
    kml = simplekml.Kml()
    
    def get_color_hex(val):
        breaks = DEFAULT_BREAKS
        # Solid Colors (Alpha ff) work best for flat maps to hide overlaps
        if val <= breaks[0]: return "ff0000ff"   # Red
        if val <= breaks[1]: return "ff0045ff"   # Red-Orange
        if val <= breaks[2]: return "ff008cff"   # Dk Orange
        if val <= breaks[3]: return "ff00a5ff"   # Orange
        if val <= breaks[4]: return "ff00ffff"   # Yellow
        if val <= breaks[5]: return "ff2fffad"   # Yellow-Green
        if val <= breaks[6]: return "ff00ff00"   # Green
        return "ff006400"                        # Dk Green

    print("Writing KML...")
    for index, row in gdf_wgs84.iterrows():
        val = row[col_name]
        
        # Handle MultiPolygons (result of complex clips)
        geoms = [row.geometry] if row.geometry.geom_type == 'Polygon' else row.geometry.geoms
        
        for geom in geoms:
            poly = kml.newpolygon(name=str(val))
            poly.outerboundaryis = list(geom.exterior.coords)
            poly.style.polystyle.color = get_color_hex(val)
            poly.style.polystyle.outline = 0  # No borders
            poly.description = f"Value: {val}"

    output_path = shp_path.replace(".shp", "_FinalMap.kmz")
    kml.savekmz(output_path)
    print(f"Success! Saved to: {output_path}")

if __name__ == "__main__":
    main()
