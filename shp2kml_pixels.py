import geopandas as gpd
import simplekml
import os
from shapely.geometry import box

# --- CONFIGURATION ---
DEFAULT_COLUMN = "EM" 
DEFAULT_PIXEL_SIZE = 10.0  # Meters (Usually your survey line spacing)
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
    print("\n--- ArcView 3.1 'Flat Shaded' Pixel Generator ---\n")

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
    
    # --- NEW: BOUNDARY INPUT ---
    print("\n(Optional) Clipping")
    bound_path = input("Full path to BOUNDARY Shapefile (Press Enter to skip): ").strip().replace('"', '')
    
    # 3. Settings
    col_name = input(f"Column Name to plot [Default: {DEFAULT_COLUMN}]: ").strip() or DEFAULT_COLUMN
    pixel_size = get_input("Pixel Size in Meters (Width of each square)", DEFAULT_PIXEL_SIZE)
    min_cut = get_input("Bottom Cutoff", DEFAULT_MIN_CUTOFF)
    max_cut = get_input("Top Cutoff", DEFAULT_MAX_CUTOFF)

    # 4. Load Data
    print("\nReading Data...")
    try:
        gdf = gpd.read_file(shp_path)
        if col_name not in gdf.columns:
            print(f"Error: Column '{col_name}' not found.")
            return
        
        # Force Projection
        gdf.set_crs(epsg=epsg, allow_override=True, inplace=True)
        
        # --- STEP 5: CLIPPING ---
        if bound_path and os.path.exists(bound_path):
            print("Reading Boundary and Clipping...")
            bound_gdf = gpd.read_file(bound_path)
            # Ensure boundary is in same projection
            if bound_gdf.crs is None:
                bound_gdf.set_crs(epsg=epsg, inplace=True)
            else:
                bound_gdf = bound_gdf.to_crs(epsg=epsg)
            
            # Perform the Clip
            gdf = gpd.clip(gdf, bound_gdf)
            print(f"Data clipped to boundary. Remaining points: {len(gdf)}")

        # Filter Values
        gdf = gdf[(gdf[col_name] >= min_cut) & (gdf[col_name] <= max_cut)].copy()
        
        # --- STEP 6: CREATE SQUARES (BUFFER) ---
        # We buffer by half the size to get the full width, cap_style=3 makes it a square
        print(f"Converting points to {pixel_size}m squares...")
        gdf['geometry'] = gdf.geometry.buffer(pixel_size / 2, cap_style=3)

        # Reproject to Lat/Long (WGS84) AFTER creating squares to maintain shape size
        gdf_wgs84 = gdf.to_crs(epsg=4326)

    except Exception as e:
        print(f"Error: {e}")
        return

    # 7. Setup KML (Polygons)
    kml = simplekml.Kml()
    
    # Helper to get color hex based on value
    def get_color_hex(val):
        breaks = DEFAULT_BREAKS
        # KML Colors are AABBGGRR (Alpha, Blue, Green, Red)
        # We use Alpha 'ff' (Solid) or 'cc' (Slightly transparent)
        if val <= breaks[0]: return "cc0000ff"   # Red
        if val <= breaks[1]: return "cc0045ff"   # Red-Orange
        if val <= breaks[2]: return "cc008cff"   # Dk Orange
        if val <= breaks[3]: return "cc00a5ff"   # Orange
        if val <= breaks[4]: return "cc00ffff"   # Yellow
        if val <= breaks[5]: return "cc2fffad"   # Yellow-Green
        if val <= breaks[6]: return "cc00ff00"   # Green
        return "cc006400"                        # Dk Green

    print("Generating Surface Map...")
    for index, row in gdf_wgs84.iterrows():
        val = row[col_name]
        
        # Create Polygon
        poly = kml.newpolygon(name=str(val))
        
        # Extract exterior coordinates of the square
        if row.geometry.geom_type == 'Polygon':
            poly.outerboundaryis = list(row.geometry.exterior.coords)
        
        # Style it
        color = get_color_hex(val)
        poly.style.polystyle.color = color
        poly.style.polystyle.outline = 0  # No border line (smooth look)
        poly.description = f"Value: {val}"

    output_path = shp_path.replace(".shp", "_PixelMap.kmz")
    kml.savekmz(output_path)
    print(f"Success! Saved to: {output_path}")

if __name__ == "__main__":
    main()
