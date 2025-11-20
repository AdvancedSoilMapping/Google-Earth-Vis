import geopandas as gpd
import simplekml
import os
import pandas as pd

# --- CONFIGURATION: DEFAULT VALUES ---
DEFAULT_COLUMN = "EM"
DEFAULT_MIN_CUTOFF = 64.0
DEFAULT_MAX_CUTOFF = 326.0

# These are the upper boundaries of your ranges.
# Range 1: Min to 180
# Range 2: 180 to 197, etc...
DEFAULT_BREAKS = [180, 197, 206, 228, 247, 260, 272, 326]

def get_input(prompt, default_val):
    """Helper to get input with a default value"""
    user_val = input(f"{prompt} [Default: {default_val}]: ").strip()
    if user_val == "":
        return default_val
    try:
        return type(default_val)(user_val)
    except ValueError:
        print(f"Invalid input. Using default: {default_val}")
        return default_val

def get_epsg(datum, zone):
    epsg_map = {
        ('GDA94', '54'): 28354, ('GDA94', '55'): 28355,
        ('GDA2020', '54'): 7854, ('GDA2020', '55'): 7855
    }
    return epsg_map.get((datum, zone))

def main():
    print("\n--- ArcView 3.1 EM Color Interpretation Script ---\n")

    # 1. Coordinate System Setup
    datum = input("Enter Datum (GDA94 / GDA2020) [Default: GDA94]: ").strip().upper() or "GDA94"
    zone = input("Enter Zone (54 / 55) [Default: 55]: ").strip() or "55"
    epsg = get_epsg(datum, zone)
    
    if not epsg:
        print("Invalid Datum/Zone combination.")
        return

    # 2. File Input
    shp_path = input("Enter full path to Shapefile (.shp): ").strip().replace('"', '')
    if not os.path.exists(shp_path):
        print("Error: File not found.")
        return

    # 3. Data Column & Cutoffs
    col_name = input(f"Enter Column Name to plot [Default: {DEFAULT_COLUMN}]: ").strip() or DEFAULT_COLUMN
    min_cut = get_input("Enter Bottom Cutoff (Discard below)", DEFAULT_MIN_CUTOFF)
    max_cut = get_input("Enter Top Cutoff (Discard above)", DEFAULT_MAX_CUTOFF)

    # 4. Define Ranges (Breaks)
    print("\n--- Define Color Ranges ---")
    print(f"Current default breaks: {DEFAULT_BREAKS}")
    use_defaults = input("Use default color breaks? (y/n) [Default: y]: ").strip().lower()
    
    breaks = DEFAULT_BREAKS
    if use_defaults == 'n':
        breaks = []
        print("Enter the UPPER limit for each of the 8 ranges.")
        for i in range(1, 9):
            val = float(input(f"Enter Upper Limit for Range {i}: "))
            breaks.append(val)

    # 5. Process Data
    print("\nReading and Filtering Data...")
    try:
        gdf = gpd.read_file(shp_path)
        
        # Check if column exists
        if col_name not in gdf.columns:
            print(f"Error: Column '{col_name}' not found. Available columns: {list(gdf.columns)}")
            return

        # Force Projection
        gdf.set_crs(epsg=epsg, allow_override=True, inplace=True)
        gdf_wgs84 = gdf.to_crs(epsg=4326)

        # Filter Data (The Cutoff)
        initial_count = len(gdf_wgs84)
        gdf_filtered = gdf_wgs84[
            (gdf_wgs84[col_name] >= min_cut) & 
            (gdf_wgs84[col_name] <= max_cut)
        ].copy() # Use copy to avoid warnings
        
        print(f"Filtered {initial_count} points down to {len(gdf_filtered)} points.")

    except Exception as e:
        print(f"Error processing data: {e}")
        return

    # 6. Setup KML Styles (Shared Styles for Performance)
    kml = simplekml.Kml()
    
    # Define styles matching your ArcView setup
    # KML Color Format: AABBGGRR (Alpha, Blue, Green, Red) - Yes, it's backwards.
    styles = []
    
    # Range 1: Red (67-180)
    s1 = simplekml.Style(); s1.iconstyle.color = simplekml.Color.red; styles.append(s1)
    # Range 2: Red-Orange (180-197) -> Blue=00, Green=45, Red=FF
    s2 = simplekml.Style(); s2.iconstyle.color = "ff0045ff"; styles.append(s2)
    # Range 3: Dark Orange (197-206)
    s3 = simplekml.Style(); s3.iconstyle.color = simplekml.Color.darkorange; styles.append(s3)
    # Range 4: Orange (206-228)
    s4 = simplekml.Style(); s4.iconstyle.color = simplekml.Color.orange; styles.append(s4)
    # Range 5: Yellow (228-247)
    s5 = simplekml.Style(); s5.iconstyle.color = simplekml.Color.yellow; styles.append(s5)
    # Range 6: Yellow-Green (247-260)
    s6 = simplekml.Style(); s6.iconstyle.color = simplekml.Color.yellowgreen; styles.append(s6)
    # Range 7: Green (260-272)
    s7 = simplekml.Style(); s7.iconstyle.color = simplekml.Color.lime; styles.append(s7)
    # Range 8: Dark Green (272-326)
    s8 = simplekml.Style(); s8.iconstyle.color = simplekml.Color.darkgreen; styles.append(s8)

    # Apply common settings to all styles (make them small circles)
    for s in styles:
        s.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png'
        s.iconstyle.scale = 0.6  # Scale down the dot size (0.6 is good for density)
        s.labelstyle.scale = 0   # Hide labels by default (so it looks like a map, not text)

    # 7. Generate Points
    print("Generating Color Interpretation Map...")
    
    # Define bins for classification
    # 0   1     2     3     4     5     6     7
    # Min -> B1 -> B2 -> B3 -> B4 -> B5 -> B6 -> B7 -> B8
    
    for index, row in gdf_filtered.iterrows():
        val = row[col_name]
        
        # Determine which bin (style) to use
        style_to_use = None
        
        if val <= breaks[0]: style_to_use = styles[0]      # Red
        elif val <= breaks[1]: style_to_use = styles[1]    # Red-Orange
        elif val <= breaks[2]: style_to_use = styles[2]    # Dk Orange
        elif val <= breaks[3]: style_to_use = styles[3]    # Orange
        elif val <= breaks[4]: style_to_use = styles[4]    # Yellow
        elif val <= breaks[5]: style_to_use = styles[5]    # Yellow-Green
        elif val <= breaks[6]: style_to_use = styles[6]    # Green
        elif val <= breaks[7]: style_to_use = styles[7]    # Dk Green
        else: style_to_use = styles[7] # Catch-all for top edge
        
        # Create Point
        pnt = kml.newpoint(coords=[(row.geometry.x, row.geometry.y)])
        pnt.style = style_to_use
        pnt.description = f"Value: {val}" # Pop-up will show value
        
    # 8. Save
    output_path = shp_path.replace(".shp", "_Interpretation.kmz")
    kml.savekmz(output_path)
    print(f"Success! Interpretation saved to: {output_path}")

if __name__ == "__main__":
    main()
