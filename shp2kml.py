import geopandas as gpd
import simplekml
import os

def get_epsg(datum, zone):
    """
    Returns the correct EPSG code based on Australian Datum and Zone.
    """
    # Dictionary mapping (Datum, Zone) to EPSG Code
    epsg_map = {
        ('GDA94', '54'): 28354,
        ('GDA94', '55'): 28355,
        ('GDA2020', '54'): 7854,
        ('GDA2020', '55'): 7855
    }
    
    return epsg_map.get((datum, zone))

def main():
    print("--- ArcView 3.1 to Google Earth Automation ---")
    
    # 1. User Inputs with Defaults
    datum_input = input("Enter Datum (GDA94 / GDA2020) [Default: GDA94]: ").strip().upper()
    if datum_input == "":
        datum_input = "GDA94"
        
    zone_input = input("Enter Zone (54 / 55) [Default: 55]: ").strip()
    if zone_input == "":
        zone_input = "55"

    # Validate EPSG
    epsg_code = get_epsg(datum_input, zone_input)
    if not epsg_code:
        print(f"Error: Combination of {datum_input} and Zone {zone_input} is not supported in this script yet.")
        return

    print(f"Selected System: {datum_input} Zone {zone_input} (EPSG:{epsg_code})")

    # 2. Get Shapefile Path
    shp_path = input("Enter full path to ArcView Shapefile (.shp): ").strip().replace('"', '')
    
    if not os.path.exists(shp_path):
        print("Error: File not found.")
        return

    # 3. Read and Reproject
    print("Reading Shapefile...")
    try:
        # Load the data
        gdf = gpd.read_file(shp_path)
        
        # Assign the projection (ArcView 3.1 .shp files often lack a .prj file, so we force it)
        gdf.set_crs(epsg=epsg_code, allow_override=True, inplace=True)
        
        # Reproject to WGS84 (Latitude/Longitude) for Google Earth
        gdf_wgs84 = gdf.to_crs(epsg=4326)
        
    except Exception as e:
        print(f"Error reading shapefile: {e}")
        return

    # 4. Create KML/KMZ
    print("Generating Google Earth KML...")
    kml = simplekml.Kml()
    
    # Iterate through rows to create points
    for index, row in gdf_wgs84.iterrows():
        # Extract coordinates
        x, y = row.geometry.x, row.geometry.y
        
        # Create point
        # Change 'Hole_ID' to whatever your actual ID field is named in ArcView
        # If you don't know the field name, you can use str(index)
        name_label = str(row.iloc[0]) # Default to taking the first column as the label
        
        pnt = kml.newpoint(name=name_label, coords=[(x, y)])
        
        # Optional: Add all attributes to the description bubble in Google Earth
        description = ""
        for col in gdf_wgs84.columns:
            if col != 'geometry':
                description += f"<b>{col}:</b> {row[col]}<br>"
        pnt.description = description

    # 5. Save Output
    output_path = shp_path.replace(".shp", "_GoogleEarth.kmz")
    kml.savekmz(output_path)
    
    print(f"Success! File saved to: {output_path}")
    
    # Optional: Auto-open Google Earth
    try:
        os.startfile(output_path)
    except AttributeError:
        # MacOS/Linux alternative
        import subprocess
        subprocess.call(('open', output_path))

if __name__ == "__main__":
    main()
