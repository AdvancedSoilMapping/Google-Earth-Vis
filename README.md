### **The Solution: Automated ArcView-to-Google Earth Script**
Here is the complete, ready-to-run Python script.

**What this script does:**
1.  Asks you for the **Datum** (defaults to GDA94).
2.  Asks you for the **Zone** (defaults to Zone 55).
3.  Asks for your **Input Shapefile** path (exported from ArcView 3.1).
4.  Automatically selects the correct Coordinate Reference System (EPSG code).
5.  Converts the coordinates to WGS84 (Lat/Long).
6.  Generates a `Drillholes.kmz` file that opens in Google Earth.

#### **Step 1: Install the required Python libraries**
You will need to install Python, then open your terminal/command prompt and run:
```bash
pip install geopandas simplekml
```
*(Note: `geopandas` can sometimes be tricky to install on Windows. If `pip install` fails, the easiest way is to install [Anaconda](https://www.anaconda.com/download) and use `conda install geopandas`)*.

#### **Step 2: The Script (`shp2kml.py`)**
Save this code into a file named `shp2kml.py` in your new project folder.

```python
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
```

### **How to set up the GitHub Repository**
Since you want to do this for R&D:

1.  Go to [GitHub.com](https://github.com) and sign up/log in.
2.  Click the **+** icon (top right) -> **New repository**.
3.  Name it: `ArcView-EM-Automation`.
4.  Select **Public** (for open source) or **Private**.
5.  Check **Add a README file**.
6.  Under ".gitignore", select **Python** (this prevents junk files from being uploaded).
7.  Click **Create repository**.

**To get the code in there:**
1.  On your computer, create a folder named `ArcView-EM-Automation`.
2.  Save the script above as `shp2kml.py` in that folder.
3.  If you have Git installed:
    ```bash
    git init
    git remote add origin https://github.com/YOUR_USERNAME/ArcView-EM-Automation.git
    git add .
    git commit -m "Initial commit of conversion script"
    git push -u origin main
    ```
    *(If you don't like command lines, you can just use the "Upload files" button on the GitHub website to drag and drop your script there).*
