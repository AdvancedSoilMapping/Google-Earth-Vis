import geopandas as gpd
import os

def main():
    print("\n--- Shapefile Column Inspector ---\n")
    
    # Get file path
    shp_path = input("Enter full path to Shapefile (.shp): ").strip().replace('"', '')
    
    if not os.path.exists(shp_path):
        print("Error: File not found.")
        return

    try:
        # Read file
        gdf = gpd.read_file(shp_path)
        
        print(f"\nSUCCESS! Loaded {len(gdf)} rows.")
        print("\nHERE ARE YOUR COLUMN NAMES:")
        print("=============================")
        for col in gdf.columns:
            # Print the name and an example value from the first row
            first_val = gdf[col].iloc[0]
            print(f"NAME: {col}   (Example data: {first_val})")
        print("=============================\n")
        
    except Exception as e:
        print(f"Error reading file: {e}")

if __name__ == "__main__":
    main()
