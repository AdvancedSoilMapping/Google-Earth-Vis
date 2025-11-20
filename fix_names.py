import os

def main():
    folder = "Test1"
    print(f"Scanning folder: {folder}...")
    
    for filename in os.listdir(folder):
        old_path = os.path.join(folder, filename)
        
        # Convert to lowercase
        new_filename = filename.lower()
        new_path = os.path.join(folder, new_filename)
        
        if old_path != new_path:
            print(f"Renaming: {filename} -> {new_filename}")
            os.rename(old_path, new_path)
            
    print("Done! All files are now lowercase.")

if __name__ == "__main__":
    main()
