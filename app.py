import streamlit as st
import geopandas as gpd
import simplekml
import os
import shutil
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.path as mpath
from shapely.ops import unary_union

# --- CONFIGURATION & SETUP ---
st.set_page_config(page_title="EM Survey Processor", layout="wide")

# Custom CSS to hide the "Made with Streamlit" footer
st.markdown("""<style>
    footer {visibility: hidden;}
    .stDeployButton {display:none;}
</style>""", unsafe_allow_html=True)

# --- HELPER FUNCTIONS ---

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
    # Avoid division by zero
    if span == 0: span = 1
    
    norm_breaks = [(b - min_val) / span for b in breaks]
    norm_breaks = [max(0, min(1, b)) for b in norm_breaks]

    cdict = [
        (0.0, "red"), (norm_breaks[0], "orangered"), (norm_breaks[1], "darkorange"),
        (norm_breaks[2], "orange"), (norm_breaks[3], "yellow"), (norm_breaks[4], "yellowgreen"),
        (norm_breaks[5], "lime"), (norm_breaks[6], "green"), (1.0, "darkgreen")
    ]
    return mcolors.LinearSegmentedColormap.from_list("em_ramp", cdict)

def save_uploaded_files(uploaded_files):
    # Streamlit keeps files in RAM. We must save them to disk for GeoPandas.
    if os.path.exists("temp_upload"):
        shutil.rmtree("temp_upload")
    os.makedirs("temp_upload")
    
    shp_path = None
    
    for uploaded_file in uploaded_files:
        file_path = os.path.join("temp_upload", uploaded_file.name)
        with open(file_path, "wb") as f:
            f.write(uploaded_file.getbuffer())
        if uploaded_file.name.lower().endswith(".shp"):
            shp_path = file_path
            
    return shp_path

# --- MAIN APP INTERFACE ---

st.title("🇦🇺 EM Survey Visualization Tool")
st.markdown("Convert ArcView 3.1 EM Surveys to Google Earth (KMZ) with separate layers and smooth clipping.")

col1, col2 = st.columns(2)

with col1:
    st.header("1. Data Inputs")
    
    # Coordinate System
    datum = st.selectbox("Datum", ["GDA94", "GDA2020"], index=0)
    zone = st.selectbox("Zone", ["54", "55"], index=1)
    
    # File Uploader
    st.info("Upload ALL shapefile parts (.shp, .shx, .dbf) for BOTH Data and Boundary.")
    uploaded_data = st.file_uploader("Upload EM Data Files", accept_multiple_files=True, type=['shp','shx','dbf'])
    uploaded_bound = st.file_uploader("Upload Boundary Files", accept_multiple_files=True, type=['shp','shx','dbf'])

with col2:
    st.header("2. Settings")
    col_name = st.text_input("Column Name", "EM")
    
    # Sliders for Cutoffs
    range_val = st.slider("Data Filter Range", 0.0, 500.0, (64.0, 326.0))
    min_cut, max_cut = range_val
    
    pixel_size = st.number_input("Interpolation Grid Size (Meters)", value=10.0, min_value=1.0)
    bound_width = st.number_input("Boundary Line Width", value=3.0)
    
    # SUPER SAMPLING
    st.markdown("---")
    quality = st.select_slider("Edge Quality (Super Sampling)", options=["Standard", "High", "Ultra"], value="High")
    dpi_map = {"Standard": 100, "High": 300, "Ultra": 600}

# --- PROCESS BUTTON ---
if st.button("Generate Map", type="primary"):
    if not uploaded_data or not uploaded_bound:
        st.error("Please upload both Data and Boundary files (ensure .shp, .shx, and .dbf are included).")
    else:
        # 1. Save files to distinct folders
        data_path = save_uploaded_files(uploaded_data)
        
        # Handle Boundary separately
        if os.path.exists("temp_bound"): shutil.rmtree("temp_bound")
        os.makedirs("temp_bound")
        bound_path = None
        for f in uploaded_bound:
            fp = os.path.join("temp_bound", f.name)
            with open(fp, "wb") as w: w.write(f.getbuffer())
            if f.name.lower().endswith(".shp"): bound_path = fp

        if not data_path or not bound_path:
            st.error("Could not find .shp file in uploads. Did you upload all extensions?")
            st.stop()
            
        with st.spinner("Processing Geometry and WGS84 Conversion..."):
            # EPSG Lookup
            epsg_map = {('GDA94', '54'): 28354, ('GDA94', '55'): 28355,
                        ('GDA2020', '54'): 7854, ('GDA2020', '55'): 7855}
            epsg = epsg_map.get((datum, zone))
            
            # LOAD DATA
            try:
                gdf = gpd.read_file(data_path)
                gdf.set_crs(epsg=epsg, allow_override=True, inplace=True)
                gdf = gdf[(gdf[col_name] >= min_cut) & (gdf[col_name] <= max_cut)]
                gdf = gdf.to_crs(epsg=4326) # Convert to Lat/Long

                bound_gdf = gpd.read_file(bound_path)
                if bound_gdf.crs is None: bound_gdf.set_crs(epsg=epsg, inplace=True)
                else: bound_gdf = bound_gdf.to_crs(epsg=epsg)
                bound_gdf = bound_gdf.to_crs(epsg=4326) # Convert to Lat/Long
            except Exception as e:
                st.error(f"Error reading shapefiles: {e}")
                st.stop()

        with st.spinner("Interpolating Heatmap..."):
            x = gdf.geometry.x.values
            y = gdf.geometry.y.values
            z = gdf[col_name].values

            # Dynamic Grid Calculation
            mean_lat = np.mean(y)
            deg_per_meter_lat = 1.0 / 111320.0
            deg_per_meter_lon = 1.0 / (111320.0 * np.cos(np.radians(mean_lat)))
            res_lat = pixel_size * deg_per_meter_lat
            res_lon = pixel_size * deg_per_meter_lon

            x_min, x_max = x.min(), x.max()
            y_min, y_max = y.min(), y.max()
            
            pad_x = res_lon * 5
            pad_y = res_lat * 5
            x_min -= pad_x; x_max += pad_x
            y_min -= pad_y; y_max += pad_y
            
            width = int(np.ceil((x_max - x_min) / res_lon))
            height = int(np.ceil((y_max - y_min) / res_lat))
            
            grid_x, grid_y = np.meshgrid(np.linspace(x_min, x_max, width), np.linspace(y_min, y_max, height))
            grid_z = griddata((x, y), z, (grid_x, grid_y), method='linear')

        with st.spinner("Rendering High-Res Layers..."):
            # Prepare Clip Path
            combined_poly = unary_union(bound_gdf.geometry)
            clip_path = pathify(combined_poly)
            
            # RENDER IMAGE (LAYER 1)
            img_output = "overlay.png"
            cmap = create_custom_cmap(min_cut, max_cut)
            
            fig = plt.figure(figsize=(10, 10 * (height/width)), frameon=False)
            ax = plt.axes([0, 0, 1, 1])
            ax.set_axis_off()
            
            # Plot Heatmap
            im = ax.imshow(grid_z, cmap=cmap, vmin=min_cut, vmax=max_cut, 
                           origin='lower', interpolation='bilinear',
                           extent=[x_min, x_max, y_min, y_max], aspect='auto')
            
            # Apply Clip (Makes outside transparent)
            clip_patch = mpatches.PathPatch(clip_path, facecolor='none', edgecolor='none', transform=ax.transData)
            ax.add_patch(clip_patch)
            im.set_clip_path(clip_patch)
            
            # IMPORTANT: We do NOT plot the boundary here. We keep them separate.
            
            ax.set_xlim(x_min, x_max)
            ax.set_ylim(y_min, y_max)
            
            # Save with selected DPI (Super Sampling)
            plt.savefig(img_output, dpi=dpi_map[quality], transparent=True, pad_inches=0)
            plt.close()

            # GENERATE KML
            kml = simplekml.Kml()
            
            # Folder 1: Interpretation
            fol_interp = kml.newfolder(name="Interpretation")
            ground = fol_interp.newgroundoverlay(name="Heatmap")
            ground.icon.href = img_output
            ground.latlonbox.north = y_max
            ground.latlonbox.south = y_min
            ground.latlonbox.east = x_max
            ground.latlonbox.west = x_min
            ground.color = "ffffffff" # Fully opaque (transparency handled by PNG)

            # Folder 2: Boundary (Vector Layer)
            fol_bound = kml.newfolder(name="Boundary")
            
            for geom in bound_gdf.geometry:
                polys = [geom] if geom.geom_type == 'Polygon' else geom.geoms
                for poly in polys:
                    # Create LineString (Outline)
                    ls = fol_bound.newlinestring(name="Boundary Line")
                    ls.coords = list(poly.exterior.coords)
                    ls.style.linestyle.color = simplekml.Color.black
                    ls.style.linestyle.width = bound_width
            
            kmz_output = "EM_Survey_Final.kmz"
            kml.savekmz(kmz_output)

        st.success("Processing Complete!")
        
        # Download Button
        with open(kmz_output, "rb") as f:
            st.download_button(
                label="🌏 Download Google Earth KMZ",
                data=f,
                file_name="EM_Survey_Interpretation.kmz",
                mime="application/vnd.google-earth.kmz"
            )
