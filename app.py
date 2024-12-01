from utility import utils
import streamlit as st

# # Streamlit App
# Tool Title
st.title("Satellite imagery-based Crop Health Monitoring")
# Sidebar section for downloading satellite data
st.sidebar.header("Download Sentinel-2 Satellite Data")

# Button to trigger satellite data download
if st.sidebar.button("Draw Area of Interest (Optional)"):
    basic_map = utils.map_display(aoi = None)
    utils.folium_static(basic_map)

# Date selection widgets in the sidebar
aoi = eval(st.sidebar.text_input("Draw AOI on map and click on the polygon drawn, copy the json and paste here", value={"type":"Feature","properties":{},"geometry":{"type":"Polygon","coordinates":[[[75.711079,29.282731],[75.711079,29.299648],[75.73245,29.299648],[75.73245,29.282731],[75.711079,29.282731]]]}}))
start_date = st.sidebar.date_input("Start Date", utils.datetime(2023, 10, 1), help="Sowing time of a season")
end_date = st.sidebar.date_input("End Date", utils.datetime(2024, 4, 12), help="Harvesting time of a season")
cloud_percentage = st.sidebar.number_input("Cloud Percentage", min_value=0, max_value=100, value=30, help="0-100")

# Button to trigger satellite data download
if st.sidebar.button("Download Satellite Data"):
    st.write("Downloading satellite data...")
    status = utils.Sentinel_2_data_download(start_date=start_date.strftime('%Y-%m-%d'),
                                            end_date=end_date.strftime('%Y-%m-%d'),
                                            cloud_percentage=cloud_percentage, aoi=aoi)
    st.write(status)  # Display download status
    # After download, calculate vegetation indices (calling function from utils)
    utils.calculate_evi_ndvi()

# Directory paths for NDVI files after download
out_dir = utils.os.getcwd()
stacked_raster_dir = utils.os.path.join(out_dir, 'Sentinel2_data')
output_dir = utils.os.path.join(stacked_raster_dir, 'Output')

# Check if output directory exists and list NDVI files
if utils.os.path.exists(output_dir):
    indices_files = [filename for filename in utils.os.listdir(output_dir) if filename.endswith('.tif')]
    indices_files.sort()
else:
    indices_files = []

# Sidebar section for plotting vegetative indices
st.sidebar.header("Plot Vegetative Indices")

if indices_files:
    RasterPlot = []
    selected_file = st.sidebar.selectbox("Select NDVI or EVI file to plot", indices_files)
    ndvi_file_path = utils.os.path.join(output_dir, selected_file)
    vegetative_index = selected_file.split("_")[0]
    date = selected_file.split("_")[1].split(".")[0]
    # Button to plot selected NDVI file
    if st.sidebar.button("Plot Vegetative Indices"):
        # st.write("Vegetative Indices Map")
        st.markdown("<h4 style='font-size: 24px;'>Vegetative Indices Map</h4>", unsafe_allow_html=True)
        if vegetative_index == 'EVI':
            st.write(f"Plotting EVI for {date}")
        else:
            st.write(f"Plotting NDVI for {date}")
        ndvi_map = utils.create_ndvi_map(ndvi_file_path)
        RasterPlot = ['Done']
        utils.folium_static(ndvi_map)
        
else:
    st.sidebar.write("Please download satellite data first to get the list of available vegitative indices during the season")
    RasterPlot = []

# Sidebar section for plotting vegetative indices
st.sidebar.header("Plot Timeseries Graph")

# Button to trigger satellite data download
coordsFarmA = eval(st.sidebar.text_input("Add Coordinates for Farm A in this form [75.716615,29.298713]", value=[75.716615,29.298713]))
coordsFarmB = eval(st.sidebar.text_input("Add Coordinates for Farm B in this form [75.716271,29.293248]", value=[75.724125,29.293174]))
if st.sidebar.button("Select Farm-A"):
    basic_map = utils.map_display(aoi)
    utils.folium.Marker(location=[coordsFarmA[1], coordsFarmA[0]], popup="Farm A").add_to(basic_map)
    utils.folium.Marker([coordsFarmA[1], coordsFarmA[0]], icon=utils.folium.DivIcon(html='<div style="font-size: 12pt; font-weight: bold; color: blue;">Farm A</div>')).add_to(basic_map)
    utils.folium.Marker(location=[coordsFarmB[1], coordsFarmB[0]], popup="Farm B").add_to(basic_map)
    utils.folium.Marker([coordsFarmB[1], coordsFarmB[0]], icon=utils.folium.DivIcon(html='<div style="font-size: 12pt; font-weight: bold; color: Orage;">Farm B</div>')).add_to(basic_map)
    utils.folium_static(basic_map)

# Button to trigger satellite data download
if st.sidebar.button("Select Farm-B"):
    basic_map = utils.map_display(aoi)
    utils.folium.Marker(location=[coordsFarmA[1], coordsFarmA[0]], popup="Farm A").add_to(basic_map)
    utils.folium.Marker([coordsFarmA[1], coordsFarmA[0]], icon=utils.folium.DivIcon(html='<div style="font-size: 12pt; font-weight: bold; color: blue;">Farm A</div>')).add_to(basic_map)
    utils.folium.Marker(location=[coordsFarmB[1], coordsFarmB[0]], popup="Farm B").add_to(basic_map)
    utils.folium.Marker([coordsFarmB[1], coordsFarmB[0]], icon=utils.folium.DivIcon(html='<div style="font-size: 12pt; font-weight: bold; color: Orage;">Farm B</div>')).add_to(basic_map)
    utils.folium_static(basic_map)

# Sidebar option to select Vegetation Index
veg_idx = st.sidebar.selectbox("Select Vegetation Index for Timeseries", ["NDVI", "EVI"])

# Button to plot selected NDVI file
if st.sidebar.button("Plot Timeseries Indices"):
    st.markdown("<h4 style='font-size: 24px;'>Vegetative indices based comparison of farms</h4>", unsafe_allow_html=True)
    farm_map = utils.create_farms_map(ndvi_file_path, coordsFarmA, coordsFarmB)
    utils.plt.figure(figsize=(10, 6))
    farmA = utils.extract_values_at_points(coordsFarmA)
    farmB = utils.extract_values_at_points(coordsFarmB)
    ndvi_df_A = farmA[farmA['index'] == veg_idx]
    ndvi_df_B = farmB[farmB['index'] == veg_idx]
    utils.plt.plot(ndvi_df_A['date'], ndvi_df_A['value'], marker='o', linestyle='-', label = 'Farm A')
    utils.plt.plot(ndvi_df_B['date'], ndvi_df_B['value'], marker='o', linestyle='-', label = 'Farm_B')
    utils.plt.title(f"{veg_idx} Time Series")
    utils.plt.xlabel("Date")
    utils.plt.ylabel(veg_idx)
    utils.plt.xticks(rotation=45)
    utils.plt.legend()
    st.pyplot(utils.plt)
    utils.folium_static(farm_map)
    RasterPlot = ['Done']