from config import settings
import ee
import geedim as gd
import geojson
import os
from datetime import datetime, timedelta
import rasterio
import numpy as np
import folium
import ee
import matplotlib.pyplot as plt
import pandas as pd
from folium.raster_layers import ImageOverlay
from folium import LayerControl
from streamlit_folium import folium_static

def Sentinel_2_data_download(
        start_date = None,
        end_date = None,
        cloud_percentage = None,
        out_dir = None,
        aoi = None
):
    """
    Downloads Sentinel-2 images from Google Earth Engine based on the specified area of interest, date range, 
    cloud percentage, and output directory. It applies a cloud mask to the images before downloading them.

    Parameters:
        start_date (str): The start date of the image collection (inclusive). The format should be 'YYYY-MM-DD'. Defaults to 6 months before the current date if None.
        end_date (str): The end date of the image collection (inclusive). The format should be 'YYYY-MM-DD'. Defaults to the current date if None.
        cloud_percentage (int): The maximum allowed cloud cover percentage for filtering the images. Images with cloud cover greater than this percentage will be excluded. Defaults to 30% if None.
        out_dir (str): The directory to save the downloaded images. Defaults to the current working directory with a sub-folder 'Sentinel2_data'.
        aoi (dict or list): The area of interest (AOI) for the image collection. It should be in GeoJSON format or a list of coordinates. Defaults to a predefined AOI if None.

    Returns:
        str: A message indicating the success of the data download process and the number of images downloaded, so that further analysis can be done.
    """
    try:
        # Initialize Google Earth Engine (GEE)
        ee.Initialize()
        print("Google Earth Engine initialized successfully.")
    except:
        print("GEE Authentication required for the first time")
        # Authenticate and initialize GEE
        ee.Authenticate()
        ee.Initialize()

    # If no output directory is specified, set the default to current working directory
    if out_dir is None:
        out_dir = os.getcwd()
        out_dir = os.path.join(out_dir, 'Sentinel2_data')

    out_dir_vi = os.path.join(out_dir, 'Output')
    
    # Clean up the 'Output' directory if it already exists by deleting all files
    if os.path.exists(out_dir_vi):
        files = os.listdir(out_dir_vi)
        for file_name in files:
            file_path = os.path.join(out_dir_vi, file_name)
            if os.path.isfile(file_path):
                os.remove(file_path)

    # Clean up the main output directory (Sentinel2_data) if it exists
    if os.path.exists(out_dir):
        files = os.listdir(out_dir)
        for file_name in files:
            file_path = os.path.join(out_dir, file_name)
            if os.path.isfile(file_path):
                os.remove(file_path)

    # If the output directory doesn't exist, create it
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Default AOI coordinates if not provided by the user
    if aoi is None:
        aoi = [[
            [75.71107959499119,29.299710037552416],
            [75.71107959499119,29.282741818313212],
            [75.73246339195231,29.282741818313212],
            [75.73246339195231,29.299710037552416],
            [75.71107959499119,29.299710037552416]
            ]]
        coordinates = aoi # Use default AOI coordinates

    else:
        # Extract coordinates from the provided AOI
        coordinates = aoi["geometry"]["coordinates"]

    # Define the region as a polygon geometry in GEE
    region = ee.Geometry.Polygon(coordinates)

    # Default cloud percentage if not specified by the user
    if cloud_percentage is None:
        cloud_percentage = 30

    # Set default start and end date for image collection if not provided
    if start_date is None:
        current_date = datetime.now().date()
        end_date = current_date.strftime('%Y-%m-%d')
        start_date = (current_date - timedelta(days=6*30)).strftime('%Y-%m-%d')

    # Load Sentinel-2 Surface Reflectance and Cloud Probability Image Collections
    s2Sr = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    s2Clouds = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')

    # Convert start and end date to GEE Date objects
    START_DATE = ee.Date(start_date)
    END_DATE = ee.Date(end_date)
    MAX_CLOUD_PROBABILITY = 65

    # Function to mask clouds based on the cloud probability
    def maskClouds(img):
        clouds = ee.Image(img.get('cloud_mask')).select('probability')
        isNotCloud = clouds.lt(MAX_CLOUD_PROBABILITY)
        return img.updateMask(isNotCloud).divide(10000)

    # Filter the Sentinel-2 Surface Reflectance collection based on the region, date range, and cloud percentage
    s2Sr = s2Sr.filterBounds(region).filterDate(start_date, end_date).filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',cloud_percentage))
    s2Clouds = s2Clouds.filterBounds(region).filterDate(start_date, end_date)

   # Join S2 SR with cloud probability dataset to add cloud mask
    joinCondition = ee.Filter.equals(leftField='system:index', rightField='system:index')
    s2SrWithCloudMask = ee.Join.saveFirst('cloud_mask').apply(
        primary=s2Sr,
        secondary=s2Clouds,
        condition=joinCondition
    )

    # Apply cloud masking function to the Sentinel-2 SR collection and select required bands
    collection = ee.ImageCollection(s2SrWithCloudMask).map(maskClouds).select(['B2', 'B3', 'B4', 'B8'])

    try:
        count = int(collection.size().getInfo())
        print(f"Total number of images: {count}\n")

        # Download each image in the collection
        for i in range(0, count):
            image = ee.Image(collection.toList(count).get(i))
            name = image.get("system:index").getInfo() + ".tif"
            filename = os.path.join(os.path.abspath(out_dir), name)
            print(f"Downloading {i + 1}/{count}: {name}")
            img = gd.download.BaseImage(image)
            img.download(filename, region = region, crs = 'EPSG:4326')
        print("All the data has been downloaded successfully")
        return f"Total {count} tiles download successfully"

    except Exception as e:
        raise Exception(f"Error in downloading image: {e}")

def calculate_ndvi(red_band, nir_band):
    """Calculate NDVI."""
    denominator = nir_band + red_band
    
    # Handle division by zero
    with np.errstate(divide='ignore', invalid='ignore'):
        ndvi = np.where(denominator == 0, np.nan, (nir_band - red_band) / denominator)
    
    return ndvi

def calculate_evi(red_band, nir_band, blue_band):
    """Calculate EVI."""
    G = 2.5
    C1 = 6
    C2 = 7.5
    L = 1

    denominator = (nir_band + C1 * red_band - C2 * blue_band + L)
    with np.errstate(divide='ignore', invalid='ignore'):
        evi = G * np.where(denominator == 0, np.nan, (nir_band - red_band) / denominator)
    return evi  


def calculate_evi_ndvi():
    try:
        # Directory containing stacked raster data
        out_dir = os.getcwd()
        stacked_raster_dir = os.path.join(out_dir, 'Sentinel2_data')

        output_dir = os.path.join(stacked_raster_dir, 'Output')
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Iterate over each raster file in the directory
        # List to store NDVI arrays and corresponding dates
        for filename in os.listdir(stacked_raster_dir):
            if filename.endswith('.tif'):
                raster_path = os.path.join(stacked_raster_dir, filename)
                
                # Read bands from the raster stack
                with rasterio.open(raster_path) as src:
                    red_band = src.read(3).astype(np.float32)  # where band 4 is red
                    nir_band = src.read(4).astype(np.float32)  # where band 8 is NIR
                    blue_band = src.read(1).astype(np.float32)  # where band 2 is blue
                
                    # Calculate NDVI and EVI
                    ndvi = calculate_ndvi(red_band, nir_band)
                    evi = calculate_evi(red_band, nir_band, blue_band)
                    
                    # Update metadata for output
                    profile = src.profile
                    profile.update(dtype=rasterio.float32, count=1)
                    
                    # Write NDVI output
                    strpath = os.path.splitext(filename)[0].split("_")[0][:8]
                    datepath = datetime.strptime(strpath, "%Y%m%d").strftime("%Y-%m-%d")
                    ndvi_output_path = os.path.join(output_dir, f'NDVI_{datepath}.tif')
                    with rasterio.open(ndvi_output_path, 'w', **profile) as dst:
                        dst.write(ndvi, 1)
                    
                    # Write EVI output
                    evi_output_path = os.path.join(output_dir, f'EVI_{datepath}.tif')
                    with rasterio.open(evi_output_path, 'w', **profile) as dst:
                        dst.write(evi, 1)
        print("Indices for all the images has been completed successfully")
    except Exception as e:
        raise Exception(f"Error in calculating indices: {e}")


# Function to extract values at points
def extract_values_at_points(point):
    # Directory containing stacked raster data
    out_dir = os.getcwd()
    stacked_raster_dir = os.path.join(out_dir, 'Sentinel2_data/Output')

    if not os.path.exists(stacked_raster_dir):
        os.makedirs(stacked_raster_dir)
    # Example filename
    values = {"index":[], "date":[], "value":[]}
    for filename in os.listdir(stacked_raster_dir):
        if filename.endswith('.tif'):
            raster_file = os.path.join(stacked_raster_dir, filename)
            # Extract date from the filename
            basename = os.path.splitext(filename)[0]  # Remove the extension
            all_str = basename.split("_")
            VI = all_str[0]
            date_str = all_str[1][:10]  # Extract the date part (YYYYMMDD)

            # Convert date string to datetime object
            date = datetime.strptime(date_str, "%Y-%m-%d")
            # Open the stacked raster file
        
            with rasterio.open(raster_file) as src:
                # Read raster data into a numpy array
                raster_data = src.read(1)

                # Extract raster values at the specified points for each band
                row, col = src.index(point[0], point[1])  # Convert coordinates to row, column indices
                value = raster_data[row, col]
                values['index'].append(VI)
                values['date'].append(date)
                values['value'].append(value)

    # Convert data to pandas DataFrame
    df = pd.DataFrame(values)

    # Convert date column to datetime
    df['date'] = pd.to_datetime(df['date'])

    # Sort DataFrame by date
    df = df.sort_values(by='date')

    # Filter out NaN values
    df = df.dropna(subset=['value'])

    return df


def plot_indices(df):

    # Create a plot for each index (evi, ndvi)
    for index, group in df.groupby('index'):
        plt.plot(group['date'], group['value'], label=index)

    # Add labels and legend
    plt.xlabel('Date')
    plt.ylabel('Value')
    plt.title('Graph showing EVI and NDVI values')
    plt.legend()

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45)

    # Show plot
    plt.tight_layout()
    plt.show()


def ndvi_evi_raster_plot(ndvi_raster_path, evi_raster_path):
    # Add NDVI raster overlay
    with rasterio.open(ndvi_raster_path) as src:
        ndvi_img = src.read(1)
        ndvi_img = ndvi_img.astype(float)
        ndvi_min = np.nanmin(ndvi_img)
        ndvi_max = np.nanmax(ndvi_img)
        ndvi_rescaled = ((ndvi_img - ndvi_min) / (ndvi_max - ndvi_min))
        # ndvi_img = ndvi_img / 10000.0  # Convert to floating point representation
        ndvi_img_overlay = ImageOverlay(
            image=ndvi_rescaled,
            bounds=[[src.bounds.bottom, src.bounds.left], [src.bounds.top, src.bounds.right]],
            colormap=lambda x: (0, int(x * 255), 0),  # Green colormap for NDVI
            opacity=1,
            name="NDVI"
        )
    extent = src.bounds  # Get the spatial extent of the raster
    
    center_lat = (extent.bottom + extent.top) / 2.0
    center_lon = (extent.left + extent.right) / 2.0

    # Creating Folium Map centered at the center of the image
    m = folium.Map(location=[center_lat, center_lon], zoom_start=15, control_scale=True)

    ndvi_img_overlay.add_to(m)

    # Add EVI raster overlay (similar to NDVI)
    # Add NDVI raster overlay
    with rasterio.open(evi_raster_path) as src:
        evi_img = src.read(1)
        evi_img = evi_img.astype(float)
        evi_min = np.nanmin(evi_img)
        evi_max = np.nanmax(evi_img)
        evi_rescaled = ((evi_img - evi_min) / (evi_max - evi_min))
        evi_img = evi_img / 10000.0  # Convert to floating point representation
        evi_img_overlay = ImageOverlay(
            image=evi_rescaled,
            bounds=[[src.bounds.bottom, src.bounds.left], [src.bounds.top, src.bounds.right]],
            colormap=lambda x: (0, int(x * 255), 0),  # Green colormap for NDVI
            opacity=1,
            name="EVI"
        )
        evi_img_overlay.add_to(m)
    # Add RGB raster overlay (similar to NDVI)

    # Add Google Satellite basemap
    folium.TileLayer(
        tiles='https://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}',
        attr='Google',
        name='Google Satellite',
        overlay=False
    ).add_to(m)

    # Add Layer Control
    layer_control = LayerControl()
    layer_control.add_to(m)
    
    # Display the map
    return folium_static(m)

# Function to rescale NDVI values to [0, 255]
def rescale_ndvi(ndvi):
    ndvi_min = -1
    ndvi_max = 1
    ndvi_rescaled = ((ndvi - ndvi_min) / (ndvi_max - ndvi_min)) * 255
    return ndvi_rescaled

# Function to create folium map with NDVI overlay
def create_ndvi_map(ndvi_file):
    # Path to your NDVI raster file
    # Open the raster file
    with rasterio.open(ndvi_file) as src:
        ndvi = src.read(1, masked=False)  # Read the first band (NDVI) and mask nodata values
        extent = src.bounds  # Get the spatial extent of the raster
    
    ndvi_rescaled = rescale_ndvi(ndvi)

    center_lat = (extent.bottom + extent.top) / 2.0
    center_lon = (extent.left + extent.right) / 2.0

    # Creating Folium Map centered at the center of the image
    m = folium.Map(location=[center_lat, center_lon], zoom_start=15, control_scale=True)
    
    # Adding Google Satellite base map
    folium.TileLayer(
        tiles='https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}',
        attr='Google Satellite with Labels',
        name='Google Satellite with Labels',
        overlay=True
    ).add_to(m)

    # Define colormap function
    def colormap(x):
        if np.isnan(x):
            return (0, 0, 0)  # Assign black color for NaN values
        else:
            return (0, int(x), 0)  # Green colormap for NDVI
    # Creating image overlay of NDVI
    ndvi_img_overlay = folium.raster_layers.ImageOverlay(
        image=ndvi_rescaled,
        bounds=[[src.bounds.bottom, src.bounds.left], [src.bounds.top, src.bounds.right]],
        opacity=1,
        colormap=colormap # Green colormap for NDVI
    )
    ndvi_img_overlay.add_to(m)
    draw = folium.plugins.Draw()
    draw.add_to(m)

    return m

# Function to create folium map with NDVI overlay
def create_farms_map(ndvi_file, coordsFarmA, coordsFarmB):
    # Path to your NDVI raster file
    # Open the raster file
    with rasterio.open(ndvi_file) as src:
        ndvi = src.read(1, masked=False)  # Read the first band (NDVI) and mask nodata values
        extent = src.bounds  # Get the spatial extent of the raster

    center_lat = (extent.bottom + extent.top) / 2.0
    center_lon = (extent.left + extent.right) / 2.0

    # Creating Folium Map centered at the center of the image
    m = folium.Map(location=[center_lat, center_lon], zoom_start=15, control_scale=True)
    
    # Adding Google Satellite base map
    folium.TileLayer(
        tiles='https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}',
        attr='Google',
        name='Google Satellite with Labels',
        overlay=True
    ).add_to(m)

    folium.Marker(location=[coordsFarmA[1], coordsFarmA[0]], popup="Farm A").add_to(m)
    folium.Marker(location=[coordsFarmB[1], coordsFarmB[0]], popup="Farm B").add_to(m)

    folium.Marker([coordsFarmA[1], coordsFarmA[0]], icon=folium.DivIcon(html='<div style="font-size: 12pt; font-weight: bold; color: blue;">Farm A</div>')).add_to(m)
    folium.Marker([coordsFarmB[1], coordsFarmB[0]], icon=folium.DivIcon(html='<div style="font-size: 12pt; font-weight: bold; color: orange;">Farm B</div>')).add_to(m)
    draw = folium.plugins.Draw()
    draw.add_to(m)
    return m

# Function to create folium map with NDVI overlay
def map_display(aoi = None):
    # Creating Folium Map centered at the center of the image
    if aoi is None:
        m = folium.Map(location=[29.298713, 75.716615], zoom_start=14, control_scale=True)
    else:
        extent = aoi["geometry"]["coordinates"]
        center_lat = (extent[0][0][1] + extent[0][1][1]) / 2.0
        center_lon = (extent[0][0][0] + extent[0][2][0]) / 2.0

        # Creating Folium Map centered at the center of the image
        m = folium.Map(location=[center_lat, center_lon], zoom_start=15, control_scale=True)
        extent_poly = extent[0]
        extent_poly = []
        for lat, lon in extent[0]:
            extent_poly.append([lon, lat])

        # Create a Polygon object using the extent coordinates
        extent_polygon = folium.Polygon(locations=extent_poly, color='blue', fill=True, fill_color='blue', fill_opacity=0.2)
        # Add the polygon to the map
        extent_polygon.add_to(m)
    # Adding Google Satellite base map
    folium.TileLayer(
        tiles='https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}',
        attr='Google',
        name='Google Satellite with Labels',
        overlay=True
    ).add_to(m)
    
    draw = folium.plugins.Draw()
    draw.add_to(m)
    return m