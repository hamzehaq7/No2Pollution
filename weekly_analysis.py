import xarray as xr
import cfgrib
import ecmwflibs
import netCDF4 as nc
import numpy as np
import os
from shapely.geometry import Point
from scipy.interpolate import griddata, interpn
import geopandas as gpd
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from geopy import distance
from matplotlib.cm import ScalarMappable

# -------------------------- Plot the variable vs. lat and long on common grid ----------------------------
def plot(common_long, common_lat, emissions, u_intp, v_intp, title):
    # Scatter plot for NO2 emissions with and without wind
    length = min(len(common_long), len(common_lat), len(emissions), len(u_intp), len(v_intp))

    geometry = [Point(lat, lon) for lat, lon in zip(common_lat[:length], common_long[:length])]
    gdf = gpd.GeoDataFrame(geometry=geometry, crs='EPSG:4326', data={'NO2_Level': emissions[:length]})

    mask = ~np.isnan(u_intp[:length]) & ~np.isnan(v_intp[:length])  
    # Plot the map with emission locations
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    fig, ax = plt.subplots(figsize=(12, 8))
    world.boundary.plot(ax=ax, linewidth=1)
    scatter = gdf[mask].plot(ax=ax, c='NO2_Level', cmap='viridis', markersize=50, label='NO2 Emission Locations')

    # Quiver plot for wind vector field
    ax.quiver(np.array(common_lat)[:length][mask], np.array(common_long)[:length][mask], np.array(u_intp)[:length][mask],
            np.array(v_intp)[:length][mask],
            scale=50, scale_units='inches', color='blue', width=0.005, label='Wind Vector Field')

    plt.title(title)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    sm = ScalarMappable(cmap="viridis", norm=plt.Normalize(vmin=0.0, vmax=0.001))
    sm.set_array(gdf[mask]['NO2_Level'])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('NO2 Level (mole/m^2)')

    ax.legend()
    plt.show()

# -------------------------- Plot the variable vs. lat and long on common grid ----------------------------
def plot2(common_long, common_lat, emissions, u_intp, v_intp, title, x1, y1, x2, y2):
    # Scatter plot for NO2 emissions with and without wind
    length = min(len(common_long), len(common_lat), len(emissions), len(u_intp), len(v_intp))

    geometry = [Point(lat, lon) for lat, lon in zip(common_lat[:length], common_long[:length])]
    gdf = gpd.GeoDataFrame(geometry=geometry, crs='EPSG:4326', data={'NO2_Level': emissions[:length]})

    mask = ~np.isnan(u_intp[:length]) & ~np.isnan(v_intp[:length])

    # Plot the map with all emission locations
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    fig, ax = plt.subplots(figsize=(12, 8))
    world.boundary.plot(ax=ax, linewidth=1)
    scatter = gdf[mask].plot(ax=ax, c='NO2_Level', cmap='viridis', markersize=50, label='NO2 Emission Locations')

    # Quiver plot for wind vector field
    ax.quiver(np.array(common_lat)[:length][mask], np.array(common_long)[:length][mask], np.array(u_intp)[:length][mask],
              np.array(v_intp)[:length][mask],
              scale=50, scale_units='inches', color='blue', width=0.005, label='Wind Vector Field')

    plt.title(title)
    plt.xlabel('Latitude')
    plt.ylabel('Longitude')

    sm = ScalarMappable(cmap="viridis", norm=plt.Normalize(vmin=0.0, vmax=0.001))
    sm.set_array(gdf[mask]['NO2_Level'])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('NO2 Level (mole/m^2)')

    ax.legend()
    plt.show()

    # Plot a more zoomed-in version around the specified coordinates
    fig, ax = plt.subplots(figsize=(8, 6))

    # Set the x-axis and y-axis limits around the specified coordinates with increased zoom
    ax.set_xlim([x1 - 0.1, x2 + 0.1])
    ax.set_ylim([y1 - 0.1, y2 + 0.1])

    world.boundary.plot(ax=ax, linewidth=1)
    scatter = gdf[mask].plot(ax=ax, c='NO2_Level', cmap='viridis', markersize=50, label='NO2 Emission Locations')

    # Quiver plot for wind vector field
    ax.quiver(np.array(common_lat)[:length][mask], np.array(common_long)[:length][mask], np.array(u_intp)[:length][mask],
              np.array(v_intp)[:length][mask],
              scale=50, scale_units='inches', color='blue', width=0.005, label='Wind Vector Field')

    plt.title('More Zoomed In: ' + title)
    plt.xlabel('Latitude')
    plt.ylabel('Longitude')

    sm = ScalarMappable(cmap="viridis", norm=plt.Normalize(vmin=0.0, vmax=0.001))
    sm.set_array(gdf[mask]['NO2_Level'])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('NO2 Level (mole/m^2)')

    ax.legend()
    plt.savefig(title + ".png")
    plt.show()

# --------------- Variables -----------------------------
no2_data, latitude_data, longitude_data, quality = [], [], [], []
path = '/Users/aboodhameed/Desktop/NO2emission/June 2023'
data_dict = []
gd_groups,data = [],[]

# ----------------- Put the files in a list --------------
for filename in os.listdir(path):
    file = os.path.join(path, filename)
    data_dict.append(file)

# ----------------- Put the values in a list -------------
for i, values in enumerate(data_dict):
    data.append(nc.Dataset(values, mode='r'))
    gd_groups.append(data[i].groups['PRODUCT'])

    no2_data.append(gd_groups[i].variables['nitrogendioxide_tropospheric_column'][:])
    latitude_data.append(gd_groups[i].variables['latitude'][:])
    longitude_data.append(gd_groups[i].variables['longitude'][:])
    quality.append(gd_groups[i].variables['qa_value'][:])


# ---------------- Flatten longitude ---------------------
flattened_longitude_list = []
for i in longitude_data:
    for j in i:
        flattened_longitude_list.append(j.flatten())
long_flat = np.concatenate(flattened_longitude_list)

# ---------------- Flatten latitude -----------------------
flattened_latitude_list = []
for i in latitude_data:
    for j in i:
        flattened_latitude_list.append(j.flatten())
lat_flat = np.concatenate(flattened_latitude_list)


# ---------------- Flatten NO2 ----------------------------
flattened_no2_list = []
for i in no2_data:
    for j in i:
        flattened_no2_list.append(j.flatten())
no2_flat = np.concatenate(flattened_no2_list)


# ---------------- Flatten quality -------------------------
flattened_quality_list = []
for i in quality:
    for j in i:
        flattened_quality_list.append(j.flatten())
quality_flat = np.concatenate(flattened_quality_list)

# ------------------------- Wind Values Extraction ------------------------------------
ds = xr.open_dataset('/Users/aboodhameed/Desktop/Hack4Climate/windy.grib', engine='cfgrib')

cams_lat = ds.latitude.values
cams_long = ds.longitude.values
u = ds.u.values
v = ds.v.values

flat_u = u.flatten()
flat_v = v.flatten()

# Filtering the data to get common positions
required_values, common_lat, common_long, common = [],[],[],[]

for i in range(len(lat_flat)):
    if ((lat_flat[i] >= 24.10) & (lat_flat[i] <= 25.46)):
        if ((long_flat[i] >= 54.75) & (long_flat[i] <= 56.11)):
            if (quality_flat[i] > 0.5):
                common.append(i)
                common_long.append(lat_flat[i])
                common_lat.append(long_flat[i])


for values in common:
    required_values.append(no2_flat[values])

required_values = no2_flat[common]

# ---------------------- Extracted Emissions from No Wind Vectors ----------------------------------
nitrogen = required_values 
emission_no_wind, emissions_with_constant_wind, emissions_with_wind_from_dataset = [],[],[]
L = 1.23
T = 14400
constant_wind = 5

for i in range(len(nitrogen)):
    emission = L*nitrogen[i]/T
    emission_no_wind.append(emission)

def uv_zero(common_lat, common_long, cams_lat, cams_long):
    u_zero = []
    v_zero = []

    for i in range(len(common_lat)):
        point = np.array([common_lat[i], common_long[i]])
        points = (np.array(cams_lat), np.array(cams_long))

        u_zero.append(0)
        v_zero.append(0)
    return u_zero, v_zero

u_zero, v_zero = uv_zero(common_lat, common_long, cams_lat, cams_long)

# ---------------------- Extracted Emissions From Constant Wind Vectors -----------------------------
for i in range(len(common_long)-2):
    distance.distance((common_lat[i],common_long[i]),(common_lat[i+1],common_long[i+1])).m

for i in range(len(nitrogen)-2):
    numerator = L*constant_wind*(nitrogen[i+1] - nitrogen[i])
    denominator_y = common_lat[i+1] - common_lat[i-1]
    denominator_x = common_long[i+1] - common_long[i-1]

    inverted_y = numerator/denominator_y
    inverted_x = numerator/denominator_x

    emissions_wind = emission + inverted_x + inverted_y
    emissions_with_constant_wind.append(emissions_wind) 

def uv_const(common_lat, common_long, cams_lat, cams_long):
    u_const = []
    v_const = []

    for i in range(len(common_lat)):
        point = np.array([common_lat[i], common_long[i]])
        points = (np.array(cams_lat), np.array(cams_long))

        u_const.append(5)
        v_const.append(5)
    return u_const, v_const

u_const, v_const = uv_const(common_lat, common_long, cams_lat, cams_long)
# ---------------------- Extracted Emissions From Wind Vectors ----------------------------------------
v_intp,u_intp = [],[]
for i in range(len(common_lat)):
    point = np.array([common_lat[i], common_long[i]])
    points = (np.array(cams_lat), np.array(cams_long))

    u_intp.append(interpn(points,u,point))
    v_intp.append(interpn(points,v,point))
       
u_intp = [item[0] for item in u_intp]
v_intp = [item[0] for item in v_intp]

for i in range(len(nitrogen)-2):    
    numerator_x = u_intp[i+1]*nitrogen[i+1] - u_intp[i]*nitrogen[i] 
    denominator_x = common_long[i+1] - common_long[i-1]
    inverted_x = numerator_x/denominator_x

    numerator_y = v_intp[i+1]*nitrogen[i+1] - v_intp[i]*nitrogen[i] 
    denominator_y = common_lat[i+1] - common_lat[i-1]
    inverted_y = numerator_y/denominator_y

    emissions_wind_new = emission + L*(inverted_x + inverted_y)
    emissions_with_wind_from_dataset.append(emissions_wind_new) 

# ---------------------- Plot The Emissions -----------------------------------
x1, y1, x2, y2 = 54.728, 24.085, 56.183, 25.476
plot2(common_long, common_lat, emissions_with_wind_from_dataset, u_intp, v_intp, 'NO2 Emissions and Wind Vector Field for June 2023', x1, y1, x2, y2)
plot2(common_long, common_lat, emissions_with_constant_wind, u_const, v_const, 'NO2 Emissions and Constant Wind for June 2023', x1, y1, x2, y2)
plot2(common_long, common_lat, emission_no_wind, u_zero, v_zero, 'NO2 Emissions and Zero Wind for June 2023', x1, y1, x2, y2)

# ---------------------- Print The Emissions -----------------------------------
print(f'''
NO2 Tropospheric Column (mol/m^2): {nitrogen}
Emissions With No Wind - (mol/s.m^2): {emission_no_wind}     
Emissions With Constant Wind of {constant_wind} m/s - (mol/s.m^2): {emissions_with_constant_wind}
Emissions With Wind from the Dataset (mol/s.m^2): {emissions_with_wind_from_dataset}''')

