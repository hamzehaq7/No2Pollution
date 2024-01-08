import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata, interpn
import math
from geopy import distance
import cfgrib
import xarray as xr
import ecmwflibs

# -------------------------- Plot the variable vs. lat and long on common grid ----------------------------
def plot(lat, long, variable, label, title, imagename):
        x, y = np.meshgrid(np.linspace(min(lat), max(lat), 100),
                        np.linspace(min(long), max(long), 100))
        z = griddata((lat, long), variable, (x, y), method='linear')
        plt.pcolormesh(x, y, z, cmap='viridis', shading='auto')
        plt.colorbar(label=label)
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        plt.title(title)
        plt.savefig(imagename, format='eps', bbox_inches='tight')
        plt.grid(True, linestyle='--', linewidth=0.5, color='black', alpha=0.5)
        plt.show()


# -------------------------- Extract the wind data from the grib file ----------------------------
dir = '/home/abuqamhs/cwd/HACK4CLIMATE/Figures/'
ds = xr.open_dataset('/home/abuqamhs/cwd/HACK4CLIMATE/windy.grib', engine='cfgrib')

# -------------------------- Extract the specific data  -------------------------------------
cams_lat = ds.latitude.values
cams_long = ds.longitude.values
u = ds.u.values
v = ds.v.values

cams_lat = cams_lat.flatten()
cams_long = cams_long.flatten()
flat_u = u.flatten()
flat_v = v.flatten()


data = nc.Dataset('/home/abuqamhs/cwd/HACK4CLIMATE/S5P_NRTI_L2__NO2____20231201T100739_20231201T101239_31780_03_020600_20231201T105021/S5P_NRTI_L2__NO2____20231201T100739_20231201T101239_31780_03_020600_20231201T105021.nc', mode='r')
gd_groups = data.groups['PRODUCT']

# The required variable in 3D arrays
no2_data = gd_groups.variables['nitrogendioxide_tropospheric_column'][:]
latitude_data = gd_groups.variables['latitude'][:]
longitude_data = gd_groups.variables['longitude'][:]
quality = gd_groups.variables['qa_value'][:]

# The flattened 3D arrays to 1D
flattened_longitude = longitude_data.flatten()
flattened_latitude = latitude_data.flatten()
flatten_no2 = no2_data.flatten()
flattened_quality = quality.flatten()

# Filtering the data to get common positions
common_lat, common_long, common = [],[],[]
for i in range(len(flattened_latitude)):
    if ((flattened_latitude[i] >= 24.10) & (flattened_latitude[i] <= 25.46)):
        if ((flattened_longitude[i] >= 54.75) & (flattened_longitude[i] <= 56.11)):
            if (flattened_quality[i] > 0.5):
                common.append(i)
                common_long.append(flattened_latitude[i])
                common_lat.append(flattened_longitude[i])


required_values = []
for values in common:
    required_values.append(flatten_no2[values])

nitrogen = required_values 
emission_no_wind, emissions_with_constant_wind, emissions_with_wind_from_dataset = [],[],[]
L = 1.23
T = 14400
constant_wind = 5

for i in range(len(nitrogen)):
    emission = L*nitrogen[i]/T
    emission_no_wind.append(emission)

plot(common_lat, common_long, nitrogen, 'NO2 Tropospheric Column (mol/m^2)','Interpolated NO2 Tropospheric Column on Grid', dir + 'no2_conc.eps')
plot(common_lat, common_long, emission_no_wind, 'Emissions With No Wind (mol/s.m^2)','Interpolated Emissions Without Wind', dir + 'em_no_wind.eps')


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


new_lat = common_lat[:364]
new_long = common_long[:364]

plot(new_lat, new_long, emissions_with_constant_wind, f'Emissions with Wind Speed of {constant_wind} m/s in mol/s.m^2', f'Interpolated Emissions - Wind Speed = {constant_wind} m/s', dir + 'em_cnst_wind.eps')

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


plot(new_lat, new_long, emissions_with_wind_from_dataset, 'Emissions with Wind from the dataset in mol/s.m^2', 'Interpolated Emissions w/wind Speed', dir + 'em_var_wind.eps')

print(f'''
NO2 Tropospheric Column (mol/m^2): {nitrogen}
Emissions With No Wind - (mol/s.m^2): {emission_no_wind}     
Emissions With Constant Wind of {constant_wind} m/s - (mol/s.m^2): {emissions_with_constant_wind}
Emissions With Wind from the Dataset (mol/s.m^2): {emissions_with_wind_from_dataset}''')

