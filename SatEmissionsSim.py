from math import pi, log10, cos, sin, sqrt
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata, interpn
from geopy import distance
import xarray as xr
import netCDF4 as nc
import ipywidgets as widgets
from IPython.display import display, HTML

# Define colors for styling
primary_color = '#3498db'  # Blue
background_color = '#ecf0f1'  # Light Gray

windpath = '/content/drive/MyDrive/HACK4CLIMATE/windy.grib'
copernicuspath = '/content/drive/MyDrive/HACK4CLIMATE/S5P_NRTI_L2__NO2____20231201T100739_20231201T101239_31780_03_020600_20231201T105021.nc'

# Styling for the input textboxes
input_style = {'description_width': 'initial'}

# Create textboxes for manual entering
min_lat_textbox = widgets.FloatText(description='Minimum Latitude (degrees):', value=24.10, style=input_style)
max_lat_textbox = widgets.FloatText(description='Maximum Latitude (degrees):', value=25.46, style=input_style)

min_long_textbox = widgets.FloatText(description='Minimum Longitude (degrees):', value=54.75, style=input_style)
max_long_textbox = widgets.FloatText(description='Maximum Longitude (degrees):', value=56.11, style=input_style)

# Output widget for header and calculation result
header_output = widgets.Output()
output_text = widgets.Output()

# -------------------------- Plot the variable vs. lat and long on a common grid ----------------------------
def plot(lat, long, variable, label, title, imagename):
    points = np.column_stack((lat, long))
    grid_x, grid_y = np.mgrid[min(lat):max(lat):100j, min(long):max(long):100j]
    grid_points = np.column_stack((grid_x.ravel(), grid_y.ravel()))

    z = griddata(points, variable, grid_points, method='linear')
    z = z.reshape(grid_x.shape)

    plt.pcolormesh(grid_x, grid_y, z, cmap='viridis', shading='auto')
    plt.colorbar(label=label)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title(title)
    # plt.savefig(imagename, format='eps', bbox_inches='tight')
    plt.grid(True, linestyle='--', linewidth=0.5, color='black', alpha=0.5)
    plt.show()

# Function to perform calculations and display results
def calculate_and_display(button):
    with output_text:
        output_text.clear_output(wait=True)
        display(get_calculated_output())
        # Plot 3D satellite position (Assumption: You need to define plot_satellite function)
        # plot_satellite(lat_textbox.value, long_textbox.value, height_textbox.value)


def get_calculated_output():
    min_lat = min_lat_textbox.value
    max_lat = max_lat_textbox.value
    min_longt = min_long_textbox.value
    max_longt = max_long_textbox.value

    # -------------------------- Extract the wind data from the grib file ----------------------------
    ds = xr.open_dataset(windpath, engine='cfgrib')

    # -------------------------- Extract the specific data  -------------------------------------
    cams_lat = ds.latitude.values
    cams_long = ds.longitude.values
    u = ds.u.values
    v = ds.v.values

    cams_lat = cams_lat.flatten()
    cams_long = cams_long.flatten()
    flat_u = u.flatten()
    flat_v = v.flatten()

    data = nc.Dataset(copernicuspath, mode='r')
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
      if ((flattened_latitude[i] >= min_lat) and (flattened_latitude[i] <= max_lat)):
        if ((flattened_longitude[i] >= min_longt) and (flattened_longitude[i] <= max_longt)):
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
        emission = L * nitrogen[i] / T
        emission_no_wind.append(emission)

    # Assuming dir is a valid directory path
    dir = '/content/drive/MyDrive/HACK4CLIMATE/'
    plot(common_lat, common_long, nitrogen, 'NO2 Tropospheric Column (mol/m^2)',
         'Interpolated NO2 Tropospheric Column on Grid', dir + 'no2_conc.eps')
    plot(common_lat, common_long, emission_no_wind, 'Emissions With No Wind (mol/s.m^2)',
         'Interpolated Emissions Without Wind', dir + 'em_no_wind.eps')

    for i in range(len(common_long)-2):
        distance.distance((common_lat[i], common_long[i]), (common_lat[i+1], common_long[i+1])).m

    for i in range(len(nitrogen)-2):
        numerator = L * constant_wind * (nitrogen[i+1] - nitrogen[i])
        denominator_y = common_lat[i+1] - common_lat[i-1]
        denominator_x = common_long[i+1] - common_long[i-1]

        inverted_y = numerator / denominator_y
        inverted_x = numerator / denominator_x

        emissions_wind = emission + inverted_x + inverted_y
        emissions_with_constant_wind.append(emissions_wind)

    new_lat = common_lat[:364]
    new_long = common_long[:364]

    plot(new_lat, new_long, emissions_with_constant_wind, f'Emissions with Wind Speed of {constant_wind} m/s in mol/s.m^2',
         f'Interpolated Emissions - Wind Speed = {constant_wind} m/s', dir + 'em_cnst_wind.eps')

    v_intp, u_intp = [], []
    for i in range(len(common_lat)):
        point = np.array([common_lat[i], common_long[i]])
        points = (np.array(cams_lat), np.array(cams_long))

        u_intp.append(interpn(points, u, point))
        v_intp.append(interpn(points, v, point))

    u_intp = [item[0] for item in u_intp]
    v_intp = [item[0] for item in v_intp]

    for i in range(len(nitrogen)-2):
        numerator_x = u_intp[i+1] * nitrogen[i+1] - u_intp[i] * nitrogen[i]
        denominator_x = common_long[i+1] - common_long[i-1]
        inverted_x = numerator_x / denominator_x

        numerator_y = v_intp[i+1] * nitrogen[i+1] - v_intp[i] * nitrogen[i]
        denominator_y = common_lat[i+1] - common_lat[i-1]
        inverted_y = numerator_y / denominator_y

        emissions_wind_new = emission + L * (inverted_x + inverted_y)
        emissions_with_wind_from_dataset.append(emissions_wind_new)

    plot(new_lat, new_long, emissions_with_wind_from_dataset, 'Emissions with Wind from the dataset in mol/s.m^2',
         'Interpolated Emissions w/wind Speed', dir + 'em_var_wind.eps')


# Create a button to trigger calculations
calculate_button = widgets.Button(description='Analyze', style=widgets.ButtonStyle(button_color=primary_color))

# Attach the function to the button click event
calculate_button.on_click(calculate_and_display)

# Header with website-like styling
header_html = HTML("""
    <div style="background-color: {}; padding: 20px; text-align: center; font-size: 24px; color: white; border-radius: 10px;">
        SatEmissionsSim
    </div>
""".format(primary_color))

# Display header and main layout
with header_output:
    display(header_html)

# Styling for the main layout
main_layout = widgets.VBox([
    widgets.HBox([min_lat_textbox, max_lat_textbox]),
    widgets.HBox([min_long_textbox, max_long_textbox]),
    calculate_button,
    output_text
])

# Apply styling to the main layout
main_layout.layout.display = 'flex'
main_layout.layout.flex_flow = 'column'
main_layout.layout.align_items = 'center'
main_layout.layout.margin = '20px'

# Display main layout
display(widgets.VBox([header_output, main_layout], layout=widgets.Layout(display='flex', flex_flow='column', align_items='center')))
