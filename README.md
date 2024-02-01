# No2Pollution
This repository is linked to our conference paper titled "No2Pollution: A Guide to Analyze Air Quality via Level-2 Satellite Imagery"

Paper link: [to be added]

## Setup
To install the relevant libraries, you can run in your terminal:

```html
pip install netCDF4 numpy matplotlib scipy geopy cfgrib xarray ecmwflibs
```

## Usage
* <kbd style="background-color: #f0f0f0; padding: 5px; border-radius: 5px;">no2pollution.py</kbd> contains the code used to analyze NO2 concentrations for a whole week (June 1st-7th, 2023) as well as the calculation of emissions (1) without wind, (2) with constant wind, and (3) with variable wind (obtained from ADS)
* <kbd style="background-color: #f0f0f0; padding: 5px; border-radius: 5px;">SatEmissionsSim.py</kbd> contains the GUI that allows users to enter the desired set of coordinates for analysis. To use the simulator, you need to input 2 different datasets, referred to as variables in the code: _windpath_ and _copernicuspath_.

   _windpath_ is the path associated with the _.grib_ file of the dataset obtained from the Copernicus Climate Data Store: https://ads.atmosphere.copernicus.eu/cdsapp#!/dataset/cams-global-atmospheric-composition-forecasts?tab=form. There is _windy.grib_ file attached, might be of use..

   _copernicuspath_ is the path associated with the _.nc_ file of the dataset obtained from the Sentinel-5P Satellite using Copernicus Browser: https://dataspace.copernicus.eu/browser/?zoom=7&lat=45.83645&lng=10.74463&demSource3D=%22MAPZEN%22&cloudCoverage=30&dateMode=SINGLE. 

## Citation

If you find this project helpful or use it in your research, please consider citing it:

```html
@inproceedings{No2Pollution: A Guide to Analyze Air Quality via Level-2 Satellite Imagery,
  author = {Abdulrahman Hameed, Adil Mahroof, Mohammed Alnuaimi, Areej Alamin, and Hamzeh Abu Qamar},
  year = {2024},
  link = {https://github.com/hamzehaq7/No2Pollution}
}
```

## Acknowledgments
This work is a student effort. We were worked on the project at a hackathon and decided that the results were good enough to contribute & publish, wish us luck!
