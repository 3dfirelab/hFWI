# hFWI - hourly Fire Weather Index Calculation

This code computes Fire Weather Index (FWI) Calculation as defined in [Von Wagner 1987](https://ostr-backend-prod.azurewebsites.net/server/api/core/bitstreams/5a865686-e097-40df-abc0-65f54c6ff379/content) for 2D fields from weather model outputs. 
It uses:
1. relative humidity at 2m
1. temperature at 2m
1. resolved wind speed at 10m
1. accumulated precipitation (mm) 
and outputs FWI as well as the intermediate indices, Fine Fuel Moistuce Code (FFMC), Duff Moisture Code (DMC), Droguht Code (DC), Initial Spread Index (ISI) and Buildingup Index (BUI).

The code can be run using hourly or daily time step.
`src/fwi.py` contains a test of the indices calculation against data provided in [Wang et al 2015](https://courses.seas.harvard.edu/climate/eli/Courses/global-change-debates/Sources/Forest-fires/aridity-indices/code-for-calculating-canadian-fire-weather-index.pdf)
