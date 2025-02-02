# hFWI - hourly Fire Weather Index Calculation

This code computes Fire Weather Index (FWI) Calculation as defined in [Von Wagner 1987](https://ostr-backend-prod.azurewebsites.net/server/api/core/bitstreams/5a865686-e097-40df-abc0-65f54c6ff379/content) for 2D fields from weather model outputs. 
It uses:
* relative humidity
* temperature 
* resolved wind speed 
* accumulated precipitation (mm)

The three first variable are extracted at model first vertical level at cell center.
The code is only working on MesoNH output for now.

It outputs FWI as well as the intermediate indices, Fine Fuel Moistuce Code (FFMC), Duff Moisture Code (DMC), Droguht Code (DC), Initial Spread Index (ISI) and Buildingup Index (BUI) in a dataarray format.

It adapts to the time steps avaiable using always -1h field variable to run the time integration. It can get as input daily (at 12h) time step or higher resolution.

In `src/fwi.py` there is a test of the indices calculation against data provided in [Wang et al 2015](https://courses.seas.harvard.edu/climate/eli/Courses/global-change-debates/Sources/Forest-fires/aridity-indices/code-for-calculating-canadian-fire-weather-index.pdf). `0.1` variation appears on some test case.


The 'mnh' branch integrates the hourly FFMC equation from Wagner et al (1977) using input data from the MesoNH model.
Below is an example:
![timeSerieFFMCFWI](/img/ffmcfwi_FCAST_000021.png)
daily DMC and DC are also calculated, and hourly FWI is output.

