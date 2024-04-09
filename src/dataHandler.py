from datetime import datetime, timedelta
import numpy as np
import pygrib
import resource
from os import path,makedirs
from scipy.ndimage import binary_dilation

SPINUP_CUTOFF=3
RUN_REFRESH_PERIOD=24

################
def get_filename_for_date(package, date):
    run = date
    while run.hour % RUN_REFRESH_PERIOD != 0:
        run = run + timedelta(hours=-1)
    echeance = (date - run).total_seconds() / 3600
    if echeance < SPINUP_CUTOFF:
        run = run + timedelta(hours=-RUN_REFRESH_PERIOD)

    run_txt = run.strftime("%Y-%m-%dT%H:%M:%SZ")

    date_txt_current = date.strftime("%Y-%m-%dT%H:%M:%SZ")

    date_previous = date + timedelta(hours=-1)
    date_txt_previous = date_previous.strftime("%Y-%m-%dT%H:%M:%SZ")

    current = f"/home/fire/data/arome_ramdisk/{run_txt}/{date_txt_current}/{package}.grib"
    previous = f"/home/fire/data/arome_ramdisk/{run_txt}/{date_txt_previous}/{package}.grib"
    
    return current, previous


################

def load_data_for_date_cached(date):
    date_txt = date.strftime("%Y-%m-%dT%H:%M:%SZ")
    cache_dir = f"/home/fire/cache/arome-{date_txt}"

    if path.exists(path.join(cache_dir, "wind.npy")):
        print("load from cache")
        wind = np.load(path.join(cache_dir, "wind.npy"))
        rh = np.load(path.join(cache_dir, "rh.npy"))
        temp = np.load(path.join(cache_dir, "temp.npy"))
        rain = np.load(path.join(cache_dir, "rain.npy"))
        mask = np.load(path.join(cache_dir, "mask.npy"))

        wind = np.ma.masked_array(wind, mask)
        rh = np.ma.masked_array(rh, mask)
        temp = np.ma.masked_array(temp, mask)
        rain = np.ma.masked_array(rain, mask)

        return wind, rh, temp, rain
    else:
        wind, rh, temp, rain = load_data_for_date(date)
    
        makedirs(cache_dir, exist_ok=True)
        
        np.save(path.join(cache_dir, "wind.npy"), wind.data)
        np.save(path.join(cache_dir, "rh.npy"), rh.data)
        np.save(path.join(cache_dir, "temp.npy"), temp.data)
        np.save(path.join(cache_dir, "rain.npy"), rain.data)
        np.save(path.join(cache_dir, "mask.npy"), temp.mask)
        return wind, rh, temp, rain

def prepare_data_for_date(date, patch_size):
    date_txt = date.strftime("%Y-%m-%dT%H:%M:%SZ")
    cache_dir = f"/home/fire/cache/arome-{date_txt}"

    wind, rh, temp, rain = load_data_for_date(date)
    
    makedirs(cache_dir, exist_ok=True)

    ny = wind.shape[0]

    for i in range(0, ny, patch_size):
        patch_file = path.join(cache_dir, f"{i}.npz")
        print(wind.data[i:i+patch_size,:].shape, np.average(wind.data[i:i+patch_size,:]))
        print(patch_file)
        np.savez(patch_file,
                 wind = wind.data[i:i+patch_size,:],
                 rh = rh.data[i:i+patch_size,:],
                 temp = temp.data[i:i+patch_size,:],
                 rain = rain.data[i:i+patch_size,:],
                 mask = temp.mask[i:i+patch_size,:]
        )

def load_patch_for_date(date, patch_index):
    date_txt = date.strftime("%Y-%m-%dT%H:%M:%SZ")
    cache_dir = f"/home/fire/cache/arome-{date_txt}"
    patch_file = path.join(cache_dir, f"{patch_index}.npz")

    data = np.load(patch_file)

    wind = np.ma.masked_array(data['wind'], data['mask'])
    rh = np.ma.masked_array(data['rh'], data['mask'])
    temp = np.ma.masked_array(data['temp'], data['mask'])
    rain = np.ma.masked_array(data['rain'], data['mask'])
    
    return wind, rh, temp, rain
    

def load_data_for_date(date):

    filename_sp1, ignore                = get_filename_for_date("SP1", date)
    filename_sp2, filename_sp2_previous = get_filename_for_date("SP2", date)

    grib_sp1 = pygrib.open(filename_sp1)
    grib_sp2 = pygrib.open(filename_sp2)
    grib_sp2_previous = pygrib.open(filename_sp2)

    wind_u = grib_sp1.select(shortName='10u')[0].values
    wind_v = grib_sp1.select(shortName='10v')[0].values

    wind = np.sqrt(wind_u**2 + wind_v**2) * 3.6

    rh = grib_sp1.select(shortName='2r')[0].values # %
    temp = grib_sp1.select(shortName='2t')[0].values -273.15 # celcius
    
    rain_acc_current = grib_sp2.select(shortName='tirf')[0].values # kg/m2 = mm/m2
    rain_acc_previous = grib_sp2_previous.select(shortName='tirf')[0].values

    rain = rain_acc_current - rain_acc_previous # acc sur 1h

    # sp1
    # 1:2 metre temperature:K (instant):regular_ll:heightAboveGround:level 2 m:fcst time 7 hrs:from 202404080000
    # 2:2 metre relative humidity:% (instant):regular_ll:heightAboveGround:level 2 m:fcst time 7 hrs:from 202404080000
    # 3:10 metre U wind component:m s**-1 (instant):regular_ll:heightAboveGround:level 10 m:fcst time 7 hrs:from 202404080000
    # 4:10 metre V wind component:m s**-1 (instant):regular_ll:heightAboveGround:level 10 m:fcst time 7 hrs:from 202404080000
    # 5:10 metre eastward wind gust since previous post-processing:m s**-1 (max):regular_ll:heightAboveGround:level 10 m:fcst time 6-7 hrs (max):from 202404080000
    # 6:10 metre northward wind gust since previous post-processing:m s**-1 (max):regular_ll:heightAboveGround:level 10 m:fcst time 6-7 hrs (max):from 202404080000
    # sp2
    # 1:Convective Available Potential Energy instantaneous:m**2 s**-2 (instant):regular_ll:unknown:levels 0-3000:fcst time 7 hrs:from 202404080000
    # 2:193:193 (instant):regular_ll:surface:level 0:fcst time 7 hrs:from 202404080000
    # 3:Graupel (snow pellets) precipitation rate:kg m-2 s-1 (accum):regular_ll:surface:level 0:fcst time 0-7 hrs (accum):from 202404080000
    # 4:Surface pressure:Pa (instant):regular_ll:surface:level 0:fcst time 7 hrs:from 202404080000
    # 5:Low cloud cover:% (instant):regular_ll:surface:level 0:fcst time 7 hrs:from 202404080000
    # 6:High cloud cover:% (instant):regular_ll:surface:level 0:fcst time 7 hrs:from 202404080000
    # 7:Medium cloud cover:% (instant):regular_ll:surface:level 0:fcst time 7 hrs:from 202404080000
    # 8:Time integral of rain flux:kg m**-2 (accum):regular_ll:surface:level 0:fcst time 0-7 hrs (accum):from 202404080000
    # 9:Snow precipitation rate:kg m**-2 s**-1 (accum):regular_ll:surface:level 0:fcst time 0-7 hrs (accum):from 202404080000


    return wind, rh, temp, rain

def load_data_for_range(start, stop):

    start_txt = start.strftime("%Y-%m-%dT%H:%M:%SZ")
    stop_txt = stop.strftime("%Y-%m-%dT%H:%M:%SZ")
    cache_dir = f"/home/fire/cache/arome-{start_txt}-{stop_txt}"
    if path.exists(path.join(cache_dir, "date.npy")):
        print("load from cache")
        date_array = np.load(path.join(cache_dir, "date.npy"), allow_pickle= True)
        wind_array = np.load(path.join(cache_dir, "wind.npy"))
        rh_array = np.load(path.join(cache_dir, "rh.npy"))
        temp_array = np.load(path.join(cache_dir, "temp.npy"))
        rain_array = np.load(path.join(cache_dir, "rain.npy"))
        mask = np.load(path.join(cache_dir, "mask.npy"))

        wind_array = np.ma.masked_array(wind_array, mask)
        rh_array = np.ma.masked_array(rh_array, mask)
        temp_array = np.ma.masked_array(temp_array, mask)
        rain_array = np.ma.masked_array(rain_array, mask)
    
        print("RAM: " + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)) + " GB")
        
        return date_array, wind_array, rh_array, temp_array, rain_array

    else:
        print("load from gribs")
        date = start
        date_list = []
        wind_list = []
        rh_list = []
        temp_list = []
        rain_list = []
        
        while date <= stop:
            print(date)
            wind, rh, temp, rain = load_data_for_date(date)
            date_list.append(date)
            wind_list.append(wind)
            rh_list.append(rh)
            temp_list.append(temp)
            rain_list.append(rain)
            date = date + timedelta(hours=1)
    
        date_array = np.array(date_list)
        wind_array = np.ma.masked_array(wind_list)
        rh_array = np.ma.masked_array(rh_list)
        temp_array = np.ma.masked_array(temp_list)
        rain_array = np.ma.masked_array(rain_list)
    
        print("RAM: " + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)) + " GB")

        makedirs(cache_dir, exist_ok=True)
        
        np.save(path.join(cache_dir, "date.npy"), date_array)
        np.save(path.join(cache_dir, "wind.npy"), wind_array.data)
        np.save(path.join(cache_dir, "rh.npy"), rh_array.data)
        np.save(path.join(cache_dir, "temp.npy"), temp_array.data)
        np.save(path.join(cache_dir, "rain.npy"), rain_array.data)
        np.save(path.join(cache_dir, "mask.npy"), temp_array.mask)
        
        return date_array, wind_array, rh_array, temp_array, rain_array

def loadSeaLandMask(path_index):
    grib = pygrib.open('./dataStatic/CONSTANT_AROME_EURW1S100_2024.grb')
    mask = grib.select()[1].values # =1 for land
    return binary_dilation(mask)