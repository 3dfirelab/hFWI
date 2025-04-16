from datetime import datetime, timedelta
import numpy as np
import pygrib
import resource
from os import path,makedirs
from scipy.ndimage import binary_dilation
import glob 
import xarray as xr 
import pdb 
import f90nml
import os 
import pandas as pd

SPINUP_CUTOFF=3
RUN_REFRESH_PERIOD=24


class MNH:
    def __init__(self, filein):
        #constant 
        c_p = 1005.   # J  / kg / K
        R   =  287.05 #  J / kg / K
        kappa = R/c_p
        pressure_ref = 1.e5 # Pa
        
        #open file
        ds =  xr.open_dataset(filein)

        #wind at cell center
        ut = np.squeeze(ds.UT)
        utc = 0.5*(np.array(ut[:,:,:-1]) + np.array(ut[:,:,1:]))
        utc = utc[1:-1,1:-1,1:]

        vt = np.squeeze(ds.VT)
        vtc = 0.5*(np.array(vt[:,:-1,:]) + np.array(vt[:,1:,:]))
        vtc = vtc[1:-1,1:,1:-1]

        wt = np.squeeze(ds.WT)
        wtc = 0.5*(np.array(wt[:-1,:,:]) + np.array(wt[1:,:,:]))
        wtc = wtc[1:,1:-1,1:-1]

        self.rvap     = np.squeeze(ds.RVT)[1:-1,1:-1,1:-1]
        theta    = np.squeeze(ds.THT)[1:-1,1:-1,1:-1] 
        self.pressure = np.squeeze(ds.PABST)[1:-1,1:-1,1:-1]
        
        self.wind     = (np.sqrt( utc**2+vtc**2+wtc**2))*1.e-3*3600
        self.temp     = np.array(theta * (self.pressure/pressure_ref)**(kappa))
        self.rh       = np.array(self.relativeHumidity())
        self.rain    = np.squeeze(ds.ACPRR)[1:-1,1:-1]
        

    def relativeHumidity(self):
        y_v,pressure,temperature = self.rvap, self.pressure, self.temp
        '''
        in:
        yv: vapor mixing ratio kg/kg
        pressure Pa
        temperature K 
        
        out:relative humidity in %
        '''
        Wvap = 18.01
        Wair = 28.85
        # pressure at saturation form sonntag http://cires.colorado.edu/~voemel/vp.html
        p_vsat =   np.exp(-6096.9385/temperature + 16.635794  - 2.711193e-2 * temperature  
                          + 1.673952e-5 * temperature**2 + 2.433502 * np.log(temperature) ) * 100 
        x_vsat = p_vsat/pressure
        y_vsat = x_vsat*Wvap / (x_vsat*Wvap+ (1-x_vsat)*Wair) # kg/kg
        return y_v/y_vsat * 100.


################
def get_filename_for_date(package, date, dir):
    
    date = pd.to_datetime(date).to_pydatetime()
    run = date
    while run.hour % RUN_REFRESH_PERIOD != 0:
        run = run + timedelta(hours=-1)
    echeance = (date - run).total_seconds() / 3600
    if echeance < SPINUP_CUTOFF:
        run = run + timedelta(hours=-RUN_REFRESH_PERIOD)

    run_txt = run.strftime("%Y-%m-%d")

    date_txt_current = date.strftime("%Y-%m-%dT%H:%M:%SZ")

    date_previous = date + timedelta(hours=-1)
    date_txt_previous = date_previous.strftime("%Y-%m-%dT%H:%M:%SZ")

    current = f"{dir}/{run_txt}.{date.hour:02d}Z.{hour:02d}.{package}.grib2"
    previous = f"{dir}/{package}.grib"
   
    pdb.set_trace()
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
    

############
def load_dates_and_filenames(flag_model, dirin, iseg):
    
    if flag_model == 'mnh':
        
        exsegnam = f90nml.read(dirin+'/EXSEG1.nam')

        files = sorted(glob.glob(dirin+'/{:s}.{:d}*.nc'.format(exsegnam['NAM_CONF']['CEXP'],iseg)))
        files_out = []
        out = None
        for ifile, file, in enumerate(files):
            if '000.nc' in file: continue
            with xr.open_dataset(file) as ds: 
                time_ = ds.time[0].data
                if out is None: 
                    out = np.array([time_])
                else:
                    if out[-1]-time_ == 0: continue
                    out = np.append(out, time_)
            files_out.append(file)
   
        cexp = exsegnam['NAM_CONF']['CEXP']

    elif flag_model == 'arome':
        
        files = sorted(glob.glob(dirin+'/*{:s}*.grib2'.format('SP1')))
        day_str =  os.path.basename(files[0]).split('.')[0]
        hour_str =  os.path.basename(files[0]).split('.')[1][:-1]
   
        base_dt = datetime.strptime(day_str + hour_str, "%Y%m%d%H") 
        
        out = None
        files_out = []
        for ifile, file, in enumerate(files):
       
            step_str =  os.path.basename(file).split('.')[2][:-1]
            if int(step_str) == 0 : continue
            time_ = base_dt + timedelta(hours=int(step_str))
            
            if out is None: 
                out = np.array([time_])
            else:
                if out[-1]-time_ == 0: continue
                out = np.append(out, time_)
            files_out.append(file)

        cexp = None 

    else:   
        print('flag_model not defined yet: flag_model=',flag_model)
        print('stop here')
        sys.exit()
    
    return np.array(out, dtype='datetime64[ns]'),np.array(files_out), cexp


##############
def load_data_for_date(flag_model,date, filein, it):

    if flag_model == 'mnh': 
        
        #3D
        mnhdata = MNH(filein)
        
        #2D
        #wind        # km/h
        #humidity    # % (<100)
        #temperature # C
        #rain        # mm
        wind        = mnhdata.wind[0,:,:]
        humidity    = mnhdata.rh[0,:,:]
        temperature = mnhdata.temp[0,:,:] - 273.14
        rain        = mnhdata.rain[:,:]
        

    elif flag_model == 'arome': 

        filename_sp1 = filein
        filename_sp2 = filein.replace('SP1','SP2') 
        filename_sp2_previous = get_previous(filename_sp2, it)

        grib_sp1 = pygrib.open(filename_sp1)
        grib_sp2 = pygrib.open(filename_sp2)
        if filename_sp2_previous is not None:
            grib_sp2_previous = pygrib.open(filename_sp2_previous)

        wind_u = grib_sp1.select(shortName='10u')[0].values
        wind_v = grib_sp1.select(shortName='10v')[0].values

        wind = np.sqrt(wind_u**2 + wind_v**2) * 3.6

        humidity = grib_sp1.select(shortName='2r')[0].values # %
        temperature = grib_sp1.select(shortName='2t')[0].values -273.15 # celcius
        
        rain_acc_current = grib_sp2.select(shortName='tirf')[0].values # kg/m2 = mm/m2
        if filename_sp2_previous is not None:
            rain_acc_previous = grib_sp2_previous.select(shortName='tirf')[0].values
            rain = rain_acc_current - rain_acc_previous # acc sur 1h
        else: 
            rain = rain_acc_current 
        
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

    else: 
        print('flag_model not defined yet: flag_model=',flag_model)
        print('stop here')
        sys.exit()
    

    return wind, humidity, temperature, rain

def get_previous(name, iname):

    day_str =  os.path.basename(name).split('.')[0]
    hour_str = os.path.basename(name).split('.')[1][:-1]
    step_str = os.path.basename(name).split('.')[2][:-1]
    try: 
        date_base = datetime.strptime(day_str + '{:02d}'.format(int(hour_str)), "%Y%m%d%H") 
        datet = datetime.strptime(day_str + '{:02d}'.format(int(hour_str)), "%Y%m%d%H") + timedelta(hours=int(step_str))
    except: 
        pdb.set_trace()

    if  iname == 0: 
        return None
    
    else: 
        datet_prev = datet + timedelta(hours=-1)
        
        day_str_prev = day_str #datet_prev.strftime("%Y%m%d")
        try:
            step_str_prev = '{:02d}'.format( int ((datet_prev - date_base).total_seconds()/3600) )
        except: 
            pdb.set_trace()
        name_prev = name.replace(day_str, day_str_prev).replace(step_str+'H',step_str_prev+'H')
    
        if not(os.path.isfile(name_prev)):pdb.set_trace()

        return name_prev



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

def loadSeaLandMask(flag_model, dirin, files):

    if flag_model == 'mnh':
        
        file = files[0]
        with xr.open_dataset(file) as ds: 
            mask = ds.Z0SEA[1:-1,1:-1] 
       
            lat2d = ds.LAT[1:-1,1:-1]
            lon2d = ds.LON[1:-1,1:-1]

        mask = np.where(mask==999,1,0)

    elif flag_model == 'arome': 
    
        grib = pygrib.open('../dataStatic/CONSTANT_AROME_EURW1S100_2024.grb')
        mask = grib.select()[1].values # =1 for land
    
        lat2d, lon2d = grib.select()[1].latlons()
        

    else: 
        print('flag_model not defined yet: flag_model=',flag_model)
        print('stop here')
        sys.exit()
    '''
    grib = pygrib.open('./dataStatic/CONSTANT_AROME_EURW1S100_2024.grb')
    mask = grib.select()[1].values # =1 for land
    return binary_dilation(mask)
    '''

    return np.array(mask), np.array(lat2d), np.array(lon2d)


def getMeanLatitdue(flag_model, dirin, files):

    if flag_model == 'mnh':
        file = files[0]
        with xr.open_dataset(file) as ds: 
            meanlat = ds.LAT[1:-1,1:-1].mean() 
    
    elif flag_model == 'arome':
        file = files[0]
        grib = pygrib.open(file)
        meanlat = grib.select()[1].latlons()[0].mean()
        
    else: 
        print('flag_model not defined yet: flag_model=',flag_model)
        print('stop here')
        sys.exit()

    return float(meanlat) 

