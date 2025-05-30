import pygrib
import numpy as np
import matplotlib.pyplot as plt
import math
import pdb 
import tools
import datetime 
import glob 
from datetime import datetime, timedelta
import pandas as pd
import sys 
import xarray as xr 
import rioxarray
import importlib 
import os 
import pickle 
import tracemalloc
import psutil
import gc 
import time
from pathlib import Path


#homebrewed
import dataHandler
import plot_fwi_arome

def get_Ed(H,T):
    return 0.942*(H**0.679) + 11*(np.exp((H - 100)/10)) + 0.18*(21.1 - T)*(1 - np.exp(-0.115*H))
    
def get_Ew(H,T):
    return 0.618*(H**0.753) + 10*(np.exp((H - 100)/10)) + 0.18*(21.1 - T)*(1 - np.exp(-0.115*H))

class FWICLASS:

    ############
    def __init__(self,temp,rhum,wind,prcp,mask):
        self.updateMeteo(temp,rhum,wind,prcp)
        if isinstance(wind, np.ndarray):
            self.mask = np.where((mask==1) & (~wind.mask), 1, 0  )
        else:
            self.mask = mask

    def updateMeteo(self,temp,rhum,wind,prcp):
        self.h = rhum
        self.t = temp
        self.w = wind
        self.p = prcp

    
    ############
    def hFFMCcalc(self,ffmc0, mode='hourly'):
        
        fmc0 = (147.2*(101.0 - ffmc0))/(59.5 + ffmc0) #*Eq. 1*#

        '''
        H HU2M_DD2_12
        W W10M_DD2_12 
        T T2M_DD2_12
        '''
       
        maskLand = self.mask
        idx = np.where(maskLand==1)
        
        m0 = fmc0[idx]
        r0 = self.p[idx]
        W  = self.w[idx]
        H  = np.where(self.h[idx]>100,100.,self.h[idx])
        T  = self.t[idx]

        #integration of rain accumulation from past 24h in initial moisture m0
        idx_ = np.where(r0>0.5)
        rf = r0[idx_] - 0.5
        #try:
        mr = np.where(m0[idx_] <= 150, 
                  m0[idx_] + 42.5*rf*(np.exp(-100/(251 - m0[idx_])))*(1 - np.exp(-6.93/rf))                                   ,
                  m0[idx_] + 42.5*rf*(np.exp(-100/(251 - m0[idx_])))*(1 - np.exp(-6.93/rf)) + 0.0015*((m0[idx_] - 150)**2)*(rf**0.5))
        mr = np.where( mr > 250, 250, mr)
        m0[idx_] = mr
        #except:   
        #    pdb.set_trace()

        # compute Equilibrium Moisture Content (EMC) for desportion (Ed) and Absportion (Ew)
        Ed = get_Ed(H,T) #0.942*(H**0.679) + 11*(np.exp((H - 100)/10)) + 0.18*(21.1 - T)*(1 - np.exp(-0.115*H))
        Ew = get_Ew(H,T) #0.618*(H**0.753) + 10*(np.exp((H - 100)/10)) + 0.18*(21.1 - T)*(1 - np.exp(-0.115*H))

        #init
        m = np.zeros_like(m0) - 999
        #try:
        kbeta = lambda nu,W,T: (0.424*(1 - nu**1.7) + 0.0694*np.sqrt(W)*(1 - nu**8) ) *0.579*np.exp(0.0365*T)
        kalpha = lambda nu,W,T: (0.424*(1 - nu**1.7) + 0.0694*np.sqrt(W)*(1 - nu**8) ) *0.581*np.exp(0.0365*T)
        if mode == 'hourly':
            kk =  kbeta
        elif mode == 'daily':
            kk =  kalpha            

        #desportion
        idx_ = np.where( m0 > Ed)
        nu = H/100
        #nu = np.where(H>0,H/100,0)
        m[idx_] = Ed[idx_] + (m0[idx_] - Ed[idx_])*(10**( -1.*kk(nu[idx_],W[idx_],T[idx_]) ))

        #absportion
        idx_ = np.where( m0 < Ew)
        #nu = np.where(H<100,(100-H)/100,0)
        nu = (100-H)/100
        if nu.min() < 0: pdb.set_trace()
        m[idx_] = Ew[idx_] - (Ew[idx_] - m0[idx_])*(10**( -1.*kk(nu[idx_],W[idx_],T[idx_]) ))
        #except:   
        #    pdb.set_trace()

        #equilibrium
        idx_ = np.where( (Ed >= m0) & (m0 >= Ew) )
        m[idx_] = m0[idx_] 

        ffmc1 = np.zeros_like(ffmc0)
        ffmc1[idx] = (59.5 * (250.0 -m)) / (147.2 + m)   
        
        ffmc1 = np.where(ffmc1>101.0,101.0, ffmc1)
        ffmc1 = np.where(ffmc1<=0.0,0.0, ffmc1)
        
        return ffmc1


    ############
    def DMCcalc(self,dmc0,day,latitude):
        #El = [6.5,7.5,9.0,12.8,13.9,13.9,12.4,10.9,9.4,8.0,7.0,6.0]
        El = tools.DayLength_(latitude, day.timetuple().tm_yday, day.month)

        maskLand = self.mask
        idx = np.where(maskLand==1)
        
        t = self.t[idx]
        p = self.p[idx]
        t = np.where(t < -1.1,-1.1,t)
        rk = 1.894*(t+1.1) * (100.0-self.h[idx]) * (El*0.0001)     #*Eqs. 16 and 17*#
        pr = np.zeros_like(rk)
        
        idx_ = np.where(p>1.5)
        if len(idx_[0])>0:
            ra    = p[idx_]
            dmc0_ = dmc0[idx][idx_]
            rw = 0.92*ra - 1.27                                        #*Eq. 11*#
            wmi = 20.0 + 280.0/np.exp(0.023*dmc0_)                     #*Eq. 12*#
            
            b = np.zeros_like(dmc0_)
            idx__ = np.where(dmc0_<= 33.0)
            b[idx__] = 100.0 /(0.5 + 0.3*dmc0_[idx__])                   #*Eq. 13a*#
            
            idx__ = np.where((dmc0_ > 33.0) & (dmc0_ <= 65.0))
            b[idx__] = 14.0 - 1.3*np.log(dmc0_[idx__])                    #*Eq. 13b*#
            idx__ = np.where((dmc0_ > 33.0) & (dmc0_ > 65.0))
            b[idx__] = 6.2 * np.log(dmc0_[idx__]) - 17.2                  #*Eq. 13c*#
                     
            wmr = wmi + (1000*rw) / (48.77+b*rw)                        #*Eq. 14*#
            pr[idx_] = 43.43 * (5.6348 - np.log(wmr-20.0))               #*Eq. 15*#
        
        idx_ = np.where(p<=1.5)
        if len(idx_[0])>0:
            pr[idx_] = dmc0[idx][idx_]
        
        pr = np.where(pr<0.0, 0.0, pr)
        
        dmc_ = pr + rk
        
        dmc_ = np.where(dmc_<=1.0, 1.0, dmc_)

        dmc = np.zeros_like(dmc0)
        dmc[idx] = dmc_
        
        return dmc


    ############
    def DCcalc(self,dc0,day,latitude):

        maskLand = self.mask
        idx = np.where(maskLand==1)

        p = self.p[idx]
        t = self.t[idx]
        dc0_ = dc0[idx]

        idx_ = np.where(p > 2.8)
        if len(idx_[0])>0:
        
            rd = 0.83 * p[idx_] - 1.27                  #*Eq. 18*#
            Qo = 800.0 * np.exp(-dc0_[idx_] / 400.0)          #*Eq. 19*#
            Qr = Qo + 3.937 * rd                            #*Eq. 20*#
            Dr = 400.0 * np.log(800.0 / Qr)                 #*Eq. 21*# 

            dc_prev_ = np.copy(dc0_)
            dc_prev_[idx_] = np.where(Dr>0, Dr, 0.0)

        else: 
            dc_prev_ = dc0_

        Lf = tools.DryingFactor_(latitude, day.month)

        V = np.where(t >= -2.8, 0.36 * (t+2.8) + Lf, Lf)
        V = np.where(V<0, 0.0, V)
        
        dc_ = dc_prev_ + 0.5 * V

        dc = np.zeros_like(dc0)
        dc[idx] = dc_
        
        return dc

    
    ############
    def ISIcalc(self,ffmc):
                
        maskLand = self.mask
        idx = np.where(maskLand==1)

        ffmc_ = ffmc[idx]
        w_    = self.w[idx]
        
        mo = 147.2*(101.0-ffmc_) / (59.5+ffmc_)                          #*Eq. 1*#
        ff = 19.115*np.exp(mo*-0.1386) * (1.0+(mo**5.31)/49300000.0)     #*Eq. 25*#
        isi_ = ff * np.exp(0.05039*w_)                                    #*Eq. 26*#
        
        isi = np.zeros_like(ffmc)
        isi[idx] = isi_
        
        return isi


    ############
    def BUIcalc(self,dmc,dc):
        maskLand = self.mask
        idx = np.where(maskLand==1)

        dmc_ = dmc[idx]
        dc_  = dc[idx]

        bui_ = np.where(dmc_<=0.4*dc_, 
                       (0.8*dc_*dmc_) / (dmc_+0.4*dc_),                             #*Eq. 27a*#
                       dmc_-(1.0-0.8*dc_/(dmc_+0.4*dc_))*(0.92+(0.0114*dmc_)**1.7)  #*Eq. 27b*#
                       )
        
        bui_ = np.where(bui_<0, 0, bui_)

        bui = np.zeros_like(dc)
        bui[idx] = bui_
        
        return bui

    
    ############
    def FWIcalc(self,isi,bui):
        maskLand = self.mask
        idx = np.where(maskLand==1)

        isi_ = isi[idx]
        bui_  = bui[idx]
        
        bb = np.where(bui_ <= 80.0, 
                      0.1 * isi_ * (0.626*bui_**0.809 + 2.0),                       #*Eq. 28a*#
                      0.1 * isi_ * (1000.0/(25. + 108.64/np.exp(0.023*bui_)))       #*Eq. 28b*#
                     )
        
        fwi_ = np.copy(bb)                                                            #*Eq. 30b*#
        idx_ = np.where(bb > 1.0) 
        fwi_[idx_] = np.exp(2.72 * (0.434*np.log(bb[idx_]))**0.647)                     #*Eq. 30a*#

        fwi = np.zeros_like(isi)
        fwi[idx] = fwi_
        
        return fwi
    
    def self_destruct(self):
        #print(f"{self.name} is destroying itself...")
        # Optional: Clean up attributes
        for attr in list(self.__dict__):
            setattr(self, attr, None)
        # Optional: Notify the GC
        gc.collect()


class Indices:

    ############
    def __init__(self, shape, latitude, patch_index):

        self.latitude = 46.45 # valeur harcodée pour arome
        self.seaLandMask = np.ones(shape) #dataHandler.loadSeaLandMask(patch_index)
        
        self.ffmc = np.zeros(shape) + 85.0 # valeurs prises dans la func initialisation plus bas. merci de valider 
        self.dmc = np.zeros(shape) + 6.0
        self.dc = np.zeros(shape) + 15.0
        self.isi = np.zeros(shape)
        self.bui = np.zeros(shape)
        self.fwi = np.zeros(shape)

        self.hour = 0
        self.rain_buffer = np.zeros((24, shape[0], shape[1]))

        self.needSpinup = True

        self.fwisystem = FWICLASS(0, 0, 0, 0, self.seaLandMask) # valeurs updatés dans la boucle
    
    def integrate(self, date, IntMode, wind, rh, temp, rain):
        
        self.hour = date.hour

        self.rain_buffer[self.hour] = rain
        rain24 = np.sum(self.rain_buffer, axis=0)
        
        self.fwisystem.updateMeteo(temp, rh, wind, rain24)

        if self.needSpinup:
            # fore ffmc
            ii = 0
            diff = 1.e6
            ffmc0 = self.ffmc
            while(diff > 1.e-6):
                ffmc1 = self.fwisystem.hFFMCcalc(ffmc0) 
                diff = ((ffmc1 - ffmc0)**2).sum()
                #print(diff)
                ffmc0 = ffmc1
                ii += 1
            self.ffmc = ffmc1
            self.needSpinup = False

        #3 first indices
        if IntMode == 'hourly':
            self.ffmc = self.fwisystem.hFFMCcalc(self.ffmc,mode=IntMode)
                
        if self.hour == 11:  # if it is 12h00 update dmc and dc:
            if IntMode == 'daily':
                self.ffmc = self.fwisystem.hFFMCcalc(self.ffmc,mode=IntMode)
            self.dmc  = self.fwisystem.DMCcalc(self.dmc,date,self.latitude)
            self.dc   = self.fwisystem.DCcalc(self.dc,date,self.latitude)

        #last three        
        if IntMode == 'hourly':
            self.isi = self.fwisystem.ISIcalc(self.ffmc)
            self.bui = self.fwisystem.BUIcalc(self.dmc,self.dc)
            self.fwi = self.fwisystem.FWIcalc(self.isi, self.bui)
        if IntMode == 'daily':
            if self.hour == 11:
                self.isi = self.fwisystem.ISIcalc(self.ffmc)
                self.bui = self.fwisystem.BUIcalc(self.dmc,self.dc)
                self.fwi = self.fwisystem.FWIcalc(self.isi, self.bui)
    
    def test(self):
        '''
        run test against data available in 
        https://courses.seas.harvard.edu/climate/eli/Courses/global-change-debates/Sources/Forest-fires/aridity-indices/code-for-calculating-canadian-fire-weather-index.pdf
        '''
        
        df_input_test  = pd.read_csv('../dataTest/data.txt', delimiter=' ')
        df_output_test  = pd.read_csv('../dataTest/outputRef.txt', delimiter=' ')

        year=2024 # not important
        
        self.__init__(shape, latitude, 0)
        self.needSpinup = False

        it =0
        for rowIn, rowOut in zip(df_input_test.itertuples(),df_output_test.itertuples()):
            date = datetime(year, rowIn.month, rowIn.day, 11, 0, 0)
            print(date.strftime('%Y-%m-%dT%H:%M:%SZ'), end = ' ' )
            
            wind = np.zeros(shape) + rowIn.wind
            rh   = np.zeros(shape) + rowIn.rh
            temp = np.zeros(shape) + rowIn.temp
            rain = np.zeros(shape) + rowIn.prep
            
            IntMode = 'hourly'
    
            self.integrate(date, IntMode, wind, rh, temp, rain)
            #print(self.fwisystem.h[0,0],
            #     self.fwisystem.t[0,0],
            #     self.fwisystem.w[0,0],
            #     self.fwisystem.p[0,0])
            #print(self.ffmc[0,0],self.dmc[0,0],self.dc[0,0],
            #      self.isi[0,0],self.bui[0,0],self.fwi[0,0])
            #print(rowOut.ffmc, rowOut.dmc,rowOut.dc,rowOut.isi,rowOut.bui,rowOut.fwi)
            print( [ '{:.1f}  '.format(xx) for xx in [
                   abs( np.round(self.ffmc[0,0],1)-rowOut.ffmc), 
                   abs(self.dmc[0,0]-rowOut.dmc),
                   abs(self.dc[0,0] -rowOut.dc),
                   abs(self.isi[0,0]-rowOut.isi),
                   abs(self.bui[0,0]-rowOut.bui),
                   abs(self.fwi[0,0]-rowOut.fwi) ]])

            it+=1


############
def timeIntegration(dirin,dirout,flag_model,iseg):

    date_array, files_array, cexp = dataHandler.load_dates_and_filenames(flag_model, dirin, iseg)
    seaMask,lat2d,lon2d = dataHandler.loadSeaLandMask(flag_model, dirin, files_array)
    latitude =  dataHandler.getMeanLatitdue(flag_model, dirin, files_array)
    
    print ('{:d} files found'.format(len(files_array)))

    times_seconds = np.array([float(xx-date_array[0])/(1.e9) for xx in date_array])
    dt = times_seconds[1]-times_seconds[0]
  
    #init loop
    rain_last24h = []
    time_seconds_last = 0
    #for output
    times_ffmc = []
    times_ffmc_s = []
    ffmc_arr = []
    #dmc_arr = []
    #dc_arr = []
    #isi_arr = []
    #bui_arr = []
    fwi_arr = []

    #check if a data save in FWI calculation from previous forecast can be used.
    #we need the indices and rain_24 at t0 of the current run, as the loop below start at t0+dt
    available_prev_fwi = sorted(glob.glob( dirout+"{:s}_previ_*.pkl".format( (pd.Timestamp(date_array[0])+timedelta(hours=-1)).strftime("%Y%m%dT%H%M")) ))
    
    if len(available_prev_fwi) > 0 : 
        # Load the pickle file
        with open(available_prev_fwi[-1], "rb") as f:
            data_loaded = pickle.load(f)

        # Access the data
        rain_last24h = data_loaded["rain24"]
        ffmc1 = data_loaded["ffmc"]
        dmc1 = data_loaded["dmc"]
        dc1 = data_loaded["dc"]
        fwi1 = data_loaded["fwi"]
        
        times_ffmc.append( pd.Timestamp(date_array[0])+timedelta(hours=-1) )  
        times_ffmc_s.append( float(0.) ) 
        ffmc_arr.append(ffmc1)
        fwi_arr.append(fwi1)

        flag_spinup = False
    else: 
        flag_spinup = True

    for it, (time_seconds, date, file) in enumerate(zip(times_seconds,date_array,files_array)):
        
        if it > 0: 
            flag_spinup = False

        print(pd.Timestamp(date).strftime("%Y%m%dT%H%M"), end=' ')
        wind, humidity, temperature, rain = dataHandler.load_data_for_date(flag_model,date, file, it)
        #wind        # km/h
        #humidity    # % (<100)
        #temperature # C
        #rain        # mm
        
        if (it ==0): 
            if (flag_spinup): 
                mask_data = np.where((seaMask==1) & (~wind.mask), 1., 0 )
                ffmc00 = mask_data * 20.0 # valeurs prises dans la func initialisation plus bas. merci de valider 
                dmc00  = mask_data *  6.0
                dc00   = mask_data * 15.0
            
        rain_last24h.append(rain)
        if len(rain_last24h)>24: 
            rain_last24h = rain_last24h[1:]
        print( ' lenRain24={:2d}'.format(len(rain_last24h)), end='' )
        rain0       = np.array(rain_last24h).sum(axis=0) # mm 
        #-- note: in the intergaration, before the end of the first day we are missing rain from the previous day
        
        date_ = datetime.strptime(str(date)[:-3], "%Y-%m-%dT%H:%M:%S.%f")
        #time_seconds_since00 = (date_.hour * 3600) + (date_.minute * 60) + date_.second
        
        humidity  = np.where(humidity>100,100.,humidity)
                
        fwisystem = FWICLASS(temperature,humidity,wind,rain0,seaMask)

        #compute time integration of ffmc, dmc, and dc
        #----------------
        
        #get previous time data if present
       
        
        if it > 0:
            if np.abs((time_seconds-3600) - np.array(times_ffmc_s)).min() != 0 : #-1h data exist, otherwise do spinup
                flag_spinup = True
        
        if not(flag_spinup) : #-1h data exist, otherwise do spinup
            ithm1 = np.abs((time_seconds-3600) - np.array(times_ffmc_s)).argmin()
            ffmc0 = ffmc_arr[ithm1]
            dmc0  = dmc1 #dmc_arr[-1] # for dmc0 and dc0 we use always the last one to update at 12h
            dc0   = dc1  #dc_arr[-1]
            if it == 0 : 
                print(' init from prev calculation: {:s}'.format(os.path.basename(available_prev_fwi[-1]).split('_')[-1].split('.pk')[0] ), end ='')
        else: 
            ffmc0 = ffmc00
            maskLand_ = fwisystem.mask
            idx = np.where(maskLand_==1)
            #H_  = np.where(fwisystem.h[idx]>100,100.,fwisystem.h[idx])
            #T_  = fwisystem.t[idx]
            #ffmc0 = np.zeros(maskLand_.shape)-999
            #ffmc0[idx] = 0.5*(get_Ed(H_,T_) + get_Ew(H_,T_))
            dmc0  = dmc00
            dc0   = dc00
            print (' spinOn, ffmc00 = ', ffmc0[idx].mean(), end ='\n')
       
        #if no previous, run spinup
        if flag_spinup: 
            # fore ffmc
            ii = 0
            diff = 1.e6
            diff_threshold = 1.e-3
            while(diff > diff_threshold):
                ffmc1 = fwisystem.hFFMCcalc(ffmc0) 
                diff = ((ffmc1 - ffmc0)**2).sum() / ffmc1.flatten().shape[0]
                #print('spinup while {:.3e} > {:.3e}  '.format(diff,diff_threshold), end='\r')
                ffmc0 = ffmc1
                ii += 1
            print('spinup done                                ', end='')
            dmc1 = dmc0
            dc1  = dc0
        
        else:
            #ffmc1 =  0.75 * fwisystem.hFFMCcalc(ffmc0) +\
            #         0.25 * ffmc_arr[-1]
            ffmc1 =  fwisystem.hFFMCcalc(ffmc0) 
            if  pd.Timestamp(date).hour == 12 : #abs( time_seconds_since00 - (12*3600) ) < 1:  # if it is 12h00 update dmc and dc:
                print(' update dmc, dc', end='')
                try: 
                    dmc1  = fwisystem.DMCcalc(dmc0,date_,latitude)
                    dc1   = fwisystem.DCcalc(dc0,date_,latitude)
                except: 
                    pdb.set_trace()
            else: 
                dmc1 = dmc0
                dc1  = dc0
        
        print('')

        isi1  = fwisystem.ISIcalc(ffmc1)
        bui1  = fwisystem.BUIcalc(dmc1,dc1)
        
        fwi1 = fwisystem.FWIcalc(isi1, bui1) 
       
        time_seconds_last = time_seconds
        
        times_ffmc.append(date) 
        times_ffmc_s.append( float(times_ffmc[-1]-date_array[0])/(1.e9) ) 

        ffmc_arr.append(np.where((seaMask==1)& (~wind.mask),ffmc1,np.nan))
        #dmc_arr.append( np.where((seaMask==1)& (~wind.mask),dmc1,np.nan))
        #dc_arr.append(  np.where((seaMask==1)& (~wind.mask),dc1,np.nan))
        #isi_arr.append( np.where((seaMask==1)& (~wind.mask),isi1,np.nan))
        #bui_arr.append( np.where((seaMask==1)& (~wind.mask),bui1,np.nan))
        fwi_arr.append( np.where((seaMask==1)& (~wind.mask),fwi1,np.nan))

   
        #at time of new previ we save rain_last24h and indices
        if pd.Timestamp(date).hour % 3 == 0: 
            # Combine everything in a dictionary (recommended for clarity)
            data_to_pickle = {
                "rain24": rain_last24h, 
                "ffmc":   ffmc_arr[-1],
                "dmc":    dmc1,
                "dc":     dc1,
                "fwi":    fwi_arr[-1]
            }

            # Dump to pickle
            with open(dirout+"{:s}_previ_{:s}.pkl".format(pd.Timestamp(date).strftime("%Y%m%dT%H%M"),cexp), "wb") as f:
                pickle.dump(data_to_pickle, f)
            del data_to_pickle

        fwisystem.self_destruct()

    # Combine all arrays into one dictionary
    data_dict = {
        "FFMC": ffmc_arr,
        #"DMC": dmc_arr,
        #"DC":  dc_arr,
        #"ISI": isi_arr,
        #"BUI": bui_arr,
        "FWI": fwi_arr
    }
    del ffmc_arr, fwi_arr

    # Create an xarray.DataArray for each variable
    data_arrays = {
        key: xr.DataArray(
            data=np.stack(value),  # Stack 2D arrays along a new dimension (time)
            dims=["time", "lat", "lon"],  # Dimensions: time, y (rows), x (cols)
            coords={"time": times_ffmc,
                    "lon": lon2d[0,:].astype(np.float32),  # Add 2D longitude
                    "lat": lat2d[:,0].astype(np.float32),  # Add 2D latitude
                    },
            name=key,
        )
        for key, value in data_dict.items()
    }
    del data_dict

    # Combine into a single xarray.Dataset
    ds = xr.Dataset(data_arrays)
    ds['FFMC'] = ds['FFMC'].astype('float32')
    ds['FWI'] = ds['FWI'].astype('float32')
    ds.rio.write_crs(4326) 
    ds.lon.attrs['long_name'] = "Longitude"
    ds.lon.attrs['units'] = "degrees_east"
    ds.lon.attrs['axis'] = "X"
    ds.lon.attrs['standard_name'] = "longitude"

    ds.lat.attrs['long_name'] = "Latitude"
    ds.lat.attrs['units'] = "degrees_north"
    ds.lat.attrs['axis'] = "Y"
    ds.lat.attrs['standard_name'] = "latitude"

    del data_arrays
    gc.collect()   # force garbage collection

    return ds.rio.write_crs(4326), cexp
  

###########################
def get_two_level_dirs(base_path):
    base = Path(base_path)
    out = sorted([
        p for p in base.rglob("*")
        if p.is_dir()
        and len(p.relative_to(base).parts) <= 2
        and not (len(p.relative_to(base).parts) == 1 and p.name == "log")
    ])
    return out


###########################
if __name__ == '__main__':
###########################
    
    importlib.reload(dataHandler)

    flag_test = False

    if flag_test:
        shape=[1,1]
        latitude = 55
        indices = Indices(shape, latitude, 0)
        indices.test()
        sys.exit() 
    
    #dirin = '/data/paugam/MNH/Cat_PdV/006_mnhSolo'
    #flag_model = 'mnh'
    #iseg = 1
    '''
    dirin = '/mnt/data3/SILEX/AROME/20250413/18Z/'
    #dirin = '/mnt/dataEstrella2/SILEX/AROME/20250413/21Z/'
    flag_model = 'arome'
    iseg = -9999
    dataset, cexp = timeIntegration(dirin,flag_model,iseg)

    os.makedirs(dirin+'/Postproc/FWI/', exist_ok=True)
    dataset.to_netcdf(dirin+'/Postproc/FWI/{:s}_fwi.nc'.format(cexp))
    '''
    
    dirArome = sys.argv[1]  #'/mnt/data3/SILEX/AROME/'
    dirout2 = dirArome+'FWI/'
    os.makedirs(dirout2, exist_ok=True)
    flag_model = 'arome'
    iseg = -9999
    
    #get thelast complete download
    #existingdir_sorted = get_two_level_dirs(dirArome+'FORECAST/')
    #for dir_ in existingdir_sorted[::-1]:
    #    last_forecastTime = datetime.strptime('T'.join(str(dir_).split('/')[-2:]), '%Y%m%dT%HZ')
    #    dirout_ = '{:s}/{:s}/'.format(dirArome+'FORECAST/',last_forecastTime.strftime('%Y%m%d/%HZ'))
    #    if len(glob.glob(dirout_+'*.grib2')) == 104:
    #        break
    with open(dirout2+'timeToCompute.txt', 'r') as f:
        time_str = f.read().strip()
    
    start = datetime.strptime(time_str, '%Y%m%dT.%H%MZ')
    step  = timedelta(hours=3)
    end   = start + step 
    if os.path.isfile(dirout2+'/{:s}_fwiffmc.nc'.format(start.strftime("%Y%m%d.%HZ"))):
        sys.exit()

    #start = datetime.strptime("20250412T0000", "%Y%m%dT%H%M")
    #end   = datetime.strptime("20250413T2100", "%Y%m%dT%H%M")
    #step  = timedelta(hours=3)



    tracemalloc.start()

    currentZ = start
    while currentZ < end:
       
        start_time = time.time()  # Record start time
        if False:# os.path.isfile(dirout2+'/{:s}_fwiffmc.nc'.format(currentZ.strftime("%Y%m%d.%HZ"))): 
            #next step
            currentZ += step
            continue
        
        print('##########')
        print('FORECAST is now: ', currentZ.strftime("%Y%m%dT%H%M"))
        print('##########')

        dirprevi = currentZ.strftime("%Y%m%d/%HZ/")
        dirin = dirArome + 'FORECAST/' + dirprevi
        
        dirout1 = dirArome+'IntermediateIndicesRainAcc/'
        os.makedirs(dirout1, exist_ok=True)
       
        #FWI Calculation
        dataset, cexp = timeIntegration(dirin,dirout1,flag_model,iseg)
        dataset.to_netcdf(dirout2+'/{:s}_fwiffmc.nc'.format(cexp))
        del dataset 

        #collection time of file in dirout1 and remove all older that currentZ
        files = glob.glob(dirout1+'*.pkl')
        date_saved = []
        for file in files:
            date_saved.append(  datetime.strptime(os.path.basename(file).split('_p')[0], "%Y%m%dT%H%M") ) 
        date_saved = np.array(date_saved)
        idx_ = np.where(date_saved<currentZ)[0]
        if len(idx_)>0:
            for file in np.array(files)[idx_]:
                os.remove(file)
        #next step
        currentZ += step

    # Get the current process
    process = psutil.Process(os.getpid())
    # Get memory usage in GB
    ram_used_gb = process.memory_info().rss / 1024**3
    print(f"RAM used by the process: {ram_used_gb:.4f} GB")
    
    end_time = time.time()  # Record end time
    elapsed_time = (end_time - start_time)/60  # Time used in iteration
    print(f"cput time by FORECAST  : {elapsed_time:.6f} minutes")      

    plot_fwi_arome.plot_latest(dirArome)
    
        
