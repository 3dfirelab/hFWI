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

import dataHandler

class FWICLASS:

    ############
    def __init__(self,temp,rhum,wind,prcp,mask):
        self.updateMeteo(temp,rhum,wind,prcp)
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
        Ed = 0.942*(H**0.679) + 11*(np.exp((H - 100)/10)) + 0.18*(21.1 - T)*(1 - np.exp(-0.115*H))
        Ew = 0.618*(H**0.753) + 10*(np.exp((H - 100)/10)) + 0.18*(21.1 - T)*(1 - np.exp(-0.115*H))

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
        El = tools.DayLength(latitude, day.month)

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

        Lf = tools.DryingFactor(latitude, day.month)

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
        
        df_input_test  = pd.read_csv('dataTest/data.txt', delimiter=' ')
        df_output_test  = pd.read_csv('dataTest/outputRef.txt', delimiter=' ')

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
            print( abs(self.ffmc[0,0]-rowOut.ffmc), 
                   abs(self.dmc[0,0]-rowOut.dmc),
                   abs(self.dc[0,0] -rowOut.dc),
                   abs(self.isi[0,0]-rowOut.isi),
                   abs(self.bui[0,0]-rowOut.bui),
                   abs(self.fwi[0,0]-rowOut.fwi) )

            it+=1


if __name__ == '__main__':
    

    shape=[1,1]
    latitude = 55
    indices = Indices(shape, latitude, 0)
    indices.test()





'''

############
def timeIntegration(day, 
                    latitude):

    date_array = load_dates(day)
    
    #shape   = wind_array.shape[1:]
    
    #out = np.zeros(timeRain0.shape, dtype=np.dtype([('ffmc',float),('dmc',float),('dc',float),
    #                                                ('isi',float),('bui',float),('fwi',float)]))
    #out = out.view(np.recarray)

    times_seconds = np.array([(xx-date_array[0]).total_seconds() for xx in date_array])
    dt = times_seconds[1]-times_seconds[0]
    
    rain_last24h = []
    
    for it, (time_seconds, date) in enumerate(zip(times_seconds,date_array)):

        print(date)
        
        wind, humidity, temperature, rain = dataHandler.load_data_for_date(date)
        #wind        # km/h
        #humidity    # % (<100)
        #temperature # C
        #rain        # mm
        
        rain_last24h.append(rain)
        if len(rain_last24h)>24: 
            rain_last24h = rain_last24[1:]
        rain0       = np.array(rain_last24h).sum(axis=0) # mm 
        #-- note: in the intergaration, before the enf of the first day we are missing rain from the previous day

        time_seconds = (date.hour * 3600) + (date.minute * 60) + date.second
        
        humidity  = np.where(humidity>100,100.,humidity)
                
        fwisystem = FWICLASS(temperature,humidity,wind,rain0)

        #compute time integration of ffmc, dmc, and dc
        #----------------
        
        #get previous time data if present
        
        if it  > 0:
            ffmc0 = ffmc1
            dmc0  = dmc1
            dc0   = dc1
            flag_spinup = False

        else: 
            ffmc0 = ffmc00
            dmc0  = dmc00
            dc0   = dc00
            flag_spinup = True
        
        #if no previous, run spinup
        if flag_spinup: 
            # fore ffmc
            ii = 0
            diff = 1.e6
            while(diff > 1.e-6):
                ffmc1 = fwisystem.hFFMCcalc(ffmc0) 
                diff = ((ffmc1 - ffmc0)**2).sum()
                ffmc0 = ffmc1
                ii += 1
            dmc1 = dmc0
            dc1  = dc0
        
        else:
            ffmc1 = fwisystem.hFFMCcalc(ffmc0)
            if abs( time_seconds - (12*3600) ) < 1:  # if it is 12h00 update dmc and dc:
                dmc1  = fwisystem.DMCcalc(dmc0,date,latitude)
                dc1   = fwisystem.DCcalc(dc0,date,latitude)
            else: 
                dmc1 = dmc0
                dc1  = dc0

        isi1  = fwisystem.ISIcalc(ffmc1)
        bui1  = fwisystem.BUIcalc(dmc1,dc1)
        
        fwi1 = fwisystem.FWIcalc(isi1, bui1) 
        
        #out['ffmc'][it] = ffmc1
        #out['dmc'][it]  = dmc1
        #out['dc'][it]   = dc1

        #out['isi'][it]  = fwisystem.ISIcalc(ffmc1)
        #out['bui'][it]  = fwisystem.BUIcalc(dmc1,dc1)
        
        #out['fwi'][it]  = fwisystem.FWIcalc(out['isi'][it], out['bui'][it])

    
    return ffmc1, dmc1, dc1, isi1, bui1, fwi1
    

def load_dates(dayStart, nbreMonth=1):
    out = []
    date = dayStart
    out.append(date)
    for i in range(int(30*nbreMonth*24)):
        date = date - timedelta(hours=1)
        out.append(date)

    return np.array(out)[::-1]
    

############
def initialization(dayStart,
                   latitude):
    
    files = glob.glob() # that depends of dayStart -- first fire is at 00h00
    dt = 3600  # files are provided every hour
    
    nbre_day_since_start = 0
    for ifile, file_ in files: 

        #form file:
        rain0       =              # mm per hour
        wind        =              # km/h
        humidity    =              # % (<100)
        temperature =              # C
        time        =              # time of the simulation - dattime

        if ifile_ == 0: 
            ffmc00 = np.zeros(shape) + 85.0
            dmc00  = np.zeros(shape) + 6.0
            dc00   = np.zeros(shape) + 15.0

            daily_rain_accumulated = np.zeros(shape)
            time00 = time
            
        #mth,day,temp,rhum,wind,prcp=[float(field) for field in line.strip().lstrip('[').rstrip(']').split()]
        humidity  = np.where(humidity>100,100.,humidity)
        nbre_day_since_start_ =  int((time-time00).total_seconds()/(3600*24))  #assume first files is at 00h00

        flag_new_day = False
        if nbre_day_since_start_ > nbre_day_since_start: 
            nbre_day_since_start = True
            nbre_day_since_start = nbre_day_since_start_

        if flag_new_day: 
            daily_rain_accumulated = np.zeros(shape)
        else: 
            daily_rain_accumulated +=  rain0 * dt
        
        fwisystem = FWICLASS(temperature,humidity,wind,daily_rain_accumulated)

        #compute time integration of ffmc
        #----------------
        
        #get previous time data if present
        if ifile_ > 0: 
            ffmc0 = ffmc_prev
            dmc0  = dmc_prev
            dc0   = dc_prev 
            flag_spinup = False

        else: 
            ffmc0 = ffmc00
            dmc0  = dmc00
            dc0   = dc00
            flag_spinup = True
                
        if flag_spinup: 
            # fore ffmc
            ii = 0
            diff = 1.e6
            while(diff > 1.e-6):
                ffmc1 = fwisystem.hFFMCcalc(ffmc0) 
                diff = ((ffmc1 - ffmc0)**2).sum()
                ffmc0 = ffmc1
                ii += 1
                    
        ffmc_prev = fwisystem.hFFMCcalc(ffmc0)
        dmc_prev = fwisystem.DMCcalc(dmc0,time,latitude)
        dc_prev = fwisystem.DCcalc(dc0,time,latitude)

    
    return ffmc_prev, dmc_prev, dc_prev



if __name__ == '__main__':


    grbs = pygrib.open('./dataTest/arome__001__HP1__01H__2024-04-08T00 00 00Z.grib2')
    wind_ = grbs.select()[4].values
    rh_   = grbs.select()[8].values
    temp_ = np.zeros_like(rh_)+24
    rain_ = np.zeros_like(rh_)
    
    wind = np.zeros([2,wind_.shape[0],wind_.shape[1]])
    rh   = np.zeros_like(wind) 
    temp = np.zeros_like(wind) 
    rain = np.zeros_like(wind) 

    wind[0] = wind_ 
    rh[0]   = rh_
    temp[0] = temp_
    rain[0] = rain_

    latitude = 45

    times = np.array([3600,7200])

    out = timeIntegration(times, 
                    rain, wind, 
                    rh, temp, 
                    latitude)



def main():
    ffmc0 = 85.0
    dmc0 = 6.0
    dc0 = 15.0
    infile = open('data.txt','r')
    outfile = open('fwioutput.txt','w')
    try:
        for line in infile:
            mth,day,temp,rhum,wind,prcp=[float(field) for field in line.strip().lstrip('[').rstrip(']').split()]
            if rhum>100.0:
                rhum = 100.0
            mth = int(mth)
            
            fwisystem = FWICLASS(temp,rhum,wind,prcp)
            
            ffmc = fwisystem.FFMCcalc(ffmc0)
            dmc  = fwisystem.DMCcalc(dmc0,mth)
            dc   = fwisystem.DCcalc(dc0,mth)
            isi  = fwisystem.ISIcalc(ffmc)
            bui  = fwisystem.BUIcalc(dmc,dc)
            fwi  = fwisystem.FWIcalc(isi,bui)
            
            ffmc0 = ffmc
            dmc0  = dmc
            dc0   = dc

            outfile.write("%s %s %s %s %s %s\n"%(str(ffmc),str(dmc),str(dc),str(isi),str(bui),str(fwi)))
    
    finally:
        infile.close()
        outfile.close()

#main()
'''
