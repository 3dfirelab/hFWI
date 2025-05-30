import xarray as xr 
import matplotlib as mpl 
import matplotlib.pyplot as plt 
import sys
import os 
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import rioxarray
import numpy as np 
from pyproj import Transformer
import glob 
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
import shutil

#############################
def plot_latest(dirArome):
    dirin = dirArome+'FWI/'

    dirout = dirin+'/PLOT_latest/png/'
    flag_model = 'arome'
    if os.path.isdir(dirout):
        shutil.rmtree(dirout)
    os.makedirs(dirout, exist_ok=True)
    
    crshere = 3035

    transformer = Transformer.from_crs("EPSG:4326", "EPSG:32631", always_xy=True)

    fwifiles = sorted(glob.glob(dirin+'*.nc'))
       
    print( 'plot png in PLOT_latest dir for :', os.path.basename(fwifiles[-1]))
    ds = xr.open_dataset(fwifiles[-1])
    ds = ds.rio.write_crs(4326)
    ds = ds.rio.reproject(crshere)
    
    projection = ccrs.epsg(ds.rio.crs.to_epsg())
    levels = [0.0, 5.2, 11.2, 21.3, 38.0, 50.0,100]
    colors = ["#008000", "#FFFF00", "#FFA500", "#FF0000", "#654321", "#000000"]
    labels = ['Very low', 'Low', 'Moderate', 'High', 'Very high', 'Extreme']

    mpl.rcdefaults()
    mpl.rcParams['text.usetex'] = False
    mpl.rcParams['figure.subplot.left'] = .05
    mpl.rcParams['figure.subplot.right'] = .83
    mpl.rcParams['figure.subplot.top'] = .92
    mpl.rcParams['figure.subplot.bottom'] = .05

    for it in range(ds.time.shape[0]):
        fig, ax = plt.subplots(figsize=(10, 9), subplot_kw={'projection': projection})
        
        # Create a divider for the existing axes instance
        #divider = make_axes_locatable(ax)

        # Append a new axes to the right of ax with a fixed size and padding
        #cax = divider.append_axes("right", size="5%", pad=0.05)

        plot = ds.FWI.isel(time=it).plot(
            levels=levels,
            colors=colors,
            add_colorbar=False,
            ax=ax
        )

        # Add custom colorbar
        cax = fig.add_axes([0.84, 0.1, 0.03, 0.76])  # adjust as needed

        cbar = fig.colorbar(plot, cax=cax, ticks=[(levels[i] + levels[i+1]) / 2 for i in range(len(levels)-1)])
        cbar.ax.set_yticklabels(labels)  # or use set_ticks + set_ticklabels for horizontal
        cbar.set_label("FWI Class")  # Optional: label for the colorbar

        ax.coastlines(resolution='10m', color='black', linewidth=1)
        ax.add_feature(cfeature.BORDERS, linestyle=':')
        ax.add_feature(cfeature.LAND, facecolor='lightgray')
        ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
        gridlines = ax.gridlines(draw_labels=True, color='gray', alpha=0.5, linestyle='--')
        gridlines.top_labels = False  # Turn off top labels
        gridlines.right_labels = False  # Turn off right labels
        gridlines.xlabel_style = {'size': 10, 'color': 'black'}
        gridlines.ylabel_style = {'size': 10, 'color': 'black'}
        ax.set_title('FWI-{:s}'.format(pd.to_datetime(ds.time[it].item()).strftime("%Y-%m-%dT%H%MZ (UTC)")))
    
        fig.savefig('{:s}/FWI-{:s}.png'.format(dirout,pd.to_datetime(ds.time[it].item()).strftime("%Y-%m-%dT%H%M")))
        plt.close(fig)


###########################
if __name__ == '__main__':
###########################
    dirArome = '/mnt/data3/SILEX/AROME/'
    #dirArome = '/mnt/dataEstrella2/SILEX/AROME/'
    
    plot_latest(dirArome)

    sys.exit()
    dirin = dirArome+'FWI/'

    dirout = dirin+'/PLOT/'
    flag_model = 'arome'
    os.makedirs(dirout, exist_ok=True)
    crshere = 3035

    transformer = Transformer.from_crs("EPSG:4326", "EPSG:32631", always_xy=True)

    
    
    #ffmcmin=np.nanpercentile(ds.FFMC,5)
    #ffmcmax=np.nanpercentile(ds.FFMC,95)
    #fwimin=np.nanpercentile(ds.FWI,20)
    #fwimax=np.nanpercentile(ds.FWI,80)
    
    tab10 = plt.get_cmap("tab10")

    fwifiles = sorted(glob.glob(dirin+'*.nc'))
    if os.path.isfile(dirout+'dsa.nc'):
        print('load existing concatenated file:')
        dsa =  xr.open_dataset(dirout+'dsa.nc')
        #dsa = dsa.rio.write_crs(crshere)

    else:
        dsa = None
        for ii, (fwifile,fwifile_next) in enumerate(zip(fwifiles[:-1],fwifiles[1:])):
       
            print( ii, os.path.basename(fwifile))
            if ii == 0 : 
                ds = xr.open_dataset(fwifile)
                ds = ds.rio.write_crs(4326)
                #ds = ds.rio.reproject(crshere)
            else:
                ds = dsn.copy()
                
            dsn = xr.open_dataset(fwifile_next)
            dsn = dsn.rio.write_crs(4326)
            #dsn = dsn.rio.reproject(crshere)

            if dsa is None:
                dsa = ds.sel(time=slice(pd.Timestamp(ds.time.min().values), pd.Timestamp(dsn.time.min().values)))
            else:
                dsa = xr.concat([dsa, ds.sel(time=slice(pd.Timestamp(ds.time[ds.time.argmin()+1].values), pd.Timestamp(dsn.time.min().values))) ] , dim='time')
        
            if ii == len(fwifiles)-2: 
                dsa = xr.concat([dsa, dsn.sel(time=slice(pd.Timestamp(dsn.time[dsn.time.argmin()+1].values), pd.Timestamp(dsn.time.max().values))) ] , dim='time')
       
        dsa.to_netcdf(dirout+'dsa.nc')

    ###############
    #plot time series 
    ###############
    ds1 = xr.open_dataset(fwifiles[0])
    ds1 = ds1.rio.write_crs(4326)

    lat_,lon_ = 45.6, 1.61

    timeSerie_dsa = dsa.sel(lat=lat_, lon=lon_)
    timeSerie_ds1 = ds1.sel(lat=lat_, lon=lon_)

    fig =plt.figure(figsize=(12,6))
    ax = plt.subplot(211)
    timeSerie_dsa.FWI.plot(ax=ax,label='concatenated')
    timeSerie_ds1.FWI.plot(ax=ax,label='{:s}'.format(os.path.basename(fwifiles[0])))
    ax.legend()
    ax = plt.subplot(212)
    timeSerie_dsa.FFMC.plot(ax=ax)
    timeSerie_ds1.FFMC.plot(ax=ax)
    fig.savefig(dirout+'timeSeries-fwi-ffmc-concateneated-firstIntgration.png')
    plt.close(fig)


    ###############
    #plot 2d
    ###############
    ds_copernicus = xr.open_dataset(dirout+'fwi_copernicus.nc')
    # Adjust longitude to the [-180, 180] range
    ds_copernicus = ds_copernicus.assign_coords(
        longitude=(((ds_copernicus.longitude + 180) % 360) - 180)
    )
    # Optional: sort by the new longitude values
    ds_copernicus = ds_copernicus.sortby('longitude')
    ds_copernicus = ds_copernicus.sortby('latitude')

    ds_copernicus_croped = ds_copernicus.sel(latitude=slice(dsa.lat.min(),dsa.lat.max()), longitude=slice(dsa.lon.min(),dsa.lon.max()))

    #projection = ccrs.epsg(dsa.rio.crs.to_epsg())
    projection = ccrs.PlateCarree()
    for ii, time in enumerate(ds_copernicus_croped.valid_time): 

        it = np.abs(dsa.time-(ds_copernicus_croped.valid_time[ii]+ np.timedelta64(12, 'h'))).argmin().item()

        #FFMC
        fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(figsize=(12, 12), ncols=2, nrows=2, subplot_kw={'projection': projection})

        
        dsa.FFMC.isel(time=it).plot(ax=ax1,vmin=10,vmax=90)
        ax1.coastlines(resolution='10m', color='black', linewidth=1)
        ax1.add_feature(cfeature.BORDERS, linestyle=':')
        ax1.add_feature(cfeature.LAND, facecolor='lightgray')
        ax1.add_feature(cfeature.OCEAN, facecolor='lightblue')
        gridlines = ax1.gridlines(draw_labels=True, color='gray', alpha=0.5, linestyle='--')
        gridlines.top_labels = False  # Turn off top labels
        gridlines.right_labels = False  # Turn off right labels
        gridlines.xlabel_style = {'size': 10, 'color': 'black'}
        gridlines.ylabel_style = {'size': 10, 'color': 'black'}
        ax1.set_title('FFMC-hFWI-{:s}'.format(pd.to_datetime(dsa.time[it].item()).strftime("%Y-%m-%dT%H%M")))

        ds_copernicus_croped.ffmcode.isel(valid_time=ii).plot(ax=ax2,vmin=10,vmax=90)
        ax2.coastlines(resolution='10m', color='black', linewidth=1)
        ax2.add_feature(cfeature.BORDERS, linestyle=':')
        ax2.add_feature(cfeature.LAND, facecolor='lightgray')
        ax2.add_feature(cfeature.OCEAN, facecolor='lightblue')
        gridlines = ax2.gridlines(draw_labels=True, color='gray', alpha=0.5, linestyle='--')
        gridlines.top_labels = False  # Turn off top labels
        gridlines.right_labels = False  # Turn off right labels
        gridlines.xlabel_style = {'size': 10, 'color': 'black'}
        gridlines.ylabel_style = {'size': 10, 'color': 'black'}
        ax2.set_title('FFMC-Copernicus-{:s}'.format(pd.to_datetime(dsa.time[it].item()).strftime("%Y-%m-%dT%H%M")))
        
        dsa.FWI.isel(time=it).plot(ax=ax3,vmin=0,vmax=15)
        ax3.coastlines(resolution='10m', color='black', linewidth=1)
        ax3.add_feature(cfeature.BORDERS, linestyle=':')
        ax3.add_feature(cfeature.LAND, facecolor='lightgray')
        ax3.add_feature(cfeature.OCEAN, facecolor='lightblue')
        gridlines = ax3.gridlines(draw_labels=True, color='gray', alpha=0.5, linestyle='--')
        gridlines.top_labels = False  # Turn off top labels
        gridlines.right_labels = False  # Turn off right labels
        gridlines.xlabel_style = {'size': 10, 'color': 'black'}
        gridlines.ylabel_style = {'size': 10, 'color': 'black'}
        ax3.set_title('FWI-hFWI-{:s}'.format(pd.to_datetime(dsa.time[it].item()).strftime("%Y-%m-%dT%H%M")))

        ds_copernicus_croped.fwinx.isel(valid_time=ii).plot(ax=ax4,vmin=0,vmax=15)
        ax4.coastlines(resolution='10m', color='black', linewidth=1)
        ax4.add_feature(cfeature.BORDERS, linestyle=':')
        ax4.add_feature(cfeature.LAND, facecolor='lightgray')
        ax4.add_feature(cfeature.OCEAN, facecolor='lightblue')
        gridlines = ax4.gridlines(draw_labels=True, color='gray', alpha=0.5, linestyle='--')
        gridlines.top_labels = False  # Turn off top labels
        gridlines.right_labels = False  # Turn off right labels
        gridlines.xlabel_style = {'size': 10, 'color': 'black'}
        gridlines.ylabel_style = {'size': 10, 'color': 'black'}
        ax4.set_title('FWI-Copernicus-{:s}'.format(pd.to_datetime(dsa.time[it].item()).strftime("%Y-%m-%dT%H%M")))

        fig.savefig(dirout+'concatenated-copernicus-{:s}.png'.format(pd.to_datetime(dsa.time[it].item()).strftime("%Y-%m-%dT%H%M")))
        plt.close(fig)
        #plt.show()

        #fig.savefig('{:s}/ffmcfwi_{:s}_{:06d}.png'.format(diroutpdV,cexp,it))
        #plt.close(fig)

