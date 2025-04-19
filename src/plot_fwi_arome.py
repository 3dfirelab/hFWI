import xarray as xr 
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

###########################
if __name__ == '__main__':
###########################
    dirArome = '/mnt/data3/SILEX/AROME/'
    dirin = dirArome+'FWI/'

    dirout = dirin+'/PNG/'
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
    dsa = None
    for ii, (fwifile,fwifile_next) in enumerate(zip(fwifiles[:-1],fwifiles[1:])):
   
        print( ii, os.path.basename(fwifile))
        if ii == 0 : 
            ds = xr.open_dataset(fwifile)
            ds = ds.rio.write_crs(4326)
            ds = ds.rio.reproject(crshere)
        else:
            ds = dsn.copy()
            
        dsn = xr.open_dataset(fwifile_next)
        dsn = dsn.rio.write_crs(4326)
        dsn = dsn.rio.reproject(crshere)

        if dsa is None:
            dsa = ds.sel(time=slice(pd.Timestamp(ds.time.min().values), pd.Timestamp(dsn.time.min().values)))
        elif ii == len(fwifile)-1: 
            dsa = xr.concat([dsa, ds.sel(time=slice(pd.Timestamp(ds.time[ds.time.argmin()+1].values), pd.Timestamp(dsn.time.max().values))) ] , dim='time')
        else:
            dsa = xr.concat([dsa, ds.sel(time=slice(pd.Timestamp(ds.time[ds.time.argmin()+1].values), pd.Timestamp(dsn.time.min().values))) ] , dim='time')
    
   
    sys.exit()
    projection = ccrs.epsg(ds.rio.crs.to_epsg())
    
    for it, time in enumerate(ds.time): 
        #FFMC
        fig, ax = plt.subplots(figsize=(9, 7), subplot_kw={'projection': projection})

        ds.FFMC.isel(time=it).plot(ax=ax,vmin=ffmcmin,vmax=ffmcmax)

        ax.coastlines(resolution='10m', color='black', linewidth=1)
        ax.add_feature(cfeature.BORDERS, linestyle=':')
        ax.add_feature(cfeature.LAND, facecolor='lightgray')
        ax.add_feature(cfeature.OCEAN, facecolor='lightblue')

        gridlines = ax.gridlines(draw_labels=True, color='gray', alpha=0.5, linestyle='--')
        gridlines.top_labels = False  # Turn off top labels
        gridlines.right_labels = False  # Turn off right labels
        gridlines.xlabel_style = {'size': 10, 'color': 'black'}
        gridlines.ylabel_style = {'size': 10, 'color': 'black'}

        ax.set_title(str(ds.time[it].data).split('.000')[0], fontsize=14)

        ax.scatter(pdV_easting, pdV_northing,c='r')

        fig.savefig('{:s}ffmc_{:s}_{:06d}.png'.format(diroutffmc,cexp,it))
        plt.close(fig)

        
        #FWI
        fig, ax = plt.subplots(figsize=(9, 7), subplot_kw={'projection': projection})

        #ds.FWI.isel(time=it).plot(ax=ax,cmap='Reds',vmin=fwimin,vmax=fwimax)

        ds.FWI.isel(time=it).plot(ax=ax,
            levels = [0.0, 5.2, 11.2, 21.3, 38.0, 50.0],
            colors = ["#008000", "#FFFF00", "#FFA500", "#FF0000", "#654321", "#000000"],
            label = ['Very low', 'Low', 'Moderate', 'High', 'Very high', 'Extreme'],
        )

        ax.coastlines(resolution='10m', color='black', linewidth=1)
        ax.add_feature(cfeature.BORDERS, linestyle=':')
        ax.add_feature(cfeature.LAND, facecolor='lightgray')
        ax.add_feature(cfeature.OCEAN, facecolor='lightblue')

        gridlines = ax.gridlines(draw_labels=True, color='gray', alpha=0.5, linestyle='--')
        gridlines.top_labels = False  # Turn off top labels
        gridlines.right_labels = False  # Turn off right labels
        gridlines.xlabel_style = {'size': 10, 'color': 'black'}
        gridlines.ylabel_style = {'size': 10, 'color': 'black'}

        ax.set_title(str(ds.time[it].data).split('.000')[0], fontsize=14)
        
        ax.scatter(pdV_easting, pdV_northing,c='k')

        fig.savefig('{:s}fwi_{:s}_{:06d}.png'.format(diroutfwi,cexp,it))
        plt.close(fig)


        #Point Value
        fig = plt.figure(figsize=(10,6))
        ax = plt.subplot(111)
        bx= ax.twinx()

        line1 = ds.FFMC.sel(x=pdV_easting,y=pdV_northing,method='nearest').plot(ax=ax)
        ax.axes.set_title("")  # Remove the title


        line2 = ds.FWI.sel(x=pdV_easting,y=pdV_northing,method='nearest').plot(ax=bx,color=tab10(1))
        bx.axes.set_title("")  # Remove the title
    
        ax.scatter(time, ds.FFMC.sel(x=pdV_easting,y=pdV_northing,time=time,method='nearest'),color='k')
        bx.scatter(time, ds.FWI.sel(x=pdV_easting,y=pdV_northing,time=time,method='nearest'),color='k')
        
        ax.set_title('time serie of FFCM and FWI at el Pont de Vilomara' , fontsize=14)
 
        ax.annotate(str(ds.time[it].data).split('.000')[0], 
             xy=(0.87, 0.05),  # Position of the annotation (in relative coordinates)
             xycoords='axes fraction',  # Coordinates are in fraction of the axes
             ha='center',  # Horizontal alignment
             va='center',  # Vertical alignment
             fontsize=12,  # Font size
             bbox=dict(facecolor='yellow', alpha=0.5))  # Optional: add a background color to the annotation


        lines = [line1[0], line2[0]]
        labels = ['FFMC', 'FWI']

        ax.legend(lines, labels, loc='upper left')

        ax.set_xlabel('time')
        ax.set_ylabel('FFMC (%)')
        bx.set_ylabel('FWI (-)')


        fig.savefig('{:s}/ffmcfwi_{:s}_{:06d}.png'.format(diroutpdV,cexp,it))
        plt.close(fig)

