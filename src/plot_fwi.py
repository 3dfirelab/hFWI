import xarray as xr 
import matplotlib.pyplot as plt 
import sys
import os 
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import rioxarray
import numpy as np 
from pyproj import Transformer

###########################
if __name__ == '__main__':
###########################

    dirin = '/data/paugam/MNH/Cat_PdV/006_mnhSolo'
    diroutfwi = dirin+'/Postproc/FWI/PNG-fwi/'
    diroutffmc = dirin+'/Postproc/FWI/PNG-ffmc/'
    diroutpdV = dirin+'/Postproc/FWI/PNG-pdV/'
    flag_model = 'mnh'
    cexp = 'FCAST'
    os.makedirs(diroutfwi, exist_ok=True)
    os.makedirs(diroutffmc, exist_ok=True)
    os.makedirs(diroutpdV, exist_ok=True)
    crshere = 32631

    #From Wikipedia, the free encyclopedia
    pdV_lat=41.702
    pdV_lon=1.873
    
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:32631", always_xy=True)
    pdV_easting, pdV_northing = transformer.transform(pdV_lon, pdV_lat)

    ds = xr.open_dataset(dirin+'/Postproc/FWI/{:s}_fwi.nc'.format(cexp))
    ds = ds.rio.write_crs(4326)
    ds = ds.rio.reproject(crshere)
    
    projection = ccrs.epsg(ds.rio.crs.to_epsg())
    
    ffmcmin=np.nanpercentile(ds.FFMC,5)
    ffmcmax=np.nanpercentile(ds.FFMC,95)
    fwimin=np.nanpercentile(ds.FWI,20)
    fwimax=np.nanpercentile(ds.FWI,80)
    
    tab10 = plt.get_cmap("tab10")

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

        ds.FWI.isel(time=it).plot(ax=ax,cmap='Reds',vmin=fwimin,vmax=fwimax)

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

