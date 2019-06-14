#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Read and plot EPA and Purple site

"""
__author__ = "Nastaran Moghimi"
__copyright__ = "Copyright 2017, UCAR/NOAA"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "nastrann.moghimi@gmail.com"



import matplotlib.pyplot as plt
import numpy as np
import os,sys
import datetime
import string
import pandas as pd



import cartopy.crs as ccrs
from cartopy.mpl.gridliner import (LONGITUDE_FORMATTER,
                                   LATITUDE_FORMATTER)
import cartopy.feature as cfeature 


def make_map(projection=ccrs.PlateCarree(), bg='m'):
    
    """
    Generate fig and ax using cartopy
    input: projection
    output: fig and ax
    """

    subplot_kw = dict(projection=projection)
    fig, ax = plt.subplots(figsize=(9, 13),
                           subplot_kw=subplot_kw)
    gl = ax.gridlines(draw_labels=True)
    gl.xlabels_top = gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    
    ax.add_feature(cfeature.LAND) #If I comment this => all ok, but I need 
    ax.add_feature(cfeature.LAKES)
    ax.add_feature(cfeature.RIVERS)
    ax.coastlines()
    
    
    
    
    return fig, ax



width, height = 650, 500
map = folium.Map(location=[45.372, -121.6972],
                 zoom_start=3, tiles='OpenStreetMap', width=width, height=height)

marker_cluster = folium.MarkerCluster().add_to(map)

for i in range (n_epa_pm10)[::5]:
    print i
    #map.simple_marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
    #                  popup='EPA PM10',
    #                  marker_icon='glyphicon-plus')


    folium.Marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
                      popup='EPA PM10 \n  Name:' + epa_pm10_locs.Name[i] ,
                      icon = epa_ic,
                      ).add_to(marker_cluster)

map.save('test.html')
for i in range (n_epa_pm10)[::5]:
    print i
    #map.simple_marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
    #                  popup='EPA PM10',
    #                  marker_icon='glyphicon-plus')


    folium.Marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
                      popup='EPA PM10 \n  Name:' + epa_pm10_locs.Name[i] ,
                      #icon = epa_ic,
                      ).add_to(marker_cluster)









map.save('test.html')
width, height = 650, 500
map = folium.Map(location=[45.372, -121.6972],
                 zoom_start=3, tiles='OpenStreetMap', width=width, height=height)

marker_cluster = folium.MarkerCluster().add_to(map)

for i in range (n_epa_pm10)[::5]:
    print i
    #map.simple_marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
    #                  popup='EPA PM10',
    #                  marker_icon='glyphicon-plus')


    folium.Marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
                      popup='EPA PM10 \n  Name:' + epa_pm10_locs.Name[i] ,
                      #icon = epa_ic,
                      ).add_to(marker_cluster)









map.save('test.html')
epa_ic = folium.Icon(color=u'blue'  )
pur_ic = folium.Icon(color=u'purple')



width, height = 650, 500
map = folium.Map(location=[42, -102],
                 zoom_start=3, tiles='OpenStreetMap', width=width, height=height)

marker_cluster = folium.MarkerCluster().add_to(map)

for i in range (n_epa_pm10)[::5]:
    print i
    #map.simple_marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
    #                  popup='EPA PM10',
    #                  marker_icon='glyphicon-plus')


    folium.Marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
                      popup='EPA PM10 \n  Name:' + epa_pm10_locs.Name[i] ,
                      icon = epa_ic,
                      ).add_to(marker_cluster)









map.save('test.html')
width, height = 650, 500
map = folium.Map(location=[42, -102],
                 zoom_start=3, tiles='OpenStreetMap', width=width, height=height)

marker_cluster = folium.MarkerCluster().add_to(map)
    folium.Marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
                      popup='EPA PM10 \n  Name:' + epa_pm10_locs.Name[i] ,
                      icon = epa_ic,
                      ).add_to(marker_cluster)

    folium.Marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
                      popup='EPA PM10 \n  Name:' + epa_pm10_locs.Name[i] ,
                      icon = epa_ic,
                      ).add_to(marker_cluster)
 folium.Marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
                      popup='EPA PM10 \n  Name:' + epa_pm10_locs.Name[i] ,
                      icon = epa_ic,
                      ).add_to(marker_cluster)

map.save('test.html')

width, height = 650, 500
map = folium.Map(location=[42, -102],
                 zoom_start=3, tiles='OpenStreetMap', width=width, height=height)

marker_cluster = folium.MarkerCluster().add_to(map)

for i in range (n_epa_pm10)[::5]:
    print i
    #map.simple_marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
    #                  popup='EPA PM10',
    #                  marker_icon='glyphicon-plus')


    folium.Marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
                      popup='EPA PM10 \n  Name:' + epa_pm10_locs.Name[i] ,
                      icon = epa_ic   ).add_to(marker_cluster)









map.save('test.html')

epa_ic = folium.Icon(color=u'blue'  , icon_color=u'white', icon=u'stop', angle=0, prefix=u'glyphicon')
pur_ic = folium.Icon(color=u'purple', icon_color=u'white' , icon=u'play', angle=0, prefix=u'glyphicon')





width, height = 650, 500
map = folium.Map(location=[42, -102],
                 zoom_start=3, tiles='OpenStreetMap', width=width, height=height)

marker_cluster = folium.MarkerCluster().add_to(map)

for i in range (n_epa_pm10)[::5]:
    print i
    #map.simple_marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
    #                  popup='EPA PM10',
    #                  marker_icon='glyphicon-plus')


    folium.Marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
                      popup='EPA PM10 \n  Name:' + epa_pm10_locs.Name[i] ,
                      icon = epa_ic   ).add_to(marker_cluster)









map.save('test.html')
folium.Marker?
width, height = 650, 500
map = folium.Map(location=[42, -102],
                 zoom_start=3, tiles='OpenStreetMap', width=width, height=height)

marker_cluster = folium.MarkerCluster().add_to(map)

for i in range (n_epa_pm10)[::5]:
    print i
    #map.simple_marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
    #                  popup='EPA PM10',
    #                  marker_icon='glyphicon-plus')


    folium.Marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
                      popup='EPA PM10 \n  Name:' + epa_pm10_locs.Name[i] ,
                      #icon = epa_ic
                      icon=folium.Icon(color=df_counters["color"][point], icon_color='white', icon='male', angle=0, prefix='fa')
                      ).add_to(marker_cluster)

map.save('test.html')
width, height = 650, 500
map = folium.Map(location=[42, -102],
                 zoom_start=3, tiles='OpenStreetMap', width=width, height=height)

marker_cluster = folium.MarkerCluster().add_to(map)

for i in range (n_epa_pm10)[::5]:
    print i
    #map.simple_marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
    #                  popup='EPA PM10',
    #                  marker_icon='glyphicon-plus')


    folium.Marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
                      popup='EPA PM10 \n  Name:' + epa_pm10_locs.Name[i] ,
                      #icon = epa_ic
                      icon=folium.Icon(color='red', icon_color='white', icon='male', angle=0, prefix='fa')
                      ).add_to(marker_cluster)

map.save('test.html')
epa_ic = folium.Icon(color='red', icon_color='white', icon='male', angle=0, prefix='fa')
pur_ic = folium.Icon(color='purple', icon_color='white', icon='play', angle=0, prefix='glyphicon')




width, height = 650, 500
map = folium.Map(location=[42, -102],
                 zoom_start=3, tiles='OpenStreetMap', width=width, height=height)

marker_cluster = folium.MarkerCluster().add_to(map)

for i in range (n_epa_pm10)[::5]:
    print i
    #map.simple_marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
    #                  popup='EPA PM10',
    #                  marker_icon='glyphicon-plus')


    folium.Marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
                      popup='EPA PM10 \n  Name:' + epa_pm10_locs.Name[i] ,
                      icon = epa_ic
                      #icon=folium.Icon(color='red', icon_color='white', icon='male', angle=0, prefix='fa')
                      ).add_to(marker_cluster)

map.save('test.html')

width, height = 650, 500
map = folium.Map(location=[42, -102],
                 zoom_start=3, tiles='OpenStreetMap', width=width, height=height)

marker_cluster = folium.MarkerCluster().add_to(map)

for i in range (n_epa_pm10)[::5]:
    print i
    #map.simple_marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
    #                  popup='EPA PM10',
    #                  marker_icon='glyphicon-plus')


    folium.Marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
                      popup='EPA PM10 \n  Name:' + epa_pm10_locs.Name[i] ,
                      #icon = epa_ic
                      icon=folium.Icon(color='red', icon_color='white', icon='male', angle=0, prefix='fa')
                      ).add_to(marker_cluster)

map.save('test.html')

width, height = 650, 500
map = folium.Map(location=[42, -102],
                 zoom_start=3, tiles='OpenStreetMap', width=width, height=height)

marker_cluster = folium.MarkerCluster().add_to(map)

for i in range (n_epa_pm10)[::5]:
    print i
    #map.simple_marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
    #                  popup='EPA PM10',
    #                  marker_icon='glyphicon-plus')


    folium.Marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
                      popup='EPA PM10 \n  Name:' + epa_pm10_locs.Name[i] ,
                      #icon = epa_ic
                      icon=folium.Icon(color='purple', icon_color='white')
                      ).add_to(marker_cluster)

map.save('test.html')
width, height = 650, 500
map = folium.Map(location=[42, -102],
                 zoom_start=3, tiles='OpenStreetMap', width=width, height=height)

marker_cluster = folium.MarkerCluster().add_to(map)

for i in range (n_epa_pm10)[::5]:
    print i
    #map.simple_marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
    #                  popup='EPA PM10',
    #                  marker_icon='glyphicon-plus')


    folium.Marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
                      popup='EPA PM10 (' + epa_pm10_locs.Name[i]+')' ,
                      #icon = epa_ic
                      icon=folium.Icon(color='purple', icon_color='white')
                      ).add_to(marker_cluster)





map.save('test.html')
width, height = 650, 500
map = folium.Map(location=[42, -102],
                 zoom_start=3, tiles='OpenStreetMap')#, width=width, height=height)

marker_cluster = folium.MarkerCluster().add_to(map)

for i in range (n_epa_pm10)[::5]:
    print i
    #map.simple_marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
    #                  popup='EPA PM10',
    #                  marker_icon='glyphicon-plus')


    folium.Marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
                      popup='EPA PM10 (' + epa_pm10_locs.Name[i]+')' ,
                      #icon = epa_ic
                      icon=folium.Icon(color='purple', icon_color='white')
                      ).add_to(marker_cluster)





map.save('test.html')
purp_infos.keys()
purp_infos.Label
purp_infos.Labelpurp_lab = purp_infos['Label'][:]
purp_lab = purp_infos['Label'][:]
purp_lab
def get_ind(lons,lats):
    lim = {}
    lim['xmin']  = -125.0
    lim['xmax']  = -55.
    lim['ymin']  =  15.0
    lim['ymax']  =  46.3
    
    [ind] = np.where((lons > lim['xmin']) & 
                    ( lons < lim['xmax']) & 
                    ( lats > lim['ymin']) & 
                    ( lats < lim['ymax']))
    
    return ind
purp_lat = purp_infos['Lat'][:]
purp_lon = purp_infos['Lon'][:]
purp_lab = purp_infos['Label'][:]

purp_ind = get_ind(lons=purp_lon,lats=purp_lat)

purp_lat = purp_lat[purp_ind]
purp_lon = purp_lon[purp_ind]
purp_lab = purp_lab[purp_ind]
purp_lat
purp_lat.shape()
len(purp_lat)

#us
lim = {}
lim['xmin']  = -125.0
lim['xmax']  = -55.
lim['ymin']  =  15.0
lim['ymax']  =  46.3




def get_ind(lim,lons,lats):

    
    [ind] = np.where((lons > lim['xmin']) & 
                    ( lons < lim['xmax']) & 
                    ( lats > lim['ymin']) & 
                    ( lats < lim['ymax']))
    
    return ind
purp_lat = purp_infos['Lat'][:]
purp_lon = purp_infos['Lon'][:]
purp_lab = purp_infos['Label'][:]

purp_ind = get_ind(lim=lim,lons=purp_lon,lats=purp_lat)

purp_lat = purp_lat[purp_ind]
purp_lon = purp_lon[purp_ind]
purp_lab = purp_lab[purp_ind]
epa_ind = get_ind(lim=lim,lons=epa_pm10_lons,lats=epa_pm10_lats)
epa_ind
purp_lat = purp_infos['Lat'][:]
purp_lon = purp_infos['Lon'][:]
purp_lab = purp_infos['Label'][:]

purp_ind = get_ind(lim=lim,lons=purp_lon,lats=purp_lat)

width, height = 1200, 800
map = folium.Map(location=[42, -102],
                 zoom_start=3, tiles='OpenStreetMap', width=width, height=height)

marker_cluster = folium.MarkerCluster().add_to(map)

for i in epa_ind:
    print i
    #map.simple_marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
    #                  popup='EPA PM10',
    #                  marker_icon='glyphicon-plus')


    folium.Marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
                      popup='EPA PM10 (' + epa_pm10_locs.Name[i]+')' ,
                      #icon = epa_ic
                      icon=folium.Icon(color='blue', icon_color='white')
                      ).add_to(marker_cluster)




purp_lat = purp_infos['Lat'][:]
purp_lon = purp_infos['Lon'][:]
purp_lab = purp_infos['Label'][:]

for i in purp_ind:
    print i
    #map.simple_marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
    #                  popup='EPA PM10',
    #                  marker_icon='glyphicon-plus')


    folium.Marker(location=[purp_lat[i],purp_lon[i]],
                      popup='Purple (' + purp_lab[i]+')' ,
                      #icon = epa_ic
                      icon=folium.Icon(color='purple', icon_color='white')
                      ).add_to(marker_cluster)





map.save('test.html')
map.save('test.html')
    import mplleaflet
    #plot
    # Set up figure and axis
    f2, ax2 = plt.subplots(1)
    ax2.scatter(x=purp_infos['Lon'][purp_ind]  ,y=purp_infos['Lat'][purp_ind] ,c='purple')
    ax2.scatter(x=epa_pm10_lons[epa_ind]       ,y=epa_pm10_lats[epa_ind]       ,c='blue')
    ax2.set_xlim(-132,-65)  #lon limits 
    ax2.set_ylim( 20 , 55)  #lat limits  
    mapfile = 'test_mpleaflet.html'
    mplleaflet.show(fig=f2,path=mapfile, tiles='mapbox bright')
    import mplleaflet
    #plot
    # Set up figure and axis
    f2, ax2 = plt.subplots(1)
    ax2.scatter(x=purp_infos['Lon'][purp_ind]  ,y=purp_infos['Lat'][purp_ind] ,c='purple',s=10)
    ax2.scatter(x=epa_pm10_lons[epa_ind]       ,y=epa_pm10_lats[epa_ind]       ,c='blue' ,s=10)
    ax2.set_xlim(-132,-65)  #lon limits 
    ax2.set_ylim( 20 , 55)  #lat limits  
    mapfile = 'test_mpleaflet.html'
    mplleaflet.show(fig=f2,path=mapfile, tiles='mapbox bright')
    import mplleaflet
    #plot
    # Set up figure and axis
    f2, ax2 = plt.subplots(1)
    ax2.scatter(x=purp_infos['Lon'][purp_ind]  ,y=purp_infos['Lat'][purp_ind] ,c='purple',s=50,alpha=0.7)
    ax2.scatter(x=epa_pm10_lons[epa_ind]       ,y=epa_pm10_lats[epa_ind]       ,c='blue' ,s=50,alpha=0.7)
    ax2.set_xlim(-132,-65)  #lon limits 
    ax2.set_ylim( 20 , 55)  #lat limits  
    mapfile = 'test_mpleaflet.html'
    mplleaflet.show(fig=f2,path=mapfile, tiles='mapbox bright')

### nicer
width, height = 1200, 800
map = folium.Map(location=[42, -102],
                 zoom_start=3, tiles='OpenStreetMap', width=width, height=height)

marker_cluster = folium.MarkerCluster().add_to(map)
print 'Map EPA pm10 obs ...'
for i in epa_ind:
    folium.Marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],
                      popup='EPA PM10 (' + epa_pm10_locs.Name[i]+')' ,
                      #icon = epa_ic
                      icon=folium.Icon(color='blue', icon_color='white')
                      ).add_to(marker_cluster)


print 'Map purple obs ...'
for i in purp_ind:
    folium.Marker(location=[purp_lat[i],purp_lon[i]],
                      popup='Purple (' + purp_lab[i]+')' ,
                      #icon = epa_ic
                      icon=folium.Icon(color='purple', icon_color='white')
                      ).add_to(marker_cluster)

map.save('test_folium.html')









