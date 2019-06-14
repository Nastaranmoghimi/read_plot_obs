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
import geopandas as gpd
import fiona


#live maps
import folium
import mplleaflet

#static maps
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import (LONGITUDE_FORMATTER,
                                   LATITUDE_FORMATTER)
import cartopy.feature as cfeature 
from matplotlib.offsetbox import AnchoredText


#us
lim = {}
lim['xmin']  = -125.0
lim['xmax']  = -55.
lim['ymin']  =  15.0
lim['ymax']  =  46.3

inp_dir = 'inp/'


def get_ind(lim,lons,lats):

    
    [ind] = np.where((lons > lim['xmin']) & 
                    ( lons < lim['xmax']) & 
                    ( lats > lim['ymin']) & 
                    ( lats < lim['ymax']))
    
    return ind


def make_map(projection=ccrs.PlateCarree()):                                                                                                                                        
                                                                                           
    """                                                                          
    Generate fig and ax using cartopy                                                                    
    input: projection                                                                                    
    output: fig and ax                             
    """                                  
    alpha = 0.5                                        
    subplot_kw = dict(projection=projection)                        
    fig, ax = plt.subplots(figsize=(9, 13),                           
                           subplot_kw=subplot_kw)   
    gl = ax.gridlines(draw_labels=True)                                 
    gl.xlabels_top = gl.ylabels_right = False 
    gl.xformatter = LONGITUDE_FORMATTER                        
    gl.yformatter = LATITUDE_FORMATTER                                
                                    
        # Put a background image on for nice sea rendering.             
    ax.stock_img()                                   
                                                          
    # Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
    states_provinces = cfeature.NaturalEarthFeature(                      
        category='cultural',                  
        name='admin_1_states_provinces_lines',
        scale='50m',           
        facecolor='none')        

    SOURCE = 'Natural Earth'
    LICENSE = 'public domain'
                                                                                                                                                                                    
    ax.add_feature(cfeature.LAND,zorder=0,alpha=alpha)          
    ax.add_feature(cfeature.COASTLINE,zorder=1,alpha=alpha)
    ax.add_feature(cfeature.BORDERS,zorder=1,alpha=2*alpha)
                       
    ax.add_feature(states_provinces, edgecolor='gray',zorder=1)
                                                          
    # Add a text annotation for the license information to the
    # the bottom right corner.                                            
    text = AnchoredText(r'$\mathcircled{{c}}$ {}; license: {}'
                        ''.format(SOURCE, LICENSE),
                        loc=4, prop={'size': 9}, frameon=True)                                    
    ax.add_artist(text)                                                                           
                                         
    ax.set_xlim(-132,-65)  #lon limits           
    ax.set_ylim( 20 , 55)  #lat limits   
    return fig, ax

#read purple exel file                                                                                                                                                              
purple_stations_info_file = inp_dir + 'PurpleAir_Locations.xlsx'
purp_infos = pd.read_excel(purple_stations_info_file)

print (purp_infos.head())
purp_lat = purp_infos['Lat'][:]
purp_lon = purp_infos['Lon'][:]
purp_lab = purp_infos['Label'][:]
purp_ind = get_ind(lim=lim,lons=purp_lon,lats=purp_lat)

data = dict(Name = purp_lab, lon = purp_lon, lat = purp_lat)
df_pm10 = pd.DataFrame(data=data)
df_pm10.to_csv  (purple_stations_info_file.replace('xlsx','csv'))

#read epa kml files
fiona.drvsupport.supported_drivers['kml'] = 'rw' # enable KML support which is disabled by default
fiona.drvsupport.supported_drivers['KML'] = 'rw' # enable KML support which is disabled by default

#epa 
epa_pm10_info_file = inp_dir + 'PM10.kml'
epa_pm10_locs = gpd.read_file(epa_pm10_info_file)
print (epa_pm10_locs.head())
 
n_epa_pm10 = len(epa_pm10_locs.centroid)
epa_pm10_lons = np.zeros((n_epa_pm10))
epa_pm10_lats = np.zeros((n_epa_pm10))

for i in range (n_epa_pm10):
    epa_pm10_lons[i] = epa_pm10_locs.centroid[i].x
    epa_pm10_lats[i] = epa_pm10_locs.centroid[i].y

epa_ind_pm10 = get_ind(lim=lim,lons=epa_pm10_lons,lats=epa_pm10_lats)

data = dict(Name = epa_pm10_locs['Name'], lon = epa_pm10_lons, lat = epa_pm10_lats)
df_pm10 = pd.DataFrame(data=data)
df_pm10.to_csv  (epa_pm10_info_file.replace('kml','csv'))

#epa 
epa_pm25_info_file = inp_dir + 'PM25.kml'
epa_pm25_locs = gpd.read_file(epa_pm25_info_file)
print (epa_pm25_locs.head())
 
n_epa_pm25 = len(epa_pm25_locs.centroid)
epa_pm25_lons = np.zeros((n_epa_pm25))
epa_pm25_lats = np.zeros((n_epa_pm25))

for i in range (n_epa_pm25):
    epa_pm25_lons[i] = epa_pm25_locs.centroid[i].x
    epa_pm25_lats[i] = epa_pm25_locs.centroid[i].y

epa_ind_pm25 = get_ind(lim=lim,lons=epa_pm10_lons,lats=epa_pm10_lats)


data = dict(Name = epa_pm25_locs['Name'], lon = epa_pm25_lons, lat = epa_pm25_lats)
df_pm25 = pd.DataFrame(data=data)
df_pm25.to_csv  (epa_pm25_info_file.replace('kml','csv'))



#static Cartopy map
fig,ax = make_map()                                             
ax.scatter(x=purp_infos['Lon'][:],y=purp_infos['Lat'][:],s=15,c='purple',lw=0,label = 'Purple' ,alpha=0.7)                                                                                           
ax.scatter(x=epa_pm10_lons       ,y=epa_pm10_lats       ,s=15,c='blue'  ,lw=0,label = 'EPA PM10',alpha=0.7)
ax.scatter(x=epa_pm25_lons       ,y=epa_pm25_lats       ,s=15,c='red'   ,lw=0,label = 'EPA PM25',alpha=0.7)
ax.legend()  
plt.savefig('test_map.png',dpi=450)
plt.close('all')



# live mplleaflet map
# Set up figure and axis
f2, ax2 = plt.subplots(1)
ax2.scatter(x=purp_infos['Lon'][purp_ind]  ,y=purp_infos['Lat'][purp_ind] ,c='purple',alpha=0.5,s=25)   
ax2.scatter(x=epa_pm10_lons[epa_ind]       ,y=epa_pm10_lats[epa_ind]       ,c='blue',alpha=0.5,s=25)   
ax2.scatter(x=epa_pm25_lons[epa_ind]       ,y=epa_pm25_lats[epa_ind]       ,c='red',alpha=0.5,s=25)   
ax2.set_xlim(-132,-65)  #lon limits 
ax2.set_ylim( 20 , 55)  #lat limits  
mapfile = 'test_mpleaflet.html'
mplleaflet.show(fig=f2,path=mapfile, tiles='mapbox bright')

mapfile = 'test_mpleaflet2.html'
mplleaflet.show(fig=f2,path=mapfile)

### nicer                                                                                                                                                                           
width, height = 1200, 800            
map = folium.Map(location=[42, -102],                                            
                 zoom_start=3, tiles='OpenStreetMap', width=width, height=height)                        
                                                                                                         
marker_cluster = folium.MarkerCluster().add_to(map)
print ('Map EPA pm10 obs ...'             )
for i in epa_ind_pm10:                                      
    folium.Marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],     
                      popup='EPA PM10 (' + epa_pm10_locs.Name[i]+')' ,
                      #icon = epa_ic                
                      icon=folium.Icon(color='blue', icon_color='white')
                      ).add_to(marker_cluster)
                                                               

for i in epa_ind_pm25:                                      
    folium.Marker(location=[epa_pm25_lats[i],epa_pm25_lons[i]],     
                      popup='EPA pm25 (' + epa_pm25_locs.Name[i]+')' ,
                      #icon = epa_ic                
                      icon=folium.Icon(color='blue', icon_color='white')
                      ).add_to(marker_cluster)


                                                                      
print  ('Map purple obs ...'          )
for i in purp_ind:                                                      
    folium.Marker(location=[purp_lat[i],purp_lon[i]],
                      popup='Purple (' + purp_lab[i]+')' ,
                      #icon = epa_ic
                      icon=folium.Icon(color='purple', icon_color='white')
                      ).add_to(marker_cluster)
                               
map.save('test_folium.html')   
                                 





"""


    ...: ### nicer                                                                                                                                                                           
    ...: width, height = 1200, 800            
    ...: map = folium.Map(location=[42, -102],                                            
    ...:                  zoom_start=3, tiles='OpenStreetMap', width=width, height=height)                        
    ...:                                                                                                          
    ...: marker_cluster = folium.MarkerCluster().add_to(map)
    ...: print 'Map EPA pm10 obs ...'             
    ...: for i in epa_ind:                                      
    ...:     folium.Marker(location=[epa_pm10_lats[i],epa_pm10_lons[i]],     
    ...:                       popup='EPA PM10 (' + epa_pm10_locs.Name[i]+')' ,
    ...:                       #icon = epa_ic                
    ...:                       icon=folium.Icon(color='blue', icon_color='white')
    ...:                       ).add_to(marker_cluster)
    ...:                                                                
    ...:                                                                       
    ...: print 'Map purple obs ...'          
    ...: for i in purp_ind:                                                      
    ...:     folium.Marker(location=[purp_lat[i],purp_lon[i]],
    ...:                       popup='Purple (' + purp_lab[i]+')' ,
    ...:                       #icon = epa_ic
    ...:                       icon=folium.Icon(color='purple', icon_color='white')
    ...:                       ).add_to(marker_cluster)
    ...:                                
    ...: map.save('test_folium.html')   
    ...:                                  

















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









"""
