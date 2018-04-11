#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Read Purple site specs

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



#### MAIN
#us
lim = {}
lim['xmin']  = -125.0
lim['xmax']  = -55.
lim['ymin']  =  15.0
lim['ymax']  =  46.3




inp_dir = 'inp/'


print 'Read ...'
#read purple exel file                                                                                                                                                              
purple_stations_info_file = inp_dir + 'PurpleAir_Locations.xlsx'
purp_infos = pd.read_excel(purple_stations_info_file)

#read EPA locations from csv                                                                                                                                                         
epa_pm10_info_file = inp_dir + 'PM10.csv'
epa_pm10_locs      = pd.read_csv(epa_pm10_info_file)

epa_pm25_info_file = inp_dir + 'PM25.csv'
epa_pm25_locs      = pd.read_csv(epa_pm25_info_file)








print 'get index ...'
#get index of obs sites inside epecific window LIM
purp_ind     = get_ind(lim=lim,lons=purp_infos   ['Lon'][:],lats=purp_infos   ['Lat'][:])
epa_ind_pm10 = get_ind(lim=lim,lons=epa_pm10_locs['lon'][:],lats=epa_pm10_locs['lat'][:])
epa_ind_pm25 = get_ind(lim=lim,lons=epa_pm25_locs['lon'][:],lats=epa_pm25_locs['lat'][:])


#plot stations
print 'Static Cartopy map ...'
fig,ax = make_map()                                             
ax.scatter(x=purp_infos['Lon'][purp_ind]          ,y=purp_infos['Lat'][purp_ind]          ,s=10,c='purple',lw=0,label = 'Purple' ,alpha=0.7)                                                                                           
ax.scatter(x=epa_pm10_locs['lon'][epa_ind_pm10]  ,y=epa_pm10_locs['lat'] [epa_ind_pm10]   ,s=10,c='blue'  ,lw=0,label = 'EPA PM10',alpha=0.7)
ax.scatter(x=epa_pm25_locs['lon'][epa_ind_pm25]  ,y=epa_pm25_locs['lat'] [epa_ind_pm25]   ,s=10,c='red'   ,lw=0,label = 'EPA PM25',alpha=0.7)
ax.legend()  
ax.set_xlim(lim['xmin'],lim['xmax'])
ax.set_ylim(lim['ymin'],lim['ymax'])

plt.savefig('test_map.png',dpi=450)
plt.show()
#plt.close('all')





if False:
    print ' Plot mplleaflet map ...'
    # live mplleaflet map
    # Set up figure and axis
    f2, ax2 = plt.subplots(1)
    ax2.scatter(x=purp_infos['Lon'][purp_ind]         ,y=purp_infos   ['Lat'][purp_ind]     ,c='purple',alpha=0.5,s=25)   
    ax2.scatter(x=epa_pm10_locs['lon'][epa_ind_pm10]  ,y=epa_pm10_locs['lat'][epa_ind_pm10] ,c='blue'  ,alpha=0.5,s=25)   
    ax2.scatter(x=epa_pm25_locs['lon'][epa_ind_pm25]  ,y=epa_pm25_locs['lat'][epa_ind_pm25] ,c='red'   ,alpha=0.5,s=25)   
    ax2.set_xlim(-132,-65)  #lon limits 
    ax2.set_ylim( 20 , 55)  #lat limits  
    #mapfile = 'test_mpleaflet.html'
    #mplleaflet.show(fig=f2,path=mapfile, tiles='mapbox bright')

    mapfile = 'test_mpleaflet2.html'
    mplleaflet.show(fig=f2,path=mapfile)

    ### nicer       
    print ' Plot folium map ...'                                                                                                                                                                    
    width, height = 1200, 800            
    map = folium.Map(location=[42, -102],                                            
                     zoom_start=3, tiles='OpenStreetMap', width=width, height=height)                        
                                                                                                             
    marker_cluster = folium.MarkerCluster().add_to(map)
    print ('Map EPA pm10 obs ...'             )
    for i in epa_ind_pm10:                                      
        folium.Marker(location=[epa_pm10_locs['lat'][i],epa_pm10_locs['lon'][i]],     
                          popup='EPA PM10 (' + str(epa_pm10_locs.Name[i])+')' ,
                          #icon = epa_ic                
                          icon=folium.Icon(color='blue', icon_color='white')
                          ).add_to(marker_cluster)
                                                                   

    for i in epa_ind_pm25:                                      
        folium.Marker(location=[epa_pm25_locs['lat'][i],epa_pm25_locs['lon'][i]],     
                          popup='EPA pm25 (' + str(epa_pm25_locs.Name[i])+')' ,
                          #icon = epa_ic                
                          icon=folium.Icon(color='blue', icon_color='white')
                          ).add_to(marker_cluster)

    print  ('Map purple obs ...'          )
    for i in purp_ind:                                                      
        folium.Marker(location=[purp_infos   ['Lat'][i],purp_infos   ['Lon'][i]],
                          popup='Purple (' + purp_infos['Label'][i]+')' ,
                          #icon = epa_ic
                          icon=folium.Icon(color='purple', icon_color='white')
                          ).add_to(marker_cluster)
                                   
    map.save('test_folium.html')   
                                     

