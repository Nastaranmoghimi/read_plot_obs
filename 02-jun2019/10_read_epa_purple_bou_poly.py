#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Read and plot EPA and Purple site

"""
__author__ = "Nastaran Moghimi"
__copyright__ = "Copyright 2017, UCAR/NOAA"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "nastarann.moghimi@gmail.com"


# Updated
# Thu 13 Jun 2019 10:12:35 PM EDT   read new json file and find stations in Cali
#
#




import matplotlib.pyplot as plt
import numpy as np
import os,sys
import datetime
import string
import pandas as pd
#import geopandas as gpd
#import fiona

from    collections import defaultdict

#live maps
import folium
#import mplleaflet

#static maps
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import (LONGITUDE_FORMATTER,
                                   LATITUDE_FORMATTER)
import cartopy.feature as cfeature 
from matplotlib.offsetbox import AnchoredText

import pl_tools

from matplotlib.path import Path

###functions
def find_purple():

    #### JSON ######
    json_file = inp_dir + 'purpleair_13june2019.json'
    purp_jsonall = pd.read_json(json_file)

    purp = defaultdict(dict)
    for il in range (len( purp_jsonall)):
        res      = purp_jsonall.iloc[il]['results']
        try:
            purp [res['ID']]['lat'  ]     = float(res['Lat'])
            purp [res['ID']]['lon'  ]     = float(res['Lon'])
            purp [res['ID']]['pm2_5']     = res['PM2_5Value']
            purp [res['ID']]['last_seen'] = datetime.datetime(1970,1,1) + datetime.timedelta(seconds = res['LastSeen'] )
            purp [res['ID']]['label'    ] = res['Label']
            purp [res['ID']]['id'    ]    = res['ID']
        except:
            print (il, '===None')

    purp_df = pd.DataFrame.from_dict(purp)
    #purp_df.dropna(inplace=True)

    return purp_df




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
        scale='10m',           
        facecolor='none')        

    #coast            = cfeature.NaturalEarthFeature(
    #    category='physical', 
    #    name='coastline',
    #    scale='10m',
    #    facecolor='none') 
    
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




def wait():
    """
    
    """
    while True:
        choice = input("Enter 1 when you are done .. > ")
        if choice == 1 :
            break


def create_boundary(name = '', region = {}):
    """
    
    """
    lim = region[name]
    filename = '../inp/'+name + '_bou.txt'

    if (not os.path.exists(filename)):
        fig,ax = make_map()       
        ax.set_xlim(lim['xmin'],lim['xmax'])  
        ax.set_ylim(lim['ymin'],lim['ymax'])

        poly = pl_tools.InteractiveLine(type='cblin')
        plt.show()
        wait()
    
        # Open filename
        f = open(filename,'w')
        for i in range(len(poly.x)):
            f.write( str(poly.x[i])+ '  ' +str(poly.y[i]) + '  \n' )
        
        f.close()    
    
        lons = poly.x
        lats = poly.y
    else:
        bou = np.loadtxt(filename)
        lons = bou[:,0]
        lats = bou[:,1] 

    return lons,lats 



def get_ind_poly(lons_bound,lats_bound,lons,lats):
    """
    Return index of points inside a polygon
    
    Input:
    bou_lons,bou_lats:  
    lons and lats: data poits coordinates
    
    Output:
    ind: index of points inside rectangle
    
    
    """

    vertices = np.c_[lons_bound,lats_bound]
    p    = Path(vertices)

    data = np.c_[lons, lats]
    mask = p.contains_points(data)

    ind = np.where(mask==True)
    return ind[0]


inp_dir = '../inp/'


#####
region = defaultdict(dict)
#####
#name = 'us'
#region[name]['xmin']  = -125.0
#region[name]['xmax']  = -55.
#region[name]['ymin']  =  15.0
#region[name]['ymax']  =  46.3

name = 'california'
region[name]['xmin']  = -125.2
region[name]['xmax']  = -112.8
region[name]['ymin']  =  31.
region[name]['ymax']  =  43.1



#### MAIN
bou_type = 'poly'

if bou_type == 'poly':
    lons_bou,lats_bou = create_boundary(name = name, region = region)



purp = find_purple()


####
lons = []
lats = []
ids  = []
lbs  = []

for key in purp.keys():
    #get purple located outside 
    try:
        
        if not purp [key]['label'].endswith(' B'):
            lons.append (purp [key]['lon']  )
            lats.append (purp [key]['lat']  )
            lbs.append  (purp [key]['label'])
            ids.append  (purp [key]['id']   )
    except:
        print (key)
        pass

lons  = np.array(lons)
lats  = np.array(lats)
lbs   = np.array(lbs)
ids   = np.array(ids)

purp_ind = get_ind_poly(lons_bound = lons_bou,lats_bound=lats_bou,lons =lons ,lats=lats)


ids_inside   = ids   [purp_ind]
lons_inside  = lons  [purp_ind]
lats_inside  = lats  [purp_ind]
lbs_inside   = lbs   [purp_ind]
#
data_sta = np.c_[ids_inside,lons_inside,lats_inside,lbs_inside]

df = pd.DataFrame(data=data_sta,columns=['id','lon','lat','label'])
df = df.dropna()
df.to_csv(name + '_purple_air_june2019.csv', sep='\t', encoding = 'utf-8',index=False)



lim =region[name]

#plot stations
print ('Static Cartopy map ...')
fig,ax = make_map()                                             
ax.set_title( name.capitalize() + ' Stations')
ax.scatter(x = lons_inside , y =lats_inside ,s=20,lw=0, c= 'purple',alpha=0.85,label = 'Purple air')
ax.legend(ncol=4)
ax.set_xlim(lim['xmin'],lim['xmax'])
ax.set_ylim(lim['ymin'],lim['ymax'])


plt.savefig(name +'_all.png',dpi=450)
plt.show()
#plt.close('all')
#######################################


































