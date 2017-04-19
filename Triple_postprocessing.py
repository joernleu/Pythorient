#!/usr/bin/env python

#!/usr/bin/env python
from matplotlib import pyplot as plt
from matplotlib import cm as cm
from matplotlib import mlab as ml
import numpy as np
import os

# my_funcs
from IO.read_data_ang_file import read_data_hex_ang,det_Step_Size,read_GB_data,read_TJ_data
from IO.write_EBSD_data import write_GB_stats, write_TJ_stats
from OIM_Mapping.get_verts import get_verts
from OIM_Mapping.plot_hex_map import plot_hex_map,plot_hex_rgb_map, plot_line_map
from OIM_Mapping.point_densities import grid_density_gaussian_filter,grid_density_boxsum
#path = 'EBSD_Data/5_Low_res_200_hexagonal_cells_rand_twin_distr_rad/'
#my_file =path + 'test_grain.ang'
#path = 'EBSD_Data/1_E6c_150_sample_holder/'
#my_file =path + '1_E6c_150_sample_holder.txt'
#path = 'EBSD_Data/52_Enano4_1rot/Center_plus_500mu/'
#my_file =path + 'Center_1500mu.txt'

path = 'EBSD_Data/52_Enano4_1rot/Center_plus_500mu/'
my_file =path + 'Center_500mu.txt'
path = 'EBSD_Data/52_Enano4_1rot/Center/'
my_file =path + 'Center.txt'
#path = 'EBSD_Data/HPTCu_150_electro/'
#my_file =path +'150_electro.txt'
#path = 'EBSD_Data/50_Enano3/Enano3_Center_plus_2000mu/'
#my_file =path + 'Enano3_Center_plus_2000mu.txt'
#my_file =path + '1_E6c_150_sample_holder.txt'
headerlines, phi1, Phi, phi2, x, y, IQ, CI, Fit, Grain_ID, edge, phase = read_data_hex_ang(my_file)
XSTEP,YSTEP, NCOLS_ODD, NCOLS_EVEN, NROWS = det_Step_Size(x,y)

verts = get_verts(x,y,YSTEP)


## Determine and plot general and CSL boundaries
#-------------------------------------------------------------
Sigma_boundary = ['General','3','9']
GB_color = ['k','y','r']
GB_linewidth =1*YSTEP**0.5
Triple_junction = ['General','3','3-3','3-3-9']
TJ_color =['b','g','y','r']
TJ_marker =['^','>','<','v']
TJ_markersize = 5*YSTEP**0.5

GB_data = []
for i in Sigma_boundary:
    GB_data.append(read_GB_data(path+'trip_out/'+i+'_GBs.txt'))

TJ_data = []
for i in Triple_junction:
    TJ_data.append(read_TJ_data(path+'trip_out/'+i+'_TJs.txt'))

## Plotting
#-------------------------------------------------------------

## IPF map
#------------------------
#IPF_r_g_b = []
#for i in range(len(phi1)):
#    r,g,b = get_IPF_color_vals([phi1[i],Phi[i],phi2[i]])
#    IPF_r_g_b.append([r,g,b])
#fig = plot_hex_rgb_map(x,y,XSTEP, YSTEP, np.array(IPF_r_g_b),NCOLS_EVEN,NCOLS_ODD,NROWS, title=str('IPF map'))

## IQ or other color coded map
#------------------------------
fig = plt.figure()
view_ymin = 0
view_ymax = max(y)
view_xmin = 0
view_xmax = max(x)

fig.set_size_inches(view_xmax/view_ymax*5, 5, forward=False)
ax = plt.Axes(fig, [0., 0., 1, 1])
ax.set_axis_off()
fig.add_axes(ax)
ax.set_xlim([0, view_xmax])
ax.set_ylim([0, view_ymax])


#IQ map
#alpha_val = 0.5
#plot_hex_map(fig,ax,alpha_val,IQ,verts, title=str(IQ))
x_TJ = []
y_TJ = []
for i in range(len(TJ_data)):

    if Triple_junction[i] in ['3']:#,'3-3','3-3-9']:
        for n in range(len(TJ_data[i][0])):
            x_TJ.append(verts[TJ_data[i][0][n]][TJ_data[i][2][n],0])
            y_TJ.append(verts[TJ_data[i][0][n]][TJ_data[i][2][n],1])
            
## view area range
density = 512
zd = grid_density_gaussian_filter(view_xmin, view_ymin, view_xmax, view_ymax, density, density, zip(x_TJ, y_TJ))
##zd = grid_density_boxsum(view_xmin, view_ymin, view_xmax, view_ymax, density, density, zip(x_g, y_g))


plt.imshow(zd ,cmap=cm.gray, origin='lower', extent=[view_xmin, view_xmax, view_ymin, view_ymax])
plt.colorbar(shrink=.5)#, pad=0.04)
#cbaxes = fig.add_axes([1.0, 0.1, 0.03, 0.8]) 
#cb = plt.colorbar( cax = cbaxes)  
#plt.savefig(path[m] + 'trip_out/Triple_junction_density.png',bbox_inches=0,pad_inches=0, dpi=900,transparent=True)
legend_names_GB = ['General GB','$\Sigma$3 GB','$\Sigma$9 GB']
for i in range(len(GB_data)):
    xy_GB = []
    for n in range(len(GB_data[i][0])):
        xy_GB.append(verts[GB_data[i][0][n]][GB_data[i][2][n],:])
        #y_GB.append(verts[GB_data[i][0][n]][GB_data[i][2][n],1])
    print(i)    
    plot_line_map(xy_GB,ax,GB_color[i],GB_linewidth,legend_names_GB[i] )

legend_names_TJ = ['General','$\Sigma$3-X-X TJ','$\Sigma$3-$\Sigma$3-X TJ','$\Sigma$3-$\Sigma$3-$\Sigma$9 TJ']
#print(x_TJ)
for i in range(len(TJ_data)):
    if Triple_junction[i] in ['3','3-3','3-3-9']:
        for n in range(len(TJ_data[i][0])):
            plt.plot(verts[TJ_data[i][0][n]][TJ_data[i][2][n],0],verts[TJ_data[i][0][n]][TJ_data[i][2][n],1],TJ_color[i]+TJ_marker[i],markeredgecolor= TJ_color[i], markersize=TJ_markersize)
        plt.plot([],[],TJ_color[i]+TJ_marker[i],markeredgecolor= TJ_color[i], markersize=TJ_markersize,label = legend_names_TJ[i])       
#else:
#    #    for p in range(len(TJ_data[i])):
#    #        for n in range(len(TJ_data[i][p][0])):
#    #          plt.plot(verts[TJ_data[i][p][0][n]][TJ_data[i][p][2][n],0],verts[TJ_data[i][p][0][n]][TJ_data[i][p][2][n],1],TJ_color[i+p]+TJ_marker[i+p],markeredgecolor= 'k', markersize=15*YSTEP)
#
#

plt.legend(loc = 1,framealpha = 1, numpoints = 1,prop={'size':6})# ,legend_names,loc=1,prop={'size':10})
plt.gca().invert_yaxis()
plt.savefig(path + 'trip_out/Triple_junction_density.png',bbox_inches=0,pad_inches=0, dpi=900)


#plt.show()

