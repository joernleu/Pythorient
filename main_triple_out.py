#!/usr/bin/env python
import numpy as np
import os
from matplotlib import pyplot as plt
# my_funcs
from IO.read_data_ang_file import read_data_hex_ang,det_Step_Size
from IO.write_EBSD_data import write_GB_data, write_TJ_data, write_GB_stats,write_TJ_stats
from OIM_Mapping.get_verts import get_verts
from OIM_Mapping.identify_GBs import identify_GBs,identify_Sigma_GBs
from OIM_Mapping.identify_TJs import identify_TJs,identify_Sigma_TJs
from OIM_Mapping.plot_hex_map import plot_hex_map,plot_hex_rgb_map, plot_line_map
from OIM_Mapping.point_densities import grid_density_gaussian_filter,grid_density_boxsum
## Import .ang file with grain ID
#-------------------------------------------------------------

path = ['Test_Files/200_hexagonal_cells_rand_twin_distr/']

ang_file =['test_grain.ang']
    
for m in range(len(path)):    
    my_file =path[m] + ang_file[m]
    
    if not os.path.exists(path[m]+'triple_out/'):
        os.makedirs(path[m]+'triple_out/')
    
    
    headerlines, phi1, Phi, phi2, x, y, IQ, CI, Fit, Grain_ID, edge, phase = read_data_hex_ang(my_file)
    XSTEP,YSTEP, NCOLS_ODD, NCOLS_EVEN, NROWS = det_Step_Size(x,y)
    
    crystal_sys = 'Cubic'
    Bravais_lattice = 'fcc'
    verts = get_verts(x,y,YSTEP)
    
    
    ## Determine and plot general and CSL boundaries
    #-------------------------------------------------------------
    Sigma_boundary = ['3','9']
    GB_color = ['y','r']
    GB_linewidth = 1 * YSTEP**0.5
    Triple_junction = [['3', '3', '9']]
    TJ_color = ['g', 'y', 'r']
    TJ_marker = ['>', '<', 'v']
    TJ_markersize = 5 * YSTEP**0.5
    
    
    Data_GB = {}
    ## General boundaries 
    #---------------------------------------
    identify_GBs(Data_GB, Grain_ID, NROWS, NCOLS_EVEN, NCOLS_ODD)
    #print(General_GB_data[0])
    
    ## CSL Sigma boundaries
    ##---------------------------------------
    
    for i in Sigma_boundary:
        identify_Sigma_GBs(Data_GB,
                            i,
                            'Brandon',
                            crystal_sys,
                            Bravais_lattice,
                            phi1, Phi, phi2)
    
    
    ## Determine and plot general and CSL triple Junctions
    #-------------------------------------------------------------
    Data_TJ = {}
    ## General triple junctions
    #---------------------------------------
    identify_TJs(Data_GB, Data_TJ, Grain_ID, NCOLS_ODD, NCOLS_EVEN)
    ## CSL Sigma triple junctions
    #---------------------------------------
    for i in Triple_junction:
        identify_Sigma_TJs(Data_GB,
                            Data_TJ,
                            i,
                            'Brandon',
                            crystal_sys,
                            Bravais_lattice,
                            phi1, Phi, phi2)

    
    for px in Data_GB:
            print(Data_GB.keys())
    for px in Data_TJ:
            print(Data_TJ.keys())
    ### Export
    no_General_GBs = len(Data_GB['General_GB'][0])
    write_GB_data(Data_GB, 'General_GB', path[m])
    write_GB_stats(Data_GB, no_General_GBs,'General_GB', path[m])
    
    for i in range(len(Sigma_boundary)):
        write_GB_data(Data_GB,Sigma_boundary[i], path[m])
        write_GB_stats(Data_GB,no_General_GBs,Sigma_boundary[i], path[m])
    
    no_General_TJs = len(Data_TJ['General_TJ'][0])
    write_TJ_data(Data_TJ, 'General_TJ', path[m])
    write_TJ_stats(Data_TJ, no_General_TJs, 'General_TJ', path[m]) 
    
    for i in range(len(Triple_junction)):
        for n in range(len(Triple_junction[i])):        
            Triple = Triple_junction[i][0:n+1]
            Triple =''.join(Triple)
            write_TJ_data(Data_TJ, Triple, path[m])
            write_TJ_stats(Data_TJ, no_General_TJs, Triple,path[m])
    

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
    plt.gca().invert_yaxis()
    # plot only 3-3-G and 3-3-9 Triple junctions
    #IQ map
    alpha_val = 0.5
    #plot_hex_map(fig,ax,alpha_val,IQ,verts, title=str(IQ))
    x_TJ = []
    y_TJ = []
    for i in ['3']:
        for n in range(len(Data_TJ[i][0])):
             x_TJ.append(verts[Data_TJ[i][0][n]][Data_TJ[i][2][n],0])
             y_TJ.append(verts[Data_TJ[i][0][n]][Data_TJ[i][2][n],1])
                        
    ## view area range
    density = 512
    zd = grid_density_gaussian_filter(view_xmin, view_ymin, view_xmax, view_ymax, density, density, zip(x_TJ, y_TJ))
    ##zd = grid_density_boxsum(view_xmin, view_ymin, view_xmax, view_ymax, density, density, zip(x_g, y_g))
    plt.imshow(zd , origin='lower', extent=[view_xmin, view_xmax, view_ymin, view_ymax])
    #
    plt.savefig(path[m] + 'triple_out/Triple_junction_density.png',bbox_inches=0,pad_inches=0, dpi=900,transparent=True)
    
    # plot general GBs
    xy_GB = []
    for n in range(len(Data_GB['General_GB'][0])):
        xy_GB.append(verts[Data_GB['General_GB'][0][n]][Data_GB['General_GB'][2][n],:])
    plot_line_map(xy_GB, ax,'k',GB_linewidth)
    
    
    legend_names_GB = ['General GB','$\Sigma$3 GB','$\Sigma$9 GB']
    col = 0
    for i in Sigma_boundary:
        xy_GB = []
        for n in range(len(Data_GB[i][0])):
            xy_GB.append(verts[Data_GB[i][0][n]][Data_GB[i][2][n],:])
        plot_line_map(xy_GB,ax,GB_color[col],GB_linewidth)
        col += 1


    #for i in range(len(General_TJ_data[0])):
    #    plt.plot(verts[General_TJ_data[0][i]][General_TJ_data[2][i],0],
    #            verts[General_TJ_data[0][i]][General_TJ_data[2][i],1],
    #            '^k',
    #            markeredgecolor= 'k',
    #            markersize=TJ_markersize)
    
    col = 0    
    for i in ['3', '33', '339']:
        for n in range(len(Data_TJ[i][0])):
            plt.plot(verts[Data_TJ[i][0][n]][Data_TJ[i][2][n],0],
                                 verts[Data_TJ[i][0][n]][Data_TJ[i][2][n],1],
                                 TJ_color[col]+TJ_marker[col],
                                 markeredgecolor= TJ_color[col],
                                 markersize=TJ_markersize)
        col += 1

    plt.savefig(path[m] + 'triple_out/Triple_junction_map.png',bbox_inches=0,pad_inches=0, dpi=1200,transparent=True)