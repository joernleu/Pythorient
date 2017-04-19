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
    GB_linewidth =1*YSTEP**0.5
    Triple_junction = [['3','3','9']]
    TJ_color =['g','y','r']
    TJ_marker =['>','<','v']
    TJ_markersize = 5*YSTEP**0.5
    
    ## General boundaries 
    #---------------------------------------
    General_GB_data = identify_GBs(Grain_ID,NROWS,NCOLS_EVEN,NCOLS_ODD)
    #print(General_GB_data[0])
    
    ## CSL Sigma boundaries
    ##---------------------------------------
    Sigma_GB_data =[]
    for i in Sigma_boundary:
        Sigma_GB_data.append(identify_Sigma_GBs(i,
                                                'Brandon',
                                                crystal_sys,
                                                Bravais_lattice,
                                                phi1, Phi, phi2,
                                                General_GB_data))
    
    
    ## Determine and plot general and CSL triple Junctions
    #-------------------------------------------------------------
    
    ## General triple junctions
    #---------------------------------------
    General_TJ_data = identify_TJs(Grain_ID,General_GB_data,NCOLS_ODD,NCOLS_EVEN)
    ## CSL Sigma triple junctions
    #---------------------------------------
    Sigma_TJ_data = []
    for i in Triple_junction:
        Sigma_TJ_data.append(identify_Sigma_TJs(i,
                                                'Brandon',
                                                crystal_sys,
                                                Bravais_lattice,
                                                phi1, Phi, phi2,
                                                General_TJ_data))
    
    ### Export
    ##-------------------------------------------------------------
    no_General_GBs = len(General_GB_data[0])
    write_GB_data(General_GB_data,'General', path[m])
    write_GB_stats(General_GB_data,no_General_GBs,'General', path[m])
    
    for i in range(len(Sigma_boundary)):
        write_GB_data(Sigma_GB_data[i],Sigma_boundary[i], path[m])
        write_GB_stats(Sigma_GB_data[i],no_General_GBs,Sigma_boundary[i], path[m])
    
    no_General_TJs = len(General_TJ_data[0])
    write_TJ_data(General_TJ_data,'General', path[m])
    write_TJ_stats(General_TJ_data,no_General_TJs,'General', path[m]) 
    
    for i in range(len(Triple_junction)):
        for n in range(len(Sigma_TJ_data[i])):        
            Triple = Triple_junction[i][0:n+1]
            Triple ='-'.join(Triple)
            write_TJ_data(Sigma_TJ_data[i][n],Triple,path[m])
            write_TJ_stats(Sigma_TJ_data[i][n],no_General_TJs,Triple,path[m])
    

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
    for i in range(len(Sigma_TJ_data)):
            for p in range(len(Sigma_TJ_data[i])):
                    for n in range(len(Sigma_TJ_data[i][p][0])):
                        x_TJ.append(verts[Sigma_TJ_data[i][p][0][n]][Sigma_TJ_data[i][p][2][n],0])
                        y_TJ.append(verts[Sigma_TJ_data[i][p][0][n]][Sigma_TJ_data[i][p][2][n],1])
                        
    ## view area range
    density = 512
    zd = grid_density_gaussian_filter(view_xmin, view_ymin, view_xmax, view_ymax, density, density, zip(x_TJ, y_TJ))
    ##zd = grid_density_boxsum(view_xmin, view_ymin, view_xmax, view_ymax, density, density, zip(x_g, y_g))
    plt.imshow(zd , origin='lower', extent=[view_xmin, view_xmax, view_ymin, view_ymax])
    #
    plt.savefig(path[m] + 'triple_out/Triple_junction_density.png',bbox_inches=0,pad_inches=0, dpi=900,transparent=True)
    
    
    xy_GB = []
    
    for n in range(len(General_GB_data[0])):
        xy_GB.append(verts[General_GB_data[0][n]][General_GB_data[2][n],:])
    plot_line_map(xy_GB,ax,'k',GB_linewidth)
    #
    legend_names_GB = ['General GB','$\Sigma$3 GB','$\Sigma$9 GB']
    for i in range(len(Sigma_GB_data)):
        xy_GB = []
        for n in range(len(Sigma_GB_data[i][0])):
            xy_GB.append(verts[Sigma_GB_data[i][0][n]][Sigma_GB_data[i][2][n],:])
        plot_line_map(xy_GB,ax,GB_color[i],GB_linewidth)
    #
    #for i in range(len(General_TJ_data)):
    # #   for n in range(len(General_TJ_data[i])):
    #for i in range(len(General_TJ_data[0])):
    #    plt.plot(verts[General_TJ_data[0][i]][General_TJ_data[2][i],0],
    #            verts[General_TJ_data[0][i]][General_TJ_data[2][i],1],
    #            '^k',
    #            markeredgecolor= 'k',
    #            markersize=TJ_markersize)
    
        
    for i in range(len(Sigma_TJ_data)):
            for p in range(len(Sigma_TJ_data[i])):
                    for n in range(len(Sigma_TJ_data[i][p][0])):
                        plt.plot(verts[Sigma_TJ_data[i][p][0][n]][Sigma_TJ_data[i][p][2][n],0],
                                 verts[Sigma_TJ_data[i][p][0][n]][Sigma_TJ_data[i][p][2][n],1],
                                 TJ_color[p]+TJ_marker[p],
                                 markeredgecolor= TJ_color[p],
                                 markersize=TJ_markersize)
    

    plt.savefig(path[m] + 'triple_out/Triple_junction_map.png',bbox_inches=0,pad_inches=0, dpi=1200,transparent=True)
