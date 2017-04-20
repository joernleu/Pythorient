#!/usr/bin/env python
import numpy as np
import os
from matplotlib import pyplot as plt
from matplotlib import cm as cm
# my_funcs
from IO.read_data_ang_file import read_data_hex_ang,det_Step_Size
from IO.write_EBSD_data import write_GB_data, write_TJ_data, write_GB_stats,write_TJ_stats
from OIM_Mapping.get_verts import get_verts
from OIM_Mapping.identify_GBs import identify_GBs,identify_Sigma_GBs
from OIM_Mapping.identify_TJs import identify_TJs,identify_Sigma_TJs
from OIM_Mapping.plot_hex_map import plot_hex_map,plot_hex_rgb_map, plot_line_map
from OIM_Mapping.point_densities import grid_density_gaussian_filter,grid_density_boxsum

def triple_out(path, ang_file, crystal_sys = 'Cubic', Bravais_lattice = 'fcc'):    
    
    '''Import OIM files and create outpu folder'''    
    if not os.path.exists(path+'triple_out/'):
        os.makedirs(path+'triple_out/')
    
    my_file =path + ang_file
    headerlines, phi1, Phi, phi2, x, y, IQ, CI, Fit, Grain_ID, edge, phase = read_data_hex_ang(my_file)
    XSTEP,YSTEP, NCOLS_ODD, NCOLS_EVEN, NROWS = det_Step_Size(x,y)
    verts = get_verts(x,y,YSTEP)

        
    '''Select CSL GBs and TJs and plot parameters'''    
    Sigma_boundary = ['3', '9']
    GB_color = ['y', 'r']
    GB_linewidth = 1 * YSTEP**0.5
    Triple_junction = [['3', '3', '9']]
    TJ_color = ['g', 'y', 'r']
    TJ_marker = ['>', '<', 'v']
    TJ_markersize = 5 * YSTEP**0.5
    
    ''' Identify General and CSL GBs'''
    Data_GB = {}
    identify_GBs(Data_GB, Grain_ID, NROWS, NCOLS_EVEN, NCOLS_ODD)
    
    for i in Sigma_boundary:
        identify_Sigma_GBs(Data_GB,
                            i,
                            'Brandon',
                            crystal_sys,
                            Bravais_lattice,
                            phi1, Phi, phi2)
    
    
    ''' Identify General and CSL TJs'''    
    Data_TJ = {}    
    identify_TJs(Data_GB, Data_TJ, Grain_ID, NCOLS_ODD, NCOLS_EVEN)

    for i in Triple_junction:
        identify_Sigma_TJs(Data_GB,
                            Data_TJ,
                            i,
                            'Brandon',
                            crystal_sys,
                            Bravais_lattice,
                            phi1, Phi, phi2)

    
    '''Export GB data and statistics'''    
    no_General_GBs = len(Data_GB['General_GB'][0])
    write_GB_data(Data_GB, 'General_GB', path)
    write_GB_stats(Data_GB, no_General_GBs,'General_GB', path)
    
    for i in range(len(Sigma_boundary)):
        write_GB_data(Data_GB,Sigma_boundary[i], path)
        write_GB_stats(Data_GB,no_General_GBs,Sigma_boundary[i], path)
    
    '''Export TJ data and statistics'''
    no_General_TJs = len(Data_TJ['General_TJ'][0])
    write_TJ_data(Data_TJ, 'General_TJ', path)
    write_TJ_stats(Data_TJ, no_General_TJs, 'General_TJ', path) 
    
    for i in range(len(Triple_junction)):
        for n in range(len(Triple_junction[i])):        
            Triple = Triple_junction[i][0:n+1]
            Triple =''.join(Triple)
            write_TJ_data(Data_TJ, Triple, path)
            write_TJ_stats(Data_TJ, no_General_TJs, Triple,path)
    

    
    '''Plot density map from selected TJs'''
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
    alpha_val = 0.5
    x_TJ = []
    y_TJ = []
    for i in ['3']:
        for n in range(len(Data_TJ[i][0])):
             x_TJ.append(verts[Data_TJ[i][0][n]][Data_TJ[i][2][n],0])
             y_TJ.append(verts[Data_TJ[i][0][n]][Data_TJ[i][2][n],1])
                        
    
    density = 512
    zd = grid_density_gaussian_filter(view_xmin, view_ymin, view_xmax, view_ymax, density, density, zip(x_TJ, y_TJ))
    #zd = grid_density_boxsum(view_xmin, view_ymin, view_xmax, view_ymax, density, density, zip(x_g, y_g))
    plt.imshow(zd ,cmap=cm.Blues, origin='lower', extent=[view_xmin, view_xmax, view_ymin, view_ymax])
    plt.colorbar(shrink=.5)
    
    '''Plot GBs'''
    xy_GB = []
    for n in range(len(Data_GB['General_GB'][0])):
        xy_GB.append(verts[Data_GB['General_GB'][0][n]][Data_GB['General_GB'][2][n],:])
    plot_line_map(xy_GB, ax,'k',GB_linewidth,'General GB')
    
    
    legend_names_GB = ['$\Sigma$'+ CSL + ' GB' for CSL in Sigma_boundary]
    col = 0
    for i in Sigma_boundary:
        xy_GB = []
        for n in range(len(Data_GB[i][0])):
            xy_GB.append(verts[Data_GB[i][0][n]][Data_GB[i][2][n],:])
        plot_line_map(xy_GB,ax,GB_color[col], GB_linewidth, legend_names_GB[col])
        col += 1

    '''Plot TJs'''
    for i in range(len(Data_TJ['General_TJ'][0])):
        plt.plot(verts[Data_TJ['General_TJ'][0][i]][Data_TJ['General_TJ'][2][i],0],
                verts[Data_TJ['General_TJ'][0][i]][Data_TJ['General_TJ'][2][i],1],
                '^k',
                markeredgecolor = 'k',
                markersize = TJ_markersize)
    
    col = 0
    legend_names_TJ = ['$\Sigma$3-X-X TJ', '$\Sigma$3-$\Sigma$3-X TJ', '$\Sigma$3-$\Sigma$3-$\Sigma$9 TJ']
    for i in ['3', '33', '339']:
        for n in range(len(Data_TJ[i][0])):
            plt.plot(verts[Data_TJ[i][0][n]][Data_TJ[i][2][n],0],
                                 verts[Data_TJ[i][0][n]][Data_TJ[i][2][n],1],
                                 TJ_color[col]+TJ_marker[col],
                                 markeredgecolor = TJ_color[col],
                                 markersize = TJ_markersize)
        plt.plot([],[],
                                 TJ_color[col]+TJ_marker[col],
                                 markeredgecolor = TJ_color[col],
                                 markersize = TJ_markersize,
                                 label = legend_names_TJ[col])        
        col += 1
    plt.legend(loc = 1,framealpha = 1, numpoints = 1,prop={'size':6})
    
    '''Save the Plot'''
    plt.savefig(path + 'triple_out/Triple_junction_map.png',bbox_inches=0,pad_inches=0, dpi=1200,transparent=True)


if __name__ == '__main__':
    '''Put folders and files here'''
    path = ['Test_Files/200_hexagonal_cells_rand_twin_distr/']
    ang_file =['test_grain.ang']
    crystal_sys = 'Cubic'
    Bravais_lattice = 'fcc'
    for i in range(len(path)):
        triple_out(path[i], ang_file[i], crystal_sys, Bravais_lattice)
    