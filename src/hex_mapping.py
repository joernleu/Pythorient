#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection,LineCollection
import matplotlib as mpl


def get_verts(x,y,YSTEP):
    numpolys = len(x)
    centers = []
    for i in range(len(x)):
        centers.append([x[i],y[i]])
    #centers = np.array(zip(x,y))
    numverts = 6
    offsets_x = []
    offsets_y = []
    for i in range(numverts):
        offsets_x.append(YSTEP/(3/2.0)*np.sin(i*60*np.pi/180))
        offsets_y.append(YSTEP/(3/2.0)*np.cos(i*60*np.pi/180))
    
    
    vert0 = np.array([offsets_x[0],offsets_y[0]])
    c0 = [vert0 for i in range(numpolys)]
    c0 = np.array(c0)+centers
    
    vert1 = np.array([offsets_x[1],offsets_y[1]])
    c1 = [vert1 for i in range(numpolys)]
    c1 = np.array(c1)+centers
    
    vert2 = np.array([offsets_x[2],offsets_y[2]])
    c2 = [vert2 for i in range(numpolys)]
    c2 = np.array(c2)+centers
    vert3 = np.array([offsets_x[3],offsets_y[3]])
    c3 = [vert3 for i in range(numpolys)]
    c3 = np.array(c3)+centers
    
    vert4 = np.array([offsets_x[4],offsets_y[4]])
    c4 = [vert4 for i in range(numpolys)]
    c4 = np.array(c4)+centers
    
    vert5 = np.array([offsets_x[5],offsets_y[5]])
    c5 = [vert5 for i in range(numpolys)]
    c5 = np.array(c5)+centers
    
    
    verts = np.array([c0,c1,c2,c3,c4,c5])
    verts = np.swapaxes(verts, 0, 1)
    return verts

def plot_hex_map(fig, ax, alpha_val, IQ,verts, title = 'IQ map'):
    z = np.array(IQ)
    coll = PolyCollection(verts, array=z, cmap=mpl.cm.gray, edgecolors='none', antialiased = 1)
    coll.set_alpha(alpha_val)
    ax.add_collection(coll)

def plot_hex_rgb_map(x, y, IPF_r_g_b, title = 'IPF map'):
    coll = PolyCollection(verts, array=None, facecolors = IPF_r_g_b, edgecolors='none', antialiased=1)
    ax.add_collection(coll)

def plot_line_map(xy_GB, ax, GB_color, GB_width, legend_name = 'GB'):
    line_segments = LineCollection(xy_GB,
                                   linewidths = GB_width,
                                   linestyles = 'solid',
                                   colors = GB_color,
                                   label = legend_name)
    ax.add_collection(line_segments)