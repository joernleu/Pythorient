#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection,LineCollection
import matplotlib as mpl

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