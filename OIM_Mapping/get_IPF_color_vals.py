#!/usr/bin/env python
import numpy as np
def get_IPF_color_vals(ang):
    h = np.sin(ang[1] * np.pi / 180) * np.sin(ang[2] * np.pi / 180)
    k = np.sin(ang[1] * np.pi / 180) * np.cos(ang[2] * np.pi/180)
    l = np.cos(ang[1] * np.pi/180)
    n = (h**2 + k**2 + l**2)**0.5
    h= h * n
    k= k * n
    l= l * n
    hkl_max = min([abs(h), abs(k), abs(l)])
    if hkl_max != 0:
        h = abs(h) / hkl_max
        k = abs(k) / hkl_max
        l = abs(l) / hkl_max
    if h < k:
        h, k = k, h
    if k > l:
        k, l = l, k
    if h > l:
        h, l = l, h
        
    c_max = max([abs(l - h), abs(h - k), abs(k)])
    
    r = (l - h) / c_max
    g = (h - k) / c_max
    b = k / c_max
    return abs(r), abs(g), abs(b)

if __name__ == '__main__':
    angle = [0, 0, 0]
    angle = [60, 0, 0]
    angle = [45, 45, 0]
    angle = [311.6, 5.8, 290]
    angle = [43.9, 137.3, 157.6]
    #angle = [190,145,10]
    #angle = [275.40583, 214.74557, 128.01273]
    h,k,l=[1., 2., 1.]
    r,g,b = get_IPF_color_vals(angle)
    
    print(r, g, b)

