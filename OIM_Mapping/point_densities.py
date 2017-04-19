# gaussian filter
import numpy as np
import scipy.ndimage as ndi

def grid_density_gaussian_filter(x0, y0, x1, y1, w, h, data):
    kx = (w - 1) / (x1 - x0)
    ky = (h - 1) / (y1 - y0)
    r = 20
    border = r
    imgw = (w + 2 * border)
    imgh = (h + 2 * border)
    img = np.zeros((imgh,imgw))
    for x, y in data:
        ix = int((x - x0) * kx) + border
        iy = int((y - y0) * ky) + border
        if 0 <= ix < imgw and 0 <= iy < imgh:
            img[iy][ix] += 1
    return ndi.gaussian_filter(img, (r,r))  ## gaussian convolution



def boxsum(img, w, h, r):
    st = [0] * (w+1) * (h+1)
    for x in xrange(w):
        st[x+1] = st[x] + img[x]
    for y in xrange(h):
        st[(y+1)*(w+1)] = st[y*(w+1)] + img[y*w]
        for x in xrange(w):
            st[(y+1)*(w+1)+(x+1)] = st[(y+1)*(w+1)+x] + st[y*(w+1)+(x+1)] - st[y*(w+1)+x] + img[y*w+x]
    for y in xrange(h):
        y0 = max(0, y - r)
        y1 = min(h, y + r + 1)
        for x in xrange(w):
            x0 = max(0, x - r)
            x1 = min(w, x + r + 1)
            img[y*w+x] = st[y0*(w+1)+x0] + st[y1*(w+1)+x1] - st[y1*(w+1)+x0] - st[y0*(w+1)+x1]

def grid_density_boxsum(x0, y0, x1, y1, w, h, data):
    kx = (w - 1) / (x1 - x0)
    ky = (h - 1) / (y1 - y0)
    r = 15
    border = r * 2
    imgw = (w + 2 * border)
    imgh = (h + 2 * border)
    img = [0] * (imgw * imgh)
    for x, y in data:
        ix = int((x - x0) * kx) + border
        iy = int((y - y0) * ky) + border
        if 0 <= ix < imgw and 0 <= iy < imgh:
            img[iy * imgw + ix] += 1
    for p in xrange(4):
        boxsum(img, imgw, imgh, r)
    a = np.array(img).reshape(imgh,imgw)
    b = a[border:(border+h),border:(border+w)]
    return b
    #n = 1000
    ## data points range
    #data_ymin = -2.
    #data_ymax = 2.
    #data_xmin = -2.
    #data_xmax = 2.
    ## view area range
    #view_ymin = -.5
    #view_ymax = .5
    #view_xmin = -.5
    #view_xmax = .5
    ## generate data
    #xl = np.random.uniform(data_xmin, data_xmax, n)    
    #yl = np.random.uniform(data_ymin, data_ymax, n)
    #zl = np.random.uniform(0, 1, n)
    #
    ## get visible data points
    #xlvis = []
    #ylvis = []
    #for i in range(0,len(xl)):
    #    if view_xmin < xl[i] < view_xmax and view_ymin < yl[i] < view_ymax:
    #        xlvis.append(xl[i])
    #        ylvis.append(yl[i])
    #
    #fig = plt.figure()

#plt3 = fig.add_subplot(223)
#plt3.set_axis_off()
#t0 = time.clock()
#zd = grid_density_gaussian_filter(view_xmin, view_ymin, view_xmax, view_ymax, 256, 256, zip(xl, yl))
##plt.title('ndi.gaussian_filter - '+str(time.clock()-t0)+"sec")
#plt.imshow(zd , origin='lower', extent=[view_xmin, view_xmax, view_ymin, view_ymax])
#plt.scatter(xlvis, ylvis)
