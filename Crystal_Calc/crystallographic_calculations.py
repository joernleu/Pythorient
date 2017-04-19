## Crystallographic calculations
# Represent Euler angles in rotation matrix
import numpy as np

#  Euler angles in matrix representation
def euler2rotmat(ang):
    sa = np.sin(ang[0])
    ca = np.cos(ang[0]) 
    sb = np.sin(ang[1])
    cb = np.cos(ang[1])
    sc = np.sin(ang[2])
    cc = np.cos(ang[2])
    
    g1 = [ ca*cc - sa*sc*cb, sa*cc + ca*sc*cb,  sc*sb]
    g2 = [ -ca*sc - sa*cc*cb, -sa*sc + ca*cc*cb, cc*sb]
    g3 = [sa*sb, -ca*sb, cb]
    
    g = np.matrix([g1, g2, g3])
    return g

# axis/angle in Euler angles representation
def axisangle2euler(axisangle):
    c = np.cos(axisangle[0]*np.pi/180)
    s = np.sin(axisangle[0]*np.pi/180)
    t =1 - c
    
    u = axisangle[1]/(axisangle[1]**2+axisangle[2]**2+axisangle[3]**2)**0.5
    v = axisangle[2]/(axisangle[1]**2+axisangle[2]**2+axisangle[3]**2)**0.5
    w = axisangle[3]/(axisangle[1]**2+axisangle[2]**2+axisangle[3]**2)**0.5
    g1 = [ t*u*u + c, t*u*v - w*s,  t*u*w + v*s]
    g2 = [ t*u*v + w*s, t*v*v + c, t*v*w - u*s]
    g3 = [t*u*w - v*s, t*v*w + u*s, t*w*w + c]
    
    phi1 = np.arctan2(-(t*u*w - v*s), (t*v*w + u*s))
    Phi =  np.arccos(t*w*w + c) 
    phi2 = np.arctan2((t*u*w + v*s), (t*v*w - u*s))
    
    return phi1, Phi, phi2


# axis/angle in rotation matrix representation
def axisangle2rotmat(angle,axis):
    c = np.cos(angle*np.pi/180)
    s = np.sin(angle*np.pi/180)
    t =1 - c
    u = axis[0]/(axis[0]**2+axis[1]**2+axis[2]**2)**0.5
    v = axis[1]/(axis[0]**2+axis[1]**2+axis[2]**2)**0.5
    w = axis[2]/(axis[0]**2+axis[1]**2+axis[2]**2)**0.5
    g1 = [ t*u*u + c, t*u*v - w*s,  t*u*w + v*s]
    g2 = [ t*u*v + w*s, t*v*v + c, t*v*w - u*s]
    g3 = [t*u*w - v*s, t*v*w + u*s, t*w*w + c]
    
    g = np.matrix([g1, g2, g3])
    return g

# Calculate ideal orientations in [hkl]<uvw> representation
def euler2miller(ang):
    h = np.sin(ang[1])*np.sin(ang[2])
    k = np.sin(ang[1])*np.cos(ang[2])
    l = np.cos(ang[1])
    n = (h**2+k**2+l**2)**0.5
    h= h*n
    k= k*n
    l= l*n
    

    hkl_max = max([abs(h),abs(k),abs(l)])
    h = h/hkl_max
    k = k/hkl_max
    l = l/hkl_max
    
    u = np.cos(ang[0])*np.cos(ang[2])-np.sin(ang[0])*np.sin(ang[2])*np.cos(ang[1])
    v = -np.cos(ang[0])*np.sin(ang[2])-np.sin(ang[0])*np.cos(ang[2])*np.cos(ang[1])
    w = np.sin(ang[0])*np.sin(ang[2])
 
    n = (u**2+v**2+w**2)**0.5
    u= u * n
    v= v * n
    w= w * n


    uvw_max = max([abs(u),abs(v),abs(w)])
        
    u = u/uvw_max
    v = v/uvw_max
    w = w/uvw_max
    
    res  = 22.
    for i in range(1,22):
        if res>(abs(round(h*i)-h*i)+abs(round(k*i)-k*i)+abs(round(l*i)-l*i)):
            res=(abs(round(h*i)-h*i)+abs(round(k*i)-k*i)+abs(round(l*i)-l*i))
            i_max = i
            if res<0.3:
                break
            
    h = h*i_max
    k = k*i_max
    l = l*i_max
    
    res  = 22.
    for i in range(1,22):
        if res>(abs(round(u*i)-u*i)+abs(round(v*i)-v*i)+abs(round(w*i)-w*i)):
            res=(abs(round(u*i)-u*i)+abs(round(v*i)-v*i)+abs(round(w*i)-w*i))
            i_max = i
            if res<0.3:
                break
  
    u = u*i_max
    v = v*i_max
    w = w*i_max   
    
    h= round(h)
    k= round(k)
    l= round(l)
    u= round(u)
    v= round(v)
    w= round(w)
    
    
    hukvlw =[h*u,k*v,l*w]
    i_max = 0
    for i in range(1,3):
        if hukvlw[i]>hukvlw[i_max]:
            i_max = i
    
    while i:
        res = (h*u + k*v + l*w)*i
        if res == 0.0:
            break
            
        if res == h and i_max==0:
            u = u*i - 1
            v = v*i
            w = w*i
            break
        if res == k and i_max == 1:
            u = u*i
            v = v*i - 1
            w = w*i
            break
        if res == l and i_max == 2:
            u = u*i
            v = v*i
            w = w*i - 1
            break
        i+=1
    
    return [h, k, l], [u, v, w]


# Calculate misorientation axis/angle and misorientation angle
def ggt(x, y):
    while y != 0:
        x, y = y, x%y
    return x


def rotmat2misor_axisangle(rotmat, Symmetry_group):
    rotmat_sym = [rotmat * x for x in Symmetry_group]
    x_trace = [x.trace() for x in rotmat_sym]
    min_idx = 0
    
    for i in range(len(x_trace)):
        if x_trace[min_idx]<x_trace[i]:
            min_idx  = i
    
    Theta = np.arccos(((x_trace[min_idx])-1)/2)
    
    omega = (1/(2* np.sin(Theta))
             * [rotmat_sym[min_idx][2,1] - rotmat_sym[min_idx][1,2],
                rotmat_sym[min_idx][0,2] - rotmat_sym[min_idx][2,0],
                rotmat_sym[min_idx][1,0] - rotmat_sym[min_idx][0,1]])
    
    omega = omega.tolist()
    
    omega = omega[0]
    
    min_idx = 1
    p_min = 1
    # The highest uvw value is chosen to be 20
    for i in range(1, 20, 1):
        omega_2 = [r * i for r in  omega]
        p = abs(omega_2[0] - round(omega_2[0]))
        + abs(omega_2[1] - round(omega_2[1]))
        + abs(omega_2[2] - round(omega_2[2]))
        if p < p_min:
            p_min = p
            min_idx = i
    
    omega = [int(round(i*min_idx)) for i in  omega]
    
    if ggt(abs(omega[0]), abs(omega[1])) == ggt(abs(omega[0]), abs(omega[2])) == ggt(abs(omega[1]), abs(omega[2])):
        omega = [x/abs(ggt(omega[0], omega[1])) for x in omega]
    
    return omega, Theta*180/np.pi


def rotmat2misor_angle(rotmat, Symmetry_group):
    rotmat_sym = [rotmat * x for x in Symmetry_group]
    x_trace = [x.trace() for x in rotmat_sym]
    min_idx = 0
    
    for i in range(len(x_trace)):
        if x_trace[min_idx]<x_trace[i]:
            min_idx  = i
    
    Theta = np.arccos(((x_trace[min_idx])-1)/2)
    
    return Theta*180/np.pi

