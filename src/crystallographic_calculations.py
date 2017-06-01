import numpy as np
from numba import jit


# Euler angles in matrix representation
@jit(nopython=True)
def euler2rotmat(ang):
    g = np.zeros((3, 3))
    sa = np.sin(ang[0])
    ca = np.cos(ang[0]) 
    sb = np.sin(ang[1])
    cb = np.cos(ang[1])
    sc = np.sin(ang[2])
    cc = np.cos(ang[2])
    
    g[0, 0] = ca*cc - sa*sc*cb
    g[0, 1] = sa*cc + ca*sc*cb
    g[0, 2] = sc*sb
    g[1, 0] = -ca*sc - sa*cc*cb
    g[1, 1] = -sa*sc + ca*cc*cb
    g[1, 2] = cc*sb
    g[2, 0] = sa*sb
    g[2, 1] = -ca*sb
    g[2, 2] = cb
    #np.array([[ca*cc - sa*sc*cb, sa*cc + ca*sc*cb,  sc*sb],
    #             [-ca*sc - sa*cc*cb, -sa*sc + ca*cc*cb, cc*sb],
    #             [sa*sb, -ca*sb, cb]])
    #g = np.array([[ca*cc - sa*sc*cb, sa*cc + ca*sc*cb,  sc*sb],
    #             [-ca*sc - sa*cc*cb, -sa*sc + ca*cc*cb, cc*sb],
    #             [sa*sb, -ca*sb, cb]])
    
    #g = np.matrix([g1, g2, g3])
    return g

# axis/angle in Euler angles representation
def axisangle2euler(axisangle):
    c = np.cos(axisangle[0]*np.pi/180)
    s = np.sin(axisangle[0]*np.pi/180)
    t =1 - c
    
    u = axisangle[1]/(axisangle[1]**2+axisangle[2]**2+axisangle[3]**2)**0.5
    v = axisangle[2]/(axisangle[1]**2+axisangle[2]**2+axisangle[3]**2)**0.5
    w = axisangle[3]/(axisangle[1]**2+axisangle[2]**2+axisangle[3]**2)**0.5
    g1 = [t*u*u + c, t*u*v - w*s,  t*u*w + v*s]
    g2 = [t*u*v + w*s, t*v*v + c, t*v*w - u*s]
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
    g1 = [t*u*u + c, t*u*v - w*s,  t*u*w + v*s]
    g2 = [t*u*v + w*s, t*v*v + c, t*v*w - u*s]
    g3 = [t*u*w - v*s, t*v*w + u*s, t*w*w + c]
    
    g = np.matrix([g1, g2, g3])
    return g

def ggt(x, y):
    while y != 0:
        x, y = y, x%y
    return x

def round_miller(omega):
    min_idx = 1
    p_min = 1
    # The highest uvw value is chosen to be 20
    for i in range(1, 20, 1):
        omega_2 = [r*i for r in  omega]
        p = abs(omega_2[0] - round(omega_2[0]))
        + abs(omega_2[1] - round(omega_2[1]))
        + abs(omega_2[2] - round(omega_2[2]))
        if p < p_min:
            p_min = p
            min_idx = i
    
    omega = [int(round(i*min_idx)) for i in  omega]
    if ggt(abs(omega[0]), abs(omega[1])) == ggt(abs(omega[0]), abs(omega[2])) == ggt(abs(omega[1]), abs(omega[2])):
        omega = [x/abs(ggt(omega[0], omega[1])) for x in omega]
    return omega


# Calculate ideal orientations in [hkl]<uvw> representation
def euler2miller(ang):
    sa = np.sin(ang[0])
    ca = np.cos(ang[0]) 
    sb = np.sin(ang[1])
    cb = np.cos(ang[1])
    sc = np.sin(ang[2])
    cc = np.cos(ang[2])
    
    g1 = [ca*cc - sa*sc*cb, sa*cc + ca*sc*cb,  sc*sb]
    g2 = [-ca*sc - sa*cc*cb, -sa*sc + ca*cc*cb, cc*sb]
    g3 = [sa*sb, -ca*sb, cb]
    uvw = [g1[0], g2[0], g3[0]]
    hkl = [g1[2], g2[2], g3[2]]
    uvw = round_miller(uvw)
    hkl = round_miller(hkl)
    return hkl, uvw 


# Calculate misorientation axis/angle and misorientation angle

def rotmat2misor_axisangle(rotmat, Symmetry_group):
    rotmat_sym = [rotmat*x for x in Symmetry_group]
    x_trace = [x.trace() for x in rotmat_sym]
    min_idx = 0
    
    for i in range(len(x_trace)):
        if x_trace[min_idx] < x_trace[i]:
            min_idx  = i
    
    Theta = np.arccos(((x_trace[min_idx]) - 1)/2)
    
    omega = (1/(2*np.sin(Theta))
             *[rotmat_sym[min_idx][2,1] - rotmat_sym[min_idx][1,2],
                rotmat_sym[min_idx][0,2] - rotmat_sym[min_idx][2,0],
                rotmat_sym[min_idx][1,0] - rotmat_sym[min_idx][0,1]])
    
    omega = omega.tolist()
    
    omega = round_miller(omega[0])
    
    return omega, Theta*180/np.pi

@jit(nopython=True)
def rotmat2misor_angle(rotmat, Symmetry_group):
    x_trace = 0
    for i in range(len(Symmetry_group)):
        x = np.dot(rotmat, Symmetry_group[i])
        x_tr = x[0, 0] + x[1, 1] + x[2, 2]
        if x_trace < x_tr:
            x_trace = x_tr    
    Theta = np.arccos(((x_trace) - 1)/2)    
    return Theta*180/np.pi


def get_IPF_color_vals(ang):
    h = np.sin(ang[1])*np.sin(ang[2])
    k = np.sin(ang[1])*np.cos(ang[2])
    l = np.cos(ang[1])
    n = (h**2 + k**2 + l**2)**0.5
    h= h * n
    k= k * n
    l= l * n
    hkl_max = min([abs(h), abs(k), abs(l)])
    if hkl_max != 0:
        h = abs(h)/hkl_max
        k = abs(k)/hkl_max
        l = abs(l)/hkl_max
    if h < k:
        h, k = k, h
    if k > l:
        k, l = l, k
    if h > l:
        h, l = l, h
        
    c_max = max([abs(l - h), abs(h - k), abs(k)])
    
    r = (l - h)/c_max
    g = (h - k)/c_max
    b = k/c_max
    return abs(r), abs(g), abs(b)

if __name__ == '__main__':
    Euler_angles1 = [0, 0, 0]
    Euler_angles2 = [90, 45, 0]
    Euler_angles3 = [149, 54, 45]
    Euler_angles_sigma_3 = [63.43, 48.18, 333.4]
    #Euler_angles = [x*np.pi/180 for x in Euler_angles]
    #r,g,b = get_IPF_color_vals(Euler_angles)
    #print(Euler_angles)
    g1 = euler2rotmat(Euler_angles2)
    g2 = euler2rotmat(Euler_angles3)
    print(g1)
    print(g2)
    print(np.dot(g1,g2))
    #ideal_or = euler2miller(Euler_angles)
    #print(ideal_or)