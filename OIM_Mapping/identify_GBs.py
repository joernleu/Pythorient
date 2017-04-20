#!/usr/bin/env python
from Crystal_Calc.crystallographic_calculations import euler2rotmat
from Crystal_Calc.get_sigma_variants_crit import get_sigma_variants_crit
from Crystal_Calc.crystallographic_calculations import rotmat2misor_angle
from Crystal_Calc.get_symmetry_group import *

def identify_GBs(Data_GB, Grain_ID, NROWS, NCOLS_EVEN, NCOLS_ODD):
    NCOLS_EVEN_ODD = NCOLS_EVEN+NCOLS_ODD
    NCOLS_EVEN_EVEN = NCOLS_EVEN+NCOLS_EVEN
    if NROWS%2 == 1:
        TOPLEFT = int(NROWS/2)*(NCOLS_EVEN+NCOLS_ODD)
        TOPRIGHT= int(NROWS/2)*(NCOLS_EVEN+NCOLS_ODD)+NCOLS_EVEN
    else:
        TOPLEFT = int(NROWS/2)*(NCOLS_EVEN+NCOLS_ODD)-NCOLS_EVEN
        TOPRIGHT = int(NROWS/2)*(NCOLS_EVEN+NCOLS_ODD)-1
    
        
    neighbor_direction = [[1, NCOLS_ODD, NCOLS_EVEN]]*len(Grain_ID)
    boundaries = [[[1,2],[0,1] ,[-1,0]]]*len(Grain_ID)
    #Bottom left
    neighbor_direction[0] = [1, NCOLS_ODD]
    boundaries[0] = [[1,2],[0,1]]
    #Bottom right
    neighbor_direction[NCOLS_EVEN]=[NCOLS_EVEN]
    boundaries[NCOLS_EVEN] = [[-1,0]]
    
    #Left side even
    LS_EVEN = [NCOLS_EVEN_ODD*x -NCOLS_EVEN for x in range(1,int(NROWS/2+1))]
    for i in LS_EVEN:
        neighbor_direction[i] = [1, NCOLS_ODD, NCOLS_EVEN]
        boundaries[i] = [[1,2],[0,1],[-1,0]]
        
    
    #Right Side even
    RS_EVEN = [NCOLS_EVEN_ODD*x-1 for x in range(1,int(NROWS/2+1))]
    for i in RS_EVEN:
        neighbor_direction[i] = [NCOLS_ODD, NCOLS_EVEN]
        boundaries[i]= [[0,1] ,[-1,0]]
    
    LS_ODD=[NCOLS_EVEN_ODD*x for x in range(1,int(NROWS/2))]
    for i in LS_ODD:
        neighbor_direction[i] = [1, NCOLS_ODD]
        boundaries[i] = [[1,2],[0,1]]
    
    RS_ODD = [NCOLS_EVEN_ODD*x+NCOLS_EVEN for x in range(1,int(NROWS/2))]
    for i in RS_ODD:
         neighbor_direction[i] = [NCOLS_EVEN]
         boundaries[i] = [[-1,0]]

    
    TOP = [x for x in range(TOPLEFT,TOPRIGHT)]
    for i in TOP:
        neighbor_direction[i] = [1]
        boundaries[i] = [[1,2]]
        
    
    neighbor_direction[TOPRIGHT] = [0]
    boundaries[TOPRIGHT] = [[0]]
    GB_point_ID = []
    GB_vert_ID = []
    GB_neighbor_ID = []
    for i in range(len(Grain_ID)):
        for k in range(len(neighbor_direction[i])):
            if Grain_ID[i] != Grain_ID[i + neighbor_direction[i][k]]:
                    GB_point_ID.append(i)
                    GB_neighbor_ID.append(i+neighbor_direction[i][k])
                    GB_vert_ID.append(boundaries[i][k])
    
    Data_GB['General_GB'] = [GB_point_ID,GB_neighbor_ID, GB_vert_ID]
    


def get_disorientation(boundary_point,neighbor_point,phi1,Phi,phi2):
    g_point  = euler2rotmat([phi1[boundary_point],Phi[boundary_point],phi2[boundary_point]])
    g_neighbor = euler2rotmat([phi1[neighbor_point],Phi[neighbor_point],phi2[neighbor_point]])
    g = g_point * g_neighbor.transpose()
    return g

def identify_Sigma_GBs(Data_GB, Sigma_CSL, crit, crystal_sys, Bravais_lattice, phi1, Phi, phi2):
    Symmetry_group = get_symmetry_group('Cubic') 
    Sigma_variants, Brandon_crit = get_sigma_variants_crit(Sigma_CSL, crit, crystal_sys, Bravais_lattice)
    Sigma_GB_point_ID = []
    Sigma_GB_neighbor_ID = []
    Sigma_GB_vert_ID = []
    GB_data = Data_GB['General_GB']
    for i in  range(len(GB_data[0])):
        g = get_disorientation(GB_data[0][i],GB_data[1][i],phi1,Phi,phi2)           
        for n in Sigma_variants:            
            #x_trace = [x.trace() for x in rotmat_sym]    
            Theta = rotmat2misor_angle( g * n.transpose(), Symmetry_group)
            #print(Theta)
            if Theta < Brandon_crit:
                Sigma_GB_point_ID.append(GB_data[0][i])                
                Sigma_GB_neighbor_ID.append(GB_data[1][i])
                Sigma_GB_vert_ID.append(GB_data[2][i])
                break
    
    Data_GB[Sigma_CSL] = [Sigma_GB_point_ID, Sigma_GB_neighbor_ID, Sigma_GB_vert_ID] 