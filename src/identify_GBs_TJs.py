#!/usr/bin/env python
from numba import jit
import numpy as np
from src.crystallographic_calculations import euler2rotmat, rotmat2misor_angle
from src.symmetry_groups import get_sigma_variants_crit, get_symmetry_group

''' Identify General and special grain boundaries '''
def get_neighbor_direction(Grain_ID, NROWS, NCOLS_EVEN, NCOLS_ODD):
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
    LS_EVEN = [NCOLS_EVEN_ODD*x -NCOLS_EVEN for x in range(1,int(NROWS/2+1))]
    for i in LS_EVEN:
        neighbor_direction[i] = [1, NCOLS_ODD, NCOLS_EVEN]
        boundaries[i] = [[1,2],[0,1],[-1,0]]
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
    return neighbor_direction, boundaries


def get_GB_IDs(Grain_ID, boundaries, neighbor_direction, GB_point_ID, GB_vert_ID1, GB_vert_ID2, GB_neighbor_ID):
    for i in range(len(Grain_ID)):
        for k in range(len(neighbor_direction[i])):
            if Grain_ID[i] != Grain_ID[i + neighbor_direction[i][k]]:                    
                    GB_point_ID.append(i)
                    GB_neighbor_ID.append(i+neighbor_direction[i][k])
                    GB_vert_ID1.append(boundaries[i][k][0])
                    GB_vert_ID2.append(boundaries[i][k][1])
    return GB_point_ID, GB_neighbor_ID, GB_vert_ID1, GB_vert_ID2

def identify_GBs(Data_GB, Grain_ID, NROWS, NCOLS_EVEN, NCOLS_ODD):
    neighbor_direction, boundaries = get_neighbor_direction(Grain_ID, NROWS, NCOLS_EVEN, NCOLS_ODD)
    GB_point_ID = []
    GB_neighbor_ID = []
    GB_vert_ID1 = []
    GB_vert_ID2 = []
    Data_GB['General_GB'] = get_GB_IDs(Grain_ID, boundaries, neighbor_direction, GB_point_ID, GB_vert_ID1, GB_vert_ID2, GB_neighbor_ID)
   
def get_disorientation(boundary_point, neighbor_point, phi1, Phi, phi2):
    g_point  = euler2rotmat([phi1[boundary_point], Phi[boundary_point], phi2[boundary_point]])
    g_neighbor = euler2rotmat([phi1[neighbor_point], Phi[neighbor_point], phi2[neighbor_point]])
    g = np.dot(g_point, g_neighbor.transpose())
    return g

def get_Sigma_GB_IDs(GB_data, phi1, Phi, phi2, Sigma_variants, Symmetry_group, Brandon_crit, Sigma_GB_point_ID, Sigma_GB_neighbor_ID, Sigma_GB_vert_ID1, Sigma_GB_vert_ID2):
    for i in  range(len(GB_data[0])):
        g = get_disorientation(GB_data[0][i], GB_data[1][i], phi1, Phi, phi2)           
        for n in Sigma_variants:             
            Theta = rotmat2misor_angle( np.dot(g, n.transpose()), Symmetry_group)
            if Theta < Brandon_crit:
                Sigma_GB_point_ID.append(GB_data[0][i])
                Sigma_GB_neighbor_ID.append(GB_data[1][i])
                Sigma_GB_vert_ID1.append(GB_data[2][i])
                Sigma_GB_vert_ID2.append(GB_data[3][i])
                break
    return Sigma_GB_point_ID, Sigma_GB_neighbor_ID, Sigma_GB_vert_ID1, Sigma_GB_vert_ID2


def identify_Sigma_GBs(Data_GB, Sigma_CSL, crit, crystal_sys, Bravais_lattice, phi1, Phi, phi2):
    Symmetry_group = get_symmetry_group(crystal_sys) 
    Sigma_variants, Brandon_crit = get_sigma_variants_crit(Sigma_CSL, crit, crystal_sys, Bravais_lattice)
    Sigma_GB_point_ID = []
    Sigma_GB_neighbor_ID = []
    Sigma_GB_vert_ID1 = []
    Sigma_GB_vert_ID2 = []
    GB_data = Data_GB['General_GB']
    Sigma_GB_point_ID, Sigma_GB_neighbor_ID, Sigma_GB_vert_ID1, Sigma_GB_vert_ID2 = get_Sigma_GB_IDs(GB_data, phi1, Phi, phi2, Sigma_variants, Symmetry_group, Brandon_crit, Sigma_GB_point_ID, Sigma_GB_neighbor_ID, Sigma_GB_vert_ID1, Sigma_GB_vert_ID2)
    Data_GB[Sigma_CSL] = [Sigma_GB_point_ID, Sigma_GB_neighbor_ID, Sigma_GB_vert_ID1, Sigma_GB_vert_ID2]


''' Identify General and special triple junctions '''

def identify_TJs(Data_GB, Data_TJ, Grain_ID, NCOLS_ODD, NCOLS_EVEN):
    TJ_point_ID = []
    TJ_neighbor_ID = []
    TJ_vert_ID = []
    GB_data = Data_GB['General_GB']
    for i in range(len(GB_data[0])-1):
        if ((GB_data[0][i] == GB_data[0][i + 1])
            and (Grain_ID[GB_data[0][i]] != Grain_ID[GB_data[1][i]])
            and (Grain_ID[GB_data[0][i+1]] != Grain_ID[GB_data[1][i + 1]])
            and (Grain_ID[GB_data[1][i]] != Grain_ID[GB_data[1][i + 1]])
            and (GB_data[1][i+1]-GB_data[1][i] == NCOLS_EVEN
                 or GB_data[1][i]-GB_data[1][i+1] == 1)):   
                    TJ_point_ID.append(GB_data[0][i])
                    TJ_neighbor_ID.append([GB_data[1][i], GB_data[1][i+1]])
                    TJ_vert_ID.append(GB_data[2][i])
    Data_TJ['General_TJ'] = [TJ_point_ID, TJ_neighbor_ID, TJ_vert_ID]

def get_disorientation_trip(triple_point, neighbor_point, phi1, Phi, phi2):
    g_point  = euler2rotmat([phi1[triple_point ],Phi[triple_point],phi2[triple_point]])
    g_neighbor1 = euler2rotmat([phi1[neighbor_point[0]],Phi[neighbor_point[0]],phi2[neighbor_point[0]]])
    g_neighbor2 = euler2rotmat([phi1[neighbor_point[1]],Phi[neighbor_point[1]],phi2[neighbor_point[1]]])
    g1 = np.dot(g_point, g_neighbor1.transpose())
    g2 = np.dot(g_point, g_neighbor2.transpose())
    g3 = np.dot(g_neighbor1, g_neighbor2.transpose())
    return g1, g2, g3



def identify_Sigma_TJs(Data_GB,
                       Data_TJ,
                       Triple_junction,
                       crit,
                       crystal_sys,
                       Bravais_lattice,
                       phi1, Phi, phi2):    
    Symmetry_group = get_symmetry_group(crystal_sys)
    Sigma_variants_crit = []
    for i in Triple_junction:
        Sigma_variants_crit.append(get_sigma_variants_crit(i,crit, crystal_sys, Bravais_lattice))

    Sigma_TJ_point_ID = []
    Sigma_TJ_neighbor_ID = []
    Sigma_TJ_vert_ID = []
    
    Sigma_Sigma_TJ_point_ID = []
    Sigma_Sigma_TJ_neighbor_ID = []
    Sigma_Sigma_TJ_vert_ID = []
    
    Sigma_Sigma_Sigma_TJ_point_ID = []
    Sigma_Sigma_Sigma_TJ_neighbor_ID = []
    Sigma_Sigma_Sigma_TJ_vert_ID = []
    
    TJ_data = Data_TJ['General_TJ']
    for i in  range(len(TJ_data[0])):
        g1, g2, g3 = get_disorientation_trip(TJ_data[0][i],TJ_data[1][i],phi1,Phi,phi2)
        g = [g1, g2, g3]
        tj = []
        boundary_no = []
        for j in range(len(Sigma_variants_crit)):
            for n in Sigma_variants_crit[j][0]:
                for k in range(len(g)):
                    rotmat_sym = g[k] * n.transpose()    
                    Theta = rotmat2misor_angle(rotmat_sym,Symmetry_group)
                    if Theta < Sigma_variants_crit[j][1] and k not in boundary_no:            
                        tj.append(Triple_junction[j])
                        boundary_no.append(k)
                        break
            if len(tj) < j:
                break
        if tj == [Triple_junction[0]]:
                Sigma_TJ_point_ID.append(TJ_data[0][i])               
                Sigma_TJ_neighbor_ID.append(TJ_data[1][i])
                Sigma_TJ_vert_ID.append(TJ_data[2][i])                
        if tj == [Triple_junction[0],Triple_junction[1]]:            
                Sigma_Sigma_TJ_point_ID.append(TJ_data[0][i])
                Sigma_Sigma_TJ_neighbor_ID.append(TJ_data[1][i])
                Sigma_Sigma_TJ_vert_ID.append(TJ_data[2][i])                
        if tj == [Triple_junction[0],Triple_junction[1],Triple_junction[2]]:
                Sigma_Sigma_Sigma_TJ_point_ID.append(TJ_data[0][i])
                Sigma_Sigma_Sigma_TJ_neighbor_ID.append(TJ_data[1][i])
                Sigma_Sigma_Sigma_TJ_vert_ID.append(TJ_data[2][i])
  
    Data_TJ[Triple_junction[0]] = [Sigma_TJ_point_ID, Sigma_TJ_neighbor_ID, Sigma_TJ_vert_ID]            
    Data_TJ[Triple_junction[0] + Triple_junction[1]] = [Sigma_Sigma_TJ_point_ID, Sigma_Sigma_TJ_neighbor_ID, Sigma_Sigma_TJ_vert_ID]
    Data_TJ[Triple_junction[0] + Triple_junction[1] + Triple_junction[2]] = [Sigma_Sigma_Sigma_TJ_point_ID,  Sigma_Sigma_Sigma_TJ_neighbor_ID, Sigma_Sigma_Sigma_TJ_vert_ID]
