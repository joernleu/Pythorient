#!/usr/bin/env python
from Crystal_Calc.crystallographic_calculations import euler2rotmat
from Crystal_Calc.get_sigma_variants_crit import get_sigma_variants_crit
from Crystal_Calc.crystallographic_calculations import rotmat2misor_angle
from Crystal_Calc.get_symmetry_group import *

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
                    TJ_vert_ID.append(GB_data[2][i][0])
    Data_TJ['General_TJ'] = [TJ_point_ID, TJ_neighbor_ID, TJ_vert_ID]


def get_disorientation_trip(triple_point,neighbor_point,phi1,Phi,phi2):
    g_point  = euler2rotmat([phi1[triple_point ],Phi[triple_point],phi2[triple_point]])
    g_neighbor1 = euler2rotmat([phi1[neighbor_point[0]],Phi[neighbor_point[0]],phi2[neighbor_point[0]]])
    g_neighbor2 = euler2rotmat([phi1[neighbor_point[1]],Phi[neighbor_point[1]],phi2[neighbor_point[1]]])
    g1 = g_point * g_neighbor1.transpose()
    g2 = g_point * g_neighbor2.transpose()
    g3 = g_neighbor1 *g_neighbor2.transpose()
    return [g1,g2,g3]

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
        g = get_disorientation_trip(TJ_data[0][i],TJ_data[1][i],phi1,Phi,phi2)
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
            if len(tj)<j:
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
