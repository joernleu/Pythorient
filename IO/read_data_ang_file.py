#!/usr/bin/env python

def read_data_hex_ang(my_file):
    f = open(my_file, 'r')
    for i,l in enumerate(f):
        pass
    f.close()
    f = open(my_file, 'r')
    headerlines=0
    phi1 = []
    Phi  = []
    phi2 = []
    x = []
    y = []
    IQ = []
    CI = []
    Fit = []
    ID = []
    edge = []
    phase = []

    while 1:
        header = f.readline()       
        if header[0]!='#':
            columns =header.strip().split()
            phi1.append(float(columns[0]))
            Phi.append(float(columns[1]))
            phi2.append(float(columns[2]))
            x.append(float(columns[3]))
            y.append(float(columns[4]))
            IQ.append(float(columns[5]))
            CI.append(float(columns[6]))
            Fit.append(float(columns[7]))
            ID.append(float(columns[8]))
            edge.append(float(columns[9]))
            phase.append(columns[10])
            for line in f:
                #line = line.strip()
                columns = line.strip().split()
                #name = columns[2]
                if '#QNAN' in str(columns[0]):
                    columns[0] = 0.0
                phi1.append(float(columns[0]))
                Phi.append(float(columns[1]))
                phi2.append(float(columns[2]))
                x.append(float(columns[3]))
                y.append(float(columns[4]))
                IQ.append(float(columns[5]))
                CI.append(float(columns[6]))
                Fit.append(float(columns[7]))
                ID.append(float(columns[8]))
                edge.append(float(columns[9]))
                phase.append(columns[10])
            f.close()
            break
        headerlines += 1    
    # Loop over lines and extract variables of interest

    return headerlines, phi1, Phi, phi2, x, y, IQ, CI, Fit, ID, edge, phase

def read_GB_data(my_file):
    f = open(my_file, 'r')
    for i,l in enumerate(f):
        pass
    f.close()
    f = open(my_file, 'r')
    headerlines=0
    GB_point_ID  = []
    GB_neighbor_ID = []
    GB_vert_ID = []
    while 1:
        header = f.readline()       
        columns =header.strip().split()
        GB_point_ID.append(int(columns[0]))
        GB_neighbor_ID.append(int(columns[1]))
        GB_vert_ID.append([int(columns[2]),int(columns[3])])
        for line in f:
            #line = line.strip()
            columns = line.strip().split()
            GB_point_ID.append(int(columns[0]))
            GB_neighbor_ID.append(int(columns[1]))
            GB_vert_ID.append([int(columns[2]),int(columns[3])])
        f.close()
        break
        headerlines += 1
    return GB_point_ID,GB_neighbor_ID, GB_vert_ID
    


def read_TJ_data(my_file):
    f = open(my_file, 'r')
    for i,l in enumerate(f):
        pass
    f.close()
    f = open(my_file, 'r')
    headerlines=0
    TJ_point_ID  = []
    TJ_neighbor_ID = []
    TJ_vert_ID = []
    while 1:
        header = f.readline()       
        columns =header.strip().split()
        TJ_point_ID.append(int(columns[0]))
        TJ_neighbor_ID.append([int(columns[1]),int(columns[2])])
        TJ_vert_ID.append(int(columns[3]))
        for line in f:
            #line = line.strip()
            columns = line.strip().split()
            TJ_point_ID.append(int(columns[0]))
            TJ_neighbor_ID.append([int(columns[1]),int(columns[2])])
            TJ_vert_ID.append(int(columns[3]))
        f.close()
        break
        headerlines += 1
    return TJ_point_ID, TJ_neighbor_ID, TJ_vert_ID
#
def det_Step_Size(x,y):
    XSTEP = x[1]
    NCOLS_ODD = 0
    while x[NCOLS_ODD+1]-x[NCOLS_ODD]>0:
        NCOLS_ODD +=1
    NCOLS_ODD = NCOLS_ODD+1
    NCOLS_EVEN = NCOLS_ODD-1
    YSTEP = y[NCOLS_ODD+1]
    NROWS = round((2*len(x)+1)/(NCOLS_ODD+NCOLS_EVEN))
    return XSTEP,YSTEP, NCOLS_ODD, NCOLS_EVEN, int(NROWS)
        