import numpy as np
def write_to_hex_ang(my_file, Conv_file,XSTEP,YSTEP, NCOLS_ODD, NCOLS_EVEN, NROWS, headerlines, phi1, Phi, phi2, x, y, IQ, CI, ID, edge, phase):
    f1 = open(my_file, 'r')
    f2 = open(Conv_file, 'w')
    for i in range(headerlines-1):        
        header = f1.readline()
        f2.write(header)
        if header[0:5] == '# Use':
                    header = f1.readline()
                    f2.write(header)
                    f2.write('# Phase 1 \n')
            
        if header[0:6] == '# GRID':
            f2.write('# XSTEP: '+str(XSTEP)+'\n'
                     '# YSTEP: '+str(YSTEP)+'\n'
                     '# NCOLS_ODD: '+str(NCOLS_ODD)+'\n'
                     '# NCOLS_EVEN: '+str(NCOLS_EVEN)+'\n'
                     '# NROWS: '+str(NROWS)+'\n')
    for n in range(len(x)):
        f2.write('%1.5f \t %1.5f \t %1.5f \t %3.5f \t %3.5f \t %3.3f \t %1.3f \t %1.0f \t %5.0f \t %1.3f \n' % (phi1[n], Phi[n], phi2[n], x[n], y[n], IQ[n], CI[n], ID[n], edge[n], phase[n]))
        
    f1.close()
    f2.close()


    
def write_GB_data(Data,Sigma_name, path):
    GB_data = Data[Sigma_name]
    GB_vert_ID_1 = [item[0] for item in GB_data[2]]
    GB_vert_ID_2 = [item[1] for item in GB_data[2]]
    data = np.column_stack((GB_data[0], GB_data[1],GB_vert_ID_1,GB_vert_ID_2))
        
    outfile = open(path + 'triple_out/'+Sigma_name+'_GBs.txt', 'w')
    for row in data:
        for column in row:
            outfile.write('%i\t' % column)
        outfile.write('\n')
    outfile.close()
    

def write_TJ_data(Data, Sigma_name, path):
    TJ_data = Data[Sigma_name]
    TJ_neighbor_ID1 = [item[0] for item in TJ_data[1]]
    TJ_neighbor_ID2 = [item[1] for item in TJ_data[1]]
    data = np.column_stack((TJ_data[0], TJ_neighbor_ID1,TJ_neighbor_ID2,TJ_data[2]))
                
    outfile = open(path + 'triple_out/'+Sigma_name +'_TJs.txt', 'w')
    for row in data:
        for column in row:
            outfile.write('%i\t' % column)
        outfile.write('\n')
    outfile.close()


def write_GB_stats(Data, no_General_GBs, Sigma_name, path):
    GB_data = Data[Sigma_name]
    no_GBs = float(len(GB_data[0]))
    frac_GBs = no_GBs/no_General_GBs
    data = [no_GBs,frac_GBs]
    outfile = open(path + 'triple_out/'+Sigma_name +'_GB_stats.txt', 'w')
    for column in data:
        outfile.write('%f\t' % column)
    outfile.write('\n')
    outfile.close()

def write_TJ_stats(Data, no_General_TJs, Sigma_name, path):
    TJ_data = Data[Sigma_name]
    no_TJs = float(len(TJ_data[0]))
    frac_TJs = no_TJs/no_General_TJs
    data = [no_TJs,frac_TJs]
    outfile = open(path + 'triple_out/'+Sigma_name +'_TJ_stats.txt', 'w')
    for column in data:
        outfile.write('%f\t' % column)
    outfile.write('\n')
    outfile.close()
    
               