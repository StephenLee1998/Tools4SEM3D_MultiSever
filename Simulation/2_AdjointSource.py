import numpy as np
import glob
import math
import os, re
from Utils.Utils import read_par, get_info, find_line_in_file, read_dat, read_semd, write_dat, __SNR__
from Utils.Signals import bpfilter, interpolate_data
from calculate_adjoint_source import cc_traveltime_misfit
from shutil import rmtree
import matplotlib as mpl
mpl.use('TkAgg')
from matplotlib import pyplot as plt


def get_seismogram(data_dir, file_suffix, fs, sta, array, period,channel):
    ### read the siesmograms either in binary format 

    if file_suffix.endswith('semd'):
        file_path = f'{data_dir}/{array}.{sta}.{channel}.{file_suffix}'
        if os.path.exists(file_path):
            [t,data] = read_semd(file_path)
        else:
            return False, False, False
    else:
        file_path = f'{data_dir}/{array}.{sta}-{file_suffix}.{channel}.dat' 
        
        if os.path.exists(file_path) == False:
            file_path = f'{data_dir}/{array}.{file_suffix}-{sta}.{channel}.dat'
        if os.path.exists(file_path):
            data = read_dat(file_path)
        else:
            return False, False, False
    
    ### bandpass the seismograms
    data = bpfilter(data, fs, period)

    T = np.arange(len(data)) / fs 
    if file_suffix.endswith('semd') == True:
        T = t
    
    return True, T, data/max(abs(data))

def calculate_adjoint_source(observed, synthetic, deltat, window, config):  # NOQA
    if config.adjsrc_type      == 'cc_traveltime':
        [adj, tshift]  =  cc_traveltime_misfit(observed, synthetic, deltat, window, config)

    return adj, tshift

def __Window__(current_index, source, loc, out_path, config):
    global noi_len
    out_path = out_path + '/' + source
    os.mkdir(out_path)
    os.mkdir(out_path+'/adj')
    os.mkdir(out_path+'/obs')
    os.mkdir(out_path+'/syn')
    

    source_id = source.split('S')[1]
    ### obtain the sampling frequency from runbase parametre file and cfs
    file   = glob.glob(config.obs_path+'/*.dat')[0]
    data   = np.loadtxt(file)
    dt     = data[3,0] - data[2,0]
    fs_obs = 1/dt

    dt     = find_line_in_file('../runbase/run0001/DATA/Par_file','DT                              =')
    fs_syn = 1/float(dt)

    fs = {'fs_obs': fs_obs, 'fs_syn': fs_syn}

    ### obtain the sources and stations of simulations
    with open('../runbase/DATA/STATIONS_FILTERED', 'r') as f1:
        stas = [line.strip().split()[:4] for line in f1]
    f1.close()

    syn_path = config.syn_path + '/run'+str(current_index+1).zfill(4)+'/OUTPUT_FILES'

    Xs = float(loc[0]); Ys = float(loc[1])

    Channels = ['FXX','FXY','FXZ']

    f2 = open(out_path+'/STATIONS_ADJOINT','w')
    
    ### loop for the period bands and stations
    for sta, array, x, y in stas:
        sta_flag = False
        
        for channel in Channels:
            #file_syn = out_path + '/syn/'+array+'.'+re.sub("[a-zA-Z]","",sta)+'.'+channel+'.syn'
            #file_obs = out_path + '/obs/'+array+'.'+re.sub("[a-zA-Z]","",sta)+'.'+channel+'.obs'
            #file_adj = out_path + '/adj/'+array+'.'+re.sub("[a-zA-Z]","",sta)+'.'+channel+'.adj'
            file_syn = out_path + '/syn/'+array+'.'+sta+'.'+channel+'.syn'
            file_obs = out_path + '/obs/'+array+'.'+sta+'.'+channel+'.obs'
            file_adj = out_path + '/adj/'+array+'.'+sta+'.'+channel+'.adj'

            for j, period in enumerate(info['periods']):               
                x = float(x); y = float(y)
                dist = math.sqrt(math.pow(Xs - x, 2) + math.pow(Ys - y, 2))/1000                
                ### read the seismograms
                [flag_syn,t_syn,data_syn]=get_seismogram(syn_path,        'semd',           fs['fs_syn'], sta, array, period,channel)
                [flag_obs,t_obs,data_obs]=get_seismogram(config.obs_path, f'NJ{source_id}', fs['fs_obs'], sta, array, period,channel)

                #print(flag_obs,flag_syn)
                if flag_obs and flag_syn:
                    data_obs = interpolate_data(t_obs,data_obs,t_syn)

                if j == 0:
                    adj_src      = np.zeros(len(data_syn))
                    data_obs_sum = np.zeros(len(data_syn))
                    data_syn_sum = np.zeros(len(data_syn))

                data_obs_sum += data_obs
                data_syn_sum += data_syn

                if flag_obs and flag_syn:
                    t1 = dist/float(info['vel_win'][j][1]) + float(info['tshift'][j][1])
                    t2 = dist/float(info['vel_win'][j][0]) + float(info['tshift'][j][0])
                    
                    print(file_obs,t1,t2)
    
                    SNR = __SNR__(t_syn, data_obs, t1, t2, config.noi_len)
    
                    tshift_threshold = period[0] * config.tshift
                    print(SNR)
    
                    if SNR > config.snr:             
                        ### apply window to seismograms
                        sta_flag = True
    
                        window   = [[t1, t2]]
    
                        [adj, tshift] = calculate_adjoint_source(data_obs, data_syn, dt, window, config)
    
                        if tshift < tshift_threshold:
                            adj_src      += adj["adjoint_source"][::-1]/max(adj["adjoint_source"])

            #adj_src = adj_src/max(adj_src)
                    
            write_dat(file_adj, t_syn, adj_src)
            write_dat(file_obs, t_syn, data_obs_sum)
            write_dat(file_syn, t_syn, data_syn_sum)
        if sta_flag:
            f2.write(sta+' '+array+' '+str(x)+' '+str(y)+' 0 0'+'\n')
    f2.close()


if __name__ == "__main__":
    ### This code is disigned for Specfem 3D Cartesian Version 4.1.0, use it after the forward simulations and the 
    ### velocity windows have been selected for your seismograms.  This script read and apply window to the synthetic
    ### and observed seismograms.  The seismograms are designed to be selected in this script by SNR(signal-noise 
    ### ratio) and tshift for the calculation of adjoint source.  The cross-correlation adjoint source can be calculated
    ### for the adjoint simulation.  

    config = read_par('./par_file')
    info   = get_info(config)

    ### construct the dir tree        
    dir_out = './scratch'
    if os.path.exists(dir_out):
        rmtree(dir_out)
    os.mkdir(dir_out)
    os.mkdir(dir_out+'/traces')

    ## read the source information 
    with open(config.sources, 'r') as f1:
        lines = [line.strip().split()[:3] for line in f1]
    for i, (source, x, y) in enumerate(lines):
        __Window__(i, source, [x, y] , dir_out+'/traces/', config)
