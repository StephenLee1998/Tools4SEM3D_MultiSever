import numpy as np
import os
import glob
import math
import matplotlib as mpl
mpl.use('TkAgg')
from matplotlib import pyplot as plt
from Utils.Utils import read_dat, read_semd, bpfilter, find_line_in_file, rms

def plot_seismogram(data_dir, loc, file_suffix, color, fs, stas, info):
    period = info[0]; vel_win = info[1]; tshift = info[2]
    Xs = float(loc[0]); Ys = float(loc[1])
    d = []
    for i, (sta, array, x, y) in enumerate(stas):
        if file_suffix.endswith('.semd'):
            if (i+1) % 5 != 0:
                continue
            file_path = f'{data_dir}/{array}.{sta}{file_suffix}'
            data = read_semd(file_path)
            dt = 0.3
        else:
            if (i+1) % 5 != 0:
                continue
            dt = 0
            file_path = f'{data_dir}/{array}.{sta}-{file_suffix}.dat' 
            print
            if os.path.exists(file_path) == False:
                file_path = f'{data_dir}/{array}.{file_suffix}-{sta}.dat'
            if os.path.exists(file_path):
                data = read_dat(file_path)
            else:
                continue
                
        
        x = float(x); y = float(y)
        dist = math.sqrt(math.pow(Xs - x, 2) + math.pow(Ys - y, 2))/1000
        if dist == 0:
            continue
        if x < Xs:
            dist = -dist
        d.append(dist)
        data = bpfilter(data, fs, period)
          
        T = np.arange(len(data)) / fs - dt
        if file_path.endswith('.dat'):
            Sig_Index = np.where((T>abs(dist)/3) & (T<abs(dist)/1.5))
            Noi_Index = np.where((T>abs(dist)/1.5+10) & (T<abs(dist)/1.5+30))
            SNR = rms(data[Sig_Index])/rms(data[Noi_Index])
            if SNR > 0:
                plt.plot(T, 0.3 * data / max(abs(data[Sig_Index])) + dist, color)
        else:
            plt.plot(T, 0.3 * data / max(abs(data)) + dist, color)

        t = np.arange(0,21)
        
        plt.plot(t+tshift[0], vel_win[0]*t,'k', linewidth=0.3)
        plt.plot(t+tshift[1], vel_win[1]*t,'k', linewidth=0.3)

        plt.plot(t+tshift[0], -vel_win[0]*t,'k',linewidth=0.3)
        plt.plot(t+tshift[1], -vel_win[1]*t,'k',linewidth=0.3)
        

    plt.xlim([0, 20])
    plt.ylim([min(d)-0.5,max(d)+0.5])

    

if __name__ == "__main__":
    ### This code is disigned for Specfem 3D Cartesian Version 4.1.0, use it after the forward simulations
    ### have been completed.  We simplify the process of window picking in ambient noise adjoint tomography
    ### with a simple velocity window applied to all seiemograms.  This program helps to select and taper
    ### the window for the calculation of adjoint source.  Another important function of this script is to
    ### check the fit condition of syn and obs seismograms
    ### Added by Chao, Li  March 11,2024
    sources  = './SOURCES'
    obs_path = '../CFs'

    ### Def the period band selected
    info = {'period_band' : [[0.5,1]  ,[1, 2]   ,[2, 3]],
            'vel_win'     : [[2.3,2.5],[2.3,2.5],[2.5,3.0]],
            'tshift'      : [[1,-1]   ,[2,-2]   ,[5,-2]]} 
    
    ### Get the fs of synthetic and observed seismograms
    dt     = find_line_in_file('../runbase/DATA/Par_file','DT                              =')
    fs_syn = 1/float(dt)

    file   = glob.glob(obs_path+'/*.dat')[0]
    data   = np.loadtxt(file)
    dt     = data[3,0] - data[2,0]
    fs_obs = 1/dt

    fs   = {'fs_obs': fs_obs, 'fs_syn': fs_syn}
    
    ### Get the stations used in forward simulation
    with open('../runbase/DATA/STATIONS_FILTERED', 'r') as f:
        stas = [line.strip().split()[:4] for line in f]

    ### Get the sources used in forward simulation
    with open(sources,'r') as f1:
        lines = [line.strip().split()[:3] for line in f1]
    
    ### Loop for thr events' forward simulations
    for i,(source, xs, ys) in enumerate(lines):
        syn_path = '../runbase/run'+str(i+1).zfill(4)+'/OUTPUT_FILES' 
        source_loc  = [xs,ys]
        source_id = source.split('S')[1]
        print(source_id)

        for j,period in enumerate(info['period_band']):
            
            plt.figure(1,figsize=(15,8))
            # Plot forward seismograms
            plt.subplot(1, 3, j+1)
            plot_seismogram(syn_path, source_loc, '.FXZ.semd',      'k', fs['fs_syn'], stas , [period, info['vel_win'][j], info['tshift'][j]])
            plot_seismogram(obs_path, source_loc, f'NJ{source_id}', 'r', fs['fs_obs'], stas , [period, info['vel_win'][j], info['tshift'][j]])

        plt.show()