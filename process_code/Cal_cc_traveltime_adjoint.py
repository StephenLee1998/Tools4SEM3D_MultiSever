#!/usr/bin/env python
# -*- encoding: utf8 -*-
from __future__ import absolute_import, division, print_function

import warnings

import numpy as np
from obspy.signal.cross_correlation import xcorr_pick_correction
from scipy.integrate import simps

from utils import window_taper,  generic_adjoint_source_plot

from obspy import read
import os
import glob
from matplotlib import pyplot as plt

def _xcorr_shift(d, s):
    cc = np.correlate(d, s, mode="full")
    time_shift = cc.argmax() - len(d) + 1
    return time_shift


def cc_error(d1, d2, deltat, cc_shift, cc_dlna, sigma_dt_min, sigma_dlna_min):
    """
    Estimate error for dt and dlna with uncorrelation assumption
    """
    nlen_t = len(d1)

    d2_cc_dt = np.zeros(nlen_t)
    d2_cc_dtdlna = np.zeros(nlen_t)

    for index in range(0, nlen_t):
        index_shift = index - cc_shift

        if 0 <= index_shift < nlen_t:
            # corrected by c.c. shift
            d2_cc_dt[index] = d2[index_shift]

            # corrected by c.c. shift and amplitude
            d2_cc_dtdlna[index] = np.exp(cc_dlna) * d2[index_shift]

    # time derivative of d2_cc (velocity)
    d2_cc_vel = np.gradient(d2_cc_dtdlna, deltat)

    # the estimated error for dt and dlna with uncorrelation assumption
    sigma_dt_top = np.sum((d1 - d2_cc_dtdlna)**2)
    sigma_dt_bot = np.sum(d2_cc_vel**2)

    sigma_dlna_top = sigma_dt_top
    sigma_dlna_bot = np.sum(d2_cc_dt**2)

    sigma_dt = np.sqrt(sigma_dt_top / sigma_dt_bot)
    sigma_dlna = np.sqrt(sigma_dlna_top / sigma_dlna_bot)

    if sigma_dt < sigma_dt_min:
        sigma_dt = sigma_dt_min

    if sigma_dlna < sigma_dlna_min:
        sigma_dlna = sigma_dlna_min

    return sigma_dt, sigma_dlna


def subsample_xcorr_shift(d, s):
    """
    Calculate the correlation time shift around the maximum amplitude of the
    synthetic trace with subsample accuracy.
    :param s:
    :param d:
    """
    # Estimate shift and use it as a guideline for the subsample accuracy
    # shift.
    time_shift = _xcorr_shift(d.data, s.data) * d.stats.delta

    # Align on the maximum amplitude of the synthetics.
    pick_time = s.stats.starttime + s.data.argmax() * s.stats.delta

    # Will raise a warning if the trace ids don't match which we don't care
    # about here.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return xcorr_pick_correction(
            pick_time, s, pick_time, d, 20.0 * time_shift,
            20.0 * time_shift, 10.0 * time_shift)[0]


def calculate_adjoint_source(observed, synthetic, config, window,
                             adjoint_src, figure):  # NOQA

    ret_val_p = {}
    ret_val_q = {}

    nlen_data = len(synthetic.data)
    deltat = synthetic.stats.delta

    fp = np.zeros(nlen_data)
    fq = np.zeros(nlen_data)

    misfit_sum_p = 0.0
    misfit_sum_q = 0.0

    # ===
    # loop over time windows
    # ===
    for wins in window:
        left_window_border = wins[0]
        right_window_border = wins[1]

        left_sample = int(np.floor(left_window_border / deltat)) + 1
        nlen = int(np.floor((right_window_border -
                             left_window_border) / deltat)) + 1
        right_sample = left_sample + nlen

        d = np.zeros(nlen)
        s = np.zeros(nlen)

        d[0:nlen] = observed.data[left_sample:right_sample]
        s[0:nlen] = synthetic.data[left_sample:right_sample]

        # All adjoint sources will need some kind of windowing taper
        d = window_taper(d, taper_percentage=config.taper_percentage,
                     taper_type=config.taper_type)
        s = window_taper(s, taper_percentage=config.taper_percentage,
                     taper_type=config.taper_type)

        i_shift = _xcorr_shift(d, s)
        t_shift = i_shift * deltat

        cc_dlna = 0.5 * np.log(sum(d[0:nlen]*d[0:nlen]) /
                               sum(s[0:nlen]*s[0:nlen]))

        sigma_dt, sigma_dlna = cc_error(d, s, deltat, i_shift, cc_dlna,
                                        config.dt_sigma_min,
                                        config.dlna_sigma_min)

        misfit_sum_p += 0.5 * (t_shift/sigma_dt) ** 2
        misfit_sum_q += 0.5 * (cc_dlna/sigma_dlna) ** 2

        dsdt = np.gradient(s, deltat)
        nnorm = simps(y=dsdt*dsdt, dx=deltat)
        fp[left_sample:right_sample] = dsdt[:] * t_shift / nnorm / sigma_dt**2

        mnorm = simps(y=s*s, dx=deltat)
        fq[left_sample:right_sample] =\
            -1.0 * s[:] * cc_dlna / mnorm / sigma_dlna ** 2

    ret_val_p["misfit"] = misfit_sum_p
    ret_val_q["misfit"] = misfit_sum_q

    if adjoint_src is True:
        ret_val_p["adjoint_source"] = fp[::-1]
        ret_val_q["adjoint_source"] = fq[::-1]

    if config.measure_type == "dt":
        if figure:
            generic_adjoint_source_plot(observed, synthetic,
                                        ret_val_p["adjoint_source"],
                                        ret_val_p["misfit"],
                                        window, VERBOSE_NAME)

        return fp,t_shift
        #ret_val_p

    if config.measure_type == "am":
        if figure:
            generic_adjoint_source_plot(observed, synthetic,
                                        ret_val_q["adjoint_source"],
                                        ret_val_q["misfit"],
                                        window, VERBOSE_NAME)

        return ret_val_q

## Cal and output the frequency dependent adj source

flag = False
if flag:
    files = ['SAC_2-6','SAC_5-9','SAC_8-12']
    Out_path = '/project/li_chao/SEM3D/specfem3d/misfit_test/ADJ'
    adjoint_src = False
    figure = False
    class item:
         def __init__(self):
             self.taper_percentage = 0.4
             self.taper_type = 'cos_p10'
             self.measure_type = 'dt'
             self.dt_sigma_min = 0.0001
             self.dlna_sigma_min = 0
    config = item()
    os.chdir(Out_path)
    D = glob.glob("30_80")
    for i in range(len(D)):
        os.chdir(Out_path + '/' + D[i])
        os.system('rm *.dat')

    for j in range(len(files)):
        data_path = '/project/li_chao/SEM3D/specfem3d/misfit_test/Window_Picking_80/' + files[j]
        os.chdir(data_path)
        
        sgfs = glob.glob("*FXZ*sac")
        t = 0.0001*np.arange(23000)
        for i in range(len(sgfs)):
            st = read(sgfs[i])
            tr = st[0]
            window = []
            window.append([int(tr.stats.sac.t1)-100,int(tr.stats.sac.t2)+100])
            #window.append(int(tr.stats.sac.t2))
            s = tr
            egfs = sgfs[i].replace('FXZ','BXZ')
            d = read(egfs)[0]
            d.data = d.data/max(abs(s.data))
            [adj,misfit] = calculate_adjoint_source(d, s, config, window, adjoint_src, figure)
            out = Out_path + '/' + '30_80' + '/JJK.' + sgfs[i].strip().split('.')[0].split('_')[1] + '.' + files[j] +'.adj'
            Out = open(out,'w')
            misfit_out = Out_path + '/' +  '30_80' + '/' + files[j] + '_misfit.dat'
            #print(misfit['misfit'])
            os.system("echo %s >> %s" % (misfit,misfit_out))
            for k in range(len(t)):
                Out.write(str(t[k])+' '+str(adj[k])+'\n')
###

Data_path = '/project/li_chao/SEM3D/specfem3d/Test_model/ADJ/'
#dirs = ['05' , '100' , '150' , '180' , '30']
dirs = ['30']
fre_band = ['2-6','5-9','8-12']
W = {} ; wf = []
for j in range(len(fre_band)):
    w = []
    for i in range(len(dirs)):
        misfitfile =  Data_path + dirs[i] + '/SAC_' + fre_band[j] + '_misfit.dat'
        Misfit = np.loadtxt(misfitfile)
        Misfit_f = np.sum(Misfit)
        w.append(Misfit_f)
        print(Misfit_f)
    W[fre_band[j]] = np.sum(w)
    wf.append(W[fre_band[j]])
misfit_sum = np.sum(wf)

for i in range(len(fre_band)):
    W[fre_band[i]] = misfit_sum/(3*W[fre_band[i]])

## Write the final adj source ##
flag = True
if flag:
    Data_path = '/project/li_chao/SEM3D/specfem3d/Test_model/ADJ/'
    #dirs = ['05' , '100' , '150' , '180' , '30']
    dirs = ['30']
    fre_band = ['2-6','5-9','8-12']

    for j in range(len(dirs)):
        d_path = Data_path + dirs[j]
        #o_path = '/project/li_chao/SEM3D/specfem3d/Jiajika_solver/' + dirs[j]
        o_path = '/project/li_chao/SEM3D/specfem3d/misfit_test'
        os.chdir(d_path)
        STA_path = '/project/li_chao/SEM3D/specfem3d/SEM_Jiajika/DATA/STATIONS'
        with open(STA_path) as f :
            info = f.readlines()
    
        sta_id = [];sta_info={}
        for i in range(len(info)):
            splitInfos = info[i].strip().split()
            tmp = splitInfos[0]
            sta_info[tmp] = info[i]
            sta_id.append(tmp)
        
        for k in range(len(sta_id)):
            files = glob.glob("*.%s.*" % (sta_id[k]))
            if len(files) == 0 :
                continue
            T =[] ; Adj = {}
            for i in range(len(fre_band)):
                adj_file = 'JJK.' + sta_id[k] + '.SAC_' + fre_band[i] + '.adj'
                adj = [] 
                with open(files[i],'r') as f:
                    lines = f.readlines()
                for j in range(len(lines)):
                    #t.append(float(lines[j].strip().split()[0]))
                    adj.append(float(lines[j].strip().split()[1]))
                Adj[fre_band[i]] = adj
            
            print(sta_id[k],files)
            Adj_sum = []
            for m in range(23000):
                Adj_sum.append(W[fre_band[0]]*Adj[fre_band[0]][m] + W[fre_band[1]]*Adj[fre_band[1]][m] + W[fre_band[2]]*Adj[fre_band[2]][m])
            adj_out = open(o_path + '/SEM/JJK.' + sta_id[k] + '.FXZ.adj','w')
            
            time = 0.0001*np.arange(0,25400)
            for m in range(len(Adj_sum)):
                adj_out.write(str(time[m]) + ' ' + str(Adj_sum[m]) + '\n')
            for m in range(len(Adj_sum),25400):
                adj_out.write(str(time[m]) + ' ' + '0.0' + '\n')
            
            adj_outX = open(o_path + '/SEM/JJK.' + sta_id[k] + '.FXX.adj','w')
            adj_outY = open(o_path + '/SEM/JJK.' + sta_id[k] + '.FXY.adj','w')
            for m in range(25400):
                adj_outX.write(str(time[m]) + ' ' + '0.0' + '\n')
                adj_outY.write(str(time[m]) + ' ' + '0.0' + '\n')
            
        
