import numpy as np
import os
from numpy import diff
from scipy.signal import butter, lfilter
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import freqz
import obspy
from scipy.signal import hilbert
from scipy import interpolate
from obspy.core import Trace
from obspy.io.sac import SACTrace

def _search(id1,id2):
    filename = sac_path + '/VV_'+id1+'-'+id2+'_14d.dat.SAC'
    if os.path.exists(filename):
        return filename
    else:
        filename = sac_path + '/VV_'+id2+'-'+id1+'_14d.dat.SAC'
        if os.path.exists(filename):
            return filename
        else:
            return False

def _stack(st):
    half_l = int(np.floor(len(st)/2))
    st_on = st[0:5001]
    st_off = st[5000:10001]
    st_on_re = st_on[::-1]
    st_out = st_on_re + st_off
    return st_out

def _readsac(sacname):
    st = obspy.read(sacname)
    tr = st[0]
    tr_filt = tr.copy()
    dist = tr.stats.sac.dist
    #tr_filt.filter('bandpass',freqmin=4.0,freqmax=6.0,corners=2, zerophase=True)
    #filted = tr_filt
    #tr_stack = _stack(filted.data)[0:200]
    tr_stack = _stack(tr.data)[0:232]
    dx = 0.0001
    tr_egf = -diff(tr_stack)/dx
    #tr_ngf = np.real(hilbert(tr_stack))
    #tr_ngf = scipy.fftpack.hilbert(tr_stack)
    return tr_egf,dist

def _readsemd(sacname):
    data = np.fromfile(sacname,dtype = 'float32')
    return data

def read_EGFs(sacs):
    ### read the EGFs and interplate ###
    [tr_stack,dist] = _readsac(sacs)
    delta1 = 0.01
    T1 = delta1*np.arange(len(tr_stack))
    f = interpolate.interp1d(T1,tr_stack,kind='linear')
    T2 = np.linspace(0,2.2999,100*(len(tr_stack)-1))
    tr_new = f(T2)
    tr_new = tr_new.astype(np.float32)
    return tr_new,dist

def _filt(data,fmin,fmax):
    tr = obspy.read()[0]
    tr.data = data
    tr.stats.sampling_rate = 10000
    filted = tr.filter('bandpass',freqmin=fmin,freqmax=fmax,corners=4, zerophase=True)
    return filted.data

def _plot(sgfs_data,egfs_data):
    plt.figure(1,figsize=(16,9))
    tr1 = _filt(sgfs_data,2,6) 
    tr2 = _filt(sgfs_data,5,9)
    tr3 = _filt(sgfs_data,8,12)
    tr4 = _filt(egfs_data,2,6) 
    tr5 = _filt(egfs_data,5,9)
    tr6 = _filt(egfs_data,8,12)
    plt.subplot(311)
    plt.plot(T,tr1/max(abs(tr1)),color='black',linewidth=1)
    plt.plot(T,tr4/max(abs(tr4)),color='red',linewidth=1)
    plt.xlim(0,2.3)
    plt.subplot(312)
    plt.plot(T,tr2/max(abs(tr2))+1,color='black',linewidth=1)
    plt.plot(T,tr5/max(abs(tr5))+1,color='red',linewidth=1)
    plt.xlim(0,2.3)
    plt.subplot(313)
    p1,=plt.plot(T,tr3/max(abs(tr3))+2,color='black',linewidth=1)
    p2,=plt.plot(T,tr6/max(abs(tr6))+2,color='red',linewidth=1)
    plt.xlim(0,2.3)
    plt.legend([p1,p2],["SGFs","EGFs"],loc = 'lower right')
    plt.show()

def write_sac(data,fmin,fmax,semd,sour_id,dist):
    tr_data = _filt(data,fmin,fmax)
    sacfile = Trace()
    sacfile.data = tr_data
    sac = SACTrace.from_obspy_trace(sacfile)
    sac_data = sac.data
    sac.dist=dist
    outid = str(fmin) + '-' + str(fmax)
    out = out_path[outid] + '/' + str(sour_id) + '_' + str(semd) + '_' + str(fmin) + '_' + str(fmax) + '.sac'
    sac.write(out)
    print(out)
    

dirs = ['test_160']

sac_path = '/project/li_chao/SEM3D/SAC_Jiajika'
file_path = '/project/li_chao/SEM3D/specfem3d/Jiajika_solver/test_30/OUTPUT_FILES'
STA_path = '/project/li_chao/SEM3D/specfem3d/SEM_Jiajika/DATA/STATIONS'
out_path = {}
out_path['2-6'] = '/project/li_chao/SEM3D/specfem3d/misfit_test/Window_Picking_160/SAC_2-6'
out_path['5-9'] = '/project/li_chao/SEM3D/specfem3d/misfit_test/Window_Picking_160/SAC_5-9'
out_path['8-12'] = '/project/li_chao/SEM3D/specfem3d/misfit_test/Window_Picking_160/SAC_8-12'

#os.system('cd /project/li_chao/SEM3D/specfem3d/misfit_test/Window_Picking_40')
#os.system('mkdir SAC_2-6 SAC_5-9 SAC_8-12')

with open(STA_path) as f :
    info = f.readlines()

sta_id = [];sta_info={}
for i in range(len(info)):
    splitInfos = info[i].strip().split()
    tmp = splitInfos[0]
    sta_info[tmp] = info[i]
    sta_id.append(tmp)

for j in range(len(dirs)):
    sour_id = '30'
    file_path = '/project/li_chao/SEM3D/specfem3d/Jiajika_solver/'+ dirs[j] +'/OUTPUT_FILES'
    #out = open('./STATIONS_ADJOINT','w') 
    df = 10000
    
    #plt.figure(1,figsize=(10,5))
    #/project/li_chao/SEM3D/SAC_Jiajika/VV_142-05_14d.dat.SAC"    
    for i in range(len(sta_id)):
        ## read the EGFs and interplate ##
        sacs = _search(sour_id,sta_id[i])
        if sacs == False :
            continue
        [egfs_data,dist] = read_EGFs(sacs)
        #out.write(sta_info[sta_id[i]])
        ## read the SGFs ##
        semd = file_path + '/JJK.'+sta_id[i]+'.FXZ.semd'
        egfs = str(sta_id[i]) + '.BXZ'
        sgfs_data = _readsemd(semd)
    
        T = 0.0001*np.arange(len(sgfs_data))
        #_plot(sgfs_data,egfs_data)
        #print(sgfs_data)
        if dist > 1.6:
            write_sac(sgfs_data,2,6,semd.split('JJK.')[1].split('.semd')[0],sour_id,dist)
            write_sac(sgfs_data,5,9,semd.split('JJK.')[1].split('.semd')[0],sour_id,dist)
            write_sac(sgfs_data,8,12,semd.split('JJK.')[1].split('.semd')[0],sour_id,dist)
            write_sac(egfs_data,2,6,egfs,sour_id,dist)
            write_sac(egfs_data,5,9,egfs,sour_id,dist)
            write_sac(egfs_data,8,12,egfs,sour_id,dist)