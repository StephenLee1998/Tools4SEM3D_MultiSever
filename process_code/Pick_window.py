import numpy as np
import os
from matplotlib import pyplot as plt
import obspy
from scipy.signal import hilbert
from scipy import interpolate
import obspy.signal.trigger as ost
from numpy import diff


sac_path = '/project/li_chao/SEM3D/SAC_Jiajika'
file_path = '/project/li_chao/SEM3D/specfem3d/Jiajika_solver/old_30/OUTPUT_FILES'
STA_path = '/project/li_chao/SEM3D/specfem3d/SEM_Jiajika/DATA/STATIONS'

with open(STA_path) as f :
    info = f.readlines()

sta_id = [];sta_info={}
for i in range(len(info)):
    splitInfos = info[i].strip().split()
    tmp = splitInfos[0]
    sta_info[tmp] = info[i]
    sta_id.append(tmp)

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
    tr_filt.filter('bandpass',freqmin=8.0,freqmax=12.0,corners=2, zerophase=True)
    filted = tr_filt
    #tr_stack = _stack(filted.data)[0:200]
    tr_stack = _stack(filted.data)[0:232]
    dx = 0.0001
    tr_egf = -diff(tr_stack)/dx
    #tr_ngf = np.real(hilbert(tr_stack))
    #tr_ngf = scipy.fftpack.hilbert(tr_stack)
    return tr_egf,dist

def _readsemd(sacname):
    data = np.fromfile(sacname,dtype = 'float32')
    return data

def _search_cft(data,T,cut1,cut2):
    for i in range(len(data)):
        if data[i] > cut1 :
            id1 = i
            break
    for i in range(id1,len(data)):
        if data[i] < cut2:
            id2 = i
            break
        else:
            id2 =19999
    return T[id1],T[id2]

def plot_single(T2,tr_new,tr_semd):
    norm = max(abs(tr_new))/max(abs(tr_semd))
    plt.plot(T2,tr_new,color = 'black')
    plt.plot(T2,norm*tr_semd,color = 'red')
    plt.xlim(0,2.3)

## PLOT SACS ##
sour_id = '30'
out = open('./STATIONS_ADJOINT','w') 
df = 10000
#/project/li_chao/SEM3D/SAC_Jiajika/VV_142-05_14d.dat.SAC"    
for i in range(117,197):
    sacs = _search(sour_id,sta_id[i])
    if sacs == False :
        continue
    [tr_egf,dist] = _readsac(sacs)
    out.write(sta_info[sta_id[i]])
    delta1 = 0.01
    T1 = delta1*np.arange(len(tr_egf))
    f = interpolate.interp1d(T1,tr_egf,kind='linear')
    T2 = np.arange(0,2.3,0.0001)
    tr_new = f(T2)
    tr_new = tr_new.astype(np.float32)
    tr_zeros = np.zeros(len(tr_new),dtype='float32')

    semd = file_path + '/JJK.'+sta_id[i]+'.FXZ.semd'
    semd1 = '/project/li_chao/SEM3D/specfem3d/Jiajika_solver/r80_30/OUTPUT_FILES' + '/JJK.'+sta_id[i]+'.FXZ.semd'
    semd2 = '/project/li_chao/SEM3D/specfem3d/Jiajika_solver/r160_30/OUTPUT_FILES' + '/JJK.'+sta_id[i]+'.FXZ.semd'
    semd3 = '/project/li_chao/SEM3D/specfem3d/Jiajika_solver/30/OUTPUT_FILES' + '/JJK.'+sta_id[i]+'.FXZ.semd'
    semd4 = '/project/li_chao/SEM3D/specfem3d/Jiajika_solver/r40_30/OUTPUT_FILES' + '/JJK.'+sta_id[i]+'.FXZ.semd'
    tr_semd = _readsemd(semd)
    tr_semd1 = _readsemd(semd1)
    tr_semd2 = _readsemd(semd2)
    tr_semd3 = _readsemd(semd3)
    tr_semd4 = _readsemd(semd4)

    #norm = max(abs(tr_new))/max(abs(tr_semd))
    #tr_semd = tr_semd * norm
    plt.figure(1,figsize=(9,16))
    plt.plot(T2,tr_new/max(abs(tr_new))/50+dist,color = 'black')
    plt.plot(T2,tr_semd/max(abs(tr_semd))/50+dist,color = 'red')
    #print(len(tr_semd),len(tr_new))
    #plt.subplot(5,1,1)
    #plot_single(T2,tr_new,tr_semd)
    #plt.subplot(5,1,2)
    #plot_single(T2,tr_new,tr_semd4)
    #plt.subplot(5,1,3)
    #plot_single(T2,tr_new,tr_semd1)
    #plt.subplot(5,1,4)
    #plot_single(T2,tr_new,tr_semd2)
    #plt.subplot(5,1,5)
    #plot_single(T2,tr_new,tr_semd3)
    
plt.xlim(0,2)
plt.show()