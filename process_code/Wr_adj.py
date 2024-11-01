#!/usr/bin/python
import numpy as np
import os
import math
from matplotlib import pyplot as plt
import obspy
from scipy.signal import hilbert
from scipy import interpolate
import obspy.signal.trigger as ost

sac_path = '/project/li_chao/SEM3D/SAC_Jiajika'
file_path = '/project/li_chao/SEM3D/specfem3d/SEM_Jiajika/OUTPUT_FILES'
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
    #tr_filt.filter('bandpass',freqmin=4.0,freqmax=6.0,corners=2, zerophase=True)
    #filted = tr_filt
    #tr_stack = _stack(filted.data)[0:200]
    tr_stack = _stack(tr.data)[0:231]
    tr_ngf = np.real(hilbert(tr_stack))
    #tr_ngf = scipy.fftpack.hilbert(tr_stack)
    return tr_stack,dist

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
    plt.subplot(2,1,1)
    plt.plot(T2,tr_new,color = 'black')
    plt.plot(T2,tr_semd,color = 'red')
    plt.plot(t1*np.ones(20),np.arange(-1,1,0.1),color = 'red')
    plt.plot(t2*np.ones(20),np.arange(-1,1,0.1),color = 'blue')
    plt.subplot(2,1,2)
    plt.plot(T2,cft1,color = 'black')
    plt.plot(T2,cft2,color = 'red')
    plt.plot(np.arange(0,2,0.1),-0.0*np.ones(20),'--',color = 'red')
    plt.plot(np.arange(0,2,0.1),-0.2*np.ones(20),'--',color = 'blue')

def cos_window(t,t1,t2):
    win = np.zeros(len(T2))
    dt = math.pi/0.1
    for i in range(len(t)):
        if t[i] < t1-0.1 or t[i] > t2+0.1:
            continue
        elif t[i] > t1 and t[i] < t2:
            win[i]= 1.0
        elif t[i] >= t1-0.1 and t[i] <= t1:
            win[i] = 0.5*math.cos(dt*(t[i]-t1)) + 0.5
        elif t[i] <= t2+0.1 and t[i] >= t2:
            win[i] = 0.5*math.cos(dt*(t[i]-t2)) + 0.5
    return win



## PLOT SACS ##
sour_id = '05'
out = open('./STATIONS_ADJOINT','w') 
df = 10000

plt.figure(1,figsize=(9,16))
#/project/li_chao/SEM3D/SAC_Jiajika/VV_142-05_14d.dat.SAC"    
for i in range(len(sta_id)):
    sacs = _search(sour_id,sta_id[i])
    if sacs == False :
        continue
    [tr_stack,dist] = _readsac(sacs)
    out.write(sta_info[sta_id[i]])
    delta1 = 0.01
    T1 = delta1*np.arange(len(tr_stack))
    f = interpolate.interp1d(T1,tr_stack,kind='linear')
    T2 = np.linspace(0,2.2999,100*(len(tr_stack)-1))
    tr_new = f(T2)
    tr_new = tr_new.astype(np.float32)
    tr_zeros = np.zeros(len(tr_new),dtype='float32')

    semd = file_path + '/JJK.'+sta_id[i]+'.FXZ.semd'
    tr_semd = _readsemd(semd)

    norm = max(abs(tr_new))/max(abs(tr_semd))
    tr_semd = tr_semd * norm
    
    #adj = tr_new - tr_semd

    cft1 = ost.z_detect(tr_new,int(0.3*df))
    cft2 = ost.z_detect(tr_semd,int(0.3*df))

    [t1,t2] = _search_cft(cft2,T2,-0.0,-0.2)
    t1 = t1 - 0.1
    if dist < 0.75:
        t1 = t2 - 0.45
    print(t1,t2,t2-t1)

    win = cos_window(T2,t1,t2)

    tr_obs = tr_new * win
    tr_sim = tr_semd * win

    adj = tr_sim - tr_obs
    plt.plot(T2,adj+int(sta_id[i]),color = 'black')
    #plt.plot(T2,tr_new/100+float(dist),color = 'red')
    #plt.plot(T2,tr_obs/100+float(dist),color = 'black')
    #plt.plot(T2,win/100+float(dist),color = 'black')
    #plt.plot(t1*np.ones(20),np.arange(-1,1,0.1)/100+dist,color = 'red')
    #plt.plot(t2*np.ones(20),np.arange(-1,1,0.1)/100+dist,color = 'blue')
    

    filenameZ =  '/project/li_chao/SEM3D/specfem3d/Forward_JIAJIKA/SEM/'"JJK." + sta_id[i] + '.FXZ.adj'
    filenameX =  '/project/li_chao/SEM3D/specfem3d/Forward_JIAJIKA/SEM/'"JJK." + sta_id[i] + '.FXX.adj'
    filenameY =  '/project/li_chao/SEM3D/specfem3d/Forward_JIAJIKA/SEM/'"JJK." + sta_id[i] + '.FXY.adj'
    out1 = open(filenameZ,'w')
    out2 = open(filenameX,'w')
    out3 = open(filenameY,'w')
    for i in range(len(T2)):
        out1.write(str(T2[i])+' '+str(adj[i])+'\n')
        out2.write(str(T2[i])+' '+str(tr_zeros[i])+'\n')
        out3.write(str(T2[i])+' '+str(tr_zeros[i])+'\n')
    #tr_new.tofile(filenameZ)
    #tr_zeros.tofile(filenameX)
    #tr_zeros.tofile(filenameY)

#plt.ylim(0,3.8)
plt.xlim(0,2.3)
plt.show()