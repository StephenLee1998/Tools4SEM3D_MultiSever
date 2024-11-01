import numpy as np
import os
import glob
import math
import matplotlib as mpl
mpl.use('TkAgg')
from matplotlib import pyplot as plt
import sys
from Utils.Utils import bpfilter

def plot_seismogram(data_dir,period,color):
    files = glob.glob(data_dir+'/*')

    id = []

    for file in files:
        isplit = len(file.split('/'))
        itrace = file.split('/')[isplit-1].split('.')[1]
        
        data = np.loadtxt(file)
        T = data[:,0]
        D = data[:,1]
        plt.plot(T, 0.7 * D / max(abs(D)) + int(itrace), color)
        id.append(int(itrace))

    plt.xlim([min(T),max(T)])
    plt.ylim([min(id)-2,max(id)+2])


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print('You shloud input the data dir at least !')
    elif len(sys.argv) == 2:
        dir = sys.argv[1]
        
        plt.figure(figsize=(5,9))
        plot_seismogram(dir,1,'k')
        plt.show()
    elif len(sys.argv) == 3:
        dir1 = sys.argv[1]
        dir2 = sys.argv[2]
        plt.figure(figsize=(5,9))
        plot_seismogram(dir1,1,'k')
        plot_seismogram(dir2,1,'r')
        plt.show()