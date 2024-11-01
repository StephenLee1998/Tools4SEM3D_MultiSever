#!/usr/bin/python
import numpy as np
import os
import glob
from obspy import read 

Out_path = '/project/li_chao/SEM3D/specfem3d/SEM_Jiajika/OUTPUT_FILES'

os.chdir(Out_path)
#sufile = glob.glob("*SU")
#
#for i in range(len(sufile)):
#    print(sufile)
#    data=read(sufile[i])
#    for j in range(len(data)):
#        data[j].data = data[j].data[2400:22400]
#    data.write(sufile[i],format='SH_ASC')

semd = glob.glob("*semd")

for i in range(len(semd)):
    data = np.fromfile(semd[i],dtype='float32')
    os.system("cp %s OLD" % (semd[i]))
    data = data[2400:25400]
    data.tofile(semd[i])
