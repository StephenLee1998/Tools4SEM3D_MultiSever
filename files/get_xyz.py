import numpy as np

f1=open('./DATA/STATIONS','r')
lines1=f1.readlines()
f1.close()

f2=open('./stations','r')
lines2=f2.readlines()
f2.close()

nlen = len(lines1)

for i in np.arange(nlen):
    print(lines1[i].strip().split()[0], lines2[i].strip().split()[0], lines2[i].strip().split()[1], lines2[i].strip().split()[2])
