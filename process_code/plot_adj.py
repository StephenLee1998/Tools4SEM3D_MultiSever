import os
import glob
from matplotlib import pyplot as plt

adj_path = '/project/li_chao/SEM3D/ADJ/05'
os.chdir(adj_path)

freq = ['2-6Hz' , '5-9Hz' , '8-12Hz']

sta = '197'
files = glob.glob('*%s*' % (sta))
print(files)
T = [];S = []
plt.figure(1,figsize = (9,16))
for i in range(len(files)):
    adj = files[i]
    t = [] ; s = []
    with open(adj) as f:
        lines = f.readlines()
    for j in range(len(lines)):
        t.append(float(lines[j].strip().split()[0]))
        s.append(float(lines[j].strip().split()[1]))
        T.append(t);S.append(s)
    plt.subplot(4,1,i+1)
    #plt.title(freq[i] + ' adjoint source')
    plt.xlim(0,2.3)
    plt.plot(t,s)

Adj_sum = []
for m in range(len(S[0])):
    Adj_sum.append(S[0][m] + S[1][m] + S[2][m])

plt.subplot(4,1,4)
#plt.title('Sum adjoint source')
plt.xlim(0,2.3)
plt.plot(t,Adj_sum)

plt.show()