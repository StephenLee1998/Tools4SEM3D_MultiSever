from matplotlib import pyplot as plt
import os
import glob

data_path = '/project/li_chao/SEM3D/ADJ/05'
os.chdir(data_path)
files = glob.glob("JJK.153*")

T =[] ; Adj = []
for i in range(len(files)):
    t = [] ; adj = [] 
    with open(files[i],'r') as f:
        lines = f.readlines()
    for j in range(len(lines)):
        t.append(float(lines[j].strip().split()[0]))
        adj.append(float(lines[j].strip().split()[1]))
    T.append(t)
    Adj.append(adj)

Adj_sum = []
for i in range(len(Adj[0])):
    Adj_sum.append(Adj[0][i] + Adj[1][i] + Adj[2][i])

plt.figure(1,figsize = (16,9))

plt.subplot(411);plt.plot(t,Adj[0]);plt.xlim(0,2.3);plt.title('2-6Hz')
plt.subplot(412);plt.plot(t,Adj[1]);plt.xlim(0,2.3);plt.title('5-9Hz')
plt.subplot(413);plt.plot(t,Adj[2]);plt.xlim(0,2.3);plt.title('8-12Hz')
plt.subplot(414);plt.plot(t,Adj_sum),plt.xlim(0,2.3);plt.title('Sum_adj_source')
plt.show()