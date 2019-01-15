import re
import pprint as pp
import numpy as np
import matplotlib.pyplot as plt

from os import listdir
from os.path import isfile, join

path = './results/stats/'
files = [f for f in listdir(path)]
files.sort()
pp.pprint(files)

plt.figure(figsize=(10, 20))

path2 = './results/graphs/'
x_all = np.zeros((len(files), 199))
y_all = np.zeros((len(files), 199))
z_all = np.zeros((len(files), 199))
t = []
for i, f in enumerate(files):
    t = []
    x = []
    y = []
    z = []
    with open(join(path, f), 'r') as fd:
        for line in fd:
            l = line.split()
            t.append(int(l[0]))
            x.append(float(l[1]))
            y.append(float(l[2]))
            z.append(float(l[3]))

    x_all[i] = x[1:]
    y_all[i] = y[1:]
    z_all[i] = z[1:]
    
#     plt.plot(t, x)
#     plt.xlabel('t')
#     plt.ylabel('x')
#     plt.savefig(join(path2, f.split('.')[0] + '_x.png'), bbox_inches='tight')
#     plt.cla()

#     plt.plot(t, y)
#     plt.xlabel('t')
#     plt.ylabel('y')
#     plt.savefig(join(path2, f.split('.')[0] + '_y.png'), bbox_inches='tight')
#     plt.cla()

#     plt.plot(t, z)
#     plt.xlabel('t')
#     plt.ylabel('z')
#     plt.savefig(join(path2, f.split('.')[0] + '_z.png'), bbox_inches='tight')
#     plt.cla()
    
    plt.plot(t, x, label='x')
    plt.plot(t, y, label='y')
    plt.plot(t, z, label='z')
    plt.xlabel('t')
    lgd = plt.legend()
    plt.savefig(join(path2, f.split('.')[0] + '_xyz.png'), bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.cla()

x_avg = np.average(x_all, axis=1)
plt.plot(x_avg)
plt.savefig(join(path2, 'avg_x.png'), bbox_inches='tight')
plt.cla()

y_avg = np.average(y_all, axis=1)
plt.plot(y_avg)
plt.savefig(join(path2, 'avg_y.png'), bbox_inches='tight')
plt.cla()

z_avg = np.average(z_all, axis=1)
plt.plot(z_avg)
plt.savefig(join(path2, 'avg_z.png'), bbox_inches='tight')
plt.cla()

