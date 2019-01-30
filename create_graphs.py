import re
import sys
import pprint as pp
import numpy as np
import matplotlib.pyplot as plt

from os import listdir
from os.path import isfile, join
from datetime import datetime

regex_pattern = ".*"
if len(sys.argv) == 2:
    regex_pattern = sys.argv[1]

path = './results/stats/'
files = [f for f in listdir(path) if re.match(regex_pattern, f)]
files.sort()
pp.pprint(files)

plt.figure(figsize=(10, 20))

graph_path = './results/graphs/'
for i, f in enumerate(files):
    t = []
    x = []
    y = []
    with open(join(path, f), 'r') as fd:
        for line in fd:
            l = line.split()
            t.append(int(l[0]))
            x.append(float(l[1]))
            y.append(float(l[2]))

        r = [np.sqrt(x[i]*x[i] + y[i]*y[i]) for i in range(len(x))]

        plt.title(f)
        # plt.plot(t, r, label=f.split('.')[0])
        plt.plot(t, x, label='x')
        plt.plot(t, y, label='y')

lgd = plt.legend()
image_name = 'distances.png'
image_path = join(graph_path, image_name)
if isfile(image_path):
    now = datetime.now()
    [ name, extension ] = image_name.split('.')
    image_name = name + '_' + now.strftime("%Y%m%d_%H%M%S") + "." + extension
    image_path = join(graph_path, image_name)

plt.savefig(image_path, bbox_extra_artists=(lgd,), bbox_inches='tight')
print("\033[1;32mSaved to file " + image_path + "\n")
plt.cla()