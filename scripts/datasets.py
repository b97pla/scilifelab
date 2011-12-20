#!/usr/bin/env python

import csv
import numpy as np
import matplotlib.pyplot as plt

dsets_lsh = csv.reader(open('datasets.csv', 'rU'))
dsets = csv.writer(open('datasets_toplot.csv', 'w'))

for row in dsets_lsh:
    dpref = row[0][-1]  #i.e 104...
    dsize = row[0][:-1] # ...G(igabytes)
    if dpref == 'T':
        dsize = int(float(dsize)*1024*1024*1024*1024)
    if dpref == 'G':
        dsize = int(float(dsize)*1024*1024*1024)
    elif dpref == 'M':
        dsize = int(float(dsize)*1024*1024)
    elif dpref == 'K':
        dsize = int(float(dsize)*1024)

    dsets.writerow([dsize, row[1]])

dsets = csv.reader(open('datasets_toplot.csv', 'rU'))

sizes = [int(row[0]) for row in dsets]

N = len(sizes)
ind = np.arange(N)      # the x locations for the groups
width = 0.35            # the width of the bars: can also be len(x) sequence

fig = plt.figure()
ax = fig.add_subplot(111)
num_values = 30
p1 = plt.bar(range(0, num_values), sizes[:num_values])

plt.ylabel('Sizes')
plt.title('Dataset sizes')
# all row[1]'s
# plt.xticks(ind+width/2., ('G1', 'G2'))

plt.show()
