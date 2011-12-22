#!/usr/bin/env python

import matplotlib.pyplot as plt
import itertools
from operator import itemgetter
import datetime
from collections import defaultdict


def convert_human_readable_to_numbers(file_size):
    """In data, convert letters which indicate file sizes to actual numbers.
    """
    dpref = file_size[-1]
    if dpref not in ["K", "M", "G", "T"]:
        return 0
    file_size = file_size[:-1]
    if dpref == 'T':
        file_size = int(float(file_size) * 1024 * 1024 * 1024 * 1024)
    elif dpref == 'G':
        file_size = int(float(file_size) * 1024 * 1024 * 1024)
    elif dpref == 'M':
        file_size = int(float(file_size) * 1024 * 1024)
    elif dpref == 'K':
        file_size = int(float(file_size) * 1024)

    return file_size

# Get the data and structure it

with open("../sizegraphing/datasets_sizes.log") as data_f:
    all_data = data_f.read()

data_lines = all_data.splitlines()
data = [line.split("\t") for line in data_lines]
for line in data:
    line[-1] = line[-1].split("/")[-1]

data = [line for line in data if "." not in line[-1]]
sizes = [int(convert_human_readable_to_numbers(line[0])) for line in data]
data = [line[-1].split("_") for line in data]

for i in range(len(sizes)):
    data[i].insert(0, sizes[i])

data = [line for line in data if len(line) == 5]

for i in range(len(data)):
    data[i][3:5] = ["_".join(data[i][3:5])]
    data[i][1] = datetime.datetime.strptime(data[i][1], "%y%m%d")

for line in data:
    line[2] = line[2].upper()
    line[3] = line[3].upper()

data = sorted(data, key=itemgetter(3))
grouped_data = [list(g) for k, g in itertools.groupby(data, itemgetter(3))]

sorted_grouped_data = []
for group in grouped_data:
    sorted_grouped_data.append(sorted(group, key=itemgetter(0)))

# Modify data for plotting (increment dates with the days past since creation)

for i in range(len(sorted_grouped_data)):
    for j in range(len(sorted_grouped_data[i])):
        sorted_grouped_data[i][j][1] += datetime.timedelta(days=j)

# Plot

fig = plt.figure()
ax = fig.add_subplot(111)
colors = ["b", "g", "r", "c", "m", "y"]
i = 0
stack_heights = defaultdict(int)
for group in sorted_grouped_data:
    clr = colors[i % 6]
    i += 1
    sizes = [line[0] for line in group]
    days = [line[1] for line in group]
    heights = []
    for day in days:
        heights.append(stack_heights[day])

    ax.bar(days, sizes, bottom=heights, color=clr, width=1.0, lw=0.01)

    for size, day in zip(sizes, days):
        stack_heights[day] += size

plt.ylabel('Sizes')
ax.set_title('Dataset sizes')

fig.savefig("../sizegraphing/stacked_sizes.pdf")
