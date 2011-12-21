#!/usr/bin/env python

import matplotlib.pyplot as plt


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

N = len(sizes)

fig = plt.figure()
ax = fig.add_subplot(111)
num_values = 75
ax.bar(range(num_values), sizes[:num_values])

plt.ylabel('Sizes')
ax.set_title('Dataset sizes')

plt.show()
