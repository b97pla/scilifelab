#!/usr/bin/env python

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import itertools
from operator import itemgetter
import datetime
from collections import defaultdict
import numpy as np


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


def get_formatted_data():
    """Returns the data, processed a bit to be more regular.
    """
    # This of course need to be pointed to wherever the log file is available.
    with open("../../sizegraphing/datasets_sizes.log") as data_f:
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

    return data


def trim_data(grouped_data, trim=True):
    """Trims away the file sizes where preprocessing has occured.
    """
    # Dictionary with project names and the corresponding number of data
    # points which should be removed from the end.
    trimmings = {
        "0044_AD035EACXX": 1,
        "0050_AD037LACXX": 3,
        "0054_AC00H7ABXX": 1,
        "0055_BB06RCABXX": 2,
        "0157_B80NNMABXX": 1,
        "0158_A80NNNABXX": 1,
        "0175_B814PKABXX": 1,
        "0177_AB02VJABXX": 1,
        "0178_BB02UMABXX": 2,
        "0192_BB051UABXX": 2,
        "0201_BB06GNABXX": 1,
        "0202_AB0AA5ABXX": 1,
        "0203_BB039NABXX": 1,
        "0204_AB028LABXX": 2,
        "0205_BB047AABXX": 2,
        "0208_AB0B0KABXX": 1,
        "0209_BB0ADVABXX": 1,
        "0212_AB0B10ABXX": 1,
        "0214_B80LRFABXX": 2,
        "0215_B80J7PABXX": 1,
        "0219_BB0B29ABXX": 1,
        "0241_AC00RAABXX": 1,
        "0242_BB08CAABXX": 1,
        "0249_A816HRABXX": 1,
        "0250_B816J0ABXX": 2,
        "0251_A8158LABXX": 3,
        "0254_AB03BEABXX": 3,
        "0255_A81BF6ABXX": 2,
        "0255_BB04BAABXX": 2,
        "0256_B81CL7ABXX": 2,
        "0257_A81BHDABXX": 3,
        "0258_B819KCABXX": 3,
        "0261_A81BFKABXX": 2,
        "0262_B819K6ABXX": 2,
        "0264_B81BELABXX": 1,
        "0268_BB02J7ABXX": 2,
        "0269_AB02EVABXX": 1,
        "0270_BD035NACXX": 2,
        "0271_AB03T6ABXX": 1,
        "0273_AC00HFABXX": 2,
        "0274_BB028HABXX": 1,
        "0277_AB02J6ABXX": 1,
        "0278_BB0B0HABXX": 2,
        "0279_AC03E3ABXX": 1,
        "0280_BB0AC9ABXX": 1,
        "0281_AD0DYGABXX": 1,
        "0282_BC04B5ABXX": 2,
        "0283_AB06PBABXX": 2,
        "0284_BB06CNABXX": 2
    }

    sorted_data = []

    for group in grouped_data:
        sizes = []
        for line in group:
            sizes.append(line[0])
        sizes = sorted(sizes)
        if trim:
            remove = trimmings.get(line[-1], None)
            if remove is not None:
                sizes = sizes[:-remove]

        sorted_data.append([sizes, [line[1]], line[2], line[3]])

    return sorted_data


def file_size_over_time_plot(plot_type="bar", given_machine=None, trim=True):
    """Makes a plot of the amount of space the project files took up at
    specific times.
    """
    # Get the data and structure it

    data = get_formatted_data()

    data = sorted(data, key=itemgetter(3))
    grouped_data = [list(g) for k, g in itertools.groupby(data, itemgetter(3))]

    sorted_data = trim_data(grouped_data, trim)

    # Add the additional days to the project tp use when plotting
    for project in sorted_data:
        for i in range(1, len(project[0])):
            project[1].append(project[1][0] + datetime.timedelta(days=i))

    # Plot

    fig = plt.figure()
    ax = fig.add_subplot(111)
    stack_heights = defaultdict(int)

    bytes_per_gigabyte = 1024. ** 3

    machines = set([])
    for group in sorted_data:
        machines.add(group[2])
    mc_val = dict(zip(machines, [float(i) / len(machines) for i in range(len(machines))]))
    if given_machine is not None:
        filtered_grouped_data = [group for group in sorted_data if group[2] == given_machine]
        sorted_data = filtered_grouped_data
    for i, group in enumerate(sorted_data):
        clr = cm.hsv(mc_val[group[2]])
        sizes = [s / bytes_per_gigabyte for s in group[0]]
        days = group[1]
        heights = []
        for day in days:
            heights.append(stack_heights[day])

        if plot_type == 'bar':
            ax.bar(days, sizes, bottom=heights, color=clr, width=1.0, lw=0.1)

        for size, day in zip(sizes, days):
            stack_heights[day] += size

    plt.ylabel('Size (GB)')
    figure_title = 'Dataset sizes'
    if given_machine is not None:
        figure_title += ' - Machine ' + given_machine
    ax.set_title(figure_title)

    if plot_type == 'plot':
        sorted_stacked = sorted(stack_heights.items(), key=itemgetter(0))
        days = [l[0] for l in sorted_stacked]
        stacks = [l[1] for l in sorted_stacked]
        zip_of_data = dict(zip(days, stacks))
        daylist = [days[0] + datetime.timedelta(days=x) for x in range(0, (days[-1] - days[0]).days)]
        stacklist = []
        for day in daylist:
            stacklist.append(zip_of_data.get(day, 0))

        plt.plot(daylist, stacklist)

    # fig.savefig("../sizegraphing/stacked_sizes.pdf")
    plt.show()


def get_average_size_per_project(plot=False):
    """Returs the average size of a project when it is finished. In Gigabytes.
    """
    data = get_formatted_data()
    data = sorted(data, key=itemgetter(-1))
    grouped_data = [list(g) for k, g in itertools.groupby(data, itemgetter(-1))]
    sorted_data = trim_data(grouped_data)
    final_sizes = []
    for project in sorted_data:
        final_sizes.append(project[0][-1])
    final_sizes = np.array([float(s) for s in final_sizes if s > 0. * 1024 ** 3])

    average = np.mean(final_sizes / (1024 ** 3))

    if plot:
        plt.scatter(final_sizes / (1024 ** 3), [0] * len(final_sizes), marker='+')
        plt.plot(average, 0, 'ro')
        plt.show()

    return average
    # Gives odd results, might need to consider the following cases:
    # * Single end
    # * Paired end
    # * Failed run
    # * Version 3 flowcell


def numer_of_projects_growth_plot():
    """Plots the growth of the number of projects over time.
    """
    data = get_formatted_data()
    data = sorted(data, key=itemgetter(-1))
    grouped_data = [list(g) for k, g in itertools.groupby(data, itemgetter(-1))]
    sorted_data = trim_data(grouped_data)
    days = [li[1][0] for li in sorted_data]
    days = sorted(days)
    fig = plt.figure()
    ax = fig.add_subplot("111")
    ax.plot(days, range(len(days)))
    ax.set_title("Total number of projects over time")
    plt.show()


## TODO: Make a plot which is easy to grasp at a glance.
## TODO: Get average size per project.
#! TODO: Get total number of basepairs generated since start.
#! TODO: Get total number of indexes. (Requires access to samplesheet files!)
#^ (TODO: Also make a graph of 'number of basepairs over time.')
## (TODO: Graph total number of projects over time.)
