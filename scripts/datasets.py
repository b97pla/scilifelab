#!/usr/bin/env python

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import itertools
from operator import itemgetter
import datetime
from collections import defaultdict
import numpy as np
import StringIO
import json

# Dictionary with project names and the corresponding number of data
# points which should be removed from the end.
TRIMMINGS = {
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


def convert_human_readable_to_numbers(file_size):
    """In data, convert letters which indicate file sizes to actual numbers.
    """
    dpref = file_size[-1]
    if dpref not in ["K", "M", "G", "T"]:
        return 0

    file_size = file_size[:-1]
    if dpref == 'G':
        file_size = int(float(file_size) * 1024 * 1024 * 1024)
    elif dpref == 'T':
        file_size = int(float(file_size) * 1024 * 1024 * 1024 * 1024)
    elif dpref == 'M':
        file_size = int(float(file_size) * 1024 * 1024)
    elif dpref == 'K':
        file_size = int(float(file_size) * 1024)

    return file_size


def get_formatted_data():
    """Returns the data, processed a bit to be more regular.
    """
    # This of course need to be pointed to wherever the log file is available.
    with open("/Users/valentinesvensson/Documents/sizegraphing/datasets_sizes.log") as data_f:
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
    trimmings = TRIMMINGS

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


def get_list_of_machines():
    """Returns a list of machine serial numbers found in the logs.
    """
    data = get_formatted_data()
    data = [list(g) for k, g in itertools.groupby(data, itemgetter(3))]
    data = trim_data(data)
    machines = set([])
    for proj in data:
        machines.add(proj[2])
    return list(machines)


def dictify(data):
    ddata = []
    for data_list in data:
        ddata.append(list_to_dict(data_list))

    return ddata


def convert_to_json(trim=True):

    trimmings = TRIMMINGS

    data = get_formatted_data()
    ddata = dictify(data)
    get_project = lambda d: d["project"]
    ddata = sorted(ddata, key=get_project)
    grouped_data = [list(g) for k, g in itertools.groupby(ddata, get_project)]
    data_with_fixed_dates = []
    for group in grouped_data:
        group = sorted(group, key=lambda d: d["size"])
        init_date = group[0]["date"]
        remove = trimmings.get(group[0]["project"], None)
        if trim and remove is not None:
            group = group[:-remove]
        for i, log in enumerate(group):
            new_date = init_date + datetime.timedelta(i)
            log["date"] = new_date
            data_with_fixed_dates.append(log)

    dthandler = lambda obj: obj.isoformat() if isinstance(obj, datetime.datetime) else None
    with open("sizelog.json", "w") as handle:
        handle.write(json.dumps(data_with_fixed_dates, default=dthandler))


def file_size_over_time_plot(plot_type="bar", given_machine=None, trim=True, \
start_date=None, return_svg_data=False):
    """Makes a plot of the amount of space the project files took up at
    specific times.
    """
    # Get the data and structure it

    data = get_formatted_data()

    data = sorted(data, key=itemgetter(3))
    grouped_data = [list(g) for k, g in itertools.groupby(data, itemgetter(3))]

    sorted_data = trim_data(grouped_data, trim)
    if start_date is not None:
        sorted_data = [p for p in sorted_data if p[1][0] > start_date]

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
    if return_svg_data:
        fig.set_size_inches(10, 2)
        svg_data = get_svg_data(fig)
        return svg_data
    else:
        plt.show()


def get_svg_data(fig):
    imgdata = StringIO.StringIO()
    fig.savefig(imgdata, format='svg')
    imgdata.seek(0)

    svg_dta = imgdata.buf

    return svg_dta


def list_to_dict(data_list):
    """Convert to a dict.
    """
    data_dict = {
    "size": data_list[0],
    "date": data_list[1],
    "machine": data_list[2],
    "project": data_list[3]
    }

    return data_dict


def get_average_size_per_project(plot=False, return_data=False, start_date=None):
    """Returs the average size of a project when it is finished. In Gigabytes.
    """
    data = get_formatted_data()
    data = sorted(data, key=itemgetter(-1))
    grouped_data = [list(g) for k, g in itertools.groupby(data, itemgetter(-1))]
    sorted_data = trim_data(grouped_data)
    special_cases = ["0201_BB06GNABXX", "0207_BC00J3ABXX", \
    "0226_AC019MACXX", "0236_AC043HACXX", "0227_BD030GACXX"]
    sorted_data = [p for p in sorted_data if p[-1] not in special_cases]
    if start_date is not None:
        sorted_data = [p for p in sorted_data if p[1][0] > start_date]
    B_data = [p for p in sorted_data if p[-1][-3] == "B"]
    C_data = [p for p in sorted_data if p[-1][-3] == "C"]
    B_final_sizes = []
    C_final_sizes = []
    for project in B_data:
        B_final_sizes.append(project[0][-1])
    for project in C_data:
        C_final_sizes.append(project[0][-1])

    B_final_sizes = np.array([float(s) for s in B_final_sizes if 10000. * 1024 ** 3 > s > 100. * 1024 ** 3])
    C_final_sizes = np.array([float(s) for s in C_final_sizes if s > 100. * 1024 ** 3])

    average = [np.mean(B_final_sizes / (1024 ** 3)), np.mean(C_final_sizes / (1024 ** 3))]

    if plot:
        plt.plot(B_final_sizes / (1024 ** 3), [0] * len(B_final_sizes), 'b+')
        plt.plot(C_final_sizes / (1024 ** 3), [1] * len(C_final_sizes), 'g+')
        plt.plot(average, [0, 1], 'ro')
        x1, x2, y1, y2 = plt.axis()
        plt.axis((x1, x2, -1, 2))
        plt.yticks([0, 1], ('v1/v1.5', 'v3'))
        plt.show()

    if return_data:
        return (B_data, C_data)

    return average

    # Gives odd results, might need to consider the following cases:
    # * Single end
    # * Paired end
    # * Failed run
    # * Version 3 flowcell

    # According to Anna, in the beginning of using version 3 flowcell there
    # were some troubles getting as much DNA in a run as specified.
    # Looking only at Q4 2011 should give results according to spec and might
    # be a good sanity check to see if the data analysis is good.

    ## Update: the final sizes are equally distributed when looking at only
    ## later dates, there are less outliers though.

    # (Also, during the v1/v1.5 times, it was possible to put twice
    # as much as specified on the flowcell, and this was something one
    # generally tried to do; this could be one reason there is twice as much
    # data than expected at the end of a v1(.5) run. (But not 3-7 times as
    # much data, as some outliers are.))


def number_of_projects_growth_plot():
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
    ax.set_title("Total number of projects over time: %i" % (len(days),))
    # plt.show()
    plt.savefig("/Users/valentinesvensson/Desktop/projects_over_time.pdf")


def get_total_number_of_gigabasepairs(plot=False):
    data = get_formatted_data()
    data = sorted(data, key=itemgetter(-1))
    grouped_data = [list(g) for k, g in itertools.groupby(data, itemgetter(-1))]
    sorted_data = trim_data(grouped_data)
    sorted_data = sorted(sorted_data, key=lambda flowcell: flowcell[1][0] + datetime.timedelta(days=len(flowcell[0])))

    # The data points in sorted_data are in terms of flowcells.
    # (Each run consists of two flowcells)
    # Each flowcell consists of 8 lanes
    # A lane yields either 12 or 28.8 Gbp depending on flowcell version.

    days = []
    total_number_of_Gbps = []
    for flowcell in sorted_data:
        days.append(flowcell[1][0] + datetime.timedelta(days=len(flowcell[0])))
        flowcell_Gbps = 8 * 28.8 if flowcell[-1][-3] == "C" else 8 * 12.0
        if len(total_number_of_Gbps) == 0:
            total_number_of_Gbps.append(flowcell_Gbps)
        else:
            total_number_of_Gbps.append(flowcell_Gbps + total_number_of_Gbps[-1])

    if plot:
        plt.fill_between(days, total_number_of_Gbps, linewidth=2)
        plt.ylabel("Total number of Gigabasepairs")
        plt.show()

    return total_number_of_Gbps[-1]


def get_total_number_of_gigabytes(plot=False):
    """In this case 'gigabytes' refers to 10^9 bytes
    """
    data = get_formatted_data()
    data = sorted(data, key=lambda l: l[-1])
    grouped_data = [list(g) for k, g in itertools.groupby(data, itemgetter(-1))]
    sorted_data = trim_data(grouped_data)
    sorted_data = sorted(sorted_data, \
        key=lambda flowcell: flowcell[1][0] + \
        datetime.timedelta(days=len(flowcell[0])))

    days = []
    total_number_of_gbs = []
    for flowcell in sorted_data:
        days.append(flowcell[1][0] + datetime.timedelta(days=len(flowcell[0])))
        flowcell_gbs = float(flowcell[0][-1]) / 1000000000000.0
        if len(total_number_of_gbs) == 0:
            total_number_of_gbs.append(flowcell_gbs)
        else:
            total_number_of_gbs.append(flowcell_gbs + total_number_of_gbs[-1])

    if plot:
        textstr = "Total amount at %s: %.2f Tbp" % \
        (days[-1].strftime("%Y-%m-%d"), total_number_of_gbs[-1])
        fig = plt.figure()
        ax = fig.add_subplot("111")
        ax.fill_between(days, total_number_of_gbs, linewidth=2, alpha=0.6)
        ax.grid(True)
        props = dict(boxstyle='round, pad=1', facecolor='white', alpha=0.8)
        ax.text(0.965, 0.2, textstr, transform=ax.transAxes, fontsize=10, \
            horizontalalignment="right", bbox=props)
        plt.title("Cumulative amount of basepairs")
        plt.ylabel("Terabasepairs")
        # plt.show()
        fig.set_size_inches(10, 2.4)
        plt.subplots_adjust(top=0.80, bottom=0.15)
        fig.savefig("/Users/valentinesvensson/Desktop/gb_over_time.pdf", \
        ppi=500, transparent=True)

    return total_number_of_gbs[-1]
