#!/usr/bin/env python
# Plots a histogram with sizes of runs over time

import math, datetime
from rpy2.robjects.vectors import DataFrame
import rpy2.robjects.lib.ggplot2 as ggplot2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
base = importr('base')
stats = importr('stats')
grdevices = importr('grDevices')

grdevices.png(file="runs.png", width=1024, height=768)

runs = DataFrame.from_csvfile("tar_stats.csv")

#gp = ggplot2.ggplot(runs)

#pp = gp + \
#     ggplot2.aes_string(x='runs', y='sizes') + \
#     ggplot2.geom_histogram(colour = "darkgreen", fill = "blue", binwidth = 0.5)

qplot = ggplot2.qplot

pp = qplot(factor(cyl), data=runs, geom="bar", fill=factor(vs))

pp.plot()
grdevices.dev_off()

#while(1): True
