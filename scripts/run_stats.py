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
grdevices.png(file="runs.png", width=512, height=512)

runs = DataFrame.from_csvfile("tar_stats.csv")
#print runs

# Set axis labels perpendicular
#ro.r.par(las=1)
gp = ggplot2.ggplot(runs)

pp = gp + \
     ggplot2.aes_string(x='run', y='size') + \
     ggplot2.geom_histogram(colour = "darkgreen", fill = "blue", binwidth = 0.5)

pp.plot()
grdevices.dev_off()

#while(1): True
