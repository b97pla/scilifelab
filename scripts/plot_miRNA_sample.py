
import rpy2.robjects as R
from miseq_data import FastQParser
import sys 
import os 
import math
import collections
from optparse import OptionParser


XRES = 1920
YRES = 1080

def main(cutadapt_metrics_file):
    cutadapt_metrics = parse_cutadapt_metrics(cutadapt_metrics_file)
    seqfile = cutadapt_metrics['trimmed_output']
    
    fileprefix = "".join(os.path.basename(seqfile).split(".")[0:-2])
    lengthfile = "%s-length.length_metrics" % fileprefix
    length_metrics = parse_length_metrics(lengthfile,seqfile)
    x_data = [l+1 for l in range(max(length_metrics['counts'].keys()))]
    
    plotfile = "%s-plot.png" % fileprefix
    setup_plot(plotfile)
    plot_adapter_removal_result(cutadapt_metrics)
    plot_length_histogram(length_metrics['counts'], x_data, [cutadapt_metrics.get('min_cutoff',0),cutadapt_metrics.get('max_cutoff',0)])
    plot_cumulative_length(length_metrics['counts'], x_data)
    
def parse_length_metrics(metrics_file, seqfile=None):
    """Parse the length metrics if available, or create and parse if a sequence file was passed"""

    counts = collections.defaultdict(int)
    if not os.path.exists(metrics_file):
        assert os.path.exists(seqfile), "The length metrics file is missing and could not find the sequence file %s" % seqfile
        fp = FastQParser(seqfile)
        for record in fp:
            counts[len(record[1])] += 1
        fp.close()
        
        with open(metrics_file,"w") as fh:
            for l in sorted(counts.keys()):
                fh.write("%d\t%d\n" % (l,counts[l]))

    else:
        with open(metrics_file) as fh:
            for row in fh:
                l,c = row.strip().split()
                counts[int(l)] = int(c)
    
    return {'counts': counts}
                

def parse_cutadapt_metrics(metrics_file):
    
    metrics = {}
    with open(metrics_file) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("cutadapt version"):
                metrics['program'] = line
            elif line.startswith("Command line parameters:"):
                parts = line.split()
                for i,p in enumerate(parts):
                    if p == "-a":
                        metrics['adapter'] = parts[i+1]
                    elif p == "-m":
                        metrics['min_cutoff'] = parts[i+1]
                    elif p == "-M":
                        metrics['max_cutoff'] = parts[i+1]
                    elif p == "-O":
                        metrics['overlap'] = parts[i+1]
                    elif p == "-o":
                        metrics['trimmed_output'] = parts[i+1]
            elif line.startswith("Processed reads:"):
                metrics['processed'] = int(line.split()[2])
            elif line.startswith("Trimmed reads:"):
                metrics['trimmed'] = int(line.split()[2])
            elif line.startswith("Too short reads:"):
                metrics['too_short'] = int(line.split()[3])
            elif line.startswith("Too long reads:"):
                metrics['too_long'] = int(line.split()[3])
    return metrics

def setup_plot(plotfile):

    R.r.png(plotfile, width=XRES, height=YRES)
    
    histplots = 4
    histrows = int(math.ceil(1.*histplots/2))
    plots = range(1,histplots+1) 
    if histplots%2 > 0:
        plots.append(0)
    plots.extend([histplots+1,histplots+1]) 
    lay = R.IntVector(plots) 
    lay = R.r['matrix'](lay,histrows+1,2,byrow="T")
    hghts = [1 for i in range(histrows)]
    hghts.extend([2])
    layout = R.r.layout(lay,heights=R.IntVector(hghts))
    

def plot_adapter_removal_result(metrics):
    
    titletext = "Adapter trimming"
    subtitletext = metrics['program']
    heights_r = R.FloatVector([100,100*float(metrics['trimmed'])/metrics['processed'],100*float(metrics['too_short'])/metrics['processed'],100*float(metrics['too_long'])/metrics['processed']])
    names_r = R.StrVector(['Total reads','Reads w/ adapter','Too short reads','Too long reads'])
    R.r.barplot(heights_r,main=titletext,sub=subtitletext,**{'names.arg': names_r})

def plot_length_histogram(length_counts, x_data, cutoffs):
    
    titletext = "Read length distribution after adapter trimming"
    #R.r.barplot(R.IntVector([length_counts.get(l,0) for l in x_data]), main=titletext, **{'names.arg': R.IntVector(x_data)})
    h = []
    for l,c in length_counts.items():
        h.extend([l for _ in range(c)])
    R.r.hist(R.IntVector(h), breaks=R.IntVector(x_data), main=titletext, xlab="Read length")
    for c in cutoffs:
        if c > 0:
            R.r.abline(v=c,lty=2)
    
def plot_cumulative_length(length_counts, x_data):
    
    titletext = "Cumulative fraction of reads vs read length"
    
    total = sum(length_counts.values())
    cumlength = [length_counts.get(x_data[0],0)]
    for i in range(1,len(x_data)):
        cumlength.append(cumlength[i-1] + length_counts.get(x_data[i],0))
    
    cumfraction = [float(l)/total for l in cumlength]    
    R.r.plot(R.IntVector(x_data),R.FloatVector(cumfraction),type="l",main=titletext,xlab="Length",ylab="Cumulative fraction of reads",xlim=R.IntVector([x_data[0],x_data[-1]]),ylim=R.IntVector([0,1]))
    
def plot_ncrna_distribution(metrics):
    
    titletext = "ncRNA classification"
    

if __name__ == "__main__":
    parser = OptionParser()
    options, args = parser.parse_args()
    main(args[0])
