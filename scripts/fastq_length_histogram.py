
#import fileinput
#
#lengths = {}
#seqs = 0
#prev_line = "-"
#for line in fileinput.input():
#    if prev_line[0] == '@':
#        l = len(line.strip())
#        if l in lengths:
#            lengths[l] += 1
#        else:
#            lengths[l] = 1
#        seqs += 1
#    prev_line = line
#    
#max_wd = 50
#max_count = max(lengths.values())
#scale = 1.*max_wd/max_count
#prev_c = 1
#
#py = []
#px = []
#cum = 0
#for i in range(max(lengths.keys())):
#    c = 0
#    py.append(1.*cum/seqs)
#    px.append(i+1)
#    if (i+1) in lengths:
#        c = int(round(scale*lengths[i+1]))
#        cum += lengths[i+1]
#        py[-1] = 1.*cum/seqs
#    if c > 0:
#        print "%d\t%s" % (i+1,"+" * c)
#    elif c == 0 and prev_c > 0:
#        print "-"
#    prev_c = c
#
#print "+ = %d reads" % int(round(1./scale))

import rpy2.robjects as robjects
from miseq_data import FastQParser
import sys 
import os 
import math
import collections

def main():
    # Plot chunks of 8
    chunksize = 8
    offset = 1
    while offset < len(sys.argv):
        upper = min(offset+chunksize,len(sys.argv))
        plot_fastq_sizes(sys.argv[offset:upper])
        offset = upper

def plot_fastq_sizes(fastq_files):
    
    length_counts = []
    maxlen = 0
    for infile in fastq_files:
        lens = collections.defaultdict(int)

        fh = FastQParser(infile)
        for record in fh:
            lens[len(record[1])] += 1
        fh.close()
        
        length_counts.append(lens)
        maxlen = max(maxlen,max(lens.keys()))    
        
    parts = os.path.basename(fastq_files[0]).split("_")
    # Lane_Date_FC
    prefix = "_".join(parts[0:3])
    # Extension and appended operation
    suffix = parts[-1]
    
    names = [os.path.basename(f)[len(prefix)+1:][:-1*(len(suffix)+1)] for f in fastq_files]
        
    ltys = robjects.IntVector(range(1,5))
    rbow = robjects.r.rainbow(len(fastq_files))
    robjects.r.png("%s_%s_cutadapt-size.png" % (prefix,"-".join(names)), width=1*1920, height=1*1080)
    
    histplots = len(fastq_files)
    histrows = int(math.ceil(1.*histplots/2))
    plots = range(1,histplots+1) 
    if histplots%2 > 0:
        plots.append(0)
    plots.extend([histplots+1,histplots+1]) 
    lay = robjects.IntVector(plots)
    lay = robjects.r['matrix'](lay,histrows+1,2,byrow="T")
    hghts = [1 for i in range(histrows)]
    hghts.append(2)
    layout = robjects.r.layout(lay,heights=robjects.IntVector(hghts))
    
    # plot histograms of lengths
    for i,lens in enumerate(length_counts):
        # create the histogram values based on the counts
        hist_values = []
        for l,c in lens.items():
            hist_values.extend([l for _ in range(c)])
            
        robjects.r.hist(robjects.IntVector(hist_values), breaks=robjects.IntVector(range(1,maxlen+1)), xlab="Length", main="%s" % (names[i]), col=rbow[i])
        
    # plot the cumulative fractions
    
    # Expand the data
    matx = robjects.r['c']()
    maty = robjects.r['c']()
    for i,lens in enumerate(length_counts):
        
        reads = float(sum(lens.values())) 
        c = [lens.get(1,0)/reads]
        for l in range(1,maxlen):
            c.append(c[l-1] + lens.get(l+1,0)/reads)
        
        matx = robjects.r['c'](matx,robjects.IntVector(range(1,maxlen+1)))
        maty = robjects.r['c'](maty,robjects.FloatVector(c))
        
    matx = robjects.r['matrix'](matx,ncol=len(fastq_files))
    maty = robjects.r['matrix'](maty,ncol=len(fastq_files))
    
    title = "Cumulative fraction of reads vs read length after trimming"
    robjects.r.matplot(matx,maty,type="l", col=rbow, lty=ltys, xlab="Length",ylab="Cumulative read fraction",main=title,xlim=robjects.r['c'](robjects.r.min(matx),robjects.r.max(matx)),ylim=robjects.IntVector([0,1]))
    #robjects.r.plot(rx[0], ry[0], type="l", col=rbow[0], lty=ltys[0], xlab="Length",ylab="Cumulative read fraction",main=title,xlim=robjects.IntVector([min(rx[0]),max(rx[0])]),ylim=robjects.IntVector([0,1]))
    #for i,x in enumerate(rx[1:]):
    #    robjects.r.lines(rx[i+1],ry[i+1],col=rbow[i+1],lty=ltys[i+1])
    robjects.r.legend(x="bottomright",legend=names,col=rbow,lty=ltys,cex=2)


if __name__ == "__main__":
    main()
