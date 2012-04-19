
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

def main():
    # Plot chunks of 8
    chunksize = 8
    offset = 1
    while offset < len(sys.argv):
        upper = min(offset+chunksize,len(sys.argv))
        plot_fastq_sizes(sys.argv[offset:upper])
        offset = upper

def plot_fastq_sizes(fastq_files):
    
    files = []
    rx = []
    ry = []
    rh = []
    for infile in fastq_files:
        lens = {}
        hdata = []
        fh = FastQParser(infile)
        for record in fh:
            l = len(record[1])
            hdata.append(l)
            if l in lens:
                lens[l] += 1
            else:
                lens[l] = 1
        tot = fh.rread()
        fh.close()
        
        files.append(os.path.basename(infile))
        py = []
        px = []
        cum = 0
        for i in range(max(lens.keys())):
            px.append(i+1)
            if (i+1) in lens:
                cum += lens[i+1]
            py.append(1.*cum/tot)
        
        rx.append(robjects.IntVector(px))
        ry.append(robjects.FloatVector(py))
        rh.append(robjects.IntVector(hdata))
    
    if len(files) > 1:
        title = os.path.commonprefix(files)[0:-1]
        
        names = [f[len(title):].strip("_") for f in files]
        title = title.strip("_")
    
        suffix = os.path.commonprefix([n[::-1] for n in names])
        id = "-".join([n[0:min(5,len(n)-len(suffix))] for n in names])
        prefix = "%s_%s" % (title,id)
    else:
        title,_ = os.path.splitext(files[0])
        names = [title]
        prefix = title
        
    ltys = robjects.IntVector(range(1,5))
    rbow = robjects.r.rainbow(len(files))
    robjects.r.png("%s_cutadapt-size.png" % (prefix), width=1*1920, height=1*1080)
    
    histplots = len(files)
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
    
    #robjects.r.par(mfrow=robjects.r['c'](len(files)+1,1))
    
    # plot histograms of lengths
    for i,h in enumerate(rh):
        robjects.r.hist(h, breaks=rx[i], xlab="Length", main="%s" % (files[i]), col=rbow[i])
        
    # plot the cumulative fractions
    matx = robjects.r['c']()
    for rvec in rx:
        matx = robjects.r['c'](matx,rvec)
    matx = robjects.r['matrix'](matx,ncol=len(files))
    maty = robjects.r['c']()
    for rvec in ry:
        maty = robjects.r['c'](maty,rvec)
    maty = robjects.r['matrix'](maty,ncol=len(files))
    
    robjects.r.matplot(matx,maty,type="l", col=rbow, lty=ltys, xlab="Length",ylab="Cumulative read fraction",main=title,xlim=robjects.r['c'](robjects.r.min(matx),robjects.r.max(matx)),ylim=robjects.IntVector([0,1]))
    #robjects.r.plot(rx[0], ry[0], type="l", col=rbow[0], lty=ltys[0], xlab="Length",ylab="Cumulative read fraction",main=title,xlim=robjects.IntVector([min(rx[0]),max(rx[0])]),ylim=robjects.IntVector([0,1]))
    #for i,x in enumerate(rx[1:]):
    #    robjects.r.lines(rx[i+1],ry[i+1],col=rbow[i+1],lty=ltys[i+1])
    robjects.r.legend(x="bottomright",legend=files,col=rbow,lty=ltys,cex=2)


if __name__ == "__main__":
    main()
