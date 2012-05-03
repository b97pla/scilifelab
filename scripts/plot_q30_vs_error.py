
import rpy2.robjects as R
import xml.etree.ElementTree

import sys 
import os 
import glob

XRES = 1920
YRES = 1080

def main(flowcells):
    
    plotfile = "pctQ30_vs_error_rate.png"
    setup_plot(plotfile,subplots=len(flowcells))
    
    for fc in flowcells:
        fc_data = []
        reportdir = os.path.join('Data','reports')
        for error_xml in glob.glob(os.path.join(fc,reportdir,'ErrorRate','Chart_*.xml')):
            qual_xml = os.path.join(fc,reportdir,'NumGT30',os.path.basename(error_xml))
            if not os.path.exists(qual_xml):
                continue
            fc_data.extend(dicts_to_list(parse_xml_data(error_xml),parse_xml_data(qual_xml)))
        plot_xy(fc_data)
    
def setup_plot(plotfile,width=XRES,height=YRES,subplots=1):
    
    R.r.png(plotfile, width=XRES, height=YRES)
    
    plotrows = int(math.ceil(1.*subplots/2))
    plots = range(1,subplots+1) 
    if subplots%2 > 0:
        plots.append(0)
    lay = R.IntVector(plots) 
    lay = R.r['matrix'](lay,plotrows,2,byrow="T")
    hghts = [1 for i in range(plotrows)]
    layout = R.r.layout(lay,heights=R.IntVector(hghts))
    
    
def parse_xml_data(xmlfile):
    et = xml.etree.ElementTree.parse(xmlfile)
    elements = et.findall('FlowCellData/TL')
    data = {}
    for e in elements:
        data[e.get('Key')] = e.get('Val')
    return data

def dicts_to_list(errs,quals):
    l = []
    for key in [k for k in errs.keys() if k in quals]:
        l.append([float(errs[k]),float(quals[k])])
    return l

def plot_xy(data):
    
    titletext = "% bases >= Q30 vs PhiX error %"
    x_data = R.FloatVector([d[0] for d in data])
    y_data = R.FloatVector([d[1] for d in data])
    
    R.r.plot(x_data,y_data,type="p",main=titletext,xlab="PhiX error %",ylab="% bases >= Q30",xlim=R.FloatVector([0,R.max(x_data)]),ylim=R.IntVector([0,100]))

if __name__ == "__main__":
    main(sys.argv[1:])
