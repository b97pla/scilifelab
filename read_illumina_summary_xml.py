import xml.etree.ElementTree as xml
import sys, os

def readSummaries(path):
    try:
        r1_tree = xml.parse(os.path.join(path, "read1.xml"))
        r2_tree = xml.parse(os.path.join(path, "read3.xml"))
    except:
        print ("Did not find expected files read1.xml and read3.xml")
        raise
    r1_root = r1_tree.getroot()
    r2_root = r2_tree.getroot()
    return([r1_root, r2_root])

def getLaneErrorRates(roots, lane):
    [r1_root, r2_root] = roots
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    for l in lanes1:
        if l.get("key")==str(lane):
            err1 = float(l.get("ErrRatePhiX"))
    for l in lanes2:
        if l.get("key")==str(lane):
            err2 = float(l.get("ErrRatePhiX"))
    #return (err1+err2)/2.0
    return [err1, err2]

def getErrorRates(summary):
    [r1_root, r2_root] = summary
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    err_rates_1 = {}
    err_rates_2 = {}
    for l in lanes1:
        err_rates_1[l.get("key")] = float(l.get("ErrRatePhiX"))
    for l in lanes2:
        err_rates_2[l.get("key")] = float(l.get("ErrRatePhiX"))
    #err_rates = {}
    #for k in err_rates_1.keys():
    #    err_rates[k] = ( err_rates_1[k] + err_rates_2[k] ) / 2.0
    #return err_rates
    return [err_rates_1, err_rates_2]

def getLaneErrorSD(roots, lane):
    [r1_root, r2_root] = roots
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    for l in lanes1:
        if l.get("key")==str(lane):
            err1 = float(l.get("ErrRatePhiXSD"))
    for l in lanes2:
        if l.get("key")==str(lane):
            err2 = float(l.get("ErrRatePhiXSD"))
    #return (err1+err2)/2.0
    return [err1, err2]

def getErrorSD(summary):
    [r1_root, r2_root] = summary
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    err_rates_1 = {}
    err_rates_2 = {}
    for l in lanes1:
        err_rates_1[l.get("key")] = float(l.get("ErrRatePhiXSD"))
    for l in lanes2:
        err_rates_2[l.get("key")] = float(l.get("ErrRatePhiXSD"))
    #err_rates = {}
    #for k in err_rates_1.keys():
    #    err_rates[k] = ( err_rates_1[k] + err_rates_2[k] ) / 2.0
    #return err_rates
    return [err_rates_1, err_rates_2]

def getLaneClustersRaw(roots, lane):
    [r1_root, r2_root] = roots
    density_1 = float(r1_root.get("densityRatio"))
    density_2 = float(r1_root.get("densityRatio"))
    assert ( density_1 == density_2 )
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    for l in lanes1:
        if l.get("key")==str(lane):
            clr1 = str ( int ( density_1 * float(l.get("ClustersRaw")) / 1000) ) + "K"
    for l in lanes2:
        if l.get("key")==str(lane):
            clr2 = str ( int ( density_2 * float(l.get("ClustersRaw")) / 1000) ) + "K"
    return [clr1, clr2]

def getClustersRaw(summary):
    [r1_root, r2_root] = summary
    density_1 = float(r1_root.get("densityRatio"))
    density_2 = float(r1_root.get("densityRatio"))
    assert ( density_1 == density_2 )
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    clr_1 = {}
    clr_2 = {}
    for l in lanes1:
        clr_1[l.get("key")] = str ( int ( density_1 * float(l.get("ClustersRaw")) / 1000) ) + "K"
    for l in lanes2:
        clr_2[l.get("key")] = str ( int ( density_2 * float(l.get("ClustersRaw")) / 1000) ) + "K"
    assert(clr_1 == clr_2)
    return clr_1

def getLaneClustersRawSD(roots, lane):
    [r1_root, r2_root] = roots
    density_1 = float(r1_root.get("densityRatio"))
    density_2 = float(r1_root.get("densityRatio"))
    assert ( density_1 == density_2 )
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    for l in lanes1:
        if l.get("key")==str(lane):
            clr1 = str ( int ( density_1 * float(l.get("ClustersRawSD")) / 1000) ) + "K"
    for l in lanes2:
        if l.get("key")==str(lane):
            clr2 = str ( int ( density_2 * float(l.get("ClustersRawSD")) / 1000) ) + "K"
    return [clr1, clr2]

def getClustersRawSD(summary):
    [r1_root, r2_root] = summary
    density_1 = float(r1_root.get("densityRatio"))
    density_2 = float(r1_root.get("densityRatio"))
    assert ( density_1 == density_2 )
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    clr_1 = {}
    clr_2 = {}
    for l in lanes1:
        clr_1[l.get("key")] = str ( int ( density_1 * float(l.get("ClustersRawSD")) / 1000) ) + "K"
    for l in lanes2:
        clr_2[l.get("key")] = str ( int ( density_2 * float(l.get("ClustersRawSD")) / 1000) ) + "K"
    assert(clr_1 == clr_2)
    return clr_1

def getLaneClustersPF(roots, lane):
    [r1_root, r2_root] = roots
    density_1 = float(r1_root.get("densityRatio"))
    density_2 = float(r1_root.get("densityRatio"))
    assert ( density_1 == density_2 )
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    for l in lanes1:
        if l.get("key")==str(lane):
            clr1 = str ( int ( density_1 * float(l.get("ClustersPF")) / 1000) ) + "K"
    for l in lanes2:
        if l.get("key")==str(lane):
            clr2 = str ( int ( density_2 * float(l.get("ClustersPF")) / 1000) ) + "K"
    return [clr1, clr2]

def getClustersPF(summary):
    [r1_root, r2_root] = summary
    density_1 = float(r1_root.get("densityRatio"))
    density_2 = float(r1_root.get("densityRatio"))
    assert ( density_1 == density_2 )
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    clr_1 = {}
    clr_2 = {}
    for l in lanes1:
        clr_1[l.get("key")] = str ( int ( density_1 * float(l.get("ClustersPF")) / 1000) ) + "K"
    for l in lanes2:
        clr_2[l.get("key")] = str ( int ( density_2 * float(l.get("ClustersPF")) / 1000) ) + "K"
    assert(clr_1 == clr_2)
    return clr_1

def getLaneClustersPFSD(roots, lane):
    [r1_root, r2_root] = roots
    density_1 = float(r1_root.get("densityRatio"))
    density_2 = float(r1_root.get("densityRatio"))
    assert ( density_1 == density_2 )
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    for l in lanes1:
        if l.get("key")==str(lane):
            clr1 = str ( int ( density_1 * float(l.get("ClustersPF")) / 1000) ) + "K"
    for l in lanes2:
        if l.get("key")==str(lane):
            clr2 = str ( int ( density_2 * float(l.get("ClustersPF")) / 1000) ) + "K"
    return [clr1, clr2]

def getClustersPFSD(summary):
    [r1_root, r2_root] = summary
    density_1 = float(r1_root.get("densityRatio"))
    density_2 = float(r1_root.get("densityRatio"))
    assert ( density_1 == density_2 )
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    clr_1 = {}
    clr_2 = {}
    for l in lanes1:
        clr_1[l.get("key")] = str ( int ( density_1 * float(l.get("ClustersPFSD")) / 1000) ) + "K"
    for l in lanes2:
        clr_2[l.get("key")] = str ( int ( density_2 * float(l.get("ClustersPFSD")) / 1000) ) + "K"
    assert(clr_1 == clr_2)
    return clr_1

def getLanePrcPF(roots, lane):
    [r1_root, r2_root] = roots
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    for l in lanes1:
        if l.get("key")==str(lane):
            p1 = float(l.get("PrcPFClusters"))
    for l in lanes2:
        if l.get("key")==str(lane):
            p2 = float(l.get("PrcPFClusters"))
    return [p1, p2]

def getPrcPF(summary):
    [r1_root, r2_root] = summary
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    p_1 = {}
    p_2 = {}
    for l in lanes1:
        p_1[l.get("key")] = float(l.get("PrcPFClusters"))
    for l in lanes2:
        p_2[l.get("key")] = float(l.get("PrcPFClusters"))
    return [p_1, p_2]

def getLanePrcPFSD(roots, lane):
    [r1_root, r2_root] = roots
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    for l in lanes1:
        if l.get("key")==str(lane):
            p1 = float(l.get("PrcPFClustersSD"))
    for l in lanes2:
        if l.get("key")==str(lane):
            p2 = float(l.get("PrcPFClustersSD"))
    return [p1, p2]

def getPrcPFSD(summary):
    [r1_root, r2_root] = summary
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    p_1 = {}
    p_2 = {}
    for l in lanes1:
        p_1[l.get("key")] = float(l.get("PrcPFClustersSD"))
    for l in lanes2:
        p_2[l.get("key")] = float(l.get("PrcPFClustersSD"))
    return [p_1, p_2]

def getLanePhasing(roots, lane):
    [r1_root, r2_root] = roots
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    for l in lanes1:
        if l.get("key")==str(lane):
            p1 = float(l.get("Phasing"))
    for l in lanes2:
        if l.get("key")==str(lane):
            p2 = float(l.get("Phasing"))
    return [p1, p2]

def getPhasing(summary):
    [r1_root, r2_root] = summary
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    p_1 = {}
    p_2 = {}
    for l in lanes1:
        p_1[l.get("key")] = float(l.get("Phasing"))
    for l in lanes2:
        p_2[l.get("key")] = float(l.get("Phasing"))
    return [p_1, p_2]


def getLanePrephasing(roots, lane):
    [r1_root, r2_root] = roots
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    for l in lanes1:
        if l.get("key")==str(lane):
            p1 = float(l.get("Prephasing"))
    for l in lanes2:
        if l.get("key")==str(lane):
            p2 = float(l.get("Prephasing"))
    return [p1, p2]

def getPrephasing(summary):
    [r1_root, r2_root] = summary
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    p_1 = {}
    p_2 = {}
    for l in lanes1:
        p_1[l.get("key")] = float(l.get("Prephasing"))
    for l in lanes2:
        p_2[l.get("key")] = float(l.get("Prephasing"))
    return [p_1, p_2]

def getLanePrcAlign(roots, lane):
    [r1_root, r2_root] = roots
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    for l in lanes1:
        if l.get("key")==str(lane):
            p1 = float(l.get("PrcAlign"))
    for l in lanes2:
        if l.get("key")==str(lane):
            p2 = float(l.get("PrcAlign"))
    return [p1, p2]

def getPrcAlign(summary):
    [r1_root, r2_root] = summary
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    p_1 = {}
    p_2 = {}
    for l in lanes1:
        p_1[l.get("key")] = float(l.get("PrcAlign"))
    for l in lanes2:
        p_2[l.get("key")] = float(l.get("PrcAlign"))
    return [p_1, p_2]

def getLanePrcAlignSD(roots, lane):
    [r1_root, r2_root] = roots
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    for l in lanes1:
        if l.get("key")==str(lane):
            p1 = float(l.get("PrcAlignSD"))
    for l in lanes2:
        if l.get("key")==str(lane):
            p2 = float(l.get("PrcAlignSD"))
    return [p1, p2]

def getPrcAlignSD(summary):
    [r1_root, r2_root] = summary
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    p_1 = {}
    p_2 = {}
    for l in lanes1:
        p_1[l.get("key")] = float(l.get("PrcAlignSD"))
    for l in lanes2:
        p_2[l.get("key")] = float(l.get("PrcAlignSD"))
    return [p_1, p_2]

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("You must supply a path to an XML summary directory")
    xmlpath = sys.argv[1]
    roots = readSummaries(xmlpath)
    print "Error rates:"
    err = getErrorRates(roots)
    print "Read 1"
    for e in sorted(err[0].keys()):
        print e, ":", err[0][e]
    print "Read 2"
    for e in sorted(err[1].keys()):
        print e, ":", err[1][e]
    print "Raw cluster densities:"
    clur = getClustersRaw(roots)
    for c in sorted(clur.keys()):
        print c, ":", clur[c]
    print "Raw cluster density SDs:"
    clursd = getClustersRawSD(roots)
    for c in sorted(clursd.keys()):
        print c, ":", clursd[c]
    print "PF cluster densities:"
    clupf = getClustersPF(roots)
    for c in sorted(clupf.keys()):
        print c, ":", clupf[c]
    print "PF cluster density SDs:"
    clupfsd = getClustersPFSD(roots)
    for c in sorted(clupfsd.keys()):
        print c, ":", clupfsd[c]
    print "Phasing:"
    phas = getPhasing(roots)
    print "Read 1:"
    for c in sorted(phas[0]):
        print c, ":", phas[0][c]
    print "Read 2:"
    for c in sorted(phas[1]):
        print c, ":", phas[1][c]
    print "Prephasing:"
    pphas = getPrephasing(roots)
    print "Read 1:"
    for c in sorted(pphas[0]):
        print c, ":", pphas[0][c]
    print "Read 2:"
    for c in sorted(pphas[1]):
        print c, ":", pphas[1][c]
    print "Prc aligned:"
    prc = getPrcAlign(roots)
    print "Read 1:"
    for c in sorted(prc[0]):
        print c, ":", prc[0][c]
    print "Read 2:"
    for c in sorted(prc[1]):
        print c, ":", prc[1][c]
