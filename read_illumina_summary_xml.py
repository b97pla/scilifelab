import xml.etree.ElementTree as xml
import sys, os

# getQCstats() is probably the method you usually want to call
def getQCstats(path):
    r = readSummaries(path)
    qc_stats = {}
    qc_stats['error_rate'] = getErrorRates(r)
    qc_stats['error_rate_sd'] = getErrorSD(r)
    qc_stats['raw_cluster_dens'] = getClustersRaw(r)
    qc_stats['raw_cluster_dens_sd'] = getClustersRawSD(r)
    qc_stats['prc_pf'] = getPrcPF(r)
    qc_stats['prc_pf_sd'] = getPrcPFSD(r)
    qc_stats['pf_cluster_dens'] = getClustersPF(r)
    qc_stats['pf_cluster_dens_sd'] = getClustersPFSD(r)
    qc_stats['phasing'] = getPhasing(r)
    qc_stats['prephasing'] = getPrephasing(r)
    qc_stats['prc_aligned'] = getPrcAlign(r)
    qc_stats['prc_aligned_sd'] = getPrcAlignSD(r)
    return qc_stats

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
    return {'read1':err_rates_1, 'read2':err_rates_2}

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
    return {'read1':err1, 'read2':err2}

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
    return {'read1':err_rates_1, 'read2':err_rates_2}

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
    return {'read1':clr1, 'read2':clr2}

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
    return {'read1':clr_1, 'read2':clr_2}

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
    return {'read1':clr1, 'read2':clr2}

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
    return {'read1':clr_1, 'read2':clr_2}

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
    return {'read1':clr1, 'read2':clr2}

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
    return {'read1':clr_1, 'read2':clr_2}

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
    return {'read1':clr1, 'read2':clr2}

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
    return {'read1':clr_1, 'read2':clr_2}

def getLanePrcPF(roots, lane):
    [r1_root, r2_root] = roots
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    for l in lanes1:
        if l.get("key")==str(lane):
            p_1 = float(l.get("PrcPFClusters"))
    for l in lanes2:
        if l.get("key")==str(lane):
            p_2 = float(l.get("PrcPFClusters"))
    return {'read1':p_1, 'read2':p_2}

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
    return {'read1':p_1, 'read2':p_2}

def getLanePrcPFSD(roots, lane):
    [r1_root, r2_root] = roots
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    for l in lanes1:
        if l.get("key")==str(lane):
            p_1 = float(l.get("PrcPFClustersSD"))
    for l in lanes2:
        if l.get("key")==str(lane):
            p_2 = float(l.get("PrcPFClustersSD"))
    return {'read1':p_1, 'read2':p_2}

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
    return {'read1':p_1, 'read2':p_2}

def getLanePhasing(roots, lane):
    [r1_root, r2_root] = roots
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    for l in lanes1:
        if l.get("key")==str(lane):
            p_1 = float(l.get("Phasing"))
    for l in lanes2:
        if l.get("key")==str(lane):
            p_2 = float(l.get("Phasing"))
    return {'read1':p_1, 'read2':p_2}

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
    return {'read1':p_1, 'read2':p_2}

def getLanePrephasing(roots, lane):
    [r1_root, r2_root] = roots
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    for l in lanes1:
        if l.get("key")==str(lane):
            p_1 = float(l.get("Prephasing"))
    for l in lanes2:
        if l.get("key")==str(lane):
            p_2 = float(l.get("Prephasing"))
    return {'read1':p_1, 'read2':p_2}

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
    return {'read1':p_1, 'read2':p_2}

def getLanePrcAlign(roots, lane):
    [r1_root, r2_root] = roots
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    for l in lanes1:
        if l.get("key")==str(lane):
            p_1 = float(l.get("PrcAlign"))
    for l in lanes2:
        if l.get("key")==str(lane):
            p_2 = float(l.get("PrcAlign"))
    return {'read1':p_1, 'read2':p_2}

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
    return {'read1':p_1, 'read2':p_2}

def getLanePrcAlignSD(roots, lane):
    [r1_root, r2_root] = roots
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    for l in lanes1:
        if l.get("key")==str(lane):
            p_1 = float(l.get("PrcAlignSD"))
    for l in lanes2:
        if l.get("key")==str(lane):
            p_2 = float(l.get("PrcAlignSD"))
    return {'read1':p_1, 'read2':p_2}

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
    return {'read1':p_1, 'read2':p_2}

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("You must supply a path to an XML summary directory")
    xmlpath = sys.argv[1]
    q = getQCstats(xmlpath)
    print q
