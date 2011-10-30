import xml.etree.ElementTree as xml
import sys, os

# getQCstats() is probably the method you usually want to call
def getQCstats(path):
    r = readSummaries(path)
    qc_stats = {}
    qc_stats['error_rate'] = getAllLaneMetrics(r, 'ErrRatePhiX', False)
    qc_stats['error_rate_sd'] = getAllLaneMetrics(r, 'ErrRatePhiXSD', False)
    qc_stats['raw_cluster_dens'] = getAllLaneMetrics(r, 'ClustersRaw', True)
    qc_stats['raw_cluster_dens_sd'] = getAllLaneMetrics(r, 'ClustersRawSD', True)
    qc_stats['prc_pf'] = getAllLaneMetrics(r, 'PrcPFClusters', False)
    qc_stats['prc_pf_sd'] = getAllLaneMetrics(r, 'PrcPFClustersSD', False)
    qc_stats['pf_cluster_dens'] = getAllLaneMetrics(r, 'ClustersPF', True)
    qc_stats['pf_cluster_dens_sd'] = getAllLaneMetrics(r, 'ClustersPFSD', True)
    qc_stats['phasing'] = getAllLaneMetrics(r, 'Phasing', False)
    qc_stats['prephasing'] = getAllLaneMetrics(r, 'Prephasing', False)
    qc_stats['prc_aligned'] = getAllLaneMetrics(r, 'PrcAlign', False)
    qc_stats['prc_aligned_sd'] = getAllLaneMetrics(r, 'PrcAlign', False)
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

def getSingleLaneMetrics(roots, metric, lane, clu_dens = False):
    [r1_root, r2_root] = roots
    if clu_dens: densRatio = float(r1_root.get("densityRatio"))
    lanes = r1_root.findall("Lane")    
    for l in lanes1:
        if l.get("key")==str(lane):
            m1 = float(l.get(metric))
            if clu_dens: m1 = str(int(round((densRatio * float(l.get(metric)))/1000))) + 'K'
    for l in lanes2:
        if l.get("key")==str(lane):
            m2 = float(l.get(metric))
            if clu_dens: m2 = str(int(round((densRatio * float(l.get(metric)))/1000))) + 'K'
    return {'read1':m1, 'read2':m2}

def getAllLaneMetrics(roots, metric, clu_dens):
    [r1_root, r2_root] = roots
    if clu_dens: densRatio = float(r1_root.get("densityRatio"))
    lanes1 = r1_root.findall("Lane")    
    lanes2 = r2_root.findall("Lane")
    m1 = {}
    m2 = {}
    for l in lanes1:
        m1[l.get("key")] = float(l.get(metric))
        if clu_dens: m1[l.get("key")] = str(int(round((densRatio * float(l.get(metric)))/1000))) + 'K'
    for l in lanes2:
        m2[l.get("key")] = float(l.get(metric))
        if clu_dens: m2[l.get("key")] = str(int(round((densRatio * float(l.get(metric)))/1000))) + 'K'
    return {'read1':m1, 'read2':m2}

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("You must supply a path to an XML summary directory")
    xmlpath = sys.argv[1]
    q = getQCstats(xmlpath)
    print q
