"""report qc module"""
import os
import yaml

from collections import OrderedDict
from cStringIO import StringIO

from scilifelab.db.statusdb import SampleRunMetricsConnection, ProjectSummaryConnection, FlowcellRunMetricsConnection
from scilifelab.bcbio.qc import SampleRunMetricsParser
from scilifelab.log import minimal_logger
from scilifelab.bcbio.run import find_samples

LOG = minimal_logger(__name__)

## QC data cutoff values
QC_CUTOFF = {
    'rnaseq':{'PCT_PF_READS_ALIGNED':70,'PERCENT_DUPLICATION':30},
    'reseq':{'PCT_PF_READS_ALIGNED':70,'PERCENT_DUPLICATION':30},
    'WG-reseq':{'PCT_PF_READS_ALIGNED':70,'PERCENT_DUPLICATION':30},
    'seqcap':{'PCT_PF_READS_ALIGNED':70,'PERCENT_ON_TARGET':60, 'PCT_TARGET_BASES_10X':90, 'PERCENT_DUPLICATION':30},
    'customcap':{'PCT_PF_READS_ALIGNED':70, 'PERCENT_DUPLICATION':30},
    'finished':{},
    }

## QC data header information
HEADER = ["sample","lane","flowcell", "date",  "TOTAL_READS",
          "MEAN_INSERT_SIZE", "GENOME_SIZE", "PERCENT_ON_TARGET",
          "PERCENT_DUPLICATION", "PCT_TARGET_BASES_10X", "PCT_PF_READS_ALIGNED",
          "dup_status", "status"]
HEADER_LABELS = ["sample","lane","flowcell", "date",  "#READS",
                 "INS_SIZE", "G_SIZE", "PCT_TGT",
                 "PCT_DUP", "PCT_10X", "PCT_ALN",
                 "dup_status", "status"]
HEADER_MAP = OrderedDict(zip(HEADER, HEADER_LABELS))

## FIXME: this should be used also in assess_qc to set column widths
COL_WIDTHS = [20, 5, 11, 8, 11, 10, 11, 10, 10, 10, 10, 12, 12]

## Mapping from genomics project list application names
APPLICATION_MAP = {'RNA-seq (Total RNA)':'rnaseq','WG re-seq':'WG-reseq','Resequencing':'reseq', 'Exome capture':'seqcap', 'Custom':'customcap', 'Finished library':'finished' , 'Custom capture':'customcap'}
APPLICATION_INV_MAP = {v:k for k, v in APPLICATION_MAP.items()}

def _srm_to_qc(srm, application=None):
    """Given a sample run metrics object, return relevant QC info"""
    qcdata = {"sample":srm.get("barcode_name", None),
              "project":srm.get("sample_prj", None),
              "lane":srm.get("lane", None),
              "flowcell":srm.get("flowcell", None),
              "date":srm.get("date", None),
              "application":application,
              "TOTAL_READS":int(srm.get("picard_metrics", {}).get("AL_PAIR", {}).get("TOTAL_READS", -1)),
              "PERCENT_DUPLICATION":srm.get("picard_metrics", {}).get("DUP_metrics", {}).get("PERCENT_DUPLICATION", "-1.0"),
              "MEAN_INSERT_SIZE":float(srm.get("picard_metrics", {}).get("INS_metrics", {}).get("MEAN_INSERT_SIZE", "-1.0").replace(",", ".")),
              "GENOME_SIZE":int(srm.get("picard_metrics", {}).get("HS_metrics", {}).get("GENOME_SIZE", -1)),
              "FOLD_ENRICHMENT":float(srm.get("picard_metrics", {}).get("HS_metrics", {}).get("FOLD_ENRICHMENT", "-1.0").replace(",", ".")),
              "PCT_USABLE_BASES_ON_TARGET":srm.get("picard_metrics", {}).get("HS_metrics", {}).get("PCT_USABLE_BASES_ON_TARGET", "-1.0"),
              "PCT_TARGET_BASES_10X":srm.get("picard_metrics", {}).get("HS_metrics", {}).get("PCT_TARGET_BASES_10X", "-1.0"),
              "PCT_PF_READS_ALIGNED":srm.get("picard_metrics", {}).get("AL_PAIR", {}).get("PCT_PF_READS_ALIGNED", "-1.0"),
              }

    target_territory = float(srm.get("picard_metrics", {}).get("HS_metrics", {}).get("TARGET_TERRITORY", -1))
    pct_labels = ["PERCENT_DUPLICATION", "PCT_USABLE_BASES_ON_TARGET", "PCT_TARGET_BASES_10X",
                  "PCT_PF_READS_ALIGNED"]
    for l in pct_labels:
        if qcdata[l]:
            qcdata[l] = float(qcdata[l].replace(",", ".")) * 100
    if qcdata["FOLD_ENRICHMENT"] and qcdata["GENOME_SIZE"] and target_territory:
        qcdata["PERCENT_ON_TARGET"] = float(qcdata["FOLD_ENRICHMENT"]/ (float(qcdata["GENOME_SIZE"]) / float(target_territory))) * 100
    return qcdata

def _get_sample_qc_data(sample_prj, application, s_con, fc_id=None):
    samples = s_con.get_samples(fc_id=fc_id, sample_prj=sample_prj)
    qcdata = {}
    for s in samples:
        qcdata[s["name"]]=_srm_to_qc(s)
    return qcdata

def assess_qc(x, application):
    status = "PASS"
    dup_status = "OK"
    for k in QC_CUTOFF[application].keys():
        LOG.debug("assessing qc metric {}".format(k))
        if k == "PERCENT_DUPLICATION":
            if float(x[k]) > QC_CUTOFF[application][k]: 
                dup_status = "HIGH"
        else:
            if float(x[k]) < QC_CUTOFF[application][k]:
                status = "FAIL"
    return ["{:20}".format(x["sample"]),
            "{:>5}".format(x["lane"]),
            "{:>11}".format(x["flowcell"]),
            "{:>8}".format(x["date"]), 
            "{:>10.2f}M".format(int(x["TOTAL_READS"])/1e6/2), 
            "{:>10.1f}".format(float(x["MEAN_INSERT_SIZE"])),
            "{:>10.1f}G".format(int(x["GENOME_SIZE"])/1e9) if x["GENOME_SIZE"]>1e9 else "{:.1f}M".format(int(x["GENOME_SIZE"])/1e6),
            "{:>10.1f}".format(float(x["PERCENT_ON_TARGET"])),
            "{:>10.1f}".format(float(x["PERCENT_DUPLICATION"])),
            "{:>10.1f}".format(float(x["PCT_TARGET_BASES_10X"])), 
            "{:>10.1f}".format(float(x["PCT_PF_READS_ALIGNED"])), 
            "{:>12}".format(dup_status), 
            "{:>12}".format(status)]

def compile_qc(path, application="seqcap", **kw):
    """Perform qc on data without access to statusdb.

    :param **kw: keyword argument

    """
    output_data = {'stdout':StringIO(), 'stderr':StringIO()}
    ### find_samples excrutiatingly slow for multi-sample projects where we can have > 100k files...
    flist = find_samples(path, **kw)
    srm_l = []
    for f in flist:
        LOG.debug("Opening config file {}".format(f))
        with open(f) as fh:
            runinfo_yaml = yaml.load(fh)
        for info in runinfo_yaml['details']:
            if info.get("multiplex", None):
                for mp in info.get("multiplex"):
                    sample_kw = dict(path=os.path.dirname(f), flowcell=runinfo_yaml.get("fc_name", None), date=runinfo_yaml.get("fc_date", None), lane=info.get("lane", None), barcode_name=mp.get("name", None), sample_prj=kw.get("project"), barcode_id=mp.get('barcode_id', None), sequence=mp.get('sequence', None))
                    obj = SampleRunMetrics(**sample_kw)
                    srm_l.append(obj)
            else:
                sample_kw = dict(path=os.path.dirname(f), flowcell=runinfo_yaml.get("fc_name", None), date=runinfo_yaml.get("fc_date", None), lane=info.get("lane", None), barcode_name=info.get("description", None), sample_prj=kw.get("project"), barcode_id=None, sequence=None)
                obj = SampleRunMetrics(**sample_kw)
                obj.read_picard_metrics()
                srm_l.append(obj)
    qcdata = []
    output_data = _qc_info_header(kw.get("project"), application, output_data)
    for s in srm_l:
        qcdata.append(_srm_to_qc(s))
    for v in qcdata:
        y = [str(x) for x in assess_qc(v, application)]
        output_data["stdout"].write("".join(y) + "\n")
    return output_data


def _qc_info_header(project, application, output_data):
    ## Add info about qc
    output_data["stdout"].write("\n\nApplication QC\n")
    output_data["stdout"].write("==============\n")
    output_data["stdout"].write("{:12}{:>16}\n".format("Project", project.rstrip()))
    output_data["stdout"].write("{:12}{:>16}\n".format("Application", APPLICATION_INV_MAP[application]))
    output_data["stdout"].write("==============\n\n")        
    output_data["stdout"].write("Application QC criteria\n")
    output_data["stdout"].write("==============\n")        
    for k,v in sorted(QC_CUTOFF[application].iteritems()):
        output_data["stdout"].write("{:<24}={:>4}\n".format(k, v))
    output_data["stdout"].write("==============\n\n") 
    output_data["stdout"].write("Header legend:\n")
    output_data["stdout"].write(", ".join(["{} = {}".format(v,k) for k, v in HEADER_MAP.iteritems()][5:11]))
    output_data["stdout"].write("\n\n")
    header_fmt = ["{" + ":>{}".format(COL_WIDTHS[i]) + "}" for i in range(0, len(COL_WIDTHS))]
    header_fmt[0] = "{" + ":<{}".format(COL_WIDTHS[0]) + "}"
    header_out = "".join([header_fmt[i].format(HEADER_MAP[HEADER_MAP.keys()[i]]) for i in range(0, len(COL_WIDTHS))])
    output_data["stdout"].write(header_out + "\n")
    return output_data

def application_qc(project=None, flowcell=None, application=None,
                   user=None, password=None, url=None,
                   use_ps_map=True, use_bc_map=False, check_consistency=False, **kw):
    """Perform application specific qc on a project.

    :param project: project identifier
    :param flowcell: flowcell identifier
    :param application: application for which to perform qc
    :param user: database username
    :param password: database password
    :param url: database url
    :param use_ps_map: use project summary mapping
    :param use_bc_map: use barcode name mapping
    :param check_consistency: check 
    """
    LOG.debug("Doing application qc for project {}, flowcell {}".format(project, flowcell))

    output_data = {'stdout':StringIO(), 'stderr':StringIO()}
    p_con = ProjectSummaryConnection(username=user, password=password, url=url)
    project = p_con.get_entry(project)
    qc_data = p_con.get_qc_data(project, flowcell)

    if not project is None:
        qc_data = p_con.get_qc_data(project, flowcell)
        if project.get("application") not in application_map.keys():
            if not application:
                LOG.warn("No such application {}. Please use the application option (available choices {})".format(application, ",".join(QC_CUTOFF.keys())))
                return output_data
            application = application
        else:
            application = application_map[project.get("application")]
    else:
        LOG.info("No such project {} in project summary. Trying to get qc data anyway.".format(project))
        if not application:
            LOG.warn("No application provided. Please use the application option (available choices {})".format(",".join(QC_CUTOFF.keys())))
            return output_data
        s_con = SampleRunMetricsConnection(username=user, password=password, url=url)
        qc_data = _get_sample_qc_data(project, application, s_con, flowcell)

    output_data = _qc_info_header(project, application, output_data)
    for k,v in sorted(qc_data.iteritems()):
        y = [str(x) for x in assess_qc(v, application)]
        output_data["stdout"].write("".join(y) + "\n")
    return output_data

def fastq_screen(project=None, flowcell=None,
                 user=None, password=None, url=None,
                 use_ps_map=True, use_bc_map=False, check_consistency=False, **kw):
    output_data = {'stdout':StringIO(), 'stderr':StringIO()}
    s_con = SampleRunMetricsConnection(username=user, password=password, url=url)
    samples = s_con.get_samples(fc_id=flowcell, sample_prj=project)
    for s in samples:
        LOG.debug("Checking fastq_screen data for sample {}, id {}, project {}".format(s.get("name", None), s.get("_id", None), s.get("sample_prj", None)))
        fqscreen_data = s.get("fastq_scr", {})
        output_data["stdout"].write(s["barcode_name"] + "\n")
        if fqscreen_data:
            header = [[x for x in v.keys()] for k, v in fqscreen_data.iteritems()]
            output_data["stdout"].write("\t\t" + "".join("{:>27}".format(x) for x in header[0]) + "\n")
            vals = ["{:>12}\t{}\n".format(k, "".join(["{:>27}".format(x) for x in v.values()])) for k, v in fqscreen_data.iteritems()]
            for v in vals:
                output_data["stdout"].write(v)
    return output_data
