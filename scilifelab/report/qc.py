"""report qc module"""

from collections import OrderedDict
from cStringIO import StringIO

from scilifelab.db.statusdb import SampleRunMetricsConnection, ProjectSummaryConnection, FlowcellRunMetricsConnection
from scilifelab.log import minimal_logger

LOG = minimal_logger(__name__)

## QC data cutoff values
qc_cutoff = {
    'rnaseq':{'PCT_PF_READS_ALIGNED':70,'PERCENT_DUPLICATION':30},
    'reseq':{'PCT_PF_READS_ALIGNED':70,'PERCENT_DUPLICATION':30},
    'WG-reseq':{'PCT_PF_READS_ALIGNED':70,'PERCENT_DUPLICATION':30},
    'seqcap':{'PCT_PF_READS_ALIGNED':70,'PERCENT_ON_TARGET':60, 'PCT_TARGET_BASES_10X':90, 'PERCENT_DUPLICATION':30},
    'customcap':{'PCT_PF_READS_ALIGNED':70, 'PERCENT_DUPLICATION':30},
    'finished':{},
    }

## Mapping from genomics project list application names
application_map = {'RNA-seq (Total RNA)':'rnaseq','WG re-seq':'WG-reseq','Resequencing':'reseq', 'Exome capture':'seqcap', 'Custom':'customcap', 'Finished library':'finished' , 'Custom capture':'customcap'}
application_inv_map = {v:k for k, v in application_map.items()}

def application_qc(project_id=None, flowcell_id=None, application=None,
                   user=None, password=None, url=None,
                   use_ps_map=True, use_bc_map=False, check_consistency=False, **kw):
    """Perform application specific qc on a project.

    :param project_id: project identifier
    :param flowcell_id: flowcell identifier
    :param application: application for which to perform qc
    :param user: database username
    :param password: database password
    :param url: database url
    :param use_ps_map: use project summary mapping
    :param use_bc_map: use barcode name mapping
    :param check_consistency: check 
    """
    LOG.debug("Doing application qc for project {}, flowcell {}".format(project_id, flowcell_id))
    header = ["sample","lane","flowcell", "date",  "TOTAL_READS",
              "MEAN_INSERT_SIZE", "GENOME_SIZE", "PERCENT_ON_TARGET",
              "PERCENT_DUPLICATION", "PCT_TARGET_BASES_10X", "PCT_PF_READS_ALIGNED",
              "dup_status", "status"]
    header_labels = ["sample","lane","flowcell", "date",  "#READS",
                     "INS_SIZE", "G_SIZE", "PCT_TGT",
                     "PCT_DUP", "PCT_10X", "PCT_ALN",
                     "dup_status", "status"]
    header_map = OrderedDict(zip(header, header_labels))

    output_data = {'stdout':StringIO(), 'stderr':StringIO()}
    p_con = ProjectSummaryConnection(username=user, password=password, url=url)
    project = p_con.get_entry(project_id)
    qc_data = p_con.get_qc_data(project_id, flowcell_id)
        
    if project.get("application") not in application_map.keys():
        if not application:
            LOG.warn("No such application {}. Please use the application option (available choices {})".format(app_label, ",".join(qc_cutoff.keys())))
            return
        application = application
    else:
        application = application_map[project.get("application")]

    def assess_qc(x, application):
        status = "PASS"
        dup_status = "OK"
        for k in qc_cutoff[application].keys():
            LOG.debug("assessing qc metric {}".format(k))
            if k == "PERCENT_DUPLICATION":
                if float(x[k]) > qc_cutoff[application][k]: 
                    dup_status = "HIGH"
            else:
                if float(x[k]) < qc_cutoff[application][k]:
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
    ## FIXME: this should be used also in assess_qc to set column widths
    col_widths = [20, 5, 11, 8, 11, 10, 11, 10, 10, 10, 10, 12, 12]
    ## Add info about qc
    output_data["stdout"].write("\n\nApplication QC\n")
    output_data["stdout"].write("==============\n")
    output_data["stdout"].write("{:12}{:>16}\n".format("Project", project_id.rstrip()))
    output_data["stdout"].write("{:12}{:>16}\n".format("Application", application_inv_map[application]))
    output_data["stdout"].write("==============\n\n")        
    output_data["stdout"].write("Application QC criteria\n")
    output_data["stdout"].write("==============\n")        
    for k,v in sorted(qc_cutoff[application].iteritems()):
        output_data["stdout"].write("{:<24}={:>4}\n".format(k, v))
    output_data["stdout"].write("==============\n\n") 
    output_data["stdout"].write("Header legend:\n")
    output_data["stdout"].write(", ".join(["{} = {}".format(v,k) for k, v in header_map.iteritems()][5:11]))
    output_data["stdout"].write("\n\n")
    header_fmt = ["{" + ":>{}".format(col_widths[i]) + "}" for i in range(0, len(col_widths))]
    header_fmt[0] = "{" + ":<{}".format(col_widths[0]) + "}"
    header_out = "".join([header_fmt[i].format(header_map[header_map.keys()[i]]) for i in range(0, len(col_widths))])
    output_data["stdout"].write(header_out + "\n")
    for k,v in sorted(qc_data.iteritems()):
        y = [str(x) for x in assess_qc(v, application)]
        output_data["stdout"].write("".join(y) + "\n")
    return output_data

def fastq_screen(project_id=None, flowcell_id=None,
                 user=None, password=None, url=None,
                 use_ps_map=True, use_bc_map=False, check_consistency=False, **kw):
    output_data = {'stdout':StringIO(), 'stderr':StringIO()}
    s_con = SampleRunMetricsConnection(username=user, password=password, url=url)
    samples = s_con.get_samples(fc_id=flowcell_id, sample_prj=project_id)
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
