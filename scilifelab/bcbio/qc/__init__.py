"""bcbio qc module. Parsers for collecting qc metrics."""
import os
import re
import yaml
import glob
import xml.parsers.expat
from uuid import uuid4
import json
import numpy as np
import csv
import collections
import xml.etree.cElementTree as ET
from bs4 import BeautifulSoup
import datetime

from cement.core import backend
LOG = backend.minimal_logger("bcbio")

from bcbio.broad.metrics import PicardMetricsParser
from bcbio.pipeline.qcsummary import FastQCParser

class MetricsParser():
    """Basic class for parsing metrics"""
    def __init__(self, log=None):
        self.log = LOG
        if log:
            self.log = log

    def parse_bc_metrics(self, in_handle):
        data = {}
        while 1:
            line = in_handle.readline()
            if not line:
                break
            vals = line.rstrip("\t\n\r").split("\t")
            data[vals[0]] = int(vals[1])
        return data
    
    def parse_filter_metrics(self, in_handle):
        data = {}
        data["reads"] = int(in_handle.readline().rstrip("\n").split(" ")[-1])
        data["reads_aligned"] = int(in_handle.readline().split(" ")[-2])
        data["reads_fail_align"] = int(in_handle.readline().split(" ")[-2])
        return data

    def parse_fastq_screen_metrics(self, in_handle):
        in_handle.readline()
        data = {}
        while 1:
            line = in_handle.readline()
            if not line:
                break
            vals = line.rstrip("\t\n").split("\t")
            data[vals[0]] = {}
            data[vals[0]]["Unmapped"] = float(vals[1])
            data[vals[0]]["Mapped_One_Library"] = float(vals[2])
            data[vals[0]]["Mapped_Multiple_Libraries"] = float(vals[3])
        return data

    def parse_undemultiplexed_barcode_metrics(self, in_handle):

        data = collections.defaultdict(list)
        for line in in_handle:
            data[line['lane']].append({c:line[c] for c in in_handle.fieldnames if c != 'lane'})
        return data

    def parse_bcbb_checkpoints(self, in_handle):
        
        TIMEFORMAT = "%Y-%m-%dT%H:%M:%S.%f"
        timestamp = []
        for line in in_handle:
            try:
                ts = datetime.datetime.strptime(line.strip(), TIMEFORMAT)
                timestamp.append(ts)
            except ValueError:
                pass 
    
        return timestamp

class ExtendedPicardMetricsParser(PicardMetricsParser):
    """Extend basic functionality and parse all picard metrics"""

    def __init__(self):
        PicardMetricsParser.__init__(self)

    def _get_command(self, in_handle):
        while 1:
            line = in_handle.readline()
            if line.startswith("# net.sf.picard.analysis") or line.startswith("# net.sf.picard.sam"):
                break
        return line.rstrip("\n")

    def _read_off_header(self, in_handle):
        while 1:
            line = in_handle.readline()
            if line.startswith("## METRICS"):
                break
        return in_handle.readline().rstrip("\n").split("\t")

    def _read_vals_of_interest(self, header, info):
        want_indexes = [header.index(w) for w in header]
        vals = dict()
        for i in want_indexes:
            vals[header[i]] = info[i]
        return vals

    def _parse_align_metrics(self, in_handle):
        command = self._get_command(in_handle)
        header = self._read_off_header(in_handle)
        d = dict([[x, []] for x in header])
        res = dict(command=command, FIRST_OF_PAIR = d, SECOND_OF_PAIR = d, PAIR = d)
        while 1:
            info = in_handle.readline().rstrip("\n").split("\t")
            category = info[0]
            if len(info) <= 1:
                break
            vals = self._read_vals_of_interest(header, info)
            res[category] = vals
        return res

    def _parse_dup_metrics(self, in_handle):
        command = self._get_command(in_handle)
        header = self._read_off_header(in_handle)
        info = in_handle.readline().rstrip("\n").split("\t")
        vals = self._read_vals_of_interest(header, info)
        histvals = self._read_histogram(in_handle)
        return dict(command=command, metrics = vals, hist = histvals)

    def _parse_insert_metrics(self, in_handle):
        command = self._get_command(in_handle)
        header = self._read_off_header(in_handle)
        info = in_handle.readline().rstrip("\n").split("\t")
        vals = self._read_vals_of_interest(header, info)
        histvals = self._read_histogram(in_handle)
        return dict(command=command, metrics = vals, hist = histvals)

    def _parse_hybrid_metrics(self, in_handle):
        command = self._get_command(in_handle)
        header = self._read_off_header(in_handle)
        info = in_handle.readline().rstrip("\n").split("\t")
        vals = self._read_vals_of_interest(header, info)
        return dict(command=command, metrics = vals)
    
    def _read_histogram(self, in_handle):
        labels = self._read_to_histogram(in_handle)
        if labels is None:
            return None
        vals = dict([[x, []] for x in labels])
        while 1:
            line = in_handle.readline()
            info = line.rstrip("\n").split("\t")
            if len(info) < len(labels):
                break
            for i in range(0, len(labels)):
                vals[labels[i]].append(info[i])
        return vals

    def _read_to_histogram(self, in_handle):
        while 1:
            line = in_handle.readline()
            if line.startswith("## HISTOGRAM"):
                break
            if not line:
                return None
        return in_handle.readline().rstrip("\n").split("\t")
          
class RunInfoParser():
    """RunInfo parser"""
    def __init__(self):
        self._data = {}
        self._element = None

    def parse(self, fp):
        self._parse_RunInfo(fp)
        return self._data

    def _start_element(self, name, attrs):
        self._element=name
        if name == "Run":
            self._data["Id"] = attrs["Id"]
            self._data["Number"] = attrs["Number"]
        elif name == "FlowcellLayout":
            self._data["FlowcellLayout"] = attrs
        elif name == "Read":
            self._data["Reads"].append(attrs)
            
    def _end_element(self, name):
        self._element=None

    def _char_data(self, data):
        want_elements = ["Flowcell", "Instrument", "Date"]
        if self._element in want_elements:
            self._data[self._element] = data
        if self._element == "Reads":
            self._data["Reads"] = []

    def _parse_RunInfo(self, fp):
        p = xml.parsers.expat.ParserCreate()
        p.StartElementHandler = self._start_element
        p.EndElementHandler = self._end_element
        p.CharacterDataHandler = self._char_data
        p.ParseFile(fp)


# Generic XML to dict parsing
# See http://code.activestate.com/recipes/410469-xml-as-dictionary/
class XmlToList(list):
    def __init__(self, aList):
        for element in aList:
            if element:
                # treat like dict
                if len(element) == 1 or element[0].tag != element[1].tag:
                    self.append(XmlToDict(element))
                # treat like list
                elif element[0].tag == element[1].tag:
                    self.append(XmlToList(element))
            elif element.text:
                text = element.text.strip()
                if text:
                    self.append(text)
            else:
                # Set dict for attributes
                self.append({k:v for k,v in element.items()})

class XmlToDict(dict):
    '''
    Example usage:

    >>> tree = ET.parse('your_file.xml')
    >>> root = tree.getroot()
    >>> xmldict = XmlToDict(root)

    Or, if you want to use an XML string:

    >>> root = ET.XML(xml_string)
    >>> xmldict = XmlToDict(root)

    And then use xmldict for what it is... a dict.
    '''
    def __init__(self, parent_element):
        if parent_element.items():
            self.update(dict(parent_element.items()))
        for element in parent_element:
            if element:
                # treat like dict - we assume that if the first two tags
                # in a series are different, then they are all different.
                if len(element) == 1 or element[0].tag != element[1].tag:
                    aDict = XmlToDict(element)
                # treat like list - we assume that if the first two tags
                # in a series are the same, then the rest are the same.
                else:
                    # here, we put the list in dictionary; the key is the
                    # tag name the list elements all share in common, and
                    # the value is the list itself 
                    aDict = {element[0].tag: XmlToList(element)}
                # if the tag has attributes, add those to the dict
                if element.items():
                    aDict.update(dict(element.items()))
                self.update({element.tag: aDict})
            # this assumes that if you've got an attribute in a tag,
            # you won't be having any text. This may or may not be a 
            # good idea -- time will tell. It works for the way we are
            # currently doing XML configuration files...
            elif element.items():
                self.update({element.tag: dict(element.items())})
                # add the following line
                self[element.tag].update({"__Content__":element.text})

            # finally, if there are no child tags and no attributes, extract
            # the text
            else:
                self.update({element.tag: element.text})

class IlluminaXMLParser():
    """Illumina xml data parser. Parses xml files in flowcell directory."""
    def __init__(self):
        self._data = {}
        self._element = None
        self._tmp = None
        self._header = None

    def _chart_start_element(self, name, attrs):
        self._element = name
        if name == "FlowCellData":
            self._header = attrs
        if name == "Layout":
            n_tiles_per_lane = int(attrs['RowsPerLane']) * int(attrs['ColsPerLane'])
            if self._tmp is None:
                self._tmp = self._header
                self._tmp.update(attrs)
                for i in range(1,int(attrs['NumLanes'])+1):
                    for j in range(1,n_tiles_per_lane+1):
                        key = "%s_%s" % (i,j)
                        self._tmp[key] = {}
                     
        if name == "TL":
            self._tmp[attrs["Key"]][self._index] = {}
            for k in attrs.keys():
                if k == "Key":
                    continue
                if attrs[k] == "NaN":
                    v = None
                else:
                    v = float(attrs[k])
                                        
                self._tmp[attrs["Key"]][self._index] = v

    def _chart_end_element(self, name):
        self._element = None
    def _chart_char_data(self, data):
        pass

    def _parse_charts(self, files):
        for f in files:
            self._index = os.path.basename(f).rstrip(".xml").lstrip("Chart_")
            p = xml.parsers.expat.ParserCreate()
            p.StartElementHandler = self._chart_start_element
            p.EndElementHandler = self._chart_end_element
            p.CharacterDataHandler = self._chart_char_data
            fp = open(f)
            p.ParseFile(fp)
            fp.close()

    def _summary_start_element(self, name, attrs):
        self._element = name
        if name == "Summary":
            self._tmp[self._index] = attrs
        if name == "Lane":
            self._tmp[self._index][attrs['key']] = attrs
    def _summary_end_element(self, name):
        self._element = None
    def _summary_char_data(self, data):
        pass

    def _parse_summary(self, files):
        for f in files:
            self._index = os.path.basename(f).rstrip(".xml")
            p = xml.parsers.expat.ParserCreate()
            p.StartElementHandler = self._summary_start_element
            p.EndElementHandler = self._summary_end_element
            p.CharacterDataHandler = self._summary_char_data
            fp = open(f)
            p.ParseFile(fp)
            fp.close()

    def _clusters_start_element(self, name, attrs):
        self._element = name
        if name == "Data":
            self._tmp[self._index] = attrs
        if name == "Lane":
            self._tmp[self._index][attrs['key']] = attrs
    def _clusters_end_element(self, name):
        self._element = None
    def _clusters_char_data(self, data):
        pass

    def _parse_clusters(self, files):
        for f in files:
            self._index = os.path.basename(f).rstrip(".xml")
            p = xml.parsers.expat.ParserCreate()
            p.StartElementHandler = self._clusters_start_element
            p.EndElementHandler = self._clusters_end_element
            p.CharacterDataHandler = self._clusters_char_data
            fp = open(f)
            p.ParseFile(fp)
            fp.close()

    ## Caution: no assert statements for file existence
    def parse(self, files, fullRTA=False):
        """Full parsing includes all RTA files"""
        if fullRTA:
            error_files = filter(lambda x: os.path.dirname(x).endswith("ErrorRate"), files)
            self._parse_charts(error_files)
            self._data["ErrorRate"] = self._tmp
            self._tmp = None
            FWHM_files = filter(lambda x: os.path.dirname(x).endswith("FWHM"), files)
            self._parse_charts(FWHM_files)
            self._data["FWHM"] = self._tmp
            self._tmp = None
            intensity_files = filter(lambda x: os.path.dirname(x).endswith("Intensity"), files)
            self._parse_charts(intensity_files)
            self._data["Intensity"] = self._tmp
            self._tmp = None
            numgt30_files = filter(lambda x: os.path.dirname(x).endswith("NumGT30"), files)
            self._parse_charts(numgt30_files)
            self._data["NumGT30"] = self._tmp
            self._tmp = None
            chart_files = filter(lambda x: os.path.basename(x).endswith("_Chart.xml"), files)
            self._parse_charts(chart_files)
            self._data["Charts"] = self._tmp
            self._tmp = None

        ## Parse Summary and clusters
        self._tmp = {}
        summary_files = filter(lambda x: os.path.dirname(x).endswith("Summary"), files)
        self._parse_summary(summary_files)
        self._data["Summary"] = self._tmp
        self._tmp = {}
        cluster_files = filter(lambda x: os.path.basename(x).startswith("NumClusters By"), files)
        self._parse_clusters(cluster_files)
        self._data["NumClusters"] = self._tmp

        return self._data

class ExtendedFastQCParser(FastQCParser):
    def __init__(self, base_dir):
        FastQCParser.__init__(self, base_dir)

    def get_fastqc_summary(self):
        metric_labels = ["Per base sequence quality", "Basic Statistics", "Per sequence quality scores",
                         "Per base sequence content", "Per base GC content", "Per sequence GC content",
                         "Per base N content", "Sequence Length Distribution", "Sequence Duplication Levels",
                         "Overrepresented sequences", "Kmer Content"]
        metrics = {x : self._to_dict(self._fastqc_data_section(x)) for x in metric_labels}
        return metrics

    def _to_dict(self, section):
        if len(section) == 0:
            return {}
        header = [x.strip("#") for x in section[0].rstrip("\t").split("\t")]
        d = []
        for l in section[1:]:
            d.append(l.split("\t"))
        data = np.array(d)
        df = {header[i]:data[:,i].tolist() for i in range(0,len(header))}
        return df
##############################
##  objects
##############################
class RunMetricsParser(dict):
    """Generic Run Parser class"""
    _metrics = []
    ## Following paths are ignored
    ignore = "|".join(["tmp", "tx", "-split", "log"])
    reignore = re.compile(ignore)

    def __init__(self, log=None):
        super(RunMetricsParser, self).__init__()
        self.files = []
        self.path=None
        self.log = LOG
        if log:
            self.log = log

    def _collect_files(self):
        if not self.path:
            return
        if not os.path.exists(self.path):
            raise IOError
        self.files = []
        for root, dirs, files in os.walk(self.path):
            if re.search(self.reignore, root):
                continue
            self.files = self.files + [os.path.join(root, x) for x in files]

    def filter_files(self, pattern, filter_fn=None):
        """Take file list and return those files that pass the filter_fn criterium"""
        def filter_function(f):
            return re.search(pattern, f) != None
        if not filter_fn:
            filter_fn = filter_function
        return filter(filter_fn, self.files)

class SampleRunMetricsParser(RunMetricsParser):
    """Sample-level class for parsing run metrics data"""

    def __init__(self, path):
        RunMetricsParser.__init__(self)
        self.path = path
        self._collect_files()
        
    def read_picard_metrics(self, barcode_name, sample_prj, lane, flowcell, barcode_id, **kw):
        self.log.debug("read_picard_metrics for sample {}, project {}, lane {} in run {}".format(barcode_name, sample_prj, lane, flowcell))
        picard_parser = ExtendedPicardMetricsParser()
        pattern = "|".join(["{}_[0-9]+_[0-9A-Za-z]+(_nophix)?_{}-.*.(align|hs|insert|dup)_metrics".format(lane, barcode_id),
                            "{}_[0-9]+_[0-9A-Za-z]+_{}(_nophix)?-.*.(align|hs|insert|dup)_metrics".format(lane, barcode_id)])
        files = self.filter_files(pattern)
        if len(files) == 0:
            self.log.warn("no picard metrics files for sample {}; pattern {}".format(barcode_name, pattern))
            return {}
        try:
            self.log.debug("files {}".format(",".join(files)))
            metrics = picard_parser.extract_metrics(files)
            return metrics
        except:
            self.log.warn("no picard metrics for sample {}".format(barcode_name))
            return {}

    def parse_fastq_screen(self, barcode_name, sample_prj, lane, flowcell, barcode_id, **kw):
        self.log.debug("parse_fastq_screen for sample {}, project {}, lane {} in run {}".format(barcode_name, sample_prj, lane, flowcell))
        parser = MetricsParser()
        pattern = "|".join(["{}_[0-9]+_[0-9A-Za-z]+(_nophix)?_{}_[12]_screen.txt".format(lane, barcode_id),
                            "{}_[0-9]+_[0-9A-Za-z]+_{}(_nophix)?_[12]_screen.txt".format(lane, barcode_id)])
        files = self.filter_files(pattern)
        self.log.debug("files {}".format(",".join(files)))
        try:
            fp = open(files[0])
            data = parser.parse_fastq_screen_metrics(fp)
            fp.close()
            return data
        except:
            self.log.warn("no fastq screen metrics for sample {}".format(barcode_name))
            return {}

    def parse_bcbb_checkpoints(self, barcode_name, sample_prj, flowcell, barcode_id, **kw):
        self.log.debug("parse_bcbb_checkpoints for sample {}, project {} in run {}".format(barcode_name, sample_prj, flowcell))
        parser = MetricsParser()
        pattern = r"[0-9][0-9]_.*\.txt"
        files = self.filter_files(pattern)
        self.log.debug("files {}".format(",".join(files)))
        
        checkpoints = {}
        for f in files:
            try: 
                with open(f) as fh:
                    checkpoints[os.path.splitext(os.path.basename(f))[0]] = parser.parse_bcbb_checkpoints(fh)
            except Exception as e:
                self.log.warn("Exception: {}".format(e))
                self.log.warn("no bcbb checkpoint for sample {} using pattern '{}'".format(barcode_name, pattern))
        
        return checkpoints

    def read_fastqc_metrics(self, barcode_name, sample_prj, lane, flowcell, barcode_id, **kw):
        self.log.debug("read_fastqc_metrics for sample {}, project {}, lane {} in run {}".format(barcode_name, sample_prj, lane, flowcell))
        if barcode_name == "unmatched":
            return
        pattern = "fastqc/{}_[0-9]+_[0-9A-Za-z]+(_nophix)?_{}-*".format(lane, barcode_id)
        files = self.filter_files(pattern)
        self.log.debug("files {}".format(",".join(files)))
        try:
            fastqc_dir = os.path.dirname(files[0])
            fqparser = ExtendedFastQCParser(fastqc_dir)
            stats = fqparser.get_fastqc_summary()
            return {'stats':stats}
        except Exception as e:
            self.log.warn("Exception: {}".format(e))
            self.log.warn("no fastqc metrics for sample {} using pattern '{}'".format(barcode_name, pattern))
            return {'stats':{}}

    def parse_filter_metrics(self, **kw):
        """CASAVA: Parse filter metrics at sample level"""
        self.log.debug("parse_filter_metrics for lane {}, project {} in flowcell {}".format(lane, sample_prj, flowcell))
        pattern = "{}_[0-9]+_[0-9A-Za-z]+_{}(_nophix)?.filter_metrics".format(lane, barcode_id)
        files = self.filter_files(pattern)
        self.log.debug("files {}".format(",".join(files)))
        try:
            fp = open(files[0])
            parser = MetricsParser()
            data = parser.parse_filter_metrics(fp)
            fp.close()
            return data
        except:
            self.log.warn("No filter nophix metrics for lane {}".format(lane))
            return {"reads":None, "reads_aligned":None, "reads_fail_align":None}

    def get_bc_count(self, barcode_name, sample_prj, flowcell, lane, barcode_id, demultiplex_stats=None, **kw):
        """Parse bc metrics at sample level and get *bc_count* for a sample run!"""
        self.log.debug("get_bc_count for sample {}, project {} in flowcell {}".format(barcode_name, sample_prj, flowcell))
        # If demultiplex_stats passed use this info instead
        if demultiplex_stats:
            demux_stats_dict = {"{}_{}".format(l.get('Sample ID', None), l.get('Lane', None)):l for l in demultiplex_stats.get('Barcode_lane_statistics', [])}
            sample_lane = "{}_{}".format(barcode_name, lane)
            if sample_lane in demux_stats_dict:
                self.log.debug("sample {}, lane {} found in demultiplex_stats - using this information".format(barcode_name, lane))
                return int(demux_stats_dict[sample_lane]["# Reads"].replace(",", ""))/2
        pattern = "{}_[0-9]+_[0-9A-Za-z]+(_nophix)?[\._]bc[\._]metrics".format(lane)
        files = self.filter_files(pattern)
        if len(files) == 0:
            self.log.debug("no bc metrics files for sample {}, lane {}; pattern {}".format(barcode_name, lane, pattern))
            return None
        self.log.debug("files {}".format(",".join(files)))
        try:
            parser = MetricsParser()
            fp = open(files[0])
            data = parser.parse_bc_metrics(fp)
            fp.close()
            return data[str(barcode_id)]
        except:
            self.log.warn("No bc_metrics info for lane {}".format(lane))
            return None


class FlowcellRunMetricsParser(RunMetricsParser):
    """Flowcell level class for parsing flowcell run metrics data."""
    _lanes = range(1,9)
    def __init__(self, path):
        RunMetricsParser.__init__(self)
        self.path = path
        self._collect_files()

    def parseRunInfo(self, fn="RunInfo.xml", **kw):
        infile = os.path.join(os.path.abspath(self.path), fn)
        self.log.debug("parseRunInfo: going to read {}".format(infile))
        if not os.path.exists(infile):
            self.log.warn("No such file {}".format(infile))
            return {}
        try:
            fp = open(infile)
            parser = RunInfoParser()
            data = parser.parse(fp)
            fp.close()
            return data
        except:
            self.log.warn("Reading file {} failed".format(os.path.join(os.path.abspath(self.path), fn)))
            return {}

    def parseRunParameters(self, fn="runParameters.xml", **kw):
        """Parse runParameters.xml from an Illumina run.

        :param fn: filename
        :param **kw: keyword argument

        :returns: parsed data structure
        """
        infile = os.path.join(os.path.abspath(self.path), fn)
        self.log.debug("parseRunParameters: going to read {}".format(infile))
        if not os.path.exists(infile):
            self.log.warn("No such files {}".format(infile))
            return {}
        try:
            with open(infile) as fh:
                tree = ET.parse(fh)
            root = tree.getroot()
            data = XmlToDict(root)
            return data
        except:
            self.log.warn("Reading file {} failed".format(os.path.join(os.path.abspath(self.path), fn)))
            return {}

    def parse_samplesheet_csv(self, runinfo_csv="SampleSheet.csv", **kw):
        infile = os.path.join(os.path.abspath(self.path), runinfo_csv)
        self.log.debug("parse_samplesheet_csv: going to read {}".format(infile))
        if not os.path.exists(infile):
            self.log.warn("No such file {}".format(infile))
            return {}
        try:
            fp = open(infile)
            runinfo = [x for x in csv.DictReader(fp)]
            fp.close()
            return runinfo
        except:
            self.log.warn("Reading file {} failed".format(infile))
            return {}
            
    def parse_run_info_yaml(self, run_info_yaml="run_info.yaml", **kw):
        infile = os.path.join(os.path.abspath(self.path), run_info_yaml)
        self.log.debug("parse_run_info_yaml: going to read {}".format(infile))
        if not os.path.exists(infile):
            self.log.warn("No such file {}".format(infile))
            return {}
        try:
            fp = open(infile)
            runinfo = yaml.load(fp)
            fp.close()
            return runinfo
            return True
        except:
            self.log.warn("No such file {}".format(infile))
            return False

    def parse_illumina_metrics(self, fullRTA=False, **kw):
        self.log.debug("parse_illumina_metrics")
        fn = []
        for root, dirs, files in os.walk(os.path.abspath(self.path)):
            for f in files:
                if f.endswith(".xml"):
                    fn.append(os.path.join(root, f))
        self.log.debug("Found {} RTA files {}...".format(len(fn), ",".join(fn[0:10])))
        parser = IlluminaXMLParser()
        metrics = parser.parse(fn, fullRTA)
        return metrics

    def parse_filter_metrics(self, fc_name, **kw):
        """pre-CASAVA: Parse filter metrics at flowcell level"""
        self.log.debug("parse_filter_metrics for flowcell {}".format(fc_name))
        lanes = {str(k):{} for k in self._lanes}
        for lane in self._lanes:
            pattern = "{}_[0-9]+_[0-9A-Za-z]+(_nophix)?.filter_metrics".format(lane)
            lanes[str(lane)]["filter_metrics"] = {"reads":None, "reads_aligned":None, "reads_fail_align":None}
            files = self.filter_files(pattern)
            self.log.debug("filter metrics files {}".format(",".join(files)))
            try:
                fp = open(files[0])
                parser = MetricsParser()
                data = parser.parse_filter_metrics(fp)
                fp.close()
                lanes[str(lane)]["filter_metrics"] = data
            except:
                self.log.warn("No filter nophix metrics for lane {}".format(lane))
        return lanes

    def parse_bc_metrics(self, fc_name, **kw):
        """Parse bc metrics at sample level"""
        self.log.debug("parse_bc_metrics for flowcell {}".format(fc_name[1:]))
        lanes = {str(k):{} for k in self._lanes}
        for lane in self._lanes:
            pattern = "{}_[0-9]+_[0-9A-Za-z]+(_nophix)?[\._]bc[\._]metrics".format(lane)
            lanes[str(lane)]["bc_metrics"] = {}
            files = self.filter_files(pattern)
            self.log.debug("bc metrics files {}".format(",".join(files)))
            try:
                parser = MetricsParser()
                fp = open(files[0])
                data = parser.parse_bc_metrics(fp)
                fp.close()
                lanes[str(lane)]["bc_metrics"] = data
            except:
                self.log.warn("No bc_metrics info for lane {}".format(lane))
        return lanes

    def parse_undemultiplexed_barcode_metrics(self, fc_name, **kw):
        """Parse the undetermined indices top barcodes materics
        """
        
        lanes = {str(k):{} for k in self._lanes}
        # Use a glob to allow for multiple fastq folders
        metrics_file_pattern = os.path.join(self.path, "Unaligned*", "Basecall_Stats_{}".format(fc_name[1:]), "Undemultiplexed_stats.metrics")
        for metrics_file in glob.glob(metrics_file_pattern):
            self.log.debug("parsing {}".format(metrics_file))
            if not os.path.exists(metrics_file):
                self.log.warn("No such file {}".format(metrics_file))
                continue
            
            with open(metrics_file) as fh:
                parser = MetricsParser()
                in_handle = csv.DictReader(fh, dialect=csv.excel_tab)
                data = parser.parse_undemultiplexed_barcode_metrics(in_handle)
                for k in lanes.keys():
                    lanes[str(k)]["undemultiplexed_barcodes"] = collections.defaultdict(list)
                    try:
                        for barcode in data[str(k)]:
                            # Warn that we are replacing previously parsed results
                            if len(lanes[str(k)]) > 0:
                                self.log.warn("Conflicting undemultiplexed barcode metrics for lane {}. Using values parsed from {}".format(k,metrics_file))
                            for key, val in barcode.items():
                                lanes[str(k)]["undemultiplexed_barcodes"][key].append(val)
                    except KeyError:
                        self.log.warn("No undemultiplexed barcode metrics for lane {}".format(k))
        return lanes
    
    def parse_demultiplex_stats_htm(self, fc_name, **kw):
        """Parse the Unaligned*/Basecall_Stats_*/Demultiplex_Stats.htm file
        generated from CASAVA demultiplexing and returns barcode metrics.
        """
        metrics = {"Barcode_lane_statistics": [],
                   "Sample_information": []}
        # Use a glob to allow for multiple fastq directories
        htm_file_pattern = os.path.join(self.path, "Unaligned*", "Basecall_Stats_{}".format(fc_name[1:]), "Demultiplex_Stats.htm")
        for htm_file in glob.glob(htm_file_pattern):
            self.log.debug("parsing {}".format(htm_file))
            if not os.path.exists(htm_file):
                self.log.warn("No such file {}".format(htm_file))
                continue
            with open(htm_file) as fh:
                htm_doc = fh.read()
            soup = BeautifulSoup(htm_doc)
            ## 
            ## Find headers
            allrows = soup.findAll("tr")
            column_gen=(row.findAll("th") for row in allrows)
            parse_row = lambda row: row 
            headers = [h for h in map(parse_row, column_gen) if h]
            bc_header = [str(x.string) for x in headers[0]]
            smp_header = [str(x.string) for x in headers[1]]
            ## 'Known' headers from a Demultiplex_Stats.htm document
            bc_header_known = ['Lane', 'Sample ID', 'Sample Ref', 'Index', 'Description', 'Control', 'Project', 'Yield (Mbases)', '% PF', '# Reads', '% of raw clusters per lane', '% Perfect Index Reads', '% One Mismatch Reads (Index)', '% of >= Q30 Bases (PF)', 'Mean Quality Score (PF)']
            smp_header_known = ['None', 'Recipe', 'Operator', 'Directory']
            if not bc_header == bc_header_known:
                self.log.warn("Barcode lane statistics header information has changed. New format?\nOld format: {}\nSaw: {}".format(",".join((["'{}'".format(x) for x in bc_header_known])), ",".join(["'{}'".format(x) for x in bc_header])))
            if not smp_header == smp_header_known:
                self.log.warn("Sample header information has changed. New format?\nOld format: {}\nSaw: {}".format(",".join((["'{}'".format(x) for x in smp_header_known])), ",".join(["'{}'".format(x) for x in smp_header])))
            ## Fix first header name in smp_header since htm document is mal-formatted: <th>Sample<p></p>ID</th>
            smp_header[0] = "Sample ID"
    
            ## Parse Barcode lane statistics
            soup = BeautifulSoup(htm_doc)
            table = soup.findAll("table")[1]
            rows = table.findAll("tr")
            column_gen = (row.findAll("td") for row in rows)
            parse_row = lambda row: {bc_header[i]:str(row[i].string) for i in range(0, len(bc_header)) if row}
            metrics["Barcode_lane_statistics"].extend(map(parse_row, column_gen))
    
            ## Parse Sample information
            soup = BeautifulSoup(htm_doc)
            table = soup.findAll("table")[3]
            rows = table.findAll("tr")
            column_gen = (row.findAll("td") for row in rows)
            parse_row = lambda row: {smp_header[i]:str(row[i].string) for i in range(0, len(smp_header)) if row}
            metrics["Sample_information"].extend(map(parse_row, column_gen))
            
        ## Set data
        return metrics
