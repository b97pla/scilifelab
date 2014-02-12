"""

  :sequencing_report.py
  
  Module for creating reports based on data, extracted from StatusDB, and render them in rst format based on mako templates.
  The module contains ReportObject classes that act as connectors on top of StatusDB documents and are used for retrieving 
  field values as well as performing aggregating operations and gathering statistics. 
  
  The general idea is that the objects are instantiated and then passed as parameters to the method that renders the reports. 
  The mako templates contain the layout and content of the reports and receive the connector objects when rendering. The retrieval 
  of fields and statistics are done from the mako templates, although these should contain a minimum of logic. The logic is in
  the connector objects.
  
"""

import copy
import datetime
import numpy as np

from scilifelab.db.statusdb import ProjectSummaryConnection, SampleRunMetricsConnection, FlowcellRunMetricsConnection, ProjectSummaryDocument
from scilifelab.report.rst import render_rest_note, write_rest_note, _make_rst_table

def report(ctrl, project_name, **kw):
    """Main method for generating a report summarizing the progress of a project. Renders rst reports, concatenates them
    and write them to a file.
    
    :param ctrl: a Controller object 
    :param project_name: name of the project to generate report for (e.g. J.Doe_11_01)
    :param **kw: additional arguments
    
    """
    
    notes = project_summary_report(ctrl,project_name,**kw)
    write_rest_note("test.rst","/Users/pontus/Downloads/",contents=notes)

def project_summary_report(ctrl, project_name, sample_name=None, flowcell_id=None, **kw):
    """Generate and render a summary report for a project, consisting of a general summary, a
    flowcell-centric summary and one sample-centric summary for each sample.
    
    Instantiates a ProjectReportObject and passes them to the rendering methods.
    
    :param ctrl: the DeliveryReportController that invoked the method
    :param project_name: the project name (e.g. J.Doe_11_01)
    
    :returns: a list of rendered summary reports in rst format 
    """
    
    # Get connections to the databases in StatusDB
    pcon = ProjectSummaryConnection(**kw)
    assert pcon, "Could not connect to {} database in StatusDB".format("project")
    fcon = FlowcellRunMetricsConnection(**kw)
    assert fcon, "Could not connect to {} database in StatusDB".format("flowcell")
    scon = SampleRunMetricsConnection(**kw)
    assert scon, "Could not connect to {} database in StatusDB".format("samples")
    
    # Fetch the database document for the project 
    project = pcon.get_entry(project_name)
    if not project:
        ctrl.log.warn("No such project '{}'".format(project_name))
        return None
    # If the data source is not lims, we cannot generate a report using this method
    if project.get('source') != 'lims':
        ctrl.log.warn("The source for data for project {} is not LIMS. This type of report can only be generated for data from LIMS. Please refer to older versions of pm".format(project_name))
        return None
    
    # Instantiate a ProjectReportObject. This will be passed to the mako render method and used within the template
    psr = ProjectReportObject(ctrl, project, fcon, scon)
    
    # Render the project overview
    project_report = render_rest_note(tables={}, 
                                      report="project_report", 
                                      **{"project": psr,
                                         "skip_samples": sample_name,
                                         "skip_flowcells": flowcell_id})
    # Render the sample-centric reports, one for each sample in the project
    sample_reports = [render_rest_note(tables={}, 
                                       report="sample_report", 
                                       **{"project": psr,
                                          "sample": sample,
                                          "skip_samples": sample_name,
                                          "skip_flowcells": flowcell_id}) for sample in psr.project_samples]
    # Render the flowcell-centric reports, one for each flowcell run within the scope of the project
    flowcell_reports = [render_rest_note(tables={}, 
                                         report="flowcell_report", 
                                         **{"project": psr,
                                            "flowcell": fc,
                                            "skip_samples": sample_name,
                                            "skip_flowcells": flowcell_id}) for fc in psr.flowcells]
    # Return the list of rendered reports
    return [project_report] + sample_reports + flowcell_reports
   
class ReportObject():
    """The base class for ReportObjects. Specialized ReportObjects inherit this.
    """

    def __init__(self, ctrl, dbobj, dbcon=None):
        self.ctrl = ctrl
        self.log = self.ctrl.log
        self.dbobj = dbobj
        self.dbcon = dbcon
        self.date_format = "%Y-%m-%d"
        
    def __getattr__(self, name):
        return self.dbobj.get(name)
        
    def __repr__(self):
        return "<ReportObject>"
    
    def _parse_date(self, val, fmt=None):
        if fmt is None:
            fmt = self.date_format
        org = val
        try:
            val = datetime.datetime.strptime(val,fmt)
        except TypeError, ValueError:
            val = org
        return val

    def _date_field(self, field, fmt=None):
        if fmt is None:
            fmt = self.date_format
        try:
            return self._parse_date(getattr(self,field)).strftime(fmt)
        except AttributeError:
            return "N/A"

class ProjectReportObject(ReportObject):
    """The ProjectReportObject mostly represents a project document in StatusDB. It contains access methods to fields and 
    convenience methods for e.g. creating rst tables. The object also serves as an entry point to ReportObjects 
    representing e.g. flowcells and samples belonging to the project.
    """
    
    def __init__(self, ctrl, dbobj, fcon=None, scon=None):
        """
        
        :param ctrl: reference to a Controller object that holds e.g. reference to log and config
        :param dbobj: the database document representing the project
        :param fcon: a database connection to the flowcell database, used to fetch objects for flowcells relevant to the project
        :param scon: a database connection to the samples database, used to fetch objects for sample runs performed in the project
        
        """ 
        ReportObject.__init__(self,ctrl,dbobj,fcon)
        self.details = self.dbobj.get("details",{})
        self.project_samples = sorted([ProjectSampleReportObject(ctrl,sample,self.project_name,scon) for sample in self.dbobj.get("samples",{}).values()], key=lambda name: name.scilife_name)
        self.project_flowcells()     
        
    def __repr__(self):
        return "<ProjectReportObject {}>".format(self.project_name)

    def __getattr__(self, name):
        return self.dbobj.get(name,self.details.get(name))

    def sequencing_units_ordered(self):
        return self.__getattr__("sequence_units_ordered_(lanes)")
    def order_received(self):
        return self._date_field("order_received")
    def contract_received(self):
        return self._date_field("contract_received")
    def samples_received(self):
        return self._date_field("samples_received")
    def queued_date(self):
        return self._date_field("queued")
    def project_sample_name_table(self):
        return _make_rst_table([["SciLife ID", "Submitted ID"]] + [sample.project_sample_name_table_data() for sample in self.project_samples])
    def project_sample_status_table(self):
        return _make_rst_table([["SciLife ID", "Arrival QC", "Library QC", "Status", "M read( pair)s sequenced"]] + [sample.sample_status_table_data() for sample in self.project_samples])
    def project_flowcells(self):
        """Fetch flowcell documents for project-related flowcells
        """
        # skip if we have already fetched the documents
        if self.flowcells and len(self.flowcells) > 0:
            return
        self.flowcells = []
        lane_fc = sorted(list(set([fc for sample in self.project_samples for fc in sample.sample_run_flowcells()])))
        self.flowcell_lanes = {}
        for lane, fc in [(lfc.split("_")[0],"_".join(lfc.split("_")[1:3])) for lfc in lane_fc]:
            if fc not in self.flowcell_lanes:
                self.flowcell_lanes[fc] = []
            self.flowcell_lanes[fc].append(lane)
             
        for fcid in self.flowcell_lanes.keys():
            fcdoc = self.dbcon.get_entry(fcid)
            if not fcdoc:
                self.log.warn("Could not find flowcell document for {}".format(fcid))
                continue
            self.flowcells.append(FlowcellReportObject(self.ctrl,fcdoc))
            
    def project_flowcell_summary_table_data(self):
        if not self.flowcells:
            self.project_flowcells()        
        return [[fc.Date,
                 fc.Barcode,
                 fc.ScannerID,
                 fc.run_setup,
                 fc.FCPosition,
                 lane,
                 fc.lane_clusters(lane)/1e6,
                 str(fc.q30(lane)),
                 str(fc.phix(lane))] for fc in self.flowcells for lane in self.flowcell_lanes[fc.name]]
    
    def project_flowcell_summary_table(self):
        return _make_rst_table([["Date","Flowcell","Instrument","Run setup","Position","Lane","Clusters (M)","Q30 (%)","PhiX (%)"]] + self.project_flowcell_summary_table_data())
        
class ProjectSampleReportObject(ReportObject):
    
    def __init__(self, ctrl, dbobj, project_name, scon=None):
        ReportObject.__init__(self,ctrl,dbobj,scon)
        self.details = self.dbobj.get("details",{})
        self.library_prep = self.dbobj.get("library_prep",{})
        self.project_name = project_name
        self._sample_run_objs()
    
    def __repr__(self):
        return "<ProjectSampleReportObject {}>".format(self.scilife_name)

    def __getattr__(self, name):
        return self.dbobj.get(name,self.details.get(name))
    
    def _sample_run_objs(self):
        """Connect to the samples database and fetch sample_run_metrics documents for sequencing runs belonging to this sample
        """
        # Get all sample runs for the project and create objects 
        self.sample_run_objs = {sample_run_doc.get("_id"):SampleRunReportObject(self.ctrl,sample_run_doc) for sample_run_doc in self.dbcon.get_project_sample(self.scilife_name, sample_prj=self.project_name)}
        ids = self.sample_run_objs.keys()
        
        # Match the document ids to the sample_run_metrics entries under the corresponding prep
        for prep in self.library_prep.keys():
            for lane_run_data in self.library_prep[prep].get("sample_run_metrics",{}).values():
                doc_id = lane_run_data.get("sample_run_metrics_id")
                try:
                    ids.remove(doc_id)
                    self.sample_run_objs[doc_id].library_prep = prep
                except ValueError:
                    self.log.warn("No sample_run_metrics document for run {} with id {} found in samples database!".format(lane_run,doc_id))
        for id in ids:
            self.log.warn("sample_run_metrics document found for run {} but no corresponding entry in projects database. Please check for inconsistencies!".format(self.sample_run_objs[id].get("name","N/A")))
            self.sample_run_objs[id].library_prep = "N/A"
        self.sample_run_objs = self.sample_run_objs.values() 
            
    def library_prep_status(self):
        """Go through all preps and check if any of them have passed"""
        return "PASS" if any([prep.get("prep_status","").lower() == "passed" for prep in self.library_prep.values()]) else "FAIL"
    def library_prep_table_data(self):
        return [[k,
                 self.library_prep_fragment_size(k),
                 ", ".join(self.library_prep[k].get("reagent_labels",[])),
                 self.library_prep[k].get("prep_status")] for k in sorted(self.library_prep.keys())]
    def library_prep_table(self):
        return _make_rst_table([["Library prep", "Average fragment size (bp)", "Index", "Library validation"]] + self.library_prep_table_data())
    def library_prep_fragment_size(self, prep):
        validation = self._latest_library_validation(prep)
        if not validation:
            return None
        return validation.get("average_size_bp")
    def sample_run_flowcells(self, prep=None):
        return list(set(["_".join(k.split("_")[0:3]) for lbl, data in self.library_prep.items() for k in data.get("sample_run_metrics",{}).keys() if not prep or prep == lbl]))
    def project_sample_name_table_data(self):
        return [self.scilife_name, self.customer_name]
    def sample_status_table_data(self):
        return [self.scilife_name, 
                self.initial_QC_status or "N/A", 
                self.library_prep_status(), 
                self.status or "N/A", 
                self.m_reads_sequenced]

    def _latest_library_validation(self, prep):
        """Determines the latest library validation based on 1) date, 2) lims id
        """
        if prep not in self.library_prep:
            self.log.warn("Library prep {} is not in the list of library preps for sample {}".format(prep,self.scilife_name))
            return None
        
        validations = self.library_prep[prep].get("library_validation")
        if not validations:
            self.log.warn("Library prep {} for sample {} does not have any library validation".format(prep,self.scilife_name))
            return None
            
        # Sort first by the lims id which will break ties when sorting on date
        keys = sorted(validations.keys(), reverse=True)
        # Then  sort by starting date
        keys = sorted(keys, key=lambda k: self._parse_date(validations[k].get("start_date")), reverse=True)
        # Return the most recent validation
        return validations[keys[0]]

    def sample_run_table(self, project, flowcell=None):
        return _make_rst_table([["Library prep", "Sequencing started", "Flowcell ID", "Flowcell position", "Lane", "Reads (M)", "Yield (Mbases)", "Average Q", "Q30 (%)", "Perfect index read (%)"]] + [summary for prep in sorted(self.library_prep.keys()) for summary in self.sample_run_table_data(project, prep, flowcell)])
        
    def sample_run_table_data(self, project, prep=None, flowcell=None):
        """Compile the results related to a sample run by combining information from the sample object and the flowcell object
        """
        sample_run_data = []
        # Loop over the preps and compile data for each of the sample runs
        for sprep in self.library_prep.keys():
            if prep and prep != sprep:
                continue
            # Get the list of flowcells for which there exists sample_run_metrics sections in the project database entry
            sample_runs = [[sr.split("_")[0],"_".join(sr.split("_")[1:3])] for sr in self.sample_run_flowcells(sprep)]
            # If we are restricting to a particular flowcell, prune the sample_runs list
            if flowcell:
                sample_runs = [sr for sr in sample_runs if sr[1] == flowcell]
            # Compile information for each flowcell in the project
            for lane, run in sorted(sample_runs, key=lambda s: (s[1],s[0])):
                [fcobj] = [fc for fc in project.flowcells if fc.name == run]
                if fcobj is None:
                    self.log.warn("Could not find flowcell object corresponding to sample run {} for sample {}, please check database/python code!".format(run,self.scilife_name))
                    continue
                splobj = fcobj.flowcell_samples.get(self.scilife_name,{}).get(lane)
                if splobj is None:
                    self.log.warn("Could not find sample flowcell object corresponding to sample run {} for sample {}, please check database/python code!".format(run,self.scilife_name))
                    continue
                for spl in splobj.values():
                    sample_run_data.append([sprep,
                                            fcobj.Date,
                                            fcobj.Flowcell,
                                            fcobj.FCPosition,
                                            spl.lane,
                                            int(spl.reads)/1e6,
                                            spl.yield_bases,
                                            spl.avgQ,
                                            spl.q30,
                                            spl.mm0])
        return sample_run_data

    def sample_delivery_table(self):
        return _make_rst_table([["Library prep", "Flowcell ID", "Lane", "Read", "File size (bytes)", "Checksum", "Path"]] + [f for obj in sorted(self.sample_run_objs, key=lambda o: "_".join([o.library_prep,o.flowcell,o.lane])) for f in obj.sample_delivery_table_data()])
        
        
class FlowcellReportObject(ReportObject):
    
    def __init__(self, ctrl, dbobj):
        ReportObject.__init__(self,ctrl,dbobj)
        self.RunParameters = self.dbobj.get("RunParameters",{})
        self.RunInfo = self.dbobj.get("RunInfo",{})
        self.run_summary = self.dbobj.get("illumina",{}).get("run_summary",{})
        self.demultiplex_stats = self.dbobj.get("illumina",{}).get("Demultiplex_Stats",{}).get("Barcode_lane_statistics",[])
        self.flowcell_samples = self._flowcell_samples(self.demultiplex_stats) 
        
    def __repr__(self):
        return "<FlowcellReportObject {}>".format(self.name)

    def __getattr__(self, name):
        return self.dbobj.get(name,
                              self.RunInfo.get(name,
                                               self.RunParameters.get(name)))
    def _flowcell_samples(self, stats):
        samples = {}
        for sample in stats:
            id = sample.get("Sample ID")
            lane = sample.get("Lane")
            index = sample.get("Index")
            if id not in samples:
                samples[id] = {}
            if lane not in samples[id]:
                samples[id][lane] = {}
            samples[id][lane][index] = DemultiplexSample(sample)
        return samples

    def _get_lane_summary(self, lane, field):
        """Get the run summary of the flowcell that was fetched from LIMS
        """
        # MiSeq lanes are referred to by letter
        lane = str(lane)
        if lane not in self.run_summary and lane == "1":
            lane = "A"
        return self.run_summary.get(lane,{}).get(field)
    
    def lanes(self):
        return int(self.RunInfo.get("FlowcellLayout",{}).get("LaneCount","0"))
    
    def template_reads(self):
        return len([r for r in self.RunInfo.get("Reads",[]) if r.get("IsIndexedRead") == "N"])
    
    def index_reads(self):
        return len([r for r in self.RunInfo.get("Reads",[]) if r.get("IsIndexedRead") == "Y"])
    
    def index_setup(self):
        return "DUAL" if self.index_reads() == 2 else "SINGLE" if self.index_reads() == 1 else "No index"
    
    def lane_clusters(self, lane):
        return self._avg_lane_read_summary_field(lane,"Clusters PF")

    def lane_clusters_raw(self, lane):
        return self._avg_lane_read_summary_field(lane,"Clusters Raw")

    def lane_cluster_density(self, lane):
        return self._avg_lane_read_summary_field(lane,"Cluster Density (K/mm^2)")

    def prephasing(self, lane):
        return self._avg_lane_read_summary_field(lane,"% Prephasing")

    def phasing(self, lane):
        return self._avg_lane_read_summary_field(lane,"% Phasing")

    def lane_yield(self, lane):
        return self._avg_lane_read_summary_field(lane,"Yield PF (Gb)")

    def q30(self, lane):
        return self._avg_lane_read_summary_field(lane,"% Bases >=Q30")

    def avgQ(self, lane):
        return self._avg_lane_read_summary_field(lane,"Avg Q Score")

    def phix(self, lane):
        return self._avg_lane_read_summary_field(lane,"% Error Rate")

    def qc(self, lane):
        return self._get_lane_summary(lane,"qc")
    
    def demultiplex_recovery(self, lane):
        total = 0
        demuxed = 0
        for sample in self.flowcell_samples.values():
            if lane not in sample:
                continue
            total += np.sum([int(s.reads) for s in sample[lane].values()])
            demuxed += np.sum([int(s.reads) for s in sample[lane].values() if s.index != "Undetermined"])
        return 100.*demuxed/total
    
    def flowcell_lane_summary_table_data(self, lane):
        return [lane, 
                self.lane_cluster_density(lane),
                self.lane_clusters_raw(lane)/1e6, 
                self.lane_clusters(lane)/1e6, 
                self.lane_yield(lane),
                self.demultiplex_recovery(lane),
                self.prephasing(lane),
                self.phasing(lane),
                self.avgQ(lane),
                self.q30(lane),
                self.phix(lane),
                self.qc(lane)]
    
    def _avg_lane_read_summary_field(self, lane, field):
        return np.average([float(v) for v in self._lane_read_summary_field(lane,field)])
    
    def _lane_read_summary_field(self,lane,field):
        return [self._get_lane_summary(lane,"{} R{}".format(field,str(r+1))) for r in xrange(self.template_reads())]
    
    def flowcell_lane_summary_table(self, lanes):
        return _make_rst_table([["Lane", "Cluster density (K/mm2)", "Clusters Raw (M)", "Clusters PF (M)", "Yield (Gb)", "Demultiplex recovery (%)", "Prephasing (%)", "Phasing (%)", "Average Q", "Q30 (%)", "PhiX(%)", "QC"]] + [self.flowcell_lane_summary_table_data(lane) for lane in lanes])
    
    def flowcell_lane_sample_summary_table(self, lanes, samples):
        return _make_rst_table([["Lane", "Sample", "Index", "Reads (M)", "Yield (Mbases)", "Perfect index read (%)", "Percentage of lane (%)", "Average Q", "Q30 (%)"]] + [summary for lane in lanes for summary in self.flowcell_lane_sample_summary_table_data(lane, samples)])
    
    def flowcell_lane_sample_summary_table_data(self, lane, samples):
        summary = []
        for id in samples + ["lane{}".format(lane)]:
            if id not in self.flowcell_samples:
                continue
            if lane not in self.flowcell_samples[id]:
                continue
            for sample in self.flowcell_samples[id][lane].values():
                summary.append([sample.lane,
                                sample.scilife_name,
                                sample.index,
                                int(sample.reads)/1e6,
                                sample.yield_bases,
                                sample.mm0,
                                sample.pctLane,
                                sample.avgQ,
                                sample.q30])
        return summary
 
class DemultiplexSample():
    
    def __init__(self, data):
        self.data = data
        self.mapping = {"scilife_name": "Sample ID",
                        "index": "Index",
                        "lane": "Lane",
                        "description": "Description",
                        "q30": "% of >= Q30 Bases (PF)",
                        "reference_genome": "Sample Ref",
                        "reads": "# Reads",
                        "avgQ": "Mean Quality Score (PF)",
                        "project_name": "Project",
                        "pctPF": "% PF",
                        "yield_bases": "Yield (Mbases)",
                        "control": "Control",
                        "pctLane": "% of raw clusters per lane",
                        "mm0": "% Perfect Index Reads",
                        "mm1": "% One Mismatch Reads (Index)"}
        
    def __getattr__(self, name): 
        return self.data.get(self.mapping[name],"").replace("__",".").replace(",","")
        
class SampleRunReportObject(ReportObject):
    
    def __init__(self, ctrl, dbobj):
        ReportObject.__init__(self,ctrl,dbobj)
    
    def sample_delivery_table_data(self):
        if self.raw_data_delivery is None:
            return [[self.library_prep,
                     self.flowcell, 
                     self.lane,
                     "","","",""]]
                     
        return [[self.library_prep,
                 self.flowcell, 
                 self.lane, 
                 "{}".format("Forward" if r == "R1" else "Reverse"), 
                 self.raw_data_delivery["files"][r].get("size_in_bytes"), 
                 self.raw_data_delivery["files"][r].get("md5"), 
                 self.raw_data_delivery["files"][r].get("path")] for r in sorted(self.raw_data_delivery.get("files",{}).keys())]
        
          