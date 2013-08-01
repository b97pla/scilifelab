"""Demultiplex a CASAVA 1.8+ FastQ file based on the information in the fastq header 
assumed to be run in Undermined_indeces
"""

import sys
import os
import argparse
import collections
from scilifelab.utils.fastq_utils import FastQParser, is_read_pair, FastQWriter, parse_header
from scilifelab.illumina.hiseq import HiSeqSampleSheet


def demultiplex(CVS):
    sdata =  HiSeqSampleSheet(CVS)
    folderStruct = collections.defaultdict(lambda : collections.defaultdict(dict))
    outfiles = {}
    counts = {}
    
    for sd in sdata:
        sampleProject = sd['SampleProject']
        sampleID = sd['SampleID']
        lane = sd["Lane"]
        index = sd["Index"]
        if lane not in outfiles:
            outfiles[lane] = {}
            counts[lane] = {}
        outfiles[lane][index] = []
        counts[lane][index] = 0
        read = 1
        read1 = "tmp_{}_{}_L00{}_R{}_001.fastq.gz".format(sampleID,
                                                              index,
                                                              lane,
                                                              read)
        read = 2
        read2 = "tmp_{}_{}_L00{}_R{}_001.fastq.gz".format(sampleID,
                                                              index,
                                                              lane,
                                                              read)
        folderStruct[sampleProject][sampleID][lane] = [read1, read2]
        projectDirName = "Project_"+sampleProject
        sampleDirName = "Sample_"+sampleID
        outfiles[lane][index].append(os.path.join(projectDirName+ "/" + sampleDirName + "/" , read1))
        outfiles[lane][index].append(os.path.join(projectDirName+ "/" + sampleDirName + "/" , read2))
        
    for lane in outfiles:
        outfiles[lane]["Undetermined"] = [];
        read1 = "lane" + lane + "_Undetermined_L00" + lane + "_R1_001.fastq"
        read2 = "lane" + lane + "_Undetermined_L00" + lane + "_R2_001.fastq"
        outfiles[lane]["Undetermined"].append("Undetermined_indices/Sample_lane" + lane + "/" + read1)
        outfiles[lane]["Undetermined"].append("Undetermined_indices/Sample_lane" + lane + "/" + read2)
        folderStruct["Undetermined_indices"]["Sample_lane" + lane][lane] = [read1,read2]
    
    
    prepareDirectories(folderStruct)
    for lane in outfiles:
        for index in outfiles[lane]:
            files = outfiles[lane][index]
            outfiles[lane][index] = [] 
            for file in files:
                outfiles[lane][index].append(FastQWriter(file))
    
    demultiplex_lanes(outfiles)
    return


def prepareDirectories(folderStruct):
    for project in folderStruct:
        projectDirName = "Project_"+project
        if project is "Undetermined_indices":
            projectDirName = "Undetermined_indices"
        assert not os.path.exists(projectDirName), "Directory already exist, I am not going to overwrite previous results" 
        os.mkdir(projectDirName)
        os.chdir(projectDirName)
        for sample in folderStruct[project]:
            sampleDirName = "Sample_" + sample
            if project is "Undetermined_indices":
                sampleDirName = sample
            assert not os.path.exists(sampleDirName), "Directory already exist, I am not going to overwrite previous results" 
            os.mkdir(sampleDirName)
            os.chdir(sampleDirName)
            os.chdir("..")
        os.chdir("..")
    

def change_index(header,new_index):
    """Changes the index field in the header with the one specified in a  FASTQ header as specified by CASAVA 1.8.2 
        and returns the new header
       @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<index sequence>
    """
    if header[0] != '@':
        return None
    
    instrument, run_number, flowcell_id, lane, tile, x_pos, y_pos_read, is_filtered, control_number, index = header[1:].split(":")
    new_header = ":".join([instrument, run_number, flowcell_id, lane, tile, x_pos, y_pos_read, is_filtered, control_number, new_index])
    return new_header

def demultiplex_lanes(outfiles):
    for lane in outfiles:
        unmultiplexedDirForLane = "Sample_lane"+lane
        if os.path.exists(unmultiplexedDirForLane):
            unmultiplexedFileRead1  = "lane"+lane+"_Undetermined_L00"+lane+"_R1_001.fastq"
            unmultiplexedFileRead2  = "lane"+lane+"_Undetermined_L00"+lane+"_R2_001.fastq"
            print "working on lane {}".format(lane)
            demultiplex_lane(unmultiplexedDirForLane, unmultiplexedFileRead1, unmultiplexedFileRead2, outfiles)
        else:
            print 'Lane {0} with folder {1} not found'.format(lane, unmultiplexedDirForLane)
    return
    
def demultiplex_lane(unmultiplexedFolder, unmultiplexedFastq1, unmultiplexedFastq2, outfiles):
    fp1 = FastQParser(unmultiplexedFolder+"/"+unmultiplexedFastq1)
    fp2 = FastQParser(unmultiplexedFolder+"/"+unmultiplexedFastq2)
    for r1 in fp1:
        r2 = fp2.next()
        assert is_read_pair(r1,r2), "Mismatching headers for expected read pair" 
        header_r1 = parse_header(r1[0])
        lane_r1 = str(header_r1['lane'])
        index_r1 = header_r1['index']
        real_index = check_index(index_r1 , outfiles[lane_r1])
        
        r1[0] = change_index(r1[0], real_index)
        r2[0] = change_index(r2[0], real_index)
        
        outfiles[lane_r1][real_index][0].write(r1)
        outfiles[lane_r1][real_index][1].write(r2)
       
       


        
def check_index(index_r1 , filter):
    num_matches = 0
    index_matched = ""
    for index in filter:
        indexLength = len(index)
        index_to_compare = index_r1[0:indexLength]
        num_mismatches =  hamdist(index_to_compare, index)
        if num_mismatches <= 1:
            index_matched = index
            num_matches += 1
    if num_matches == 1:
        return index_matched
    elif num_matches == 0:
        return "Undetermined"
    else:
        return "ambiguity"
    

def hamdist(str1, str2):
    """Count the # of differences between equal length strings str1 and str2"""
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
       if ch1 != ch2:
           diffs += 1
           
    return diffs
    
def main():
    
    parser = argparse.ArgumentParser(description="Demultiplex a CASAVA 1.8+ FastQ file based on the information in the fastq header")
    
    
    parser.add_argument('CVS', action='store',
                        help="CVS file describing the experiment")
    
    args = parser.parse_args()
    
    demultiplex(args.CVS)



if __name__ == "__main__":
    main()
        