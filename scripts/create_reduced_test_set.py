"""
Creates a minimal raw sequencing data test set from a flowcell directory by copying a
subset of the tiles and cycles and just the files needed for downstream analysis, as 
well as modifying the relevant config files to reflect this.
    
Usage:
    create_reduced_test_set.py <flow cell dir> [options]
    
Options:
    
    -t, --tile_subset=<FLOAT | INT>       Specify the size of the subset of tiles to use. 
                                           0 < FLOAT < 1 or 1 <= INT <= NUM_TILES. If < 1, 
                                           the specified fraction of tiles (rounded down) 
                                           will be used. If >= 1, the specified number of 
                                           tiles will be used. 
    
    -c, --cycle_subset=<FLOAT | INT>      Specify the size of the subset of cycles to use. 
                                           0 < FLOAT < 1 or 1 <= INT <= NUM_CYCLES. If < 1, 
                                           the specified fraction of cycles (rounded down) 
                                           will be used. If >= 1, the specified number of 
                                           cycles will be used. 
                                  
    -o, --output_dir=<output directory>    The destination directory where the minimal test
                                           set will be written to a sub-directory with the
                                           same name as the source directory. There must not
                                           already exist a folder for the flowcell within the 
                                           target directory.
                    

"""

import sys
import glob
import os
import re
import math
from optparse import OptionParser
import xml.etree.ElementTree as xml
from random import randint
from shutil import (copy2,copytree,rmtree)

def main(run_dir,tile_subset,cycle_subset,target_dir):
    
    # Get the name of the run folder
    _,run_name = os.path.split(os.path.normpath(run_dir))

    # Change to the target directory if specified
    cwd = os.getcwd()
    if target_dir is not None:
        target_dir = os.path.abspath(target_dir)
        assert os.path.exists(target_dir)
        os.chdir(target_dir)
        
    # Create the folder hierarchy to the intensity files
    assert not os.path.exists(run_name)
    os.mkdir(run_name,0770)
    os.chdir(run_name)

    # Copy the data files
    copy_data_files(run_dir,tile_subset,cycle_subset)
    # Copy the meta files
    print "Copying meta files"
    copy_meta_files(run_dir)
    
    # Change back to the original cwd
    os.chdir(cwd)
    
def copy_data_files(run_dir,tile_subset,cycle_subset):

    # Copy and parse the config.xml file to get the available tiles common to all lanes
    config_file = os.path.join(run_dir,'Data','Intensities','config.xml')
    all_tiles, cycle_range = _parse_config(config_file)
    
    # Randomly pick a subset of tiles to copy
    subset = []
    if tile_subset < 1.0:
        desired = max(1,int(math.floor(tile_subset*len(all_tiles))))
    else:
        desired = min(math.floor(tile_subset),len(all_tiles))
    
    while len(subset) < desired:
        subset.append(all_tiles.pop(randint(1,len(all_tiles))-1))

    print "%s tiles available, will copy a subset of %s tiles" % (len(all_tiles),len(subset))
    
    if cycle_subset > 0:
        for range in cycle_range:
            range[1] = min(range[1],int(range[0] + cycle_subset - 1))
    
    print "Will copy cycles in the ranges: %s" % cycle_range
    
    to_copy = []
    size = 0
    idir = os.path.join(run_dir,'Data','Intensities')
    bcdir = os.path.join(run_dir,'Data','Intensities','BaseCalls')
    for tile in subset:
        globs = [os.path.join(idir,"*%s*" % tile),
                 os.path.join(idir,"L*","*%s*" % tile),
                 os.path.join(idir,"L*","C*","*%s*" % tile),
                 os.path.join(idir,"Offsets","*%s*" % tile),
                 os.path.join(bcdir,"*%s*" % tile),
                 os.path.join(bcdir,"L*","*%s*" % tile),
                 os.path.join(bcdir,"L*","C*","*%s*" % tile),
                 os.path.join(bcdir,"Matrix", "*%s*" % tile),
                 os.path.join(bcdir,"Phasing", "*%s*" % tile)]
        for g in globs:
            for item in glob.glob(g):
                if _file_in_cycle(item,cycle_range):
                    to_copy.append([item,os.path.relpath(os.path.dirname(item),run_dir)])
                    size += os.path.getsize(to_copy[-1][0])
                
    print "Will copy %s files, totalling %s Mbytes" % (len(to_copy),size/(1024*1024))
    
    # Copy the files, creating subdirectories as necessary
    for src_file, dst_dir in to_copy:
        if not os.path.exists(dst_dir):
            print "Creating directory %s" % dst_dir
            os.makedirs(dst_dir,0770)
        print "Copying file %s to %s" % (src_file,dst_dir)
        copy2(src_file,dst_dir)
    
    config_file_cpy = os.path.join('Data','Intensities',os.path.basename(config_file))
    print "Copying and updating config file %s" % config_file
    _update_config(config_file,config_file_cpy,subset,cycle_range)

    config_file = os.path.join(run_dir,'Data','Intensities','BaseCalls','config.xml')
    config_file_cpy = os.path.join('Data','Intensities','BaseCalls',os.path.basename(config_file))
    print "Copying and updating config file %s" % config_file
    _update_config(config_file,config_file_cpy,subset,cycle_range)

def copy_meta_files(run_dir):
    
    bcdir = os.path.join('Data','Intensities','BaseCalls')
    globs = ['*.xml',
             '*.csv',
             '*.txt',
             os.path.join(bcdir,'*.xsl'),
             os.path.join(bcdir,'*.htm'),
             os.path.join(bcdir,'Plots'),
             os.path.join('Data','reports'),
             os.path.join('Data','Status_Files'),
             os.path.join('Data','Status.htm'),
             'InterOp']
    
    for g in globs:
        path = os.path.join(run_dir,g)
        print "Copying %s..." % path
        for item in glob.glob(path):
            relpath = os.path.relpath(item,run_dir)
            if os.path.isdir(item):
                # Remove the directory if it already exists
                if os.path.exists(relpath):
                    rmtree(relpath)
                copytree(item,relpath)
            elif os.path.isfile(item):
                copy2(item,relpath)
    
def _file_in_subset(file,subset):
    return _file_tile(file) in subset

def _file_tile(file):
    m = re.search(r"_(\d{4})",file)
    if m is not None and len(m.groups()):
        return int(m.group(1))
    return -1

def _file_in_cycle(file,cycle_range):
    m = re.search(r"C(\d+)\.1",os.path.dirname(file))
    if m is None or len(m.groups()) == 0:
        return True
    cycle = int(m.group(1))
    for range in cycle_range:
        if cycle >= range[0] and cycle <= range[1]:
            return True
    return False

def _parse_config(config_file):
    all_tiles = []
    cycle_range = []
    cfg = xml.parse(config_file)
    root = cfg.getroot()
    reads = root.findall('Run/RunParameters/Reads')
    for read in reads:
        c1 = int(read.find('FirstCycle').text)
        c2 = int(read.find('LastCycle').text)
        cycle_range.append([c1,c2])
        
    tiles = root.find('Run/TileSelection')
    if tiles is None:
        return all_tiles
    for lane in tiles.findall('Lane'):
        lane_tiles = []
        for tile in lane.findall('Tile'):
            lane_tiles.append(int(tile.text))
        if len(all_tiles) == 0:
            all_tiles = lane_tiles
        else:
            all_tiles = list(set(all_tiles).intersection(lane_tiles))
    return all_tiles, cycle_range

def _update_config(config_file_src,config_file_dst,subset,cycle_range):
    cfg = xml.parse(config_file_src)
    root = cfg.getroot()
    
    # Since we're not using an xml module with xpath capability, just traverse the entire tree of this (relatively small) xml file
    root = _update_cycle_range(root,cycle_range)
    cycles = root.find('Run/Cycles')
    if cycles is not None:
        num_cycles = 0
        last = 0
        for range in cycle_range:
            num_cycles += range[1] - range[0] + 1
            last = max(last,range[1])
            
        cycles.set('Last',str(last))
        cycles.set('Number',str(num_cycles))
        
    tiles = root.find('Run/TileSelection')
    if tiles is None:
        return
    for lane in tiles.findall('Lane'):
        for tile in lane.findall('Tile'):
            if int(tile.text) not in subset:
                lane.remove(tile)
    cfg.write(config_file_dst,encoding='UTF8')
    
def _update_cycle_range(elem,cycle_range):
    
    # Check if a child element matches 'FirstCycle'
    fc = elem.find('FirstCycle')
    if fc is not None:
        lc = elem.find('LastCycle')
        if lc is not None:
            for range in cycle_range:
                if int(fc.text) == range[0]:
                    lc.text = str(range[1])
                    
    # Check if the element has an attribute 'FirstCycle' and an element 'LastCycle'
    fc = elem.get('FirstCycle')
    if fc is not None and elem.get('LastCycle') is not None:
        for range in cycle_range:
            if int(fc) == range[0]:
                elem.set('LastCycle',str(range[1]))
                
    # Recurse into subelements
    for child in list(elem.getiterator(tag="*")):
        if child is None or child == elem:
            continue
        child = _update_cycle_range(child, cycle_range)
    
    return elem

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-t", "--tile_subset", dest="tile_subset", default=0.05)
    parser.add_option("-c", "--cycle_subset", dest="cycle_subset", default=100)
    parser.add_option("-o", "--output_dir", dest="output_dir", default=None)
    parser.add_option("-v", "--verbose", dest="verbose", default=False, \
                                                        action="store_true")
    options, args = parser.parse_args()
    if len(args) == 1:
        run_dir, = args
    else:
        print __doc__
        sys.exit()
    main(os.path.abspath(run_dir),float(options.tile_subset),float(options.cycle_subset),options.output_dir)
    