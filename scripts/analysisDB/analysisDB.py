#!/usr/bin/env python

import commands

"""A module for building up the best practice analysis objects that 
build up the analysis database on statusdb.

Maya Brandi, Science for Life Laboratory, Stockholm, Sweden.
"""

class BP_RNA():
    def __init__(self, project_name):
        self.obj = {'entity_type': 'RNA_BP_analysis',
                    'samples': {},
                    'project_name': project_name}
        self.get_stat()
        self.get_RSeQC()
        self.get_rRNA()
        self.get_picardDup()
        self.get_picardIns()
        self.get_top_dups()

    def get_proj_db_inf(self, proj_db):
        for key in ['application', 'project_id', 'no_of_samples']:
            if proj_db.has_key(key):
                self.obj[key] = proj_db[key]

    def get_stat(self):
        try:stat = open('stat.json','r')
        except: return
        stat_dict = eval(stat.read())
        scilife_names,preps = self.strip_scilife_name(stat_dict.keys())
        for samp in scilife_names:
            samp_name_stripped = scilife_names[samp]
            if self.obj['samples'].has_key(samp_name_stripped):
                self.obj['samples'][samp_name_stripped]['mapping_statistics'] = stat_dict[samp]
                self.obj['samples'][samp_name_stripped]['prep'] = preps[samp]
            else:
                self.obj['samples'][samp_name_stripped] = {'mapping_statistics' : stat_dict[samp],
                                                           'prep': preps[samp]}
    def get_rRNA(self):
        try:rRNA = open('rRNA.quantification','r')
        except: return
        for line in rRNA:
            samp = line.split()[0]
            scilife_names, preps = self.strip_scilife_name([samp])
            samp_name_stripped  = scilife_names[samp]
            if self.obj['samples'].has_key(samp_name_stripped):
                self.obj['samples'][samp_name_stripped]['percent_rRNA'] = float(line.split()[1].strip('%'))
            else:
                self.obj['samples'][samp_name_stripped] = {'percent_rRNA' : float(line.split()[1].strip('%'))}

    def get_RSeQC(self):
        try:RSeQC = open('RSeQC_rd.json','r')
        except: return
        RSeQC_dict = eval(RSeQC.read())
        scilife_names,preps = self.strip_scilife_name(RSeQC_dict.keys())
        for samp in scilife_names:
            samp_name_stripped = scilife_names[samp]
            if self.obj['samples'].has_key(samp_name_stripped):
                self.obj['samples'][samp_name_stripped]['read_distribution'] = RSeQC_dict[samp]
            else:
                self.obj['samples'][samp_name_stripped] = {'read_distribution' : RSeQC_dict[samp]}

    def get_top_dups(self):
        try:top_dups = open('top_dups.json','r')
        except: return 
        top_dups_dict = eval(top_dups.read())
        scilife_names,preps = self.strip_scilife_name(top_dups_dict.keys())
        for samp in scilife_names:
            samp_name_stripped = scilife_names[samp]
            if self.obj['samples'].has_key(samp_name_stripped):
                self.obj['samples'][samp_name_stripped]['top_dups'] = top_dups_dict[samp]
            else:
                self.obj['samples'][samp_name_stripped] = {'top_dups' : top_dups_dict[samp]}


    def get_picardDup(self):
        names = commands.getoutput("ls -d tophat_out_*|sed 's/tophat_out_//g'").split('\n')
        scilife_names, preps = self.strip_scilife_name(names)
        for samp in scilife_names:
            samp_name_stripped = scilife_names[samp]
            try:picardDup = self.pars_picard_metrics('tophat_out_' + samp + '/' + samp + '_picardDup_metrics')
            except: return
            if self.obj['samples'].has_key(samp_name_stripped):
                self.obj['samples'][samp_name_stripped]['picard_dup'] = picardDup
            else:
                self.obj['samples'][samp_name_stripped] = {'picard_dup' : picard_dup} 

    def get_picardIns(self):
        names = commands.getoutput("ls -d tophat_out_*|sed 's/tophat_out_//g'").split('\n')
        scilife_names, preps = self.strip_scilife_name(names)
        for samp in scilife_names:
            samp_name_stripped = scilife_names[samp]
            try:picardEstInSize = self.pars_picard_metrics('tophat_out_' + samp + '/' + samp + '.picard_estimated_insert_size')
            except:return
            if self.obj['samples'].has_key(samp_name_stripped):
                self.obj['samples'][samp_name_stripped]['picard_estimated_insert_size'] = picardEstInSize
            else:
                self.obj['samples'][samp_name_stripped] = {'picard_estimated_insert_size' : picardEstInSize}


    def pars_picard_metrics(self,file):
        f=open(file,'r')
        lines=f.readlines()
        for i,line in enumerate(lines):
            if len(line.split('\t'))>0 and line.split('\t')[0]=='## METRICS CLASS':
                keys=lines[i+1].strip().split('\t')
                vals=lines[i+2].strip().split('\t')
                return dict(zip(keys,vals))

    def pars_picard_histogram(self,file):
        f=open(file,'r')
        lines=f.readlines()
        histogram = False
        for i,line in enumerate(lines):
            if histogram:
                values = lines[i+1].strip().split('\t')
                if values[0] !='':
                    histogram[keys[0]].append(values[0])
                    histogram[keys[1]].append(values[1])
            if len(line.split('\t'))>0 and line.split('\t')[0]=='## HISTOGRAM':
                keys=lines[i+1].strip().split('\t')
                histogram = {keys[0]:[], keys[1]:[]}
        return histogram

    def strip_scilife_name(self, names): 
        N = {}
        P = {}
        preps = 'F_BCDE'
        for name_init in names:
            prep = 'A'
            name = name_init.replace('-', '_').replace(' ', '').split("_index")[0].split("_ss")[0].split("_dual")[0].strip()
            if name != '':
                while name[-1] in preps:
                    prep = name[-1] + prep
                    name = name[0: -1]
                if name != '':
                    N[name_init] = name
                    P[name_init] = prep.replace('_', '')
        return N, P
        
