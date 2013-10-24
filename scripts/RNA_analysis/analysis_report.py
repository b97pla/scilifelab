#!/usr/bin/env python
import os
import sys
import commands
from optparse import OptionParser
from operator import itemgetter
from numpy import *
from pylab import *
import json
import yaml
import glob
import re
from datetime import date
from mako.template import Template
from mako.lookup import TemplateLookup
from texttable import Texttable
from bcbio.pipeline.config_loader import load_config
from scilifelab.google.project_metadata import ProjectMetaData
#from scilifelab.report.rst import make_logo_table
import operator

def image(fp, width):
   res = ".. figure:: %s\n    :width: %s\n\n" % (fp, width)
   return res

def indent_texttable_for_rst(ttab, indent=4, add_spacing=True):
    """Texttable needs to be indented for rst.

    :param ttab: texttable object
    :param indent: indentation (should be 4 *spaces* for rst documents)
    :param add_spacing_row: add additional empty row below class directives

    :returns: reformatted texttable object as string
    """
    output = ttab.draw()
    new_output = []
    for row in output.split("\n"):
        new_output.append(" " * indent + row)
        if re.search('.. class::', row):
            new_row = [" " if x != "|" else x for x in row]
            new_output.append(" " * indent + "".join(new_row))
    return "\n".join(new_output)

def make_template(Map_Stat, FPKM, GBC, Read_Dist, rRNA_table, strandness_table, complexity):    
#    TEMPLATE=make_logo_table+"""\
    TEMPLATE="""\
<%
    tablenum=0
    figurenum=0
%>
.. image:: sll_logo.gif
   :align: right

=======================
RNA-seq analysis report
=======================

SciLifeLab Stockholm
--------------------

${date}
-----------------   

Project name
^^^^^^^^^^^^

${project_id} 

UPPNEX project id
^^^^^^^^^^^^^^^^^

${uppnex}

Access the results on UPPMAX
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Please refer to the `Sequencing FAQ`_ for detailed instructions on how
to log in and access your files. If you have problems to access your
data, please contact genomics_support@scilifelab.se. If you have
questions regarding UPPNEX, please contact support@uppmax.uu.se.


.. _Sequencing FAQ: http://www.scilifelab.se/archive/pdf/tmp/SciLifeLab_Sequencing_FAQ.pdf

Specification of RNA-seq analysis delivery
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The RNA-seq analysis is delivered to /proj/${uppnex}/INBOX/${project_id}/analysis. The delivery consists of:

- **Analysis Report (this document):** Contains reports about: Mapping statistics, Read distribution, Correlation between FPKM values of samples, Gene body covarage, Library complexity and Strand specificity.

- **BAM files:** Two alignment files are placed under the alignments folder. We deliver two files per sample, one with and one without duplicates.

- **FPKM files:** The quantification folder contains FPKM values calculated using the Cufflinks program using ENSEMBL annotation of genes and transcripts for each sample. These files also contain the upper and lower limits of the confidence interval for the FPKM estimate. FPKM values are the paired-end equivalent of RPKM (Reads Per Kilobase per Million mapped reads; the standard measure for gene expression in RNA-seq.)

- **fpkm_table.txt:** A single file containing all of the FPKM values per gene and sample (duplicates included). This can be opened in Excel or a regular text processing application.

- **fpkm_table_isoforms.txt:**  A single file containing all of the FPKM values per isoform and sample (duplicates included). This can be opened in Excel or a regular text processing application.

- **count_table.txt:** Counts per gene and sample (duplicates included). For analyzing differential expression of genes or transcripts, it may be useful to have the raw read counts (the number of sequences that map to each gene/transcript). These are calculated using the HTSeq software and are collected into the table count_table.txt.

.. raw:: pdf

    PageBreak


The RNA-seq pipeline step by step
---------------------------------------

Your samples were analysed as follows:

1. **Mapping:** Reads were mapped with ${mapping} to the ${species} genome assembly, build ${genombuild}. 

2. **Merging:** Bamfiles from samples run on different lanes were merged with samtools.

3. **Sorting and Marking duplicates:** Merged bam files were sorted and duplicates removed using ${dup_rem}.

4. **Counts for genes:** Gene counts were generated using ${read_count} on bam files with duplicates included.

5. **FPKMs for genes and transcripts:** FPKMs for genes and transcripts were generated using ${quantifyer} on bamfiles with duplicates included.

6. **Mapping statistics:** Mapping statistics were calculated from numbers obtained by running bam_stat.py (included in rseqc/${rseqc_version}) on bam files with and without duplicates.

7. **Readdistribution:** read_distribution.py (included in rseqc/${rseqc_version}) was run on bam files with duplicates.

8. **Genebodycovarage:** geneBody_coverage.py (included in rseqc/${rseqc_version}) was run on bam files with duplicates.

9. **Strandedness:** infer_experiment.py (included in rseqc/${rseqc_version}) was run on bam files with duplicates.  

10. **Library complexity:** The complexity plot is generated by running ${Preseq} on bam files with duplicates.

11. **Corelation heatmap:** The correlation heatmap is generated using the R-package pheatmap on the FPKM table for genes.

.. raw:: pdf

    PageBreak

Results
--------
"""
    
    if Map_Stat:
        TEMPLATE=TEMPLATE+"""
<%
    tablenum=tablenum+1
%>
Mapping statistics
^^^^^^^^^^^^^^^^^^^^^

.. table:: **Table ${tablenum}.**

${Mapping_statistics}
    
**Tot NO Reads:** If paired-end reads, the total number of reads indicates the total number of sequenced paired-end reads. Since a paired-end read is made up of two sequenced fragments (mates), the total number of sequenced 100-bp regions is twice the number shown in this column. If single-end reads, this column reflects the total numer of sequences.

**UniqMapped:** The number of fragments that are mapped relative to the total number of sequenced fragments. 
    
**UniqMapped DuplRem:** We remove duplicate reads i.e. reads that map to the same genomic location (if paired-end, duplicates are defined as the paired end reads where both mates map to the same genomic positions). If the sequenced reads/pairs will be used to call sequence variants, it is recommended to use duplicate-removed data to obtain unbiased variant frequency calls. It is not uncommon to observe high rates of read duplication rates in RNA-seq libraries and they may not necessarily be the result of libraries with low complexity and/or over-amplification of libraries. The duplication rates can be gene specific i.e. a small number of genes are being highly expressed and comprise most of the duplicated reads/pairs, hence reflecting the underlying RNA distribution of the sample being investigated. In the case of small number of highly abundant genes consuming majority of the sequenced reads, one might need to sequence the samples deeper to obtain a better distribution of the whole transcriptome. When performing differential expression analyses, it is not recommended to use duplicate-removed reads/pairs. This is due to the fact that removing duplicate reads may effect expression values of transcripts. Since different samples have inherently different RNA distributions with different levels of duplicate reads/pairs for different transcripts, there could be a risk of removing biologically relevant information hence leading to incorrect comparisons.

.. raw:: pdf

    PageBreak
"""

    if Read_Dist:
        TEMPLATE=TEMPLATE+"""
<%
    tablenum=tablenum+1
%>
Read distribution
^^^^^^^^^^^^^^^^^^^
Table ${tablenum} contain information about the extent to which sequences from each sample mapped to different structural parts of genes, like coding exons, untranslated regions, and transcription start sites. The actual number itself is less important than the relative values for the different kinds of regions. For a normal RNA-seq experiment you should have a higher value in the CDS Exon column than in the others, for example. Perhaps the most easily interpretable column is the mRNA column, which gives the percentage of sequences that mapped to ENSEMBL-annotated mRNA (including coding regions and UTRs). While this number is not completely accurate (because ENSEMBL doe not completely describe the transcriptome), it is a useful summary statistic which should be relatively high for an mRNA-seq experiment, typically above 80%. 

**CDS:** Coding sequence exons.

**UTR:** Untranslated exon region.

**TES:** Transcription end site Down 1kb.

**TSS:** Transcription start site Up 1kb. 

**Intron:** Intronic or intergenic regions.

**mRNA:** Percentage of sequences that mapped to ENSEMBL-annotated mRNA.

.. table:: **Table ${tablenum}.**

${Read_Distribution}


.. raw:: pdf

    PageBreak
"""
    if GBC:
        TEMPLATE=TEMPLATE+"""
<%
    figurenum=figurenum+1
%>
Gene body covarage
^^^^^^^^^^^^^^^^^^^^^
Read coverage over gene body. To check if reads coverage is uniform and if there is any 5'/3' bias. All transcripts are scaled to 100 nucleotides and the read number is then calculated as the number of reads covering each nucleotide position. Figure ${figurenum} shows the average gene body coverage from all samples.

${GBC}

    **Figure ${figurenum}.**

.. raw:: pdf

    PageBreak
"""

    if FPKM:
        TEMPLATE=TEMPLATE+"""
<%
    figurenum=figurenum+1
%>
FPKM heatmap
^^^^^^^^^^^^^^^^^
Figure ${figurenum} shows the (Pearson) correlation between FPKM values of samples. 
    

${FPKM_heatmap}
    
    **Figure ${figurenum}.**
    
.. raw:: pdf

    PageBreak
""" 
    if complexity:
        TEMPLATE=TEMPLATE+"""
<%
    figurenum=figurenum+1
%>
Library complexity 
^^^^^^^^^^^^^^^^^^^^
Library complexity is an important measure when assessing the quality of samples and libraries. It refers to the amount of unique molecules that are present in a sample/library. Complex libraries, i.e. those with high number of unique molecules can be explored by sequencing with high coverage, whereas those with low complexity are probably not worthy of sequencing further or deeper. Figure ${figurenum} shows the number of unique molecules detected as a function of number of reads sequenced. It is obtained by extrapolating from the actual sequencing run and can be regarded as a predictor of the complexity of sample/library and assess its potential for further sequencing. For further information, please refer to: http://www.nature.com/nmeth/journal/v10/n4/full/nmeth.2375.html)


${complexity_plot}

    **Figure ${figurenum}.**

.. raw:: pdf

    PageBreak
                """

    if rRNA_table:
        TEMPLATE=TEMPLATE+"""
<%
    tablenum=tablenum+1
%>
Quantification of rRNA present in the samples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Presented percentage rRNA per sample. The numbers are calculated as the number of reads mapping to Ensembl genes that have their <source> field set to "rRNA" in the gtf file, divided by the total number of mapped reads.

.. table:: **Table ${tablenum}.**

${rRNA_table}

.. raw:: pdf

    PageBreak
"""
    if strandness_table:
        TEMPLATE=TEMPLATE+"""
<%
    tablenum=tablenum+1
%>
Percentage of reads/pairs mapped to the expected strand
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                     
The script infer_experiment.py from RSeQC package is used to evaluate strand-specificity of the sequencing. Table ${tablenum}  shows the percentage of read pairs that were detected on the strand the gene is expressed.

.. table:: **Table ${tablenum}.**

${strandness_table}


.. raw:: pdf

    PageBreak
"""
    TEMPLATE=TEMPLATE+"""
Tools and references
---------------------------

**Reference genome assembly:**
${species}, ${genombuild}

http://www.ensembl.org

**Ensembl annotation:** Release ${anotation_version}

**Mapping:** 
${mapping}

http://tophat.cbcb.umd.edu/
    
**Duplicate removal:**
${dup_rem}

http://picard.sourceforge.net/command-line-overview.shtml
    
**Read count:**
${read_count}

http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html
    
**RPKM/FPKM values:**
${quantifyer}

http://cufflinks.cbcb.umd.edu/

**Gene body coverage, Read distribution and Strandspecifisity:**
rseqc/${rseqc_version}

http://rseqc.sourceforge.net/


**Library Complexity:**    
${Preseq}

http://smithlab.usc.edu/plone/software/librarycomplexity

.. raw:: pdf

    PageBreak
"""
    return TEMPLATE


def generate_report(proj_conf,single_end,stranded):

    d = {
        'project_id': proj_conf['id'],
        'samplenames': ' '.join(proj_conf['samples']),
        'latex_opt' : "",
        'uppnex': "",
        'mapping':"",
        'dup_rem':"",
        'read_count':"",
        'quantifyer':"",
        'gene_body_cov':"",
        'FPKM_heatmap':"",
        'Mapping_statistics': "",
        'Read_Distribution':"",
        'rRNA_table':"",
        'GBC':"",
        'strandness_table':"",
        'complexity_plot': "",
        'species':"",
        'genombuild':"",
        'rseqc_version':'',
        'Preseq':'',
        'date':date.today(),
        'anotation_version':''
        }

    if stranded:
        analysis_type='RNA-seq stranded analysis'
    #the name in custom_algorithms section in opt/config/post_process.yaml
    else:
        analysis_type='RNA-seq analysis'

    ## Latex option (no of floats per page)
    floats_per_page = '.. raw:: latex\n\n   \setcounter{totalnumber}{8}'
    d['latex_opt'] = floats_per_page


    ## Metadata fetched from the 'Genomics project list' on Google Docs 
    try:
        proj_data = ProjectMetaData(proj_conf['id'], proj_conf['config'])
        uppnex_proj = proj_data.uppnex_id
        if proj_data.ref_genome == "hg19":
            d['species'] = 'Human'
        elif proj_data.ref_genome == "mm9":
            d['species'] = 'Mouse'
        else:
            d['species'] = proj_data.ref_genome
    except:
        uppnex_proj = "b201YYXX"
        print "No uppnex ID fetched"
        pass
    d['uppnex'] = uppnex_proj 


     ## RNA-seq tools fetched from config file post_process.yaml
    try:
        tools = proj_conf['config']['custom_algorithms'][analysis_type]
        d['mapping'] = os.path.join(tools['aligner'],tools['aligner_version'])
        d['dup_rem'] = os.path.join(tools['dup_remover'],tools['dup_remover_version'])
        d['read_count'] = os.path.join(tools['counts'],tools['counts_version'])
        d['quantifyer'] = os.path.join(tools['quantifyer'],tools['quantifyer_version'])
        d['genombuild'] = tools[proj_data.ref_genome]['name']
        d['rseqc_version'] = tools['rseqc_version']
        d['Preseq'] = tools['preseq']
        d['anotation_version'] = tools[proj_data.ref_genome]['annotation_release']
    except:
        print "Could not fetched RNA-seq tools from config file post_process.yaml"
        pass


    ## Mapping Statistics
    tab = Texttable()
    tab.set_cols_dtype(['t','t','t','t'])
    tab.header(['Sample','Tot NO Reads','UniqMapped','UniqMapped DuplRem'])
    statistics={}
    try:
        for sample_name in proj_conf['samples']:
            try:
                f = open('tophat_out_'+sample_name+'/logs/prep_reads.log', 'r')
                tot_NO_read_pairs = f.readlines()[2].split()[3]
                f.close()
                f = open('tophat_out_'+sample_name+'/stat'+sample_name, 'r')
                dict = make_stat(f,tot_NO_read_pairs,single_end)
                tab.add_row([sample_name,tot_NO_read_pairs,str(dict['bef_dup_rem']['%uniq_mapped'])+'%',str(dict['aft_dup_rem']['%uniq_mapped'])+'%'])
                statistics[sample_name] = dict
            except:
                print 'Could not make mapping statistics for sample '+sample_name

        d['Mapping_statistics'] = indent_texttable_for_rst(tab)
        stat_json = open('stat.json','w')
        print>> stat_json, statistics
        stat_json.close()
    except:
        print "Could not make Mapping Statistics table"
        pass

    ## Read Distribution 
    try:
        tab = Texttable()
        tab.set_cols_dtype(['t','t','t','t','t','t','t','t'])
        print dir(tab)
        tab.header(["Sample","CDS","5'UTR","3'UTR","Intron","TSS","TES","mRNA"])
        read_dist = {}
        for i in range(len(proj_conf['samples'])):
            sample_name = proj_conf['samples'][i]
            dict = {}
            try:
                f = open('RSeQC_rd_'+sample_name+'.out','r')
                dict = read_RSeQC_rd233(f)
                row = [sample_name,dict['CDS_Exons']['Tags/Kb'], dict["5'UTR_Exons"]['Tags/Kb'],
                dict["3'UTR_Exons"]['Tags/Kb'],dict['Introns']['Tags/Kb'], 
                dict['TSS_up_1kb']['Tags/Kb'],dict['TES_down_1kb']['Tags/Kb'], dict['mRNA_frac']]
                tab.add_row(row)
                read_dist[sample_name] = dict
            except:
                print "Could not make read distribution for sample "+sample_name
                pass
        RSeQC_rd_json = open('RSeQC_rd.json','w')
        print >> RSeQC_rd_json, read_dist
        RSeQC_rd_json.close()
        d['Read_Distribution'] = indent_texttable_for_rst(tab)
    except:
        print "Could not make Read Distribution table"
        pass

    ## Gene Body Coverage
    try:
        figure()
        x = range(0,101)
        for i in range(len(proj_conf['samples'])):
            y = zeros(101)
            sample_name = proj_conf['samples'][i]
            f = open(sample_name + '.geneBodyCoverage.txt','r')
            for line in f.readlines():
                try:
                    key = int(line.split()[0])
                    val = int(line.split()[1])
                    y[key] = val
                except:
                    pass
            plot(x,y)#,label=proj_conf['samples'][i])
        #legend(loc='upper left',fontsize='xx-small')
        ylabel("read number")
        xlabel("percentile of gene body (5'->3')")
        savefig('gbc.pdf')
        d['GBC'] = image("gbc.pdf", width="100%")
    except:
        print "could not make GBC plot"

    ##  FPKM_heatmap
    if os.path.exists("FPKM_heatmap.pdf"):
        d['FPKM_heatmap'] = image("FPKM_heatmap.pdf", width="100%")
    else:
        print "could not make FPKM heatmap"

    ## complexity plot
    if os.path.exists("Rplots.pdf"):
        d['complexity_plot'] = image("Rplots.pdf", width="100%")
    else:
        complexity=False
        print "could not make complexity plot"

    ## rRNA_table
    try:
        tab = Texttable()
        tab.set_cols_dtype(['t','t'])
        tab.header(["Sample","rRNA"])
        f=open('rRNA.quantification','r')
        D={}
        for line in f:
            D[str(line.split('\t')[0].strip())]=str(line.split('\t')[1].strip())
        for sample_name in proj_conf['samples']:
            if D.has_key(sample_name):
                tab.add_row([sample_name,D[sample_name]])
        d['rRNA_table']=indent_texttable_for_rst(tab)
        f.close()

    except:
        print "could not generate rRNA table"
        pass   
 
    ## strandness_table
    try:
        tab = Texttable()
        tab.set_cols_dtype(['t','t'])
        tab.header(["Sample","strand-specific reads"])
        try:
            f=open('infer_experiment.json', 'rb')
            data=json.load(f)        
        except:
            print "can't open infer_experiment.json\n"
        D=data
        for sample_name in proj_conf['samples']:
            if D.has_key(sample_name):
                tab.add_row([sample_name,str(float(D[sample_name])*100)+'%'])
                print str(float(D[sample_name])*100)

        d['strandness_table']=indent_texttable_for_rst(tab)
        f.close()
    except:
        print "could not generate strandness_table"
        pass

    return d

def make_stat(f,counts,single_end):
    aft_dup_rem={}
    bef_dup_rem={}
    stat={'Total_No_reads':counts}
    step=0
    file_content = f.readlines()
    for i, line in enumerate(file_content):
        if line.strip() != '':
            if line.split()[0]=='bam_stat.py':
                stat['version']=line.strip()
            if line[0]=='#':
                step = step+1
            elif step==1:
                try:
                    line = line.split(':')
                    bef_dup_rem[line[0].strip()] = line[1].strip()
                except:
                    pass
            elif step==2:
                try:
                    line = line.split(':')
                    aft_dup_rem[line[0].strip()] = line[1].strip()
                except:
                    pass
    f.close()
    if single_end:
        bef_dup_rem['%uniq_mapped'] = round(100*(float(bef_dup_rem['Uniquely mapped']))/(float(counts)),2)
        aft_dup_rem['%uniq_mapped'] = round(100*(float(aft_dup_rem['Uniquely mapped']))/(float(counts)),2)
        aft_dup_rem['%spliced'] = round(100*float(aft_dup_rem['spliced'])/float(aft_dup_rem['Uniquely mapped']))
    else:
        if float(counts) > 0:
            bef_dup_rem['%uniq_mapped'] = round(100*(float(bef_dup_rem['Read-1'])+float(bef_dup_rem['Read-2']))/(2*float(counts)),2)
            aft_dup_rem['%uniq_mapped'] = round(100*(float(aft_dup_rem['Read-1'])+float(aft_dup_rem['Read-2']))/(2*float(counts)),2)
        else:
            aft_dup_rem['%uniq_mapped'] = 0
            bef_dup_rem['%uniq_mapped'] = 0
        if (float(aft_dup_rem['Read-1']) + float(aft_dup_rem['Read-2'])) > 0:
            aft_dup_rem['%spliced'] = round(100*float(aft_dup_rem['spliced'])/(float(aft_dup_rem['Read-1'])+float(aft_dup_rem['Read-2'])))
        else:
            aft_dup_rem['%spliced'] = 0
    stat['bef_dup_rem']=bef_dup_rem
    stat['aft_dup_rem']=aft_dup_rem
    return stat

def read_RSeQC_rd233(f):
    dict={}
    delim = 0
    for i, l in enumerate(f):
        if l[0] != '=':
            line = l.split()
            if delim == 1:
                dict[line[0]] = {'Total_bases':line[1],'Tag_count':line[2],'Tags/Kb':line[3]}
            elif line[0] == "read_distribution.py":
                dict['version'] = l.strip()
            elif line[1] == "Tags":
                Total_Tags = int(line[2])
        else:
            delim = delim + 1
    f.close()
    mRNA_frac = (float(dict['CDS_Exons']['Tag_count']) + float(dict["5'UTR_Exons"]['Tag_count'])
                        + float(dict["3'UTR_Exons"]['Tag_count'])) / Total_Tags
    dict['mRNA_frac'] = str(round(mRNA_frac*100,2))+'%'
    return dict


def main(project_id,sample_names,single_end,config_file,Map_Stat,Read_Dist,FPKM,rRNA_table,GBC,stranded, strandness_table,complexity):
    if not sample_names:
        sample_names = commands.getoutput("ls -d tophat_out_*|sed 's/tophat_out_//g'").split('\n')
    else:
        sample_names = sample_names.split(',')
    TEMPLATE = make_template(Map_Stat, FPKM, GBC, Read_Dist, rRNA_table, strandness_table, complexity)
    if config_file:
        config = load_config(config_file)
    else:
        config = {}

    projectfile = "%s.mako" % (project_id)
    fp = open(projectfile, "w")
    fp.write(TEMPLATE)
    fp.close()
    tmpl = Template(filename=projectfile)

    proj_conf = {
        'id' : project_id,
        'config' : config,
        'samples': sample_names
         } 
    d = generate_report(proj_conf,single_end,stranded)
    rstfile = "%s.rst" % (project_id)
    fp = open(rstfile, "w")
    fp.write(tmpl.render(**d))
    fp.close()


if __name__ == "__main__":
    parser = OptionParser(usage = "analysis_reports.py <project id> [Options]")
    parser.add_option("-e", "--single_end", dest="single_end", action="store_true",default=False,
    help = "For single end runs")
    parser.add_option("-n", "--dry_run", dest="dry_run", action="store_true",default=False)
    parser.add_option("-c", "--config-file", dest="config_file", default=None,
    help = "FILE should be a config file. (post_process.yaml)")
    parser.add_option("-r", "--rRNA_table", dest="rRNA_table",action="store_true", default=False,
    help = "to include an rRNA table in the report. Requires rRNA.quantification from quantify_rRNA.sh")
    parser.add_option("-s", "--Map_Stat", dest="Map_Stat",action="store_true", default=False,
    help = "to include mapping statistics table in the report. Requires stat from get_stat.sh")
    parser.add_option("-d", "--Read_Dist", dest="Read_Dist",action="store_true", default=False,
    help = "to include read distribution table in the report. Requires stderr output files from read_distribution.py")
    parser.add_option("-f", "--FPKM", dest="FPKM",action="store_true", default=False,
    help = "to include fpkm-heatmap. Requires FPKM_heatmap.pdf from correl.R")
    parser.add_option("-g", "--gbc", dest="GBC",action="store_true", default=False,
    help = "to include gene body coverage plot in the report.")
    parser.add_option( "-w", "--stranded", dest="stranded",action="store_true", default=False,
    help = "to include strandness table with % of strand-specifically mapped reads")
    parser.add_option( "-a", "--sample_names", dest="sample_names", default=None,
    help = "Samplenames should be given as a coma delimited string. Default will be the samples given by the tophat_out - directories")
    parser.add_option( "-b", "--complexity", dest="complexity",action="store_true", default=False,
    help = "to include library complexity plot")
    (options, args) = parser.parse_args()
    if len(args) < 1:
        print __doc__
        sys.exit()
    kwargs = dict(sample_names = options.sample_names,
        single_end = options.single_end,
        config_file = options.config_file,
        Map_Stat = options.Map_Stat,
        Read_Dist = options.Read_Dist,
        FPKM = options.FPKM,
        rRNA_table = options.rRNA_table,
        GBC = options.GBC,
        stranded = options.stranded,
        strandness_table = options.stranded,
        complexity = options.complexity) 
    main(*args, **kwargs)
