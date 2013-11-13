import numpy as np

import matplotlib.pyplot as plt
import sys
import os

NAME = sys.argv[1]
out_dir='Check_duplicates'
try:
	os.makedirs(out_dir)
	os.makedirs(out_dir + '/Counts_log')
	os.makedirs(out_dir + '/Fraction_dup')
	os.makedirs(out_dir + '/batch_Fraction_dup_of_all')
	os.makedirs(out_dir + '/Counts_per_1%batch')
	os.makedirs(out_dir + '/batch_Fraction_dup')
	os.makedirs(out_dir + '/Fraction_on_topX%')
except:
	pass
p = 20

withdup = open(str('tophat_out_'+NAME+'/cufflinks_out_'+NAME+'/genes.fpkm_tracking'), 'r')
duprem = open(str('tophat_out_'+NAME+'/cufflinks_out_dupRemoved_'+NAME+'/genes.fpkm_tracking'), 'r')

duprem_dict = {}
withdup_dict = {}

for line in withdup:
        l = line.split()
	key=l[0]
	val=l[9]
        if key[0]=='E':
                withdup_dict[key] = float(val)

for line in duprem:
        l = line.split()
        key=l[0]
        val=l[9]
        if key[0]=='E':
                duprem_dict[key] = float(val)

tupl = sorted([(value,key) for (key,value) in duprem_dict.items()],reverse=True)
sorted_genes_duprem = list( j for i,j in tupl)

NO_genes = len(duprem_dict)
total_withdup = sum(i for i in withdup_dict.values())
total_duprem = sum(i for i in duprem_dict.values())

### Top 3 ###
f=open(str(out_dir + '/top3.txt'),'a')
g1=sorted_genes_duprem[0]
g2=sorted_genes_duprem[1]
g3=sorted_genes_duprem[2]
d1=withdup_dict[g1]-duprem_dict[g1]
d2=withdup_dict[g2]-duprem_dict[g2]
d3=withdup_dict[g3]-duprem_dict[g3]
f1=round(float(d1)/float(total_withdup),4)
f2=round(float(d2)/float(total_withdup),4)
f3=round(float(d3)/float(total_withdup),4)
print >> f , """===========================================
Sample	{17}

Total fpkm				{16}
Total fpkm duplicates removed		{15}

Top 3 genes, (Duplicates remooved)
----------------------------------
        gene			fpkm(duprem)		fpkm(with dup)		diff		fract dupl/(total fpkm)
1.      {0}		{3}		{6}			{9}			{12}		
2.      {1}		{4}		{7}			{10}			{13}
3.      {2}		{5}		{8}			{11}			{14}
==========================================
""".format(g1,g2,g3,str(duprem_dict[g1]),str(duprem_dict[g2]),str(duprem_dict[g3]),str(withdup_dict[g1]),str(withdup_dict[g2]),str(withdup_dict[g3]),str(d1),str(d2),str(d3),str(f1),str(f2),str(f3),str(total_duprem),str(total_withdup),NAME)
#######################################################

topX_duprem = []
topX_withdup = []
group_duprem = []
group_withdup = []

N=0
for X in range(1,p + 1):  #np.arange(0,0.21,0.01)
	top_X_genes = sorted_genes_duprem[0:int(NO_genes*X/100.0)]
	topX_withdup.append(float(sum(withdup_dict[i] for i in top_X_genes)))
	topX_duprem.append(float(sum(duprem_dict[i] for i in top_X_genes)))
	group = sorted_genes_duprem[N:int(NO_genes*X/100.0)]
	group_duprem.append(sum(duprem_dict[i] for i in group))
	group_withdup.append(sum(withdup_dict[i] for i in group))
	N=int(NO_genes*X/100.0)


### 	plotting    ###
plt.rc('legend',**{'fontsize':11})

##      Plot fraction topX%
plt.figure(1)
#plt.subplot(221) 
plt.plot(range(1,p + 1) ,np.true_divide(topX_duprem, total_withdup) , label=r'(counts on topX dup rem)/(all counts with dup)')
plt.plot(range(1,p + 1) ,np.true_divide(topX_withdup, total_withdup) , label=r'(counts on topX with dup)/(all counts with dup)')
plt.xlabel(r'Top X% of genes, (sorted by counts after dup rem)')
plt.ylabel(r'Fraction on topX%')
plt.title(NAME)
plt.legend(loc='lower right')
#plt.show()
plt.savefig(out_dir+'/Fraction_on_topX%/'+NAME+'_Fraction_on_topX%.pdf')

##      Plot counts per 1%-batch of genes sorted. Only top 20% shown (Not cumulative)
plt.figure(2)
#plt.subplot(222) 
plt.plot(range(1,p + 1),group_duprem,label=r'dup rem',)
plt.plot(range(1,p + 1),group_withdup,label=r'with dup')
plt.xlabel(r'Top X% of genes, (sorted by counts after dup rem)')
plt.ylabel(r'Counts per 1%-batch')
plt.title(NAME)
plt.legend(loc='upper right')
#plt.show()
plt.savefig(out_dir+'/Counts_per_1%batch/'+NAME+'_Counts_per_1%batch.pdf')

##	Plot fraction duplicates per 1%-batch of genes sorted. Only top 20% shown (Not cumulative) (/All)
plt.figure(3)
#plt.subplot(223)
alphab = range(1,p + 1)
frecuencies = np.true_divide(np.array(group_withdup) - np.array(group_duprem), total_withdup)
pos = np.arange(len(alphab))
width = 1.0     
ax = plt.axes()
ax.set_xticks(pos + (width / 2))
ax.set_xticklabels(alphab)
plt.bar(pos, frecuencies, width, color='r')
plt.xlabel(r'Top X% - top (X-1)% of genes, (sorted by counts after dup rem)')
plt.ylabel(r'(Duplicates per 1%-batch)/(total counts with dup)')
plt.title(NAME)
#plt.show()
plt.savefig(out_dir+'/batch_Fraction_dup_of_all/'+NAME+'_batch_Fraction_dup_of_all.pdf')

##1      Plot fraction duplicates per 1%-batch of genes sorted. Only top 20% shown (Not cumulative) (/All per 1%-batch)
plt.figure(4)
#plt.subplot(224)
alphab = range(1,p + 1)
frecuencies = np.true_divide(np.array(group_withdup) - np.array(group_duprem), group_withdup)
pos = np.arange(len(alphab))
width = 1.0     
ax = plt.axes()
ax.set_xticks(pos + (width / 2))
ax.set_xticklabels(alphab)
plt.bar(pos, frecuencies, width, color='r')
plt.xlabel(r'Top X% - top (X-1)% of genes, (sorted by counts after dup rem)')
plt.ylabel(r'(Duplicates per 1%-batch)/(total counts with dup per 1%-batch)')
plt.title(NAME)
#plt.show()
plt.savefig(out_dir+'/batch_Fraction_dup/'+NAME+'_batch_Fraction_dup.pdf')

P1_counts_duprem=[]
P1_counts_withdup=[]
P1_fract_dupoly=[]
for gene in sorted_genes_duprem[0:int(NO_genes*0.2)]:
        if withdup_dict[gene]>5:
                P1_counts_duprem.append(duprem_dict[gene])
                P1_counts_withdup.append(withdup_dict[gene])
                P1_fract_dupoly.append(float(withdup_dict[gene]-duprem_dict[gene])/withdup_dict[gene])


##	Plot fract dupl (no bach) - corespond to ##1 
plt.figure(5)
#plt.subplot(211)
plt.plot(range(len(P1_fract_dupoly)),P1_fract_dupoly)
plt.xlabel(r'Genes sorted by counts after dup rem, (only top 20% genes shown)')
plt.ylabel(r'dupl/(counts with dupl)')
plt.title(NAME)
#plt.show()
plt.savefig(out_dir+'/Fraction_dup/'+NAME+'_Fraction_dup.pdf')

##      Plot counts log shale
plt.figure(6)
#plt.subplot(212)
plt.plot(range(len(P1_counts_duprem)),np.log(P1_counts_duprem),label=r'dup rem',)
plt.plot(range(len(P1_counts_duprem)),np.log(P1_counts_withdup),label=r'with dup')
plt.xlabel(r'Genes sorted by counts after dup rem, (only top 20% genes shown)')
plt.ylabel(r'Counts log scale')
plt.title(r'P352_120B')
plt.legend(loc='upper right')
#plt.show()
plt.savefig(out_dir+'/Counts_log/'+NAME+'_Counts_log.pdf')

