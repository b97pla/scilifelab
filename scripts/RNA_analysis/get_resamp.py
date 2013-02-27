import sys

if len(sys.argv) ==1:
        print """
To be run from the analysis directory: eg /proj/a2012043/private/nobackup/projects/e_lundberg_11_02/intermediate/merged
It takes as argument <M reads ordered>.  Find this number in genomics project list.

Usage:

get_resamp.py <M reads ordered> 	

"""
	sys.exit()

ordered=sys.argv[1]

n_passed=0
n_notpassed=0

dict = eval(open('stat.json').read())
samps=dict.keys()
print 'sample    M_reads       M_reads_(dup_rem)         Status'
for sample in samps:
	Total_M_read=float(dict[sample]['Total_No_read-pairs'])/1000000
	R1=dict[sample]['aft_dup_rem']['Read-1']
	R2=dict[sample]['aft_dup_rem']['Read-2']
	M_reads_aft_dup_rem=(float(R2)+float(R1))/2000000
	if M_reads_aft_dup_rem<float(ordered)/2:
		status='NP'
		n_notpassed=n_notpassed+1
	else:
		status='P'
		n_passed=n_passed+1
	print sample+'	'+str(round(Total_M_read,2))+'	'+str(round(M_reads_aft_dup_rem,2))+'	'+status

print 'Number of P samples:	'+ str(n_passed)
print 'Number of NP samples:	'+ str(n_notpassed)


