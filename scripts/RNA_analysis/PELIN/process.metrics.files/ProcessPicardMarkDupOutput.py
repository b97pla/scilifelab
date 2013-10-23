#!/usr/bin/python

import sys
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument list:', str(sys.argv)

#Open the dupmetrics file
filename = sys.argv[1] + '_alignments/' + sys.argv[1] + '.dupmetrics.txt'

infile = open(filename,'r');
outfile = open('dupmetrics.txt','a');

for line in infile:
    columns = line.split('\t')
    if columns[0] == 'Unknown Library':
	optdup = float(columns[6])/float(columns[2])
	outfile.write(sys.argv[1])
	outfile.write('\t')
	outfile.write(columns[2])
	outfile.write('\t')
	outfile.write(str(optdup))
	outfile.write('\t')
	outfile.write(columns[7])
	outfile.write('\n')

outfile.close()   
infile.close()
