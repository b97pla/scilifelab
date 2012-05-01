
import os
import sys

PREFIX_LENGTH = 1

infiles = sys.argv[1:]

to_merge = {}
for i in range(len(infiles)):
        if len(infiles[i]) == 0: continue 
	name_suffix = os.path.basename(infiles[i])[PREFIX_LENGTH:]
	to_merge[name_suffix] = [infiles[i]]
	for j in range(i+1,len(infiles)):
		if os.path.basename(infiles[j])[PREFIX_LENGTH:] != name_suffix:
			continue
		to_merge[name_suffix].append(infiles[j])
		infiles[j] = ""

for suffix, files in to_merge.items():
	if len(files) < 2: continue
	prefix = "-".join([os.path.basename(f)[0:PREFIX_LENGTH] for f in files])
	concatfile = "%s%s" % (prefix,suffix)
	with open(concatfile,"w") as outh:
		for f in files:
			with open(f,"r") as inh:
				for row in inh:
					outh.write(row)
