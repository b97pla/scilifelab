#!/bin/sh 
# Poor man's rsync copy from sequencer raw dumps to pre-processing machines

rsync -craz -e ssh --exclude "*_qseq.txt" --exclude "*.tmp" --exclude .AppleDouble --exclude .DS_Store \
	--progress hiseq@comicbookguy.scilifelab.se:/home/hiseq.hiseq/illumina/ /mnt/illumina

rsync -craz -e ssh --delete --exclude .AppleDouble --exclude "*.old" --exclude .DS_Store \
	--progress hiseq@comicbookguy.scilifelab.se:/home/hiseq.hiseq/samplesheets_illumina /mnt
