#!/bin/bash -l
#SBATCH -A a2010002
#SBATCH -p node
#SBATCH -t 40:00:00
#SBATCH -J merge
#SBATCH -e merge.err
#SBATCH -o merge.out
#SBATCH --mail-user maya.brandi@scilifelab.se
#SBATCH --mail-type=ALL
module load bioinfo-tools
module load samtools
WP=/bubo/home/h24/mayabr/glob/RNA_analysis
while getopts ":p:g:m:c:" option; do
        case ${option} in
                p) path=${OPTARG};;
                g) gtf_file=${OPTARG};;
                m) mail=${OPTARG};;
		c) config_file=${OPTARG};;
        esac
done
shift $(( OPTIND - 1 ))

flowcells=''
for dir in "$@"; do
        flowcells=$flowcells" "$dir
done
flowcells=($flowcells)
nr_fc=${#flowcells[*]}
name_list=`for dir in ${flowcells[*]};do ls -d $dir/tophat_out_*|cut -f 2 -d '/'|sed 's/tophat_out_//g';done|sort|uniq`
echo ${flowcells[*]}

if [ $nr_fc  -lt 2 -o $# -lt 1 ]; then
echo "
Usage:
	merge.sh -p <path> -g [gtf file] -m [mail] -c [config file]  <run dir 1> <run dir 2> ... <run dir N>

Arguments:
	<path> 
	- Path to the intermediate directory containing the runs to be merged
       	- eg: /proj/a2010002/projects/m_muurinen_11_01a/intermediate, containing 
	the three runs: 20120127B_hiseq2000 20120228A_hiseq2000 20120328A_hiseq2000

	[mail] 
	- if you want to run cufflincs and HTseq on merged samples give your mail for SLURM
	messages

	[config_file]
	- if you want to run cufflincs and HTseq on merged samples give config_file
	specifying HTseq version

       	[gtf fie]
	- Optional! If not given, the script will only merge, but not run cufflinks 
	and HTseq on the merged data.
        - Reference annotation in gtf format, used by cufflinks and HTseq

	<run dir i>
	name of directories to merge, eg
	120727_AD14B4ACXX 120727_BD1591ACXX 120823_BD16GRACXX
	

Output:
	- A new directory for the merged data, called merged and placed in the intermediate directory.
	"
else
	# merge first two flowcells
	mkdir $path/merged
	for name in ${name_list[*]};do
		samp_dir=tophat_out_$name 
                bam_file_dupRem=accepted_hits_sorted_dupRemoved_${name}.bam
                bam_file=accepted_hits_${name}.bam
		if [[ -e ${flowcells[0]}/$samp_dir && -e ${flowcells[1]}/$samp_dir ]];then
			## merge read 
			echo 'merge'
			mkdir $path/merged/$samp_dir $path/merged/$samp_dir/logs
			samtools merge $path/merged/$samp_dir/$bam_file_dupRem $path/${flowcells[0]}/$samp_dir/$bam_file_dupRem $path/${flowcells[1]}/$samp_dir/$bam_file_dupRem
                        samtools merge $path/merged/$samp_dir/$bam_file $path/${flowcells[0]}/$samp_dir/$bam_file $path/${flowcells[1]}/$samp_dir/$bam_file
			## sum read counts
                        counts_0=(`grep 'reads have been filtered out' $path/${flowcells[0]}/$samp_dir/logs/prep_reads.log|cut -f 1,4 -d ' '`)
                        counts_1=(`grep 'reads have been filtered out' $path/${flowcells[1]}/$samp_dir/logs/prep_reads.log|cut -f 1,4 -d ' '`)
                        counts=$((${counts_0[1]}+${counts_1[1]}))
                        sorted=$((${counts_0[0]}+${counts_1[0]}))

			sed -e "s/\b${counts_0[1]}\b/$counts/g" $path/${flowcells[0]}/$samp_dir/logs/prep_reads.log|sed -e "s/\b${counts_0[0]}\b/$sorted/g" > $path/merged/$samp_dir/logs/prep_reads.log
		else if [ -e ${flowcells[0]}/tophat_out_$name ];then 
			echo 'cp'
			## copy reads and counts
                        cp -r $path/${flowcells[0]}/$samp_dir $path/merged
		else if [ -e ${flowcells[1]}/tophat_out_$name ];then 
                        ## copy reads and counts
                        cp -r $path/${flowcells[1]}/$samp_dir $path/merged
		fi;fi;fi 
	done

	# merge all remaining flowcells with "merged"
	if [ $nr_fc -gt 2 ]; then
	echo "merge all remaining flowcells with merged"
	for ((i=2; i<$nr_fc; i++));do
		for name in ${name_list[*]};do  
			echo ${flowcells[i]}               
			samp_dir=tophat_out_$name
			bam_file_dupRem=accepted_hits_sorted_dupRemoved_${name}.bam
			bam_file=accepted_hits_${name}.bam
			if [[ -e ${flowcells[i]}/$samp_dir && -e merged/$samp_dir ]];then
                                ## merge read
				echo 'merge'
				ls $path/merged/
                                samtools merge $path/merged/$samp_dir/temp_${bam_file_dupRem} $path/merged/$samp_dir/$bam_file_dupRem $path/${flowcells[$i]}/$samp_dir/$bam_file_dupRem
                                echo 'hej§'
				samtools merge $path/merged/$samp_dir/temp_${bam_file} $path/merged/$samp_dir/$bam_file $path/${flowcells[$i]}/$samp_dir/$bam_file
                                mv $path/merged/$samp_dir/temp_${bam_file_dupRem} $path/merged/$samp_dir/$bam_file_dupRem
                                mv $path/merged/$samp_dir/temp_${bam_file} $path/merged/$samp_dir/$bam_file
                                ## sum read counts
                                counts_0=(`grep 'reads have been filtered out' $path/merged/$samp_dir/logs/prep_reads.log|cut -f 1,4 -d ' '`)
                                counts_1=(`grep 'reads have been filtered out' $path/${flowcells[$i]}/$samp_dir/logs/prep_reads.log|cut -f 1,4 -d ' '`)
                                counts=$((${counts_0[1]}+${counts_1[1]}))
                                sorted=$((${counts_0[0]}+${counts_1[0]}))
                                sed -e "s/\b${counts_0[1]}\b/$counts/g" $path/merged/$samp_dir/logs/prep_reads.log|sed -e "s/\b${counts_0[0]}\b/$sorted/g" > $path/merged/$samp_dir/logs/temp_prep_reads.log
                                mv $path/merged/$samp_dir/logs/temp_prep_reads.log $path/merged/$samp_dir/logs/prep_reads.log
                	else if [ -e ${flowcells[i]}/$samp_dir ];then
                        	## copy reads and counts
				echo 'cp'
                        	cp -r $path/${flowcells[i]}/$samp_dir $path/merged
                	fi;fi
        	done
        done	
	fi
	if [ "$gtf_file" != "" ]; then
        ## get names of samples to be merged
	rerun=`for dir in ${flowcells[*]};do ls -d $dir/tophat_out_*|cut -f 2 -d '/'|sed 's/tophat_out_//g';done|sort|uniq -d`
	## run HTseq and cufflinks on meged samples
	for i in ${rerun[*]};do
		make_HT_cuff.py $i $gtf_file $mail $path/merged $config_file
        done
	fi
fi
