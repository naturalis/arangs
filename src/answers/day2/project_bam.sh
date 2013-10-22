#!/bin/bash

# set up our root directories, which anchor all other paths
project_root="${HOME}/arangs"
result_root="${project_root}/results"

# usage string for use when arguments missing
usage=$0" sam_file_name [original_run_date]"

sam_name=$1
if [ -z $sam_name ]
then
  echo $usage 1>&2
  exit 1
fi

date_stamp=`date "+%Y-%m-%d"`
this_result_dir="${result_root}/${date_stamp}"

orig_run_date=$2
if [ -n $orig_run_date ]
then
  this_result_dir="${result_root}/${orig_run_date}"
  if [ ! -d $this_result_dir ]
  then
    echo "${orig_run_date} does not appear to exist in ${result_root}" 1>&2
    exit 1
  fi
fi
    
sam_path="${this_result_dir}/${sam_name}"
if [ ! -f $sam_path ]
then
  echo "${ sam_name } not found in ${ this_result_dir }" 1>&2
  exit 1
fi

#samtools sort will kindly add .bam to each of these for us
unaligned_bam="${this_result_dir}/unaligned_reads"
aligned_bam="${this_result_dir}/aligned_reads"

echo "producing aligned_reads.bam" 1>&2
# -bS   = input is SAM, output is BAM
# -F 4  = remove unmapped reads
samtools view -bS -F 4 $sam_path | samtools sort - $aligned_bam
if [ 0 != $? ]
then
  echo "Problem encountered producing aligned_reads.bam" 1>&2
  exit 1
fi

echo "producing unaligned_reads.bam" 1>&2
#use the perl script to get the unaligned reads
perl ${project_root}/src/filter_sam.pl -u $sam_path | samtools -bS - | samtools sort - $unaligned_bam
if [ 0 != $? ]
then
  echo "Problem encountered producing unaligned_bam" 1>&2
  exit 1
fi

echo "ALL COMPLETE" 1>&2
exit
