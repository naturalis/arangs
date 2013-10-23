#!/bin/bash

# usage string for use when arguments missing
usage=$0" project_name reference pair1 pair2 [orig_run_date]"

$project_name=$1
if [ -z $project_name ]
then
  echo $usage 1>&2
  exit 1
fi

# set up our root directories, which anchor all other paths
project_root="${HOME}/${project_name}"
if [ ! -d $project_root ]
then
  echo "${project_root} directory does not exist!"
  exit 1
fi

data_root="${project_root}/data"
if [ ! -d $data_root ]
then
  echo "${project_root} does not appear to be in the right structure, it needs a data directory" 1>&2
  exit 1
fi

fasta_root="${data_root}/fasta"
fastq_root="${data_root}/fastq"

result_root="${project_root}/results"
if [ ! -d $result_root ]
then
  echo "${project_root} does not appear to be in the right structure, it needs a results directory" 1>&2
  exit 1
fi

reference_name=$2
if [ -z $reference_name ]
then
  echo $usage 1>&2
  exit 1
fi

reference_path="${fasta_root}/${reference_name}"
if [ ! -f $reference_path ]
then
  echo "${ reference_name } not found in ${ fasta_root }" 1>&2
  exit 1
fi

pair1_name=$3
if [ -z $pair1_name ]
then
  echo $usage 1>&2
  exit 1
fi

pair1_path="${fastq_root}/${pair1_name}"
if [ ! -f $pair1_path ]
then
  echo "${ pair1_name } not found in ${ fastq_root }" 1>&2
  exit 1
fi

pair2_name=$4
if [ -z $pair2_name ]
then
  echo $usage 1>&2
  exit 1
fi

pair2_path="${fastq_root}/${pair2_name}"
if [ ! -f $pair2_path ]
then
  echo "${ pair2_name } not found in ${ fastq_root }" 1>&2
  exit 1
fi

date_stamp=`date "+%Y-%m-%d"`
this_result_dir="${result_root}/${date_stamp}"

orig_run_date=$5
if [ -n $orig_run_date ]
then
  this_result_dir="${result_root}/${orig_run_date}"
  if [ ! -d $this_result_dir ]
  then
    echo "${orig_run_date} does not appear to exist in ${result_root}" 1>&2
    exit 1
  fi
fi

sai1="${this_result_dir}/${pari1_name}.sai"
sai2="${this_result_dir}/${pari1_name}.sai"
sam="${this_result_dir}/all_reads.sam"

#mkdir -p will create a directory only if it does not already exist, otherwise it exits cleanly
#without -p it will exit with an error if the directory already exists
mkdir -p $this_result_dir

echo "Starting alignment" 1>&2

# create BWA index if not exists
if [ ! -e $reference_path.bwt ]; then
	echo "indexing ${reference_name}" 1>&2

	# Warning: "-a bwtsw" does not work for short genomes, 
	# while "-a is" and "-a div" do not work not for long 
	# genomes. Please choose "-a" according to the length 
	# of the genome.
	bwa index -a bwtsw $reference_path

        #always check the exit status and exit immediately if something goes wrong
	if [ 0 != $? ]
	then
	    echo "Problem encountered indexing ${reference_name}" 1>&2
	    exit 1
	fi
fi

#align pair1 against the reference
bwa aln -t 8 $reference_path $pair1_path > $
if [ 0 != $? ]
then
  echo "Problem encountered aligning ${pair1_name}" 1>&2
  exit 1
fi

#align pair2 against the reference
bwa aln -t 8 $reference_path $pair2_path > $sai2
if [ 0 != $? ]
then
  echo "Problem encountered aligning ${pair2_name}" 1>&2
  exit 1
fi

# use the sampe bwa command to produce sam output
bwa sampe $reference_path $sai1 $sai2 $pair1_path $pair2_path > $sam
if [ 0 != $? ]
then
  echo "Problem encountered producing sampe output" 1>&2
  exit 1
fi
echo "ALL COMPLETE" 1>&2
exit
