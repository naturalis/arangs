#!/bin/bash

# all shell scripts (and perl, ruby, python, etc) get $0, which
# stores the quoted path to the executable being run.
# use $0 with a list of arguments required.  You should change this to
# fit your scripts argument needs.  Change arg1 to be something more
# explanatory, e.g. path_to_reference_file
usage=$0" coordinates bam_file"


# shell scripts give you commandline arguments in numeric order from 1 onward
arg1=$1

# -z tests for zero bytes in a string
if [ -z $arg1 ]
then
  # print usage to STDERR
  echo $usage 1>&2
  echo "example: sh ex_day3.sh SL2.40ch01:1-1000 /gsd/To_Participant/ARANGS13/U0015717.bam "
  exit 1
fi

# you may only need one argument, or you may need more than 2.
# remove the next 7 lines, or copy them and increment the number
# to fit your scripts needs
arg2=$2
if [ -z $arg2 ]
then
  # print usage to STDERR
  echo $usage 1>&2
  exit 1
fi

#bai=`echo $arg2 | sed -e "s/bam/bai/"`


# do bwa sampe if needed
if [ ! -e $arg2.bai ]; then

	# create paired-end SAM file
	echo "you need to do the indexing - indexed bam file missing"
	#bwa sampe $REFERENCE $SAIS $FASTQS > $SAM
else
	echo "indexed bam file present"
fi
echo $arg1 $arg2
samtools mpileup -u -r $arg1 $arg2 | bcftools view -cg - | /usr/share/samtools/vcfutils.pl vcf2fq > consensus.fastq 
perl src/fastq2fastaRegion.pl --infile consensus.fastq --region $arg1 > consensus.fasta
#python remove_ends.py consensus.fasta > final.fasta
python src/remove_ends.py consensus.fasta



# do things with arg1 (and arg2, etc)
# for each command that you run, check its exit status to make sure it ran ok
if [ 0 != $? ]
then
  echo "Problem encountered" 1>&2
  exit 1
fi

# exit with a zero exit status, to signal shell, or wrapper scripts that it ran to completion
exit
