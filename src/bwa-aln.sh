#!/bin/bash

# this is a simple shell script that goes through the steps of a reference-guided
# assembly of paired-end illumina data

# this is the reference genome in FASTA format. Typically there would be as many
# sequences in this as there are chromosomes
REFERENCE=reference.fasta

# these are the mate pairs in FASTQ format as they are produced by the illumina
# pipeline
FASTQS="foo_R1_001.fastq foo_R2_001.fastq"

# threads for BWA align
CORES=4

# recreate BWA index if not exists
if [ ! -e $REFERENCE.bwt ]; then
	echo "going to index $REFERENCE"

	# Warning: "-a bwtsw" does not work for short genomes, 
	# while "-a is" and "-a div" do not work not for long 
	# genomes. Please choose "-a" according to the length 
	# of the genome.
	bwa index -a bwtsw $REFERENCE
else 
	echo "$REFERENCE already indexed"
fi

# this will hold the two *.sai files that are produced
SAIS=""

# this will hold the on SAM file that will be produced
SAM=""
for FASTQ in $FASTQS; do

	# create new names
	OUTFILE=$FASTQ-$REFERENCE.sai
	SAIS="$SAIS $OUTFILE"
	SAM=`echo $OUTFILE | sed -e "s/_R.*/-$REFERENCE.sam/"`

	# note: we don't do basic QC here, because that might mean
	# that the mate pairs in the FASTQ files go out of order,
	# which will result in the bwa sampe step taking an inordinate
	# amount of time

	# do bwa aln if needed
	if [ ! -e $OUTFILE ]; then
		echo "going to align $FASTQ against $REFERENCE"

		# use $CORES threads
		bwa aln -t $CORES $REFERENCE $FASTQ > $OUTFILE
	else
		echo "alignment $OUTFILE already created"
	fi
done

# do bwa sampe if needed
if [ ! -e $SAM ]; then

	# create paired-end SAM file
	echo "going to run bwa sampe $FASTA $SAIS $FASTQS > $SAM"
	bwa sampe $REFERENCE $SAIS $FASTQS > $SAM
else
	echo "sam file $SAM already created"
fi		

# do samtools filter if needed
if [ ! -e $SAM.filtered ]; then
	# -bS   = input is SAM, output is BAM
	# -F 4  = remove unmapped reads
	# -q 50 = remove reads with mapping qual < 50
	echo "going to run samtools view -bS -F 4 -q 50 $SAM > $SAM.filtered"
	samtools view -bS -F 4 -q 50 $SAM > $SAM.filtered
else
	echo "sam file $SAM.filtered already created"
fi

# do samtools sorting if needed
if [ ! -e $SAM.sorted.bam ]; then

	# sorting is needed for indexing
	echo "going to run samtools sort $SAM.filtered $SAM.sorted"
	samtools sort $SAM.filtered $SAM.sorted
else
	echo "sam file $SAM.sorted already created"
fi

# created index for BAM file if needed
if [ ! -e $SAM.sorted.bam.bai ]; then

	# this should result in faster processing
	echo "going to run samtools index $SAM.sorted.bam"
	samtools index $SAM.sorted.bam
else
	echo "BAM file index $SAM.sorted.bam.bai already created"
fi
