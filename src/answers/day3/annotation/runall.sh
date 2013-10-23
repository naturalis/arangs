#!/bin/bash

# variables for this script
DATA=data
SRC=src
DATABASE=$DATA/ITAG2_3_genomic.fasta

# command line argument 1
query=$1

# get the genomic coordinates
REGION=`bash $SRC/ketchup.sh $query $DATABASE`

# parse the coordinates
CHROMO=`echo $REGION | cut -f 1 -d ':'`
START=`echo $REGION | cut -f 1 -d '-' | sed -e 's/.*://'`
STOP=`echo $REGION | cut -f 2 -d '-'`

# get the gene ID
ID=`python $SRC/chromo_coords.py $DATA/ITAG2.3_gene_models.gff3 $START $STOP $CHROMO | cut -f 2 -d ' '`

# get the consensus sequence
sh $SRC/ex_day3.sh $CHROMO:$STOP-$START $DATA/U0015717.bam
 


