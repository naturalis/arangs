#!/bin/bash
# this script extracts the seqid and start and stop locations of CDS features whose
# parent ID (an mRNA feature) is provided on the command line. Usage:
# mrna2regions.sh <gff3 file> <mrna ID>

GFF3=$1
MRNA=$2

# check if it's a gzipped file, adjust CAT accordingly
CAT=cat
if [[ `echo $GFF3 | grep -e '\.gz$'` != '' ]]; then
	CAT=zcat
fi

# parse out the CDS features, translate to region
CDSs=`$CAT $GFF3 | grep Parent=$MRNA | grep CDS | cut -f 9 | cut -f 1 -d ';' | sed -e 's/ID=//'`
for cds in $CDSs; do

	# column 1 is seq ID
	chromo=`$CAT $GFF3 | grep ID=$cds | cut -f 1`
	
	# column 4 is start location
	startloc=`$CAT $GFF3 | grep ID=$cds | cut -f 4`
	
	# column 5 is stop location
	stoploc=`$CAT $GFF3 | grep ID=$cds | cut -f 5`
	
	# print result
	echo "${chromo}:${startloc}-${stoploc}"
done