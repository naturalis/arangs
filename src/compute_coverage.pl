#!/usr/bin/perl
use strict;
use warnings;
use Bio::DB::Sam;

# these are our assumptions from the assignment.
# it is more flexible of course if we can specify
# these on the command line.
my $chromosome = 'SL2.40ch00';
my $bamfile    = 'U0015717.bam';
my $start      = 1;
my $end        = 100_000;

# here we create a Bio::DB::Sam object, which 
# presents the contents of a BAM file as a 
# database which we can navigate and query.
my $sam = Bio::DB::Sam->new( -bam => $bamfile );

# here we query the sam object to get the aligned
# reads within our genomic region
my @reads = $sam->get_features_by_location(
	-seq_id => $chromosome,
	-start  => $start,
	-end    => $end,
);

# here we iterate over all the reads within our
# genomic region and sum their lengths
my $total_reads_length;
for my $read ( @reads ) {
	$total_reads_length += $read->length;
} 

# the coverage within the region is the total 
# lengths of the reads divided by the length
# of the region
my $coverage = $total_reads_length / ( $end - $start );
print "coverage on $chromosome between $start and $end is $coverage\n";


