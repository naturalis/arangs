#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use lib '../lib';
use Fasterq;

# process command line arguments
my $freq     = 10000; # print progress every $freq reads
my $minqual  = 30;    # minimum phred score
my $p        = 0.01;  # sampling probability
my ( $in1, $in2 );    # mate pair files
GetOptions(
	'in1=s'     => \$in1,
	'in2=s'     => \$in2,
	'minqual=i' => \$minqual,
	'p=f'       => \$p,
);

# instantiate readers
my $fq1 = Fasterq->new($in1);
my $fq2 = Fasterq->new($in2);

# create out handles
my ( $outFH1, $outFH2 ) = ( make_handle( $in1 ), make_handle( $in2 ) );

my $counter = 1;
my %seen;
my %args = (
	'-minqual'     => $minqual,
	'-probability' => $p,
	'-encoding'    => Fasterq::ILLUMINA, # ASCII+64, unused 0, 1
	'-seen'        => \%seen,
);
while ( my $seq1 = $fq1->next_seq( %args ) ) {

	# fetch the corresponding mate
	my $name = $seq1->get_name;
	$name =~ s|\/1$|/2|;
	my $seq2 = $fq2->next_seq( '-name' => $name );
	
	# write sequences
	$seq1->to_string($outFH1);
	$seq2->to_string($outFH2);
	
	# report progress
	print "Done pair $counter\n" unless $counter++ % $freq;
}

sub make_handle {
	my $file = shift;
	$file =~ s/\.[^\.]+$/-sampled.fastq/;
	open my $fh, '>', $file or die "Can't open $file: $!";
	return $fh;
}
