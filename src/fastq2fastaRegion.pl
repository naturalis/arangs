#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# process command line arguments. both are optional.
my ( $infile, $region );
GetOptions(
	'infile=s' => \$infile,
	'region=s' => \$region,
);

# create a handle to read from. this is either from an input file or from stdin;
my $handle;
if ( $infile ) {
	open $handle, '<', $infile;
}
else {
	$handle = \*STDIN;
}

# parse the region argument, if provided
my ( $seqid, $start, $stop );
if ( $region and $region =~ /^(.+?):(\d+)-(\d+)$/ ) {
	( $seqid, $start, $stop ) = ( $1, $2, $3 );
}

# start reading
my ( $seq, $id, $is_seq );
while(<$handle>) {
	chomp;
	if ( /^\@(.+)/ and not $is_seq ) {
		my $new_id = $1;
		
		# print the sequence
		if ( $id and $seq ) {
			if ( not $seqid or $id eq $seqid ) {
				print ">$id\n";
				print substr $seq, $start || 0, ( $stop - $start ) || length $seq;
			}
		}
		
		# start a new record
		$id       = $new_id;
		$is_seq   = 1;
		$seq      = '';
		next;
	}
	
	# reached the separator between sequence and phred, switch the is_seq flag
	$is_seq = 0 if /^\+$/ or /^\+$id$/;
	
	# grow the sequence
	$seq .= $_ if $is_seq;
}

# print the sequence
if ( $id and $seq ) {
	if ( not $seqid or $id eq $seqid ) {
		print $region ? ">$region\n" : ">$id\n";
		print substr $seq, $start || 0, ( $stop - $start ) || length $seq;
		print "\n";
	}
}