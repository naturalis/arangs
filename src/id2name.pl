#!/usr/bin/perl
use strict;
use LWP::UserAgent;

# this script accepts a gene ID as per the ITAG2.3 naming conventions
# in their GFF3 and will produce the associated external gene name
# (if it exists) by accessing BioMart. Gene ID's should be formatted
# such as 'Solyc01g005000.2'
my $gene_id = shift || die "Usage: $0 <gene ID>\n";

# this is a template that, given an EnsEMBL gene ID, returns the
# gene name, if any
my $xml = <<"TEMPLATE";
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  
	virtualSchemaName="default2" 
	formatter="CSV" 
	header="0" 
	uniqueRows="0" 
	count="" 
	datasetConfigVersion="0.7">
	<Dataset 
		name="slycopersicum_eg_gene" interface="default">
		<Filter name="ensembl_gene_id" value="$gene_id"/>
		<Attribute name="external_gene_id" />
	</Dataset>
</Query>
TEMPLATE

# base URL for the biomart service
my $path="http://www.biomart.org/biomart/martservice?";

# create POST request object that includes the XML
my $request = HTTP::Request->new("POST",$path,HTTP::Headers->new(),'query='.$xml."\n");

# create client object that will execute the request
my $ua = LWP::UserAgent->new;

# do the request
$ua->request(
	$request, 
	sub{   
		my ( $data, $response ) = @_;
		if ($response->is_success) {
			print "$data";
		 }
		 else {
		     warn ("Problems with the web server: ".$response->status_line);
		 }
	},
	1000
);
