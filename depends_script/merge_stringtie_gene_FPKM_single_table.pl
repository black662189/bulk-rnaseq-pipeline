#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

my $usage = "usage: $0 1/1_abundance.txt 2/2_abundance.txt 3/3_abundance.txt ...\n\n";

unless (@ARGV) {
	die $usage;
}

my @files = @ARGV;

if (scalar @files == 1) {
	if (-s $files[0]) {
		# allow for a file listing the various files.
		@files = `cat $files[0]`;
		chomp @files;
	} else {
		die $usage;
	}
}

=header_format

0	Gene ID
1	Gene Name	
2	Reference
3	Strand
4	Start	
5	End
6	Coverage	
7	FPKM
8	TPM

=cut

main: {
	my %data;
	
	foreach my $file (@files) {	
		open (my $fh, $file) or die "Error, cannot open file $file";
		my $header = <$fh>; # ignore it
		while (<$fh>) {
			chomp;
			my @x = split(/\t/);
			my $acc = $x[0];
			my $count = $x[7];
			$data{$acc}->{$file} = $count;
		}
		close $fh;
	}
	
	my @filenames = @files;
	foreach my $file (@filenames) {
		$file = dirname($file);
	}

	print join("\t", "gene_id", @filenames) . "\n";
	foreach my $acc (keys %data) {
		print "$acc";

		foreach my $file (@files) {
		    my $count = $data{$acc}->{$file};
		    unless (defined $count) {
		        $count = "NA";
		    }
		    print "\t$count";
		}
		
		print "\n";
	}

	exit(0);
}
