#!/usr/bin/perl -w 

if (@ARGV < 2) {
	print STDERR "<CNV input file>  <bed-type output file> \n";
	print STDERR "converts PLN data into bed file\n";
	exit;
}

open IN, $ARGV[0] or die "Could not open file $ARGV[0]\n";

my $outputFile = $ARGV[1];
my $chr = "";
my $name = "";

open OUT, ">$outputFile" ;	
while (<IN>) {
	chomp;
	s/\r//g;
	my @line= split /\t/;
	
	$chr = "chr" . $line[0];
	$name = $line[4] . " to " . $line[3] . " alleles";

	print OUT $chr;
	print OUT "\t$line[1]\t$line[2]\t$name\t0\t+\n";
	}

close IN;

exit;

