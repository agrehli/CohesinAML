#!/usr/bin/env perl
use warnings;
use lib "/misc/software/ngs/homer/v4.11/.//bin";

my $tmpdir = "/loctmp";

if (@ARGV < 3) {
	print STDERR "Script to either add a constant value or scale bedGrap/bigWig files\n";
	print STDERR "\n";
	print STDERR "myScaleBedGraph <input bedGraph file> <options> -o <outputfile w/o ext.>\n";
	print STDERR "\toptions:\n";
	print STDERR "\t\t-S <#> add factor (default: 1)\n";
	print STDERR "\t\t-M <#> multiply by factor\n";
	print STDERR "\t\t-bigWig <genome> (IN and OUT are bigWig; available genomes: hg19, hg38, mm10))\n";
	exit;
}

my $input = $ARGV[0];
my $scale = 1;
my $factor = 1;
my $rand = rand();

my $job = "add" ;
my $type = "bedGraph";
my $outputname = "";
my $output = "";
my $genome = "";
my $chrsizes = "";

for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-S') {
		$scale = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-M') {
		$factor = $ARGV[++$i];
		$job = "multiply";
	} elsif ($ARGV[$i] eq '-o') {
		$outputname = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-bigWig') {
		$type = "bigWig";
		$genome = $ARGV[++$i];
	}
}

	print STDERR "\n";
	print STDERR "input file: $input\n";
	print STDERR "file type: $type ";
	print STDERR "genome: $genome \n";
	print STDERR "output name: $outputname\n";
	if ($job eq "add"){
		print STDERR "value to add: $scale\n";
	} else {
	    print STDERR "scaling factor: $factor\n";
	}


if ($type eq "bedGraph"){
    $output = $outputname . ".bedGraph";

	print STDERR "\nscaling bedGraph\n\n";


	open OUT, ">$output" ;		
	open IN, $input or die "Could not open file $ARGV[0]\n";
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line= split /\t/;
		if ($job eq "add"){
		$adjustedcount = $line[3]+$scale;
		} else {
		$adjustedcount = $line[3]*$factor;
		}
		print OUT $line[0];
		print OUT "\t$line[1]\t$line[2]\t$adjustedcount\n";
	}
	close IN;
	close OUT; 

} elsif ($type eq "bigWig"){
    $output = $outputname . ".bigWig";

    if ($genome eq "hg19") {
    $chrsizes = "/misc/software/viewer/IGV/IGVTools_2.3.98/genomes/hg19.chrom.sizes";
    } elsif ($genome eq "hg38") {
    $chrsizes = "/misc/software/viewer/IGV/IGVTools_2.3.98/genomes/GRCh38.PRI_p10.chrom.sizes";
   } elsif ($genome eq "mm10") {
    $chrsizes = "/misc/software/viewer/IGV/IGVTools_2.3.98/genomes/mm10.chrom.sizes";
    } else {
	print STDERR "\tgenome $genome not available))\n";
	exit;
    
    }

    my $tmpin = $tmpdir . "/" . $rand . ".in.bg" ;
    my $tmpout = $tmpdir . "/" . $rand . ".out.bg" ;

	print STDERR "\nconverting bigWig to bedGraph\n\n";

	`bigWigToBedGraph $input $tmpin`;

	print STDERR "\nscaling bedGraph\n\n";

	open OUT, ">$tmpout" ;		
	open IN, $tmpin or die "Could not open file $ARGV[0]\n";
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line= split /\t/;
		if ($job eq "add"){
		$adjustedcount = $line[3]+$scale;
		} else {
		$adjustedcount = $line[3]*$factor;
		}
		print OUT $line[0];
		print OUT "\t$line[1]\t$line[2]\t$adjustedcount\n";
	}
	close IN;
	close OUT; 
	
    print STDERR "\nconverting bedGraph to bigWig\n\n";

	`bedGraphToBigWig $tmpout $chrsizes $output`;
	`rm $tmpin $tmpout`;

	
} else {
	print STDERR "\twrong file format?))\n";
	exit;

}






exit;

