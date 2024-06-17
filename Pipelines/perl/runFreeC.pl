#!/usr/bin/perl -w

if (@ARGV < 5) {
	print STDERR "runFreeC.pl <output directory> <sex (XY/XX)> <genome> <sam-file directory> <file name> <Path to chrom.sizes> <Path to gem file>\n";
	exit;
}
my $outputDir = $ARGV[0];
my $gender = $ARGV[1];
my $genome = $ARGV[2];
my $samDir = $ARGV[3];
my $samFileName = $ARGV[4];
my $chromsizes = $ARGV[5];
my $gempath = $ARGV[6];

my $r = rand();
my $tmpConfig = "/loctmp/" . $r . ".txt";

open OUT, ">$tmpConfig";
print OUT "[general]\n";
print OUT "chrFiles = /misc/software/ngs/homer/v4.9/data/genomes/$genome\n";
print OUT "chrLenFile = $chromsizes\n";
print OUT "coefficientOfVariation = 0.05\n";
print OUT "gemMappabilityFile = $gempath\n";
print OUT "maxThreads = 4\n";
print OUT "outputDir = $outputDir\n";
print OUT "ploidy = 2,3,4\n";
print OUT "sex = $gender\n";
print OUT "telocentromeric = 50000\n";
print OUT "uniqueMatch = TRUE\n";
print OUT "#contaminationAdjustment = TRUE\n";
print OUT "#contamination .5\n\n";
print OUT "[sample]\n";
print OUT "mateFile = $samDir/$samFileName.sam\n";
print OUT "inputFormat = sam\n";
print OUT "mateOrientation = 0\n";
print OUT "[control]\n\n";
print OUT "[BAF]\n\n";
print OUT "[target]\n";
close OUT;

	print STDERR "\nrunning FreeC for CNV analysis \n";

`/misc/software/ngs/freec/FREEC-11.0/src/freec -conf "$tmpConfig"`;

	print STDERR "creating Plot with FreeChistogram_$genome.pl, output: $outputDir/$samFileName.pdf \n";
	
`FreeChistogram_$genome.pl "$outputDir"/"$samFileName".sam_ratio.txt "$outputDir" "$samFileName" "$genome" "$gender"`;
	
	
	print STDERR "removing temporary scripts and files \n\n";


 `rm -f \"$tmpConfig\" `;
