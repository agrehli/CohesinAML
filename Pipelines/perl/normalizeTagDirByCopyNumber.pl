#!/usr/bin/env perl
use warnings;
use lib "/misc/software/ngs/homer/v4.9/.//bin";

use POSIX;
if (@ARGV < 1) {
	print STDERR "\n\tThis program is used to normalize tagDirectories by CNV\n";
	print STDERR "\t   normalizeTagDirByCopyNumber.pl <tagDir> -cnv <CNV file>\n";
	print STDERR "option: -remove (removes scaffolds from tagDir)\n";
	print STDERR "\n";
	exit;
}

my $tagDirectory = $ARGV[0] ;
my $CNVfile = '' ;
my $CNV = 0;
my $remove = 0;
my $rand = rand();
my $tmpFreeCBedFile = $rand . ".FreeC.bed";
my $tmpBedFile = $rand . ".tagDir.bed";
my $tmpOverFile = $rand . ".over.bed";
my $tmpReformatedBedFile = $rand . ".reformated.bed";
my $tagDirectoryNorm = $ARGV[0] . "_norm";
my $BEDTOOLS = "/misc/software/ngs/bedtools/bedtools2-2.27.1/bin/bedtools";

for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-cnv') {
		$CNVfile = $ARGV[++$i];
		$CNV = 1;
	} elsif ($ARGV[$i] eq '-remove') {
		$remove = 1;
	}
}

if ($CNV == 0) {	
	$tagDirectoryNorm = $ARGV[0] . "_refChr";
} elsif ($remove == 1) {
	$tagDirectoryNorm = $ARGV[0] . "_CNVnormRefChr";
} else {
	$tagDirectoryNorm = $ARGV[0] . "_CNVnorm";
}
	
	
my $value = 0 ;

if ($CNV == 0){

print STDERR "\n\tremoving scaffolds from tagDir $tagDirectory\n\n";


	`tagDir2bed.pl $tagDirectory >$tmpBedFile`;

open OUT, ">$tmpReformatedBedFile" ;		
open IN, $tmpBedFile ;		
while (<IN>) {
	chomp;
	s/\r//g;
	my @line= split /\t/;
	
	if ($line[0] =~ /_gl/) {
	next;
	} else {
			print OUT $line[0];
			print OUT "\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\n";
			}
	}
close IN;		
close OUT;

} else {

print STDERR "\n\treformating and annotating tagDir $tagDirectory\n";
print STDERR "\twith CNV info $CNVfile\n\n";

	`myConvertFreeC2Bed.pl $CNVfile > $tmpFreeCBedFile`;

	`tagDir2bed.pl $tagDirectory >$tmpBedFile`;

	`$BEDTOOLS intersect -a $tmpBedFile -b $tmpFreeCBedFile -wao -f 0.5 > $tmpOverFile`;

print STDERR "\n\treformating file\n\n";

open OUT, ">$tmpReformatedBedFile" ;		
open IN, $tmpOverFile ;		
while (<IN>) {
	chomp;
	s/\r//g;
	my @line= split /\t/;
	
	if ($remove == 1 && $line[0] =~ /_gl/) {
	next;
	}
	
		if ($line[7] == -1 ) {
			
			print OUT $line[0];
			print OUT "\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\n";

		} else {
			
			if ($line[10] == 0) {
			$value = 0;
			} else {
			$value = $line[4] * ( 2/ $line[10] );
			}
		
			print OUT $line[0];
			print OUT "\t$line[1]\t$line[2]\t$line[3]\t$value\t$line[5]\n";
		
			}
	}
close IN;		
close OUT;

}


print STDERR "\n\tcreating (CNV normalized) tagDir without scaffolds \n\n";

	`makeTagDirectory $tagDirectoryNorm $tmpReformatedBedFile -format bed -force5th -precision 3 `;

print STDERR "\n\tremoving temporary files...\n\n\n";		
				
	`rm "$rand".*`;

exit;

