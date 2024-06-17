#!/usr/bin/perl 

#########################################################################
#																		#
#   Perl script to generate an averaged bigWig from replicate samples	#  
#																		#
#   Note that the script just adds up counts and divides each signal 	#
#	by the number of samples!											#
#																		#
#	If you need to combine negative strands please add the minus 		#
#	option. This will invert each bigWig before processing and invert	#
#	again after processing. 									 		#
#																		#
#########################################################################

if (@ARGV < 3) {
	print STDERR "myAverageBigWig.pl -bw <bw-file1> <bw-file2> -chr <chrom.sizes> -o\n";
	print STDERR "\n";
	print STDERR "option\n";
	print STDERR "-minus inverts all values before merging\n";
	print STDERR "\n";
	exit;
}

my @bwFiles = ();
my @invbwFiles = ();
my $output = '';
my $rand = rand();
my $chrsizes = '';
my $inv = 0;
my $tmp = '';

for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-bw') {
		my $bail = 0;
		print STDERR "\tmerging and averaging the following bigWig files:\n";
		while ($ARGV[++$i] !~ /^\-/) {
			push(@bwFiles, $ARGV[$i]);
			print STDERR "\t\t$ARGV[$i]\n";
			if ($i>=@ARGV-1) {
				$bail = 1;
				last;
			}
		}
		last if ($bail == 1);
		$i--;
	} elsif ($ARGV[$i] eq '-o') {
		$output = $ARGV[++$i];
		print STDERR "\n\toutput file: $output\n";

	} elsif ($ARGV[$i] eq '-chr') {
		$chrsizes = $ARGV[++$i];
		print STDERR "\n\tchromosome size from: $chrsizes\n";

	} elsif ($ARGV[$i] eq '-minus') {
		$inv = 1;
	} 
}

my $fileCount = scalar(@bwFiles);
my $files = join(" ", @bwFiles);
my $tmp1 = $rand . ".bedGraph";
my $tmp2 = $rand . ".ave.bedGraph";
my $tmp3 = $rand . ".sorted.bedGraph";
my $tmp4 = $rand . ".inv.bedGraph";

if ($inv == 1) {

	for (my $i=0;$i<@bwFiles;$i++) {
		print STDERR "\n\tinverting bigWig $i...\n";		
		$invbwFiles[$i] = $rand . $i . ".bw";
		`bigWigToBedGraph $bwFiles[$i] $tmp1 `;
		`LC_NUMERIC=\"en_US.UTF-8\" awk \'BEGIN{OFS=\"\t\"}{print \$1,\$2,\$3,\$4*-1}\' $tmp1 > $tmp2`;
		`LC_COLLATE=C sort -k1,1 -k2,2n $tmp2 > $tmp3`;
		`bedGraphToBigWig $tmp3 $chrsizes $invbwFiles[$i]`;
	}
	$files = join(" ", @invbwFiles);
}

print STDERR "\n\tmerging bigWigs...\n";		

	`bigWigMerge $files $tmp1`;

print STDERR "\n\taveraging counts...\n";		

	`LC_NUMERIC=\"en_US.UTF-8\" awk \'BEGIN{OFS=\"\t\"}{print \$1,\$2,\$3,\$4/$fileCount}\' $tmp1 > $tmp2`;
	
print STDERR "\n\tsorting bedgraph...\n";		

	`LC_COLLATE=C sort -k1,1 -k2,2n $tmp2 > $tmp3`;
	
if ($inv == 1) {

	print STDERR "\n\tinverting counts and creating bigWig...\n";		
	`LC_NUMERIC=\"en_US.UTF-8\" awk \'BEGIN{OFS=\"\t\"}{print \$1,\$2,\$3,\$4*-1}\' $tmp3 > $tmp4`;
	`bedGraphToBigWig $tmp4 $chrsizes $output`;

} else {
	
	print STDERR "\n\tcreating bigWig...\n";		
	`bedGraphToBigWig $tmp3 $chrsizes $output`;
	
}
	
print STDERR "\n\tremoving temporary files...\n";		
				
		`rm $rand.*`;

exit;

