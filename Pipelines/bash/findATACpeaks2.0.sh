#!/bin/bash
# by Michael Rehli, March 9, 2018
# updated with additional options , November 6, 2018

###################
# script for peak finding in ATACseq data
# will return three files - one with centered peak, one with broad region and one intersect

# update: can return a fourth file which contains extended "small" peaks for motif searches
###################

#setting homer environment

DIR_PKG="/misc/software/ngs"

PATH_PERL=/misc/software/package/perl/perl-5.26.1/bin
PATH_SAMTOOLS=${DIR_PKG}/samtools/samtools-1.6/bin
PATH_HOMER=${DIR_PKG}/homer/v4.9/bin
PATH_R=/misc/software/package/RBioC/3.4.3/bin

export PATH=${PATH_R}:${PATH_PERL}:${PATH_SAMTOOLS}:${PATH_HOMER}:${PATH}
export PATH

BEDTOOLS="/misc/software/ngs/bedtools/bedtools2-2.27.1/bin/bedtools"
TMPDIR="/loctmp"


# minimal required input:
# 
# -t <tagdir>
# -o <outputdir>
#
# optional:
# -l local fold change passed to findPeaks
# -f fdr passed to findPeaks


# Set Script Name variable
SCRIPT=`basename ${BASH_SOURCE[0]}`

# Set fonts for Help.
NORM=`tput sgr0`
BOLD=`tput bold`
REV=`tput smso`

# USAGE function
function USAGE {
  echo -e \\n"Help documentation for ${BOLD}${SCRIPT}.${NORM}"\\n
  echo -e "${REV}Basic usage:${NORM}\n${BOLD}$SCRIPT -t <tagdir> -o <outputdir> -g <genome> ${NORM}"\\n

  echo "Required:"
  echo "${BOLD}-t${NORM} <tagDir> Homer-style tag directory for ATAC data"
  echo "${BOLD}-o${NORM} <outputDir> target folder for resulting peak files"
  echo "${BOLD}-g${NORM} <genome> available genomes: hg19 hg38 mm9 mm10"
  echo -e \\n"Optional:"
  echo "${BOLD}-f${NORM} <FDR value> passed to findpeaks for peak calling (optional/default 0.00001)"
  echo "${BOLD}-l${NORM} <fold change> local fold change passed to findPeaks (optional/default 2)"
  echo "${BOLD}-s${NORM} also determine small peaks based on fragLength 1 (optional/default: not done)"

  echo -e "Example: ${BOLD}$SCRIPT -t MO_18h -o /misc/data/analysis/someProject/ATACpeaks -g hg19 -l 3 ${NORM}"\\n
  exit 1
}

#Check the number of arguments. If none are passed, print help and exit.
NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
  USAGE
fi

# Set defaults
FDR=0.00001
LOCAL=2
SMALLFLAG=0

# Defining available genomes
GENOMES=(hg19 hg38 mm9 mm10)


# Parse command line options
while getopts :g:t:o:l:f:sh OPTIONS; do
  case $OPTIONS in
    g)  #set option "g" - Bowtie will will only accept GRCh38 while HOMER needs hg38
      GENOME=$OPTARG
      ;;
    t)  #set option "t"
      TAGDIR=$OPTARG
      SAMPLEBASE=${TAGDIR##*/}     
      ;;
    o)  #set option "o"
      PEAKDIR=$OPTARG
      ;;
    f)  #set option "f"
      FDR=$OPTARG
      ;;
    l)  #set option "l"
      LOCAL=$OPTARG
      ;;
    s)  #show help
      SMALLFLAG=1
      ;;
    h)  #show help
      USAGE
      ;;
    \?) #unrecognized option - show USAGE
      echo -e \\n"Option -${BOLD}$OPTARG${NORM} not allowed."
      USAGE
      ;;
  esac
done
shift $((OPTIND-1)) 


# Check whether the tagDir exists

if [ ! -d "$TAGDIR" ]; then
	echo "TagDir does not exist!"
    echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
fi

# Check whether the output directory exists - - if not available, asks whether to make it

if [ ! -d "$PEAKDIR" ]; then
	echo "$PEAKDIR"
	read -p "Do you want to create this output folder (y/n)?" choice
	case "$choice" in 
  		y|Y ) 	mkdir -p "$PEAKDIR"
		;;
  		n|N )   echo "Didn't create the folder and stopped!"
    	        echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
    	        exit 8
		;;
  		* ) echo "invalid response";;
	esac
fi	

# Check whether GENOME is available

INARRAY=$(echo ${GENOMES[@]} | grep -o "$GENOME" | wc -w)

if [ $INARRAY == 0 ] ; then
      echo -e \\n"Genome -${BOLD}$GENOME${NORM} not available."
      echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
      exit 2
fi

# path to chromosome sizes 
CHROMSIZES="/misc/software/viewer/IGV/IGVTools_2.3.98/genomes/$GENOME.chrom.sizes"

case $GENOME in
	hg38)	  CHROMSIZES="/misc/software/viewer/IGV/IGVTools_2.3.98/genomes/GRCh38.PRI_p10.chrom.sizes"
	;;
esac


# Stitching to regions with nucleosome size peak sizes with stringent fdr cut-off

findPeaks ${TAGDIR} -size 150 -minDist 250 -L ${LOCAL} -fdr ${FDR} -region -o "${PEAKDIR}/${SAMPLEBASE}.s150.md250.L2.fdr-5.region.peaks.txt"
pos2bed.pl "${PEAKDIR}/${SAMPLEBASE}.s150.md250.L2.fdr-5.region.peaks.txt" > "${PEAKDIR}/${SAMPLEBASE}.s150.md250.L2.fdr-5.region.peaks.bed"

# centered factor peaks with 250bp peak sizes with stringent fdr cut-off and peak size equal minDist

findPeaks ${TAGDIR} -size 250 -minDist 250 -L ${LOCAL} -fdr ${FDR} -style factor -o "${PEAKDIR}/${SAMPLEBASE}.s250.md250.L2.fdr-5.factor.peaks.txt"
pos2bed.pl "${PEAKDIR}/${SAMPLEBASE}.s250.md250.L2.fdr-5.factor.peaks.txt" > "${PEAKDIR}/${SAMPLEBASE}.s250.md250.L2.fdr-5.factor.peaks.bed"

# intersecting bed files

$BEDTOOLS intersect -a "${PEAKDIR}/${SAMPLEBASE}.s150.md250.L2.fdr-5.region.peaks.bed" -b "${PEAKDIR}/${SAMPLEBASE}.s250.md250.L2.fdr-5.factor.peaks.bed" -u > "${PEAKDIR}/${SAMPLEBASE}.intersect.peaks.bed"
bed2pos.pl "${PEAKDIR}/${SAMPLEBASE}.intersect.peaks.bed" > "${PEAKDIR}/${SAMPLEBASE}.intersect.peaks.txt"

# finding smaler peaks (using fraglength 1)


if [ $SMALLFLAG==1 ] ; then

findPeaks ${TAGDIR} -size 12 -fragLength 1 -minDist 16 -L 0 -region -o "${PEAKDIR}/${SAMPLEBASE}.fl1.s12.md16.L0.region.peaks.txt"
pos2bed.pl "${PEAKDIR}/${SAMPLEBASE}.fl1.s12.md16.L0.region.peaks.txt" > "${PEAKDIR}/${SAMPLEBASE}.fl1.s12.md16.L0.region.peaks.bed"
$BEDTOOLS intersect -a "${PEAKDIR}/${SAMPLEBASE}.fl1.s12.md16.L0.region.peaks.bed" -b "${PEAKDIR}/${SAMPLEBASE}.intersect.peaks.bed" -u > ${TMPDIR}/tmp.inter.bed
$BEDTOOLS slop -i ${TMPDIR}/tmp.inter.bed -g $CHROMSIZES -b 48 > ${TMPDIR}/tmp.ext.bed
sort -k1,1 -k2,2n ${TMPDIR}/tmp.ext.bed > ${TMPDIR}/tmp.ext.sort.bed
$BEDTOOLS merge -i ${TMPDIR}/tmp.ext.sort.bed > ${TMPDIR}/tmp.ext.merge.bed
$BEDTOOLS slop -i ${TMPDIR}/tmp.ext.merge.bed -g $CHROMSIZES -b -24 > "${PEAKDIR}/${SAMPLEBASE}.extended.smallpeaks.bed"
bed2pos.pl "${PEAKDIR}/${SAMPLEBASE}.extended.smallpeaks.bed" > "${PEAKDIR}/${SAMPLEBASE}.extended.smallpeaks.pos.txt"

fi

