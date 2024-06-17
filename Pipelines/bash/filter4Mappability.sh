#!/bin/bash
# by Michael Rehli, March 17, 2018
# bash script to filter-out regions with low mappability 
# 

#setting homer environment

DIR_PKG="/misc/software/ngs"

PATH_PERL=/misc/software/package/perl/perl-5.26.1/bin
PATH_SAMTOOLS=${DIR_PKG}/samtools/samtools-1.6/bin
PATH_HOMER=${DIR_PKG}/homer/v4.9/bin
PATH_R=/misc/software/package/RBioC/3.4.3/bin

export PATH=${PATH_R}:${PATH_PERL}:${PATH_SAMTOOLS}:${PATH_HOMER}:${PATH}
export PATH

# Defining the program versions (needs to be adjusted to the server/workstation used)
# rhskl17 version

BOWTIE="/misc/software/ngs/bowtie/bowtie2-2.3.4-linux-x86_64/bowtie2"
SAMTOOLS="/misc/software/ngs/samtools/samtools-1.6/bin/samtools"
SAMBAMBA="/misc/software/ngs/sambamba/v0.6.7/sambamba"
BEDTOOLS="/misc/software/ngs/bedtools/bedtools2-2.27.1/bin/bedtools"
PICARD="/usr/bin/java -jar /misc/software/ngs/picard/src/v2.17.3/picard.jar"
IGVTOOLS="/usr/bin/java -Xmx4g -Djava.awt.headless=true -jar /misc/software/viewer/IGV/IGVTools_2.3.98/igvtools.jar"

# required input:
# 
# -p <peak file> accepts position or bed files
# -g <genome> available genomes include hg19, GRCh38, mm9, mm10, default GRCh38
# -s <mapping size> currently accepts 50, 75, 100 (available mappability traacks)
# -f <cut-off> filter for mappability score of regions [range:0-1], default 0.8
#

# Set Script Name variable
SCRIPT=`basename ${BASH_SOURCE[0]}`

# Set fonts for Help.
NORM=`tput sgr0`
BOLD=`tput bold`
REV=`tput smso`

# USAGE function
function USAGE {
  echo -e \\n"Help documentation for ${BOLD}${SCRIPT}.${NORM}"\\n
  echo -e "${REV}Basic usage:${NORM}\n${BOLD}$SCRIPT -p <peak file> -g <genome> -f <cut-off> -s <mapping size>"\\n
  echo "Required:"
  echo "${BOLD}-p${NORM} <peak/region file> accepts position or bed files"
  echo "${BOLD}-s${NORM} <mapping size> currently accepts 50, 75, 100 (available mappability tracks), default:50"
  echo -e \\n"Optional:"
  echo "${BOLD}-g${NORM} <genome>  available genomes include hg19, GRCh38, mm9, mm10, GRCm38 default:GRCh38"
  echo "${BOLD}-f${NORM} <cut-off> filter for mappability score of regions [range:0-1], default:0.8"
  echo -e \\n"Example: ${BOLD}$SCRIPT -p peaks.txt -g hg19 -f 0.9 -s 50 ${NORM}"\\n
  exit 1
}

#Check the number of arguments. If none are passed, print help and exit.
NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
  USAGE
fi

# Defining available data sets
GENOMES=(hg19 GRCh38 hg38 mm9 mm10 GRCm38)
MAPPABILITY=(50 75 100)

# Set defaults
GENOME="GRCh38"
CUTOFF=0.8
SIZE=50

# set path
MAPDIR="/misc/data/analysis/generalStuff/mappability"
SIZES_DIR="/misc/software/viewer/IGV/IGVTools_2.3.98/genomes"
TEMPFILES="/loctmp/$RANDOM.tmp"

# Parse command line options
while getopts :g:f:p:s:h OPTIONS; do
  case $OPTIONS in
    p)  #set option "p"
      PEAKS=$OPTARG
      ;;
    g)  #set option "g" - Bowtie will will only accept GRCh38 while HOMER needs hg38
      GENOME=$OPTARG
      ;;
    f)  #set option "f"
      CUTOFF=$OPTARG
      ;;
    s)  #set option "s"
      SIZE=$OPTARG
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

# Sanity checks

# Check whether GENOME is available

INARRAY=$(echo ${GENOMES[@]} | grep -o "$GENOME" | wc -w)

if [ $INARRAY == 0 ] ; then
      echo -e \\n"Genome -${BOLD}$GENOME${NORM} not available."
      echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
      exit 2
fi

INMAP=$(echo ${MAPPABILITY[@]} | grep -o "$SIZE" | wc -w)

if [ $INMAP == 0 ] ; then
      echo -e \\n"Mappability for size ${BOLD}$SIZE${NORM} not available."
      echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
      exit 2
fi

# Define mappability track and chrom.sizes and genome version for available genomes

case $GENOME in
	hg19)   MAPTRACK=${MAPDIR}/hg19.mappability_${SIZE}.bedGraph
			HOMERGENOME=hg19
	;;
	mm9)   MAPTRACK=${MAPDIR}/mm9.mappability_${SIZE}.bedGraph
			HOMERGENOME=mm9
	;;
	hg38|GRCh38)   MAPTRACK=${MAPDIR}/hg38.mappability_${SIZE}.bedGraph
			HOMERGENOME=hg38
	;;
	mm10|GRCm38)   MAPTRACK=${MAPDIR}/mm10.mappability_${SIZE}.bedGraph
			HOMERGENOME=mm10
	;;
esac

# Run annotatePeaks.pl

annotatePeaks.pl ${PEAKS} ${HOMERGENOME} -size given -bedGraph ${MAPTRACK} -nogene -noann -cpu 3 > ${TEMPFILES}.ann.peaks.txt

# Filter out mappability below cut-off

FILTEREDPEAKS=${PEAKS%.*}.mapScoreFiltered.txt

myFilterFile.pl ${TEMPFILES}.ann.peaks.txt -column 8 -lowerlimit ${CUTOFF} > ${FILTEREDPEAKS}

LENGTHA="$(wc -l <"$PEAKS")"
LENGTHB="$(wc -l <"$FILTEREDPEAKS")"
DIFFERENCE=$((LENGTHA-LENGTHB))

echo -e \\n"Done. Filtered out ${DIFFERENCE} regions with mappability scores below ${CUTOFF}, leaving ${LENGTHB} regions."\\n

# removing tmp files

rm ${TEMPFILES}.*



