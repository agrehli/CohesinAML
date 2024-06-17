#!/bin/bash
# by Michael Rehli, Febr 23, 2018
# bash script to map ATAC-seq data, generate tagDirectories, coverage bigWigs, 
# and to find ATAC-seq peaks
# 

# Nov 20, 2023 changed script to allow adding multiple fastq files, updated R & Samtools versions





#setting homer environment

DIR_PKG="/misc/software/ngs"

PATH_PERL=/misc/software/package/perl/perl-5.26.1/bin
PATH_SAMTOOLS=${DIR_PKG}/samtools/samtools-1.16.1/bin
PATH_HOMER=${DIR_PKG}/homer/v4.11/bin
PATH_R=${DIR_SOFT}/package/RBioC/buster_R-4.3.0_Bioc-3.17_intel-mkl-64bit-2020.1-102/lib/R/bin

export PATH=${PATH_R}:${PATH_PERL}:${PATH_SAMTOOLS}:${PATH_HOMER}:${PATH}
export PATH

# Defining the program versions (needs to be adjusted to the server/workstation used)
# rhskl17 version

SKEWER="${DIR_PKG}/skewer/skewer-0.2.2/skewer"
BOWTIE="${DIR_PKG}/bowtie/bowtie2-2.3.4-linux-x86_64/bowtie2"
SAMTOOLS="${DIR_PKG}/samtools/samtools-1.16.1/bin/samtools"
PICARD="/usr/bin/java -jar /misc/software/ngs/picard/src/v2.17.3/picard.jar"
IGVTOOLS="/usr/bin/java -Xmx4g -Djava.awt.headless=true -jar /misc/software/viewer/IGV/IGVTools_2.3.98/igvtools.jar"

# required input:
# 
# -f <fastq or fastq.gz> forward read 
# -r <fastq or fastq.gz> reverse read (if not available, will run as single end)
# -g <genome>  available genomes include hg19, GRCh38, mm9, mm10
# -o <output subdirectory> only subdirectories above ATAC-folder are given
# -n <sample name> used for sam, bigwig and tag directory
#
# optional:
#
# -t <number of threads>
# -s <path to adapter sequences> (or type NEXTERA to use Nextera-adapters)

# Set Script Name variable
SCRIPT=`basename ${BASH_SOURCE[0]}`

# Set fonts for Help.
NORM=`tput sgr0`
BOLD=`tput bold`
REV=`tput smso`

# USAGE function
function USAGE {
  echo -e \\n"Help documentation for ${BOLD}${SCRIPT}.${NORM}"\\n
  echo -e "${REV}Basic usage:${NORM}\n${BOLD}$SCRIPT -g <genome> -f <fastq> -r <fastq> -o <subdirectory> -n <name> -t <threads>${NORM}"\\n
  echo "Required:"
  echo "${BOLD}-g${NORM} <genome>  available genomes include hg19, GRCh38, mm9, mm10"
  echo "${BOLD}-f${NORM} <fastq or fastq.gz> forward read (several files can be given, separated by komma w/o space)"
  echo "${BOLD}-r${NORM} <fastq or fastq.gz> reverse read (several files can be given, separated by komma w/o space. if not available, will run as single end)"
  echo "${BOLD}-d${NORM} <output subdirectory> only subdirectories above ATAC-folder are given"
  echo "${BOLD}-n${NORM} <sample name> used for sam, bigwig and tag directory"
  echo -e "${BOLD}-t${NORM} <number of threads> (optional/default 16)"\\n
  echo -e "${BOLD}-s${NORM} <path to adapter sequences> (or type NEXTERA to use Nextera-adapters)"\\n
  echo -e "Example: ${BOLD}$SCRIPT -g hg19 -f readA.fastq -r readB.fastq -d MOMACDC -n test -t 32 -s NEXTERA${NORM}"\\n
  exit 1
}

# Note:
echo -e \\n"Please note"\\n
echo -e "This is an updated version of the mapATAC.sh script!"\\n
echo "Changes:"
echo "- updated R & Samtools versions"
echo "- now allows the simultanous mapping of several fastq files"
echo -e "Example: ${BOLD}$SCRIPT -g hg38 -f readA1.fastq,readA2.fastq -r readB1.fastq,readB2.fastq -d MOMACDC -n test -t 32 -s NEXTERA${NORM}"\\n
echo "${BOLD}If you notice inconsistencies, please report to Michael!${NORM}"
echo -e \\n"If you like to use the old version refer to mapATAC.v1.sh"\\n

#Check the number of arguments. If none are passed, print help and exit.
NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
  USAGE
fi

# Defining available genomes
#GENOMES=(hg19 GRCh38 hg38 mm9 mm10 GRCm39)
GENOMES=(hg19 GRCh38 hg38 mm9 mm10)

# Set default number of threads to 16
THREADS=16
TRIM=0
MULTI=0

# Parse command line options
while getopts :g:f:r:d:n:t:s:h OPTIONS; do
  case $OPTIONS in
    g)  #set option "g" - Bowtie will will only accept GRCh38 while HOMER needs hg38
      GENOME=$OPTARG
      HOMERGENOME=$OPTARG
      ;;
    f)  #set option "f"
      FASTQR1=$OPTARG
      ;;
    r)  #set option "r"
      FASTQR2=$OPTARG
      ;;
    d)  #set option "d"
      DIRECTORY=$OPTARG
      ;;
    n)  #set option "n"
      SAMPLENAME=$OPTARG
      ;;
    t)  #set option "t"
      THREADS=$OPTARG
      SAMTHREADS=$((THREADS-1))
      ;;
    s)  #set option "s"
      ADAPTER=$OPTARG
      TRIM=1
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

# Check whether a reasonable number of cores is given

if [ $THREADS -lt 1 ] || [ "$THREADS" -gt 63 ]; then
	echo "Number of threads should be in the range of 1 - 63!"
    echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
    exit 9
fi


# Check whether GENOME is available

INARRAY=$(echo ${GENOMES[@]} | grep -o "$GENOME" | wc -w)

if [ $INARRAY == 0 ] ; then
      echo -e \\n"Genome -${BOLD}$GENOME${NORM} not available."
      echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
      exit 2
fi

GENOMEPATH=$GENOME
HOMERGENOME=$GENOME
BIGWIGGENOME=$GENOME

# path to chromosome sizes 
CHROMSIZES="/misc/software/viewer/IGV/IGVTools_2.3.98/genomes/$GENOME.chrom.sizes"

if [ $GENOME == "hg38" ] || [ $GENOME == "GRCh38" ] ; then
      GENOME="hg38"
      GENOMEPATH="GRCh38.PRI_p10"
      BIGWIGGENOME="GRCh38"
      HOMERGENOME="hg38"
      CHROMSIZES="/misc/software/viewer/IGV/IGVTools_2.3.98/genomes/GRCh38.PRI_p10.chrom.sizes"
fi

# if [ $GENOME == "GRCm39" ] ; then
#       GENOME="hg38"
#       GENOMEPATH="GRCm39.PRI_M33"
#       BIGWIGGENOME="GRCm39"
#       HOMERGENOME="hg38"
#       CHROMSIZES="${DIR_PKG}/genome/sequence/GRCm39.PRI_M33/chrom.sizes"
# fi

# Bowtie index
INDEX="$GENOMEPATH/$GENOME"

 
# Check whether:
# – more than one fastq is provided (added Nov 20, 2023 by MR), 
# – whether one or two reads are given, 
# – whether files are available, non-identical and have the right extensions

## Read1

if [[ "$FASTQR1" =~ "," ]]; then
    echo -e \\n"... multiple Fastqs given:"\\n
    MULTI=1
    R1ARRAY=($(echo $FASTQR1 | tr "," "\n"))
    for R1 in "${R1ARRAY[@]}"; do
        if [ ! -f "$R1" ]; then
            echo "Fastq file $R1 (read 1) not found!"
            echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
            exit 3
        fi
        case "$R1" in
            *.fastq.gz) ;;
            *.fastq)    ;;
            *.fq.gz)    ;;
            *.fq)       ;;
            *)          echo "Fastq file (read 1) has wrong extension (should be fastq, fq, fastq.gz or fq.gz!"
            echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
            exit 4
            ;;
        esac
    done
elif [ ! -f "$FASTQR1" ]; then
    echo "Fastq file (read 1) not found!"
    echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
    exit 3
else
    case "$FASTQR1" in
        *.fastq.gz) ;;
        *.fastq)    ;;
        *.fq.gz)    ;;
        *.fq)       ;;
        *)          echo "Fastq file (read 1) has wrong extension (should be fastq, fq, fastq.gz or fq.gz!"
        echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
        exit 4
        ;;
    esac
fi

## Read2

if [ -z "$FASTQR2" ]; then
    SEQTYPE=1
    TYPE="SE"
    SKEWTYPE="any"
    echo -e \\n"Data processing for ${BOLD}single-end${NORM} reads"
elif [ "$FASTQR1" == "$FASTQR2" ]; then
    echo "Fastq files for read 1 and 2 are identical!"
    echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
    exit 6
elif [ $MULTI == 1 ]; then
    R2ARRAY=($(echo $FASTQR2 | tr "," "\n"))
    for R2 in "${R2ARRAY[@]}"; do
        if [ ! -f "$R2" ]; then
            echo "Fastq file $R2 (read 1) not found!"
            echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
            exit 3
        fi
        case "$R2" in
            *.fastq.gz) ;;
            *.fastq)    ;;
            *.fq.gz)    ;;
            *.fq)       ;;
            *)          echo "Fastq file (read 1) has wrong extension (should be fastq, fq, fastq.gz or fq.gz!"
            echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
            exit 4
            ;;
        esac
        # Check whether equal numbers of fastq files are given for both reads
        if [ ${#R2ARRAY[@]} != ${#R1ARRAY[@]} ] ; then 
            echo "File counts for read 1 (${#R1ARRAY[@]}) and 2 (${#R2ARRAY[@]}) don't match"
            echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
            exit 4
        fi
    done
    SEQTYPE=2
    TYPE="PE"
    SKEWTYPE="pe"
    echo -e \\n"Data Processing for ${BOLD}paired-end${NORM} reads"
else
    SEQTYPE=2
    TYPE="PE"
    SKEWTYPE="pe"
    echo -e \\n"Data Processing for ${BOLD}paired-end${NORM} reads"
    case "$FASTQR2" in
        *.fastq.gz) ;;
        *.fastq)    ;;
        *.fq.gz)    ;;
        *.fq)       ;;
        *)          echo "Fastq file (read 1) has wrong extension (should be fastq, fq, fastq.gz or fq.gz!"
        echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
        exit 4
        ;;
    esac
fi


# Check wether the subfolders are already there - if not ask whether path should be made
 
if [ -z "$DIRECTORY" ]; then
	echo "Subdirectory is missing!"
    echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
    exit 7
fi

# rhskl17 version
MAPPINGDIR="/misc/data/processedData/mapping/chromatin/${BIGWIGGENOME}/ATAC/${DIRECTORY}"
TAGDIRDIR="/misc/data/processedData/tagDir/chromatin/${BIGWIGGENOME}/ATAC/${DIRECTORY}"
BIGWIGDIR="/misc/data/processedData/bigWig/chromatin/${BIGWIGGENOME}/ATAC/${DIRECTORY}"

if [ ! -d "$MAPPINGDIR" ]; then			#looks for the mapping output folder - if not available, askes whether to make it
	echo "$MAPPINGDIR"
	read -p "Do you want to create this mapping folder (y/n)?" choice
	case "$choice" in 
  		y|Y ) 	mkdir -p "$MAPPINGDIR"
		;;
  		n|N )   echo "Didn't create the folder and stopped!"
    	        echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
    	        exit 8
		;;
  		* ) echo "invalid response";;
	esac
fi	
	
if [ ! -d "$TAGDIRDIR" ]; then			#looks for the tagDir output folder - if not available, askes whether to make it
	echo "$TAGDIRDIR"
	read -p "Do you want to create this tagDirectory folder (y/n)?" choice
	case "$choice" in 
  		y|Y ) 	mkdir -p "$TAGDIRDIR"
		;;
  		n|N )   echo "Didn't create the folder and stopped!"
    	        echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
    	        exit 8
		;;
  		* ) echo "invalid response";;
	esac
fi	
	
if [ ! -d "$BIGWIGDIR" ]; then			#looks for the bigWig output folder - if not available, askes whether to make it
	echo "$BIGWIGDIR"
	read -p "Do you want to create this bigWig folder (y/n)?" choice
	case "$choice" in 
  		y|Y ) 	mkdir -p "$BIGWIGDIR"
		;;
  		n|N )   echo "Didn't create the folder and stopped!"
    	        echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
    	        exit 8
		;;
  		* ) echo "invalid response";;
	esac
fi	
 
# Check whether sample name is given 
 
if [ -z "$SAMPLENAME" ]; then
	echo "Sample name is missing!"
    echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
    exit 9
fi

SAMPLENAME="${SAMPLENAME// /_}"  # replaces spaces with underscores in sample name, just in case...

# Check if path for trimming is given or set for NEXTERA

if [ $TRIM == 1 ] && [ $ADAPTER == "NEXTERA" ]; then
    ADAPTER="/misc/software/ngs/skewer/adapterFiles/nextera_adapters.fa"
elif [ $TRIM == 1 ] && [ ! -f "$ADAPTER" ]; then
    echo "Adapter file not found!"
    echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
    exit 11   
fi

if [ ! -d "${MAPPINGDIR}/logs" ]; then			#looks for the logs output folder - if not available, makes it
   mkdir "${MAPPINGDIR}/logs"
fi

# If everything passed, create a log file for all the commands that are executed during the processing

LOGFILE="${MAPPINGDIR}/logs/${SAMPLENAME}.commands.txt"

touch $LOGFILE

# Function to append command lines into logfile and then execute

exe() { echo -e "\$ $@"\\n >> $LOGFILE; "$@" ; }


# Generating a random number file name for temporary files

TEMPFILES="/loctmp/$RANDOM.tmp"

# Optionally starting the pipeline by trimming the reads with skewer

if [ $TRIM == 1 ] && [ $MULTI == 0 ] ; then
   case "$SEQTYPE" in
   1 )  exe $SKEWER -f sanger -t $THREADS -m any -x $ADAPTER -o $TEMPFILES $FASTQR1 
        FASTQR1="$TEMPFILES-trimmed.fastq"
        mv "$TEMPFILES-trimmed.log" "${MAPPINGDIR}/logs/${SAMPLENAME}.trimming.log"
        ;;
   2 )  exe $SKEWER -f sanger -t $THREADS -m pe -x $ADAPTER -o $TEMPFILES $FASTQR1 $FASTQR2
        FASTQR1="$TEMPFILES-trimmed-pair1.fastq"
        FASTQR2="$TEMPFILES-trimmed-pair2.fastq"
        mv "$TEMPFILES-trimmed.log" "${MAPPINGDIR}/logs/${SAMPLENAME}.trimming.log"
        ;;       
   esac     

elif [ $TRIM == 1 ] && [ $MULTI == 1 ] ; then
   case "$SEQTYPE" in
   1 )  COUNTER=0
        for R1 in "${R1ARRAY[@]}"; do
            exe $SKEWER -f sanger -t $THREADS -m pe -x $ADAPTER -o "$TEMPFILES.$COUNTER" $R1
            mv "$TEMPFILES.$COUNTER-trimmed.log" "${MAPPINGDIR}/logs/${SAMPLENAME}.$COUNTER.trimming.log"
            COUNTER=$((COUNTER+=1))
        done
        FASTQR1="$TEMPFILES.0-trimmed-pair1.fastq"
        for (( i=1; i<${#R1ARRAY[@]}; i++ )); do
            FASTQR1="$FASTQR1,$TEMPFILES.$i-trimmed-pair1.fastq"
        done
        ;;
   2 )  COUNTER=0
        for R1 in "${R1ARRAY[@]}"; do
            R2="${R2ARRAY[${COUNTER}]}"
            exe $SKEWER -f sanger -t $THREADS -m pe -x $ADAPTER -o "$TEMPFILES.$COUNTER" $R1 $R2
            mv "$TEMPFILES.$COUNTER-trimmed.log" "${MAPPINGDIR}/logs/${SAMPLENAME}.$COUNTER.trimming.log"
            COUNTER=$((COUNTER+=1))
        done
        FASTQR1="$TEMPFILES.0-trimmed-pair1.fastq"
        FASTQR2="$TEMPFILES.0-trimmed-pair2.fastq"
        for (( i=1; i<${#R1ARRAY[@]}; i++ )); do
            FASTQR1="$FASTQR1,$TEMPFILES.$i-trimmed-pair1.fastq"
            FASTQR2="$FASTQR2,$TEMPFILES.$i-trimmed-pair2.fastq"
        done
        ;;
    esac
fi


# Running Bowtie2

echo -e \\n"Running Bowtie2 on $SAMPLENAME."\\n

case "$SEQTYPE" in
	1 )  exe $BOWTIE --very-sensitive -p $THREADS -x $INDEX --met-file $MAPPINGDIR/logs/$SAMPLENAME.aln_metrics.txt -U $FASTQR1 2> $MAPPINGDIR/logs/$SAMPLENAME.aln_rates.txt -S $TEMPFILES.sam
	;;
	2 )  exe $BOWTIE --very-sensitive --no-discordant -p $THREADS -x $INDEX --met-file $MAPPINGDIR/logs/$SAMPLENAME.aln_metrics.txt -1 $FASTQR1 -2 $FASTQR2 2> $MAPPINGDIR/logs/$SAMPLENAME.aln_rates.txt -S $TEMPFILES.sam
    ;;
esac

# Sorting and reducing the SAM-file to q>10

echo -e "Sorting and reducing $SAMPLENAME to q>10."\\n


exe $SAMTOOLS view -h -S -q 10 $TEMPFILES.sam -@ $SAMTHREADS -o $TEMPFILES.red.sam 

exe $SAMTOOLS sort -O sam -@ $SAMTHREADS -T $TEMPFILES -o $TEMPFILES.sorted.sam $TEMPFILES.red.sam 


# Shifting the read positions for ATAC

echo -e "Shifting the read positions."\\n

exe myATACshiftSAM.pl "$TEMPFILES.sorted.sam" "$TEMPFILES.shifted" "$TYPE"


echo -e "Generating TagDir, BigWig, and TDF and finding Peaks."\\n

# Generating a HOMER tag directory 

exe makeTagDirectory $TAGDIRDIR/$SAMPLENAME $TEMPFILES.shifted.sam -genome $HOMERGENOME -checkGC

# Making a UCSC browser track

exe makeUCSCfile $TAGDIRDIR/$SAMPLENAME -bigWig $CHROMSIZES -fragLength 73 -o $BIGWIGDIR/$SAMPLENAME.bigwig

# Finding peaks

exe findPeaks $TAGDIRDIR/$SAMPLENAME -region -size 150 -o auto

if [ "$SEQTYPE" == 2 ]; then
exe $PICARD CollectInsertSizeMetrics \
      I="$TEMPFILES.shifted.sam" \
      O="$TAGDIRDIR/$SAMPLENAME/insert_size_metrics.txt" \
      H="$TAGDIRDIR/$SAMPLENAME/insert_size_histogram.pdf" \
      W=500
fi

exe $SAMTOOLS view -bS -@ $SAMTHREADS -o "${MAPPINGDIR}/${SAMPLENAME}.shifted.bam" $TEMPFILES.shifted.sam 

# removing tmp files

rm ${TEMPFILES}*

echo -e \\n"Done."\\n


