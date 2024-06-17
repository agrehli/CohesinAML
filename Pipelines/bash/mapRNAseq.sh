#!/bin/bash
# by Michael Rehli, Febr 23, 2018
# bash script to map RNA-seq data, generate a count table, and coverage bigWigs
# will only keep unique.aligned.bam, bigwig and tdf files as well as gene count data
# IF YOU REQUIRE MULTIMAPPERS, OR ANY OTHER ORIGINAL FILES, PLEASE RUN STAR MANUALLY

# modified by MR, March16, 2018 to incorporate genome GRCm38

# modified by MR, November9, 2018 to incorporate alternative strand orientation

#setting homer environment

DIR_PKG="/misc/software/ngs"

PATH_PERL=/misc/software/package/perl/perl-5.26.1/bin
PATH_SAMTOOLS=${DIR_PKG}/samtools/samtools-1.6/bin
PATH_HOMER=${DIR_PKG}/homer/v4.11/bin
PATH_R=/misc/software/package/RBioC/3.4.3/bin

export PATH=${PATH_R}:${PATH_PERL}:${PATH_SAMTOOLS}:${PATH_HOMER}:${PATH}
export PATH

# Defining the program versions (needs to be adjusted to the server/workstation used)
# rhskl17 version

STAR="/misc/software/ngs/STAR/v2.5.3a/STAR/source/STAR"
IGVTOOLS="/usr/bin/java -Xmx4g -Djava.awt.headless=true -jar /misc/software/viewer/IGV/IGVTools_2.3.98/igvtools.jar"
INDEXSDIR="/misc/software/ngs/genome/index/STAR"


# required input:
# 
# -f <fastq or fastq.gz> forward read 
# -r <fastq or fastq.gz> reverse read (if not available, or set "NA" will run as single end)
# -g <genome>  available genomes include hg19, GRCh38, GRCm38
# -l <read length>
# -o <output subdirectory> only subdirectories above RNA-folder are given
# -n <sample name> used for sam, bigwig and tag directory
#
# optional:
#
# -t <number of threads>
# -c <read length> trim reads to the given length
# -s revert strands for wiggle output (e.g. for Illumina's stranded Tru-Seq)


# Set Script Name variable
SCRIPT=`basename ${BASH_SOURCE[0]}`

# Set fonts for Help.
NORM=`tput sgr0`
BOLD=`tput bold`
REV=`tput smso`

# USAGE function
function USAGE {
  echo -e \\n"Help documentation for ${BOLD}${SCRIPT}.${NORM}"\\n
  echo -e "${REV}Basic usage:${NORM}\n${BOLD}$SCRIPT -g <genome> -f <fastq> -r <fastq> -o <subdirectory> -n <name> -t <threads> -c <trim length>${NORM}"\\n
  echo "Required:"
  echo "${BOLD}-g${NORM} <genome>  available genomes include hg19, GRCh38, GRCm38"
  echo "${BOLD}-f${NORM} <fastq or fastq.gz> forward read"
  echo "${BOLD}-r${NORM} <fastq or fastq.gz> reverse read (if not available, or set \"NA\" will run as single end)"
  echo "${BOLD}-d${NORM} <output subdirectory> only subdirectories above RNA-folder are given"
  echo "${BOLD}-n${NORM} <sample name> used for sam, bigwig and tag directory"
  echo "${BOLD}-t${NORM} <number of threads> (optional/default 12)"
  echo "${BOLD}-l${NORM} <read length> used for index"
  echo "${BOLD}-s${NORM} revert strands for wiggle output (e.g. for Illumina's stranded Tru-Seq)"
  echo "${BOLD}-c${NORM} <length> trim reads to the given length"
  echo "Available indices:"
  for entry in $INDEXSDIR/*
	do
    echo $entry
  done
 
  
  
  echo -e \\n"Example: ${BOLD}$SCRIPT -g hg19 -l 50 -f read1.fastq.gz -r read2.fastq.gz -d test -n testrun -t 24${NORM} -c 43"\\n
  exit 1
}

#Check the number of arguments. If none are passed, print help and exit.
NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
  USAGE
fi

# Defining available genomes
GENOMES=(hg19 GRCh38 hg38 GRCm38)

# Set default number of threads to 16
THREADS=12
TRIM=0
REVERT=0
COMP=0

# Parse command line options
while getopts :g:f:r:d:n:t:l:c:sh OPTIONS; do
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
    l)  #set option "l"
      LENGTH=$OPTARG
      ;;
    c)  #set option "c"
      TRIMLENGTH=$OPTARG
      TRIM=1
      ;;
    h)  #show help
      USAGE
      ;;
    s)  #set option "s"
      REVERT=1
      ;;
    \?) #unrecognized option - show USAGE
      echo -e \\n"Option -${BOLD}$OPTARG${NORM} not allowed."
      USAGE
      ;;
  esac
done
shift $((OPTIND-1)) 

                                  #################
                                  # Sanity checks #
                                  #################


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
BIGWIGGENOME=$GENOME

# path to chromosome sizes 
CHROMSIZES="/misc/software/viewer/IGV/IGVTools_2.3.98/genomes/$GENOME.chrom.sizes"
# STAR index
INDEX="$INDEXSDIR/index_${GENOME}_${LENGTH}"

if [ $GENOME == "hg38" ] || [ $GENOME == "GRCh38" ] ; then
      GENOME="hg38"
      INDEX="$INDEXSDIR/index_GRCh38.PRI.p10_$LENGTH"
      BIGWIGGENOME="GRCh38"
      CHROMSIZES="/misc/software/viewer/IGV/IGVTools_2.3.98/genomes/GRCh38.PRI_p10.chrom.sizes"
fi

if [ $GENOME == "mm10" ] || [ $GENOME == "GRCm38" ] ; then
      GENOME="GRCm38"
      INDEX="$INDEXSDIR/index_GRCm38.PRI.p5_$LENGTH"
      BIGWIGGENOME="GRCm38"
      CHROMSIZES="/misc/software/viewer/IGV/IGVTools_2.3.98/genomes/GRCm38.PRI_p5.chrom.sizes"
fi


#check if index is available
if [ ! -d "$INDEX" ]; then
	echo "Index is missing!"
    echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
    exit 7
fi
 
# Check whether one or two reads are given, whether files are available, non-identical and have the right extensions

if [ ! -f "$FASTQR1" ]; then
    echo "Fastq file (read 1) not found!"
    echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
    exit 3
fi
case "$FASTQR1" in
  *.fastq.gz) ZIP="--readFilesCommand zcat " 
              COMP=1
              ;;
  *.fastq)    ZIP=""
              ;;
  *.fq.gz)    ZIP="--readFilesCommand zcat "
              COMP=1
              ;;
  *.fq)       ZIP=""
              ;;
  *)          echo "Fastq file (read 1) has wrong extension (should be fastq, fq, fastq.gz or fq.gz!"
              echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
              exit 4
              ;;
esac

if [ -z "$FASTQR2" ] || [ "$FASTQR2" == "NA" ] || [ "$FASTQR2" == "na" ] ; then
    SEQTYPE=1
    TYPE="SE"
    SKEWTYPE="any"
    echo -e \\n"Data processing for ${BOLD}single-end${NORM} reads"
elif [ ! -f "$FASTQR2" ]; then
    echo "Fastq file (read 2) not found!"
    echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
    exit 5   
elif [ "$FASTQR1" == "$FASTQR2" ]; then
    echo "Fastq files for read 1 and 2 are identical!"
    echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
    exit 6
else
	SEQTYPE=2
	TYPE="PE"
	SKEWTYPE="pe"
	echo -e \\n"Data Processing for ${BOLD}paired-end${NORM} reads"
	case "$FASTQR2" in
       *.fastq.gz) ZIP="--readFilesCommand zcat " 
       ;;
       *.fastq)    ZIP=""
       ;;
       *.fq.gz)    ZIP="--readFilesCommand zcat "
       ;;
       *.fq)       ZIP=""
       ;;
       *)          echo "Fastq file (read 2) has wrong extension (should be fastq, fq, fastq.gz or fq.gz!"
       	           echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
                   exit 7
       ;;
    esac
fi

# Check whether the subfolders are already there - if not ask whether path should be made
 
if [ -z "$DIRECTORY" ]; then
	echo "Subdirectory is missing!"
    echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
    exit 7
fi

# rhskl17 version
MAPPINGDIR="/misc/data/processedData/mapping/RNA/${BIGWIGGENOME}/RNAseq/${DIRECTORY}"
BIGWIGDIR="/misc/data/processedData/bigWig/RNA/${BIGWIGGENOME}/RNAseq/${DIRECTORY}"

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

# Creating the folder for STAR

if [ ! -d "${MAPPINGDIR}/${SAMPLENAME}" ]; then			#looks for the mapping output folder - if not available, askes whether to make it
   mkdir "${MAPPINGDIR}/${SAMPLENAME}"
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


# Optionally starting the pipeline by trimming the reads with homerTools
# needs unzipped files

if [ $TRIM == 1 && $COMP == 0 ]; then
   ZIP=""
   case "$SEQTYPE" in
   1 )  exe homerTools trim -len $TRIMLENGTH $FASTQR1
        exe mv "$FASTQR1.trimmed" "$TEMPFILES.trimmed.fastq" 
        mv "$FASTQR1.lengths" "${MAPPINGDIR}/logs/${SAMPLENAME}.trimming.log"
        FASTQR1="$TEMPFILES.trimmed.fastq"
        ;;
   2 )  exe homerTools trim -len $TRIMLENGTH $FASTQR1
        exe mv "$FASTQR1.trimmed" "$TEMPFILES.trimmed-pair1.fastq" 
        mv "$FASTQR1.lengths" "${MAPPINGDIR}/logs/${SAMPLENAME}.R1.trimming.log"     
        FASTQR1="$TEMPFILES.trimmed-pair1.fastq"
        exe homerTools trim -len $TRIMLENGTH $FASTQR2
        exe mv "$FASTQR2.trimmed" "$TEMPFILES.trimmed-pair2.fastq" 
        mv "$FASTQR2.lengths" "${MAPPINGDIR}/logs/${SAMPLENAME}.R2.trimming.log"
        FASTQR2="$TEMPFILES.trimmed-pair2.fastq"
        ;;        
   esac     
elif [ $TRIM == 1 ] && [ $COMP == 1 ]; then
   ZIP=""
   case "$SEQTYPE" in
   1 )  exe gunzip -c $FASTQR1 > $TEMPFILES.fastq
        exe homerTools trim -len $TRIMLENGTH $TEMPFILES.fastq
        FASTQR1="$TEMPFILES.fastq.trimmed"
        mv "$TEMPFILES.fastq.lengths" "${MAPPINGDIR}/logs/${SAMPLENAME}.trimming.log"
        ;;
   2 )  exe gunzip -c $FASTQR1 > $TEMPFILES.R1.fastq
        exe homerTools trim -len $TRIMLENGTH $TEMPFILES.R1.fastq
        FASTQR1="$TEMPFILES.R1.fastq.trimmed"
        mv "$TEMPFILES.R1.fastq.lengths" "${MAPPINGDIR}/logs/${SAMPLENAME}.R1.trimming.log"     
        exe gunzip -c $FASTQR2 > $TEMPFILES.R2.fastq
        exe homerTools trim -len $TRIMLENGTH $TEMPFILES.R2.fastq
        FASTQR1="$TEMPFILES.R1.fastq.trimmed"
        mv "$TEMPFILES.R2.fastq.lengths" "${MAPPINGDIR}/logs/${SAMPLENAME}.R2.trimming.log"
        ;;        
   esac     
fi


# Running STAR

echo -e \\n"Running STAR on $SAMPLENAME."\\n

cd "${MAPPINGDIR}/${SAMPLENAME}"

case "$SEQTYPE" in
	1 )  exe $STAR --runThreadN $THREADS \
	--genomeDir $INDEX \
	--readFilesIn $FASTQR1 \
	$ZIP\
	--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
	--outReadsUnmapped Fastq \
	--alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 \
	--alignMatesGapMax 1000000 \
	--alignIntronMax 1000000 \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode GeneCounts \
	--outWigType bedGraph \
	--outWigStrand Stranded
	;;
	2 )  exe $STAR --runThreadN $THREADS \
	--genomeDir $INDEX \
	--readFilesIn $FASTQR1 $FASTQR2 \
	$ZIP\
	--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
	--outReadsUnmapped Fastq \
	--alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 \
	--alignMatesGapMax 1000000 \
	--alignIntronMax 1000000 \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode GeneCounts \
	--outWigType bedGraph \
	--outWigStrand Stranded
	;;
esac


# Sorting bedgraph files

echo -e \\n"Sorting the bedGraph files for $SAMPLENAME."\\n

LC_COLLATE=C sort -k1,1 -k2,2n $MAPPINGDIR/$SAMPLENAME/Signal.Unique.str1.out.bg -o $TEMPFILES.str1.sorted.bg 
LC_COLLATE=C sort -k1,1 -k2,2n $MAPPINGDIR/$SAMPLENAME/Signal.Unique.str2.out.bg -o $TEMPFILES.str2.sorted.bg 

echo -e "\$ LC_COLLATE=C sort -k1,1 -k2,2n $MAPPINGDIR/$SAMPLENAME/Signal.Unique.str1.out.bg -o $TEMPFILES.str1.sorted.bg"\\n >> $LOGFILE
echo -e "\$ LC_COLLATE=C sort -k1,1 -k2,2n $MAPPINGDIR/$SAMPLENAME/Signal.Unique.str2.out.bg -o $TEMPFILES.str2.sorted.bg"\\n >> $LOGFILE


# inverting values for the second strand
case "$REVERT" in
     0)  exe awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4*(-1)}' $TEMPFILES.str2.sorted.bg > $TEMPFILES.str2.inv.bg
         mv $TEMPFILES.str2.inv.bg $TEMPFILES.str2.sorted.bg
         ;;
     1)  exe awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4*(-1)}' $TEMPFILES.str1.sorted.bg > $TEMPFILES.str1.inv.bg
         mv $TEMPFILES.str1.inv.bg $TEMPFILES.str1.sorted.bg
	     ;;
esac
        
         
# generating bigWigs for both strands

echo -e \\n"Generating bigWig files for $SAMPLENAME."\\n

exe bedGraphToBigWig $TEMPFILES.str1.sorted.bg $CHROMSIZES $BIGWIGDIR/$SAMPLENAME.Unique.for.bigwig
exe bedGraphToBigWig $TEMPFILES.str2.sorted.bg $CHROMSIZES $BIGWIGDIR/$SAMPLENAME.Unique.rev.bigwig

# generating corresponding tdf files for IGV browser

echo -e "Generating TDF files for $SAMPLENAME."\\n

exe bigWigToWig $BIGWIGDIR/$SAMPLENAME.Unique.for.bigwig $TEMPFILES.str1.sorted.wig -udcDir=/loctmp
exe bigWigToWig $BIGWIGDIR/$SAMPLENAME.Unique.rev.bigwig $TEMPFILES.str2.sorted.wig -udcDir=/loctmp

exe $IGVTOOLS toTDF $TEMPFILES.str1.sorted.wig $BIGWIGDIR/$SAMPLENAME.Unique.for.cov.tdf $CHROMSIZES
exe $IGVTOOLS toTDF $TEMPFILES.str2.sorted.wig $BIGWIGDIR/$SAMPLENAME.Unique.rev.cov.tdf $CHROMSIZES

# move and rename gene count table and bam file

case "$TRIM" in
	1 )  exe mv ${MAPPINGDIR}/${SAMPLENAME}/Aligned.sortedByCoord.out.bam ${MAPPINGDIR}/${SAMPLENAME}.trimmed.sorted.bam
         exe mv ${MAPPINGDIR}/${SAMPLENAME}/ReadsPerGene.out.tab ${MAPPINGDIR}/${SAMPLENAME}.trimmed.ReadsPerGene.txt
         exe mv ${MAPPINGDIR}/${SAMPLENAME}/Log.final.out ${MAPPINGDIR}/logs/${SAMPLENAME}.trimmed.summary.log.txt
         exe mv ${MAPPINGDIR}/${SAMPLENAME}/Log.out ${MAPPINGDIR}/logs/${SAMPLENAME}.trimmed.processing.log.txt
         ;;
    0 )  exe mv ${MAPPINGDIR}/${SAMPLENAME}/Aligned.sortedByCoord.out.bam ${MAPPINGDIR}/${SAMPLENAME}.sorted.bam
         exe mv ${MAPPINGDIR}/${SAMPLENAME}/ReadsPerGene.out.tab ${MAPPINGDIR}/${SAMPLENAME}.ReadsPerGene.txt
         exe mv ${MAPPINGDIR}/${SAMPLENAME}/Log.final.out ${MAPPINGDIR}/logs/${SAMPLENAME}.summary.log.txt
         exe mv ${MAPPINGDIR}/${SAMPLENAME}/Log.out ${MAPPINGDIR}/logs/${SAMPLENAME}.processing.log.txt
         ;;
esac


# removing tmp files

exe rm -r "${MAPPINGDIR}/${SAMPLENAME}"

exe rm ${TEMPFILES}.*



