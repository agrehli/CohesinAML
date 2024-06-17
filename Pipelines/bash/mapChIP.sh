#!/bin/bash
# by Michael Rehli, Febr 23, 2018, mod. March 5, 2018, mod. March 11, 2018, mod. March 26, 2016
# bash script to map ChIPseq-like data (SE- only), generate tagDirectories, coverage bigWigs and TDF files 
# optional for input data: copy number alterations 
# 

# June 12, 2023 changed script to allow adding multiple fastq files




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

BOWTIE="${DIR_PKG}/bowtie/bowtie2-2.3.4-linux-x86_64/bowtie2"
SAMTOOLS="${DIR_PKG}/samtools/samtools-1.16.1/bin/samtools"
BEDTOOLS="${DIR_PKG}/bedtools/bedtools2-2.27.1/bin/bedtools"
IGVTOOLS="/usr/bin/java -Xmx4g -Djava.awt.headless=true -jar /misc/software/viewer/IGV/IGVTools_2.3.98/igvtools.jar"

# required input:
# 
# -f <fastq or fastq.gz>
# -g <genome>  available genomes include hg19, GRCh38, mm9, mm10
# -o <output subdirectory> only subdirectories above ChIP- or Input-folder are given
# -n <sample name> used for sam, bigwig and tag directory
#
# optional:
# -i indicates that the sample is input
# -m indicates that the sample is MCIp
# -t <number of threads>

# options specific for input samples:
# -c <read length for mappability> generates copy number profiles and annotates CNVs with Cancer Genes from Cosmic Database (Available mappability data for 50, 75, and 100 bp)
# -s <XX/XY> gender XX female, XY male (default) - only necessary if -c is given
# -a <tissue> E (epithelial), L (leukaemia/lymphoma), M (mesenchymal), O (other) will output Cancer Genes only for given tissues (only human!)

# Set Script Name variable
SCRIPT=`basename ${BASH_SOURCE[0]}`

# Set fonts for Help.
NORM=`tput sgr0`
BOLD=`tput bold`
REV=`tput smso`

# USAGE function
function USAGE {
  echo -e \\n"Help documentation for ${BOLD}${SCRIPT}.${NORM}"\\n
  echo -e "${REV}Basic usage:${NORM}\n${BOLD}$SCRIPT -g <genome> -f <fastq> -o <subdirectory> -n <name> -t <threads>${NORM} "\\n
  # -c <read length> -s <gender>"\\n
  echo "Required:"
  echo "${BOLD}-g${NORM} <genome>  available genomes include hg19, GRCh38, mm9, mm10"
  echo "${BOLD}-f${NORM} <fastq or fastq.gz> forward read (several files can be given, separated by komma w/o space)"
  echo "${BOLD}-d${NORM} <output subdirectory> only subdirectories above ChIP- or Input-folder are given"
  echo "${BOLD}-n${NORM} <sample name> used for bam, bigwig and tag directory, etc."
  echo -e \\n"Optional:"
  echo "${BOLD}-t${NORM} <number of threads> (optional/default 16)"
  echo "${BOLD}-i${NORM} treat sample as genomic input (stored in DNA directory)"
  echo "${BOLD}-m${NORM} treat sample as MCIp (stored in DNA directory, TagDir with -keepOne)"
  echo "${BOLD}-c${NORM} <read length> will use FreeC to generate CNV profiles of Input files. Read length refers to the mappability files available for a given genome,"
  echo "   which are found in /misc/software/ngs/genome/sequence/XXXX (XXXX = genome)"
  echo "${BOLD}-s${NORM} <gender> accepts either XX or XY (default) - only necessary if -c is given"
  echo "${BOLD}-a${NORM} <tissue> E (epithelial), L (leukaemia/lymphoma), M (mesenchymal), O (other) will output COSMIC Cancer Genes for given tissues ${BOLD}HUMAN ONLY!${NORM}"
  echo -e \\n"Example: ${BOLD}$SCRIPT -g hg19 -f TFread.fastq.gz -d MOMACDC -n TF -t 32 ${NORM}"\\n
  exit 1
}

#Check the number of arguments. If none are passed, print help and exit.
NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
  USAGE
fi

# Defining available genomes
GENOMES=(hg19 GRCh38 hg38 mm9 mm10)
TISSUETYPES=(E L M O)
KEEP=""

# Set defaults
GENDER="XY"
REF="ref"
THREADS=16
RUNFREEC=0
ANNOTATE=0
INPUT=0
MCIP=0
FOLDER="ChIP"
TOPFOLDER="chromatin"

# Parse command line options
while getopts :g:f:d:n:t:c:s:a:imh OPTIONS; do
  case $OPTIONS in
    g)  #set option "g" - Bowtie will will only accept GRCh38 while HOMER needs hg38
      GENOME=$OPTARG
      HOMERGENOME=$OPTARG
      ;;
    f)  #set option "f"
      FASTQR1=$OPTARG
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
    c)  #set option "c"
      MAPPAB=$OPTARG
      RUNFREEC=1
      ;;
    s)  #set option "s"
      GENDER=$OPTARG
      ;;
    h)  #show help
      USAGE
      ;;
    i)  #set option "i"
      echo -e \\n"Sample is genomic input."\\n
      INPUT=1
      FOLDER="Input"
      TOPFOLDER="DNA"
      ;;
    m)  #set option "m"
      echo -e \\n"Sample is MCIp."\\n
      MCIP=1
      FOLDER="MCIp"
      TOPFOLDER="DNA"
      KEEP="-keepOne "
      ;;
    a)  #set option "a"
	  ANNOTATE=1
	  TISSUE=$OPTARG
      ;;

    \?) #unrecognized option - show USAGE
      echo -e \\n"Option -${BOLD}$OPTARG${NORM} not allowed."
      USAGE
      ;;
  esac
done
shift $((OPTIND-1)) 

# Sanity checks

# excude using -m and -i together
if [ $INPUT == 1 ] && [ $MCIP == 1 ]; then
      echo -e \\n"${BOLD}Use either -i or -m.{NORM}"
      echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
      exit 2
fi

# Check whether Tissues are selected and correct

if [ $ANNOTATE == 1 ] && [ $INPUT == 1 ]; then

	INTISSUE=$(echo ${TISSUETYPES[@]} | grep -o "$TISSUE" | wc -w)

	if [ $INTISSUE == 0 ] ; then
      echo -e \\n"Tissue -${BOLD}$TISSUE${NORM} not available (Options include E, M, L & O)."
      echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
      exit 2
	fi
fi



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

# path to GEM files (mappability info), and chromosome sizes for FreeC
# needs two different one's due to a bug in FreeC
GEMPATH="/misc/software/ngs/genome/sequence/${GENOME}/mappability"
GEM="${GEMPATH}/${GENOME}_${MAPPAB}.mappability"
case $GENDER in
	XY)   REF="ref"
	;;
	XX)	  REF="ref-Y"
	;;
esac
REDCHROMSIZES="/misc/software/ngs/freec/FREEC-11.0/data/${GENOME}.${REF}.chrom.sizes"

if [ $GENOME == "hg38" ] || [ $GENOME == "GRCh38" ] ; then
      GENOME="hg38"
      GENOMEPATH="GRCh38.PRI_p10"
      BIGWIGGENOME="GRCh38"
      HOMERGENOME="hg38"
      CHROMSIZES="/misc/software/viewer/IGV/IGVTools_2.3.98/genomes/GRCh38.PRI_p10.chrom.sizes"
      GEMPATH="/misc/software/ngs/genome/sequence/GRCh38.PRI_p10/mappability"
      GEM="${GEMPATH}/GRCh38_${MAPPAB}.mappability"
      REDCHROMSIZES="/misc/software/ngs/freec/FREEC-11.0/data/${GENOME}.${REF}.chrom.sizes"
fi

#bowtie index
INDEX="$GENOMEPATH/$GENOME"

# if -c is given, check whether mappability track is available

if [ "$RUNFREEC" == 1 ] && [ ! -f "$GEM" ]; then
    echo "Gem file for ${GENOME} not found!"
    echo -e "Available gem files for ${GENOME} are:"
    for entry in $GEMPATH/*
	   do
       echo $entry
    done
    echo -e \\n"Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
    exit 3
fi
 
# Check whether fastq file is available and has the right extension

## First check if more than one fastq is provided (added June 12th, 2023 by MR)

if [[ "$FASTQR1" =~ "," ]]; then
    echo "... multiple Fastqs given"
    R1ARRAY=($(echo $FASTQR1 | tr "," "\n"))
    for R1 in "${R1ARRAY[@]}"; do
      if [ ! -f "$R1" ]; then
        echo "Fastq file $R1 (read 1) not found!"
       echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
       exit 3
       fi
    done
elif [ ! -f "$FASTQR1" ]; then
    echo "Fastq file (read 1) not found!"
    echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
    exit 3
fi
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

# Check whether the subfolders are given and already there - if not ask whether path should be made
 
if [ -z "$DIRECTORY" ]; then
	echo "Subdirectory is missing!"
    echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
    exit 7
fi

# rhskl17 version
MAPPINGDIR="/misc/data/processedData/mapping/${TOPFOLDER}/${BIGWIGGENOME}/${FOLDER}/${DIRECTORY}"
TAGDIRDIR="/misc/data/processedData/tagDir/${TOPFOLDER}/${BIGWIGGENOME}/${FOLDER}/${DIRECTORY}"
BIGWIGDIR="/misc/data/processedData/bigWig/${TOPFOLDER}/${BIGWIGGENOME}/${FOLDER}/${DIRECTORY}"

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
	
if [ ! -d "$TAGDIRDIR" ]; then			#looks for the tagDir output folder - if not available, asks whether to make it
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


# Running Bowtie2

echo -e \\n"Running Bowtie2 on $SAMPLENAME."\\n

exe $BOWTIE --very-sensitive -p $THREADS -x $INDEX --met-file "$MAPPINGDIR/logs/$SAMPLENAME.aln_metrics.txt" -U $FASTQR1 -S $TEMPFILES.sam 2> "$MAPPINGDIR/logs/$SAMPLENAME.aln_rates.txt"

####################################################
#modified on March 26, 2018
#MCIp should not be filtered for mapping quality

case "$MCIP" in 
	0 ) 	# Sorting and reducing the SAM-file to q>10
		echo -e "Sorting and reducing $SAMPLENAME to q>10."\\n
		exe $SAMTOOLS view -h -S -q 10 $TEMPFILES.sam -@ $SAMTHREADS -o $TEMPFILES.red.sam 
		exe $SAMTOOLS sort -O sam -@ $SAMTHREADS -T $TEMPFILES -o $TEMPFILES.sorted.sam $TEMPFILES.red.sam 
	;;
	1 )  # Sorting
		echo -e "Sorting $SAMPLENAME."\\n
		exe $SAMTOOLS sort -O sam -@ $SAMTHREADS -T $TEMPFILES -o $TEMPFILES.sorted.sam $TEMPFILES.sam 
	;;
esac
####################################################


echo -e "Generating TagDir, BigWig, and TDF and finding Peaks."\\n

# Generating a HOMER tag directory 

exe makeTagDirectory $TAGDIRDIR/$SAMPLENAME $TEMPFILES.sorted.sam $KEEP -genome $HOMERGENOME -checkGC

# Making a UCSC browser track

exe makeUCSCfile $TAGDIRDIR/$SAMPLENAME -bigWig $CHROMSIZES -o $BIGWIGDIR/$SAMPLENAME.bigwig

# Making an IGV browser track

#exe $IGVTOOLS count $TEMPFILES.sorted.sam $BIGWIGDIR/$SAMPLENAME.cov.tdf $CHROMSIZES


# Run FreeC if requested

if [ $RUNFREEC == 1 ]; then
    if [ ! -d "${MAPPINGDIR}/CNVdata" ]; then			    #looks for the CNVdata output folder - if not available, makes it
       mkdir "${MAPPINGDIR}/CNVdata"
    fi
    if [ ! -d "${MAPPINGDIR}/CNVprofiles" ]; then			#looks for the CNVprofile output folder - if not available, makes it
       mkdir "${MAPPINGDIR}/CNVprofiles"
    fi

    exe mv $TEMPFILES.sorted.sam "${MAPPINGDIR}/${SAMPLENAME}.sam"
    
    case "$GENOME" in 
  		hg19 ) 	COSMIC="/misc/software/ngs/freec/FREEC-11.0/data/COSMIC_Cancer_Gene_Census_GRCh37v84_2018.tsv"
		;;
  		hg38 )  COSMIC="/misc/software/ngs/freec/FREEC-11.0/data/COSMIC_Cancer_Gene_Census_GRCh38v84_2018.tsv"
		;;
	esac

    exe runFreeC.pl ${MAPPINGDIR}/CNVdata ${GENDER} ${GENOME} ${MAPPINGDIR} ${SAMPLENAME} ${REDCHROMSIZES} ${GEM}
	
	# annotation of CNVs with COSMIC genes
	exe FreeC2bed.pl ${MAPPINGDIR}/CNVdata/${SAMPLENAME}.sam_CNVs ${TEMPFILES}.CNV.bed
	exe mv ${MAPPINGDIR}/CNVdata/${SAMPLENAME}.pdf ${MAPPINGDIR}/CNVprofiles
	
	if [ $GENOME == "hg38" ] || [ $GENOME == "GRCh38" ] || [ $GENOME == "hg19" ]; then
		case "$ANNOTATE" in 
			0 )  echo -e "\$ ConvertCOSMIC.pl ${COSMIC} >${TEMPFILES}.cosmic.bed"\\n  >> $LOGFILE
				 ConvertCOSMIC.pl ${COSMIC} > ${TEMPFILES}.cosmic.bed
				 echo -e "\$ ${BEDTOOLS} intersect -a ${TEMPFILES}.CNV.bed -b ${TEMPFILES}.cosmic.bed -wa -wb >${TEMPFILES}.inter.bed"\\n >> $LOGFILE
				 ${BEDTOOLS} intersect -a ${TEMPFILES}.CNV.bed -b ${TEMPFILES}.cosmic.bed -wa -wb >${TEMPFILES}.inter.bed
				 exe reformatCOSMICannotation.pl ${TEMPFILES}.inter.bed ${MAPPINGDIR}/CNVdata/${SAMPLENAME}.annotated_CNVs.txt
			;;
			1 )  echo -e "\$ ConvertCOSMIC.pl ${COSMIC} -tissues ${TISSUE} >${TEMPFILES}.cosmic.bed"\\n >> $LOGFILE
				 ConvertCOSMIC.pl ${COSMIC} -tissues ${TISSUE} > ${TEMPFILES}.cosmic.bed
				 echo -e "\$ ${BEDTOOLS} intersect -a ${TEMPFILES}.CNV.bed -b ${TEMPFILES}.cosmic.bed -wa -wb >${TEMPFILES}.inter.bed"\\n >> $LOGFILE
				 ${BEDTOOLS} intersect -a ${TEMPFILES}.CNV.bed -b ${TEMPFILES}.cosmic.bed -wa -wb >${TEMPFILES}.inter.bed
				 exe reformatCOSMICannotation.pl ${TEMPFILES}.inter.bed ${MAPPINGDIR}/CNVdata/${SAMPLENAME}.${TISSUE}.annotated_CNVs.txt
			;;
		esac
	fi
    #convert sam to bam
    exe $SAMTOOLS view -bS -@ $SAMTHREADS -o "${MAPPINGDIR}/${SAMPLENAME}.bam" "${MAPPINGDIR}/${SAMPLENAME}.sam" 
    exe rm "${MAPPINGDIR}/${SAMPLENAME}.sam" 

else

   # just convert sam to bam
   exe $SAMTOOLS view -bS -@ $SAMTHREADS -o "${MAPPINGDIR}/${SAMPLENAME}.bam" $TEMPFILES.sorted.sam 

fi



# removing tmp files

rm ${TEMPFILES}.*



