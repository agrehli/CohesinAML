#!/bin/bash

#--------------------------------------------------------------------------------------------------
# Name        : DNAMappingBWA.sh
# Author      : Remco Hoogenboezem
# Version     :
# Copyright   :
# Description : DNAMappingBWA.sh -r referenceSequence -1 R1.fastq [-2 R2.fastq][-o outputPrefix][-t threads][-R readGroup][-h]
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#Get input arguments
#--------------------------------------------------------------------------------------------------

display_usage() {
	echo -e 'Usage: DNAMappingBWA.sh -r referenceSequence -1 R1.fastq [-2 R2.fastq][-o outputPrefix][-t threads][-R readGroup][-b bedFile][-h]\n'
	echo -e 'Optional arguments:'
	echo -e '-t\t\tNumber of threads to use in the pipeline. Please note that not all the steps are parallelizable [1]';
	echo -e '-b\t\tPath to the BED file with the target regions of NGS, which will be extracted from the BAM file [NULL]'
	echo -e '-R\t\tName of the set of reads that were generated from a single run of a sequencing instrument [NULL]'
 }


THREADS='1'
PICARD='/tools/picard-tools/picard-2.17.10.jar'

while getopts ':r:1:2:o:t:R:c:h' ARG
do
	case $ARG in

		r)

			REFERENCE_SEQUENCE=$OPTARG

		;;

		1)

			R1_FASTQ=$OPTARG

		;;

		2)
			
			R2_FASTQ=$OPTARG

		;;

		o)
	
			OUTPUT_PREFIX=$OPTARG

		;;

		t)

			THREADS=$OPTARG
		;;

		R)

			READ_GROUP=$OPTARG

		;;		

		b)

			BED_FILE=$OPTARG

		;;	

		h)
			
			echo 'Usage: DNAMappingBWA.sh -r referenceSequence -1 R1.fastq [-2 R2.fastq][-o outputPrefix][-t threads][-R readGroup][-h]'
			exit 0			

		;;			

		\?)

			echo "Error: Invalid option -$OPTARG (-h for help)"
			exit -1
		;;

	esac

done

#--------------------------------------------------------------------------------------------------
#Check input arguments
#--------------------------------------------------------------------------------------------------

if [[ -z $REFERENCE_SEQUENCE ]]
then

	echo 'Error: Please specify a reference sequence (-h for help)'
	exit -1

fi

if [[ -z $R1_FASTQ ]]
then

	echo 'Error: Please specify a fastq file (-h for help)'
	exit -1

fi

if [[ ! -f $R1_FASTQ ]]
then

	echo "Warning: Specified fastq file \"$R1_FASTQ\" is not a valid file or does not exist"
#	exit -1

fi

if [[ ! -z $R2_FASTQ ]]
then

	if [[ ! -f $R2_FASTQ ]]
	then

		echo "Warning: Specified fastq file \"$R2_FASTQ\" is not a valid file or does not exist"
	#	exit -1

	fi

fi

if [[ -z $OUTPUT_PREFIX ]]
then

	OUTPUT_PREFIX=`echo $R1_FASTQ | sed 's:\(_R[12]\)*\.fastq\(.gz\)\{0,1\}$::'`

fi

#--------------------------------------------------------------------------------------------------
#Extract readgroup from fastq file (assume illumina fastq file!)
#--------------------------------------------------------------------------------------------------

if [[ -z $READ_GROUP ]]
then

	SAMPLE=`echo $R1_FASTQ | grep -o '[^/]*$' | sed 's:\(_R[12]\)*\.fastq\(.gz\)\{0,1\}$::'`

	if [[ -z $SAMPLE ]]
	then

		echo 'Error: Empty sample name after removing extension'
		exit -1

	fi

	TOKENS=(`zcat -f $R1_FASTQ | head -n1 | sed 's;:;\n;g'`)

	if [[ ${#TOKENS[@]} -lt 4 ]]
	then

		echo "Error: Invalid illumina sequence identifier line (expected at least 4 colon seperated fields): $R1_FASTQ"
#exit -1

	fi

	FLOWCELL=${TOKENS[2]}
	LANE=${TOKENS[3]}

	ID=$RANDOM
	SM=$SAMPLE
	PL='ILLUMINA'
	LB="LIB-$SAMPLE"
	PU="${FLOWCELL}.${LANE}"
	
	READ_GROUP="@RG\tID:$ID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU"

fi

#--------------------------------------------------------------------------------------------------
#Create statistics folder
#--------------------------------------------------------------------------------------------------

STATISTICS_FOLDER=`echo $OUTPUT_PREFIX | sed 's:[^/]*$:mapping_statistics:'`

mkdir -p $STATISTICS_FOLDER

if [[ $? -ne 0 ]]
then

	echo "Error: Could not create statistics folder: $STATISTICS_FOLDER"
	exit -1

fi

STATISTICS_PREFIX=`echo $OUTPUT_PREFIX | sed 's:\([^/]*$\):mapping_statistics/\1:'`

#--------------------------------------------------------------------------------------------------
#Alignment using bwa
#--------------------------------------------------------------------------------------------------

echo 'Info: run bwa mem'

if [[ ! -e ${OUTPUT_PREFIX}_unsorted.bam && ! -e ${OUTPUT_PREFIX}.bam ]]; then

    bwa mem -t $THREADS -R $READ_GROUP -M $REFERENCE_SEQUENCE $R1_FASTQ $R2_FASTQ | samtools view -b - 2> $STATISTICS_PREFIX.bwa_log > ${OUTPUT_PREFIX}_unsorted.bam

    if [[ $? -ne 0 ]]
    then

	    echo 'Error: bwa mem failed!'
	    exit -1

    fi

else

    echo "File ${OUTPUT_PREFIX}_unsorted.bam found. Skipping..."

fi

echo ''

#--------------------------------------------------------------------------------------------------
#run SortSam
#--------------------------------------------------------------------------------------------------

if [[ ! -e ${OUTPUT_PREFIX}.bam || $(samtools view -H ${OUTPUT_PREFIX}.bam | grep 'SO:coordinate') == "" ]]; then

    echo 'Info: run sambamba sort'

    sambamba sort -m 8GB -t $THREADS -o ${OUTPUT_PREFIX}.bam ${OUTPUT_PREFIX}_unsorted.bam

    if [ $? -ne 0 ]
    then

	    echo 'Error: sambamba sort failed!'
	    exit -1

    fi

    rm ${OUTPUT_PREFIX}_unsorted.bam

else 

    echo 'Info: BAM file already sorted. Skipping...'

fi

#--------------------------------------------------------------------------------------------------
#run MarkDuplicates
#--------------------------------------------------------------------------------------------------

if [[ ! -e ${OUTPUT_PREFIX}.bam || $(samtools view -H ${OUTPUT_PREFIX}.bam | grep 'ID:MarkDuplicates') == "" ]]; then

    echo 'Info: run MarkDuplicates'

    FLOWCELL=$(samtools view ${OUTPUT_PREFIX}.bam | head -1 | cut -f1 | awk -F':' '$3~/H[A-Z0-9]{4}DMXX$/{print "PATTERNED";exit}{print "UNPATTERNED"}')
    echo -e "Flowcell type is $FLOWCELL"

    if [[ $FLOWCELL == "PATTERNED" ]]; then
        PIXEL_DISTANCE=2500 # Recommended by https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates
    else
        PIXEL_DISTANCE=100
    fi
    
    java -Xmx24g -jar $PICARD MarkDuplicates INPUT=${OUTPUT_PREFIX}.bam OUTPUT=${OUTPUT_PREFIX}_markdup.bam METRICS_FILE=${STATISTICS_PREFIX}.duplicate_metrics REMOVE_DUPLICATES=false ASSUME_SORTED=true CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT OPTICAL_DUPLICATE_PIXEL_DISTANCE=$PIXEL_DISTANCE

    if [[ $? -ne 0 ]]
    then
	    
	    echo 'Error: MarkDuplicates failed!'
	    exit -1

    fi

    mv ${OUTPUT_PREFIX}_markdup.bam ${OUTPUT_PREFIX}.bam
    mv ${OUTPUT_PREFIX}_markdup.bai ${OUTPUT_PREFIX}.bam.bai # bam.bai is widely accepted

else

    echo 'Info: Duplicates found to be marked. Skipping...'

fi

echo ''

#--------------------------------------------------------------------------------------------------
#run CollectMultipleMetrics
#--------------------------------------------------------------------------------------------------

echo 'Info: run CollectMultipleMetrics'

java -Xmx16g -jar $PICARD CollectMultipleMetrics INPUT=${OUTPUT_PREFIX}.bam OUTPUT=$STATISTICS_PREFIX PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution PROGRAM=MeanQualityByCycle PROGRAM=CollectBaseDistributionByCycle PROGRAM=CollectGcBiasMetrics PROGRAM=CollectSequencingArtifactMetrics PROGRAM=CollectQualityYieldMetrics REFERENCE_SEQUENCE=$REFERENCE_SEQUENCE ASSUME_SORTED=true

if [ $? -ne 0 ]
then

	echo 'Error: CollectMultipleMetrics failed!'
	exit -1

fi

echo ''

#--------------------------------------------------------------------------------------------------
#Remove phix sequence (necessary for some downstream analysis)
#--------------------------------------------------------------------------------------------------

phix_reads=$(samtools idxstats ${OUTPUT_PREFIX}.bam | grep -i phix | awk 'BEGIN{sum=0}NR>0{sum=sum+$3+$4}END{print(sum)}')
total_mapping=$(samtools idxstats ${OUTPUT_PREFIX}.bam | awk '{sum = sum + $3 + $4}END{print(sum-($3+$4))}')
total_reads=$(samtools idxstats ${OUTPUT_PREFIX}.bam | awk '{sum = sum + $3 + $4}END{print(sum)}')
echo -e "${OUTPUT_PREFIX}\t${phix_reads}\t${total_mapping}\t${total_reads}\t$(echo $phix_reads $total_mapping $total_reads | awk -v OFS='\t' '{print($1*100/$2,$1*100/$3)}')" > ${STATISTICS_PREFIX}.phix

samtools view -h ${OUTPUT_PREFIX}.bam | grep -v phix | grep -v "decoy" | samtools view -b > ${OUTPUT_PREFIX}_nophix.bam

mv ${OUTPUT_PREFIX}_nophix.bam ${OUTPUT_PREFIX}.bam
sambamba index ${OUTPUT_PREFIX}.bam


# We remove phix, which is used as a sort of decoy to capture contamination| grep -v "phix" -- we do it after quality controls because 
# otherwise the program fail due to missing sequences in the BAM file (ie phix). Alternatively, reference genome without phix can be used downstream.

#--------------------------------------------------------------------------------------------------
#Coverage statistics
#--------------------------------------------------------------------------------------------------

mkdir -p coverage_statistics

if [[ ! -z $BED_FILE ]]; then
       /tools/depth_of_coverage/depth_of_coverage -a ${OUTPUT_PREFIX}.bam -e $BED_FILE -o coverage_statistics/${OUTPUT_PREFIX##*/}
else
       /tools/depth_of_coverage/depth_of_coverage -a ${OUTPUT_PREFIX}.bam -o coverage_statistics/${OUTPUT_PREFIX##*/}
fi

if [ $? -ne 0 ]
then

	echo 'Error: depth_of_coverage failed!'
	exit -1

fi

#--------------------------------------------------------------------------------------------------
#Generate bigwig (normalized?)
#--------------------------------------------------------------------------------------------------

mkdir -p bigwig
/tools/pcap_core/bin/bam2bw -i ${OUTPUT_PREFIX}.bam -o bigwig/${OUTPUT_PREFIX}.bigwig # Twice as fast as BamToBigwig.sh

if [ $? -ne 0 ]
then

	echo 'Error: bam2bw failed!'
	exit -1

fi
#--------------------------------------------------------------------------------------------------
#done
#--------------------------------------------------------------------------------------------------

echo -e "sample\tphix_reads\tmapping_reads\ttotal_reads\tphix/mapping\tphix/total" > summary_phix_counts.txt
cat $STATISTICS_FOLDER/*.phix >> summary_phix_counts.txt

exit 0

