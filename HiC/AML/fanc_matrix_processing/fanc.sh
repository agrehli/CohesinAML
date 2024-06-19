#!/bin/bash
# fanc.sh

### Running fanc auto on trimmed FASTQ files to obtain processed Hi-C matrices

# FAN-C version: 0.9.25
# BWA version: 0.7.17
# Mapping quality cutoff q: 3
# Restriction enzyme DpnII
# Genome hg38

WDIR='/private/'
INPUTDIR='fastq_trimmed/'
SAMPLES=(ctr_2236 ctr_1551 ctr_1747 ctr_3323 ctr_3488 ctr_6246 ctr_5285 RAD21mut_UKR186 RAD21mut_7314 RAD21mut_12557 RAD21mut_41580a SA2mut_12514 SA2mut_12567 SA2mut_24603 SA2mut_2193 siCtrl_14 siCtrl_17 siCtrl_18 siCtrl_20 siCtrl_21 siCtrl_22 siCtrl_27 siCtrl_28 RAD21KD_18 RAD21KD_20 RAD21KD_22 RAD21KD_27 RAD21KD_28 SA1KD_14 SA1KD_17 SA1KD_20 SA1KD_21 SA1KD_27 SA1KD_28 SA2KD_14 SA2KD_17 SA2KD_20 SA2KD_21 SA2KD_22 SA2KD_28)

# the hg38 genome was downloaded from UCSC using
# wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/*'
# the hg38 genome fasta file was modified to have the order established in the negspy repository 
# in order to allow visualisation with HiGlass: https://docs.higlass.io/data_preparation.html
# hg38 negspy coordinates https://github.com/pkerpedjiev/negspy/blob/master/negspy/data/hg38/chromInfo.txt

GENOME='/private/Hsap/UCSC/hg38/hg38.negspy.fa'
GENOMEIDX='/private/Hsap/UCSC/hg38/BWAIndexNegspy/hg38.negspy.fa'

# Submit trimming jobs using a custom script for Sun Grid Engine's qsub command
# 'q32b' submits jobs with Grid Engine's qsub command allocating 32 cores and giving a custom jobid
# Usage: q32b <jobid> <command>

cd $WDIR
for i in $(seq 1 ${#SAMPLES[*]}); do \
	OUTPUTDIR=hic_BWA/$SAMPLES[${i}]/
	q32b fanc.$SAMPLES[${i}] \
		fanc.$SAMPLES[${i}-2],trim_galore.$SAMPLES[${i}] \
			fanc auto \
				$INPUTDIR/$SAMPLES[${i}]_1_val_1.fq.gz \
				$INPUTDIR/$SAMPLES[${i}]_2_val_2.fq.gz \
				$OUTPUTDIR \
				--genome $GENOME \
				--restriction-enzyme DpnII \
				--genome-index $GENOMEIDX \
				--basename $SAMPLES[${i}] \
				-q 3 \
				--threads 32 \
				-tmp --split-ligation-junction -f ;
done
	