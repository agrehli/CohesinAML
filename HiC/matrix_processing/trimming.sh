#!/bin/bash
# trimming.sh

### Running Trim Galore on FASTQ files to filter out illumina adaptors

# Trim Galore version: 0.6.7
# Cutadapt version: 3.4
# Quality Phred score cutoff: 30
# Quality encoding type selected: ASCII+33

WDIR='/private/'
INPUTDIR='fastq/'
OUTPUTDIR='fastq_trimmed/'
SAMPLES=(ctr_2236 ctr_1551 ctr_1747 ctr_3323 ctr_3488 ctr_6246 ctr_5285 RAD21mut_UKR186 RAD21mut_7314 RAD21mut_12557 RAD21mut_41580a SA2mut_12514 SA2mut_12567 SA2mut_24603 SA2mut_2193 siCtrl_14 siCtrl_17 siCtrl_18 siCtrl_20 siCtrl_21 siCtrl_22 siCtrl_27 siCtrl_28 RAD21KD_18 RAD21KD_20 RAD21KD_22 RAD21KD_27 RAD21KD_28 SA1KD_14 SA1KD_17 SA1KD_20 SA1KD_21 SA1KD_27 SA1KD_28 SA2KD_14 SA2KD_17 SA2KD_20 SA2KD_21 SA2KD_22 SA2KD_28)

# Submit trimming jobs using a custom script for Sun Grid Engine's qsub command
# 'q16c' submits jobs with Grid Engine's qsub command allocating 16 cores and waiting for a prior job to finish
# Usage: q16c <jobid> <jobtowaitfor1,jobtowaitfor2> <command>

cd $WDIR
for i in $(seq 1 ${#SAMPLES[*]}); do \
	q16c trim_galore.$SAMPLES[${i}] \
		trim_galore.$SAMPLES[${i}-2],wget.$SAMPLES[${i}] \
			trim_galore -j 16 \
			--paired --illumina --fastqc \
			--output_dir $OUTPUTDIR \
			$INPUTDIR/$SAMPLES[${i}]_R1.fastq.gz \
			$INPUTDIR/$SAMPLES[${i}]_R2.fastq.gz ;
done