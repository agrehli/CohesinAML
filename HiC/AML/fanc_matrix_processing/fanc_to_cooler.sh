#!/bin/bash
# fanc_to_cooler.sh

### Convert Hi-C matrices from FAN-C's '.hic' format into multi-resolution cooler files '.cool'

# FAN-C version: 0.9.25
# cooler version: 0.8.9
# initial resolution for .hic object: 10kb

WDIR='/private/'
SAMPLES=(ctr_2236 ctr_1551 ctr_1747 ctr_3323 ctr_3488 ctr_6246 ctr_5285 RAD21mut_UKR186 RAD21mut_7314 RAD21mut_12557 RAD21mut_41580a SA2mut_12514 SA2mut_12567 SA2mut_24603 SA2mut_2193 siCtrl_14 siCtrl_17 siCtrl_18 siCtrl_20 siCtrl_21 siCtrl_22 siCtrl_27 siCtrl_28 RAD21KD_18 RAD21KD_20 RAD21KD_22 RAD21KD_27 RAD21KD_28 SA1KD_14 SA1KD_17 SA1KD_20 SA1KD_21 SA1KD_27 SA1KD_28 SA2KD_14 SA2KD_17 SA2KD_20 SA2KD_21 SA2KD_22 SA2KD_28)

# Submit trimming jobs using a custom script for Sun Grid Engine's qsub command
# 'q8c' submits jobs with Grid Engine's qsub command allocating 8 cores and waiting for a prior job to finish
# Usage: q8c <jobid> <jobtowaitfor1,jobtowaitfor2> <command>

cd $WDIR
for i in $(seq 1 ${#SAMPLES[*]}); do \
	INPUTDIR=hic_BWA/$SAMPLES[${i}]/hic/binned/
	q8c to_cooler.10kb.$SAMPLES[${i}] \
		to_cooler.10kb.$SAMPLES[${i}-4],fanc.$SAMPLES[${i}] \
		fanc to-cooler \
			$INPUTDIR/$SAMPLES[${i}]_10kb.hic \
			$INPUTDIR/$SAMPLES[${i}]_10kb.cool \
			--threads 8 -tmp;
done

SAMPLES=(merged_sa2_mut merged_rad21_mut merged_aml_ctrl merged_CD34_rad21 merged_CD34_siCtrl merged_CD34_stag1 merged_CD34_stag2)
cd $WDIR
for i in $(seq 1 ${#SAMPLES[*]}); do \
	INPUTDIR=fanc/hic/binned/$SAMPLES[${i}]
	q8c to_cooler.10kb.$SAMPLES[${i}] \
		to_cooler.10kb.$SAMPLES[${i}-4],fanc.$SAMPLES[${i}] \
		fanc to-cooler \
			$INPUTDIR/$SAMPLES[${i}]_10kb.hic \
			$INPUTDIR/$SAMPLES[${i}]_10kb.cool \
			--threads 8 -tmp;
done
