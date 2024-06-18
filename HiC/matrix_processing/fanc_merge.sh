#!/bin/bash
# fanc_merge.sh

### Merge Hi-C matrices of biological replicates

# FAN-C version: 0.9.25
# chromosomes included: chr1...chr25,chrX,chrY,chrM
# normalisation: Knight-Ruiz
# Filter bins with a coverage lower than the 0.1% of the median contact count


### Merging of hic objects

WDIR='/private/'

# Submit trimming jobs using a custom script for Sun Grid Engine's qsub command
# 'qb' submits jobs with Grid Engine's qsub command allocating 1 core and and giving a custom jobid
# Usage: qb <jobid> <command>

cd $WDIR
OUTPUTDIR=fanc/hic/merged_rad21_mut/
mkdir -p $OUTPUTDIR
qb fanc_merge.merged_rad21_mut \
	fanc hic \
		hic_BWA/$SAMPLES[${i}]/hic/RAD21mut_UKR186.hic \
		hic_BWA/$SAMPLES[${i}]/hic/RAD21mut_7314.hic \
		hic_BWA/$SAMPLES[${i}]/hic/RAD21mut_12557.hic \
		hic_BWA/$SAMPLES[${i}]/hic/RAD21mut_41580a.hic \
		$OUTPUTDIR/merged_rad21_mut.hic

OUTPUTDIR=fanc/hic/merged_sa2_mut/
mkdir -p $OUTPUTDIR
qb fanc_merge.merged_sa2_mut \
	fanc hic \
		hic_BWA/$SAMPLES[${i}]/hic/SA2mut_12514.hic \
		hic_BWA/$SAMPLES[${i}]/hic/SA2mut_12567.hic \
		hic_BWA/$SAMPLES[${i}]/hic/SA2mut_24603.hic \
		hic_BWA/$SAMPLES[${i}]/hic/SA2mut_2193.hic \
		$OUTPUTDIR/merged_sa2_mut.hic

OUTPUTDIR=fanc/hic/merged_aml_ctrl/
mkdir -p $OUTPUTDIR
qb fanc_merge.merged_aml_ctrl \
	fanc hic \
		hic_BWA/$SAMPLES[${i}]/hic/ctr_2236.hic \
		hic_BWA/$SAMPLES[${i}]/hic/ctr_1551.hic \
		hic_BWA/$SAMPLES[${i}]/hic/ctr_1747.hic \
		hic_BWA/$SAMPLES[${i}]/hic/ctr_3323.hic \
		hic_BWA/$SAMPLES[${i}]/hic/ctr_3488.hic \
		hic_BWA/$SAMPLES[${i}]/hic/ctr_6246.hic \
		hic_BWA/$SAMPLES[${i}]/hic/ctr_5285.hic \
		$OUTPUTDIR/merged_aml_ctrl.hic

OUTPUTDIR=fanc/hic/merged_CD34_rad21/
mkdir -p $OUTPUTDIR
qb fanc_merge.merged_CD34_rad21 \
	fanc hic \
		hic_BWA/$SAMPLES[${i}]/hic/RAD21KD_18.hic \
		hic_BWA/$SAMPLES[${i}]/hic/RAD21KD_20.hic \
		hic_BWA/$SAMPLES[${i}]/hic/RAD21KD_22.hic \
		hic_BWA/$SAMPLES[${i}]/hic/RAD21KD_27.hic \
		hic_BWA/$SAMPLES[${i}]/hic/RAD21KD_28.hic \
		$OUTPUTDIR/merged_CD34_rad21.hic

OUTPUTDIR=fanc/hic/merged_CD34_siCtrl/
mkdir -p $OUTPUTDIR
qb fanc_merge.merged_CD34_siCtrl \
	fanc hic \
		hic_BWA/$SAMPLES[${i}]/hic/siCtrl_14.hic \
		hic_BWA/$SAMPLES[${i}]/hic/siCtrl_17.hic \
		hic_BWA/$SAMPLES[${i}]/hic/siCtrl_18.hic \
		hic_BWA/$SAMPLES[${i}]/hic/siCtrl_20.hic \
		hic_BWA/$SAMPLES[${i}]/hic/siCtrl_21.hic \
		hic_BWA/$SAMPLES[${i}]/hic/siCtrl_22.hic \
		hic_BWA/$SAMPLES[${i}]/hic/siCtrl_27.hic \
		hic_BWA/$SAMPLES[${i}]/hic/siCtrl_28.hic \
		$OUTPUTDIR/merged_CD34_siCtrl.hic

OUTPUTDIR=fanc/hic/merged_CD34_stag1/
mkdir -p $OUTPUTDIR
qb fanc_merge.merged_CD34_stag1 \
	fanc hic \
		hic_BWA/$SAMPLES[${i}]/hic/SA1KD_14.hic \
		hic_BWA/$SAMPLES[${i}]/hic/SA1KD_17.hic \
		hic_BWA/$SAMPLES[${i}]/hic/SA1KD_20.hic \
		hic_BWA/$SAMPLES[${i}]/hic/SA1KD_21.hic \
		hic_BWA/$SAMPLES[${i}]/hic/SA1KD_27.hic \
		hic_BWA/$SAMPLES[${i}]/hic/SA1KD_28.hic \
		$OUTPUTDIR/merged_CD34_stag1.hic

OUTPUTDIR=fanc/hic/merged_CD34_stag2/
mkdir -p $OUTPUTDIR
qb fanc_merge.merged_CD34_stag2 \
	fanc hic \
		hic_BWA/$SAMPLES[${i}]/hic/SA2KD_14.hic \
		hic_BWA/$SAMPLES[${i}]/hic/SA2KD_17.hic \
		hic_BWA/$SAMPLES[${i}]/hic/SA2KD_20.hic \
		hic_BWA/$SAMPLES[${i}]/hic/SA2KD_21.hic \
		hic_BWA/$SAMPLES[${i}]/hic/SA2KD_22.hic \
		hic_BWA/$SAMPLES[${i}]/hic/SA2KD_28.hic \
		$OUTPUTDIR/merged_CD34_stag2.hic


### Binning merged '.hic' files into specific resolutions

SAMPLES=(merged_sa2_mut merged_rad21_mut merged_aml_ctrl merged_CD34_rad21 merged_CD34_siCtrl merged_CD34_stag1 merged_CD34_stag2)
RESOLUTIONS=(1mb 500kb 100kb 50kb 10kb)

# Submit trimming jobs using a custom script for Sun Grid Engine's qsub command
# 'q8c' submits jobs with Grid Engine's qsub command allocating 8 cores and waiting for a prior job to finish
# Usage: q8c <jobid> <jobtowaitfor1,jobtowaitfor2> <command>

cd $WDIR
for i in $(seq 1 ${#SAMPLES[*]}); do \
	mkdir -p fanc/hic/binned/$SAMPLES[${i}]/
	mkdir -p hic_BWA/$SAMPLES[${i}]/plots/stats/
	for j in ${RESOLUTIONS[*]}; do \
		q8c fanc.bin.${j}.SAMPLES.$SAMPLES[${i}] \
			fanc.bin.${j}.SAMPLES.$SAMPLES[${i}-4],fanc_merge.$SAMPLES[${i}] \
				fanc hic \
				hic_BWA/$SAMPLES[${i}]/hic/$SAMPLES[${i}].hic \
				fanc/hic/binned/$SAMPLES[${i}]/$SAMPLES[${i}]_${j}.hic \
				-b ${j} -f \
				--filter-low-coverage-relative 0.1 \
				--normalise --norm-method KR \
				--statistics hic_BWA/$SAMPLES[${i}]/plots/stats/$SAMPLES[${i}].bincorrectfilter.stats.${j}.txt \
				--statistics-plot hic_BWA/$SAMPLES[${i}]/plots/stats/$SAMPLES[${i}].bincorrectfilter.stats.${j}.pdf \
				--chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM \
				--threads 8 \
				-tmp ;
    done ;
done