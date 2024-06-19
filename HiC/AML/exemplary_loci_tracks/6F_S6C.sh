#!/bin/bash
# 6F_S6C.sh

### Plot 2 Mb window around ILDR2 and in an example locus using fancplot. 
### Control matrices in the top right corner and Mutant/Patient matrices in the bottom left.

# FAN-C version: 0.9.25
# .hic object resolution: 10kb
# coordinates: 
# 6F: ILDR2 : chr1:165975540-167975539
# S6C: example : chr5:156000000-160000000
# view: 'split'

WDIR='/private/'
cd $WDIR
OUTDIR=plots/Fischer2023/FoldChange
mkdir -p $OUTDIR

# Plotting parameters
DATE='2038_01_19' #23_03_08,23_02_21
VMAXS=(0.015 0.02) # maximum saturation for Hi-C plots

# ILDR2 coordinates, and example coordinates
NAME=(ILDR2_2Mb S6C)
CHR=(chr1 chr1)
START=(165975540 156000000)
END=(167975539 160000000)

COMPARISONS=(merged_sa2_mut_vs_siCtrl merged_rad21_mut_vs_siCtrl merged_CD34_stag2_vs_siCtrl merged_CD34_stag1_vs_siCtrl merged_CD34_rad21_vs_siCtrl)
MATRIXA=(merged_aml_ctrl merged_aml_ctrl merged_CD34_siCtrl merged_CD34_siCtrl merged_CD34_siCtrl)
MATRIXB=(merged_sa2_mut merged_rad21_mut merged_CD34_stag2 merged_CD34_stag1 merged_CD34_rad21)
LABELA=(CTRL-AML CTRL-AML CTRL-HSPCs CTRL-HSPCs CTRL-HSPCs)
LABELB=(STAG2-mut RAD21-mut STAG2-KD STAG1-KD RAD21-KD)

# Submit trimming jobs using a custom script for Sun Grid Engine's qsub command
# 'qc' submits jobs with Grid Engine's qsub command allocating 1 core and waiting for a prior job to finish
# Usage: qc <jobid> <jobtowaitfor1,jobtowaitfor2> <command>

# split-view
for i in $(seq 1 ${#NAME[*]}); do \
	for l in $(seq 1 ${#COMPARISONS[*]}); do \
		for m in $(seq 1 ${#VMAXS[*]}); do \
		qc fancplot.fanc.compare.fc.split.6F.$COMPARISONS[${l}].$NAME[${i}].thredholds.VMAX.$VMAXS[${m}].pdf \
			fancplot.fanc.compare.fc.split.6F.$COMPARISONS[${l}].$NAME[${i}].thredholds.VMAX.$VMAXS[${m}-1].pdf \
			fancplot $CHR[${i}]:$START[${i}]-$END[${i}] \
				--width 2 -o $OUTDIR/${DATE}.$NAME[${i}].comp.split.6F.$COMPARISONS[${l}]_10kb.VMAX.$VMAXS[${m}].pdf --pdf-text-as-font --tick-locations $START[${i}] $END[${i}] \
				-p split fanc/hic/binned/$MATRIXA[${l}]/$MATRIXA[${l}]_10kb.hic fanc/hic/binned/$MATRIXB[${l}]/$MATRIXB[${l}]_10kb.hic -vmin 0 -vmax $VMAXS[${m}] --hide-major-ticks --title $LABELB[${l}]_under_$LABELA[${l}] \
				-p gene -s 8 --group-by gene_name --squash -cf dodgerblue -cr slategrey gencode.v33.annotation.edit.gtf;
		done;
	done;
done
