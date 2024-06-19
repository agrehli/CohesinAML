#!/bin/bash
# S6B.sh

### Plot 60 Mb window using fancplot. 
### Control matrices in the top right corner and Mutant/Patient matrices in the bottom left.

# FAN-C version: 0.9.25
# .hic object resolution: 50kb
# coordinates: chr1:160000000-220000000
# view: 'split'

WDIR='/private/'
cd $WDIR
OUTDIR=plots/Fischer2023/FoldChange
mkdir -p $OUTDIR

# Plotting parameters
DATE='2038_01_19' #23_03_09
VMAXS=(0.001) # maximum saturation for Hi-C plots
RESS=(50kb) #matrix resolution

# coordinates : 60 Mb
NAME=(S6B)
CHR=(chr1)
START=(160000000)
END=(220000000)

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
			for n in $(seq 1 ${#RESS[*]}); do \
				qc fancplot.fanc.compare.fc.split.$COMPARISONS[${l}].$NAME[${i}].thredholds.VMAX.$VMAXS[${m}].RES.$RESS[${n}].pdf \
					fanc.bin.AMLSAMPLES.$MATRIXA[${l}].RESOLUTION.$RESS[${n}],fanc.bin.AMLSAMPLES.$MATRIXB[${l}].RESOLUTION.$RESS[${n}],fancplot.fanc.compare.fc.split.$COMPARISONS[${l}].$NAME[${i}].thredholds.VMAX.$VMAXS[${m}-1].RES.$RESS[${n}].pdf \
					fancplot $CHR[${i}]:$START[${i}]-$END[${i}] \
						--width 2 -o $OUTDIR/${DATE}.$NAME[${i}].comp.split.$COMPARISONS[${l}]_$RESS[${n}].diff.$DIFF[${j}].VMAX.$VMAXS[${m}].pdf --pdf-text-as-font --tick-locations $START[${i}] $END[${i}] \
						-p split fanc/hic/binned/$MATRIXA[${l}]/$MATRIXA[${l}]_$RESS[${n}].hic fanc/hic/binned/$MATRIXB[${l}]/$MATRIXB[${l}]_$RESS[${n}].hic -vmin 0 -vmax $VMAXS[${m}] --hide-major-ticks --title $LABELB[${l}]_under_$LABELA[${l}];
			done;
		done;
	done;
done