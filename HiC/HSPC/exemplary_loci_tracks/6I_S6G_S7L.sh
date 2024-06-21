#!/bin/bash
# 6I_S6G_S7L.sh

### Plot 1.5 Mb window around SOCS2, SH2D2A, PAWR, ITGA9 and DACT1 using fancplot. 
### Control matrices on the top row, Mutant/Patient matrices in middle row, and 'difference' matrices in the bottom row.

# FAN-C version: 0.9.25
# .hic object resolution: 10kb
# coordinates: 
# 6I: SOCS2 : chr12:92855000-94355000
# S6G: SH2D2A : chr1:156066848-157566847
# S6G: PAWR : chr12:78940964-80440963
# S7L: ITGA9 : chr3:36702115-38202114
# S7L: DACT1 : chr14:57883967-59383966
# view: 'triangular'

WDIR='/private/'
cd $WDIR
OUTDIR=plots/Fischer2023/FoldChange
mkdir -p $OUTDIR

# Plotting parameters
DATE='2038_01_19' #23_03_10
DIFF=(0.005) # lower (-#) and maximum (+#) saturation for 'difference 'Hi-C plots
COLS=(RdBu_r) # colorscale for 'difference 'Hi-C plots
VMAXS=(0.015 0.02) # maximum saturation for Hi-C plots
MAX_DIST=(700kb) # cap of maximum distance away from diagonal to plot

# 1.5 Mb coordinates
NAME=(SOCS2_15Mb SH2D2A_15Mb PAWR_15Mb ITGA9_15Mb DACT1_15Mb)
CHR=(chr12 chr1 chr12 chr3 chr14)
START=(92855000 156066848 78940964 36702115 57883967)
END=(94355000 157566847 80440963 38202114 59383966)

COMPARISONS=(merged_sa2_mut_vs_siCtrl merged_rad21_mut_vs_siCtrl merged_CD34_stag2_vs_siCtrl merged_CD34_stag1_vs_siCtrl merged_CD34_rad21_vs_siCtrl)
MATRIXA=(merged_aml_ctrl merged_aml_ctrl merged_CD34_siCtrl merged_CD34_siCtrl merged_CD34_siCtrl)
MATRIXB=(merged_sa2_mut merged_rad21_mut merged_CD34_stag2 merged_CD34_stag1 merged_CD34_rad21)
LABELA=(CTRL-AML CTRL-AML CTRL-HSPCs CTRL-HSPCs CTRL-HSPCs)
LABELB=(STAG2-mut RAD21-mut STAG2-KD STAG1-KD RAD21-KD)

# Submit trimming jobs using a custom script for Sun Grid Engine's qsub command
# 'qc' submits jobs with Grid Engine's qsub command allocating 1 core and waiting for a prior job to finish
# Usage: qc <jobid> <jobtowaitfor1,jobtowaitfor2> <command>

# triangular + difference style plot; duplicating each matrix
for i in $(seq 1 ${#NAME[*]}); do \
	for j in $(seq 1 ${#DIFF[*]}); do \
		for k in $(seq 1 ${#COLS[*]}); do \
			for l in $(seq 1 ${#COMPARISONS[*]}); do \
				for m in $(seq 1 ${#VMAXS[*]}); do \
					qc fancplot.fanc.compare.fc.triangular_plus_diff.$COMPARISONS[${l}].$NAME[${i}].thredholds.col.$COLS[${k}].VMAX.$VMAXS[${m}].pdf \
						fancplot.fanc.compare.fc.triangular_plus_diff.$COMPARISONS[${l}].$NAME[${i}].thredholds.col.$COLS[${k}].VMAX.$VMAXS[${m}-1].pdf \
						fancplot $CHR[${i}]:$START[${i}]-$END[${i}] \
							--width 2 -o $OUTDIR/${DATE}.$NAME[${i}].comp.triangular_plus_diff.$COMPARISONS[${l}]_10kb.diff.$DIFF[${j}].col.$COLS[${k}].VMAX.$VMAXS[${m}].pdf --tick-locations $START[${i}] $END[${i}] \
							-p triangular --maximum-distance $MAX_DIST \
									-vmin 0 -vmax $VMAXS[${m}] --hide-x --hide-major-ticks --title $LABELA[${l}] \
									fanc/hic/binned/$MATRIXA[${l}]/$MATRIXA[${l}]_10kb.hic \
							-p triangular --maximum-distance $MAX_DIST \
									-vmin 0 -vmax $VMAXS[${m}] --hide-x --hide-major-ticks --title $LABELB[${l}] \
									fanc/hic/binned/$MATRIXB[${l}]/$MATRIXB[${l}]_10kb.hic \
							-p triangular --maximum-distance $MAX_DIST \
									-vmin -$DIFF[${j}] -vmax $DIFF[${j}] -c $COLS[${k}] --hide-major-ticks --title "Difference" \
									fanc/hic/binned/$MATRIXA[${l}]/$COMPARISONS[${l}]_10kb.diff.hic \
							-p gene -s 8 --group-by gene_name --squash -cf dodgerblue -cr slategrey gencode.v33.annotation.edit.gtf;
				done;
			done;
		done;
	done;
done