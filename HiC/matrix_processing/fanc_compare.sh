#!/bin/bash
# fanc_compare.sh

### Get contrast Hi-C matrices with fanc compare: KD vs CTRL / MUTANT vs CTRL

# FAN-C version: 0.9.25
# contrast: difference
# .hic object resolution: 10kb

WDIR='/private/'

cd $WDIR

# Submit trimming jobs using a custom script for Sun Grid Engine's qsub command
# 'qb' submits jobs with Grid Engine's qsub command allocating 1 core and and giving a custom jobid
# Usage: qb <jobid> <command>

# 1. CD34 model

# SA2 KD vs Ctrl
qb fanc.compare.diff.stag2_vs_siCtrl \
	fanc compare \
		-c difference \
		fanc/hic/binned/merged_CD34_stag2/merged_CD34_stag2_10kb.hic \
		fanc/hic/binned/merged_CD34_siCtrl/merged_CD34_siCtrl_10kb.hic \
		fanc/hic/binned/merged_CD34_siCtrl/merged_CD34_stag2_vs_siCtrl_10kb.diff.hic
# SA1 KD vs Ctrl
qb fanc.compare.diff.stag1_vs_siCtrl \
	fanc compare \
		-c difference \
		fanc/hic/binned/merged_CD34_stag1/merged_CD34_stag1_10kb.hic \
		fanc/hic/binned/merged_CD34_siCtrl/merged_CD34_siCtrl_10kb.hic \
		fanc/hic/binned/merged_CD34_siCtrl/merged_CD34_stag1_vs_siCtrl_10kb.diff.hic
# RAD21 KD vs Ctrl
qb fanc.compare.diff.rad21_vs_siCtrl \
	fanc compare \
		-c difference \
		fanc/hic/binned/merged_CD34_rad21/merged_CD34_rad21_10kb.hic \
		fanc/hic/binned/merged_CD34_siCtrl/merged_CD34_siCtrl_10kb.hic \
		fanc/hic/binned/merged_CD34_siCtrl/merged_CD34_rad21_vs_siCtrl_10kb.diff.hic

# 2. AML patients

#SA2 mut vs Ctrl
qb fanc.compare.diff.stag2mut_vs_Ctrl \
	fanc compare \
		-c difference \
		fanc/hic/binned/merged_sa2_mut/merged_sa2_mut.10kb.hic \
		fanc/hic/binned/merged_aml_ctrl/merged_aml_ctrl.10kb.hic \
		fanc/hic/binned/merged_aml_ctrl/merged_sa2_mut_vs_siCtrl_10kb.diff.hic
#RAD21 mut vs Ctrl
qb fanc.compare.diff.rad21mut_vs_Ctrl \
	fanc compare \
		-c difference \
		fanc/hic/binned/merged_rad21_mut/merged_rad21_mut.10kb.hic \
		fanc/hic/binned/merged_aml_ctrl/merged_aml_ctrl.10kb.hic \
		fanc/hic/binned/merged_aml_ctrl/merged_rad21_mut_vs_siCtrl_10kb.diff.hic
