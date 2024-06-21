#!/bin/bash
# S7B_coolpuppy.sh

### Perform Aggregate Peak Analysis (APA) with coolpup.py on indivdual patient or replicate Hi-C matrices of Ctrl and Cohesin deficient cells/patients.
### Enrichment is calculated against 100 shifted windows as controls

### Hi-C signal was aggregated over loop anchors overlapping 'weakened' or 'strengthened' loops
### between patients with STAG2 mutations compared to control patients
### Besides being differential loops, their anchors or 'edges' have to overlap Cohesin-associated enhancers
### The underlying coordinates were defined in 7B_cohesin_e_p_diff_loops_size_split.R
### The Hi-C matrices in cooler format were produced in ../fanc_matrix_processing/fanc_to_cooler.sh

# cooler version: 0.9.25
# .hic object resolution: 10kb
# excluded chromosomes: chrX,chrY,chr17_GL000205v2_random

# APA on merged Cohesin patients and merged CD34 Hi-C matrices at 10kb resolution
WDIR='/private/'
cd $WDIR
OUTDIR=$WDIR/APA/DifferentialLoopAnalysis/
mkdir -p $OUTDIR

COMPARISONS=(SA2mutvsCTRL)
EDGES=(both minus plus either)
DIRECTIONS=(strenghtened weakened)
SIZEGROUPS=(All 1 2 3 4 5 6 7)
TYPE=(CohesinAssEnhancersOverlap H3K27acOverlap)
AMLSAMPLES=(RAD21_UKR186_Rep1_R1 AML_SA2_9708_Rep2_S56_L002_R1 AML_SA2_27396_Rep1_R1 AML_ctr_18136_Rep1_R1 AML_ctr_19405_Rep1_R1 AML_ctr_21047_Rep1_R1 AML_SA2_29728_R1 AML_SA2_24743_R1 AML_RAD21_38455_R1 AML_RAD21_26830_R1 AML_RAD21_23039_R1 AML_ctr_21290_R1 AML_ctr_19416_R1 AML_ctr_18519_R1 AML_ctr_16911_R1 CD34_28_5_SA2_KD_R1 CD34_27_3_SA1_KD_R1 CD34_22_2_SA2_KD_R1 CD34_21_4_siCtrl_Rep1_R1 CD34_21_3_SA2_KD_R1 CD34_21_2_SA1_KD_R1 CD34_20_6_siCtrl_Rep1_R1 CD34_20_5_SA2_KD_R1 CD34_20_4_SA1_KD_R1 CD34_20_1_RAD21_KD_R1 CD34_18_4_siCtrl_R1 CD34_18_1_RAD21_KD_R1 CD34_17_2_SA2_KD_R1 CD34_17_1_SA1_KD_R1 CD34_14_3_siCtrl_R1 CD34_14_2_SA2_KD_R1 CD34_28_6_siCtrl_R1 CD34_28_4_SA1_KD_R1 CD34_28_4_SA1_KD_R1 CD34_27_4_siCtrl_R1 CD34_17_3_siCtrl_R1 CD34_27_1_RAD21_KD_R1 CD34_14_1_SA1_KD_R1 CD34_22_3_siCtrl_R1 CD34_28_1_RAD21_KD_R1 CD34_22_1_RAD21_KD_R1)
RESOLUTIONS=(10kb)

# Submit trimming jobs using a custom script for Sun Grid Engine's qsub command
# 'q8b' submits jobs with Grid Engine's qsub command allocating 8 cores and giving a custom jobid
# Usage: q8b <jobid> <command>

for o in ${TYPE[*]}; do \
	for k in ${COMPARISONS[*]}; do \
		for l in ${EDGES[*]}; do \
			for m in ${DIRECTIONS[*]}; do \
				for n in ${SIZEGROUPS[*]}; do \
					for i in $(seq 1 ${#AMLSAMPLES[*]}); do \
						for j in ${RESOLUTIONS[*]}; do \
							q8b coolpup.py.APA.AML2022.$AMLSAMPLES[${i}].${j}.over.22_12_05.${k}.loopAnchors.${m}.${o}.bed.source_loops.${k}.loops.AssociatedAnchorIs.${l}.SizeGroup.${n}.bedpe \
								coolpup.py \
									--n_proc 8 \
									--nshifts 100 \
									--excl_chrs chrX,chrY,chr17_GL000205v2_random \
									--pad 200 \
									--outdir $OUTDIR \
									--outname AML2022.$AMLSAMPLES[${i}].${j}.over.22_12_05.${k}.loopAnchors.${m}.${o}.bed.source_loops.${k}.loops.AssociatedAnchorIs.${l}.SizeGroup.${n}.bedpe \
									fanc/hic/binned/$AMLSAMPLES[${i}]/10kb_HiC_$AMLSAMPLES[${i}]_val_1.fq.gzorderedchr.cool::/resolutions/10000 \
									${WDIR}/LoopCoordinates/DifferentialAnchorsEnhancer/22_12_05.${k}.loopAnchors.${m}.${o}.bed.source_loops.${k}.loops.AssociatedAnchorIs.${l}.SizeGroup.${n}.bedpe;
						done;
					done;
				done;
			done;
		done;
	done;
done
