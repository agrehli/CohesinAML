#!/bin/bash
# 6B_6G_coolpuppy.sh


### Perform Aggregate Peak Analysis (APA) with coolpup.py on merged Hi-C matrices of Ctrl and Cohesin deficient cells/patients.
### Enrichment is calculated against 100 shifted windows as controls

# cooler version: 0.9.25
# .hic object resolution: 10kb
# coordinates: 
# excluded chromosomes: chrX,chrY,chr17_GL000205v2_random
# S6C: example : chr5:156000000-160000000
# view: 'split'

# APA on merged Cohesin patients and merged CD34 Hi-C matrices at 10kb resolution
WDIR='/private/'
cd $WDIR
OUTDIR=$WDIR/APA/DifferentialLoopAnalysis/
mkdir -p $OUTDIR

COMPARISONS=(RAD21KDvsCTRL SA1KDvsCTRL SA2KDvsCTRL RAD21mutvsCTRL SA2mutvsCTRL)
FDRS=(0.001 0.01 0.05) # we tested 3 different significance thresholds for 'differential looping' with very similar results
DIRECTIONS=(Strengthened Weakened)
SIZEGROUPS=(1 2 3 4 5 6 7 All) # defined in 6B_diff_loops_size_split.R # '<50kb','50kb-100kb','100kb-500kb','500kb-1Mb','1Mb-2Mb','2Mb-5Mb','>5Mb'
AMLSAMPLES=(merged_sa2_mut merged_rad21_mut merged_aml_ctrl merged_CD34_rad21 merged_CD34_siCtrl merged_CD34_stag1 merged_CD34_stag2)
RESOLUTIONS=(10kb)

# Submit trimming jobs using a custom script for Sun Grid Engine's qsub command
# 'q8b' submits jobs with Grid Engine's qsub command allocating 8 cores and giving a custom jobid
# Usage: q8b <jobid> <command>

for k in ${COMPARISONS[*]}; do \
	for l in ${FDRS[*]}; do \
		for m in ${DIRECTIONS[*]}; do \
			for n in ${SIZEGROUPS[*]}; do \
				for i in $(seq 1 ${#AMLSAMPLES[*]}); do \
					for j in ${RESOLUTIONS[*]}; do \
						q8b coolpup.py.APA.AML2022.$AMLSAMPLES[${i}].${j}.over.22_12_05.${k}.loops.${m}.FDR.${l}.SizeGroup.${n}.bedpe \
							coolpup.py \
								--n_proc 8 \
								--nshifts 100 \
								--excl_chrs chrX,chrY,chr17_GL000205v2_random \
								--pad 200 \
								--outdir $OUTDIR \
								--outname AML2022.$AMLSAMPLES[${i}].${j}.over.22_12_05.${k}.loops.${m}.FDR.${l}.SizeGroup.${n}.bedpe.txt \
								fanc/hic/binned/$AMLSAMPLES[${i}]/$AMLSAMPLES[${i}].cool::/resolutions/10000 \
								LoopCoordinates/DifferentialLoopAnalysis/22_12_05.${k}.loops.${m}.FDR.${l}.SizeGroup.${n}.bedpe;
					done;
				done;
			done;
		done;
	done;
done
