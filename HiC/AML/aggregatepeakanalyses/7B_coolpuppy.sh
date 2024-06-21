#!/bin/bash
# 7B_coolpuppy.sh

### Perform Aggregate Peak Analysis (APA) with coolpup.py on merged Hi-C matrices of Ctrl and Cohesin deficient cells/patients.
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
EDGES=(both minus plus either) # we tested whether the results differed when considering loops whose anchors overlap with cohesin only on one 'edge' or both 'edges' with similar results, we kept the 'either' category, as it was the biggest
DIRECTIONS=(Strengthened Weakened)
SIZEGROUPS=(1 2 3 4 5 6 7 All) # defined in 7B_cohesin_e_p_diff_loops_size_split.R # '<50kb','50kb-100kb','100kb-500kb','500kb-1Mb','1Mb-2Mb','2Mb-5Mb','>5Mb'
TYPE=(CohesinAssEnhancersOverlap H3K27acOverlap)
AMLSAMPLES=(merged_sa2_mut merged_rad21_mut merged_aml_ctrl merged_CD34_rad21 merged_CD34_siCtrl merged_CD34_stag1 merged_CD34_stag2)
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
									fanc/hic/binned/$AMLSAMPLES[${i}]/$AMLSAMPLES[${i}]_${j}.cool::/resolutions/10000 \
									${WDIR}/LoopCoordinates/DifferentialAnchorsEnhancer/22_12_05.${k}.loopAnchors.${m}.${o}.bed.source_loops.${k}.loops.AssociatedAnchorIs.${l}.SizeGroup.${n}.bedpe;
						done;
					done;
				done;
			done;
		done;
	done;
done
