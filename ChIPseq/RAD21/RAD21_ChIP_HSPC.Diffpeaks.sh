#!/bin/bash
# by AF, APR 2022

###############################################################################
###############################################################################
##                   Analysis of differential RAD21 Peaks                    ##
##                  detected in Cohesin subunit KD HSPCs                     ##
###############################################################################
###############################################################################

# setting basic path and linux OS
DIR_SOFT="/misc/software"
DIR_PKG="${DIR_SOFT}/ngs"
DIR_DATA="/misc/data"
OS=$(lsb_release -c |grep "^Codename" | awk -F: '{print $2}' | sed 's/[[:blank:]]//g')

# setting homer environment
PATH_PERL=${DIR_SOFT}/package/perl/perl-5.26.1/bin
PATH_SAMTOOLS=${DIR_PKG}/samtools/samtools-1.6/bin
PATH_HOMER=${DIR_PKG}/homer/v4.11/bin
PATH_R=${DIR_SOFT}/package/RBioC/3.4.3/bin
export PATH=${PATH_R}:${PATH_PERL}:${PATH_SAMTOOLS}:${PATH_HOMER}:${PATH}
export PATH

# defining the program versions
BEDTOOLS="${DIR_PKG}/bedtools/bedtools2-2.27.1/bin/bedtools"

# general files
CHROMSIZES_HG38="/misc/software/viewer/IGV/IGVTools_2.3.98/genomes/GRCh38.PRI_p10.chrom.sizes"
BLACKLIST_HG38="/misc/data/analysis/generalStuff/annotation/GRCh38/hg38.blacklist.bed"

# general directories
TMPDIR="/loctmp"

## directories with data created by mapChIP.sh pipeline
RAWDIR="${DIR_DATA}/rawData/chromatin/ChIP/CD34/"
FASTQCDIR="${DIR_DATA}/rawData/chromatin/ChIP/CD34/FastQC"
TAGDIR="${DIR_DATA}/processedData/tagDir/chromatin/GRCh38/ChIP/CD34/siRNA_KD" #tag directories ChIP (HOMER)
BWDIR="${DIR_DATA}/processedData/bigWig/chromatin/GRCh38/ChIP/CD34/siRNA_KD" #bigwigs
INPUTDIR="${DIR_DATA}/processedData/tagDir/DNA/GRCh38/Input/CD34" #tag directories DNA-input (HOMER)

##  directories for analysis
PROJECTDIR="${DIR_DATA}/analysis/project_cohesin"
WORKDIR="${PROJECTDIR}/CD34/ChIP_KD_analysis/Cohesin_CTCF_MED12"
PEAKDIR="${WORKDIR}/peaks"
DIFFPEAKS="${WORKDIR}/diffPeaks"
FIGURESDIR="${WORKDIR}/figures"
RAD21SCALEDBW="${WORKDIR}/RAD21/scaledBigWigs"
MERGEDBW="${WORKDIR}/mergedBigWigs"
MOTIFDIR="${WORKDIR}/motifs"
DTmatrixdir="${WORKDIR}/deeptoolsMatrix"

## Other ChIP/ATAC coverage tracks/peaks to show/compare
RAD21SCALEDBW="${WORKDIR}/RAD21/scaledBigWigs"
MED12SCALEDBW="${WORKDIR}/MED12/scaledBigWigs"
CTCFSCALEDBW="${WORKDIR}/CTCF/scaledBigWigs"
ATACSCALEDBW="${PROJECTDIR}/CD34/ATAC/scaledBigWigs"
H3K27acSCALEDBW="${PROJECTDIR}/CD34/ChIP_KD_analysis/H3K27ac/scaledBigWigs"
PEAKDIRH3K="${PROJECTDIR}/CD34/ChIP_KD_analysis/H3K27ac/peaks"

## AML data for comparisons
AMLRAD21DIR="${PROJECTDIR}/Cohesin_AML/ChIP_analysis"
AMLPEAKS="${AMLRAD21DIR}/peaks"
AMLDIFFPEAKS="${AMLRAD21DIR}/diffPeaks"

#Transcription Start Site positions (FANTOM CAGE seq data derived)
TSSpos="${PROJDIR}/TSS.bed"






            ############################################################
            #          diff RAD21 peakset TSS overlap                  #
            ############################################################
mkdir ${DIFFPEAKS}/TSSOverlap/
#------------------------------------------------------------------------------------ 

#same naming convention of diffpeak pos files:
cp ${DIFFPEAKS}/qstat_SA2KDvsCTRL.RAD21.peaks_edgeR.2foldDown.txt ${DIFFPEAKS}/SA2KDvsCTRL.RAD21.Peaks_edgeR.model.2folddown.txt 
cp ${DIFFPEAKS}/qstat_SA2KDvsCTRL.RAD21.peaks_edgeR.2foldUp.txt ${DIFFPEAKS}/SA2KDvsCTRL.RAD21.Peaks_edgeR.model.2foldup.txt 

##direct diffpeak overlaps
PEAKSETS="SA1KDvsCTRL.RAD21.Peaks_edgeR.model.2foldup SA1KDvsCTRL.RAD21.Peaks_edgeR.model.2folddown SA2KDvsCTRL.RAD21.Peaks_edgeR.model.2foldup SA2KDvsCTRL.RAD21.Peaks_edgeR.model.2folddown"
#convert to bedfiles and intersect
for PEAKSET in ${PEAKSETS}; do
pos2bed.pl ${DIFFPEAKS}/${PEAKSET}.txt > ${DIFFPEAKS}/${PEAKSET}.bed
$BEDTOOLS intersect -b ${DIFFPEAKS}/${PEAKSET}.bed  -a ${TSSpos} -u > ${DIFFPEAKS}/TSSOverlap/${PEAKSET}.TSS.bed
done

for PEAKSET in ${PEAKSETS}; do
wc -l ${DIFFPEAKS}/TSSOverlap/${PEAKSET}.TSS.bed 
done

###same for all RAD21 peaks could be used for RAD21 FC - RNA FC correlation analysis
pos2bed.pl ${PEAKDIR}/CD34_RAD21.filtered.peaks.txt > ${PEAKDIR}/CD34_RAD21.filtered.peaks.bed
#intersect with loj option (left outer join reports both A and B), needs to be filtered for overlaps later
$BEDTOOLS intersect -b ${PEAKDIR}/CD34_RAD21.filtered.peaks.bed -a ${TSSpos} -loj > ${DIFFPEAKS}/TSSOverlap/All.CD34.RAD21peaks.TSS.bed
#------------------------------------------------------------------------------------ 

            ############################################################
            #          diff peakset H3K27ac overlap                    #
            ############################################################
##cohesin ass enhancers at diffepeaks -> "indirect" + direct TSS overlap
#------------------------------------------------------------------------------------ 
cd ${DIFFPEAKS}
mkdir diffPeakAssEnh
pos2bed.pl ${PEAKDIRH3K}/mergedCD34_H3K27ac.filtered.peaks.txt > diffPeakAssEnh/CD34_H3K27ac.tmp.peaks.bed
#------------------------------------------------------------------------------------ 
###intersect diffpeaks with H3K27acpeak set centered on H3K27ac! then intersect with TSS
for PEAKSET in ${PEAKSETS}; do
$BEDTOOLS intersect -a diffPeakAssEnh/CD34_H3K27ac.tmp.peaks.bed -b ${DIFFPEAKS}/${PEAKSET}.bed -u > diffPeakAssEnh/${PEAKSET}.AssEnhancers.bed
$BEDTOOLS intersect -b diffPeakAssEnh/${PEAKSET}.AssEnhancers.bed -a ${TSSpos} -u > diffPeakAssEnh/${PEAKSET}.AssEnhancers.TSS.bed
done
#------------------------------------------------------------------------------------ 
for PEAKSET in ${PEAKSETS}; do
wc -l ${DIFFPEAKS}/diffPeakAssEnh/${PEAKSET}.AssEnhancers.TSS.bed
done
#------------------------------------------------------------------------------------ 

                     ######################################
                     #       CTCF/RAD21 intersections     #
                     ######################################

#------------------------------------------------------------------------------------ 
CTCFpeaks=${PEAKDIR}/CD34_CTCF.filtered.peaks.bed #60648
#for association analyses: which RAD21 peak is associated with which CTCF? -> 2way overlaps
$BEDTOOLS intersect -a ${PEAKDIR}/CD34_RAD21.filtered.peaks.bed -b $CTCFpeaks  -wao > ${PEAKDIR}/CD34_RAD21.filtered.peaks.CTCF.2wayOverlap.bed
#for diffpeaks:
######this makes mostly really sense for the downregulated peaks since upregulated ones seem to be almost fully CTCF associated
####Intersect diffpeak sets with full CTCF peakset
mkdir ${DIFFPEAKS}/subsetByCTCF
cd ${DIFFPEAKS}
for set in ${PEAKSETS};do
$BEDTOOLS intersect -a  ${DIFFPEAKS}/${set}.bed -b $CTCFpeaks -u > ${DIFFPEAKS}/subsetByCTCF/${set}.withCTCF.bed
$BEDTOOLS intersect -a  ${DIFFPEAKS}/${set}.bed -b $CTCFpeaks -v > ${DIFFPEAKS}/subsetByCTCF/${set}.noCTCF.bed
done
#------------------------------------------------------------------------------------ 



                         #################################################################
                         #   ChIP and ATAC signals at RAD21 differential positions       #
                         #################################################################
conda activate deepTools.3.3.0
#deeptools matrices/heatmaps
mkdir ${FIGURESDIR}/Heatmaps/
mkdir ${FIGURESDIR}/Heatmaps/RAD21diffpos
#------------------------------------------------------------------------------------ 
######SA2 KD vs CTRL
#without showing PU1 MED12 tracks (same layout as for AML patients)
for set in ${PEAKSETS};do
computeMatrix reference-point -S \
${MERGEDBW}/merged.ChIP_CTRL_SA1.bigWig \
${MERGEDBW}/merged.ChIP_SA2KD_SA1.bigWig \
${MERGEDBW}/merged.ChIP_CTRL_SA2.bigWig \
${MERGEDBW}/merged.ChIP_SA2KD_SA2.bigWig \
${RAD21SCALEDBW}/ave.CTRL.RAD21.scaled.bigWig \
${RAD21SCALEDBW}/ave.SA2KD.RAD21.scaled.bigWig \
${CTCFSCALEDBW}/ave.CTRL.CTCF.scaled.bigWig \
${CTCFSCALEDBW}/ave.SA2KD.CTCF.scaled.bigWig \
${H3K27acSCALEDBW}/ave.CD34_CTRL_H3K27ac.scaled.bigWig \
${H3K27acSCALEDBW}/ave.CD34_SA2KD_H3K27ac.scaled.bigWig \
${ATACSCALEDBW}/ave.ctrl.scaled.bigWig \
${ATACSCALEDBW}/ave.SA2KD.scaled.bigWig \
-R ${DIFFPEAKS}/${set}.bed \
-b 2000 -a 2000 -bs 10 -p 30 --referencePoint center \
-o ${DTmatrixdir}/ChIP_CD34.4000kbp_Matrix.${set}.gz
done

#------------------------------------------------------------------------------------ 

#plot heatmap
mkdir ${FIGURESDIR}/Heatmaps/smallversion/
set_scale_small=(
'SA2KDvsCTRL.RAD21.Peaks_edgeR.model.2foldup' 4
'SA2KDvsCTRL.RAD21.Peaks_edgeR.model.2folddown' 6.3
)
for (( idx=0 ; idx<${#set_scale_small[@]} ; idx+=2 )) ; do
    set=${set_scale_small[idx]}
    scale=${set_scale_small[idx+1]}
    echo "set=$set"
    echo "scale=$scale"
    echo
    echo "creating matrix for loop anchors ${set} with height ${scale} with"
plotHeatmap -m ${DTmatrixdir}/ChIP_CD34.4000kbp_Matrix.${set}.gz --sortRegions descend --sortUsing mean \
--colorList 'white,darkgoldenrod' 'white,darkgoldenrod' \
'white,lime' 'white,lime' \
'white,mediumvioletred' 'white,mediumvioletred' \
'white,darkblue' 'white,darkblue' \
'white,orange' 'white,orange' \
'white,lightsalmon' 'white,lightsalmon' \
--colorNumber 10 --refPointLabel Peak \
--heatmapHeight $scale \
--yMax 5 \
--samplesLabel \
CTRL SA2KD CTRL SA2KD CTRL SA2KD CTRL SA2KD CTRL SA2KD CTRL SA2KD \
--whatToShow 'plot and heatmap' \
--legendLocation none \
--regionsLabel '' --refPointLabel '' --xAxisLabel '' \
-o ${FIGURESDIR}/Heatmaps/RAD21diffpos/ChIP_CD34.4000kbp_Matrix.${set}.pdf
done
done
#------------------------------------------------------------------------------------ 
##for STAG1 KD diffpeaks - version without PU1./MED12
STAG1sets="SA1KDvsCTRL.RAD21.peaks_edgeR.2foldDown SA1KDvsCTRL.RAD21.peaks_edgeR.2foldUp"
for set in ${STAG1sets};do
computeMatrix reference-point -S \
${MERGEDBW}/merged.ChIP_CTRL_SA1.bigWig \
${MERGEDBW}/merged.ChIP_SA1KD_SA1.bigWig \
${MERGEDBW}/merged.ChIP_CTRL_SA2.bigWig \
${MERGEDBW}/merged.ChIP_SA1KD_SA2.bigWig \
${RAD21SCALEDBW}/ave.CTRL.RAD21.scaled.bigWig \
${RAD21SCALEDBW}/ave.SA1KD.RAD21.scaled.bigWig \
${CTCFSCALEDBW}/ave.CTRL.CTCF.scaled.bigWig \
${CTCFSCALEDBW}/ave.SA1KD.CTCF.scaled.bigWig \
${H3K27acSCALEDBW}/ave.CD34_CTRL_H3K27ac.scaled.bigWig \
${H3K27acSCALEDBW}/ave.CD34_SA1KD_H3K27ac.scaled.bigWig \
${ATACSCALEDBW}/ave.ctrl.scaled.bigWig \
${ATACSCALEDBW}/ave.SA1KD.scaled.bigWig \
-R ${DIFFPEAKS}/qstat_${set}.bed \
-b 2000 -a 2000 -bs 10 -p 30 --referencePoint center \
-o ${DTmatrixdir}/RAD21diffpeak_revised/ChIP_CD34.4000kbp_Matrix.${set}.wo.PU1.MED12.gz
done
##cave needs different ymax scales for up and down sets
set_scale=(
'SA1KDvsCTRL.RAD21.peaks_edgeR.2foldDown' 4 10
'SA1KDvsCTRL.RAD21.peaks_edgeR.2foldUp' 4 38
)

for (( idx=0 ; idx<${#set_scale[@]} ; idx+=3 )) ; do
    set=${set_scale[idx]}
    scale=${set_scale[idx+1]}
    ymaxval=${set_scale[idx+2]}
    echo "set=$set"
    echo "scale=$scale"
    echo "scale=$ymaxval"
    echo
    echo "creating matrix for loop anchors ${set} with height ${scale}"
plotHeatmap -m ${DTmatrixdir}/ChIP_CD34.4000kbp_Matrix.${set}.wo.PU1.MED12.gz --sortRegions descend --sortUsing mean \
--colorList 'white,darkgoldenrod' 'white,darkgoldenrod' \
'white,lime' 'white,lime' \
'white,mediumvioletred' 'white,mediumvioletred' \
'white,darkblue' 'white,darkblue' \
'white,orange' 'white,orange' \
'white,lightsalmon' 'white,lightsalmon' \
--colorNumber 10 --refPointLabel Peak \
--heatmapHeight $scale \
--yMax $ymaxval \
--samplesLabel \
CTRL SA1KD CTRL SA1KD CTRL SA1KD CTRL SA1KD CTRL SA1KD CTRL SA1KD \
--whatToShow 'plot and heatmap' \
--legendLocation none \
--regionsLabel '' --refPointLabel '' --xAxisLabel '' \
-o ${FIGURESDIR}/Heatmaps/RAD21diffpos/ChIP_CD34.4000kbp_Matrix.${set}.pdf
done
done
#------------------------------------------------------------------------------------ 

            ############################################################
            #          PU.1 coverage at differential RAD21 peaks       #
            ############################################################

##PU1 at RAD21 differential positions in STAG2 KD
mkdir ${FIGURESDIR}/ghist/PU1
CONDITIONS="CTRL SA1KD SA2KD"
PEAKSETS="SA2KDvsCTRL.RAD21.Peaks_edgeR.model.2foldup SA2KDvsCTRL.RAD21.Peaks_edgeR.model.2folddown"
#annotate PU1 bedgraphs to RAD21 diff peak positions
for PEAKSET in ${PEAKSETS};do
for COND in ${CONDITIONS};do
annotatePeaks.pl ${DIFFPEAKSRAD21}/${PEAKSET}.txt hg38 -size 2000 -hist 25 -ghist -bedGraph "${PU1SCALEDBW}/ave.CD34_${COND}_PU1.bedGraph" > "${FIGURESDIR}/ghist/PU1/${COND}_PU1.${PEAKSET}.ghist.txt"
done
done
#plot histograms SA2 vs CTRL
PEAKSET="SA2KDvsCTRL.RAD21.Peaks_edgeR.model"
plotHIST.sh -g "${FIGURESDIR}/ghist/PU1/SA2KD_PU1.${PEAKSET}.2folddown.ghist.txt ${FIGURESDIR}/ghist/PU1/CTRL_PU1.${PEAKSET}.2folddown.ghist.txt" \
-s "SA2KD CTRL" -c "seagreen3 firebrick" -x 1000 -y "0 16" -d ${FIGURESDIR}/ghist/PU1/RAD21diffpeaks -n PU1.scaled.${PEAKSET}.2folddown.SA2vsCTRLonly
plotHIST.sh -g "${FIGURESDIR}/ghist/PU1/SA2KD_PU1.${PEAKSET}.2foldup.ghist.txt ${FIGURESDIR}/ghist/PU1/CTRL_PU1.${PEAKSET}.2foldup.ghist.txt" \
-s "SA2KD CTRL" -c "seagreen3 firebrick" -x 1000 -y "0 1" -d ${FIGURESDIR}/ghist/PU1/RAD21diffpeaks -n PU1.scaled.${PEAKSET}.2foldup.SA2vsCTRLonly






            ############################################################
            #          Motif signatures in differential peaks          #
            ############################################################

# regular peakset edgeR results STAG2/STAG1 vs CTRL
#------------------------------------------------------------------------------------ 
PEAKSETS="SA2KDvsCTRL.RAD21.Peaks_edgeR.2foldUp SA2KDvsCTRL.RAD21.Peaks_edgeR.2foldDown SA1KDvsCTRL.RAD21.Peaks_edgeR.2foldUp SA1KDvsCTRL.RAD21.Peaks_edgeR.2foldDown"
cd ${DIFFPEAKS}
_DATE=$(date +%s)
for PEAKSET in ${PEAKSETS}; do
	cat >"${TMPDIR}/motifs.${PEAKSET}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd ${DIFFPEAKS}
findMotifsGenome.pl ${DIFFPEAKS}/qstat_${PEAKSET}.txt hg38 ${MOTIFDIR}/${PEAKSET} -size given -len 7,8,9,10,11,12,13,14 -p 6 -h
compareMotifs.pl ${MOTIFDIR}/${PEAKSET}/homerMotifs.all.motifs ${MOTIFDIR}/${PEAKSET}/final -reduceThresh .75 -matchThresh .6 -pvalue 1e-12 -info 1.5 -cpu 6
EOF
	chmod 750 "${TMPDIR}/motifs.${PEAKSET}.${_DATE}.sh"
	echo "finding motifs for ${PEAKSET}"
	screen -dm -S motifs.${PEAKSET} bash -c "bash ${TMPDIR}/motifs.${PEAKSET}.${_DATE}.sh"
done
#------------------------------------------------------------------------------------ 




                     ######################################################################
                     #     CTCF coverage ghist plots comparing up/down RAD21 peaks        #
                     ######################################################################

## look at CTCF signals at RAD21 diffpeaks for SA2mut SA2KD SA1KD diffpeaks in CTRL,SA1KD,SA2KD
#------------------------------------------------------------------------------------ 
mkdir ${FIGURESDIR}/ghist/RAD21diffpeaks/CTCF

conds="CTRL SA1KD SA2KD"
sets="2folddown 2foldup"

for cond in ${conds};do
for SET in ${sets};do
	cat >"${TMPDIR}/ghist.${SET}.${cond}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd /loctmp
annotatePeaks.pl ${AMLDIFFPEAKS}/subsetByCTCF/SA2mutvsCTRL.RAD21.stringentPeaks_DESEQ.model.${SET}.withCTCF.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${CTCFSCALEDBW}/ave.${cond}.CTCF.scaled.bedGraph > "${FIGURESDIR}/ghist/RAD21diffpeaks/CTCF/${cond}.RAD21diffpeaks_SA2mutAML.str.peaks.DESEQmodel.${SET}.withCTCF.av.scaled.CTCF.ghist.txt"
annotatePeaks.pl ${DIFFPEAKS}/SA2KDvsCTRL.RAD21.Peaks_edgeR.model.${SET}.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${CTCFSCALEDBW}/ave.${cond}.CTCF.scaled.bedGraph > "${FIGURESDIR}/ghist/RAD21diffpeaks/CTCF/${cond}.RAD21diffpeaks_SA2KD.edgeR.${SET}.av.scaled.CTCF.ghist.txt"
annotatePeaks.pl ${DIFFPEAKS}/SA1KDvsCTRL.RAD21.Peaks_edgeR.model.${SET}.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${CTCFSCALEDBW}/ave.${cond}.CTCF.scaled.bedGraph > "${FIGURESDIR}/ghist/RAD21diffpeaks/CTCF/${cond}.RAD21diffpeaks_SA1KD.edgeR.${SET}.av.scaled.CTCF.ghist.txt"
EOF
	chmod 750 "${TMPDIR}/ghist.${SET}.${cond}.${_DATE}.sh"
	echo "plotting ghists for ${SET} in ${cond} HSPCs"
	screen -dm -S ${SET}.${cond}2 bash -c "bash ${TMPDIR}/ghist.${SET}.${cond}.${_DATE}.sh"
done
done
#------------------------------------------------------------------------------------ 

##plot hists: compare up vs down within condition
for cond in ${conds};do
plotHIST.sh -g "${cond}.RAD21diffpeaks_SA2mutAML.str.peaks.DESEQmodel.2foldup.withCTCF.av.scaled.CTCF.ghist.txt ${cond}.RAD21diffpeaks_SA2mutAML.str.peaks.DESEQmodel.2folddown.withCTCF.av.scaled.CTCF.ghist.txt" \
-s "CTCF_increased_RAD21 CTCF_decreased_RAD21" -c "darkblue steelblue1" -x 2000 -y "0 8" -d ${FIGURESDIR}/ghist/RAD21diffpeaks/CTCF/ -n SA2mut.RAD21diffpeaks_up.vs.down_in.${cond}.HSPCs.CTCF

plotHIST.sh -g "${cond}.RAD21diffpeaks_SA2KD.edgeR.2foldup.av.scaled.CTCF.ghist.txt ${cond}.RAD21diffpeaks_SA2KD.edgeR.2folddown.av.scaled.CTCF.ghist.txt" \
-s "CTCF_increased_RAD21 CTCF_decreased_RAD21" -c "darkblue steelblue1" -x 2000 -y "0 6" -d ${FIGURESDIR}/ghist/RAD21diffpeaks/CTCF/ -n SA2KD.RAD21diffpeaks_up.vs.down_in.${cond}.HSPCs.CTCF

plotHIST.sh -g "${cond}.RAD21diffpeaks_SA1KD.edgeR.2foldup.av.scaled.CTCF.ghist.txt ${cond}.RAD21diffpeaks_SA1KD.edgeR.2folddown.av.scaled.CTCF.ghist.txt" \
-s "CTCF_increased_RAD21 CTCF_decreased_RAD21" -c "darkblue steelblue1" -x 2000 -y "0 12" -d ${FIGURESDIR}/ghist/RAD21diffpeaks/CTCF/ -n SA1KD.RAD21diffpeaks_up.vs.down_in.${cond}.HSPCs.CTCF
done
#------------------------------------------------------------------------------------ 



                     ########################################################
                     #     CTCF scores comparing up/down RAD21 peaks        #
                     ########################################################
mkdir $MOTIFDIR/CTCFscores
#------------------------------------------------------------------------------------ 
#calculate CTCF motif log odds score distribution
CTCFmotif="/misc/data/analysis/project_BoundariesMOMACDC/ChIP/CTCF/motifs/MAC_CTRL_CTCF.20bp.consensus1/homerResults/motif1.motif"
for SET in ${sets};do
annotatePeaks.pl ${DIFFPEAKS}/SA2KDvsCTRL.RAD21.Peaks_edgeR.model.${SET}.txt hg38 -m $CTCFmotif -mscore > "${MOTIFDIR}/CTCFscores/SA2KDvsCTRL.RAD21.Peaks_edgeR.model.${SET}.CTCFmotifscored.txt"
annotatePeaks.pl ${DIFFPEAKS}/SA1KDvsCTRL.RAD21.Peaks_edgeR.model.${SET}.txt hg38 -m $CTCFmotif -mscore > "${MOTIFDIR}/CTCFscores/SA1KDvsCTRL.RAD21.Peaks_edgeR.model.${SET}.CTCFmotifscored.txt"
done
#------------------------------------------------------------------------------------ 



