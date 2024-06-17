#!/bin/bash
# by AF, APR 2022

###############################################################################
###############################################################################
##                   Analysis of differential RAD21 Peaks                    ##
##                  detected in STAG2-mutant AML patients                    ##
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
BOWTIE2="${DIR_PKG}/bowtie/bowtie2-2.3.4-linux-x86_64/bowtie2"


# general files
CHROMSIZES_HG38="${DIR_SOFT}/viewer/IGV/IGVTools_2.3.98/genomes/GRCh38.PRI_p10.chrom.sizes"
CHAINFILE="${DIR_PKG}/genome/chainFiles/hg19ToHg38.over.chain"
BLACKLIST_HG38="${DIR_DATA}/analysis/generalStuff/annotation/GRCh38/hg38.blacklist.bed"

# general directories
TMPDIR="/loctmp"
## directories with data created by mapChIP.sh pipeline
TAGDIR="${DIR_DATA}/processedData/tagDir/chromatin/GRCh38/ChIP/Cohesin_AML" #Tag directories (HOMER)
BWDIR="${DIR_DATA}/processedData/bigWig/chromatin/GRCh38/ChIP/Cohesin_AML" #bigwigs
INPUTDIR="${DIR_DATA}/processedData/tagDir/DNA/GRCh38/Input/AML" # Tag directorie (HOMER) of input DNA-sequencing data
CNVDIR="${DIR_DATA}/processedData/mapping/DNA/GRCh38/Input/AML/CNVdata" # processed input DNA-sequencing data CNV FreeC output

## new directories for analysis
PROJECTDIR="${DIR_DATA}/analysis/project_cohesin"
WORKDIR="${PROJECTDIR}/Cohesin_AML/ChIP_analysis"
PEAKDIR="${WORKDIR}/peaks"
PEAKDIRH3K="${WORKDIR}/H3K27ac/peaks"
ANNDIR="${WORKDIR}/annotation_tables"
DIFFPEAKS="${WORKDIR}/diffPeaks"
FIGURESDIR="${WORKDIR}/figures"
MERGEDBW="${WORKDIR}/mergedBigWigs"
MERGEDBW_AML_H3K="${WORKDIR}/H3K27ac/mergedBigWigs"
ATACBW="${PROJECTDIR}/Cohesin_AML/ATAC/mergedBigWigs"

DTmatrixdir=${WORKDIR}/deeptoolsMatrix
MOTIFDIR="${WORKDIR}/motifs"

mkdir $DTmatrixdir #directory for deeptools output

#Transcription Start Site positions (FANTOM CAGE seq data derived)
TSSpos="${PROJDIR}/TSS.bed"

#HSPC RAD21-Coverage tracks for comparison
HSPCscaledRAD21="${PROJECTDIR}/CD34/ChIP_KD_analysis/Cohesin_CTCF_MED12/RAD21/scaledBigWigs"

#variables for looping
FILTS="2folddown 2foldup"
MUTS="SA2mut RAD21mut"
#------------------------------------------------------------------------------------ 



            ############################################################
            #          diff peakset TSS overlap                        #
            ############################################################

##direct TSS overlap
mkdir ${DIFFPEAKS}/TSSOverlap

##for diffpeaks
for MUT in ${MUTS};do
for FILT in ${FILTS};do
pos2bed.pl ${DIFFPEAKS}/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.${FILT}.txt > ${DIFFPEAKS}/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.${FILT}.bed
$BEDTOOLS intersect -b ${DIFFPEAKS}/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.${FILT}.bed -a ${TSSpos} -u > ${DIFFPEAKS}/TSSOverlap/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.${FILT}.TSS.bed
done done

##same for all (stringent) RAD21 peaks .. could be used for RAD21 FC - RNA FC correlation analysis
#intersect with loj option (left outer join reports both A and B), needs to be filtered for overlaps later
pos2bed.pl ${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.stringent.txt  > ${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.stringent.bed
$BEDTOOLS intersect -b ${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.stringent.bed -a ${TSSpos} -loj > ${DIFFPEAKS}/TSSOverlap/All.CoAML.RAD21stringentPeaks.TSS.bed




#------------------------------------------------------------------------------------ 

            ############################################################
            #          diff peakset H3K27ac overlap                    #
            ############################################################
cd ${DIFFPEAKS}
mkdir diffPeakAssEnh
pos2bed.pl ${PEAKDIRH3K}/Allpat_mergePeaks_H3K27ac.CNVnorm.filtered.XYrem.peaks.txt > diffPeakAssEnh/CoAML_H3K27ac.tmp.peaks.bed
###intersect diffpeaks with H3K27acpeak set centered on H3K27ac! then intersect with TSS

for MUT in ${MUTS};do
for FILT in ${FILTS};do
$BEDTOOLS intersect -a diffPeakAssEnh/CoAML_H3K27ac.tmp.peaks.bed -b  ${DIFFPEAKS}/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.${FILT}.bed -u > diffPeakAssEnh/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.${FILT}.AssEnhancers.bed
$BEDTOOLS intersect -b diffPeakAssEnh/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.${FILT}.AssEnhancers.bed -a ${TSSpos} -u > diffPeakAssEnh/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.${FILT}.AssEnhancers.TSS.bed
done
done
#------------------------------------------------------------------------------------ 

                     ######################################
                     #       CTCF/RAD21 intersections     #
                     ######################################

#------------------------------------------------------------------------------------ 

#for association analyses: which RAD21 peak is associated with which CTCF? -> 2way overlaps
$BEDTOOLS intersect -a ${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.bed -b ${WORKDIRotherChIP}/peaks/Allpat_mergePeaks_CTCF.filtered.peaks.bed  -wao > ${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.CTCF.2wayOverlap.bed
$BEDTOOLS intersect -a ${PEAKDIR}/Allpat_mergePeaks_RAD21_stringent.filtered.peaks.bed -b ${WORKDIRotherChIP}/peaks/Allpat_mergePeaks_CTCF_stringent.filtered.peaks.bed  -wao > ${PEAKDIR}/Allpat_mergePeaks_str.RAD21.filtered.peaks.str.CTCF.2wayOverlap.bed


#eparate differential regions by CTCF overlap status
##divide sets by CTCF overlap: separates "enhancer" from "structural" cohesin pos
CTCFpeaks=${PROJECTDIR}/ChIP_analysis/peaks/Allpat_mergePeaks_CTCF.filtered.peaks.bed #63087
#for all peaks
$BEDTOOLS intersect -a  ${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.bed -b $CTCFpeaks -u > ${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.withCTCF.bed
$BEDTOOLS intersect -a  ${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.bed -b $CTCFpeaks -v > ${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.noCTCF.bed
#for diffpeaks
mkdir ${DIFFPEAKS}/subsetByCTCF
cd ${DIFFPEAKS}
for FILT in ${FILTS};do
for MUT in ${MUTS};do
$BEDTOOLS intersect -a  ${DIFFPEAKS}/${MUT}tvsCTRL.RAD21.stringentPeaks_DESEQ.model.${FILT}.bed -b $CTCFpeaks -u > ${DIFFPEAKS}/subsetByCTCF/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.${FILT}.withCTCF.bed
$BEDTOOLS intersect -a  ${DIFFPEAKS}/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.${FILT}.bed -b $CTCFpeaks -v > ${DIFFPEAKS}/subsetByCTCF/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.${FILT}.noCTCF.bed
done
#------------------------------------------------------------------------------------ 

#look at the numbers
for MUT in ${MUTS};do
for FILT in ${FILTS};do
wc -l ${DIFFPEAKS}/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.${FILT}.bed
wc -l ${DIFFPEAKS}/subsetByCTCF/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.${FILT}.withCTCF.bed
wc -l ${DIFFPEAKS}/subsetByCTCF/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.${FILT}.noCTCF.bed
done done
#------------------------------------------------------------------------------------ 




                         #################################################################
                         #   ChIP and ATAC signals at RAD21 differential positions       #
                         #################################################################
conda activate deepTools.3.3.0
#------------------------------------------------------------------------------------ 

#CTCF separated deeptools matrices/heatmaps
mkdir ${FIGURESDIR}/Heatmaps/
mkdir ${FIGURESDIR}/Heatmaps/RAD21diffpos

##using average bigwigs create deeptoolsmatrix with CTCF separated bedfiles (both at once)
for MUT in ${MUTS};do
for FILT in ${FILTS};do
computeMatrix reference-point -S \
${MERGEDBW}/ave.AML_CTRL_SA1_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_${MUT}_SA1_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_CTRL_SA2_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_${MUT}_SA2_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_CTRL_RAD21_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_${MUT}_RAD21_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_CTRL_CTCF_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_${MUT}_CTCF_CNVnorm.bigwig \
${MERGEDBW_AML_H3K}/ave.AML_CTRL_H3K27ac_CNVnorm.bigwig \
${MERGEDBW_AML_H3K}/ave.AML_${MUT}_H3K27ac_CNVnorm.bigwig \
${ATACBW}/ave.CTRL.ATAC.bigWig \
${ATACBW}/ave.${MUT}.ATAC.bigWig \
-R ${DIFFPEAKS}/subsetByCTCF/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.${FILT}.withCTCF.bed ${DIFFPEAKS}/subsetByCTCF/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.${FILT}.noCTCF.bed \
-b 2000 -a 2000 -bs 10 -p 30 --referencePoint center \
-o ${DTmatrixdir}/ChIP_${MUT}vsCTRL_4000kbp_Matrix_RAD21.stringentPeaks_DESEQ.${FILT}.complete.byCTCF.gz
done done
#------------------------------------------------------------------------------------ 

#generate the heatmaps for SA2mut vs CTRL
set_scale_small=(
'stringentPeaks_DESEQ.2folddown' 12.9
'stringentPeaks_DESEQ.2foldup' 10.6
)
for (( idx=0 ; idx<${#set_scale_small[@]} ; idx+=2 )) ; do
    set=${set_scale_small[idx]}
    scale=${set_scale_small[idx+1]}
plotHeatmap -m ${DTmatrixdir}/ChIP_SA2mutvsCTRL_4000kbp_Matrix_RAD21.${set}.complete.byCTCF.gz --sortRegions descend --sortUsing mean \
--colorList 'white,darkgoldenrod' 'white,darkgoldenrod' \
'white,lime' 'white,lime' \
'white,mediumvioletred' 'white,mediumvioletred' \
'white,darkblue' 'white,darkblue' \
'white,orange' 'white,orange' \
'white,salmon' 'white,salmon' \
--colorNumber 10 --refPointLabel '' \
--heatmapHeight $scale \
--yMax 6.8 \
--samplesLabel \
CTRL SA2mut CTRL SA2mut CTRL SA2mut CTRL SA2mut CTRL SA2mut CTRL SA2mut \
--whatToShow 'plot and heatmap' \
--legendLocation upper-left \
--missingDataColor 1 \
--regionsLabel withCTCF noCTCF \
--xAxisLabel '' \
-o ${FIGURESDIR}/Heatmaps/RAD21diffpos/SA2mutvsCTRL.RAD21.${set}.complete.byCTCF.pdf
done
#------------------------------------------------------------------------------------ 





                         #########################################################
                         #     Motif analysis  in RAD21 Diffpeaks by CTCF status #
                         #########################################################


MUT="SA2mut"
CTCFcats="withCTCF noCTCF"
MOTIFDIRbyCTCF=${MOTIFDIR}/byCTCFoverlap
mkdir ${MOTIFDIRbyCTCF}
#------------------------------------------------------------------------------------ 
for CTCFcat in ${CTCFcats};do
cat >"${TMPDIR}/motifs.${CTCFcat}.${_DATE}.sh" <<EOF
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd ${DIFFPEAKS}
bed2pos.pl ${DIFFPEAKS}/subsetByCTCF/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.2foldup.${CTCFcat}.bed > ${DIFFPEAKS}/subsetByCTCF/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.2foldup.${CTCFcat}.txt
bed2pos.pl ${DIFFPEAKS}/subsetByCTCF/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.2folddown.${CTCFcat}.bed > ${DIFFPEAKS}/subsetByCTCF/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.2folddown.${CTCFcat}.txt
#up RAD21 peaks
findMotifsGenome.pl ${DIFFPEAKS}/subsetByCTCF/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.2foldup.${CTCFcat}.txt hg38 ${MOTIFDIRbyCTCF}/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.2foldup.${CTCFcat} -size given -len 7,8,9,10,11,12,13,14 -p 6 -h
compareMotifs.pl ${MOTIFDIRbyCTCF}/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.2foldup.${CTCFcat}/homerMotifs.all.motifs ${MOTIFDIRbyCTCF}/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.2foldup.${CTCFcat}/final -reduceThresh .75 -matchThresh .6 -pvalue 1e-12 -info 1.5 -cpu 6
#down RAD21 peaks
findMotifsGenome.pl ${DIFFPEAKS}/subsetByCTCF/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.2folddown.${CTCFcat}.txt hg38 ${MOTIFDIRbyCTCF}/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.2folddown.${CTCFcat} -size given -len 7,8,9,10,11,12,13,14 -p 6 -h
compareMotifs.pl ${MOTIFDIRbyCTCF}/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.2folddown.${CTCFcat}/homerMotifs.all.motifs ${MOTIFDIRbyCTCF}/${MUT}vsCTRL.RAD21.stringentPeaks_DESEQ.model.2folddown.${CTCFcat}/final -reduceThresh .75 -matchThresh .6 -pvalue 1e-12 -info 1.5 -cpu 6
EOF
	chmod 750 "${TMPDIR}/motifs.${CTCFcat}.${_DATE}.sh"
	echo "motif analysis ${CTCFcat} "
	screen -dm -S ${CTCFcat}${_DATE} bash -c "bash ${TMPDIR}/motifs.${CTCFcat}.${_DATE}.sh"
done
#plotting top motif enrichments subset by CTCF status: see RAD21diff.peaks.motifplots.AML.rmd script!
#------------------------------------------------------------------------------------ 



                         #########################################################
                         #     RAD21 Coverage histograms at diff RAD21 pos        #
                         #########################################################

####annotate to total and differential RAD21 stringent peak merged positions of all AML pat RAD21 tracks
patgroups="CTRL SA2mut Rad21mut"
####run in annotate peaks on peak files using the averaged bedgraphs --> launch paralell screen sessions
#------------------------------------------------------------------------------------
_DATE=$(date +%s)
for group in ${patgroups};do
	cat >"${TMPDIR}/ghist.${group}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd /loctmp
#averaged bedgraphs
#annotatePeaks.pl ${PEAKDIR}/Allpat_mergePeaks_RAD21_stringent.filtered.peaks.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_${group}_RAD21.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/RAD21/averageNonScaled/RAD21_${group}.stringentRAD21pos.all.ghist.txt"
annotatePeaks.pl ${DIFFPEAKS}/SA2mutvsCTRL.RAD21.stringentPeaks_DESEQ.model.2folddown.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_${group}_RAD21.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/RAD21/averageNonScaled/RAD21_${group}.stringentRAD21pos.SA2mutdown.ghist.txt"
annotatePeaks.pl ${DIFFPEAKS}/SA2mutvsCTRL.RAD21.stringentPeaks_DESEQ.model.2foldup.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_${group}_RAD21.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/RAD21/averageNonScaled/RAD21_${group}.stringentRAD21pos.SA2mutup.ghist.txt"
annotatePeaks.pl ${DIFFPEAKS}/RAD21mutvsCTRL.RAD21.stringentPeaks_DESEQ.model.2folddown.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_${group}_RAD21.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/RAD21/averageNonScaled/RAD21_${group}.stringentRAD21pos.RAD21mutdown.ghist.txt"
annotatePeaks.pl ${DIFFPEAKS}/RAD21mutvsCTRL.RAD21.stringentPeaks_DESEQ.model.2foldup.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_${group}_RAD21.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/RAD21/averageNonScaled/RAD21_${group}.stringentRAD21pos.RAD21mutup.ghist.txt"
EOF
	chmod 750 "${TMPDIR}/ghist.${group}.${_DATE}.sh"
	echo "plotting ghists for ${group} "
	screen -dm -S ${group}${_DATE} bash -c "bash ${TMPDIR}/ghist.${group}.${_DATE}.sh"
done
#------------------------------------------------------------------------------------


#### plot the histograms
#------------------------------------------------------------------------------------

#ghist with average CNVnorm bedgraphs at RAD21 pos comparing: CTRL SA2mut Rad21mut
mkdir -p ${FIGURESDIR}/ghist/RAD21/averageNonScaled

cd ${FIGURESDIR}/ghist/RAD21/averageNonScaled
plotHIST.sh -g "RAD21_CTRL.RAD21pos.all.ghist.txt RAD21_SA2mut.RAD21pos.all.ghist.txt RAD21_Rad21mut.RAD21pos.all.ghist.txt" \
-s "AML_CTRL AML_SA2mut AML_Rad21mut" -c "firebrick mediumspringgreen magenta1" -x 1000 -y "0 12" -d ${FIGURESDIR}/ghist/RAD21/averageNonScaled -n RAD21.Peaks.all.CNVnorm.mut.vs.CTRL
plotHIST.sh -g "RAD21_CTRL.RAD21pos.all.ghist.txt RAD21_SA2mut.RAD21pos.all.ghist.txt RAD21_Rad21mut.RAD21pos.all.ghist.txt" \
-s "AML_CTRL AML_SA2mut AML_Rad21mut" -c "firebrick mediumspringgreen magenta1" -x 500 -y "1 8" -d ${FIGURESDIR}/ghist/RAD21/averageNonScaled -n RAD21.Peaks.all.CNVnorm.mut.vs.CTRL.zoomedIn

plotHIST.sh -g "RAD21_CTRL.RAD21pos.all.ghist.txt RAD21_SA2mut.RAD21pos.all.ghist.txt" \
-s "CTRL-AML STAG2-mut" -c "firebrick seagreen3" -x 1000 -y "0 12" -d ${FIGURESDIR}/ghist/RAD21/averageNonScaled -n RAD21.Peaks.all.CNVnorm.SA2mut.vs.CTRL
plotHIST.sh -g "RAD21_CTRL.RAD21pos.all.ghist.txt RAD21_Rad21mut.RAD21pos.all.ghist.txt" \
-s "CTRL-AML RAD21-mut" -c "firebrick magenta1" -x 1000 -y "0 12" -d ${FIGURESDIR}/ghist/RAD21/averageNonScaled -n RAD21.Peaks.all.CNVnorm.RAD21mut.vs.CTRL
plotHIST.sh -g "RAD21_CTRL.RAD21pos.all.ghist.txt RAD21_Rad21mut.RAD21pos.all.ghist.txt" \
-s "CTRL-AML RAD21-mut" -c "firebrick magenta1" -x 1000 -y "0 15" -d ${FIGURESDIR}/ghist/RAD21/averageNonScaled -n RAD21.Peaks.all.CNVnorm.RAD21mut.vs.CTRL.max15
#------------------------------------------------------------------------------------



                     ######################################################################
                     #     CTCF coverage ghist plots comparing up/down RAD21 peaks        #
                     ######################################################################
##visualize average RAD21 coverage across all TSS of differentially expressed genes (DEGs) defined in RNAseq analysis
###using either all TSS of the FC1.5.CPM1 DEGs or only those of DEGSs that correlate with RAD21 by FC direction
## DEG TSS pos

FIGDIRCTCF=${WORKDIRotherChIP}/figures/ghist/CTCF/RAD21diffpos/
mkdir $FIGDIRCTCF

#------------------------------------------------------------------------------------
patgroups="CTRL SA2mut Rad21mut"

_DATE=$(date +%s)
for group in ${patgroups};do
for FILT in ${FILTS};do
	cat >"${TMPDIR}/ghist.${group}${FILT}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd /loctmp
#use CNVnormalized averaged bedgraphs
annotatePeaks.pl ${DIFFPEAKS}/subsetByCTCF/SA2mutvsCTRL.RAD21.stringentPeaks_DESEQ.model.${FILT}.withCTCF.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_${group}_CTCF.CNVnorm.bedGraph > "${FIGDIRCTCF}/AML_${group}.SA2mut_diffRAD21.${FILT}.CNVnorm.CTCF.ghist.txt"
EOF
	chmod 750 "${TMPDIR}/ghist.${group}${FILT}.${_DATE}.sh"
	echo "plotting ghists for ${group}${FILT} "
	screen -dm -S ${group}${FILT}${_DATE} bash -c "bash ${TMPDIR}/ghist.${group}${FILT}.${_DATE}.sh"
done done

#------------------------------------------------------------------------------------
#plot using plotHIST.sh script (R-based plot)
###SA2 mut vs CTRL in one plot
for FILT in ${FILTS};do
plotHIST.sh -g "AML_CTRL.SA2mut_diffRAD21.${FILT}.CNVnorm.CTCF.ghist.txt AML_SA2mut.SA2mut_diffRAD21.${FILT}.CNVnorm.CTCF.ghist.txt" \
-s "CTRL-AML STAG2-mut" -c "firebrick3 seagreen3" -x 2000 -y "0 6" -d ${FIGDIRCTCF} -n CTCF.RAD21diffpeaks.${FILT}.SA2mut.vs.CTRL
done

###up vs down in one plot
for group in ${patgroups};do
plotHIST.sh -g "AML_${group}.SA2mut_diffRAD21.2foldup.CNVnorm.CTCF.ghist.txt AML_${group}.SA2mut_diffRAD21.2folddown.CNVnorm.CTCF.ghist.txt" \
-s "CTCF_increased_RAD21 CTCF_decreased_RAD21" -c "darkblue steelblue1" -x 2000 -y "0 6" -d ${FIGDIRCTCF} -n CTCF.RAD21diffpeaks.SA2mut.up.vs.down.${group}
done
#------------------------------------------------------------------------------------


                     ######################################################################
                     #     CTCF scores comparing up/down differenetial RAD21 peaks        #
                     ######################################################################
mkdir $MOTIFDIR/CTCFscores
#------------------------------------------------------------------------------------
#calculate CTCF motif log odds score distribution
CTCFmotif="/misc/data/analysis/project_BoundariesMOMACDC/ChIP/CTCF/motifs/MAC_CTRL_CTCF.20bp.consensus1/homerResults/motif1.motif"
for FILT in ${FILTS};do
annotatePeaks.pl ${DIFFPEAKS}/SA2mutvsCTRL.RAD21.stringentPeaks_DESEQ.model.${FILT}.txt hg38 -m $CTCFmotif -mscore > "${MOTIFDIR}/CTCFscores/SA2mutvsCTRL.RAD21.stringentPeaks_DESEQ.model.${FILT}.CTCFmotifscored.txt"
done
#------------------------------------------------------------------------------------






                         ###################################################################
                         #    Individual patient RAD21 tracks at RAD21 diffpos             #
                         ###################################################################
#------------------------------------------------------------------------------------ 
####calculate matrix
for PT in ${Peaktypes};do
for FILT in ${FILTS};do
computeMatrix reference-point -S \
${BWDIR}/ChIP_AML_ctr_10156_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_ctr_11374_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_ctr_14094_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_ctr_16911_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_ctr_18136_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_ctr_18367_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_ctr_18519_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_ctr_19049_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_ctr_19284_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_ctr_19405_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_ctr_19416_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_ctr_19504_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_ctr_20180_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_ctr_20458_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_ctr_20899_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_ctr_21047_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_ctr_21290_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_ctr_22023_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_ctr_22759_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_ctr_22874_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_ctr_9873_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_SA2_15365_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_SA2_15640_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_SA2_19142_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_SA2_22464_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_SA2_22525_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_SA2_24743_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_SA2_25458_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_SA2_27396_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_SA2_27493_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_SA2_28041_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_SA2_29728_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_SA2_31501_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_SA2_32635_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_SA2_9708_Rad21_CNVnorm.bigWig \
${BWDIR}/ChIP_AML_SA2_BB007_Rad21_CNVnorm.bigWig \
-R ${DIFFPEAKS}/SA2mutvsCTRL.RAD21.stringentPeaks_DESEQ.model.${FILT}.bed \
-b 2000 -a 2000 -bs 10 -p 30 --referencePoint center \
-o ${DTmatrixdir}/ChIP_SA2mutvsCTRL_4000kbp_Matrix_RAD21.stringentPeaks_DESEQ.${FILT}.RAD21individual.gz
done
done
#------------------------------------------------------------------------------------ 
set_scale=(
'stringentPeaks_DESEQ.2folddown' 12.5
'stringentPeaks_DESEQ.2foldup' 10.5
)

KMEANS="1"
#generate the heatmaps
for km in $KMEANS;do
for (( idx=0 ; idx<${#set_scale[@]} ; idx+=2 )) ; do
    set=${set_scale[idx]}
    scale=${set_scale[idx+1]}
    echo "set=$set"
    echo "scale=$scale"
    echo
    echo "creating matrix for regions ${set} with height ${scale} with ${km} kmeans" 
plotHeatmap -m ${DTmatrixdir}/ChIP_SA2mutvsCTRL_4000kbp_Matrix_RAD21.${set}.RAD21individual.gz --kmeans ${km} --sortRegions descend --sortUsing mean \
--colorList 'white,mediumvioletred' 'white,mediumvioletred' \
--colorNumber 10 \
--heatmapHeight $scale \
--heatmapWidth 1 \
--yMax 5 \
--samplesLabel \
C C C C C C C C C C C C C C C C C C C C C M M M M M M M M M M M M M M M \
--whatToShow 'plot and heatmap' \
--legendLocation none \
--labelRotation 90 \
--regionsLabel '' --refPointLabel '' --xAxisLabel '' \
-o ${FIGURESDIR}/Heatmaps/RAD21diffpos/ChIP_SA2mutvsCTRL_4000kbp_Matrix_RAD21.${set}.RAD21individual.pdf
done
done
#------------------------------------------------------------------------------------ 

## average tracks HSPCs and AMLs at the differential RAD21 positions determined in STAG2-mutant AMLs
for PT in ${Peaktypes};do
for FILT in ${FILTS};do
computeMatrix reference-point -S \
${MERGEDBW}/ave.AML_CTRL_RAD21_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_SA2mut_RAD21_CNVnorm.bigwig \
${HSPCscaledRAD21}/ave.CTRL.RAD21.scaled.bigWig \
${HSPCscaledRAD21}/ave.SA2KD.RAD21.scaled.bigWig \
${HSPCscaledRAD21}/ave.SA1KD.RAD21.scaled.bigWig \
-R ${DIFFPEAKS}/SA2mutvsCTRL.RAD21.stringentPeaks_DESEQ.model.${FILT}.bed \
-b 2000 -a 2000 -bs 10 -p 30 --referencePoint center \
-o ${DTmatrixdir}/ChIP_SA2mutvsCTRL_4000kbp_Matrix_RAD21.stringentPeaks_DESEQ.${FILT}.averageRAD21.gz
done
done

for (( idx=0 ; idx<${#set_scale[@]} ; idx+=2 )) ; do
    set=${set_scale[idx]}
    scale=${set_scale[idx+1]}
    echo "set=$set"
    echo "scale=$scale"
    echo
    echo "creating matrix for regions ${set} with height ${scale} with ${km} kmeans" 
plotHeatmap -m ${DTmatrixdir}/ChIP_SA2mutvsCTRL_4000kbp_Matrix_RAD21.${set}.averageRAD21.gz --kmeans ${km} --sortRegions descend --sortUsing mean \
--colorList 'white,mediumvioletred' 'white,mediumvioletred' \
--colorNumber 10 \
--heatmapHeight $scale \
--heatmapWidth 1 \
--yMax 5 \
--samplesLabel \
Ca Ma CHa S2KHa S1KHa \
--whatToShow 'plot and heatmap' \
--legendLocation none \
--labelRotation 90 \
--regionsLabel '' --refPointLabel '' --xAxisLabel '' \
-o ${FIGURESDIR}/Heatmaps/RAD21diffpos/ChIP_SA2mutvsCTRLvsHSPCs_4000kbp_Matrix_RAD21.${set}.averageRAD21.pdf
done
#------------------------------------------------------------------------------------ 


