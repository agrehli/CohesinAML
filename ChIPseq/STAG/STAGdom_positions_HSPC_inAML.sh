#!/bin/bash
# by Alexander Fischer, APR 2022

###############################################################################
###############################################################################
##                                                                           ##
## Analysis of STAG1 or STAG2 dominant or shared positions                   ##
##              (determined in HSPCs) in AMLs                                ##
##                                                                           ##
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
JUICERTOOLS="${DIR_PKG}/juicertools/run_juicertools-1.8.9.sh"


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

##  directories for analysis
PROJECTDIR="${DIR_DATA}/analysis/project_cohesin"
WORKDIR="${PROJECTDIR}/Cohesin_AML/ChIP_analysis"
PEAKDIR="${WORKDIR}/peaks"
ANNDIR="${WORKDIR}/annotation_tables"
DIFFPEAKS="${WORKDIR}/diffPeaks"
FIGURESDIR="${WORKDIR}/figures"
MOTIFDIR="${WORKDIR}/motifs"
MERGEDBW="${WORKDIR}/mergedBigWigs"
MERGEDBW_AML_H3K="${WORKDIR}/H3K27ac/mergedBigWigs"
MERGEDBW_AML_ATAC="${PROJECTDIR}/Cohesin_AML/ATAC/mergedBigWigs"
DTmatrixdir=${WORKDIR}/deeptoolsMatrix
MOTIFDIR="${WORKDIR}/motifs"

## HSPC data for comparisons
PEAKDIRCD34="${PROJECTDIR}/CD34/ChIP_KD_analysis/Cohesin_CTCF_MED12/peaks"
RAD21peaksSTAGdominant="${PROJECTDIR}/CD34/ChIP_KD_analysis/Cohesin_CTCF_MED12/diffPeaks/SA1vsSA2"
RAD21peaksSTAGdominantbyCTCF="${RAD21peaksSTAGdominant}/subsetByCTCF"
HSPCscaledCTCF="${PROJECTDIR}/CD34/ChIP_KD_analysis/Cohesin_CTCF_MED12/CTCF/scaledBigWigs"
HSPCscaledRAD21="${PROJECTDIR}/CD34/ChIP_KD_analysis/Cohesin_CTCF_MED12/RAD21/scaledBigWigs"

## DEG TSS pos
RNADEGDIR="${PROJECTDIR}/Cohesin_AML/RNAseq/Analysis/Resulttables/RNAseq_Cohesin_AML_AF3/STAG2pat_vs_CTRL"



                         ##########################################################################
                         # Histogram plots  at STAG1/2 dominant sites determined in HSPCs levels  #
                         ##########################################################################
                         
#annotate to positions with STAG2 or STAG1 dominace detected in HSPCs
mkdir -p ${FIGURESDIR}/ghist/STAGdomSites_inHPSCs/
PEAKsetsSTAGdom="allCD34.RAD21.peaks_edgeR.2folddown allCD34.RAD21.peaks_edgeR.2foldup allCD34.RAD21.peaks_edgeR.common"
######################################### SA1 and SA2 signals
for SET in ${PEAKsetsSTAGdom};do
	cat >"${TMPDIR}/ghist.${SET}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd /loctmp
annotatePeaks.pl ${RAD21peaksSTAGdominant}/SA2vsSA1.${SET}.txt  hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_CTRL_SA1.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/STAGdomSites_inHPSCs/AML_CTRL.SA2vs1.${SET}.CNVnorm.SA1.ghist.txt"
annotatePeaks.pl ${RAD21peaksSTAGdominant}/SA2vsSA1.${SET}.txt  hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_SA2mut_SA1.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/STAGdomSites_inHPSCs/AML_SA2mut.SA2vs1.${SET}.CNVnorm.SA1.ghist.txt"
annotatePeaks.pl ${RAD21peaksSTAGdominant}/SA2vsSA1.${SET}.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_CTRL_SA2.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/STAGdomSites_inHPSCs/AML_CTRL.SA2vs1.${SET}.CNVnorm.SA2.ghist.txt"
annotatePeaks.pl ${RAD21peaksSTAGdominant}/SA2vsSA1.${SET}.txt  hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_SA2mut_SA2.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/STAGdomSites_inHPSCs/AML_SA2mut.SA2vs1.${SET}.CNVnorm.SA2.ghist.txt"
EOF
	chmod 750 "${TMPDIR}/ghist.${SET}.${_DATE}.sh"
	echo "plotting ghists for ${SET}"
	screen -dm -S ${SET} bash -c "bash ${TMPDIR}/ghist.${SET}.${_DATE}.sh"
done


#make ghist plots showing average ChiPs at SA1vsSA2 diff pos. comparing control and one KD at a time for the same ChiP
chips="SA1 SA2"
for chip in ${chips};do
for SET in ${PEAKsetsSTAGdom};do
cd ${FIGURESDIR}/ghist/STAGdomSites_inHPSCs
plotHIST.sh -g "AML_CTRL.SA2vs1.${SET}.CNVnorm.${chip}.ghist.txt AML_SA2mut.SA2vs1.${SET}.CNVnorm.${chip}.ghist.txt" \
-s "CTRL_AML-${chip} STAG2mut-${chip}" -c "firebrick seagreen3" -x 1000 -y "0 16" -d ${FIGURESDIR}/ghist/SA1vsSA2/STAGdomSites_inHPSCs -n SA2vs1.${SET}.${chip}.SA2mutvsCTRL.AML
done
done

#make ghist plots showing average ChiPs at SA1vsSA2 diff pos. SA1 to SA2 ChiP within the conditions separately
for SET in ${PEAKsetsSTAGdom};do
cd ${FIGURESDIR}/ghist/STAGdomSites_inHPSCs
plotHIST.sh -g "AML_CTRL.SA2vs1.${SET}.CNVnorm.SA2.ghist.txt AML_CTRL.SA2vs1.${SET}.CNVnorm.topQC.SA1.ghist.txt" \
-s "CTRL_AML-SA2 CTRL_AML-SA1" -c "springgreen2 darkgoldenrod1" -x 1000 -y "0 16" -d ${FIGURESDIR}/ghist/SA1vsSA2/STAGdomSites_inHPSCs -n SA2vs1.${SET}.CTRL.AML.SA2vsSA1
plotHIST.sh -g "AML_SA2mut.SA2vs1.${SET}.CNVnorm.SA2.ghist.txt AML_SA2mut.SA2vs1.${SET}.CNVnorm.topQC.SA1.ghist.txt" \
-s "SA2mut-SA2 SA2mut-SA1" -c "springgreen2 darkgoldenrod1" -x 1000 -y "0 16" -d ${FIGURESDIR}/ghist/SA1vsSA2/STAGdomSites_inHPSCs -n SA2vs1.${SET}.SA2mut.SA2vsSA1
done


######################################### RAD21 CTCF H3K27ac signals
conds="CTRL SA2mut RAD21mut"

for cond in ${conds};do
for SET in ${PEAKsetsSTAGdom};do
	cat >"${TMPDIR}/ghist.${SET}.${cond}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd /loctmp
annotatePeaks.pl ${RAD21peaksSTAGdominant}/SA2vsSA1.${SET}.txt  hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_${cond}_CTCF.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/STAGdomSites_inHPSCs/AML_${cond}.SA2vs1.${SET}.CNVnorm.av.RAD21.ghist.txt"
annotatePeaks.pl ${RAD21peaksSTAGdominant}/SA2vsSA1.${SET}.txt  hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_${cond}_RAD21.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/STAGdomSites_inHPSCs/AML_${cond}.SA2vs1.${SET}.CNVnorm.av.CTCF.ghist.txt"
annotatePeaks.pl ${RAD21peaksSTAGdominant}/SA2vsSA1.${SET}.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW_AML_H3K}/ave.AML_${cond}_H3K27ac.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/STAGdomSites_inHPSCs/AML_${cond}.SA2vs1.${SET}.CNVnorm.av.H3K27ac.ghist.txt"
annotatePeaks.pl ${RAD21peaksSTAGdominant}/SA2vsSA1.${SET}.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW_AML_ATAC}/ave.${cond}.ATAC.bedGraph > "${FIGURESDIR}/ghist/STAGdomSites_inHPSCs/AML_${cond}.SA2vs1.${SET}.CNVnorm.av.ATAC.ghist.txt"
EOF
	chmod 750 "${TMPDIR}/ghist.${SET}.${cond}.${_DATE}.sh"
	echo "plotting ghists for ${SET} in ${cond} AML"
	screen -dm -S ${SET}.${cond}2 bash -c "bash ${TMPDIR}/ghist.${SET}.${cond}.${_DATE}.sh"
done
done



#SA1vsSA2domvscommon ghists
chromfeatures="RAD21 CTCF H3K27ac ATAC"
for cond in ${conds};do
for feat in ${chromfeatures};do
plotHIST.sh -g "AML_${cond}.SA2vs1.allCD34.RAD21.peaks_edgeR.2folddown.CNVnorm.av.${feat}.ghist.txt AML_${cond}.SA2vs1.allCD34.RAD21.peaks_edgeR.2foldup.CNVnorm.av.${feat}.ghist.txt AML_${cond}.SA2vs1.allCD34.RAD21.peaks_edgeR.common.CNVnorm.av.${feat}.ghist.txt" \
-s "STAG1dom STAG2dom common" -c "darkgoldenrod3 seagreen3 darkgrey" -x 1000 -y "0 13" -d ${FIGURESDIR}/ghist/STAGdomSites_inHPSCs -n SA1domvsSA2domvsCommon.${cond}.${feat}
done
done
#higher scale for CTCF
for cond in ${conds};do
feat="CTCF"
plotHIST.sh -g "AML_${cond}.SA2vs1.allCD34.RAD21.peaks_edgeR.2folddown.CNVnorm.av.${feat}.ghist.txt AML_${cond}.SA2vs1.allCD34.RAD21.peaks_edgeR.2foldup.CNVnorm.av.${feat}.ghist.txt AML_${cond}.SA2vs1.allCD34.RAD21.peaks_edgeR.common.CNVnorm.av.${feat}.ghist.txt" \
-s "STAG1dom STAG2dom common" -c "darkgoldenrod3 seagreen3 darkgrey" -x 1000 -y "0 30" -d ${FIGURESDIR}/ghist/STAGdomSites_inHPSCs -n SA1domvsSA2domvsCommon.${cond}.${feat}.max30
done
#larger region for H3K27ac
for cond in ${conds};do
feat="H3K27ac"
plotHIST.sh -g "AML_${cond}.SA2vs1.allCD34.RAD21.peaks_edgeR.2folddown.CNVnorm.av.${feat}.ghist.txt AML_${cond}.SA2vs1.allCD34.RAD21.peaks_edgeR.2foldup.CNVnorm.av.${feat}.ghist.txt AML_${cond}.SA2vs1.allCD34.RAD21.peaks_edgeR.common.CNVnorm.av.${feat}.ghist.txt" \
-s "STAG1dom STAG2dom common" -c "darkgoldenrod3 seagreen3 darkgrey" -x 2000 -y "0 13" -d ${FIGURESDIR}/ghist/STAGdomSites_inHPSCs -n SA1domvsSA2domvsCommon.${cond}.${feat}.2MB
done



                         ##########################################################################
                         # RAD21/CTCF coverage  at STAG1/2 dominant sites determined in HSPCs     #
                         ##########################################################################
                         
#Positions with SA1/SA2 dominace/equality in HSPCs: show RAD21 and CTCF in AMLs and HSPCs
mkdir ${FIGURESDIR}/Heatmaps/SA1_SA2dom_pos/


STAGdomsets="SA2dom SA1dom common"
for set in ${STAGdomsets};do
computeMatrix reference-point -S \
${HSPCscaledRAD21}/ave.CTRL.RAD21.scaled.bigWig \
${HSPCscaledRAD21}/ave.SA2KD.RAD21.scaled.bigWig \
${HSPCscaledRAD21}/ave.SA1KD.RAD21.scaled.bigWig \
${MERGEDBW}/ave.AML_CTRL_RAD21_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_SA2mut_RAD21_CNVnorm.bigwig \
${HSPCscaledCTCF}/ave.CTRL.CTCF.scaled.bigWig \
${HSPCscaledCTCF}/ave.SA2KD.CTCF.scaled.bigWig \
${HSPCscaledCTCF}/ave.SA1KD.CTCF.scaled.bigWig \
${MERGEDBW}/ave.AML_CTRL_CTCF_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_SA2mut_CTCF_CNVnorm.bigwig \
-R ${RAD21peaksSTAGdominantbyCTCF}/${set}_withCTCF.bed ${RAD21peaksSTAGdominantbyCTCF}/${set}_noCTCF.bed \
-b 2000 -a 2000 -bs 10 -p 30 --referencePoint center \
-o ${DTmatrixdir}/ChIP_SA2mutvsCTRLvsHSPCs_4000kbp_Matrix.${set}.RAD21.CTCF.mbw.gz
done

set_scale=(
'SA1dom' 10
'SA2dom' 6
'common' 20)
for (( idx=0 ; idx<${#set_scale[@]} ; idx+=2 )) ; do
    SET=${set_scale[idx]}
    scale=${set_scale[idx+1]}
plotHeatmap -m ${DTmatrixdir}/ChIP_SA2mutvsCTRLvsHSPCs_4000kbp_Matrix.${SET}.RAD21.CTCF.mbw.gz --sortRegions descend --sortUsing mean \
--colorList \
'white,mediumvioletred' 'white,mediumvioletred' 'white,mediumvioletred' 'white,mediumvioletred' 'white,mediumvioletred' \
'white,darkblue' 'white,darkblue' 'white,darkblue' 'white,darkblue' 'white,darkblue' \
--colorNumber 10 --refPointLabel '' \
--heatmapHeight $scale \
--yMax 26 \
--samplesLabel \
CTRLH SA2KD SA1KD CTRLA SA2mut CTRLH SA2KD SA1KD CTRLA SA2mut \
--whatToShow 'plot and heatmap' \
--legendLocation upper-left \
--missingDataColor 1 \
--regionsLabel withCTCF noCTCF \
--xAxisLabel '' \
-o ${FIGURESDIR}/Heatmaps/SA1_SA2dom_pos/ChIP_SA2mutvsCTRLvsHSPCs_4000kbp_Matrix.${SET}.RAD21.CTCF.mbw.pdf
done






                         ###################################################################################
                         # Histogram plots of STAG1/STAG2 signals at TSS of differentiall expressed genes  #
                         ###################################################################################
                         
#annotate to positions
mkdir -p ${FIGURESDIR}/ghist/DEG_TSSpos/
DEGsets="up down"
for SET in ${DEGsets};do
	cat >"${TMPDIR}/ghist.${SET}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd /loctmp
#annotatePeaks.pl ${RNADEGDIR}/STAG2mut.${SET}.FC1.5.CPM1.TSS.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_CTRL_SA1.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/DEG_TSSpos/AML_CTRL.DEG.${SET}.FC1.5.TSS.CNVnorm.SA1.ghist.txt"
#annotatePeaks.pl ${RNADEGDIR}/STAG2mut.${SET}.FC1.5.CPM1.TSS.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_SA2mut_SA1.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/DEG_TSSpos/AML_SA2mut.DEG.${SET}.FC1.5.TSS.CNVnorm.SA1.ghist.txt"
#annotatePeaks.pl ${RNADEGDIR}/STAG2mut.${SET}.FC1.5.CPM1.TSS.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_CTRL_SA2.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/DEG_TSSpos/AML_CTRL.DEG.${SET}.FC1.5.TSS.CNVnorm.SA2.ghist.txt"
#annotatePeaks.pl ${RNADEGDIR}/STAG2mut.${SET}.FC1.5.CPM1.TSS.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_SA2mut_SA2.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/DEG_TSSpos/AML_SA2mut.DEG.${SET}.FC1.5.TSS.CNVnorm.SA2.ghist.txt"

#annotatePeaks.pl ${RNADEGDIR}/DEGs${SET}.FC1.5.STAG2mut.RAD21${SET}.TSS.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_CTRL_SA1.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/DEG_TSSpos/AML_CTRL.DEGs${SET}.FC1.5.STAG2mut.RAD21${SET}.TSS.CNVnorm.SA1.ghist.txt"
#annotatePeaks.pl ${RNADEGDIR}/DEGs${SET}.FC1.5.STAG2mut.RAD21${SET}.TSS.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_CTRL_SA2.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/DEG_TSSpos/AML_CTRL.DEGs${SET}.FC1.5.STAG2mut.RAD21${SET}.TSS.CNVnorm.SA2.ghist.txt"
#annotatePeaks.pl ${RNADEGDIR}/DEGs${SET}.FC1.5.STAG2mut.RAD21${SET}.TSS.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_SA2mut_SA1.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/DEG_TSSpos/AML_SA2mut.DEGs${SET}.FC1.5.STAG2mut.RAD21${SET}.TSS.CNVnorm.SA1.ghist.txt"
#annotatePeaks.pl ${RNADEGDIR}/DEGs${SET}.FC1.5.STAG2mut.RAD21${SET}.TSS.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_SA2mut_SA2.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/DEG_TSSpos/AML_SA2mut.DEGs${SET}.FC1.5.STAG2mut.RAD21${SET}.TSS.CNVnorm.SA2.ghist.txt"

annotatePeaks.pl ${RNADEGDIR}/DEGs${SET}.FC1.5.STAG2mut.RAD21stable.TSS.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_CTRL_SA1.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/DEG_TSSpos/AML_CTRL.DEGs${SET}.FC1.5.STAG2mut.RAD21stable.TSS.CNVnorm.SA1.ghist.txt"
annotatePeaks.pl ${RNADEGDIR}/DEGs${SET}.FC1.5.STAG2mut.RAD21stable.TSS.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_CTRL_SA2.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/DEG_TSSpos/AML_CTRL.DEGs${SET}.FC1.5.STAG2mut.RAD21stable.TSS.CNVnorm.SA2.ghist.txt"
annotatePeaks.pl ${RNADEGDIR}/DEGs${SET}.FC1.5.STAG2mut.RAD21stable.TSS.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_SA2mut_SA1.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/DEG_TSSpos/AML_SA2mut.DEGs${SET}.FC1.5.STAG2mut.RAD21stable.TSS.CNVnorm.SA1.ghist.txt"
annotatePeaks.pl ${RNADEGDIR}/DEGs${SET}.FC1.5.STAG2mut.RAD21stable.TSS.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_SA2mut_SA2.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/DEG_TSSpos/AML_SA2mut.DEGs${SET}.FC1.5.STAG2mut.RAD21stable.TSS.CNVnorm.SA2.ghist.txt"
EOF
	chmod 750 "${TMPDIR}/ghist.${SET}.${_DATE}.sh"
	echo "plotting ghists for ${SET} DEG TSS"
	screen -dm -S ${SET} bash -c "bash ${TMPDIR}/ghist.${SET}.${_DATE}.sh"
done


#make ghist plots showing average ChiP signals at DEG TSS within the conditions separately

#all DEG TSS pos
for SET in ${DEGsets};do
plotHIST.sh -g "AML_CTRL.DEG.${SET}.FC1.5.TSS.CNVnorm.SA2.ghist.txt AML_CTRL.DEG.${SET}.FC1.5.TSS.CNVnorm.SA1.ghist.txt" \
-s "CTRL_AML-SA2 CTRL_AML-SA1" -c "springgreen2 darkgoldenrod1" -x 2000 -y "0 4.5" -d ${FIGURESDIR}/ghist/DEG_TSSpos -n DEGs.${SET}.FC1.5.CTRL.AML.SA2vsSA1
plotHIST.sh -g "AML_SA2mut.DEG.${SET}.FC1.5.TSS.CNVnorm.SA2.ghist.txt AML_SA2mut.DEG.${SET}.FC1.5.TSS.CNVnorm.SA1.ghist.txt" \
-s "SA2mut-SA2 SA2mut-SA1" -c "springgreen2 darkgoldenrod1" -x 2000 -y "0 4.5" -d ${FIGURESDIR}/ghist/DEG_TSSpos -n DEGs.${SET}.FC1.5.SA2mut.SA2vsSA1
done

#RAD21 accordingly changed
for SET in ${DEGsets};do
plotHIST.sh -g "AML_CTRL.DEGs${SET}.FC1.5.STAG2mut.RAD21${SET}.TSS.CNVnorm.SA2.ghist.txt AML_CTRL.DEGs${SET}.FC1.5.STAG2mut.RAD21${SET}.TSS.CNVnorm.SA1.ghist.txt" \
-s "CTRL_AML-SA2 CTRL_AML-SA1" -c "springgreen2 darkgoldenrod1" -x 2000 -y "0 2.5" -d ${FIGURESDIR}/ghist/DEG_TSSpos -n DEGs${SET}.FC1.5.STAG2mut.RAD21${SET}.CTRL.AML.SA2vsSA1
plotHIST.sh -g "AML_SA2mut.DEGs${SET}.FC1.5.STAG2mut.RAD21${SET}.TSS.CNVnorm.SA2.ghist.txt AML_SA2mut.DEGs${SET}.FC1.5.STAG2mut.RAD21${SET}.TSS.CNVnorm.SA1.ghist.txt" \
-s "SA2mut-SA2 SA2mut-SA1" -c "springgreen2 darkgoldenrod1" -x 2000 -y "0 2.5" -d ${FIGURESDIR}/ghist/DEG_TSSpos -n DEGs${SET}.FC1.5.STAG2mut.RAD21${SET}.SA2mut.SA2vsSA1
done
#####higher ylim
plotHIST.sh -g "AML_CTRL.DEGs${SET}.FC1.5.STAG2mut.RAD21${SET}.TSS.CNVnorm.SA2.ghist.txt AML_CTRL.DEGs${SET}.FC1.5.STAG2mut.RAD21${SET}.TSS.CNVnorm.SA1.ghist.txt" \
-s "CTRL_AML-SA2 CTRL_AML-SA1" -c "springgreen2 darkgoldenrod1" -x 2000 -y "0 4.5" -d ${FIGURESDIR}/ghist/DEG_TSSpos -n DEGs${SET}.FC1.5.STAG2mut.RAD21${SET}.CTRL.AML.SA2vsSA1.max5
plotHIST.sh -g "AML_SA2mut.DEGs${SET}.FC1.5.STAG2mut.RAD21${SET}.TSS.CNVnorm.SA2.ghist.txt AML_SA2mut.DEGs${SET}.FC1.5.STAG2mut.RAD21${SET}.TSS.CNVnorm.SA1.ghist.txt" \
-s "SA2mut-SA2 SA2mut-SA1" -c "springgreen2 darkgoldenrod1" -x 2000 -y "0 4.5" -d ${FIGURESDIR}/ghist/DEG_TSSpos -n DEGs${SET}.FC1.5.STAG2mut.RAD21${SET}.SA2mut.SA2vsSA1.max5
done

#RAD21stable
for SET in ${DEGsets};do
plotHIST.sh -g "AML_CTRL.DEGs${SET}.FC1.5.STAG2mut.RAD21stable.TSS.CNVnorm.SA2.ghist.txt AML_CTRL.DEGs${SET}.FC1.5.STAG2mut.RAD21stable.TSS.CNVnorm.SA1.ghist.txt" \
-s "CTRL_AML-SA2 CTRL_AML-SA1" -c "springgreen2 darkgoldenrod1" -x 2000 -y "0 2.5" -d ${FIGURESDIR}/ghist/DEG_TSSpos -n DEGs${SET}.FC1.5.STAG2mut.RAD21stable.CTRL.AML.SA2vsSA1
plotHIST.sh -g "AML_SA2mut.DEGs${SET}.FC1.5.STAG2mut.RAD21stable.TSS.CNVnorm.SA2.ghist.txt AML_SA2mut.DEGs${SET}.FC1.5.STAG2mut.RAD21stable.TSS.CNVnorm.SA1.ghist.txt" \
-s "SA2mut-SA2 SA2mut-SA1" -c "springgreen2 darkgoldenrod1" -x 2000 -y "0 2.5" -d ${FIGURESDIR}/ghist/DEG_TSSpos -n DEGs${SET}.FC1.5.STAG2mut.RAD21stable.SA2mut.SA2vsSA1
done

