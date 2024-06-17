#!/bin/bash
# by Alexander Fischer, APR 2022

###############################################################################
###############################################################################
##                                                                           ##
## Analysis of STAG1 or STAG2 dominant or shared position in HSPCs           ##
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
DIFFPEAKS="${WORKDIR}/diffPeaks/SA1vsSA2"
FIGURESDIR="${WORKDIR}/figures"
RAD21SCALEDBW="${WORKDIR}/RAD21/scaledBigWigs"
MERGEDBW="${WORKDIR}/mergedBigWigs"
MOTIFDIR="${WORKDIR}/motifs"
DTmatrixdir="${WORKDIR}/deeptoolsMatrix"

## Other ChIP/ATAC coverage tracks/peaks to show/compare
RAD21SCALEDBW="${WORKDIR}/RAD21/scaledBigWigs"
MED12SCALEDBW="${WORKDIR}/MED12/scaledBigWigs"
PU1SCALEDBW="${WORKDIR}/PU1/scaledBigWigs"
CTCFSCALEDBW="${WORKDIR}/CTCF/scaledBigWigs"
ATACSCALEDBW="${PROJECTDIR}/CD34/ATAC/scaledBigWigs"
H3K27acSCALEDBW="${PROJECTDIR}/CD34/ChIP_KD_analysis/H3K27ac/scaledBigWigs"
PEAKDIRH3K="${PROJECTDIR}/CD34/ChIP_KD_analysis/H3K27ac/peaks"

#Transcription Start Site positions (FANTOM CAGE seq data derived)
TSSpos="${PROJDIR}/TSS.bed"


                         #####################################
                         #        CTCF association           #
                         #####################################

##divide sets by CTCF overlap:
CTCFpeaks=${PEAKDIR}/CD34_CTCF.filtered.peaks.bed #60648
mkdir ${DIFFPEAKS}/subsetByCTCF
#------------------------------------------------------------------------------------ 
#common
$BEDTOOLS intersect -a  ${DIFFPEAKS}/SA2vsSA1.allCD34.RAD21.peaks_edgeR.common.bed -b $CTCFpeaks -u > ${DIFFPEAKS}/subsetByCTCF/common_withCTCF.bed  #48247/75492
$BEDTOOLS intersect -a  ${DIFFPEAKS}/SA2vsSA1.allCD34.RAD21.peaks_edgeR.common.bed -b $CTCFpeaks -v > ${DIFFPEAKS}/subsetByCTCF/common_noCTCF.bed  #27245/75492
#SA2dom (=2foldup)
$BEDTOOLS intersect -a  ${DIFFPEAKS}/SA2vsSA1.allCD34.RAD21.peaks_edgeR.2foldup.bed -b $CTCFpeaks -u > ${DIFFPEAKS}/subsetByCTCF/SA2dom_withCTCF.bed #969/1319
$BEDTOOLS intersect -a  ${DIFFPEAKS}/SA2vsSA1.allCD34.RAD21.peaks_edgeR.2foldup.bed -b $CTCFpeaks -v > ${DIFFPEAKS}/subsetByCTCF/SA2dom_noCTCF.bed #350/1319
#SA1dom (=2folddown)
$BEDTOOLS intersect -a  ${DIFFPEAKS}/SA2vsSA1.allCD34.RAD21.peaks_edgeR.2folddown.bed -b $CTCFpeaks -u > ${DIFFPEAKS}/subsetByCTCF/SA1dom_withCTCF.bed #927/2056
$BEDTOOLS intersect -a  ${DIFFPEAKS}/SA2vsSA1.allCD34.RAD21.peaks_edgeR.2folddown.bed -b $CTCFpeaks -v > ${DIFFPEAKS}/subsetByCTCF/SA1dom_noCTCF.bed #1129/2056
#All
$BEDTOOLS intersect -a  ${PEAKDIR}/CD34_RAD21.filtered.peaks.bed -b $CTCFpeaks -u > ${DIFFPEAKS}/subsetByCTCF/All_withCTCF.bed #50143
$BEDTOOLS intersect -a  ${PEAKDIR}/CD34_RAD21.filtered.peaks.bed -b $CTCFpeaks -v > ${DIFFPEAKS}/subsetByCTCF/All_noCTCF.bed #28724
#------------------------------------------------------------------------------------ 



                         ##############################################################################
                         #       All types of ChIP/ATAC coverage at dominant/common positions         #
                         ##############################################################################

conda activate deepTools.3.3.0
mkdir ${FIGURESDIR}/Heatmaps/SA1vsSA2/
#------------------------------------------------------------------------------------ 
#CTRL HSPCs only: separate clustering by CTCF categories
##run withCTCF and noCTCF region bed files at the same time
sets="SA2dom SA1dom common"
for set in ${sets};do
computeMatrix reference-point -S \
${MERGEDBW}/merged.ChIP_CTRL_SA1.bigWig \
${MERGEDBW}/merged.ChIP_CTRL_SA2.bigWig \
${RAD21SCALEDBW}/ave.CTRL.RAD21.scaled.bigWig \
${CTCFSCALEDBW}/ave.CTRL.CTCF.scaled.bigWig \
${MED12SCALEDBW}/ave.CD34_CTRL_MED12.scaled.bigWig \
${PU1SCALEDBW}/ave.CD34_CTRL_PU1.scaled.bigWig \
${H3K27acSCALEDBW}/ave.CD34_CTRL_H3K27ac.scaled.bigWig \
${ATACSCALEDBW}/ave.ctrl.scaled.bigWig \
-R ${DIFFPEAKS}/subsetByCTCF/${set}_withCTCF.bed ${DIFFPEAKS}/subsetByCTCF/${set}_noCTCF.bed \
-b 2000 -a 2000 -bs 10 -p 30 --referencePoint center \
-o ${DTmatrixdir}/ChIP_CD34_CTRLs_4000kbp_Matrix_SA1vsSA2_${set}.byCTCFstat.complete.gz
done


###plot without kmeans option to keep the two region sets separated
set_scale=(
'SA1dom' 10
'SA2dom' 6
'common' 20)
for (( idx=0 ; idx<${#set_scale[@]} ; idx+=2 )) ; do
    set=${set_scale[idx]}
    scale=${set_scale[idx+1]}
plotHeatmap -m ${DTmatrixdir}/ChIP_CD34_CTRLs_4000kbp_Matrix_SA1vsSA2_${set}.byCTCFstat.complete.gz --sortRegions descend --sortUsing mean \
--colorList 'white,darkgoldenrod' \
'white,lime' \
'white,mediumvioletred' \
'white,darkblue' \
'white,darkgrey' \
'white,mediumturquoise' \
'white,orange' \
'white,lightsalmon' \
--colorNumber 10 --refPointLabel '' \
--heatmapHeight $scale \
--yMax 26 \
--samplesLabel \
SA1 SA2 RAD21 CTCF MED12 PU1 H3K27ac ATAC \
--whatToShow 'plot and heatmap' \
--legendLocation upper-left \
--missingDataColor 1 \
--regionsLabel withCTCF noCTCF \
--xAxisLabel '' \
-o ${FIGURESDIR}/Heatmaps/SA1vsSA2/CTRLs_${set}.byCTCFstat.pdf
done
#------------------------------------------------------------------------------------ 


            ############################################################
            #                        ghist plots					   #
            ############################################################
#------------------------------------------------------------------------------------ 
#Average coverage historgrams at STAG dom/common pos in CTRLs and KDs
##annotate to edgeR diffpos (RAD21pos annotated data using total Peaks)


#STAG1 vs STAG2 signals
PEAKsets="allCD34.RAD21.peaks_edgeR.2folddown allCD34.RAD21.peaks_edgeR.2foldup allCD34.RAD21.peaks_edgeR.common"
groups="SA1KD SA2KD CTRL"
chips="SA1 SA2"
mkdir ${FIGURESDIR}/ghist/SA1vsSA2/
#------------------------------------------------------------------------------------ 

for SET in ${PEAKsets};do
for grp in ${groups};do
for chip in ${chips};do
	cat >"${TMPDIR}/ghist.${SET}.${grp}.${chip}${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd /loctmp
annotatePeaks.pl ${DIFFPEAKS}/SA2vsSA1.${SET}.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/merged.${grp}_${chip}.bedGraph > "${FIGURESDIR}/ghist/SA1vsSA2/mbw/${grp}.SA2vsSA1.${SET}.${chip}.ghist.txt"
EOF
	chmod 750 "${TMPDIR}/ghist.${SET}.${grp}.${chip}${_DATE}.sh"
	echo "plotting ghists for ${SET} in ${grp} for ${chip} ChIP"
	screen -dm -S ${SET}.${grp}.${chip} bash -c "bash ${TMPDIR}/ghist.${SET}.${grp}.${chip}${_DATE}.sh"
done
done
done
#------------------------------------------------------------------------------------ 


#make ghist plots showing average ChiPs at SA1vsSA2 diff pos. comparing SA1 to SA2 chip signal within same condition
for SET in ${PEAKsets};do
for grp in ${groups};do
cd ${FIGURESDIR}/ghist/SA1vsSA2/
plotHIST.sh -g "${grp}.SA2vsSA1.${SET}.SA2.ghist.txt ${grp}.SA2vsSA1.${SET}.SA1.ghist.txt" \
-s "${grp}-SA2 ${grp}-SA1" -c "springgreen2 darkgoldenrod1" -x 1000 -y "0 16" -d ${FIGURESDIR}/ghist/SA1vsSA2 -n SA2vsSA1.${SET}.${grp}vsCTRL.SA2vsSA1
done
done
#------------------------------------------------------------------------------------ 



#RAD21
#------------------------------------------------------------------------------------ 
#####for all dom. region types and all RAD21peaks
cp ${PEAKDIR}/CD34_RAD21.filtered.peaks.txt ${DIFFPEAKS}/SA2vsSA1.all.txt
for SET in ${PEAKsets};do
for grp in ${groups};do
	cat >"${TMPDIR}/ghist.${SET}.${grp}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd /loctmp
annotatePeaks.pl ${DIFFPEAKS}/SA2vsSA1.${SET}.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${RAD21SCALEDBW}/ave.${grp}.RAD21.scaled.bedGraph > "${FIGURESDIR}/ghist/SA1vsSA2/${grp}.SA2vsSA1.${SET}.RAD21scaled.ghist.txt"
EOF
	chmod 750 "${TMPDIR}/ghist.${SET}.${grp}.${_DATE}.sh"
	echo "plotting ghists for ${SET} in ${grp}"
	screen -dm -S ${SET}.${grp} bash -c "bash ${TMPDIR}/ghist.${SET}.${grp}.${_DATE}.sh"
done
done
#------------------------------------------------------------------------------------ 

#all conditions compared
for SET in ${PEAKsets};do
plotHIST.sh -g "CTRL.SA2vsSA1.${SET}.RAD21scaled.ghist.txt SA1KD.SA2vsSA1.${SET}.RAD21scaled.ghist.txt SA2KD.SA2vsSA1.${SET}.RAD21scaled.ghist.txt" \
-s "CTRL-HSPCs STAG1KD STAG2KD" -c "firebrick darkgoldenrod seagreen2" -x 1000 -y "0 15" -d ${FIGURESDIR}/ghist/SA1vsSA2/ -n SA2vsSA1.${SET}.RAD21scaled.KDs.vs.CTRL
done
#SA1dom vs SA2dom vs common in CTRLs
plotHIST.sh -g "CTRL.SA2vsSA1.allCD34.RAD21.peaks_edgeR.2folddown.RAD21scaled.ghist.txt CTRL.SA2vsSA1.allCD34.RAD21.peaks_edgeR.2foldup.RAD21scaled.ghist.txt CTRL.SA2vsSA1.allCD34.RAD21.peaks_edgeR.common.RAD21scaled.ghist.txt" \
-s "STAG1dom STAG2dom common" -c "darkgoldenrod3 seagreen3 darkgrey" -x 1000 -y "0 16" -d ${FIGURESDIR}/ghist/SA1vsSA2 -n SA1domvsSA2domvsCommon.RAD21scaled.CTRL
#------------------------------------------------------------------------------------ 


### CTCF ATAC and H3K27ac ghists for CTRLS
#------------------------------------------------------------------------------------ 
for SET in ${PEAKsets};do
	cat >"${TMPDIR}/ghist.${SET}.${grp}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd /loctmp
annotatePeaks.pl ${DIFFPEAKS}/SA2vsSA1.${SET}.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${CTCFSCALEDBW}/ave.CTRL.CTCF.scaled.bedGraph > "${FIGURESDIR}/ghist/SA1vsSA2/CTRL.SA2vsSA1.${SET}.CTCFscaled.ghist.txt"
annotatePeaks.pl ${DIFFPEAKS}/SA2vsSA1.${SET}.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${H3K27acSCALEDBW}/ave.CD34_CTRL_H3K27ac.bedGraph > "${FIGURESDIR}/ghist/SA1vsSA2/CTRL.SA2vsSA1.${SET}.H3K27acscaled.ghist.txt"
annotatePeaks.pl ${DIFFPEAKS}/SA2vsSA1.${SET}.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${ATACSCALEDBW}/ave.ctrl.scaled.bedGraph > "${FIGURESDIR}/ghist/SA1vsSA2/CTRL.SA2vsSA1.${SET}.ATACscaled.ghist.txt"
EOF
	chmod 750 "${TMPDIR}/ghist.${SET}.${grp}.${_DATE}.sh"
	echo "plotting ghists for ${SET} in ${grp}"
	screen -dm -S ${SET}.${grp} bash -c "bash ${TMPDIR}/ghist.${SET}.${grp}.${_DATE}.sh"
done

#SA1vsSA2domvscommon ghists
#CTCF
plotHIST.sh -g "CTRL.SA2vsSA1.allCD34.RAD21.peaks_edgeR.2folddown.CTCFscaled.ghist.txt CTRL.SA2vsSA1.allCD34.RAD21.peaks_edgeR.2foldup.CTCFscaled.ghist.txt CTRL.SA2vsSA1.allCD34.RAD21.peaks_edgeR.common.CTCFscaled.ghist.txt" \
-s "STAG1dom STAG2dom common" -c "darkgoldenrod3 seagreen3 darkgrey" -x 1000 -y "0 30" -d ${FIGURESDIR}/ghist/SA1vsSA2 -n SA1domvsSA2domvsCommon.CTCFscaled.CTRL
#H3K27ac
plotHIST.sh -g "CTRL.SA2vsSA1.allCD34.RAD21.peaks_edgeR.2folddown.H3K27acscaled.ghist.txt CTRL.SA2vsSA1.allCD34.RAD21.peaks_edgeR.2foldup.H3K27acscaled.ghist.txt CTRL.SA2vsSA1.allCD34.RAD21.peaks_edgeR.common.H3K27acscaled.ghist.txt" \
-s "STAG1dom STAG2dom common" -c "darkgoldenrod3 seagreen3 darkgrey" -x 2000 -y "0 13" -d ${FIGURESDIR}/ghist/SA1vsSA2 -n SA1domvsSA2domvsCommon.H3K27acscaled.CTRL
#ATAC
plotHIST.sh -g "CTRL.SA2vsSA1.allCD34.RAD21.peaks_edgeR.2folddown.ATACscaled.ghist.txt CTRL.SA2vsSA1.allCD34.RAD21.peaks_edgeR.2foldup.ATACscaled.ghist.txt CTRL.SA2vsSA1.allCD34.RAD21.peaks_edgeR.common.ATACscaled.ghist.txt" \
-s "STAG1dom STAG2dom common" -c "darkgoldenrod3 seagreen3 darkgrey" -x 1000 -y "0 13" -d ${FIGURESDIR}/ghist/SA1vsSA2 -n SA1domvsSA2domvsCommon.ATACscaled.CTRL
#------------------------------------------------------------------------------------ 




#Subset dominant/common regions by CTCF presence: RAD21 ghists
mkdir ${FIGURESDIR}/ghist/SA1vsSA2/subsetByCTCF
#------------------------------------------------------------------------------------ 
CTCFsubsets="SA2dom_withCTCF SA2dom_noCTCF SA1dom_withCTCF SA1dom_noCTCF common_withCTCF common_noCTCF"
for SET in ${CTCFsubsets};do
for grp in ${groups};do
	cat >"${TMPDIR}/ghist.${SET}.${grp}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd /loctmp
annotatePeaks.pl ${DIFFPEAKS}/subsetByCTCF/${SET}.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${RAD21SCALEDBW}/ave.${grp}.RAD21.scaled.bedGraph > "${FIGURESDIR}/ghist/SA1vsSA2/subsetByCTCF/${grp}.${SET}.RAD21scaled.ghist.txt"
EOF
	chmod 750 "${TMPDIR}/ghist.${SET}.${grp}.${_DATE}.sh"
	echo "plotting ghists for ${SET} in ${grp}"
	screen -dm -S ${SET}.${grp} bash -c "bash ${TMPDIR}/ghist.${SET}.${grp}.${_DATE}.sh"
done
done
#------------------------------------------------------------------------------------ 
withCTCFsubsets="SA2dom_withCTCF SA1dom_withCTCF common_withCTCF"
for SET in ${withCTCFsubsets};do
cd ${FIGURESDIR}/ghist/SA1vsSA2/subsetByCTCF/
plotHIST.sh -g "CTRL.${SET}.RAD21scaled.ghist.txt SA1KD.${SET}.RAD21scaled.ghist.txt SA2KD.${SET}.RAD21scaled.ghist.txt" \
-s "CTRL-HSPCs STAG1KD STAG2KD" -c "firebrick darkgoldenrod seagreen2" -x 1000 -y "0 25" -d ${FIGURESDIR}/ghist/SA1vsSA2/subsetByCTCF -n ${SET}.RAD21scaled.KDs.vs.CTRL
done

noCTCFsubsets="SA2dom_noCTCF SA1dom_noCTCF common_noCTCF"
for SET in ${noCTCFsubsets};do
cd ${FIGURESDIR}/ghist/SA1vsSA2/subsetByCTCF/
plotHIST.sh -g "CTRL.${SET}.RAD21scaled.ghist.txt SA1KD.${SET}.RAD21scaled.ghist.txt SA2KD.${SET}.RAD21scaled.ghist.txt" \
-s "CTRL-HSPCs STAG1KD STAG2KD" -c "firebrick darkgoldenrod seagreen2" -x 1000 -y "0 8" -d ${FIGURESDIR}/ghist/SA1vsSA2/subsetByCTCF -n ${SET}.RAD21scaled.KDs.vs.CTRL
done
#------------------------------------------------------------------------------------ 




            ############################################################
            #      Motif signatures in STAG dom/commmon peaks          #
            ############################################################


########################divide by CTCF subset: motifs based based on all CD34 cond RAD21 pos using only the edgeR results
##regular motif analyis
MOTIFDIR_SA1vs2ctcf=${MOTIFDIR}/SA1vsSA2byctcf
mkdir $MOTIFDIR_SA1vs2ctcf
sets="SA2dom SA1dom common"
#------------------------------------------------------------------------------------ 

_DATE=$(date +%s)
for PEAKSET in ${sets};do
	cat >"${TMPDIR}/motifs.${PEAKSET}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd ${DIFFPEAKS}/subsetByCTCF/
bed2pos.pl ${PEAKSET}_withCTCF.bed > ${PEAKSET}_withCTCF.txt
bed2pos.pl ${PEAKSET}_noCTCF.bed > ${PEAKSET}_noCTCF.txt
findMotifsGenome.pl ${PEAKSET}_noCTCF.txt hg38 ${MOTIFDIR_SA1vs2ctcf}/${PEAKSET}_noCTCF -size given -len 7,8,9,10,11,12,13,14 -p 6 -h
compareMotifs.pl ${MOTIFDIR_SA1vs2ctcf}/${PEAKSET}_noCTCF/homerMotifs.all.motifs ${MOTIFDIR_SA1vs2ctcf}/${PEAKSET}_noCTCF/final -reduceThresh .75 -matchThresh .6 -pvalue 1e-12 -info 1.5 -cpu 6
findMotifsGenome.pl ${PEAKSET}_withCTCF.txt hg38 ${MOTIFDIR_SA1vs2ctcf}/${PEAKSET}_withCTCF -size given -len 7,8,9,10,11,12,13,14 -p 6 -h
compareMotifs.pl ${MOTIFDIR_SA1vs2ctcf}/${PEAKSET}_withCTCF/homerMotifs.all.motifs ${MOTIFDIR_SA1vs2ctcf}/${PEAKSET}_withCTCF/final -reduceThresh .75 -matchThresh .6 -pvalue 1e-12 -info 1.5 -cpu 6
EOF
	chmod 750 "${TMPDIR}/motifs.${PEAKSET}.${_DATE}.sh"
	echo "finding motifs for ${PEAKSET}"
	screen -dm -S motifs.${PEAKSET} bash -c "bash ${TMPDIR}/motifs.${PEAKSET}.${_DATE}.sh"
done
#------------------------------------------------------------------------------------ 



            ############################################################
            #      Gene association by TSS overlap                     #
            ############################################################
###by TSS overlap --> all promoters
mkdir ${DIFFPEAKS}/TSSOverlap/
#------------------------------------------------------------------------------------ 
for PEAKSET in ${PEAKsets}; do
$BEDTOOLS intersect -b ${DIFFPEAKS}/SA2vsSA1.${PEAKSET}.bed  -a $TSSpos -u > ${DIFFPEAKS}/TSSOverlap/SA2vsSA1.${PEAKSET}.TSS.bed
done
#check number of TSS
for PEAKSET in ${PEAKsets}; do
wc -l ${DIFFPEAKS}/SA2vsSA1.${PEAKSET}.bed
wc -l ${DIFFPEAKS}/TSSOverlap/SA2vsSA1.${PEAKSET}.TSS.bed
done
#------------------------------------------------------------------------------------ 
