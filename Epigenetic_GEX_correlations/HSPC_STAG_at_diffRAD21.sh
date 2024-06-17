#!/bin/bash
# by Alexander Fischer, APR 2024


                         ##############################################################
                         #  STAG1/STAG2 signals at differential RAD21 peaks in HSPCs  #
                         ##############################################################
#####NOTE: for differential regulatory loop anchors referr to HiC_HSPC.diff.LoopAnchors.sh

# setting basic path and linux OS
TMPDIR="/loctmp"
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


# general directories

TAGDIR="${DIR_DATA}/processedData/tagDir/chromatin/GRCh38/ChIP/CD34/siRNA_KD"
INPUTDIR="${DIR_DATA}/processedData/tagDir/DNA/GRCh38/Input/CD34"
PROJECTDIR="${DIR_DATA}/analysis/project_cohesin/CD34"
WORKDIR="${PROJECTDIR}/ChIP_KD_analysis/Cohesin_CTCF_MED12"
PEAKDIR="${WORKDIR}/peaks"
FIGURESDIR="${WORKDIR}/figures"
MERGEDBW="${WORKDIR}/mergedBigWigs"
MOTIFDIR="${WORKDIR}/motifs"
DIFFPEAKS="${WORKDIR}/diffPeaks/"



            #####################################################################
            #         Coverage of STAG1/STAG2 at differential RAD21 peaks       #
            #####################################################################



mkdir ${FIGURESDIR}/ghist/SA1vsSA2/RAD21diffpeaks/

##annotate all types of signals also to all AML SE
filts="2foldDown 2foldUp"
conds="CTRL SA2KD SA1KD"
for filt in $filts;do
for cond in $conds;do
annotatePeaks.pl ${DIFFPEAKS}/qstat_SA2KDvsCTRL.RAD21.Peaks_edgeR.${filt}.txt hg38 -size 4000 -hist 25 -ghist -bedGraph ${MERGEDBW}/merged.ChIP_${cond}_SA1.bedGraph > "${FIGURESDIR}/ghist/SA1vsSA2/RAD21diffpeaks/${cond}_SA2KD.RAD21diff.${filt}.SA1.ghist.txt"
annotatePeaks.pl ${DIFFPEAKS}/qstat_SA2KDvsCTRL.RAD21.Peaks_edgeR.${filt}.txt hg38 -size 4000 -hist 25 -ghist -bedGraph ${MERGEDBW}/merged.ChIP_${cond}_SA2.bedGraph > "${FIGURESDIR}/ghist/SA1vsSA2/RAD21diffpeaks/${cond}_SA2KD.RAD21diff.${filt}.SA2.ghist.txt"
annotatePeaks.pl ${DIFFPEAKS}/qstat_SA1KDvsCTRL.RAD21.Peaks_edgeR.${filt}.txt hg38 -size 4000 -hist 25 -ghist -bedGraph ${MERGEDBW}/merged.ChIP_${cond}_SA1.bedGraph > "${FIGURESDIR}/ghist/SA1vsSA2/RAD21diffpeaks/${cond}_SA1KD.RAD21diff.${filt}.SA1.ghist.txt"
annotatePeaks.pl ${DIFFPEAKS}/qstat_SA1KDvsCTRL.RAD21.Peaks_edgeR.${filt}.txt hg38 -size 4000 -hist 25 -ghist -bedGraph ${MERGEDBW}/merged.ChIP_${cond}_SA2.bedGraph > "${FIGURESDIR}/ghist/SA1vsSA2/RAD21diffpeaks/${cond}_SA1KD.RAD21diff.${filt}.SA2.ghist.txt"
done done

for filt in $filts;do
for cond in $conds;do
plotHIST.sh -g "${cond}_SA2KD.RAD21diff.${filt}.SA1.ghist.txt ${cond}_SA2KD.RAD21diff.${filt}.SA2.ghist.txt" \
-s "STAG1-ChIP STAG2-ChIP" -c "darkgoldenrod1 seagreen1" -x 1000 -y "0 3" -d ${FIGURESDIR}/ghist/SA1vsSA2/RAD21diffpeaks -n ${cond}_SA2KD.RAD21diff.${filt}_SA1vsSA2
done done


for filt in $filts;do
for cond in $conds;do
plotHIST.sh -g "${cond}_SA1KD.RAD21diff.${filt}.SA1.ghist.txt ${cond}_SA1KD.RAD21diff.${filt}.SA2.ghist.txt" \
-s "STAG1-ChIP STAG2-ChIP" -c "darkgoldenrod1 seagreen1" -x 1000 -y "0 3" -d ${FIGURESDIR}/ghist/SA1vsSA2/RAD21diffpeaks -n ${cond}_SA1KD.RAD21diff.${filt}_SA1vsSA2
done done


for cond in $conds;do
plotHIST.sh -g "${cond}_SA1KD.RAD21diff.2foldUp.SA1.ghist.txt ${cond}_SA1KD.RAD21diff.2foldUp.SA2.ghist.txt" \
-s "STAG1-ChIP STAG2-ChIP" -c "darkgoldenrod1 seagreen1" -x 1000 -y "0 5" -d ${FIGURESDIR}/ghist/SA1vsSA2/RAD21diffpeaks -n ${cond}_SA1KD.RAD21diff.2foldUp_SA1vsSA2.max5
done

