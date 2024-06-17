#!/bin/bash
# by Alexander Fischer, APR 2024


                         ######################################################################
                         #  STAG1/STAG2 signals at differential loci (DEGs/Diff RAD) in AMLs  #
                         ######################################################################
###coverage histograms of STAG1 vs STAG2 at sites with DEG TSS, differential RAD21 peaks, 
#####NOTE: for differential regulatory loop anchors referr to HiC_AML.diff.LoopAnchors.sh

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

# directories for analysis
PROJECTDIR="${DIR_DATA}/analysis/project_cohesin"
WORKDIR="${PROJECTDIR}/Cohesin_AML/ChIP_analysis"
PEAKDIR="${WORKDIR}/peaks"
DIFFPEAKS="${WORKDIR}/diffPeaks"
FIGURESDIR="${WORKDIR}/figures"
MERGEDBW="${WORKDIR}/mergedBigWigs"
MERGEDBW_AML_H3K="${WORKDIR}/H3K27ac/mergedBigWigs"
MERGEDBW_AML_ATAC="${PROJECTDIR}/Cohesin_AML/ATAC/mergedBigWigs"
DTmatrixdir=${WORKDIR}/deeptoolsMatrix


## DEG TSS pos
RNADEGDIR="${PROJECTDIR}/Cohesin_AML/RNAseq/Analysis/Resulttables/RNAseq_Cohesin_AML_AF3/STAG2pat_vs_CTRL"



                         ###################################################################################
                         # Histogram plots of STAG1/STAG2 signals at TSS of differentially expressed genes #
                         #                      split by overcompensation status of RAD21                  #
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
#####using higher ylim
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


                         ###################################################################################
                         # Histogram plots of STAG1/STAG2 signals at differential RAD21 positions          #
                         #                           split by CTCF presence                                #
                         ###################################################################################


#subset by CTCF presence
mkdir -p ${FIGURESDIR}/ghist/SA1vsSA2/RAD21diffpeaks/subsetbyCTCF/
CTCFcats="withCTCF noCTCF"
groups="CTRL SA2mut"
filtsRAD="2folddown 2foldup"

for SET in ${CTCFcats};do
for grp in ${groups};do
for FILT in ${filtsRAD};do
	cat >"${TMPDIR}/ghist.${SET}.${grp}.${FILT}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd /loctmp
annotatePeaks.pl ${DIFFPEAKS}/subsetByCTCF/SA2mutvsCTRL.RAD21.stringentPeaks_DESEQ.model.${FILT}.${SET}.txt hg38 -size 4000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_${grp}_SA1.CNVnorm.topQC.bedGraph  > "${FIGURESDIR}/ghist/SA1vsSA2/RAD21diffpeaks/subsetbyCTCF/RAD21diffpeaks.SA2mutvsCTRL.str.RAD21.${SET}.${FILT}.${grp}.av.SA1.ghist.txt"
annotatePeaks.pl ${DIFFPEAKS}/subsetByCTCF/SA2mutvsCTRL.RAD21.stringentPeaks_DESEQ.model.${FILT}.${SET}.txt hg38 -size 4000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_${grp}_SA2.CNVnorm.bedGraph  > "${FIGURESDIR}/ghist/SA1vsSA2/RAD21diffpeaks/subsetbyCTCF/RAD21diffpeaks.SA2mutvsCTRL.str.RAD21.${SET}.${FILT}.${grp}.av.SA2.ghist.txt"
EOF
	chmod 750 "${TMPDIR}/ghist.${SET}.${grp}.${FILT}.${_DATE}.sh"
	echo "plotting ghists for ${SET} ${FILT} ${grp}"
	screen -dm -S ${SET}.${grp}.${FILT} bash -c "bash ${TMPDIR}/ghist.${SET}.${grp}.${FILT}.${_DATE}.sh"
done done done

#STAG2 mut diff peaks separately in CTRL and SA2mut (with 3 or 5 as max value)
for grp in ${groups};do
for SET in ${CTCFcats};do
for FILT in ${filtsRAD};do
plotHIST.sh -g "RAD21diffpeaks.SA2mutvsCTRL.str.RAD21.${SET}.${FILT}.${grp}.av.SA1.ghist.txt RAD21diffpeaks.SA2mutvsCTRL.str.RAD21.${SET}.${FILT}.${grp}.av.SA2.ghist.txt" \
-s "STAG1-ChIP STAG2-ChIP" -c "darkgoldenrod1 seagreen1" -x 1000 -y "0 5" -d ${FIGURESDIR}/ghist/SA1vsSA2/RAD21diffpeaks/subsetbyCTCF -n ${grp}_AML_RAD21diffpeaks.SA2mutvsCTRL.${SET}.${FILT}.SA1vsSA2.max5
plotHIST.sh -g "RAD21diffpeaks.SA2mutvsCTRL.str.RAD21.${SET}.${FILT}.${grp}.av.SA1.ghist.txt RAD21diffpeaks.SA2mutvsCTRL.str.RAD21.${SET}.${FILT}.${grp}.av.SA2.ghist.txt" \
-s "STAG1-ChIP STAG2-ChIP" -c "darkgoldenrod1 seagreen1" -x 1000 -y "0 3" -d ${FIGURESDIR}/ghist/SA1vsSA2/RAD21diffpeaks/subsetbyCTCF -n ${grp}_AML_RAD21diffpeaks.SA2mutvsCTRL.${SET}.${FILT}.SA1vsSA2
done done done
