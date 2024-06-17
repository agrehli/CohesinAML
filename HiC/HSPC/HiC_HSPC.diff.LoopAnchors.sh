#!/bin/bash
#by Alexander Fischer, Jan 2022

#################################################################################
#         Analysis of loop anchors in Cohesin KD HSPCs                          #
#################################################################################

#setting homer environment
DIR_PKG="/misc/software/ngs"
DIR_SOFT="/misc/software"
DIR_PKG="${DIR_SOFT}/ngs"
DIR_DATA="/misc/data"
OS=$(lsb_release -c |grep "^Codename" | awk -F: '{print $2}' | sed 's/[[:blank:]]//g')
PATH_PERL=/misc/software/package/perl/perl-5.26.1/bin
PATH_SAMTOOLS=${DIR_PKG}/samtools/samtools-1.6/bin
PATH_HOMER=${DIR_PKG}/homer/v4.11/bin
PATH_R=/misc/software/package/RBioC/3.4.3/bin
export PATH=${PATH_R}:${PATH_PERL}:${PATH_SAMTOOLS}:${PATH_HOMER}:${PATH}
export PATH
#software versions and general files
BEDTOOLS="${DIR_PKG}/bedtools/bedtools2-2.27.1/bin/bedtools"
BOWTIE2="${DIR_PKG}/bowtie/bowtie2-2.3.4-linux-x86_64/bowtie2"
CHROMSIZES_HG38="${DIR_SOFT}/viewer/IGV/IGVTools_2.3.98/genomes/GRCh38.PRI_p10.chrom.sizes"
BAD_HG38="${DIR_DATA}/analysis/generalStuff/annotation/GRCh38/hg38.badRegions.bed"


#directories and files
TMPDIR="/loctmp"
MAPHICDIR_CD34="${DIR_DATA}/processedData/mapping/3Dchromatin/HiC/GRCh38/CD34"
TAGHICDIR_CD34="${DIR_DATA}/processedData/tagDir/3Dchromatin/HiC/GRCh38/CD34"
PROJDIR="${DIR_DATA}/analysis/project_cohesin"
WORKDIR_CD34="${PROJDIR}/CD34/HiC_KDs/DeepSeq_Analysis"
LOOPDIR="${WORKDIR_CD34}/loops"
LOOPDIR_CD34="${LOOPDIR}/withoutXY"
MATRIXDIR_CD34="${WORKDIR_CD34}/Interactionmatrices"
FIGURESDIR="${WORKDIR_CD34}/figures"
RAWHICDATA_CD34="${DIR_DATA}/rawData/3Dchromatin/HiC/CD34"
diffanchorsn2t="${LOOPDIR_CD34}/differentialanchorsn2t"
PairedAnch=${diffanchorsn2t}/AnchorsOverlapPairs
overlapsplitdir="${diffanchorsn2t}/splitbyoverlap"
DTmatrixdir=${FIGURESDIR}/Heatmaps/DeeptoolsMatrices

mkdir ${FIGURESDIR}/Heatmaps/
mkdir ${DTmatrixdir}

#other data types for integrative analyses
#TF ChIP
CHIPTAGDIR="${DIR_DATA}/processedData/tagDir/chromatin/GRCh38/ChIP/CD34"
WORKDIR_CHIP="${PROJDIR}/CD34/ChIP_KD_analysis/Cohesin_CTCF_MED12"
CHIPPEAKDIR="${WORKDIR_CHIP}/peaks"
PEAKDIRH3K="${PROJDIR}/CD34/ChIP_KD_analysis/H3K27ac/peaks"
CHIPDIFFDIR="${WORKDIR_CHIP}/diffPeaks"
CHIPMERGEBW="${WORKDIR_CHIP}/mergedBigWigs"
CHIPSCALEDBWR="${WORKDIR_CHIP}/RAD21/scaledBigWigs"
RAD21AVEBW="${WORKDIR_CHIP}/RAD21/averageBigWigs"
CHIPSCALEDBWC="${WORKDIR_CHIP}/CTCF/scaledBigWigs"
CHIPSCALEDBWM="${WORKDIR_CHIP}/MED12/scaledBigWigs"
CHIPSCALEDBWP="${WORKDIR_CHIP}/PU1/scaledBigWigs"
#Histone ChIP
H3K27acSCALEDBW="${DIR_DATA}/analysis/project_cohesin/CD34/ChIP_KD_analysis/H3K27ac/scaledBigWigs"
H3K27acMERGEDBW="${DIR_DATA}/analysis/project_cohesin/CD34/ChIP_KD_analysis/H3K27ac/mergedBigWigs"
#ATAC
WORKDIR_ATAC="${DIR_DATA}/analysis/project_cohesin/CD34/ATAC"
MOTIFDIR_ATAC="${WORKDIR_ATAC}/motifs"
SCALEDBW_ATAC="${WORKDIR_ATAC}/scaledBigWigs"
ATACSCALEDBW="${DIR_DATA}/analysis/project_cohesin/CD34/ATAC/scaledBigWigs"
FIGURESDIR_ATAC="${WORKDIR_ATAC}/figures"
ATAC_DIFFDIR="${WORKDIR_ATAC}/diffPeaksPeaks"
PEAKDIR_ATAC="${WORKDIR_ATAC}/peaks"
#gene expression
RNAsResults="${DIR_DATA}/analysis/project_cohesin/CD34/RNAseq/"
#general TSS pos from FANTOM consortium cage seq data
TSSpos="${PROJDIR}/TSS.bed"
#AML diff. loop anchors peak coordinates
AMLdiffanchorsH3KEP=${PROJDIR}/Cohesin_AML/HiC/loops/differentialanchorsn2t/H3K27overlap/
AMLdiffanchorsEP=${AMLdiffanchorsH3KEP}/TSSanchorOverlap/peaksSplitByAnchorTSSstatus

#Define sample names
HIC_CTRL='HiC_CD34_14_3_siCtrl HiC_CD34_17_3_siCtrl HiC_CD34_18_4_siCtrl HiC_CD34_20_6_siCtrl_Rep1 HiC_CD34_21_4_siCtrl_Rep1 HiC_CD34_22_3_siCtrl HiC_CD34_27_4_siCtrl HiC_CD34_28_6_siCtrl ' 
HIC_SA1KD='HiC_CD34_14_1_SA1_KD HiC_CD34_17_1_SA1_KD HiC_CD34_20_4_SA1_KD HiC_CD34_21_2_SA1_KD HiC_CD34_27_3_SA1_KD HiC_CD34_28_4_SA1_KD '
HIC_SA2KD='HiC_CD34_14_2_SA2_KD HiC_CD34_17_2_SA2_KD HiC_CD34_20_5_SA2_KD HiC_CD34_21_3_SA2_KD HiC_CD34_22_2_SA2_KD HiC_CD34_28_5_SA2_KD '
HIC_RAD21KD='HiC_CD34_18_1_RAD21_KD HiC_CD34_20_1_RAD21_KD HiC_CD34_22_1_RAD21_KD HiC_CD34_27_1_RAD21_KD HiC_CD34_28_1_RAD21_KD '
HIC_CD34=$HIC_CTRL$HIC_SA1KD$HIC_SA2KD$HIC_RAD21KD

#Define important variables
KDs="SA1KD SA2KD RAD21KD"
STAGKDs="SA2KD SA1KD"
SETS="SA1KD SA2KD RAD21KD"
Filts="up.FC1 down.FC1 up.FC2 down.FC2"
FiltsFC1="up.FC1 down.FC1"



                   ###########################################################################################
                   #           splitting loop anchors top FC loops from DESEQn2t-differential analysis       #
                   ###########################################################################################

# "differential"/top changed loops
# split loop anchors of differential loops
mkdir $diffanchorsn2t
declare -a DIFFS=("SA2KDvsCTRL" "SA1KDvsCTRL" "RAD21KDvsCTRL")

for DIFF in "${DIFFS[@]}";do
cd ${LOOPDIR_CD34}
#SA2 and SA1 KD vs CTRL (Loops homer normalized)
#sort by FC
sort -g -k9,9 DESEQn2total/${DIFF}.merged.loop.scores.DESeq.norm2Total.loopcoords2.txt | tail -n +2 > ${diffanchorsn2t}/tmp.loopn2t.txt
#filter for logFC>0.558 <-0.558 && padj < 0.05, print anchor point start and stop coordinates plus/minus 5000 bases #DSEQ nt2 data
awk 'BEGIN {FS=OFS="\t"} {if( $9>0.558 && $13<.05 ) print $2, $3-5000, $4+5000, $1"_plus", "1", "+" ; if($9>0.558 && $13<.05) print $5, $6-5000, $7+5000, $1"_minus", "1", "-"}' ${diffanchorsn2t}/tmp.loopn2t.txt > ${diffanchorsn2t}/${DIFF}.loopAnchors.up.FC1sig.n2t.bed
awk 'BEGIN {FS=OFS="\t"} {if( $9<-0.558 && $13<.05 ) print $2, $3-5000, $4+5000, $1"_plus", "1", "+" ; if( $9<-0.558 && $13<.05) print $5, $6-5000, $7+5000, $1"_minus", "1", "-"}' ${diffanchorsn2t}/tmp.loopn2t.txt > ${diffanchorsn2t}/${DIFF}.loopAnchors.down.FC1sig.n2t.bed
#filter for logFC>1 <-1&& padj < 0.05, print anchor point start and stop coordinates plus/minus 5000 bases #DSEQ nt2 data
awk 'BEGIN {FS=OFS="\t"} {if( $9>1 && $13<.05 ) print $2, $3-5000, $4+5000, $1"_plus", "1", "+" ; if($9>1 && $13<.05) print $5, $6-5000, $7+5000, $1"_minus", "1", "-"}' ${diffanchorsn2t}/tmp.loopn2t.txt > ${diffanchorsn2t}/${DIFF}.loopAnchors.up.FC2sig.n2t.bed
awk 'BEGIN {FS=OFS="\t"} {if( $9<-1 && $13<.05 ) print $2, $3-5000, $4+5000, $1"_plus", "1", "+" ; if( $9<-1 && $13<.05) print $5, $6-5000, $7+5000, $1"_minus", "1", "-"}' ${diffanchorsn2t}/tmp.loopn2t.txt > ${diffanchorsn2t}/${DIFF}.loopAnchors.down.FC2sig.n2t.bed
#all loops no filter
awk 'BEGIN {FS=OFS="\t"} { print $2, $3-5000, $4+5000, $1"_plus", "1", "+" ; print $5, $6-5000, $7+5000, $1"_minus", "1", "-"}' ${diffanchorsn2t}/tmp.loopn2t.txt > ${diffanchorsn2t}/tmp.bed
head -n -2 ${diffanchorsn2t}/tmp.bed > ${diffanchorsn2t}/tmp.txt > ${diffanchorsn2t}/${DIFF}.loopAnchors.n2t.all.bed
done

# remove anchors that overlap between up and down (keep as separate class = "*.overlap.*")
FCs="FC1 FC2"

cd ${diffanchorsn2t}
for DIFF in "${DIFFS[@]}";do
for FC in ${FCs};do
$BEDTOOLS merge -i <(sort -k1,1 -k2,2n ${DIFF}.loopAnchors.up.${FC}sig.n2t.bed) -s -d -7500 -c 4,5,6 -o distinct,distinct,distinct > tmp.sorted.up.bed
$BEDTOOLS merge -i <(sort -k1,1 -k2,2n ${DIFF}.loopAnchors.down.${FC}sig.n2t.bed) -s -d -7500  -c 4,5,6 -o distinct,distinct,distinct > tmp.sorted.down.bed
$BEDTOOLS intersect -a  tmp.sorted.up.bed -b ${DIFF}.loopAnchors.down.${FC}sig.n2t.bed -v -f 0.75 > ${DIFF}.loopAnchors.red.up.${FC}sig.n2t.bed
$BEDTOOLS intersect -a  tmp.sorted.down.bed -b ${DIFF}.loopAnchors.up.${FC}sig.n2t.bed -v -f 0.75 > ${DIFF}.loopAnchors.red.down.${FC}sig.n2t.bed
$BEDTOOLS intersect -a  tmp.sorted.up.bed -b tmp.sorted.down.bed -u -f 0.75 > tmp.loop1.txt
$BEDTOOLS intersect -b  tmp.sorted.up.bed -a tmp.sorted.down.bed -u -f 0.75 > tmp.loop2.txt
cat tmp.loop1.txt tmp.loop2.txt > tmp.loop3.txt
$BEDTOOLS merge -i <(sort -k1,1 -k2,2n tmp.loop3.txt) -s -c 4,5,6 -o distinct,distinct,distinct > ${DIFF}.loopAnchors.overlap.${FC}.n2t.bed
$BEDTOOLS intersect -a  tmp.sorted.up.bed -b ${DIFF}.loopAnchors.down.${FC}sig.n2t.bed -v -f 0.75 -s > ${DIFF}.loopAnchors.red.up.${FC}sig.n2t.stranded.bed 
$BEDTOOLS intersect -a  tmp.sorted.down.bed -b ${DIFF}.loopAnchors.up.${FC}sig.n2t.bed -v -f 0.75 -s > ${DIFF}.loopAnchors.red.down.${FC}sig.n2t.stranded.bed 
$BEDTOOLS intersect -a  tmp.sorted.up.bed -b tmp.sorted.down.bed -u -f 0.75 -s > tmp.loop1.txt
$BEDTOOLS intersect -b  tmp.sorted.up.bed -a tmp.sorted.down.bed -u -f 0.75 -s > tmp.loop2.txt
cat tmp.loop1.txt tmp.loop2.txt > tmp.loop3.txt
$BEDTOOLS merge -i <(sort -k1,1 -k2,2n tmp.loop3.txt) -s -c 4,5,6 -o distinct,distinct,distinct > ${DIFF}.loopAnchors.overlap.${FC}.n2t.overlap.stranded.bed
done
done

# convert bed to pos
cd ${diffanchorsn2t}
for DIFF in "${DIFFS[@]}";do
for FC in ${FCs};do
bed2pos.pl ${DIFF}.loopAnchors.up.${FC}sig.n2t.bed > ${DIFF}.loopAnchors.up.${FC}sig.n2t.txt
bed2pos.pl ${DIFF}.loopAnchors.down.${FC}sig.n2t.bed > ${DIFF}.loopAnchors.down.${FC}sig.n2t.txt
bed2pos.pl ${DIFF}.loopAnchors.n2t.all.bed >${DIFF}.loopAnchors.n2t.all.txt
$BEDTOOLS merge -i <(sort -k1,1 -k2,2n ${DIFF}.loopAnchors.n2t.all.bed) -s -d -7500 -c 4,5,6 -o distinct,distinct,distinct > tmp.sorted.all.bed
bed2pos.pl tmp.sorted.all.bed > tmp.${DIFF}.sorted.all.txt
done
done

                   ################################################################
                   #  insulation and DI ratios at differential loop anchors       #
                   ################################################################
mkdir ${FIGURESDIR}/ghist/
mkdir ${FIGURESDIR}/ghist/sigLoopAnchors/
###annotate insulation and directionality bedgraphs to differntial loop anchor regions
_DATE=$(date +%s)
for SET in ${SETS};
do
for filt in ${Filts};
do
	cat >"${TMPDIR}/ghist.${SET}${filt}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd /loctmp
annotatePeaks.pl ${diffanchorsn2t}/${SET}vsCTRL.loopAnchors.${filt}sig.n2t.txt hg38 -size 200000 -hist 250 -bedGraph "${LOOPDIR}/HiC_${SET}_combined.2.5K.10K.Insulation.bedGraph" > "${FIGURESDIR}/ghist/sigLoopAnchors/${SET}.loopAnchors.${SET}.${filt}.Insulation.hist.txt"
annotatePeaks.pl ${diffanchorsn2t}/${SET}vsCTRL.loopAnchors.${filt}sig.n2t.txt hg38 -size 200000 -hist 250 -bedGraph "${LOOPDIR}/HiC_CD34_siCtrl_combined.2.5K.10K.Insulation.bedGraph" > "${FIGURESDIR}/ghist/sigLoopAnchors/CTRL.loopAnchors.${SET}.${filt}.Insulation.hist.txt"
annotatePeaks.pl ${diffanchorsn2t}/${SET}vsCTRL.loopAnchors.${filt}sig.n2t.txt hg38 -size 200000 -hist 250 -bedGraph "${LOOPDIR}/HiC_${SET}_combined.2.5K.10K.Insulation.bedGraph" > "${FIGURESDIR}/ghist/sigLoopAnchors/${SET}.loopAnchors.${SET}.${filt}.Insulation.hist.txt"
annotatePeaks.pl ${diffanchorsn2t}/${SET}vsCTRL.loopAnchors.${filt}sig.n2t.txt hg38 -size 200000 -hist 250 -bedGraph "${LOOPDIR}/HiC_CD34_siCtrl_combined.2.5K.10K.Insulation.bedGraph" > "${FIGURESDIR}/ghist/sigLoopAnchors/CTRL.loopAnchors.${SET}.${filt}.Insulation.hist.txt"
annotatePeaks.pl ${diffanchorsn2t}/${SET}vsCTRL.loopAnchors.${filt}sig.n2t.txt hg38 -size 200000 -hist 250 -bedGraph "${LOOPDIR}/HiC_${SET}_combined.2.5K.10K.DI.bedGraph" > "${FIGURESDIR}/ghist/sigLoopAnchors/${SET}.loopAnchors.${SET}.${filt}.DI.hist.txt"
annotatePeaks.pl ${diffanchorsn2t}/${SET}vsCTRL.loopAnchors.${filt}sig.n2t.txt hg38 -size 200000 -hist 250 -bedGraph "${LOOPDIR}/HiC_CD34_siCtrl_combined.2.5K.10K.DI.bedGraph" > "${FIGURESDIR}/ghist/sigLoopAnchors/CTRL.loopAnchors.${SET}.${filt}.DI.hist.txt"
annotatePeaks.pl ${diffanchorsn2t}/${SET}vsCTRL.loopAnchors.${filt}sig.n2t.txt hg38 -size 200000 -hist 250 -bedGraph "${LOOPDIR}/HiC_${SET}_combined.2.5K.10K.DI.bedGraph" > "${FIGURESDIR}/ghist/sigLoopAnchors/${SET}.loopAnchors.${SET}.${filt}.DI.hist.txt"
annotatePeaks.pl ${diffanchorsn2t}/${SET}vsCTRL.loopAnchors.${filt}sig.n2t.txt hg38 -size 200000 -hist 250 -bedGraph "${LOOPDIR}/HiC_CD34_siCtrl_combined.2.5K.10K.DI.bedGraph" > "${FIGURESDIR}/ghist/sigLoopAnchors/CTRL.loopAnchors.${SET}.${filt}.DI.hist.txt"
EOF
	chmod 750 "${TMPDIR}/ghist.${SET}${filt}.${_DATE}.sh"
	echo "plotting ghists for ${SET}${filt}"
	screen -dm -S ${SET}${filt} bash -c "bash ${TMPDIR}/ghist.${SET}${filt}.${_DATE}.sh"
done
done


# plot histograms comparison of siCTRL and KDs
mkdir ${FIGURESDIR}/Insulation_DI_plots/

#paste histograms of insulation and DI into one table and then into one plot
for SET in ${SETS};
do
for filt in ${Filts};
do
cd ${FIGURESDIR}/ghist/sigLoopAnchors/
paste <(cut -f 1,2 ${SET}.loopAnchors.${SET}.${filt}.Insulation.hist.txt) <(cut -f 2 CTRL.loopAnchors.${SET}.${filt}.Insulation.hist.txt) \
<(cut -f 2 ${SET}.loopAnchors.${SET}.${filt}.DI.hist.txt) <(cut -f 2 CTRL.loopAnchors.${SET}.${filt}.DI.hist.txt) | tail -n +2 > ${TMPDIR}/tmp.txt
echo $'Dist\tKD\tsiCTRL\tKD_DI\tsiCTRL_DI' | cat - ${TMPDIR}/tmp.txt > tmp.${SET}.${filt}.txt
done
done

for filt in ${Filts};
do
_DATE=$(date +%s)
cat >"${TMPDIR}/R.plot.P.${_DATE}.R" <<EOF
dSA1 <- read.delim("${FIGURESDIR}/ghist/sigLoopAnchors/tmp.SA1KD.${filt}.txt", header=T)
dSA2 <- read.delim("${FIGURESDIR}/ghist/sigLoopAnchors/tmp.SA2KD.${filt}.txt", header=T)
dRAD21 <- read.delim("${FIGURESDIR}/ghist/sigLoopAnchors/tmp.RAD21KD.${filt}.txt", header=T)

pdf(file="${FIGURESDIR}/Insulation_DI_plots/SA1vsCTRL.${filt}.loopAnchors.DI.I.hist.pdf", height=8, width=8)
plot(dSA1\$Dist, dSA1\$siCTRL, type="l", col="darkred", lwd=5, xlab="loop anchor center", ylab="DI & Insulation", ylim=c(-.25,.5))
lines(dSA1\$Dist, dSA1\$KD, col="darkgoldenrod2", lwd=5)
lines(dSA1\$Dist, dSA1\$siCTRL_DI, col="firebrick", lwd=5, lty=3)
lines(dSA1\$Dist, dSA1\$KD_DI, col="darkgoldenrod2", lwd=5, lty=3)
abline(h = 0, v = 0, col = "gray60")
dev.off()


pdf(file="${FIGURESDIR}/Insulation_DI_plots/SA2vsCTRL.${filt}.loopAnchors.DI.I.hist.pdf", height=8, width=8)
plot(dSA2\$Dist, dSA2\$siCTRL, type="l", col="darkred", lwd=5, xlab="loop anchor center", ylab="DI & Insulation", ylim=c(-.25,.5))
lines(dSA2\$Dist, dSA2\$KD, col="seagreen", lwd=5)
lines(dSA2\$Dist, dSA2\$siCTRL_DI, col="firebrick", lwd=5, lty=3)
lines(dSA2\$Dist, dSA2\$KD_DI, col="seagreen", lwd=5, lty=3)
abline(h = 0, v = 0, col = "gray60")
dev.off()


pdf(file="${FIGURESDIR}/Insulation_DI_plots/RAD21vsCTRL.${filt}.loopAnchors.DI.I.hist.pdf", height=8, width=8)
plot(dRAD21\$Dist, dRAD21\$siCTRL, type="l", col="darkred", lwd=5, xlab="loop anchor center", ylab="DI & Insulation", ylim=c(-.25,.5))
lines(dRAD21\$Dist, dRAD21\$KD, col="mediumorchid", lwd=5)
lines(dRAD21\$Dist, dRAD21\$siCTRL_DI, col="firebrick", lwd=5, lty=3)
lines(dRAD21\$Dist, dRAD21\$KD_DI, col="mediumorchid", lwd=5, lty=3)
abline(h = 0, v = 0, col = "gray60")
dev.off()
EOF
chmod 750 "${TMPDIR}/R.plot.P.${_DATE}.R"
R < ${TMPDIR}/R.plot.P.${_DATE}.R  --no-save
rm ${TMPDIR}/R.plot.P.${_DATE}.R
done



                          ################################################
                          # CTCF,RAD21,MED12 - loop anchor overlap       #
                          ################################################

###overlap of Loopanchors (identified in DESEQ2 norm2 total loop anchor analysis) with ChIP-seq peaks
ChIPs="CTCF RAD21 MED12"

mkdir ${diffanchorsn2t}/ChIPoverlap
cd ${diffanchorsn2t}
for SET in ${SETS};do
for filt in ${Filts};do
for ChIP in ${ChIPs};do
#diff loop anchors with all CHIPseqpeaks (condition + ctrl peaks)
mergePeaks ${CHIPPEAKDIR}/CTRL_${ChIP}.filtered.peaks.txt ${CHIPPEAKDIR}/${SET}_${ChIP}.filtered.peaks.txt -code > tmpPeaks.txt
pos2bed.pl tmpPeaks.txt > tmpPeaks.bed
$BEDTOOLS intersect -b tmpPeaks.bed -a ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed -u > ChIPoverlap/${SET}vsCTRL.loopAnchors.${filt}.${ChIP}over.bed
#all loop anchors intersected with all CHIPseqpeaks
$BEDTOOLS intersect -b tmpPeaks.bed -a ${SET}vsCTRL.loopAnchors.n2t.all.bed -u > ChIPoverlap/${SET}vsCTRL.loopAnchors.n2t.all.${ChIP}over.bed
#diff loops intersected with differential CHIPseqpeaks
pos2bed.pl ${CHIPDIFFDIR}/qstat_${SET}vsCTRL.${ChIP}.Peaks_edgeR.2foldDown.txt > tmpPeaks.bed
$BEDTOOLS intersect -b tmpPeaks.bed -a ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed -u > ChIPoverlap/${SET}vsCTRL.loopAnchors.${filt}.${ChIP}.diff.down.over.bed
pos2bed.pl ${CHIPDIFFDIR}/qstat_${SET}vsCTRL.${ChIP}.Peaks_edgeR.2foldUp.txt > tmpPeaks.bed
$BEDTOOLS intersect -b tmpPeaks.bed -a ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed -u > ChIPoverlap/${SET}vsCTRL.loopAnchors.${filt}.${ChIP}.diff.up.over.bed
done done done

##summarize:
cd ${diffanchorsn2t}/ChIPoverlap/
for SET in ${SETS};
do
for filt in ${Filts};
do
echo "# ${SET} loops ${filt} #" $(wc -l < ${diffanchorsn2t}/${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed)
echo "# loops ${filt} CTCF overlap #" $(wc -l < ${SET}vsCTRL.loopAnchors.${filt}.CTCFover.bed)
echo "# loops ${filt} RAD21 overlap #" $(wc -l < ${SET}vsCTRL.loopAnchors.${filt}.RAD21over.bed)
echo "# loops ${filt} MED12 overlap #" $(wc -l < ${SET}vsCTRL.loopAnchors.${filt}.MED12over.bed)
done
done

#for SA2KD check also PU1 ChIP peaks
cd ${diffanchorsn2t}
for filt in ${Filts};
do
#diff loop anchors with all CHIPseqpeaks
mergePeaks ${CHIPPEAKDIR}/CTRL_PU1.filtered.peaks.txt ${CHIPPEAKDIR}/SA2KD_PU1.filtered.peaks.txt -code > tmpPeaks.txt
pos2bed.pl tmpPeaks.txt > tmpPeaks.bed
$BEDTOOLS intersect -b tmpPeaks.bed -a SA2KDvsCTRL.loopAnchors.red.${filt}sig.n2t.bed -u > ChIPoverlap/SA2KDvsCTRL.loopAnchors.${filt}.PU1over.bed
#all loop anchors intersected with all CHIPseqpeaks
$BEDTOOLS intersect -b tmpPeaks.bed -a SA2KDvsCTRL.loopAnchors.n2t.all.bed -u > ChIPoverlap/SA2KDvsCTRL.loopAnchors.n2t.all.PU1over.bed
done

cd ${diffanchorsn2t}/ChIPoverlap/
for filt in ${Filts};
do
echo "# SA2KD loops ${filt} #" $(wc -l < ${diffanchorsn2t}/SA2KDvsCTRL.loopAnchors.red.${filt}sig.n2t.bed)
echo "# loops ${filt} PU1 overlap #" $(wc -l < SA2KDvsCTRL.loopAnchors.${filt}.PU1over.bed)
done

                     #####################################################
                     #   ChIP peaks  at differential anchors             #
                     #####################################################
#repeat intersections with a and b options exchanged to get coordinates of peaks overapping the differential loop anchor regions
cd ${diffanchorsn2t}
for SET in ${SETS};
do
for filt in ${Filts};
do
 for ChIP in ${ChIPs}; do
pos2bed.pl ${CHIPPEAKDIR}/CD34_${ChIP}.filtered.peaks.txt > CD34_${ChIP}.tmpPeaks.bed 
$BEDTOOLS intersect -a CD34_${ChIP}.tmpPeaks.bed -b ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed -u > ChIPoverlap/${SET}vsCTRL.loopAnchors.${filt}.${ChIP}centered.bed
$BEDTOOLS intersect -a CD34_${ChIP}.tmpPeaks.bed -b ${SET}vsCTRL.loopAnchors.n2t.all.bed -u > ChIPoverlap/${SET}vsCTRL.loopAnchors.n2t.all.${ChIP}centered.bed
bed2pos.pl ChIPoverlap/${SET}vsCTRL.loopAnchors.${filt}.${ChIP}centered.bed > ChIPoverlap/${SET}vsCTRL.loopAnchors.${filt}.${ChIP}centered.pos.txt
bed2pos.pl ChIPoverlap/${SET}vsCTRL.loopAnchors.n2t.all.${ChIP}centered.bed > ChIPoverlap/${SET}vsCTRL.loopAnchors.n2t.all.${ChIP}centered.pos.txt
 done
done
done
rm *.tmpPeaks.bed

#for SA2KD check also PU1 ChIP peaks
for filt in ${Filts};
do
pos2bed.pl ${CHIPPEAKDIR}/CD34_PU1.filtered.peaks.txt > CD34_PU1.tmpPeaks.bed 
$BEDTOOLS intersect -a CD34_PU1.tmpPeaks.bed -b SA2KDvsCTRL.loopAnchors.red.${filt}sig.n2t.bed -u > ChIPoverlap/SA2KDvsCTRL.loopAnchors.${filt}.PU1centered.bed
$BEDTOOLS intersect -a CD34_PU1.tmpPeaks.bed -b SA2KDvsCTRL.loopAnchors.n2t.all.bed -u > ChIPoverlap/SA2KDvsCTRL.loopAnchors.n2t.all.PU1centered.bed
bed2pos.pl ChIPoverlap/SA2KDvsCTRL.loopAnchors.${filt}.PU1centered.bed > ChIPoverlap/SA2KDvsCTRL.loopAnchors.${filt}.PU1centered.pos.txt
bed2pos.pl ChIPoverlap/SA2KDvsCTRL.loopAnchors.n2t.all.PU1centered.bed > ChIPoverlap/SA2KDvsCTRL.loopAnchors.n2t.all.PU1centered.pos.txt
done


##check CTCF and RAD21 co-associated loop anchor overlap
cd ${diffanchorsn2t}/ChIPoverlap/
$BEDTOOLS intersect -a CD34_RAD21.tmpPeaks.bed -b CD34_CTCF.tmpPeaks.bed -u > tmp.RAD21_CTCF_PEAKoverlap.bed #50143 
$BEDTOOLS intersect -a CD34_CTCF.tmpPeaks.bed -b CD34_RAD21.tmpPeaks.bed -v > tmp.CTCF_wo_RAD21_PEAKoverlap.bed #10336
$BEDTOOLS intersect -b CD34_CTCF.tmpPeaks.bed -a CD34_RAD21.tmpPeaks.bed -v > tmp.RAD21_woCTCF_PEAKoverlap.bed #28724 

for SET in ${SETS}; do
for filt in ${Filts}; do
$BEDTOOLS intersect -b tmp.RAD21_CTCF_PEAKoverlap.bed -a ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed   -u > ChIPoverlap/${SET}vsCTRL.loopAnchors.red.${filt}.RAD21andCTCF.overlap.bed
$BEDTOOLS intersect -a tmp.RAD21_CTCF_PEAKoverlap.bed -b ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed   -u > ChIPoverlap/${SET}vsCTRL.loopAnchors.red.${filt}.RAD21andCTCF.centred.bed
$BEDTOOLS intersect -b tmp.CTCF_wo_RAD21_PEAKoverlap.bed -a ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed   -u > ChIPoverlap/${SET}vsCTRL.loopAnchors.red.${filt}.CTCFwoRAD21.overlap.bed
done
done

for SET in ${SETS}; do
for filt in ${FiltsFC1}; do
echo "# ${SET} loop anchors ${filt} with cohesin #" $(wc -l < ChIPoverlap/${SET}vsCTRL.loopAnchors.${filt}.RAD21over.bed)
echo "# ${SET} loop anchors ${filt} with cohesin and CTCF #" $(wc -l < ChIPoverlap/${SET}vsCTRL.loopAnchors.red.${filt}.RAD21andCTCF.overlap.bed)
echo "# ${SET} loop anchors ${filt} with CTCF and no Cohesin #" $(wc -l < ChIPoverlap/${SET}vsCTRL.loopAnchors.red.${filt}.CTCFwoRAD21.overlap.bed)
echo "# ${SET} cohesin with CTCF at loop anchors ${filt} #" $(wc -l < ChIPoverlap/${SET}vsCTRL.loopAnchors.red.${filt}.RAD21andCTCF.centred.bed)
done
done

cd ${diffanchorsn2t}
mkdir ${diffanchorsn2t}/ChIPoverlap/subsetByCTCF/
for filt in ${Filts};do
for KD in ${KDs};do
$BEDTOOLS intersect -a  ChIPoverlap/${KD}vsCTRL.loopAnchors.${filt}.RAD21centered.bed -b CD34_CTCF.tmpPeaks.bed -u > ChIPoverlap/subsetByCTCF/${KD}vsCTRL.loopAnchors.${filt}.RAD21centered.withCTCF.bed
$BEDTOOLS intersect -a  ChIPoverlap/${KD}vsCTRL.loopAnchors.${filt}.RAD21centered.bed -b CD34_CTCF.tmpPeaks.bed -v > ChIPoverlap/subsetByCTCF/${KD}vsCTRL.loopAnchors.${filt}.RAD21centered.noCTCF.bed
done
done

for filt in ${FiltsFC1};do
for KD in ${KDs};do
wc -l ${diffanchorsn2t}/ChIPoverlap/subsetByCTCF/${KD}vsCTRL.loopAnchors.${filt}.RAD21centered.withCTCF.bed
wc -l ${diffanchorsn2t}/ChIPoverlap/subsetByCTCF/${KD}vsCTRL.loopAnchors.${filt}.RAD21centered.noCTCF.bed
cd ${diffanchorsn2t}/ChIPoverlap/subsetByCTCF
bed2pos.pl ${KD}vsCTRL.loopAnchors.${filt}.RAD21centered.withCTCF.bed > ${KD}vsCTRL.loopAnchors.${filt}.RAD21centered.withCTCF.pos.txt
bed2pos.pl ${KD}vsCTRL.loopAnchors.${filt}.RAD21centered.noCTCF.bed > ${KD}vsCTRL.loopAnchors.${filt}.RAD21centered.noCTCF.pos.txt
done
done
                          #########################################################################
                          #                  H3K27ac  at  diff loop anchors                       #
                          #########################################################################

#check which H3K27ac peaks are in diffloops
##allenhancers used for differential H3K27ac analysis
pos2bed.pl ${PEAKDIRH3K}/mergedCD34_H3K27ac.filtered.peaks.txt > ${TMPDIR}/tmp.peaks.bed

##overlap with sig. loop anchor regions ##look from both sides (a/b options)
mkdir ${diffanchorsn2t}/H3K27overlap
cd ${diffanchorsn2t}
for SET in ${SETS}; do
for filt in ${Filts}; do
$BEDTOOLS intersect -a ${TMPDIR}/tmp.peaks.bed -b ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed  -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27accentred.bed
$BEDTOOLS intersect -b ${TMPDIR}/tmp.peaks.bed -a ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed  -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27acOverlap.bed
done
done

for SET in ${SETS}; do
for filt in ${Filts}; do
echo "# ${SET} loop anchors ${filt} with H3K27ac peaks #" $(wc -l < H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27acOverlap.bed)
echo "# ${SET} H3K27ac peaks at loop anchors ${filt} #" $(wc -l < H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27accentred.bed)
echo "# ${SET} loop anchors ${filt} EP contacts #" $(wc -l < H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27acOverlap.bed)/$(wc -l < ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed)
done
done

#cohesin associated enhancers at loop anchors
cd ${diffanchorsn2t}
pos2bed.pl ${PEAKDIRH3K}/mergedCD34_H3K27ac.filtered.peaks.txt > CD34_H3K27ac.tmp.peaks.bed
pos2bed.pl ${CHIPPEAKDIR}/CD34_RAD21.filtered.peaks.txt > CD34_RAD21.tmpPeaks.bed 
pos2bed.pl ${CHIPPEAKDIR}/CD34_CTCF.filtered.peaks.txt > CD34_CTCF.tmpPeaks.bed 
$BEDTOOLS intersect -a CD34_H3K27ac.tmp.peaks.bed -b CD34_RAD21.tmpPeaks.bed  -u > H3K27overlap/CohesinAssEnhancers.bed
$BEDTOOLS intersect -b CD34_H3K27ac.tmp.peaks.bed -a CD34_RAD21.tmpPeaks.bed  -u > H3K27overlap/CohesinAssEnhancers.RAD21centred.bed
###at diff loop anchors
for SET in ${SETS}; do
for filt in ${Filts}; do
$BEDTOOLS intersect -b H3K27overlap/CohesinAssEnhancers.bed -a ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed  -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.Overlap.bed  #Anchors with Cohesin Ass. enhancer
$BEDTOOLS intersect -b H3K27overlap/CohesinAssEnhancers.bed -a ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed  -v > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.Overlap.bed #Anchors without Cohesin Ass. enhancer
$BEDTOOLS intersect -a H3K27overlap/CohesinAssEnhancers.bed -b ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed  -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.centred.bed
done
done

#### structural anchors (=wo CohesinAssEnhancers): 2 way overlaps with RAD21 (stringent) and CTCF
STRUCLOOPDIR=${diffanchorsn2t}/ChIPoverlap/structural_Loops
mkdir ${STRUCLOOPDIR}
cp ${diffanchorsn2t}/H3K27overlap/*.wo.CohesinAssEnhancers.Overlap.bed ${STRUCLOOPDIR}/
cd ${diffanchorsn2t}
for SET in ${SETS}; do
for filt in ${Filts}; do
$BEDTOOLS intersect -a ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.Overlap.bed -b CD34_RAD21.tmpPeaks.bed -wao > ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.structural.RAD21.2wayoverlap.bed
$BEDTOOLS intersect -a ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.Overlap.bed -b CD34_CTCF.tmpPeaks.bed -wao > ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.structural.CTCF.2wayoverlap.bed
$BEDTOOLS intersect -a ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.Overlap.bed -b CD34_H3K27ac.tmp.peaks.bed -wao > ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.structural.H3K27ac.2wayoverlap.bed
done
done

#center structual anchors on RAD21
cd ${diffanchorsn2t}
for SET in ${SETS}; do
for filt in ${Filts}; do
$BEDTOOLS intersect -b ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.Overlap.bed -a CD34_RAD21.tmpPeaks.bed -u > ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.structural.RAD21.centred.bed
done
done

###report both loop and peak ID
cd ${diffanchorsn2t}
mkdir H3K27overlap/2wayOverlap
for SET in ${SETS}; do
for filt in ${FiltsFC1}; do
$BEDTOOLS intersect -b CD34_H3K27ac.tmp.peaks.bed -a ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed  -wao > H3K27overlap/2wayOverlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27ac.2wayOverlap.bed
$BEDTOOLS intersect -b H3K27overlap/CohesinAssEnhancers.bed -a ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed  -wao > H3K27overlap/2wayOverlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.2wayOverlap.bed #Anchors with Cohesin Ass. enhancer
$BEDTOOLS intersect -b H3K27overlap/CohesinAssEnhancers.RAD21centred.bed -a ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed  -wao > H3K27overlap/2wayOverlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.RAD21centred.2wayOverlap.bed #Anchors with Cohesin Ass. enhancer centred on RAD21
done
done

###all loop anchors
for SET in ${SETS}; do
$BEDTOOLS intersect -b H3K27overlap/CohesinAssEnhancers.bed -a ${SET}vsCTRL.loopAnchors.n2t.all.bed -u > H3K27overlap/${SET}vsCTRL.loopAnchors.all.CohesinAssEnhancers.centred.bed
done

###diffloop anchor overlaps without cohesin ass. enhancers
for SET in ${SETS}; do
for filt in ${Filts}; do
##without cohesin ass. enhancers, but cohesin
$BEDTOOLS intersect -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.Overlap.bed -b CD34_RAD21.tmpPeaks.bed -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.RAD21.Overlap.bed
##without cohesin ass. enhancers, without cohesin
$BEDTOOLS intersect -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.Overlap.bed -b CD34_RAD21.tmpPeaks.bed -v > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.wo.RAD21.Overlap.bed
##without cohesin ass. enhancers, but cohesin with CTCF
$BEDTOOLS intersect -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.RAD21.Overlap.bed -b CD34_CTCF.tmpPeaks.bed -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.RAD21.CTCF.Overlap.bed
##without cohesin ass. enhancers, but cohesin without CTCF
$BEDTOOLS intersect -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.RAD21.Overlap.bed -b CD34_CTCF.tmpPeaks.bed -v > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.RAD21.woCTCF.Overlap.bed
##without cohesin ass. enhancers, without cohesin
$BEDTOOLS intersect -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.wo.RAD21.Overlap.bed -b CD34_CTCF.tmpPeaks.bed -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.wo.RAD21.CTCF.Overlap.bed
$BEDTOOLS intersect -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.wo.RAD21.Overlap.bed -b CD34_CTCF.tmpPeaks.bed -v > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.wo.RAD21.woCTCF.Overlap.bed
$BEDTOOLS intersect -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.wo.RAD21.woCTCF.Overlap.bed -b CD34_H3K27ac.tmp.peaks.bed -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.wo.RAD21.woCTCF.H3K27ac.Overlap.bed
$BEDTOOLS intersect -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.wo.RAD21.woCTCF.Overlap.bed -b CD34_H3K27ac.tmp.peaks.bed -v > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.wo.RAD21.woCTCF.woH3K27ac.Overlap.bed
done
done

#cohesin associated enhancers association with CTCF?
for SET in ${SETS}; do
for filt in ${Filts}; do
$BEDTOOLS intersect -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.Overlap.bed -b CD34_CTCF.tmpPeaks.bed -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.CTCF.Overlap.bed
$BEDTOOLS intersect -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.Overlap.bed -b CD34_CTCF.tmpPeaks.bed -v > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.wo.CTCF.Overlap.bed
$BEDTOOLS intersect -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.centred.bed -b CD34_CTCF.tmpPeaks.bed   -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancersCTCF.centred.bed
$BEDTOOLS intersect -b H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancersCTCF.centred.bed -a ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed   -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancersCTCF.overlap.bed
$BEDTOOLS intersect -b H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancersCTCF.centred.bed -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.centred.bed   -v > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.wo.CTCF.overlap.bed
done
done
FiltsFC1="up.FC1 down.FC1"
for SET in ${SETS}; do
for filt in ${FiltsFC1}; do
echo "# ${SET} loop anchors ${filt} with cohesin-H3K27ac #" $(wc -l < H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.Overlap.bed)
echo "# ${SET} loop anchors ${filt} without cohesin-H3K27ac #" $(wc -l < H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.Overlap.bed)
echo "# ${SET} loop anchors ${filt} with cohesin-H3K27ac and CTCF#" $(wc -l < H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.CTCF.Overlap.bed)
echo "# ${SET} loop anchors ${filt} with cohesin-H3K27ac without CTCF#" $(wc -l < H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.wo.CTCF.Overlap.bed)
done
done

##focus on H3K27 centred anchors centred on cohesin --> "Cohesin-dependent-enhancers at loop anchors"
cd ${diffanchorsn2t}
pos2bed.pl ${CHIPPEAKDIR}/CD34_RAD21.filtered.peaks.txt > CD34_RAD21.tmpPeaks.bed 
for SET in ${SETS}; do
for filt in ${Filts}; do
$BEDTOOLS intersect -a CD34_RAD21.tmpPeaks.bed -b H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27accentred.bed -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27accentred.RAD21centred.bed
done
done

for SET in ${SETS}; do
for filt in ${Filts}; do
echo "# ${SET} loop anchors ${filt} with H3K27ac peaks #" $(wc -l < H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27acOverlap.bed)
echo "# ${SET} RAD21 peaks in H3K27ac marks at loop anchors ${filt} #" $(wc -l < H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27accentred.RAD21centred.bed)
done
done


                          ##############################################################
                          #            Overlaps of diff. anchors with TSS              #
                          ##############################################################

mkdir ${diffanchorsn2t}/H3K27overlap/TSSanchorOverlap/
mkdir ${PairedAnch}/TSSoverlap/

#report all regulatory loops: this can be used for enhancer promoter distinction of anchors:
cd ${diffanchorsn2t}/H3K27overlap
for KD in ${KDs}; do
for filt in ${FiltsFC1}; do
$BEDTOOLS intersect -a ${KD}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.Overlap.bed -b ${TSSpos} -wao > TSSanchorOverlap/${KD}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.Overlap.TSS.2wayoverlap.bed
done
done

#also test for structural loops
for SET in ${KDs}; do
for filt in ${FiltsFC1}; do
##report all structural loops:
$BEDTOOLS intersect -a ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.Overlap.bed -b ${TSSpos} -wao > ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.structural.TSS.2wayoverlap.bed
##report only TSS:
$BEDTOOLS intersect -b ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.Overlap.bed -a ${TSSpos} -u > ${diffanchorsn2t}/H3K27overlap/TSSanchorOverlap/${SET}vsCTRL.loopAnchors.red.${filt}.structural.TSS.bed
done
done



                          #########################################################################
                          #      HSPCs    ChIP signals at SA2mut diff loop anchors                #
                          #########################################################################
# ChIPcoverage of RAD21/H3K27ac/STAG split by Enhancer of Promoter ID  (ignoring paired/unpaired status) 

mkdir ${FIGURESDIR}/ghist/AMLanchors/ChIPbyEP
mkdir ${FIGURESDIR}/ghist/AMLanchors/ChIPbyEP/scbw/
mkdir ${FIGURESDIR}/ghist/AMLanchors/ChIPbyEP/STAG/
################################################################################for SA2 mut anchors ################################################################################
################################################################################RAD21/H3K27ac
types="Enhancer Promoter"
for filt in ${FiltsFC1}; do
 for type in ${types}; do
     cat >"${TMPDIR}/ghist.${filt}.${type}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd ${AMLdiffanchorsEP}/
annotatePeaks.pl ${type}_RAD21peaks.${filt}_SA2mut.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${CHIPSCALEDBWR}/ave.CTRL.RAD21.scaled.bedGraph > "${FIGURESDIR}/ghist/AMLanchors/ChIPbyEP/scbw/CTRL.${type}_RAD21peaks.${filt}_SA2mut.RAD21.ghist.txt"
annotatePeaks.pl ${type}_RAD21peaks.${filt}_SA2mut.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${CHIPSCALEDBWR}/ave.SA2KD.RAD21.scaled.bedGraph > "${FIGURESDIR}/ghist/AMLanchors/ChIPbyEP/scbw/SA2KD.${type}_RAD21peaks.${filt}_SA2mut.RAD21.ghist.txt"
annotatePeaks.pl ${type}_RAD21peaks.${filt}_SA2mut.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${CHIPSCALEDBWR}/ave.SA1KD.RAD21.scaled.bedGraph > "${FIGURESDIR}/ghist/AMLanchors/ChIPbyEP/scbw/SA1KD.${type}_RAD21peaks.${filt}_SA2mut.RAD21.ghist.txt"
EOF
    chmod 750 "${TMPDIR}/ghist.${filt}.${type}.${_DATE}.sh"
	echo "plotting ghists for ${type} SA2mut ${filt} loop anchors in HSPCs"
	screen -dm -S ${filt}.${type} bash -c "bash ${TMPDIR}/ghist.${filt}.${type}.${_DATE}.sh"
 done
done


for filt in ${FiltsFC1};do
for type in ${types}; do
plotHIST.sh -g "SA2KD.${type}_RAD21peaks.${filt}_SA2mut.RAD21.ghist.txt CTRL.${type}_RAD21peaks.${filt}_SA2mut.RAD21.ghist.txt" \
-s "STAG2-KD CTRL-HSPCs" -c "seagreen1 firebrick1" -x 1000 -y "0 20" -d ${FIGURESDIR}/ghist/AMLanchors/ChIPbyEP/scbw/ -n RAD21.SA2KDvsCTRL.SA2mut.${type}_RAD21peaks.${filt}

plotHIST.sh -g "SA1KD.${type}_RAD21peaks.${filt}_SA2mut.RAD21.ghist.txt CTRL.${type}_RAD21peaks.${filt}_SA2mut.RAD21.ghist.txt" \
-s "STAG1-KD CTRL-HSPCs" -c "darkgoldenrod3 firebrick1" -x 1000 -y "0 20" -d ${FIGURESDIR}/ghist/AMLanchors/ChIPbyEP/scbw/ -n RAD21.SA1KDvsCTRL.SA2mut.${type}_RAD21peaks.${filt}
done done


################################################################################STAG1/STAG2
#############split by E-P
STAGs="SA1 SA2"
for STAG in ${STAGs}; do
for filt in ${FiltsFC1}; do
 for type in ${types}; do
     cat >"${TMPDIR}/ghist.${filt}.${type}.${STAG}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd ${AMLdiffanchorsEP}/
annotatePeaks.pl ${type}_RAD21peaks.${filt}_SA2mut.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${CHIPMERGEBW}/merged.ChIP_CTRL_${STAG}.bedGraph > "${FIGURESDIR}/ghist/AMLanchors/ChIPbyEP/STAG/CTRL.${type}_RAD21peaks.${filt}_SA2mut.${STAG}.ghist.txt"
annotatePeaks.pl ${type}_RAD21peaks.${filt}_SA2mut.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${CHIPMERGEBW}/merged.ChIP_SA2KD_${STAG}.bedGraph > "${FIGURESDIR}/ghist/AMLanchors/ChIPbyEP/STAG/SA2KD.${type}_RAD21peaks.${filt}_SA2mut.${STAG}.ghist.txt"
annotatePeaks.pl ${type}_RAD21peaks.${filt}_SA2mut.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${CHIPMERGEBW}/merged.ChIP_SA1KD_${STAG}.bedGraph > "${FIGURESDIR}/ghist/AMLanchors/ChIPbyEP/STAG/SA1KD.${type}_RAD21peaks.${filt}_SA2mut.${STAG}.ghist.txt"
EOF
    chmod 750 "${TMPDIR}/ghist.${filt}.${type}.${STAG}.${_DATE}.sh"
	echo "plotting ${STAG} ghists for ${type} SA2mut ${filt} loop anchors in HSPCs"
	screen -dm -S ${filt}.${type}.${STAG} bash -c "bash ${TMPDIR}/ghist.${filt}.${type}.${STAG}.${_DATE}.sh"
 done
done
done
#############combined E-P
for STAG in ${STAGs}; do
for filt in ${FiltsFC1}; do
     cat >"${TMPDIR}/ghist.${filt}.${STAG}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd ${AMLdiffanchorsH3KEP}/
annotatePeaks.pl SA2mutvsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.RAD21centred.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${CHIPMERGEBW}/merged.ChIP_CTRL_${STAG}.bedGraph > "${FIGURESDIR}/ghist/AMLanchors/ChIPbyEP/STAG/CTRL.allEP_RAD21peaks.${filt}_SA2mut.${STAG}.ghist.txt"
annotatePeaks.pl SA2mutvsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.RAD21centred.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${CHIPMERGEBW}/merged.ChIP_SA2KD_${STAG}.bedGraph > "${FIGURESDIR}/ghist/AMLanchors/ChIPbyEP/STAG/SA2KD.allEP_RAD21peaks.${filt}_SA2mut.${STAG}.ghist.txt"
annotatePeaks.pl SA2mutvsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.RAD21centred.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${CHIPMERGEBW}/merged.ChIP_SA1KD_${STAG}.bedGraph > "${FIGURESDIR}/ghist/AMLanchors/ChIPbyEP/STAG/SA1KD.allEP_RAD21peaks.${filt}_SA2mut.${STAG}.ghist.txt"
EOF
    chmod 750 "${TMPDIR}/ghist.${filt}.${STAG}.${_DATE}.sh"
	echo "plotting ${STAG} ghists for all SA2mut ${filt} EP-loop anchors in HSPCs"
	screen -dm -S ${filt}.${type}.${STAG} bash -c "bash ${TMPDIR}/ghist.${filt}.${STAG}.${_DATE}.sh"
 done
done

#############SA1/SA2 plots for SA2mut specific loops by condition
conds="CTRL SA2KD SA1KD"
types2="Enhancer Promoter allEP"
for cond in ${conds};do
for filt in ${FiltsFC1};do
for type in ${types2}; do
plotHIST.sh -g "${cond}.${type}_RAD21peaks.${filt}_SA2mut.SA2.ghist.txt ${cond}.${type}_RAD21peaks.${filt}_SA2mut.SA1.ghist.txt" \
-s "STAG2 STAG1" -c "springgreen2 darkgoldenrod2" -x 1000 -y "0 10" -d ${FIGURESDIR}/ghist/AMLanchors/ChIPbyEP/STAG -n ${cond}.SA1vsSA2.${type}_RAD21peaks.${filt}
done done done



                          #########################################################################
                          #      HSPCs    ChIP signals at KD diff loop anchors                    #
                          #########################################################################

#RAD21/H3K27ac for diff. EP anchors directly found in HSPCs (ignoring paired/unpaired)
mkdir ${FIGURESDIR}/ghist/sigLoopAnchors/ChIPallEP/
mkdir ${FIGURESDIR}/ghist/sigLoopAnchors/ChIPallEP/scbw/
types="Enhancer Promoter"
for filt in ${FiltsFC1}; do
for KD in ${STAGKDs}; do
 for type in ${types}; do
     cat >"${TMPDIR}/ghist.${KD}.${filt}.${type}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd ${diffanchorsn2t}/H3K27overlap/TSSanchorOverlap/peaksSplitByAnchorTSSstatus
annotatePeaks.pl ${type}_RAD21peaks.${filt}_${KD}.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${CHIPSCALEDBWR}/ave.CTRL.RAD21.scaled.bedGraph > "${FIGURESDIR}/ghist/sigLoopAnchors/ChIPallEP/scbw/CTRL.${type}_RAD21peaks.${filt}_${KD}.RAD21.ghist.txt"
annotatePeaks.pl ${type}_RAD21peaks.${filt}_${KD}.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${CHIPSCALEDBWR}/ave.SA2KD.RAD21.scaled.bedGraph > "${FIGURESDIR}/ghist/sigLoopAnchors/ChIPallEP/scbw/SA2KD.${type}_RAD21peaks.${filt}_${KD}.RAD21.ghist.txt"
annotatePeaks.pl ${type}_RAD21peaks.${filt}_${KD}.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${CHIPSCALEDBWR}/ave.SA1KD.RAD21.scaled.bedGraph > "${FIGURESDIR}/ghist/sigLoopAnchors/ChIPallEP/scbw/SA1KD.${type}_RAD21peaks.${filt}_${KD}.RAD21.ghist.txt"
EOF
    chmod 750 "${TMPDIR}/ghist.${KD}.${filt}.${type}.${_DATE}.sh"
	echo "plotting ghists for ${type} in ${KD} ${filt} loop anchors"
	screen -dm -S ${KD}.${filt}.${type} bash -c "bash ${TMPDIR}/ghist.${KD}.${filt}.${type}.${_DATE}.sh"
 done
done
done

for filt in ${FiltsFC1};do
for type in ${types}; do
#STAG2 diff EP
plotHIST.sh -g "SA2KD.${type}_RAD21peaks.${filt}_SA2KD.RAD21.ghist.txt CTRL.${type}_RAD21peaks.${filt}_SA2KD.RAD21.ghist.txt" \
-s "STAG2-KD CTRL-HSPCs" -c "seagreen1 firebrick1" -x 1000 -y "0 30" -d ${FIGURESDIR}/ghist/sigLoopAnchors/ChIPallEP/scbw/ -n RAD21.SA2KDvsCTRL.SA2KD.${type}_RAD21peaks.${filt}.2
plotHIST.sh -g "SA1KD.${type}_RAD21peaks.${filt}_SA2KD.RAD21.ghist.txt CTRL.${type}_RAD21peaks.${filt}_SA2KD.RAD21.ghist.txt" \
-s "STAG1-KD CTRL-HSPCs" -c "darkgoldenrod3 firebrick1" -x 1000 -y "0 30" -d ${FIGURESDIR}/ghist/sigLoopAnchors/ChIPallEP/scbw/ -n RAD21.SA1KDvsCTRL.SA2KD.${type}_RAD21peaks.${filt}.2

#STAG1 diff EP
plotHIST.sh -g "SA2KD.${type}_RAD21peaks.${filt}_SA1KD.RAD21.ghist.txt CTRL.${type}_RAD21peaks.${filt}_SA1KD.RAD21.ghist.txt" \
-s "STAG2-KD CTRL-HSPCs" -c "seagreen1 firebrick1" -x 1000 -y "0 30" -d ${FIGURESDIR}/ghist/sigLoopAnchors/ChIPallEP/scbw/ -n RAD21.SA2KDvsCTRL.SA1KD.${type}_RAD21peaks.${filt}.2
plotHIST.sh -g "SA1KD.${type}_RAD21peaks.${filt}_SA1KD.RAD21.ghist.txt CTRL.${type}_RAD21peaks.${filt}_SA1KD.RAD21.ghist.txt" \
-s "STAG1-KD CTRL-HSPCs" -c "darkgoldenrod3 firebrick1" -x 1000 -y "0 30" -d ${FIGURESDIR}/ghist/sigLoopAnchors/ChIPallEP/scbw/ -n RAD21.SA1KDvsCTRL.SA1KD.${type}_RAD21peaks.${filt}.2
done done

cd ${FIGURESDIR}/ghist/sigLoopAnchors/ChIPallEP/scbw/
mkdir STAG2diffEP
mv *.SA2KD.*_RAD21peaks.*hist.pdf STAG2diffEP/

mkdir STAG1diffEP
mv *.SA1KD.*_RAD21peaks.*hist.pdf STAG1diffEP/

################################################################################STAG1/STAG2
mkdir ${FIGURESDIR}/ghist/sigLoopAnchors/ChIPallEP/STAG

#############split by E-P
STAGs="SA1 SA2"
for KD in ${STAGKDs}; do
for STAG in ${STAGs}; do
for filt in ${FiltsFC1}; do
 for type in ${types}; do
     cat >"${TMPDIR}/ghist.${KD}.${filt}.${type}.${STAG}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd ${diffanchorsn2t}/H3K27overlap/TSSanchorOverlap/peaksSplitByAnchorTSSstatus
annotatePeaks.pl ${type}_RAD21peaks.${filt}_${KD}.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${CHIPMERGEBW}/merged.ChIP_CTRL_${STAG}.bedGraph > "${FIGURESDIR}/ghist/sigLoopAnchors/ChIPallEP/STAG/CTRL.${type}_RAD21peaks.${filt}_${KD}.${STAG}.ghist.txt"
annotatePeaks.pl ${type}_RAD21peaks.${filt}_${KD}.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${CHIPMERGEBW}/merged.ChIP_SA2KD_${STAG}.bedGraph > "${FIGURESDIR}/ghist/sigLoopAnchors/ChIPallEP/STAG/SA2KD.${type}_RAD21peaks.${filt}_${KD}.${STAG}.ghist.txt"
annotatePeaks.pl ${type}_RAD21peaks.${filt}_${KD}.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${CHIPMERGEBW}/merged.ChIP_SA1KD_${STAG}.bedGraph > "${FIGURESDIR}/ghist/sigLoopAnchors/ChIPallEP/STAG/SA1KD.${type}_RAD21peaks.${filt}_${KD}.${STAG}.ghist.txt"
EOF
    chmod 750 "${TMPDIR}/ghist.${KD}.${filt}.${type}.${STAG}.${_DATE}.sh"
	echo "plotting ${STAG} ghists for ${type} ${KD} ${filt} loop anchors in HSPCs"
	screen -dm -S ${filt}.${type}.${STAG}.${KD} bash -c "bash ${TMPDIR}/ghist.${KD}.${filt}.${type}.${STAG}.${_DATE}.sh"
 done done done
done

#############combined E-P
for KD in ${STAGKDs}; do
for STAG in ${STAGs}; do
for filt in ${FiltsFC1}; do
     cat >"${TMPDIR}/ghist.${KD}.${filt}.${STAG}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd ${diffanchorsn2t}/H3K27overlap/
annotatePeaks.pl ${KD}vsCTRL.loopAnchors.red.${filt}.H3K27accentred.RAD21centred.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${CHIPMERGEBW}/merged.ChIP_CTRL_${STAG}.bedGraph > "${FIGURESDIR}/ghist/sigLoopAnchors/ChIPallEP/STAG/CTRL.allEP_RAD21peaks.${filt}_${KD}.${STAG}.ghist.txt"
annotatePeaks.pl ${KD}vsCTRL.loopAnchors.red.${filt}.H3K27accentred.RAD21centred.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${CHIPMERGEBW}/merged.ChIP_SA2KD_${STAG}.bedGraph > "${FIGURESDIR}/ghist/sigLoopAnchors/ChIPallEP/STAG/SA2KD.allEP_RAD21peaks.${filt}_${KD}.${STAG}.ghist.txt"
annotatePeaks.pl ${KD}vsCTRL.loopAnchors.red.${filt}.H3K27accentred.RAD21centred.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${CHIPMERGEBW}/merged.ChIP_SA1KD_${STAG}.bedGraph > "${FIGURESDIR}/ghist/sigLoopAnchors/ChIPallEP/STAG/SA1KD.allEP_RAD21peaks.${filt}_${KD}.${STAG}.ghist.txt"
EOF
    chmod 750 "${TMPDIR}/ghist.${KD}.${filt}.${STAG}.${_DATE}.sh"
	echo "plotting ${STAG} ghists for all ${KD} ${filt} EP-loop anchors in HSPCs"
	screen -dm -S ${filt}.${KD}.${STAG}.all bash -c "bash ${TMPDIR}/ghist.${KD}.${filt}.${STAG}.${_DATE}.sh"
 done done
done

for STAG in ${STAGs}; do
for KD in ${STAGKDs}; do
for filt in ${FiltsFC1}; do
wc -l ${FIGURESDIR}/ghist/sigLoopAnchors/ChIPallEP/STAG/SA2KD.allEP_RAD21peaks.${filt}_${KD}.${STAG}.ghist.txt
done done done

#############SA1/SA2 plots for KD specific loops by condition
conds="CTRL SA2KD SA1KD"
types2="Promoter Enhancer allEP"
for cond in ${conds};do
for filt in ${FiltsFC1};do
for type in ${types2}; do
plotHIST.sh -g "${cond}.${type}_RAD21peaks.${filt}_SA2KD.SA2.ghist.txt ${cond}.${type}_RAD21peaks.${filt}_SA2KD.SA1.ghist.txt" \
-s "STAG2 STAG1" -c "springgreen2 darkgoldenrod2" -x 1000 -y "0 14" -d ${FIGURESDIR}/ghist/sigLoopAnchors/ChIPallEP/STAG/ -n ${cond}.SA1vsSA2.SA2KD.${type}_RAD21peaks.${filt}

plotHIST.sh -g "${cond}.${type}_RAD21peaks.${filt}_SA1KD.SA2.ghist.txt ${cond}.${type}_RAD21peaks.${filt}_SA1KD.SA1.ghist.txt" \
-s "STAG2 STAG1" -c "springgreen2 darkgoldenrod2" -x 1000 -y "0 14" -d ${FIGURESDIR}/ghist/sigLoopAnchors/ChIPallEP/STAG/ -n ${cond}.SA1vsSA2.SA1KD.${type}_RAD21peaks.${filt}
done done done

cd ${FIGURESDIR}/ghist/sigLoopAnchors/ChIPallEP/STAG/
#mkdir STAG2diffEP
#mkdir STAG1diffEP
mv *.SA2KD.*_RAD21peaks.*hist.pdf STAG2diffEP/
mv *.SA1KD.*_RAD21peaks.*hist.pdf STAG1diffEP/