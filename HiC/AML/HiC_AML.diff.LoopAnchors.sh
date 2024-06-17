#!/bin/bash
#by Alexander Fischer, Jan 2022

#################################################################################
#         Analysis of loop anchors in Cohesin mutated AML samples               #
#################################################################################

#general paths
DIR_PKG="/misc/software/ngs"
DIR_SOFT="/misc/software"
DIR_PKG="${DIR_SOFT}/ngs"
DIR_DATA="/misc/data"
OS=$(lsb_release -c |grep "^Codename" | awk -F: '{print $2}' | sed 's/[[:blank:]]//g')

#setting homer environment
DIR_PKG="/misc/software/ngs"
PATH_PERL=${DIR_SOFT}/package/perl/perl-5.26.1/bin
PATH_SAMTOOLS=${DIR_PKG}/samtools/samtools-1.6/bin
PATH_HOMER=${DIR_PKG}/homer/v4.10/bin
PATH_R=${DIR_SOFT}/package/RBioC/3.4.3/bin
PATH_SKEWER=${DIR_PKG}/skewer/skewer-0.2.2
export PATH=${PATH_R}:${PATH_PERL}:${PATH_SAMTOOLS}:${PATH_HOMER}:${PATH_SKEWER}:${PATH}
export PATH

BEDTOOLS="${DIR_PKG}/bedtools/bedtools2-2.27.1/bin/bedtools"
CHROMSIZES_HG38="${DIR_SOFT}/viewer/IGV/IGVTools_2.3.98/genomes/GRCh38.PRI_p10.chrom.sizes"
R=rbioc_3-12 #R version on server 

#directories and files
TMPDIR="/loctmp"
TAGHICDIR="${DIR_DATA}/processedData/tagDir/3Dchromatin/HiC/GRCh38/AML" #HOMER tagdirs
WORKDIR="${DIR_DATA}/analysis/project_cohesin/Cohesin_AML/HiC"
LOOPDIR="${WORKDIR}/loops"
FIGURESDIR="${WORKDIR}/figures"
DTmatrixdir=${FIGURESDIR}/Heatmaps/DeeptoolsMatrices #output from deeptools comput matrix
MATRIXDIR="${WORKDIR}/Interactionmatrices"

#other data types for integrative analyses
#ChIP
WORKDIR_CHIP="${DIR_DATA}/analysis/project_cohesin/Cohesin_AML/ChIP_analysis"
CHIPPEAKDIR="${WORKDIR_CHIP}/peaks"
PEAKDIRH3K="${WORKDIR_CHIP}/H3K27ac/peaks"
CHIPDIFFDIR="${WORKDIR_CHIP}/diffPeaks"
CHIPMERGEBW="${WORKDIR_CHIP}/mergedBigWigs"
MERGEDBW="${WORKDIR_CHIP}/mergedBigWigs"
MERGEDBW_AML_H3K="${WORKDIR_CHIP}/H3K27ac/mergedBigWigs"
#RNAseqResults
RNAsResults="${DIR_DATA}/analysis/project_cohesin/Cohesin_AML/RNAseq"
#TSS pos from cage seq data
TSSpos="${PROJDIR}/TSS.bed"

#important variables for looping
MUTS="SA2mut RAD21mut"
Filts="up.FC1 down.FC1 up.FC2 down.FC2"
FiltsFC1="up.FC1 down.FC1"

#new dirs for Loop Anchor Anlysis
diffanchorsn2t="${LOOPDIR}/differentialanchorsn2t"
PairedAnch="${diffanchorsn2t}/AnchorsOverlapPairs"
mkdir $diffanchorsn2t
mkdir $PairedAnch


                   #####################################################
                   #           splitting loops into  anchors           #
                   #####################################################


# split all loops into anchors in individual patient groups
###extend interacting pos by 5kb in both directions
declare -a SAMPLES=("HiC_AML_CTRL_combined" "HiC_AML_SA2mut_combined" "HiC_AML_RAD21mut_combined")
cd ${LOOPDIR}
for SAMPLE in "${SAMPLES[@]}";do
awk 'BEGIN {FS=OFS="\t"} { print $1, $2-5000, $3+5000, $1":"$2"-"$6"_plus", "1", "+" ; print $4, $5-5000, $6+5000, $1":"$2"-"$6"_minus", "1", "-"}' <(tail -n +2 ${SAMPLE}.merged.loop.2D.bed) > ${SAMPLE}.merged.loopAnchors.bed
bed2pos.pl ${SAMPLE}.merged.loopAnchors.bed > ${SAMPLE}.merged.loopAnchors.pos.txt
done

for SET in "${SAMPLES[@]}";
do
	cat >"${TMPDIR}/ghist.${SET}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd /loctmp
annotatePeaks.pl ${LOOPDIR}/${SET}.merged.loopAnchors.pos.txt hg38 -size 200000 -hist 250 -bedGraph "${LOOPDIR}/${SET}.2.5K.10K.Insulation.bedGraph" > "${FIGURESDIR}/ghist/${SET}.loopAnchors.Insulation.hist.txt"
annotatePeaks.pl ${LOOPDIR}/${SET}.merged.loopAnchors.pos.txt hg38 -size 200000 -hist 250 -bedGraph "${LOOPDIR}/${SET}.2.5K.10K.DI.bedGraph" > "${FIGURESDIR}/ghist/${SET}.loopAnchors.DI.hist.txt"
EOF
	chmod 750 "${TMPDIR}/ghist.${SET}.${_DATE}.sh"
	echo "plotting ghists for ${SET}"
	screen -dm -S ${SET} bash -c "bash ${TMPDIR}/ghist.${SET}.${_DATE}.sh"
done

# differential/changed loops
# split diff loops into anchors ###extend interacting pos by 5kb in both directions
declare -a DIFFS=("SA2mutvsCTRL" "RAD21mutvsCTRL")

for DIFF in "${DIFFS[@]}";do
cd ${LOOPDIR}
#sort by FC
sort -g -k22,22 ${DIFF}.merged.diff.XYrm.loop.txt | tail -n +2 > ${diffanchorsHOMER}/tmp.loop.txt
sort -g -k9,9 DESEQn2total/${DIFF}.merged.loop.scores.DESeq.norm2Total.loopcoords2.txt | tail -n +2 > ${diffanchorsn2t}/tmp.loop.txt
sort -g -k9,9 DESEQ/${DIFF}.merged.loop.scores.DESeq.loopcoords2.txt | tail -n +2 > ${diffanchorsDS}/tmp.loop.txt
#Loops determined by DESEQ norm2 total
#filter for logFC>0.558 <-0.558 && padj < 0.05, print anchor point start and stop coordinates plus/minus 5000 bases #DSEQ nt2 data
awk 'BEGIN {FS=OFS="\t"} {if( $9>0.558 && $13<.05 ) print $2, $3-5000, $4+5000, $1"_plus", "1", "+" ; if($9>0.558 && $13<.05) print $5, $6-5000, $7+5000, $1"_minus", "1", "-"}' ${diffanchorsn2t}/tmp.loop.txt > ${diffanchorsn2t}/${DIFF}.loopAnchors.up.FC1sig.n2t.bed
awk 'BEGIN {FS=OFS="\t"} {if( $9<-0.558 && $13<.05 ) print $2, $3-5000, $4+5000, $1"_plus", "1", "+" ; if( $9<-0.558 && $13<.05) print $5, $6-5000, $7+5000, $1"_minus", "1", "-"}' ${diffanchorsn2t}/tmp.loop.txt > ${diffanchorsn2t}/${DIFF}.loopAnchors.down.FC1sig.n2t.bed
#filter for logFC>1 <-1&& padj < 0.05, print anchor point start and stop coordinates plus/minus 5000 bases #DSEQ nt2 data
awk 'BEGIN {FS=OFS="\t"} {if( $9>1 && $13<.05 ) print $2, $3-5000, $4+5000, $1"_plus", "1", "+" ; if($9>1 && $13<.05) print $5, $6-5000, $7+5000, $1"_minus", "1", "-"}' ${diffanchorsn2t}/tmp.loop.txt > ${diffanchorsn2t}/${DIFF}.loopAnchors.up.FC2sig.n2t.bed
awk 'BEGIN {FS=OFS="\t"} {if( $9<-1 && $13<.05 ) print $2, $3-5000, $4+5000, $1"_plus", "1", "+" ; if( $9<-1 && $13<.05) print $5, $6-5000, $7+5000, $1"_minus", "1", "-"}' ${diffanchorsn2t}/tmp.loop.txt > ${diffanchorsn2t}/${DIFF}.loopAnchors.down.FC2sig.n2t.bed
#all loops no filter
awk 'BEGIN {FS=OFS="\t"} { print $2, $3-5000, $4+5000, $1"_plus", "1", "+" ; print $5, $6-5000, $7+5000, $1"_minus", "1", "-"}' ${diffanchorsn2t}/tmp.loop.txt > ${diffanchorsn2t}/tmp.bed
head -n -2 ${diffanchorsn2t}/tmp.bed > ${diffanchorsn2t}/tmp.txt > ${diffanchorsn2t}/${DIFF}.loopAnchors.all.bed
####alternative approach: combined anchors file: loop name 
awk 'BEGIN {FS=OFS="\t"} {if( $9>0.558 && $13<.05 ) print $1, $2, $3-5000, $4+5000, $1"_plus", "1", "+", $5, $6-5000, $7+5000, $1"_minus", "1", "-"}' ${diffanchorsn2t}/tmp.loop.txt > ${diffanchorsn2t}/${DIFF}vsCTRL.loopAnchors.paired.up.FC1sig.n2t.txt
awk 'BEGIN {FS=OFS="\t"} {if( $9<-0.558 && $13<.05 ) print $1, $2, $3-5000, $4+5000, $1"_plus", "1", "+", $5, $6-5000, $7+5000, $1"_minus", "1", "-"}' ${diffanchorsn2t}/tmp.loop.txt > ${diffanchorsn2t}/${DIFF}vsCTRL.loopAnchors.paired.down.FC1sig.n2t.txt
done
#check numbers:
cd ${diffanchorsn2t}
for DIFF in "${DIFFS[@]}";do
wc -l ${DIFF}.loopAnchors.up.FC1sig.n2t.bed #
wc -l ${DIFF}.loopAnchors.down.FC1sig.n2t.bed #
wc -l ${DIFF}.loopAnchors.up.FC2sig.n2t.bed #
wc -l ${DIFF}.loopAnchors.down.FC2sig.n2t.bed #
wc -l ${DIFF}.loopAnchors.all.bed  #
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

#check numbers:
cd ${diffanchorsn2t}
for DIFF in "${DIFFS[@]}";do
wc -l ${DIFF}.loopAnchors.red.up.FC1sig.n2t.bed 
wc -l ${DIFF}.loopAnchors.red.down.FC1sig.n2t.bed
wc -l ${DIFF}.loopAnchors.overlap.FC1.n2t.bed
wc -l ${DIFF}.loopAnchors.red.up.FC2sig.n2t.bed 
wc -l ${DIFF}.loopAnchors.red.down.FC2sig.n2t.bed
wc -l ${DIFF}.loopAnchors.overlap.FC2.n2t.bed
done

# convert bed to pos
for DIFF in "${DIFFS[@]}";do
for FC in ${FCs};do
cd ${diffanchorsn2t}
bed2pos.pl ${DIFF}.loopAnchors.up.${FC}sig.n2t.bed > ${DIFF}.loopAnchors.up.${FC}sig.n2t.txt
bed2pos.pl ${DIFF}.loopAnchors.down.${FC}sig.n2t.bed > ${DIFF}.loopAnchors.down.${FC}sig.n2t.txt
bed2pos.pl ${DIFF}.loopAnchors.all.bed >${DIFF}.loopAnchors.all.txt
$BEDTOOLS merge -i <(sort -k1,1 -k2,2n ${DIFF}.loopAnchors.all.bed) -s -d -7500 -c 4,5,6 -o distinct,distinct,distinct > tmp.sorted.all.bed
bed2pos.pl tmp.sorted.all.bed > tmp.${DIFF}.sorted.all.txt
done
done

                     #####################################################
                     #                                                   #
                     #   DI and Insulation at differential anchors       #
                     #                                                   #
                     #####################################################
mkdir -p ${FIGURESDIR}/ghist/Ins_DI/
## annotate Insulation (Ins) and Directionality index (DI) tracks to differential anchor postions
MUTS="SA2mut RAD21mut"
Filts="up.FC1 down.FC1 up.FC2 down.FC2"
_DATE=$(date +%s)
for filt in ${Filts};
do
for SET in ${MUTS};
do
	cat >"${TMPDIR}/ghist.${SET}${filt}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd /loctmp
###INS
annotatePeaks.pl ${diffanchorsn2t}/SA2mutvsCTRL.loopAnchors.${filt}sig.n2t.txt hg38 -size 200000 -hist 250 -bedGraph "${LOOPDIR}/HiC_AML_${SET}_combined.2.5K.10K.Insulation.bedGraph" > "${FIGURESDIR}/ghist/Ins_DI/${SET}.loopAnchors.SA2mut.${filt}.n2t.Insulation.hist.txt"
annotatePeaks.pl ${diffanchorsn2t}/SA2mutvsCTRL.loopAnchors.${filt}sig.n2t.txt hg38 -size 200000 -hist 250 -bedGraph "${LOOPDIR}/HiC_AML_CTRL_combined.2.5K.10K.Insulation.bedGraph" > "${FIGURESDIR}/ghist/Ins_DI/CTRL.loopAnchors.SA2mut.${filt}.n2t.Insulation.hist.txt"
annotatePeaks.pl ${diffanchorsn2t}/RAD21mutvsCTRL.loopAnchors.${filt}sig.n2t.txt hg38 -size 200000 -hist 250 -bedGraph "${LOOPDIR}/HiC_AML_${SET}_combined.2.5K.10K.Insulation.bedGraph" > "${FIGURESDIR}/ghist/Ins_DI/${SET}.loopAnchors.RAD21mut.${filt}.n2t.Insulation.hist.txt"
annotatePeaks.pl ${diffanchorsn2t}/RAD21mutvsCTRL.loopAnchors.${filt}sig.n2t.txt hg38 -size 200000 -hist 250 -bedGraph "${LOOPDIR}/HiC_AML_CTRL_combined.2.5K.10K.Insulation.bedGraph" > "${FIGURESDIR}/ghist/Ins_DI/CTRL.loopAnchors.RAD21mut.${filt}.n2t.Insulation.hist.txt"
###DI
annotatePeaks.pl ${diffanchorsn2t}/SA2mutvsCTRL.loopAnchors.${filt}sig.n2t.txt hg38 -size 200000 -hist 250 -bedGraph "${LOOPDIR}/HiC_AML_${SET}_combined.2.5K.10K.DI.bedGraph" > "${FIGURESDIR}/ghist/Ins_DI/${SET}.loopAnchors.SA2mut.${filt}.n2t.DI.hist.txt"
annotatePeaks.pl ${diffanchorsn2t}/SA2mutvsCTRL.loopAnchors.${filt}sig.n2t.txt hg38 -size 200000 -hist 250 -bedGraph "${LOOPDIR}/HiC_AML_CTRL_combined.2.5K.10K.DI.bedGraph" > "${FIGURESDIR}/ghist/Ins_DI/CTRL.loopAnchors.SA2mut.${filt}.n2t.DI.hist.txt"
annotatePeaks.pl ${diffanchorsn2t}/RAD21mutvsCTRL.loopAnchors.${filt}sig.n2t.txt hg38 -size 200000 -hist 250 -bedGraph "${LOOPDIR}/HiC_AML_${SET}_combined.2.5K.10K.DI.bedGraph" > "${FIGURESDIR}/ghist/Ins_DI/${SET}.loopAnchors.RAD21mut.${filt}.n2t.DI.hist.txt"
annotatePeaks.pl ${diffanchorsn2t}/RAD21mutvsCTRL.loopAnchors.${filt}sig.n2t.txt hg38 -size 200000 -hist 250 -bedGraph "${LOOPDIR}/HiC_AML_CTRL_combined.2.5K.10K.DI.bedGraph" > "${FIGURESDIR}/ghist/Ins_DI/CTRL.loopAnchors.RAD21mut.${filt}.n2t.DI.hist.txt"
EOF
	chmod 750 "${TMPDIR}/ghist.${SET}${filt}.${_DATE}.sh"
	echo "plotting ghists for ${SET}${filt} "
	screen -dm -S ${SET}${filt} bash -c "bash ${TMPDIR}/ghist.${SET}${filt}.${_DATE}.sh"
done
done

#SA2mut: paste histograms of insulation and DI into one table and then into two plots
#plot for both diffloop normalizations
NORMs="n2t"
cd ${FIGURESDIR}/ghist/Ins_DI/
for FC in ${FCs};do
for norm in ${NORMs};do
paste <(cut -f 1,2 CTRL.loopAnchors.SA2mut.up.${FC}.${norm}.Insulation.hist.txt) <(cut -f 2 SA2mut.loopAnchors.SA2mut.up.${FC}.${norm}.Insulation.hist.txt) <(cut -f 2 CTRL.loopAnchors.SA2mut.down.${FC}.${norm}.Insulation.hist.txt) <(cut -f 2 SA2mut.loopAnchors.SA2mut.down.${FC}.${norm}.Insulation.hist.txt) \
<(cut -f 2 CTRL.loopAnchors.SA2mut.up.${FC}.${norm}.DI.hist.txt) <(cut -f 2 SA2mut.loopAnchors.SA2mut.up.${FC}.${norm}.DI.hist.txt) <(cut -f 2 CTRL.loopAnchors.SA2mut.down.${FC}.${norm}.DI.hist.txt) <(cut -f 2 SA2mut.loopAnchors.SA2mut.down.${FC}.${norm}.DI.hist.txt) | tail -n +2 > ${TMPDIR}/tmp.txt
echo $'Dist\tCTRLUP\tSA2mutUP\tCTRLDown\tSA2mutDown\tCTRLUP_DI\tSA2mutUP_DI\tCTRLDown_DI\tSA2mutDown_DI' | cat - ${TMPDIR}/tmp.txt > ${FIGURESDIR}/ghist/Ins_DI/tmp.AML.txt

_DATE=$(date +%s)
cat >"${TMPDIR}/R.plot.P.${_DATE}.R" <<EOF
d <- read.delim("${FIGURESDIR}/ghist/Ins_DI/tmp.AML.txt", header=T)
pdf(file="${FIGURESDIR}/Insulation_DI_plots/CohAML.SA2mut.diff.loopAnchors.up.${FC}.${norm}.DI.I.hist.pdf", height=8, width=8)
plot(d\$Dist, d\$CTRLUP, type="l", col="firebrick", lwd=5, xlab="loop anchor center", ylab="DI & Insulation", ylim=c(-.25,.5))
lines(d\$Dist, d\$SA2mutUP, col="seagreen", lwd=5)
lines(d\$Dist, d\$CTRLUP_DI, col="firebrick", lwd=5, lty=3)
lines(d\$Dist, d\$SA2mutUP_DI, col="seagreen", lwd=5, lty=3)
abline(h = 0, v = 0, col = "gray60")
dev.off()
pdf(file="${FIGURESDIR}/Insulation_DI_plots/CohAML.SA2mut.diff.loopAnchors.down.${FC}.${norm}.DI.I.hist.pdf", height=8, width=8)
plot(d\$Dist, d\$CTRLDown, type="l", col="firebrick", lwd=5, xlab="loop anchor center", ylab="DI & Insulation", ylim=c(-.25,.5))
lines(d\$Dist, d\$SA2mutDown, col="seagreen", lwd=5)
lines(d\$Dist, d\$CTRLDown_DI, col="firebrick", lwd=5, lty=3)
lines(d\$Dist, d\$SA2mutDown_DI, col="seagreen", lwd=5, lty=3)
abline(h = 0, v = 0, col = "gray60")
dev.off()
EOF
chmod 750 "${TMPDIR}/R.plot.P.${_DATE}.R"
R < ${TMPDIR}/R.plot.P.${_DATE}.R  --no-save
rm ${TMPDIR}/R.plot.P.${_DATE}.R
done
done



#RAD21mut: paste histograms of insulation and DI into one table and then into two plots
cd ${FIGURESDIR}/ghist/Ins_DI/
for FC in ${FCs};do
for norm in ${NORMs};do
paste <(cut -f 1,2 CTRL.loopAnchors.RAD21mut.up.${FC}.${norm}.Insulation.hist.txt) <(cut -f 2 RAD21mut.loopAnchors.RAD21mut.up.${FC}.${norm}.Insulation.hist.txt) <(cut -f 2 CTRL.loopAnchors.RAD21mut.down.${FC}.${norm}.Insulation.hist.txt) <(cut -f 2 RAD21mut.loopAnchors.RAD21mut.down.${FC}.${norm}.Insulation.hist.txt) \
<(cut -f 2 CTRL.loopAnchors.RAD21mut.up.${FC}.${norm}.DI.hist.txt) <(cut -f 2 RAD21mut.loopAnchors.RAD21mut.up.${FC}.${norm}.DI.hist.txt) <(cut -f 2 CTRL.loopAnchors.RAD21mut.down.${FC}.${norm}.DI.hist.txt) <(cut -f 2 RAD21mut.loopAnchors.RAD21mut.down.${FC}.${norm}.DI.hist.txt) | tail -n +2 > ${TMPDIR}/tmp.txt
echo $'Dist\tCTRLUP\tRAD21mutUP\tCTRLDown\tRAD21mutDown\tCTRLUP_DI\tRAD21mutUP_DI\tCTRLDown_DI\tRAD21mutDown_DI' | cat - ${TMPDIR}/tmp.txt > ${FIGURESDIR}/ghist/Ins_DI/tmp.AML.txt

_DATE=$(date +%s)
cat >"${TMPDIR}/R.plot.P.${_DATE}.R" <<EOF
d <- read.delim("${FIGURESDIR}/ghist/Ins_DI/tmp.AML.txt", header=T)
pdf(file="${FIGURESDIR}/Insulation_DI_plots/CohAML.RAD21mut.diff.loopAnchors.up.${FC}.${norm}.DI.I.hist.pdf", height=8, width=8)
plot(d\$Dist, d\$CTRLUP, type="l", col="firebrick", lwd=5, xlab="loop anchor center", ylab="DI & Insulation", ylim=c(-.25,.5))
lines(d\$Dist, d\$RAD21mutUP, col="orchid4", lwd=5)
lines(d\$Dist, d\$CTRLUP_DI, col="firebrick", lwd=5, lty=3)
lines(d\$Dist, d\$RAD21mutUP_DI, col="orchid4", lwd=5, lty=3)
abline(h = 0, v = 0, col = "gray60")
dev.off()
pdf(file="${FIGURESDIR}/Insulation_DI_plots/CohAML.RAD21mut.diff.loopAnchors.down.${FC}.${norm}.DI.I.hist.pdf", height=8, width=8)
plot(d\$Dist, d\$CTRLDown, type="l", col="firebrick", lwd=5, xlab="loop anchor center", ylab="DI & Insulation", ylim=c(-.25,.5))
lines(d\$Dist, d\$RAD21mutDown, col="orchid4", lwd=5)
lines(d\$Dist, d\$CTRLDown_DI, col="firebrick", lwd=5, lty=3)
lines(d\$Dist, d\$RAD21mutDown_DI, col="orchid4", lwd=5, lty=3)
abline(h = 0, v = 0, col = "gray60")
dev.off()
EOF
chmod 750 "${TMPDIR}/R.plot.P.${_DATE}.R"
R < ${TMPDIR}/R.plot.P.${_DATE}.R  --no-save
rm ${TMPDIR}/R.plot.P.${_DATE}.R
done
done


                          ################################################
                          #         loop anchor RAD21/CTCF overlap       #
                          ################################################

CHIPs="RAD21 CTCF"

for MUT in ${MUTS}; do
for filt in ${Filts}; do
#overlap with all CTCF peaks
cd ${diffanchorsn2t}
$BEDTOOLS intersect -b ${CHIPPEAKDIR}/Allpat_mergePeaks_CTCF.filtered.peaks.bed -a ${MUT}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed -u > AnchorChIPoverlap/${MUT}vsCTRL.loopAnchors.red.${filt}.n2t.CTCFover.bed
$BEDTOOLS intersect -b ${CHIPPEAKDIR}/Allpat_mergePeaks_CTCF.filtered.peaks.bed -a ${MUT}vsCTRL.loopAnchors.all.bed -u > AnchorChIPoverlap/${MUT}vsCTRL.loopAnchors.all.n2t.CTCFover.bed
#overlap with all RAD21 peaks
cd ${diffanchorsn2t}
$BEDTOOLS intersect -b ${CHIPPEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.bed -a ${MUT}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed -u > AnchorChIPoverlap/${MUT}vsCTRL.loopAnchors.red.${filt}.n2t.RAD21over.bed
$BEDTOOLS intersect -b ${CHIPPEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.bed -a ${MUT}vsCTRL.loopAnchors.all.bed -u > AnchorChIPoverlap/${MUT}vsCTRL.loopAnchors.all.n2t.RAD21over.bed
done
done

#summarize loop anchor overlaps found (DESEQ n2t)
cd ${diffanchorsn2t}
for MUT in ${MUTS}; do
for ChIP in ${CHIPs}; do
for filt in ${Filts}; do
echo "# ${MUT} all loopanchors ${ChIP} overlap#" $(wc -l < AnchorChIPoverlap/${MUT}vsCTRL.loopAnchors.all.n2t.${ChIP}over.bed)
echo "# ${MUT} loopanchors ${filt} ${ChIP} overlap#" $(wc -l < AnchorChIPoverlap/${MUT}vsCTRL.loopAnchors.red.${filt}.n2t.${ChIP}over.bed)
done
done
done


                     #####################################################
                     #   ChIP peaks  at differential anchors             #
                     #####################################################

#repeat intersection with a and b exchanged


#MUT vs CTRL sets: differential or complete (all) anchor sets
for MUT in ${MUTS}; do
for filt in ${Filts}; do
#centered to all CTCF peaks
cd ${diffanchorsn2t}
$BEDTOOLS intersect -a ${CHIPPEAKDIR}/Allpat_mergePeaks_CTCF.filtered.peaks.bed -b ${MUT}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed -u > AnchorChIPoverlap/${MUT}vsCTRL.loopAnchors.red.${filt}.n2t.CTCFcentered.bed
$BEDTOOLS intersect -a ${CHIPPEAKDIR}/Allpat_mergePeaks_CTCF.filtered.peaks.bed -b ${MUT}vsCTRL.loopAnchors.all.bed -u > AnchorChIPoverlap/${MUT}vsCTRL.loopAnchors.all.n2t.CTCFcentered.bed
#centered to RAD21 peaks
$BEDTOOLS intersect -a ${CHIPPEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.bed -b ${MUT}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed -u > AnchorChIPoverlap/${MUT}vsCTRL.loopAnchors.red.${filt}.n2t.RAD21centered.bed
$BEDTOOLS intersect -a ${CHIPPEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.bed -b ${MUT}vsCTRL.loopAnchors.all.bed -u > AnchorChIPoverlap/${MUT}vsCTRL.loopAnchors.all.n2t.RAD21centered.bed
$BEDTOOLS intersect -a ${CHIPPEAKDIR}/Allpat_mergePeaks_RAD21_stringent.filtered.peaks.bed -b ${MUT}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed -u > AnchorChIPoverlap/${MUT}vsCTRL.loopAnchors.red.${filt}.n2t.RAD21.stringent.centered.bed
$BEDTOOLS intersect -a ${CHIPPEAKDIR}/Allpat_mergePeaks_RAD21_stringent.filtered.peaks.bed -b ${MUT}vsCTRL.loopAnchors.all.bed -u > AnchorChIPoverlap/${MUT}vsCTRL.loopAnchors.all.n2t.RAD21.stringent.centered.bed
done
done
#summarize loop anchor overlaps found (n2t)
cd ${diffanchorsn2t}
for MUT in ${MUTS}; do
for ChIP in ${CHIPs}; do
for filt in ${FiltsFC1}; do
echo "# ${MUT} all loopanchors ${ChIP} overlap#" $(wc -l < AnchorChIPoverlap/${MUT}vsCTRL.loopAnchors.all.n2t.${ChIP}centered.bed)
echo "# ${MUT} loopanchors ${filt} ${ChIP} overlap#" $(wc -l < AnchorChIPoverlap/${MUT}vsCTRL.loopAnchors.red.${filt}.n2t.${ChIP}centered.bed)
done
done
done

#convert to pos
for MUT in ${MUTS}; do
for ChIP in ${ChIPs}; do
for filt in ${Filts}; do
cd ${diffanchorsn2t}
bed2pos.pl AnchorChIPoverlap/${MUT}vsCTRL.loopAnchors.red.${filt}.n2t.${ChIP}centered.bed > AnchorChIPoverlap/${MUT}vsCTRL.loopAnchors.red.${filt}.n2t.${ChIP}centered.pos.txt
bed2pos.pl AnchorChIPoverlap/${MUT}vsCTRL.loopAnchors.all.n2t.${ChIP}centered.bed > AnchorChIPoverlap/${MUT}vsCTRL.loopAnchors.all.n2t.${ChIP}centered.pos.txt
bed2pos.pl AnchorChIPoverlap/${MUT}vsCTRL.loopAnchors.red.${filt}.n2t.RAD21.stringent.centered.bed > AnchorChIPoverlap/${MUT}vsCTRL.loopAnchors.red.${filt}.n2t.RAD21.stringent.centered.pos.txt
done
done
done


##check CTCF and RAD21 co-associated loop anchor overlap
$BEDTOOLS intersect -a CoAML_RAD21.tmpPeaks.bed -b CoAML_CTCF.tmpPeaks.bed  -u > tmp.RAD21_CTCF_PEAKoverlap.bed #56749
$BEDTOOLS intersect -a CoAML_CTCF.tmpPeaks.bed  -b CoAML_RAD21.tmpPeaks.bed -v > tmp.CTCF_wo_RAD21_PEAKoverlap.bed #6282
$BEDTOOLS intersect -b CoAML_CTCF.tmpPeaks.bed  -a CoAML_RAD21.tmpPeaks.bed -v > tmp.RAD21_woCTCF_PEAKoverlap.bed #124779 

for SET in ${MUTS}; do
for filt in ${Filts}; do
$BEDTOOLS intersect -b tmp.RAD21_CTCF_PEAKoverlap.bed -a ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed -u > AnchorChIPoverlap/${SET}vsCTRL.loopAnchors.red.${filt}.RAD21andCTCF.overlap.bed
$BEDTOOLS intersect -a tmp.RAD21_CTCF_PEAKoverlap.bed -b ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed -u > AnchorChIPoverlap/${SET}vsCTRL.loopAnchors.red.${filt}.RAD21andCTCF.centred.bed
$BEDTOOLS intersect -b tmp.CTCF_wo_RAD21_PEAKoverlap.bed -a ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed  -u > AnchorChIPoverlap/${SET}vsCTRL.loopAnchors.red.${filt}.CTCFwoRAD21.overlap.bed
done
done

for SET in ${MUTS}; do
for filt in ${FiltsFC1}; do
echo "# ${SET} loop anchors ${filt} with cohesin #" $(wc -l < AnchorChIPoverlap/${SET}vsCTRL.loopAnchors.red.${filt}.n2t.RAD21over.bed)
echo "# ${SET} loop anchors ${filt} with cohesin and CTCF #" $(wc -l < AnchorChIPoverlap/${SET}vsCTRL.loopAnchors.red.${filt}.RAD21andCTCF.overlap.bed)
echo "# ${SET} loop anchors ${filt} with CTCF and no Cohesin #" $(wc -l < AnchorChIPoverlap/${SET}vsCTRL.loopAnchors.red.${filt}.CTCFwoRAD21.overlap.bed)
echo "# ${SET} cohesin with CTCF at loop anchors ${filt} #" $(wc -l < AnchorChIPoverlap/${SET}vsCTRL.loopAnchors.red.${filt}.RAD21andCTCF.centred.bed)
done
done
                     #####################################################
                     #   H3K27ac  at differential Loopanchors            #
                     #####################################################
mkdir ${diffanchorsn2t}/H3K27overlap
#check which H3K27ac peaks are at diffloop anchors
##all enhancers/promoters used for the differential H3K27ac analysis
pos2bed.pl ${PEAKDIRH3K}/Allpat_mergePeaks_H3K27ac.CNVnorm.filtered.XYrem.peaks.txt > ${TMPDIR}/tmp.peaks.bed

##overlap with sig. loop anchor regions n2t results ##look from both sides (a/b options)
cd ${diffanchorsn2t}
for SET in ${MUTS}; do
for filt in ${Filts}; do
$BEDTOOLS intersect -a ${TMPDIR}/tmp.peaks.bed -b ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed  -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27accentred.bed
$BEDTOOLS intersect -b ${TMPDIR}/tmp.peaks.bed -a ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed  -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27acOverlap.bed
done

for SET in ${MUTS}; do
for filt in ${Filts}; do
echo "# ${SET} loop anchors ${filt} with H3K27ac peaks #" $(wc -l < H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27acOverlap.bed)
echo "# ${SET} H3K27ac peaks at loop anchors ${filt} #" $(wc -l < H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27accentred.bed)
echo "# ${SET} loop anchors ${filt} EP contacts #" $(wc -l < H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27acOverlap.bed)/$(wc -l < ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed)
done
done

#cohesin associated enhancers at loop anchors
cd ${diffanchorsn2t}
pos2bed.pl ${PEAKDIRH3K}/Allpat_mergePeaks_H3K27ac.CNVnorm.filtered.XYrem.peaks.txt   > CoAML_H3K27ac.tmp.peaks.bed
pos2bed.pl ${CHIPPEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.txt > CoAML_RAD21.tmpPeaks.bed 
pos2bed.pl ${CHIPPEAKDIR}/Allpat_mergePeaks_CTCF.filtered.peaks.txt > CoAML_CTCF.tmpPeaks.bed 
$BEDTOOLS intersect -a CoAML_H3K27ac.tmp.peaks.bed -b CoAML_RAD21.tmpPeaks.bed  -u > H3K27overlap/CohesinAssEnhancers.bed #57043
$BEDTOOLS intersect -b CoAML_H3K27ac.tmp.peaks.bed -a CoAML_RAD21.tmpPeaks.bed  -u > H3K27overlap/CohesinAssEnhancers.RAD21centred.bed #72087


###diffloop anchor overlaps with cohesin ass. enhancers
for SET in ${MUTS}; do
for filt in ${Filts}; do
$BEDTOOLS intersect -b H3K27overlap/CohesinAssEnhancers.bed -a ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed  -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.Overlap.bed #Anchors with Cohesin Ass. enhancer
$BEDTOOLS intersect -b H3K27overlap/CohesinAssEnhancers.bed -a ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed  -v > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.Overlap.bed #Anchors without Cohesin Ass. enhancer
$BEDTOOLS intersect -a H3K27overlap/CohesinAssEnhancers.bed -b ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed  -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.centred.bed  #Anchors with Cohesin Ass. enhancer centred on H3K27ac
$BEDTOOLS intersect -a H3K27overlap/CohesinAssEnhancers.RAD21centred.bed -b ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed  -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.RAD21centred.bed #Anchors with Cohesin Ass. enhancer centred on RAD21
done
done

###report both loop and peak ID
cd ${diffanchorsn2t}
mkdir H3K27overlap/2wayOverlap
for SET in ${MUTS}; do
for filt in ${FiltsFC1}; do
$BEDTOOLS intersect -b CoAML_H3K27ac.tmp.peaks.bed -a ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed  -wao > H3K27overlap/2wayOverlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27ac.2wayOverlap.bed
$BEDTOOLS intersect -b H3K27overlap/CohesinAssEnhancers.bed -a ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed  -wao > H3K27overlap/2wayOverlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.2wayOverlap.bed #Anchors with Cohesin Ass. enhancer
$BEDTOOLS intersect -b H3K27overlap/CohesinAssEnhancers.RAD21centred.bed -a ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed  -wao > H3K27overlap/2wayOverlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.RAD21centred.2wayOverlap.bed #Anchors with Cohesin Ass. enhancer centred on RAD21
done
done

##same for complete loop anchors sets
cd ${diffanchorsn2t}
for SET in ${MUTS}; do
$BEDTOOLS intersect -b H3K27overlap/CohesinAssEnhancers.bed -a ${SET}vsCTRL.loopAnchors.all.bed -u > H3K27overlap/${SET}vsCTRL.loopAnchors.all.CohesinAssEnhancers.centred.bed #14527
$BEDTOOLS intersect -b CoAML_H3K27ac.tmp.peaks.bed -a ${SET}vsCTRL.loopAnchors.all.bed  -wao > H3K27overlap/2wayOverlap/${SET}vsCTRL.loopAnchors.all.H3K27ac.2wayOverlap.bed
$BEDTOOLS intersect -b CoAML_RAD21.tmpPeaks.bed -a ${SET}vsCTRL.loopAnchors.all.bed  -wao > AnchorChIPoverlap/${SET}vsCTRL.loopAnchors.all.RAD21.2wayOverlap.bed
done

### diffloop anchor overlaps without cohesin ass. enhancers i.e. Structural loops
for SET in ${MUTS}; do
for filt in ${Filts}; do
##without cohesin ass. enhancers, but cohesin
$BEDTOOLS intersect -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.Overlap.bed -b CoAML_RAD21.tmpPeaks.bed -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.RAD21.Overlap.bed
##without cohesin ass. enhancers, without cohesin
$BEDTOOLS intersect -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.Overlap.bed -b CoAML_RAD21.tmpPeaks.bed -v > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.wo.RAD21.Overlap.bed
##without cohesin ass. enhancers, but cohesin with CTCF
$BEDTOOLS intersect -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.RAD21.Overlap.bed -b CoAML_CTCF.tmpPeaks.bed -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.RAD21.CTCF.Overlap.bed
##without cohesin ass. enhancers, but cohesin without CTCF
$BEDTOOLS intersect -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.RAD21.Overlap.bed -b CoAML_CTCF.tmpPeaks.bed -v > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.RAD21.woCTCF.Overlap.bed
##without cohesin ass. enhancers, without cohesin
$BEDTOOLS intersect -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.wo.RAD21.Overlap.bed -b CoAML_CTCF.tmpPeaks.bed -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.wo.RAD21.CTCF.Overlap.bed
$BEDTOOLS intersect -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.wo.RAD21.Overlap.bed -b CoAML_CTCF.tmpPeaks.bed -v > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.wo.RAD21.woCTCF.Overlap.bed
$BEDTOOLS intersect -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.wo.RAD21.woCTCF.Overlap.bed -b CoAML_H3K27ac.tmp.peaks.bed -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.wo.RAD21.woCTCF.H3K27ac.Overlap.bed
$BEDTOOLS intersect -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.wo.RAD21.woCTCF.Overlap.bed -b CoAML_H3K27ac.tmp.peaks.bed -v > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.wo.RAD21.woCTCF.woH3K27ac.Overlap.bed
done
done

#### structural anchors: 2 way overlaps with RAD21 (stringent) and CTCF
STRUCLOOPDIR=${diffanchorsn2t}/AnchorChIPoverlap/structural_Loops
mkdir ${STRUCLOOPDIR}
cp ${diffanchorsn2t}/H3K27overlap/*.wo.CohesinAssEnhancers.Overlap.bed ${STRUCLOOPDIR}/
cd ${diffanchorsn2t}
for SET in ${MUTS}; do
for filt in ${Filts}; do
$BEDTOOLS intersect -a ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.Overlap.bed -b CoAML_RAD21.tmpPeaks.bed -wao > ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.structural.RAD21.2wayoverlap.bed
$BEDTOOLS intersect -a ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.Overlap.bed -b ${CHIPPEAKDIR}/Allpat_mergePeaks_RAD21_stringent.filtered.peaks.bed -wao > ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.structural.strRAD21.2wayoverlap.bed
$BEDTOOLS intersect -a ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.Overlap.bed -b CoAML_CTCF.tmpPeaks.bed -wao > ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.structural.CTCF.2wayoverlap.bed
$BEDTOOLS intersect -a ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.Overlap.bed -b CoAML_H3K27ac.tmp.peaks.bed -wao > ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.structural.H3K27ac.2wayoverlap.bed
done
done

#center structual anchors on RAD21
cd ${diffanchorsn2t}
for SET in ${MUTS}; do
for filt in ${Filts}; do
$BEDTOOLS intersect -b ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.Overlap.bed -a CoAML_RAD21.tmpPeaks.bed -u > ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.structural.RAD21.centred.bed
$BEDTOOLS intersect -a ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.Overlap.bed -b ${CHIPPEAKDIR}/Allpat_mergePeaks_RAD21_stringent.filtered.peaks.bed -u > ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.structural.strRAD21.centred.bed
done
done


#cohesin associated enhancers association with CTCF
for SET in ${MUTS}; do
for filt in ${Filts}; do
#cohesin ass. loop anchors with or without CTCF
$BEDTOOLS intersect -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.Overlap.bed -b CoAML_CTCF.tmpPeaks.bed -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.CTCF.Overlap.bed
$BEDTOOLS intersect -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.Overlap.bed -b CoAML_CTCF.tmpPeaks.bed -v > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.wo.CTCF.Overlap.bed
#coh enhancers with CTCF
$BEDTOOLS intersect -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.centred.bed -b CoAML_CTCF.tmpPeaks.bed   -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancersCTCF.centred.bed
#coh enhancers without CTCF
$BEDTOOLS intersect -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.centred.bed -b CoAML_CTCF.tmpPeaks.bed   -v > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.woCTCF.centred.bed
#loops with coh enhancers with CTCF in coh enhancer# Problem: loops can have both types of coh enhancers!
$BEDTOOLS intersect -b H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancersCTCF.centred.bed -a ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed   -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CTCF-CohesinAssEnhancers.Overlap.bed
#loops with coh enhancers without CTCF in coh enhancer # Problem: loops can have both types of coh enhancers!
$BEDTOOLS intersect -b H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.woCTCF.centred.bed -a ${SET}vsCTRL.loopAnchors.red.${filt}sig.n2t.bed   -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.non-CTCF-CohesinAssEnhancers.Overlap.bed
done
done

#look at values #these are plotted in the HiC_AML.diff.LoopAnchors.sh script
FiltsFC1="up.FC1 down.FC1"
for SET in ${MUTS}; do
for filt in ${FiltsFC1}; do
echo "# ${SET} loop anchors ${filt} with cohesin-H3K27ac #" $(wc -l < H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.Overlap.bed)
echo "# ${SET} loop anchors ${filt} without cohesin-H3K27ac #" $(wc -l < H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.Overlap.bed)
echo "# ${SET} loop anchors ${filt} with cohesin-H3K27ac and CTCF#" $(wc -l < H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.CTCF.Overlap.bed)
echo "# ${SET} loop anchors ${filt} with cohesin-H3K27ac without CTCF#" $(wc -l < H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.wo.CTCF.Overlap.bed)
done
done





# enter Anchors on Cohesin associated enhancer (could be used for clustering of regions etc)
## H3K27 centred anchors centred on cohesin --> "Cohesin-dependent-enhancers at loop anchors" #use for EP anchor regions to display/cluster and for RAD21-PSEA
###also centre on H3K27ac regions for H3K27ac-PSEA
cd ${diffanchorsn2t}
for SET in ${MUTS}; do
for filt in ${Filts}; do
$BEDTOOLS intersect -a CoAML_RAD21.tmpPeaks.bed -b H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27accentred.bed -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27accentred.RAD21centred.bed
$BEDTOOLS intersect -b CoAML_RAD21.tmpPeaks.bed -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27accentred.bed -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27accentred.RAD21overlap.bed

$BEDTOOLS intersect -a ${CHIPPEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.bed -b H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27accentred.bed -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27accentred.RAD21revcentred.bed
$BEDTOOLS intersect -b ${CHIPPEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.bed -a H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27accentred.bed -u > H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27accentred.RAD21revoverlap.bed
done
done

for SET in ${MUTS}; do
for filt in ${Filts}; do
echo "# ${SET} loop anchors ${filt} with H3K27ac peaks #" $(wc -l < H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27acOverlap.bed)
echo "# ${SET} RAD21 peaks in H3K27ac marks at loop anchors ${filt} #" $(wc -l < H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27accentred.RAD21centred.bed)
echo "# ${SET} H3K27ac marks at loop anchors ${filt} with coehsin#" $(wc -l < H3K27overlap/${SET}vsCTRL.loopAnchors.red.${filt}.H3K27accentred.RAD21overlap.bed)
done
done


                     #############################################################################
                     #   Transcription Start sites (TSS)  at differential Loopanchors            #
                     #############################################################################

##overlap the potential E-P loop anchors with TSS postitions (CAGE seq derived)
cd ${diffanchorsn2t}/H3K27overlap
mkdir TSSanchorOverlap

for SET in ${MUTS}; do
for filt in ${Filts}; do
#report all EP loop: this can be used for enhancer promoter distinction:
$BEDTOOLS intersect -a ${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.Overlap.bed -b ${TSSpos} -wao > TSSanchorOverlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.Overlap.TSS.2wayoverlap.bed
#report all structural loops:
$BEDTOOLS intersect -a ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.Overlap.bed -b ${TSSpos} -wao > ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.structural.TSS.2wayoverlap.bed
#report only TSS:
$BEDTOOLS intersect -b ${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.Overlap.bed -a ${TSSpos} -u > TSSanchorOverlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.Overlap.TSS.bed
$BEDTOOLS intersect -b ${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.wo.CTCF.overlap.bed -a ${TSSpos} -u > TSSanchorOverlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.wo.CTCF.Overlap.TSS.bed
$BEDTOOLS intersect -b ${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancersCTCF.overlap.bed -a ${TSSpos} -u > TSSanchorOverlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.CTCF.Overlap.TSS.bed
$BEDTOOLS intersect -b ${STRUCLOOPDIR}/${SET}vsCTRL.loopAnchors.red.${filt}.wo.CohesinAssEnhancers.Overlap.bed -a ${TSSpos} -u > TSSanchorOverlap/${SET}vsCTRL.loopAnchors.red.${filt}.structural.TSS.bed
done
done

for SET in ${MUTS}; do
for filt in ${Filts}; do
wc -l TSSanchorOverlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.Overlap.TSS.bed
wc -l TSSanchorOverlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.wo.CTCF.Overlap.TSS.bed
wc -l TSSanchorOverlap/${SET}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.CTCF.Overlap.TSS.bed
wc -l TSSanchorOverlap/${SET}vsCTRL.loopAnchors.red.${filt}.structural.TSS.bed
done
done


                     ########################################################################
                     #                                                                      #
                     #                paired differential loop anchors                      #
                     #                                                                      #
                     ########################################################################

#define which loop anchors are pairs (i.e. belong to the same loop)
#anchor pairs will be divided by H3K27ac presence on one or both sides
#run the HIC_AML.Diff_LoopAnchor.stats.rmd script:
cd ${WORKDIR}/scripts
rbioc_3-12 -e "rmarkdown::render('HIC_AML.Diff_LoopAnchor.stats.rmd',output_file='HIC_AML.Diff_LoopAnchor.stats.html')" #generates html report and output files

                          ##############################################################
                          #             TSS overlaps for anchor pairs                  #
                          ##############################################################


#TSS overlaps separately by paired / unpaired status 
mkdir ${PairedAnch}/TSSoverlap
###which TSS are in loop anchor categories?
for MUT in ${MUTS}; do
for filt in ${FiltsFC1}; do
#paired H3K27ac anchors  (both sides have the H3K27ac)
$BEDTOOLS intersect -b ${PairedAnch}/H3K27acOverlap_${filt}_anchors_${MUT}_paired_plus.bed -a ${diffanchorsn2t}/H3K27overlap/${TSSpos} -u > ${PairedAnch}/TSSoverlap/H3K27acOverlap_${filt}_anchors_${MUT}_paired_plus.TSS.bed
$BEDTOOLS intersect -b ${PairedAnch}/H3K27acOverlap_${filt}_anchors_${MUT}_paired_minus.bed -a ${diffanchorsn2t}/H3K27overlap/${TSSpos} -u > ${PairedAnch}/TSSoverlap/H3K27acOverlap_${filt}_anchors_${MUT}_paired_minus.TSS.bed
#"unpaired" H3K27ac anchors (only one side has the H3K27ac)
$BEDTOOLS intersect -b ${PairedAnch}/H3K27acOverlap_${filt}_anchors_${MUT}_unpairedPLUS_plus.bed -a ${diffanchorsn2t}/H3K27overlap/${TSSpos} -u > ${PairedAnch}/TSSoverlap/H3K27acOverlap_${filt}_anchors_${MUT}_unpairedPLUS_plus.TSS.bed
$BEDTOOLS intersect -b ${PairedAnch}/H3K27acOverlap_${filt}_anchors_${MUT}_unpairedMINUS_minus.bed -a ${diffanchorsn2t}/H3K27overlap/${TSSpos} -u > ${PairedAnch}/TSSoverlap/H3K27acOverlap_${filt}_anchors_${MUT}_unpairedMINUS_minus.TSS.bed
$BEDTOOLS intersect -b ${PairedAnch}/H3K27acOverlap_${filt}_anchors_${MUT}_unpairedPLUS_minus.bed -a ${diffanchorsn2t}/H3K27overlap/${TSSpos} -u > ${PairedAnch}/TSSoverlap/H3K27acOverlap_${filt}_anchors_${MUT}_unpairedPLUS_minus.TSS.bed
$BEDTOOLS intersect -b ${PairedAnch}/H3K27acOverlap_${filt}_anchors_${MUT}_unpairedMINUS_plus.bed -a ${diffanchorsn2t}/H3K27overlap/${TSSpos} -u > ${PairedAnch}/TSSoverlap/H3K27acOverlap_${filt}_anchors_${MUT}_unpairedMINUS_plus.TSS.bed
done
done

###which TSS are in which loop anchor? -->two sided overlap; could be used for sorting anhcors by promoter/Enhancer 
###important: use wao to alos report anchors without TSS!
for MUT in ${MUTS}; do
for filt in ${FiltsFC1}; do
for strand in ${strands}; do
$BEDTOOLS intersect -a ${PairedAnch}/H3K27acOverlap_${filt}_anchors_${MUT}_paired_${strand}.bed -b ${diffanchorsn2t}/H3K27overlap/${TSSpos} -wao > ${PairedAnch}/TSSoverlap/H3K27acOverlap_${filt}_anchors_${MUT}_paired_${strand}.TSS.2wayoverlap.bed
$BEDTOOLS intersect -a ${PairedAnch}/H3K27acOverlap_${filt}_anchors_${MUT}_unpairedPLUS_${strand}.bed -b ${diffanchorsn2t}/H3K27overlap/${TSSpos} -wao > ${PairedAnch}/TSSoverlap/H3K27acOverlap_${filt}_anchors_${MUT}_unpairedPLUS_${strand}.TSS.2wayoverlap.bed
$BEDTOOLS intersect -a ${PairedAnch}/H3K27acOverlap_${filt}_anchors_${MUT}_unpairedMINUS_${strand}.bed -b ${diffanchorsn2t}/H3K27overlap/${TSSpos} -wao > ${PairedAnch}/TSSoverlap/H3K27acOverlap_${filt}_anchors_${MUT}_unpairedMINUS_${strand}.TSS.2wayoverlap.bed
done done done


#define loop anchors as enhancer or promoter based on these overlaps
##run the HIC_AML.diffEP_LoopAnchor_TSS_filtering.rmd rmd script
cd ${WORKDIR}/scripts
rbioc_3-12 -e "rmarkdown::render('HIC_AML.diffEP_LoopAnchor_TSS_filtering.rmd',output_file='HIC_AML.diffEP_LoopAnchor_TSS_filtering.html')" #generates html report and output files


                          ##############################################################
                          #       RAD21 centering of TSS-filtered anchor pairs         #
                          ##############################################################
mkdir ${PairedAnch}/TSSoverlap/RAD21centred
###The files are now split by EP status in R: "true" EP, EE, PP etc
####intersect them again with RAD21 --> useful for centering tracks and PSEA
for SET in ${MUTS}; do
for filt in ${FiltsFC1}; do
$BEDTOOLS intersect -b ${diffanchorsn2t}/H3K27overlap/CohesinAssEnhancers.RAD21centred.bed -a ${PairedAnch}/TSSoverlap/${SET}_${filt}paired_EP_anchors_${SET}_enhancers.bed -wao > ${PairedAnch}/TSSoverlap/RAD21centred/${SET}_${filt}_paired_EP_anchors_enhancers.RAD21.2wayoverlap.bed
$BEDTOOLS intersect -b ${diffanchorsn2t}/H3K27overlap/CohesinAssEnhancers.RAD21centred.bed -a ${PairedAnch}/TSSoverlap/${SET}_${filt}paired_EP_anchors_${SET}_promoters.bed -wao > ${PairedAnch}/TSSoverlap/RAD21centred/${SET}_${filt}_paired_EP_anchors_promoters.RAD21.2wayoverlap.bed
$BEDTOOLS intersect -b ${diffanchorsn2t}/H3K27overlap/CohesinAssEnhancers.RAD21centred.bed -a ${PairedAnch}/TSSoverlap/${SET}_${filt}paired_PP_anchors_${SET}_promoters.plus.bed -wao > ${PairedAnch}/TSSoverlap/RAD21centred/${SET}_${filt}_paired_PP_anchors_promoters.plus.RAD21.2wayoverlap.bed
$BEDTOOLS intersect -b ${diffanchorsn2t}/H3K27overlap/CohesinAssEnhancers.RAD21centred.bed -a ${PairedAnch}/TSSoverlap/${SET}_${filt}paired_PP_anchors_${SET}_promoters.minus.bed -wao > ${PairedAnch}/TSSoverlap/RAD21centred/${SET}_${filt}_paired_PP_anchors_promoters.minus.RAD21.2wayoverlap.bed
$BEDTOOLS intersect -b ${diffanchorsn2t}/H3K27overlap/CohesinAssEnhancers.RAD21centred.bed -a ${PairedAnch}/TSSoverlap/${SET}_${filt}paired_EE_anchors_${SET}_enhancers.plus.bed -wao > ${PairedAnch}/TSSoverlap/RAD21centred/${SET}_${filt}_paired_EE_anchors_enhancers.plus.RAD21.2wayoverlap.bed
$BEDTOOLS intersect -b ${diffanchorsn2t}/H3K27overlap/CohesinAssEnhancers.RAD21centred.bed -a ${PairedAnch}/TSSoverlap/${SET}_${filt}paired_EE_anchors_${SET}_enhancers.minus.bed -wao > ${PairedAnch}/TSSoverlap/RAD21centred/${SET}_${filt}_paired_EE_anchors_enhancers.minus.RAD21.2wayoverlap.bed
done
done
#same for unpaired H3K27ac anchors
mkdir ${PairedAnch}/TSSoverlap/RAD21centred/reg_struc_pairs/
for SET in ${MUTS}; do
for filt in ${FiltsFC1}; do
$BEDTOOLS intersect -b ${diffanchorsn2t}/H3K27overlap/CohesinAssEnhancers.RAD21centred.bed -a ${PairedAnch}/TSSoverlap/reg_struc_pairs/${SET}_${filt}paired_PsP_anchors_promoters.bed -wao > ${PairedAnch}/TSSoverlap/RAD21centred/reg_struc_pairs/${SET}_${filt}_paired_PsP_anchors_promoters.RAD21.2wayoverlap.bed
$BEDTOOLS intersect -b ${diffanchorsn2t}/H3K27overlap/CohesinAssEnhancers.RAD21centred.bed -a ${PairedAnch}/TSSoverlap/reg_struc_pairs/${SET}_${filt}paired_EsP_anchors_enhancers.bed -wao > ${PairedAnch}/TSSoverlap/RAD21centred/reg_struc_pairs/${SET}_${filt}_paired_EsP_anchors_enhancers.RAD21.2wayoverlap.bed
$BEDTOOLS intersect -b ${diffanchorsn2t}/H3K27overlap/CohesinAssEnhancers.RAD21centred.bed -a ${PairedAnch}/TSSoverlap/reg_struc_pairs/${SET}_${filt}paired_ES_anchors_enhancers.bed -wao > ${PairedAnch}/TSSoverlap/RAD21centred/reg_struc_pairs/${SET}_${filt}_paired_ES_anchors_enhancers.RAD21.2wayoverlap.bed
$BEDTOOLS intersect -b ${diffanchorsn2t}/H3K27overlap/CohesinAssEnhancers.RAD21centred.bed -a ${PairedAnch}/TSSoverlap/reg_struc_pairs/${SET}_${filt}paired_PS_anchors_promoters.bed -wao > ${PairedAnch}/TSSoverlap/RAD21centred/reg_struc_pairs/${SET}_${filt}_paired_PS_anchors_promoters.RAD21.2wayoverlap.bed
done
done
#for their partners without H3K27ac use all RAD21 pos.
for SET in ${MUTS}; do
for filt in ${FiltsFC1}; do
$BEDTOOLS intersect -b ${CHIPPEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.bed -a ${PairedAnch}/TSSoverlap/reg_struc_pairs/${SET}_${filt}paired_PsP_anchors_structuralPromoter.bed -wao > ${PairedAnch}/TSSoverlap/RAD21centred/reg_struc_pairs/${SET}_${filt}_paired_PsP_anchors_structuralPromoter.RAD21.2wayoverlap.bed
$BEDTOOLS intersect -b ${CHIPPEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.bed -a ${PairedAnch}/TSSoverlap/reg_struc_pairs/${SET}_${filt}paired_EsP_anchors_structuralPromoter.bed -wao > ${PairedAnch}/TSSoverlap/RAD21centred/reg_struc_pairs/${SET}_${filt}_paired_EsP_anchors_structuralPromoter.RAD21.2wayoverlap.bed
$BEDTOOLS intersect -b ${CHIPPEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.bed -a ${PairedAnch}/TSSoverlap/reg_struc_pairs/${SET}_${filt}paired_ES_anchors_structural.bed -wao > ${PairedAnch}/TSSoverlap/RAD21centred/reg_struc_pairs/${SET}_${filt}_paired_ES_anchors_structural.RAD21.2wayoverlap.bed
$BEDTOOLS intersect -b ${CHIPPEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.bed -a ${PairedAnch}/TSSoverlap/reg_struc_pairs/${SET}_${filt}paired_PS_anchors_structural.bed -wao > ${PairedAnch}/TSSoverlap/RAD21centred/reg_struc_pairs/${SET}_${filt}_paired_PS_anchors_structural.RAD21.2wayoverlap.bed
done
done
#### also intersect EP/EE/PP  again with H3K27ac --> useful for PSEA
mkdir ${PairedAnch}/TSSoverlap/H3K27accentred
for SET in ${MUTS}; do
for filt in ${FiltsFC1}; do
$BEDTOOLS intersect -b ${diffanchorsn2t}/H3K27overlap/CohesinAssEnhancers.bed -a ${PairedAnch}/TSSoverlap/${SET}_${filt}paired_EP_anchors_${SET}_enhancers.bed -wao > ${PairedAnch}/TSSoverlap/H3K27accentred/${SET}_${filt}_paired_EP_anchors_enhancers.H3K27ac.2wayoverlap.bed
$BEDTOOLS intersect -b ${diffanchorsn2t}/H3K27overlap/CohesinAssEnhancers.bed -a ${PairedAnch}/TSSoverlap/${SET}_${filt}paired_EP_anchors_${SET}_promoters.bed -wao > ${PairedAnch}/TSSoverlap/H3K27accentred/${SET}_${filt}_paired_EP_anchors_promoters.H3K27ac.2wayoverlap.bed
$BEDTOOLS intersect -b ${diffanchorsn2t}/H3K27overlap/CohesinAssEnhancers.bed -a ${PairedAnch}/TSSoverlap/${SET}_${filt}paired_PP_anchors_${SET}_promoters.plus.bed -wao > ${PairedAnch}/TSSoverlap/H3K27accentred/${SET}_${filt}_paired_PP_anchors_promoters.plus.H3K27ac.2wayoverlap.bed
$BEDTOOLS intersect -b ${diffanchorsn2t}/H3K27overlap/CohesinAssEnhancers.bed -a ${PairedAnch}/TSSoverlap/${SET}_${filt}paired_PP_anchors_${SET}_promoters.minus.bed -wao > ${PairedAnch}/TSSoverlap/H3K27accentred/${SET}_${filt}_paired_PP_anchors_promoters.minus.H3K27ac.2wayoverlap.bed
$BEDTOOLS intersect -b ${diffanchorsn2t}/H3K27overlap/CohesinAssEnhancers.bed -a ${PairedAnch}/TSSoverlap/${SET}_${filt}paired_EE_anchors_${SET}_enhancers.plus.bed -wao > ${PairedAnch}/TSSoverlap/H3K27accentred/${SET}_${filt}_paired_EE_anchors_enhancers.plus.H3K27ac.2wayoverlap.bed
$BEDTOOLS intersect -b ${diffanchorsn2t}/H3K27overlap/CohesinAssEnhancers.bed -a ${PairedAnch}/TSSoverlap/${SET}_${filt}paired_EE_anchors_${SET}_enhancers.minus.bed -wao > ${PairedAnch}/TSSoverlap/H3K27accentred/${SET}_${filt}_paired_EE_anchors_enhancers.minus.H3K27ac.2wayoverlap.bed
done
done
####also intersect TSS themselves with RAD21
$BEDTOOLS intersect -b ${diffanchorsn2t}/H3K27overlap/CohesinAssEnhancers.RAD21centred.bed -a ${diffanchorsn2t}/H3K27overlap/${TSSpos} -wao > ${PairedAnch}/TSSoverlap/RAD21centred/TSS.RAD21.2wayoverlap.bed
$BEDTOOLS intersect -a ${diffanchorsn2t}/H3K27overlap/CohesinAssEnhancers.RAD21centred.bed -b ${diffanchorsn2t}/H3K27overlap/${TSSpos} -wao > ${PairedAnch}/TSSoverlap/RAD21centred/TSS.RAD21.2wayoverlap.RAD.bed

#filter on the top RAD21 peak present in each anchor - this way each anchor appears only once in the plots later
##run the rmd script
cd ${WORKDIR}/scripts
rbioc_3-12 -e "rmarkdown::render('HIC_AML.diffEP_LoopAnchor_topRAD21peak.rmd',output_file='HIC_AML.diffEP_LoopAnchor_topRAD21peak.html')" #generates html report and output files




                          #####################################################################################################################
                          #    Deeptools heatmap clustering of paired anchors by TSS status EP EE PP - centred on topRAD21 peak               #
                          #####################################################################################################################

mkdir -p ${DTmatrixdir}/n2tresults/PairedAnch/TSSoverlap/
conda activate deepTools.3.3.0
#clustering of ChIP coverage at anchors - by TSS status EP EE PP - centred on topRAD21 peak - using deeptools 
###The bed files were filtered and rearranged in the HIC_AML.diffEP_LoopAnchor_topRAD21peak.rmd script 
#####top RAD21 peak in each anchor as center luster using deeptools with average bigwigs for RAD21/CTCF/H3K27ac

#"Left" anchors with sorting (descending) using reference-point mode
##compute for all "left" anchors separately and generate a sorted bed file as output
leftanchors="EP_anchors_promoters PP_anchors_promoters.plus EE_anchors_enhancers.plus"
for filt in ${FiltsFC1};do
for MUT in ${MUTS};do
for anchor in ${leftanchors};do
computeMatrix reference-point -S \
${MERGEDBW}/ave.AML_CTRL_SA1_CNVnorm.topQC.bigwig \
${MERGEDBW}/ave.AML_${MUT}_SA1_CNVnorm.topQC.bigwig \
${MERGEDBW}/ave.AML_CTRL_SA2_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_${MUT}_SA2_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_CTRL_Rad21_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_${MUT}_Rad21_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_CTRL_CTCF_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_${MUT}_CTCF_CNVnorm.bigwig \
${MERGEDBW_AML_H3K}/ave.AML_CTRL_H3K27ac_CNVnorm.bigwig \
${MERGEDBW_AML_H3K}/ave.AML_${MUT}_H3K27ac_CNVnorm.bigwig \
${ATACBW}/ave.CTRL.ATAC.goodQC.scaled.bigWig \
${ATACBW}/ave.${MUT}.ATAC.goodQC.scaled.bigWig \
-R ${PairedAnch}/TSSoverlap/TopRAD21centred/${MUT}_${filt}_${anchor}_RAD21top.bed \
--sortRegions descend --sortUsing mean \
-b 2000 -a 2000 -bs 10 -p 30 --referencePoint center --sortRegions descend --sortUsing mean \
--outFileSortedRegions ${PairedAnch}/TSSoverlap/TopRAD21centred/${MUT}_${filt}_${anchor}_RAD21top.sorted.mbw.bed \
-o ${DTmatrixdir}/n2tresults/PairedAnch/TSSoverlap/Matrix_2KBregions_Top_RAD21peak_${MUT}_${filt}_${anchor}.sorted.mbw.gz
done
done
done


##sort corresponding right side anchors according to left side clustering order in R
### run short R script to sort the anchors:
_DATE=$(date +%s)
cat >"${TMPDIR}/R.sort.bed.${_DATE}.R" <<EOF
DIR_DATA="/misc/data"
WORKDIR=file.path(DIR_DATA,"analysis/project_cohesin/Cohesin_AML/HiC/loops/differentialanchorsn2t/AnchorsOverlapPairs/TSSoverlap/TopRAD21centred")
MUTs<-c("SA2mut","RAD21mut")
Filts<-c("up.FC1", "down.FC1")
for (MUT in MUTs) {
  for (FILT in Filts){
#EP merged/av bw
leftsorted<-read.table(file=file.path(WORKDIR,paste0(MUT,"_",FILT,"_EP_anchors_promoters_RAD21top.sorted.mbw.bed")))
rightunsorted<-read.table(file=file.path(WORKDIR,paste0(MUT,"_",FILT,"_EP_anchors_enhancers_RAD21top.bed")))
rightsorted<-rightunsorted[match(leftsorted$V4 , rightunsorted$V4),]
write.table(rightsorted, file.path(WORKDIR,paste0(MUT,"_",FILT,"_EP_anchors_enhancers_RAD21top.sorted.mbw.bed")),row.names=F,quote=F,sep="\t",col.names=F)
#PP merged/av bw
leftsorted<-read.table(file=file.path(WORKDIR,paste0(MUT,"_",FILT,"_PP_anchors_promoters.plus_RAD21top.sorted.mbw.bed")))
rightunsorted<-read.table(file=file.path(WORKDIR,paste0(MUT,"_",FILT,"_PP_anchors_promoters.minus_RAD21top.bed")))
rightsorted<-rightunsorted[match(leftsorted$V4 , rightunsorted$V4),]
write.table(rightsorted, file.path(WORKDIR,paste0(MUT,"_",FILT,"_PP_anchors_promoters.minus_RAD21top.sorted.mbw.bed")),row.names=F,quote=F,sep="\t",col.names=F)
#EE merged/av bw
leftsorted<-read.table(file=file.path(WORKDIR,paste0(MUT,"_",FILT,"_EE_anchors_enhancers.plus_RAD21top.sorted.mbw.bed")))
rightunsorted<-read.table(file=file.path(WORKDIR,paste0(MUT,"_",FILT,"_EE_anchors_enhancers.minus_RAD21top.bed")))
rightsorted<-rightunsorted[match(leftsorted$V4 , rightunsorted$V4),]
write.table(rightsorted, file.path(WORKDIR,paste0(MUT,"_",FILT,"_EE_anchors_enhancers.minus_RAD21top.sorted.mbw.bed")),row.names=F,quote=F,sep="\t",col.names=F)
}}
EOF
chmod 750 "${TMPDIR}/R.sort.bed.${_DATE}.R"
R < ${TMPDIR}/R.sort.bed.${_DATE}.R  --no-save
rm ${TMPDIR}/R.sort.bed.${_DATE}.R



## Right side anchors, cluster according to order determined for left side
rightanchors="EP_anchors_enhancers PP_anchors_promoters.minus EE_anchors_enhancers.minus"

for filt in ${FiltsFC1};do
for MUT in ${MUTS};do
for anchor in ${rightanchors};do
computeMatrix reference-point -S \
${MERGEDBW}/ave.AML_CTRL_SA1_CNVnorm.topQC.bigwig \
${MERGEDBW}/ave.AML_${MUT}_SA1_CNVnorm.topQC.bigwig \
${MERGEDBW}/ave.AML_CTRL_SA2_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_${MUT}_SA2_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_CTRL_Rad21_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_${MUT}_Rad21_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_CTRL_CTCF_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_${MUT}_CTCF_CNVnorm.bigwig \
/ave.AML_CTRL_H3K27ac_CNVnorm.bigwig \
${MERGEDBW_AML_H3K}/ave.AML_${MUT}_H3K27ac_CNVnorm.bigwig \
${ATACBW}/ave.CTRL.ATAC.goodQC.scaled.bigWig \
${ATACBW}/ave.${MUT}.ATAC.goodQC.scaled.bigWig \
-R ${PairedAnch}/TSSoverlap/TopRAD21centred/${MUT}_${filt}_${anchor}_RAD21top.sorted.mbw.bed \
-b 2000 -a 2000 -bs 10 -p 30 --referencePoint center \
-o ${DTmatrixdir}/n2tresults/PairedAnch/TSSoverlap/Matrix_2KBregions_Top_RAD21peak_${MUT}_${filt}_${anchor}.sorted.mbw.gz
done
done
done

for filt in ${FiltsFC1};do
for MUT in ${MUTS};do
for anchor in ${rightanchors};do
wc -l ${PairedAnch}/TSSoverlap/TopRAD21centred/${MUT}_${filt}_${anchor}_RAD21top.sorted.bed
done
done
done

#plot
set_scale=(
'up.FC1' 6 'SA2mut' 'EP_anchors_enhancers'
'up.FC1' 6 'SA2mut' 'EP_anchors_promoters'
'down.FC1' 11.5 'SA2mut' 'EP_anchors_enhancers'
'down.FC1' 11.5 'SA2mut' 'EP_anchors_promoters'
'up.FC1' 5 'SA2mut' 'EE_anchors_enhancers.plus'
'up.FC1' 5 'SA2mut' 'EE_anchors_enhancers.minus'
'down.FC1' 7.4 'SA2mut' 'EE_anchors_enhancers.plus'
'down.FC1' 7.4 'SA2mut' 'EE_anchors_enhancers.minus'
'up.FC1' 4 'SA2mut' 'PP_anchors_promoters.plus'
'up.FC1' 4 'SA2mut' 'PP_anchors_promoters.minus'
'down.FC1' 7 'SA2mut' 'PP_anchors_promoters.plus'
'down.FC1' 7 'SA2mut' 'PP_anchors_promoters.minus'
'up.FC1' 4 'RAD21mut' 'EP_anchors_enhancers'
'up.FC1' 4 'RAD21mut' 'EP_anchors_promoters'
'down.FC1' 9 'RAD21mut' 'EP_anchors_enhancers'
'down.FC1' 9 'RAD21mut' 'EP_anchors_promoters'
'up.FC1' 4 'RAD21mut' 'EE_anchors_enhancers.plus'
'up.FC1' 4 'RAD21mut' 'EE_anchors_enhancers.minus'
'down.FC1' 7 'RAD21mut' 'EE_anchors_enhancers.plus'
'down.FC1' 7 'RAD21mut' 'EE_anchors_enhancers.minus'
'up.FC1' 4 'RAD21mut' 'PP_anchors_promoters.plus'
'up.FC1' 4 'RAD21mut' 'PP_anchors_promoters.minus'
'down.FC1' 4 'RAD21mut' 'PP_anchors_promoters.plus'
'down.FC1' 4 'RAD21mut' 'PP_anchors_promoters.minus'
)
#generate the heatmaps
mkdir -p ${FIGURESDIR}/Heatmaps/n2tresults/PairedAnchors/RAD21centred/byTSSstatus/
ymaxvals="26 30 32" #try different maximum y axis values
for ymax in ${ymaxvals};do
mkdir ${FIGURESDIR}/Heatmaps/n2tresults/PairedAnchors/RAD21centred/byTSSstatus/ymax${ymax}
done

for ymax in ${ymaxvals};do
for (( idx=0 ; idx<${#set_scale[@]} ; idx+=4 )) ; do
    filt=${set_scale[idx]}
    scale=${set_scale[idx+1]}
    MUT=${set_scale[idx+2]}
    anchor=${set_scale[idx+3]}
    echo "filt=$filt"
    echo "scale=$scale"
    echo "MUT=$MUT"
    echo "anchor=$anchor"
    echo
    echo "creating matrix for loop ${anchor} ${filt} ${MUT} with height ${scale}"
Rlab="${anchor}"
Xlab=""
plotHeatmap -m ${DTmatrixdir}/n2tresults/PairedAnch/TSSoverlap/Matrix_2KBregions_Top_RAD21peak_${MUT}_${filt}_${anchor}.sorted.mbw.gz --sortRegions keep \
--colorList 'white,darkgoldenrod' 'white,darkgoldenrod' \
'white,lime' 'white,lime' \
'white,mediumvioletred' 'white,mediumvioletred' \
'white,darkblue' 'white,darkblue' \
'white,orange' 'white,orange' \
'white,salmon' 'white,salmon' \
--colorNumber 10 \
--heatmapHeight $scale \
--yMax ${ymax} \
--samplesLabel \
CTRL ${MUT} CTRL ${MUT} CTRL ${MUT} CTRL ${MUT} CTRL ${MUT} CTRL ${MUT} \
--whatToShow 'plot and heatmap' \
--missingDataColor 1 \
--legendLocation none \
--regionsLabel "${Rlab}" \
--xAxisLabel "${Xlab}" --startLabel '' --endLabel '' \
-o ${FIGURESDIR}/Heatmaps/n2tresults/PairedAnchors/RAD21centred/byTSSstatus/ymax${ymax}/ChIP_AML_CTRLvs${MUT}_Matrix_2KBregions_Top_RAD21peak_${filt}_${anchor}.mbw.pdf
done
done
done


                          ######################################################################################################################################
                          #    Deeptools heatmap clustering of paired regulatory-structural (RS) diff. anchors  - centred on topRAD21 peak                          #
                          ######################################################################################################################################

###################################clustering of ChIP coverage at mixed anchor pairs - by TSS/H3K27ac status ES PS PsP EsP - centred on topRAD21 peak - using deeptools 
###The files were filtered and rearranged in the rmd script 

#################################################top RAD21 peak in each anchor as center
####now cluster using deeptools with scaled bigwigs for RAD21/CTCF/H3K27ac
#"Left" anchors with sorting (descending) using reference-point mode
##approach 1: all left anchors separately
leftanchorsRS="ES_anchors_enhancers PS_anchors_promoters PsP_anchors_promoters EsP_anchors_enhancers"

for filt in ${FiltsFC1};do
for MUT in ${MUTS};do
for anchor in ${leftanchorsRS};do
computeMatrix reference-point -S \
${MERGEDBW}/ave.AML_CTRL_SA1_CNVnorm.topQC.bigwig \
${MERGEDBW}/ave.AML_${MUT}_SA1_CNVnorm.topQC.bigwig \
${MERGEDBW}/ave.AML_CTRL_SA2_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_${MUT}_SA2_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_CTRL_Rad21_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_${MUT}_Rad21_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_CTRL_CTCF_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_${MUT}_CTCF_CNVnorm.bigwig \
${MERGEDBW_AML_H3K}/ave.AML_CTRL_H3K27ac_CNVnorm.bigwig \
${MERGEDBW_AML_H3K}/ave.AML_${MUT}_H3K27ac_CNVnorm.bigwig \
${ATACBW}/ave.CTRL.ATAC.goodQC.scaled.bigWig \
${ATACBW}/ave.${MUT}.ATAC.goodQC.scaled.bigWig \
-R ${PairedAnch}/TSSoverlap/TopRAD21centred/reg_struc_pairs/${MUT}_${filt}_${anchor}_RAD21top.bed \
--sortRegions descend --sortUsing mean \
-b 2000 -a 2000 -bs 10 -p 30 --referencePoint center --sortRegions descend --sortUsing mean \
--outFileSortedRegions ${PairedAnch}/TSSoverlap/TopRAD21centred/reg_struc_pairs/${MUT}_${filt}_${anchor}_RAD21top.sorted.mbw.bed \
-o ${DTmatrixdir}/n2tresults/PairedAnch/TSSoverlap/Matrix_2KBregions_Top_RAD21peak_${MUT}_${filt}_${anchor}.sorted.mbw.gz
done done done



####sort corresponding right side anchors according to left side clustering order in R

_DATE=$(date +%s)
cat >"${TMPDIR}/R.sort.bed.${_DATE}.R" <<EOF
DIR_DATA="/misc/data"
WORKDIR=file.path(DIR_DATA,"analysis/project_cohesin/Cohesin_AML/HiC/loops/differentialanchorsn2t/AnchorsOverlapPairs/TSSoverlap/TopRAD21centred/reg_struc_pairs")
MUTs<-c("SA2mut","RAD21mut")
Filts<-c("up.FC1", "down.FC1")
for (MUT in MUTs) {
  for (FILT in Filts){
#ES merged/av bw
leftsorted<-read.table(file=file.path(WORKDIR,paste0(MUT,"_",FILT,"_ES_anchors_enhancers_RAD21top.sorted.mbw.bed")))
rightunsorted<-read.table(file=file.path(WORKDIR,paste0(MUT,"_",FILT,"_ES_anchors_structural_RAD21top.bed")))
rightsorted<-rightunsorted[match(leftsorted$V4 , rightunsorted$V4),]
write.table(rightsorted, file.path(WORKDIR,paste0(MUT,"_",FILT,"_ES_anchors_structural_RAD21top.sorted.mbw.bed")),row.names=F,quote=F,sep="\t",col.names=F)
#PS merged/av bw
leftsorted<-read.table(file=file.path(WORKDIR,paste0(MUT,"_",FILT,"_PS_anchors_promoters_RAD21top.sorted.mbw.bed")))
rightunsorted<-read.table(file=file.path(WORKDIR,paste0(MUT,"_",FILT,"_PS_anchors_structural_RAD21top.bed")))
rightsorted<-rightunsorted[match(leftsorted$V4 , rightunsorted$V4),]
write.table(rightsorted, file.path(WORKDIR,paste0(MUT,"_",FILT,"_PS_anchors_structural_RAD21top.sorted.mbw.bed")),row.names=F,quote=F,sep="\t",col.names=F)
#EsP merged/av bw
leftsorted<-read.table(file=file.path(WORKDIR,paste0(MUT,"_",FILT,"_EsP_anchors_enhancers_RAD21top.sorted.mbw.bed")))
rightunsorted<-read.table(file=file.path(WORKDIR,paste0(MUT,"_",FILT,"_EsP_anchors_structuralPromoter_RAD21top.bed")))
rightsorted<-rightunsorted[match(leftsorted$V4 , rightunsorted$V4),]
write.table(rightsorted, file.path(WORKDIR,paste0(MUT,"_",FILT,"_EsP_anchors_structuralPromoter_RAD21top.sorted.mbw.bed")),row.names=F,quote=F,sep="\t",col.names=F)
#PsP merged/av bw
leftsorted<-read.table(file=file.path(WORKDIR,paste0(MUT,"_",FILT,"_PsP_anchors_promoters_RAD21top.sorted.mbw.bed")))
rightunsorted<-read.table(file=file.path(WORKDIR,paste0(MUT,"_",FILT,"_PsP_anchors_structuralPromoter_RAD21top.bed")))
rightsorted<-rightunsorted[match(leftsorted$V4 , rightunsorted$V4),]
write.table(rightsorted, file.path(WORKDIR,paste0(MUT,"_",FILT,"_PsP_anchors_structuralPromoter_RAD21top.sorted.mbw.bed")),row.names=F,quote=F,sep="\t",col.names=F)
}}
EOF
chmod 750 "${TMPDIR}/R.sort.bed.${_DATE}.R"
R < ${TMPDIR}/R.sort.bed.${_DATE}.R  --no-save
rm ${TMPDIR}/R.sort.bed.${_DATE}.R



##Right side anchors, cluster according to order determined for left side - approach 1: all right anchors separately
rightanchorsRS="ES_anchors_structural PS_anchors_structural PsP_anchors_structuralPromoter EsP_anchors_structuralPromoter"
for filt in ${FiltsFC1};do
for MUT in ${MUTS};do
for anchor in ${rightanchorsRS};do
computeMatrix reference-point -S \
${MERGEDBW}/ave.AML_CTRL_SA1_CNVnorm.topQC.bigwig \
${MERGEDBW}/ave.AML_${MUT}_SA1_CNVnorm.topQC.bigwig \
${MERGEDBW}/ave.AML_CTRL_SA2_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_${MUT}_SA2_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_CTRL_Rad21_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_${MUT}_Rad21_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_CTRL_CTCF_CNVnorm.bigwig \
${MERGEDBW}/ave.AML_${MUT}_CTCF_CNVnorm.bigwig \
${MERGEDBW_AML_H3K}/ave.AML_CTRL_H3K27ac_CNVnorm.bigwig \
${MERGEDBW_AML_H3K}/ave.AML_${MUT}_H3K27ac_CNVnorm.bigwig \
${ATACBW}/ave.CTRL.ATAC.goodQC.scaled.bigWig \
${ATACBW}/ave.${MUT}.ATAC.goodQC.scaled.bigWig \
-R ${PairedAnch}/TSSoverlap/TopRAD21centred/reg_struc_pairs/${MUT}_${filt}_${anchor}_RAD21top.sorted.mbw.bed \
-b 2000 -a 2000 -bs 10 -p 30 --referencePoint center \
-o ${DTmatrixdir}/n2tresults/PairedAnch/TSSoverlap/Matrix_2KBregions_Top_RAD21peak_${MUT}_${filt}_${anchor}.sorted.mbw.gz
done
done
done

cd ${PairedAnch}/TSSoverlap/TopRAD21centred/reg_struc_pairs/
for filt in ${FiltsFC1};do
for MUT in ${MUTS};do
for anchor in ${leftanchorsRS};do
wc -l ${MUT}_${filt}_${anchor}_RAD21top.sorted.mbw.bed
done done done


#plot
set_scale=(
'up.FC1' 13.4 'SA2mut' 'ES_anchors_enhancers'
'up.FC1' 13.4 'SA2mut' 'ES_anchors_structural'
'up.FC1' 8.1 'SA2mut' 'PS_anchors_promoters'
'up.FC1' 8.1 'SA2mut' 'PS_anchors_structural'
'up.FC1' 4 'SA2mut' 'PsP_anchors_promoters'
'up.FC1' 4 'SA2mut' 'PsP_anchors_structuralPromoter'
'up.FC1' 4 'SA2mut' 'EsP_anchors_enhancers'
'up.FC1' 4 'SA2mut' 'EsP_anchors_structuralPromoter'
'down.FC1' 25 'SA2mut' 'ES_anchors_enhancers'
'down.FC1' 25 'SA2mut' 'ES_anchors_structural'
'down.FC1' 18 'SA2mut' 'PS_anchors_promoters'
'down.FC1' 18 'SA2mut' 'PS_anchors_structural'
'down.FC1' 4 'SA2mut' 'PsP_anchors_promoters'
'down.FC1' 4 'SA2mut' 'PsP_anchors_structuralPromoter'
'down.FC1' 4 'SA2mut' 'EsP_anchors_enhancers'
'down.FC1' 4 'SA2mut' 'EsP_anchors_structuralPromoter'
'up.FC1' 6 'RAD21mut' 'ES_anchors_enhancers'
'up.FC1' 6 'RAD21mut' 'ES_anchors_structural'
'up.FC1' 4 'RAD21mut' 'PS_anchors_promoters'
'up.FC1' 4 'RAD21mut' 'PS_anchors_structural'
'up.FC1' 4 'RAD21mut' 'PsP_anchors_promoters'
'up.FC1' 4 'RAD21mut' 'PsP_anchors_structuralPromoter'
'up.FC1' 4 'RAD21mut' 'EsP_anchors_enhancers'
'up.FC1' 4 'RAD21mut' 'EsP_anchors_structuralPromoter'
'down.FC1' 24 'RAD21mut' 'ES_anchors_enhancers'
'down.FC1' 24 'RAD21mut' 'ES_anchors_structural'
'down.FC1' 13 'RAD21mut' 'PS_anchors_promoters'
'down.FC1' 13 'RAD21mut' 'PS_anchors_structural'
'down.FC1' 4 'RAD21mut' 'PsP_anchors_promoters'
'down.FC1' 4 'RAD21mut' 'PsP_anchors_structuralPromoter'
'down.FC1' 4 'RAD21mut' 'EsP_anchors_enhancers'
'down.FC1' 4 'RAD21mut' 'EsP_anchors_structuralPromoter'
)
#generate the heatmaps

mkdir ${FIGURESDIR}/Heatmaps/n2tresults/PairedAnchors/RAD21centred/byTSSstatus/struct_reg_diffloops/
mkdir ${FIGURESDIR}/Heatmaps/n2tresults/PairedAnchors/RAD21centred/byTSSstatus/struct_reg_diffloops/mbw
mkdir ${FIGURESDIR}/Heatmaps/n2tresults/PairedAnchors/RAD21centred/byTSSstatus/struct_reg_diffloops/mbw/ymax26/
mkdir ${FIGURESDIR}/Heatmaps/n2tresults/PairedAnchors/RAD21centred/byTSSstatus/struct_reg_diffloops/mbw/ymax30/
mkdir ${FIGURESDIR}/Heatmaps/n2tresults/PairedAnchors/RAD21centred/byTSSstatus/struct_reg_diffloops/mbw/ymax35/

ymaxvals="26 30 35"
for ymax in ${ymaxvals};do
for tracktype in ${tracktypes};do
for (( idx=0 ; idx<${#set_scale[@]} ; idx+=4 )) ; do
    filt=${set_scale[idx]}
    scale=${set_scale[idx+1]}
    MUT=${set_scale[idx+2]}
    anchor=${set_scale[idx+3]}
    echo "filt=$filt"
    echo "scale=$scale"
    echo "MUT=$MUT"
    echo "anchor=$anchor"
    echo
    echo "creating matrix for loop ${anchor} ${filt} ${MUT} with height ${scale}"
Rlab="${anchor}"
Xlab=""
plotHeatmap -m ${DTmatrixdir}/n2tresults/PairedAnch/TSSoverlap/Matrix_2KBregions_Top_RAD21peak_${MUT}_${filt}_${anchor}.sorted.mbw.gz --sortRegions keep \
--colorList 'white,darkgoldenrod' 'white,darkgoldenrod' \
'white,lime' 'white,lime' \
'white,mediumvioletred' 'white,mediumvioletred' \
'white,darkblue' 'white,darkblue' \
'white,orange' 'white,orange' \
'white,salmon' 'white,salmon' \
--colorNumber 10 \
--heatmapHeight $scale \
--yMax ${ymax} \
--samplesLabel \
CTRL ${MUT} CTRL ${MUT} CTRL ${MUT} CTRL ${MUT} CTRL ${MUT} CTRL ${MUT} \
--whatToShow 'plot and heatmap' \
--missingDataColor 1 \
--legendLocation none \
--regionsLabel "${Rlab}" \
--xAxisLabel "${Xlab}" --startLabel '' --endLabel '' \
-o ${FIGURESDIR}/Heatmaps/n2tresults/PairedAnchors/RAD21centred/byTSSstatus/struct_reg_diffloops/mbw/ymax${ymax}/ChIP_AML_CTRLvs${MUT}_Matrix_2KBregions_Top_RAD21peak_${filt}_${anchor}.mbw.pdf
done
done
done



                        ################################################################################
                        # ChIPcoverage of RAD21/H3K27ac/STAG split by Enhancer of Promoter ID
                        ### ignoring paired/unpaired status: this will result in more peaks in total
                        ################################################################################
mkdir ${FIGURESDIR}/ghist/n2tresults/EPanchors
mkdir ${FIGURESDIR}/ghist/n2tresults/EPanchors/avbw
mkdir ${FIGURESDIR}/ghist/n2tresults/EPanchors/avbw/STAG
types="Enhancer Promoter"
for filt in ${FiltsFC1}; do
for MUT in ${MUTS}; do
 for type in ${types}; do
     cat >"${TMPDIR}/ghist.${MUT}.${filt}.${type}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd ${diffanchorsn2t}/H3K27overlap/TSSanchorOverlap/peaksSplitByAnchorTSSstatus
annotatePeaks.pl ${type}_RAD21peaks.${filt}_${MUT}.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_CTRL_RAD21.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/n2tresults/EPanchors/avbw/CTRL.${type}_RAD21peaks.${filt}_${MUT}.RAD21.ghist.txt"
annotatePeaks.pl ${type}_RAD21peaks.${filt}_${MUT}.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_${MUT}_RAD21.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/n2tresults/EPanchors/avbw/${MUT}.${type}_RAD21peaks.${filt}_${MUT}.RAD21.ghist.txt"
annotatePeaks.pl ${type}_RAD21peaks.${filt}_${MUT}.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW_AML_H3K}/ave.AML_CTRL_H3K27ac.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/n2tresults/EPanchors/avbw/CTRL.${type}_RAD21peaks.${filt}_${MUT}.H3K27ac.ghist.txt"
annotatePeaks.pl ${type}_RAD21peaks.${filt}_${MUT}.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW_AML_H3K}/ave.AML_${MUT}_H3K27ac.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/n2tresults/EPanchors/avbw/${MUT}.${type}_RAD21peaks.${filt}_${MUT}.H3K27ac.ghist.txt"

annotatePeaks.pl ${type}_RAD21peaks.${filt}_${MUT}.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_CTRL_SA1.CNVnorm.topQC.bedGraph > "${FIGURESDIR}/ghist/n2tresults/EPanchors/avbw/STAG/CTRL.${type}_RAD21peaks.${filt}_${MUT}.SA1.ghist.txt"
annotatePeaks.pl ${type}_RAD21peaks.${filt}_${MUT}.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_${MUT}_SA1.CNVnorm.topQC.bedGraph > "${FIGURESDIR}/ghist/n2tresults/EPanchors/avbw/STAG/${MUT}.${type}_RAD21peaks.${filt}_${MUT}.SA1.ghist.txt"
annotatePeaks.pl ${type}_RAD21peaks.${filt}_${MUT}.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_CTRL_SA2.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/n2tresults/EPanchors/avbw/STAG/CTRL.${type}_RAD21peaks.${filt}_${MUT}.SA2.ghist.txt"
annotatePeaks.pl ${type}_RAD21peaks.${filt}_${MUT}.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_${MUT}_SA2.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/n2tresults/EPanchors/avbw/STAG/${MUT}.${type}_RAD21peaks.${filt}_${MUT}.SA2.ghist.txt"
EOF
    chmod 750 "${TMPDIR}/ghist.${MUT}.${filt}.${type}.${_DATE}.sh"
	echo "plotting ghists for ${type} in ${MUT} ${filt} loop anchors"
	screen -dm -S ${MUT}.${filt}.${type} bash -c "bash ${TMPDIR}/ghist.${MUT}.${filt}.${type}.${_DATE}.sh"
 done
done
done

#RAD21/H3K27ac plots
for filt in ${FiltsFC1};do
for type in ${types}; do
plotHIST.sh -g "SA2mut.${type}_RAD21peaks.${filt}_SA2mut.RAD21.ghist.txt CTRL.${type}_RAD21peaks.${filt}_SA2mut.RAD21.ghist.txt" \
-s "STAG2-mut CTRL-AML" -c "seagreen3 firebrick" -x 1000 -y "0 20" -d ${FIGURESDIR}/ghist/n2tresults/EPanchors/avbw/ -n RAD21.SA2mutvsCTRL.${type}_RAD21peaks.${filt}
plotHIST.sh -g "SA2mut.${type}_RAD21peaks.${filt}_SA2mut.H3K27ac.ghist.txt CTRL.${type}_RAD21peaks.${filt}_SA2mut.H3K27ac.ghist.txt" \
-s "STAG2-mut CTRL-AML" -c "seagreen3 firebrick" -x 2000 -y "0 25" -d ${FIGURESDIR}/ghist/n2tresults/EPanchors/avbw/ -n H3K27ac.SA2mutvsCTRL.${type}_RAD21peaks.${filt}

plotHIST.sh -g "RAD21mut.${type}_RAD21peaks.${filt}_RAD21mut.RAD21.ghist.txt CTRL.${type}_RAD21peaks.${filt}_RAD21mut.RAD21.ghist.txt" \
-s "RAD21-mut CTRL-AML" -c "mediumvioletred firebrick" -x 1000 -y "0 20" -d ${FIGURESDIR}/ghist/n2tresults/EPanchors/avbw/ -n RAD21.RAD21mutvsCTRL.${type}_RAD21peaks.${filt}
plotHIST.sh -g "RAD21mut.${type}_RAD21peaks.${filt}_RAD21mut.H3K27ac.ghist.txt CTRL.${type}_RAD21peaks.${filt}_RAD21mut.H3K27ac.ghist.txt" \
-s "RAD21-mut CTRL-AML" -c "mediumvioletred firebrick" -x 2000 -y "0 25" -d ${FIGURESDIR}/ghist/n2tresults/EPanchors/avbw/ -n H3K27ac.RAD21mutvsCTRL.${type}_RAD21peaks.${filt}
done done

#SA1/SA2 plots
for filt in ${FiltsFC1};do
for type in ${types}; do
plotHIST.sh -g "SA2mut.${type}_RAD21peaks.${filt}_SA2mut.SA2.ghist.txt SA2mut.${type}_RAD21peaks.${filt}_SA2mut.SA1.ghist.txt" \
-s "STAG2 STAG1" -c "springgreen2 darkgoldenrod2" -x 1000 -y "0 12" -d ${FIGURESDIR}/ghist/n2tresults/EPanchors/avbw/STAG -n SA2mut.SA1vsSA2.${type}_RAD21peaks.${filt}
plotHIST.sh -g "CTRL.${type}_RAD21peaks.${filt}_SA2mut.SA2.ghist.txt CTRL.${type}_RAD21peaks.${filt}_SA2mut.SA1.ghist.txt" \
-s "STAG2 STAG1" -c "springgreen2 darkgoldenrod2" -x 1000 -y "0 12" -d ${FIGURESDIR}/ghist/n2tresults/EPanchors/avbw/STAG -n CTRL.SA1vsSA2.${type}_RAD21peaks.${filt}
done done


                        ################################################################################
                        #            Motif analyist all Promoter/Enhancer diffloop RAD21 pos combined  #
                        #                     ignoring paired/unpaired status                          #
                        ################################################################################

###############################################run motif analysis for EP-anchors with cohesin
motifanchdir=${diffanchorsn2t}/AnchorChIPoverlap/motifs
for filt in ${FiltsFC1}; do
for MUT in ${MUTS}; do
 for type in ${types}; do
 cat >"${TMPDIR}/motifs.${type}${MUT}${filt}.${_DATE}.sh" <<EOF
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
#find and filter motifs
cd ${diffanchorsn2t}/H3K27overlap/TSSanchorOverlap/peaksSplitByAnchorTSSstatus
findMotifsGenome.pl ${type}_RAD21peaks.${filt}_${MUT}.bed hg38 ${motifanchdir}/${type}_RAD21peaks.${filt}_${MUT} -size given -len 7,8,9,10,11,12,13,14 -p 6 -h
compareMotifs.pl ${motifanchdir}/${type}_RAD21peaks.${filt}_${MUT}/homerMotifs.all.motifs ${motifanchdir}/${type}_RAD21peaks.${filt}_${MUT}/final -reduceThresh .75 -matchThresh .6 -pvalue 1e-12 -info 1.5 -cpu 6
EOF
chmod 750 "${TMPDIR}/motifs.${type}${MUT}${filt}.${_DATE}.sh"
	echo "motifs for ${type} ${MUT} ${filt}"
	screen -dm -S ${type}${MUT}${filt} bash -c "bash ${TMPDIR}/motifs.${type}${MUT}${filt}.${_DATE}.sh"
done
done
done



################################################################################
# ChIPcoverage ofSTAG at all Promoter/Enhancer diffloop pos combined
### ignoring paired/unpaired status: this will result in more peaks in total
################################################################################

for filt in ${FiltsFC1}; do
for MUT in ${MUTS}; do
     cat >"${TMPDIR}/ghist.${MUT}.${filt}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd ${diffanchorsn2t}/H3K27overlap/
annotatePeaks.pl ${MUT}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.RAD21centred.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_CTRL_SA1.CNVnorm.topQC.bedGraph > "${FIGURESDIR}/ghist/n2tresults/EPanchors/avbw/STAG/CTRL.allEP_RAD21peaks.${filt}_${MUT}.SA1.ghist.txt"
annotatePeaks.pl ${MUT}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.RAD21centred.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_${MUT}_SA1.CNVnorm.topQC.bedGraph > "${FIGURESDIR}/ghist/n2tresults/EPanchors/avbw/STAG/${MUT}.allEP_RAD21peaks.${filt}_${MUT}.SA1.ghist.txt"
annotatePeaks.pl ${MUT}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.RAD21centred.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_CTRL_SA2.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/n2tresults/EPanchors/avbw/STAG/CTRL.allEP_RAD21peaks.${filt}_${MUT}.SA2.ghist.txt"
annotatePeaks.pl ${MUT}vsCTRL.loopAnchors.red.${filt}.CohesinAssEnhancers.RAD21centred.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_${MUT}_SA2.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/n2tresults/EPanchors/avbw/STAG/${MUT}.allEP_RAD21peaks.${filt}_${MUT}.SA2.ghist.txt"
EOF
    chmod 750 "${TMPDIR}/ghist.${MUT}.${filt}.${_DATE}.sh"
	echo "plotting ghists for in ${MUT} ${filt} all EP loop anchors"
	screen -dm -S ${MUT}.${filt} bash -c "bash ${TMPDIR}/ghist.${MUT}.${filt}.${_DATE}.sh"
 done
done
done


#SA1/SA2 plots
for filt in ${FiltsFC1};do
plotHIST.sh -g "SA2mut.allEP_RAD21peaks.${filt}_SA2mut.SA2.ghist.txt SA2mut.allEP_RAD21peaks.${filt}_SA2mut.SA1.ghist.txt" \
-s "STAG2 STAG1" -c "springgreen2 darkgoldenrod2" -x 1000 -y "0 10" -d ${FIGURESDIR}/ghist/n2tresults/EPanchors/avbw/STAG -n SA2mut.SA1vsSA2.allEP_RAD21peaks.${filt}
plotHIST.sh -g "CTRL.allEP_RAD21peaks.${filt}_SA2mut.SA2.ghist.txt CTRL.allEP_RAD21peaks.${filt}_SA2mut.SA1.ghist.txt" \
-s "STAG2 STAG1" -c "springgreen2 darkgoldenrod2" -x 1000 -y "0 10" -d ${FIGURESDIR}/ghist/n2tresults/EPanchors/avbw/STAG -n CTRL.SA1vsSA2.allEP_RAD21peaks.${filt}
done



################################################################################
# ChIPcoverage ofSTAG at all structural diffloop pos combined
### ignoring paired/unpaired status: this will result in more peaks in total
################################################################################
mkdir ${FIGURESDIR}/ghist/n2tresults/strucANCH/
mkdir ${FIGURESDIR}/ghist/n2tresults/strucANCH/avbw
mkdir ${FIGURESDIR}/ghist/n2tresults/strucANCH/avbw/STAG
for filt in ${FiltsFC1}; do
for MUT in ${MUTS}; do
     cat >"${TMPDIR}/ghist.${MUT}.${filt}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
#RAD21
annotatePeaks.pl ${STRUCLOOPDIR}/${MUT}vsCTRL.loopAnchors.red.${filt}.structural.RAD21.centred.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_CTRL_RAD21.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/n2tresults/strucANCH/avbw/CTRL.RAD21peaks.strucLoop.${filt}_${MUT}.RAD21.ghist.txt"
annotatePeaks.pl ${STRUCLOOPDIR}/${MUT}vsCTRL.loopAnchors.red.${filt}.structural.RAD21.centred.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_${MUT}_RAD21.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/n2tresults/strucANCH/avbw/${MUT}.RAD21peaks.strucLoop.${filt}_${MUT}.RAD21.ghist.txt"
#STAGS
annotatePeaks.pl ${STRUCLOOPDIR}/${MUT}vsCTRL.loopAnchors.red.${filt}.structural.RAD21.centred.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_CTRL_SA1.CNVnorm.topQC.bedGraph > "${FIGURESDIR}/ghist/n2tresults/strucANCH/avbw/STAG/CTRL.RAD21peaks.strucLoop.${filt}_${MUT}.SA1.ghist.txt"
annotatePeaks.pl ${STRUCLOOPDIR}/${MUT}vsCTRL.loopAnchors.red.${filt}.structural.RAD21.centred.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_${MUT}_SA1.CNVnorm.topQC.bedGraph > "${FIGURESDIR}/ghist/n2tresults/strucANCH/avbw/STAG/${MUT}.RAD21peaks.strucLoop.${filt}_${MUT}.SA1.ghist.txt"
annotatePeaks.pl ${STRUCLOOPDIR}/${MUT}vsCTRL.loopAnchors.red.${filt}.structural.RAD21.centred.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_CTRL_SA2.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/n2tresults/strucANCH/avbw/STAG/CTRL.RAD21peaks.strucLoop.${filt}_${MUT}.SA2.ghist.txt"
annotatePeaks.pl ${STRUCLOOPDIR}/${MUT}vsCTRL.loopAnchors.red.${filt}.structural.RAD21.centred.bed hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_${MUT}_SA2.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/n2tresults/strucANCH/avbw/STAG/${MUT}.RAD21peaks.strucLoop.${filt}_${MUT}.SA2.ghist.txt"
EOF
    chmod 750 "${TMPDIR}/ghist.${MUT}.${filt}.${_DATE}.sh"
	echo "plotting ghists for in ${MUT} ${filt} all structural loop anchors"
	screen -dm -S ${MUT}.${filt} bash -c "bash ${TMPDIR}/ghist.${MUT}.${filt}.${_DATE}.sh"
 done
done


#RAD21 plots: CTRL vs MUT
for MUT in ${MUTS}; do
for filt in ${FiltsFC1};do
plotHIST.sh -g "CTRL.RAD21peaks.strucLoop.${filt}_SA2mut.RAD21.ghist.txt SA2mut.RAD21peaks.strucLoop.${filt}_SA2mut.RAD21.ghist.txt" \
-s "CTRL-AML STAG2-mut" -c "firebrick seagreen3" -x 1000 -y "0 28" -d ${FIGURESDIR}/ghist/n2tresults/strucANCH/avbw -n SA2mutvsCTRL.RAD21.structural_RAD21peaks.${filt}

plotHIST.sh -g "CTRL.RAD21peaks.strucLoop.${filt}_RAD21mut.RAD21.ghist.txt RAD21mut.RAD21peaks.strucLoop.${filt}_RAD21mut.RAD21.ghist.txt" \
-s "CTRL-AML RAD21-mut" -c "firebrick mediumvioletred" -x 1000 -y "0 28" -d ${FIGURESDIR}/ghist/n2tresults/strucANCH/avbw -n RAD21mutvsCTRL.RAD21.structural_RAD21peaks.${filt}
done
done


#SA1/SA2 plots: within group
for filt in ${FiltsFC1};do
plotHIST.sh -g "SA2mut.RAD21peaks.strucLoop.${filt}_SA2mut.SA2.ghist.txt SA2mut.RAD21peaks.strucLoop.${filt}_SA2mut.SA1.ghist.txt" \
-s "STAG2 STAG1" -c "springgreen2 darkgoldenrod2" -x 1000 -y "0 18" -d ${FIGURESDIR}/ghist/n2tresults/strucANCH/avbw/STAG -n SA2mut.SA1vsSA2.structural_RAD21peaks.${filt}
plotHIST.sh -g "RAD21mut.RAD21peaks.strucLoop.${filt}_RAD21mut.SA2.ghist.txt RAD21mut.RAD21peaks.strucLoop.${filt}_RAD21mut.SA1.ghist.txt" \
-s "STAG2 STAG1" -c "springgreen2 darkgoldenrod2" -x 1000 -y "0 18" -d ${FIGURESDIR}/ghist/n2tresults/strucANCH/avbw/STAG -n RAD21mut.SA1vsSA2.structural_RAD21peaks.${filt}
plotHIST.sh -g "CTRL.RAD21peaks.strucLoop.${filt}_SA2mut.SA2.ghist.txt CTRL.RAD21peaks.strucLoop.${filt}_SA2mut.SA1.ghist.txt" \
-s "STAG2 STAG1" -c "springgreen2 darkgoldenrod2" -x 1000 -y "0 18" -d ${FIGURESDIR}/ghist/n2tresults/strucANCH/avbw/STAG -n CTRL.SA1vsSA2.structural_RAD21peaks.${filt}
done



################################################################################
#            Motif analyist all structural diffloop RAD21 pos combined         #
#                     ignoring paired/unpaired status                          #
################################################################################

###############################################run motif analysis for EP-anchors with cohesin

for filt in ${FiltsFC1}; do
for MUT in ${MUTS}; do
type="structuralLoops"
 cat >"${TMPDIR}/motifs.${type}${MUT}${filt}.${_DATE}.sh" <<EOF
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
#find and filter motifs
cd ${STRUCLOOPDIR}/
findMotifsGenome.pl ${MUT}vsCTRL.loopAnchors.red.${filt}.structural.RAD21.centred.bed hg38 ${motifanchdir}/${type}_RAD21peaks.${filt}_${MUT} -size given -len 7,8,9,10,11,12,13,14 -p 6 -h
compareMotifs.pl ${motifanchdir}/${type}_RAD21peaks.${filt}_${MUT}/homerMotifs.all.motifs ${motifanchdir}/${type}_RAD21peaks.${filt}_${MUT}/final -reduceThresh .75 -matchThresh .6 -pvalue 1e-12 -info 1.5 -cpu 6
EOF
chmod 750 "${TMPDIR}/motifs.${type}${MUT}${filt}.${_DATE}.sh"
	echo "motifs for ${type} ${MUT} ${filt}"
	screen -dm -S ${type}${MUT}${filt} bash -c "bash ${TMPDIR}/motifs.${type}${MUT}${filt}.${_DATE}.sh"
done
done

