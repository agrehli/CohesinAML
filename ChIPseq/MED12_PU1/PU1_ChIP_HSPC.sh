#!/bin/bash
## Alexander Fischer Nov 2021

###############################################################################
###############################################################################
##                                                                           ##
##       Analysis of PU1  ChIPseq  in CD34 HSPC                              ##
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

CHROMSIZES_HG38="/misc/software/viewer/IGV/IGVTools_2.3.98/genomes/GRCh38.PRI_p10.chrom.sizes"
CHAINFILE="/misc/software/ngs/genome/chainFiles/hg19ToHg38.over.chain"
BLACKLIST_HG38="/misc/data/analysis/generalStuff/annotation/GRCh38/hg38.blacklist.bed"
HG38GTF="${DIR_PKG}/genome/annotation/GRCh38.PRI_p10/gencode.v27.primary_assembly.annotation.gtf"



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
PU1SCALEDBW="${WORKDIR}/PU1/scaledBigWigs"
MERGEDBW="${WORKDIR}/mergedBigWigs"
MOTIFDIR="${WORKDIR}/motifs"
DTmatrixdir="${WORKDIR}/deeptoolsMatrix"

mkdir $PU1SCALEDBW



                         #####################################
                         #        peak finding for QC        #
                         #####################################

##Define Tagdir names
CTRL_PU1='ChIP_CD34_14_3_siCtrl_PU1 ChIP_CD34_22_3_siCtrl_PU1 ChIP_CD34_24_2_siCtrl_PU1 '
SA2KD_PU1='ChIP_CD34_14_2_SA2_529_1252_PU1 ChIP_CD34_17_2_SA2_529_1252_PU1 ChIP_CD34_20_5_SA2_529_1252_PU1 ChIP_CD34_22_2_SA2_529_1252_PU1  '
SA1KD_PU1='ChIP_CD34_17_1_SA1KD_PU1poly ChIP_CD34_21_2_SA1KD_PU1poly ChIP_CD34_27_3_SA1KD_PU1poly'
PU1TAGDIRS=$CTRL_PU1$SA2KD_PU1$SA1KD_PU1

##QC of PEAKS for individual tagdirs
cd ${TAGDIR}
for SAMPLE in ${PU1TAGDIRS};do
findPeaks ${SAMPLE} -i ${INPUTDIR}/Input_CD34_merged -style factor -o auto
done

#save peak stats to file
echo -e "Sample_ID\tTotalPeaks\tIPefficiency\tTotalTags" > $PEAKDIR/PU1peaks.default.txt
for SAMPLE in ${MED12TAGDIRS}; do
    sample_name="${SAMPLE}"
    TotalPeaks=$(grep -w "# total peaks =" $sample_name/peaks.txt | awk '{print $5}')
    IPefficiency=$(grep -w "# Approximate IP efficiency =" $sample_name/peaks.txt | awk '{print $6}')
    TotalTags=$(grep -w "# Total tags =" $sample_name/peaks.txt | awk '{print $5}')
    echo -e "${sample_name}\t${TotalPeaks}\t${IPefficiency}\t${TotalTags}" >> $PEAKDIR/PU1peaks.default.txt
done



                         #####################################
                         #       merged tagDirectories       #
                         #####################################

# all replicates combined in tagdir and merged bigwig 
cd ${TAGDIR}
makeTagDirectory ChIP_merged_CD34_CTRL_PU1 -d ${CTRL_PU1}
makeUCSCfile ChIP_merged_CD34_CTRL_PU1 -o ${MERGEDBW}/merged.ChIP_CTRL_PU1.bigWig -fragLength 200 -bigWig ${CHROMSIZES_HG38} -fsize 1e20

makeTagDirectory ChIP_merged_CD34_SA2KD_PU1 -d ${SA2KD_PU1}
makeUCSCfile ChIP_merged_CD34_SA2KD_PU1 -o ${MERGEDBW}/merged.ChIP_SA2KD_PU1.bigWig -fragLength 200 -bigWig ${CHROMSIZES_HG38} -fsize 1e20

makeTagDirectory ChIP_merged_CD34_SA1KD_PU1 -d ${SA1KD_PU1}
makeUCSCfile ChIP_merged_CD34_SA1KD_PU1 -o ${MERGEDBW}/merged.ChIP_SA1KD_PU1.bigWig -fragLength 200 -bigWig ${CHROMSIZES_HG38} -fsize 1e20


                     ############################################
                     #      peak finding in merged tagDirs      #
                     ############################################
# peak finding
findPeaks ChIP_merged_CD34_CTRL_PU1 -i ${INPUTDIR}/Input_CD34_merged -style factor -tbp 1 -fdr 0.000001 -o ${PEAKDIR}/CTRL_PU1.peaks.txt
findPeaks ChIP_merged_CD34_CTRL_PU1 -i ${INPUTDIR}/Input_CD34_merged -style factor -tbp 1 -fdr 0.000001 -ntagThreshold 10 -o ${PEAKDIR}/CTRL_PU1.stringentPeaks.txt

findPeaks ChIP_merged_CD34_SA2KD_PU1 -i ${INPUTDIR}/Input_CD34_merged -style factor -tbp 1 -fdr 0.000001 -o ${PEAKDIR}/SA2KD_PU1.peaks.txt
findPeaks ChIP_merged_CD34_SA2KD_PU1 -i ${INPUTDIR}/Input_CD34_merged -style factor -tbp 1 -fdr 0.000001 -ntagThreshold 10 -o ${PEAKDIR}/SA2KD_PU1.stringentPeaks.txt

findPeaks ChIP_merged_CD34_SA1KD_PU1 -i ${INPUTDIR}/Input_CD34_merged -style factor -tbp 1 -fdr 0.000001 -o ${PEAKDIR}/SA1KD_PU1.peaks.txt
findPeaks ChIP_merged_CD34_SA1KD_PU1 -i ${INPUTDIR}/Input_CD34_merged -style factor -tbp 1 -fdr 0.000001 -ntagThreshold 10 -o ${PEAKDIR}/SA1KD_PU1.stringentPeaks.txt

# peak filtering

declare -a PEAKSETS=("CTRL_PU1" "SA2KD_PU1"  "SA1KD_PU1")

cd ${PEAKDIR}
for PEAKSET in "${PEAKSETS[@]}";do
	pos2bed.pl ${PEAKSET}.peaks.txt > ${TMPDIR}/tmp.1.bed
	$BEDTOOLS intersect -a ${TMPDIR}/tmp.1.bed -b $BLACKLIST_HG38 -v > ${TMPDIR}/tmp.2.bed
	bed2pos.pl ${TMPDIR}/tmp.2.bed > ${TMPDIR}/tmp.1.txt
	filter4Mappability.sh -p ${TMPDIR}/tmp.1.txt -g hg38 -f 0.8 -s 75
	pos2bed.pl ${TMPDIR}/tmp.1.mapScoreFiltered.txt > ${PEAKSET}.filtered.peaks.bed
	bed2pos.pl ${PEAKSET}.filtered.peaks.bed > ${PEAKSET}.filtered.peaks.txt
done

# peak filtering stringent peaks
cd ${PEAKDIR}
for PEAKSET in "${PEAKSETS[@]}";do
	pos2bed.pl ${PEAKSET}.stringentPeaks.txt > ${TMPDIR}/tmp.1.bed
	$BEDTOOLS intersect -a ${TMPDIR}/tmp.1.bed -b $BLACKLIST_HG38 -v > ${TMPDIR}/tmp.2.bed
	bed2pos.pl ${TMPDIR}/tmp.2.bed > ${TMPDIR}/tmp.1.txt
	filter4Mappability.sh -p ${TMPDIR}/tmp.1.txt -g hg38 -f 0.8 -s 75
	pos2bed.pl ${TMPDIR}/tmp.1.mapScoreFiltered.txt > ${PEAKSET}.filtered.stringentPeaks.bed
	bed2pos.pl ${PEAKSET}.filtered.stringentPeaks.bed > ${PEAKSET}.filtered.stringentPeaks.txt
done

# merged peaksets of all conditions
mergePeaks ${PEAKDIR}/CTRL_PU1.filtered.peaks.txt ${PEAKDIR}/SA2KD_PU1.filtered.peaks.txt ${PEAKDIR}/SA1KD_PU1.filtered.peaks.txt \
-code > ${PEAKDIR}/CD34_PU1.filtered.peaks.txt

mergePeaks ${PEAKDIR}/CTRL_PU1.filtered.stringentPeaks.txt ${PEAKDIR}/SA2KD_PU1.filtered.stringentPeaks.txt ${PEAKDIR}/SA1KD_PU1.filtered.stringentPeaks.txt \
-code > ${PEAKDIR}/CD34_PU1.filtered.stringentPeaks.txt

wc -l ${PEAKDIR}/CD34_PU1.filtered.peaks.txt
#147194
wc -l ${PEAKDIR}/CD34_PU1.filtered.stringentPeaks.txt
#68181



                     ############################################
                     #          scaled averaged bigwigs         #
                     ############################################

# define top peaks for scaling
awk '$2 ~ /^chr/&& $2 ~ !/^chrM/' ${PEAKDIR}/CTRL_PU1.peaks.txt > ${PEAKDIR}/CTRL_PU1.peaks.refChr.txt
awk '$2 ~ /^chr/&& $2 ~ !/^chrM/' ${PEAKDIR}/SA1KD_PU1.peaks.txt > ${PEAKDIR}/SA1KD_PU1.peaks.refChr.txt
awk '$2 ~ /^chr/&& $2 ~ !/^chrM/' ${PEAKDIR}/SA2KD_PU1.peaks.txt > ${PEAKDIR}/SA2KD_PU1.peaks.refChr.txt

head -n 3000 ${PEAKDIR}/CTRL_PU1.peaks.refChr.txt | cut -f 1-6 > ${PEAKDIR}/tmpChIP_merged_CD34_CTRL_PU1.peaks.refChr.top.txt
head -n 3000 ${PEAKDIR}/SA1KD_PU1.peaks.refChr.txt | cut -f 1-6 > ${PEAKDIR}/tmpChIP_merged_CD34_SA1KD_PU1.peaks.refChr.top.txt
head -n 3000 ${PEAKDIR}/SA2KD_PU1.peaks.refChr.txt | cut -f 1-6 > ${PEAKDIR}/tmpChIP_merged_CD34_SA2KD_PU1.peaks.refChr.top.txt


mergePeaks ${PEAKDIR}/tmpChIP_merged_CD34_CTRL_PU1.peaks.refChr.top.txt ${PEAKDIR}/tmpChIP_merged_CD34_SA1KD_PU1.peaks.refChr.top.txt \
${PEAKDIR}/tmpChIP_merged_CD34_SA2KD_PU1.peaks.refChr.top.txt \
-code >${TMPDIR}/merge1.txt
awk -v "key=111" '$7 == key' ${TMPDIR}/merge1.txt >${PEAKDIR}/all_commonCD34.PU1.peaks.refChr.top2250.txt
rm ${PEAKDIR}/tmp*

# calculate median normalization factor
cd ${TAGDIR}
annotatePeaks.pl ${PEAKDIR}/all_commonCD34.PU1.peaks.refChr.top2250.txt hg38 -size given -d ${PU1TAGDIRS} -noann -nogene -cpu 12 > "${PEAKDIR}/all_commonCD34.PU1.peaks.refChr.top.allAnn.txt"
awk 'BEGIN{OFS="\t"} !/chrM/ {sum=0; for(i=8; i<=NF; i++) sum += $i; ave=sum/10; for(j=1;j<=NF;j++){printf "%s\t", $j} print ave }' "${PEAKDIR}/all_commonCD34.PU1.peaks.refChr.top.allAnn.txt" | tail -n +2 | cut -f 8-18 > ${TMPDIR}/test2.txt
awk 'BEGIN{OFS="\t"} { for(i=1;i<NF;i++){printf "%s\t", $i/$11}  {print ""} }' ${TMPDIR}/test2.txt >  ${TMPDIR}/test3.txt
allMEDIANS=()
for ((i=1;i<=10;i++));do
allMEDIANS[$i]=$(cut -f ${i} ${TMPDIR}/test3.txt | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }')
done
echo ${allMEDIANS[@]}
rm ${TMPDIR}/test2.txt ${TMPDIR}/test3.txt


# scale bigWigs based on the top peaks
cd ${BWDIR}
_DATE=$(date +%s)
COUNT=0
for SAMPLE in ${PU1TAGDIRS}; do
    COUNT2=$((COUNT+1))
    FACTOR=$(awk -v "key=${allMEDIANS[$COUNT2]}" 'BEGIN {print 1 / key}')
    echo $FACTOR
	cat >"${TMPDIR}/track.${SAMPLE}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.9/bin:${PATH}
export PATH
cd ${TMPDIR}
scaleBedGraph.pl ${BWDIR}/${SAMPLE}.bigwig -M ${FACTOR} -bigWig hg38 -o ${PU1SCALEDBW}/${SAMPLE}.scaled
EOF
	COUNT=$((COUNT+=1))
	chmod 750 "${TMPDIR}/track.${SAMPLE}.${_DATE}.sh"
	echo "converting track for ${SAMPLE}"
	screen -dm -S track${SAMPLE} bash -c "bash ${TMPDIR}/track.${SAMPLE}.${_DATE}.sh"
done

# generate average BigWigs
#scaled BW file names by condition
normCTRLPU1bw=()
for i in ${CTRL_PU1}; do
    bGname=${i}.scaled.bigWig
    normCTRLPU1bw+=("$bGname")
done
echo "${normCTRLPU1bw[@]}"

normSA1KDPU1bw=()
for i in ${SA1KD_PU1}; do
    bGname=${i}.scaled.bigWig
    normSA1KDPU1bw+=("$bGname")
done
echo "${normSA1KDPU1bw[@]}"

normSA2KDPU1bw=()
for i in ${SA2KD_PU1}; do
    bGname=${i}.scaled.bigWig
    normSA2KDPU1bw+=("$bGname")
done
echo "${normSA2KDPU1bw[@]}"

cd ${PU1SCALEDBW}
myAverageBigWig.pl -bw ${normCTRLPU1bw[@]} -chr ${CHROMSIZES_HG38} -o ave.CD34_CTRL_PU1.scaled.bigWig
myAverageBigWig.pl -bw ${normSA1KDPU1bw[@]} -chr ${CHROMSIZES_HG38} -o ave.CD34_SA1KD_PU1.scaled.bigWig
myAverageBigWig.pl -bw ${normSA2KDPU1bw[@]} -chr ${CHROMSIZES_HG38}  -o ave.CD34_SA2KD_PU1.scaled.bigWig

bigWigToBedGraph ave.CD34_CTRL_PU1.scaled.bigWig ave.CD34_CTRL_PU1.bedGraph
bigWigToBedGraph ave.CD34_SA1KD_PU1.scaled.bigWig ave.CD34_SA1KD_PU1.bedGraph
bigWigToBedGraph ave.CD34_SA2KD_PU1.scaled.bigWig ave.CD34_SA2KD_PU1.bedGraph



            ############################################################
            #          PU1. RAD21 intersections	                       #
            ############################################################

pos2bed.pl ${PEAKDIR}/CD34_PU1.filtered.stringentPeaks.txt > ${PEAKDIR}/CD34_PU1.filtered.stringentPeaks.bed
pos2bed.pl ${PEAKDIR}/CD34_PU1.filtered.peaks.txt > ${PEAKDIR}/CD34_PU1.filtered.peaks.bed
#overlaps with RAD21 diffpeaks
KDs="SA2KD SA1KD"
peaktypes="stringentPeaks peaks"
for pt in $peaktypes;do
for KD in $KDs;do
$BEDTOOLS intersect -b ${DIFFPEAKS}/${KD}vsCTRL.RAD21.Peaks_edgeR.model.2foldup.bed -a ${PEAKDIR}/CD34_PU1.filtered.${pt}.bed -u > ${DIFFPEAKS}/PU1.${pt}.at.increasedRAD21.${KD}vsCTRL.bed
$BEDTOOLS intersect -b ${DIFFPEAKS}/${KD}vsCTRL.RAD21.Peaks_edgeR.model.2folddown.bed -a ${PEAKDIR}/CD34_PU1.filtered.${pt}.bed -u > ${DIFFPEAKS}/PU1.${pt}.at.decreasedRAD21.${KD}vsCTRL.bed
done done
for KD in $KDs;do
for pt in $peaktypes;do
wc -l ${DIFFPEAKS}/PU1.${pt}.at.increasedRAD21.${KD}vsCTRL.bed
wc -l ${DIFFPEAKS}/PU1.${pt}.at.decreasedRAD21.${KD}vsCTRL.bed
done done





            ############################################################
            #          PU.1  peak annotation                           #
            ############################################################

# annotation to merged peak set for differnetial analyis etc.

cd ${TAGDIR}
annotatePeaks.pl "${PEAKDIR}/CD34_PU1.filtered.stringentPeaks.txt" hg38 -size 250 -d ${PU1TAGDIRS} -raw -cpu 12 > ${WORKDIR}/CD34_PU1.stringentPeaks.ann.txt
annotatePeaks.pl "${PEAKDIR}/CD34_PU1.filtered.peaks.txt" hg38 -size 250 -d ${PU1TAGDIRS} -raw -cpu 12 > ${WORKDIR}/CD34_PU1.peaks.ann.txt

tail -n +2 ${WORKDIR}/CD34_PU1.stringentPeaks.ann.txt | cut -f1,20-29 > ${TMPDIR}/tmp.1.txt
echo $'ID\t14_CTRL\t22_CTRL\t24_CTRL\t14_SA2KD\t17_SA2KD\t20_SA2KD\t22_SA2KD\t17_SA1KD\t21_SA1KD\t27_SA1KD' | cat - ${TMPDIR}/tmp.1.txt > ${TMPDIR}/CD34_PU1.stringentPeaks.ann.Rinput.txt





