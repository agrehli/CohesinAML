#!/bin/bash
# by AF, APR 2022


###############################################################################
###############################################################################
##                                                                           ##
##       Analysis of MED12 ChIPseq in Cohesin KDs in CD34+ HSPCs             ##
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
MED12SCALEDBW="${WORKDIR}/MED12/scaledBigWigs"
MERGEDBW="${WORKDIR}/mergedBigWigs"
MOTIFDIR="${WORKDIR}/motifs"
DTmatrixdir="${WORKDIR}/deeptoolsMatrix"

mkdir $MED12SCALEDBW



                         #####################################
                         #        peak finding for QC        #
                         #####################################

###Declare MED12 ChIP Samples
CTRL_MED12='ChIP_CD34_14_3_siCtrl_MED12 ChIP_CD34_14_4_Mock_MED12 ChIP_CD34_18_4_siCtrl_MED12 ChIP_CD34_19_2_siCtrl_MED12 ChIP_CD34_21_4_siCtrl_MED12 ChIP_CD34_24_2_siCtrl_MED12 ChIP_CD34_28_6_siCtrl_MED12 ChIP_CD34_21_4_siCtrl_MED12_Rep1 ChIP_CD34_22_3_siCtrl_MED12 '
SA1KD_MED12='ChIP_CD34_14_1_SA1_2259_4094_MED12 ChIP_CD34_17_1_SA1_2259_4094_MED12 ChIP_CD34_21_2_SA1_2259_4094_MED12 ChIP_CD34_27_3_SA1_2259_4094_MED12 ChIP_CD34_28_4_SA1_2259_4094_MED12 '
SA2KD_MED12='ChIP_CD34_14_2_SA2_529_1252_MED12 ChIP_CD34_17_2_SA2_529_1252_MED12 ChIP_CD34_21_3_SA2_529_1252_MED12 ChIP_CD34_22_2_SA2_529_1252_MED12 ChIP_CD34_28_5_SA2_529_1252_MED12 '
RAD21KD_MED12='ChIP_CD34_22_1_RAD21_467_2031_MED12 ChIP_CD34_27_1_RAD21_467_2031_MED12 ChIP_CD34_28_1_RAD21_467_2031_MED12 '
MED12TAGDIRS=$CTRL_MED12$SA1KD_MED12$SA2KD_MED12$RAD21KD_MED12

# MED12 PEAKS
for SAMPLE in ${MED12TAGDIRS};do
findPeaks ${SAMPLE} -i ${INPUTDIR}/Input_CD34_merged -style factor -o auto
done

#save peak stats to file
echo -e "Sample_ID\tTotalPeaks\tIPefficiency\tTotalTags" > $PEAKDIR/MED12peaks.default.txt
for SAMPLE in ${MED12TAGDIRS}; do
    sample_name="${SAMPLE}"
    TotalPeaks=$(grep -w "# total peaks =" $sample_name/peaks.txt | awk '{print $5}')
    IPefficiency=$(grep -w "# Approximate IP efficiency =" $sample_name/peaks.txt | awk '{print $6}')
    TotalTags=$(grep -w "# Total tags =" $sample_name/peaks.txt | awk '{print $5}')
    echo -e "${sample_name}\t${TotalPeaks}\t${IPefficiency}\t${TotalTags}" >> $PEAKDIR/MED12peaks.default.txt
done



                         ##################################################
                         #       merged tagDirectories by condition       #
                         ##################################################


# all replicates combined by group
cd ${TAGDIR}
makeTagDirectory ChIP_merged_CD34_CTRL_MED12 -d ${CTRL_MED12}
makeTagDirectory ChIP_merged_CD34_SA1KD_MED12 -d ${SA1KD_MED12}
makeTagDirectory ChIP_merged_CD34_SA2KD_MED12 -d ${SA2KD_MED12}
makeTagDirectory ChIP_merged_CD34_RAD21KD_MED12 -d ${RAD21KD_MED12}
#  bigWigs from merged tagdirs (meged bigwigs)
cd ${TAGDIR}
declare -a conditions=("CTRL_MED12" "SA1KD_MED12" "SA1KD_MED12" "SA2KD_MED12" "RAD21KD_MED12")
for SAMPLE in "${conditions[@]}";do
makeUCSCfile ChIP_merged_CD34_${SAMPLE} -o ${MERGEDBW}/merged.ChIP_${SAMPLE}.bigWig -fragLength 200 -bigWig ${CHROMSIZES_HG38} -fsize 1e20
done


                     ############################################
                     #      peak finding in merged tagDirs      #
                     ############################################
# peak finding
cd ${TAGDIR}
for SAMPLE in "${conditions[@]}";
do
findPeaks ChIP_merged_CD34_${SAMPLE} -i ${INPUTDIR}/Input_CD34_merged -style factor -tbp 1 -fdr 0.000001 -o ${PEAKDIR}/${SAMPLE}.peaks.txt
done
# peak filtering
cd ${PEAKDIR}
for PEAKSET in "${PEAKSETS[@]}";do
	pos2bed.pl ${PEAKSET}.peaks.txt > ${TMPDIR}/tmp.1.bed
	$BEDTOOLS intersect -a ${TMPDIR}/tmp.1.bed -b $BLACKLIST_HG38 -v > ${TMPDIR}/tmp.2.bed
	bed2pos.pl ${TMPDIR}/tmp.2.bed > ${TMPDIR}/tmp.1.txt
	filter4Mappability.sh -p ${TMPDIR}/tmp.1.txt -g hg38 -f 0.8 -s 75
	pos2bed.pl ${TMPDIR}/tmp.1.mapScoreFiltered.txt > ${PEAKSET}.filtered.peaks.bed
	bed2pos.pl ${PEAKSET}.filtered.peaks.bed > ${PEAKSET}.filtered.peaks.txt
done
#look at peak numbers
cd ${PEAKDIR}
for PEAKSET in "${PEAKSETS[@]}";do
wc -l ${PEAKSET}.filtered.peaks.txt
done

                     ##################################################
                     #             merged MED12 peak sets             #
                     ##################################################
###combine MED12 peaks from all conditions
mergePeaks ${PEAKDIR}/CTRL_MED12.filtered.peaks.txt ${PEAKDIR}/SA1KD_MED12.filtered.peaks.txt \
${PEAKDIR}/SA2KD_MED12.filtered.peaks.txt ${PEAKDIR}/RAD21KD_MED12.filtered.peaks.txt -code > ${PEAKDIR}/CD34_MED12.filtered.peaks.txt
wc -l ${PEAKDIR}/CD34_MED12.filtered.peaks.txt

                         #####################################
                         #        scaling of bigWigs         #
                         #####################################
#to adjust tracks for technical variance HSPC data is scaled based on the top peaks found in all samples from the merged peakset
# define top ~1500 peaks for scaling
cd ${PEAKDIR}
for PEAKSET in "${PEAKSETS[@]}";do
awk '$2 ~ /^chr/&& $2 ~ !/^chrM/' ${PEAKDIR}/${PEAKSET}.peaks.txt > ${PEAKDIR}/${PEAKSET}.peaks.refChr.txt
head -n 3000 ${PEAKDIR}/${PEAKSET}.peaks.refChr.txt | cut -f 1-6 > ${PEAKDIR}/${PEAKSET}.peaks.refChr.top.txt
done


###Create a merged peak set of the top peaks
mergePeaks ${PEAKDIR}/CTRL_MED12.peaks.refChr.top.txt ${PEAKDIR}/SA1KD_MED12.peaks.refChr.top.txt \
${PEAKDIR}/SA2KD_MED12.peaks.refChr.top.txt ${PEAKDIR}/RAD21KD_MED12.peaks.refChr.top.txt -code >${TMPDIR}/merge1.txt
awk -v "key=1111" '$7 == key' ${TMPDIR}/merge1.txt >${PEAKDIR}/all_commonCD34.MED12.peaks.refChr.top.txt

# calculate median normalization factor for MED12
wc -l ${PEAKDIR}/all_commonCD34.MED12.peaks.refChr.top.txt #2399
# calculate median normalization factor
cd ${TAGDIR}
#number of samples: 22
annotatePeaks.pl ${PEAKDIR}/all_commonCD34.MED12.peaks.refChr.top.txt hg38 -size given -d ${MED12TAGDIRS} -noann -nogene -cpu 12 > "${PEAKDIR}/all_commonCD34.MED12.peaks.refChr.top.allAnn.txt"
awk 'BEGIN{OFS="\t"} !/chrM/ {sum=0; for(i=8; i<=NF; i++) sum += $i; ave=sum/22; for(j=1;j<=NF;j++){printf "%s\t", $j} print ave }' "${PEAKDIR}/all_commonCD34.MED12.peaks.refChr.top.allAnn.txt" | tail -n +2 | cut -f 8-31 > ${TMPDIR}/test2.txt
awk 'BEGIN{OFS="\t"} { for(i=1;i<NF;i++){printf "%s\t", $i/$23}  {print ""} }' ${TMPDIR}/test2.txt >  ${TMPDIR}/test3.txt
allMEDIANS=()
for ((i=1;i<=22;i++));do
allMEDIANS[$i]=$(cut -f ${i} ${TMPDIR}/test3.txt | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }')
done
echo ${allMEDIANS[@]}
rm ${TMPDIR}/test2.txt ${TMPDIR}/test3.txt
#

# scale bigWigs based on the top peaks
cd ${BWDIR}
_DATE=$(date +%s)
COUNT=0
for SAMPLE in ${MED12TAGDIRS}; do
    COUNT2=$((COUNT+1))
    FACTOR=$(awk -v "key=${allMEDIANS[$COUNT2]}" 'BEGIN {print 1 / key}')
    echo $FACTOR
	cat >"${TMPDIR}/track.${SAMPLE}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.9/bin:${PATH}
export PATH
cd ${TMPDIR}
scaleBedGraph.pl ${BWDIR}/${SAMPLE}.bigwig -M ${FACTOR} -bigWig hg38 -o ${MED12SCALEDBW}/${SAMPLE}.scaled
EOF
	COUNT=$((COUNT+=1))
	chmod 750 "${TMPDIR}/track.${SAMPLE}.${_DATE}.sh"
	echo "converting track for ${SAMPLE}"
	screen -dm -S track${SAMPLE} bash -c "bash ${TMPDIR}/track.${SAMPLE}.${_DATE}.sh"
done


# generate average BigWigs
#scaled BW file names by condition
normCTRLMED12bw=()
for i in ${CTRL_MED12}; do
    bGname=${i}.scaled.bigWig
    normCTRLMED12bw+=("$bGname")
done
echo "${normCTRLMED12bw[@]}"

normSA1KDMED12bw=()
for i in ${SA1KD_MED12}; do
    bGname=${i}.scaled.bigWig
    normSA1KDMED12bw+=("$bGname")
done
echo "${normSA1KDMED12bw[@]}"

normSA2KDMED12bw=()
for i in ${SA2KD_MED12}; do
    bGname=${i}.scaled.bigWig
    normSA2KDMED12bw+=("$bGname")
done
echo "${normSA2KDMED12bw[@]}"

cd ${MED12SCALEDBW}
myAverageBigWig.pl -bw ${normCTRLMED12bw[@]} -chr ${CHROMSIZES_HG38} -o ave.CD34_CTRL_MED12.scaled.bigWig
myAverageBigWig.pl -bw ${normSA1KDMED12bw[@]} -chr ${CHROMSIZES_HG38} -o ave.CD34_SA1KD_MED12.scaled.bigWig
myAverageBigWig.pl -bw ${normSA2KDMED12bw[@]} -chr ${CHROMSIZES_HG38}  -o ave.CD34_SA2KD_MED12.scaled.bigWig

bigWigToBedGraph ave.CD34_CTRL_MED12.scaled.bigWig ave.CD34_CTRL_MED12.bedGraph
bigWigToBedGraph ave.CD34_SA1KD_MED12.scaled.bigWig ave.CD34_SA1KD_MED12.bedGraph
bigWigToBedGraph ave.CD34_SA2KD_MED12.scaled.bigWig ave.CD34_SA2KD_MED12.bedGraph







