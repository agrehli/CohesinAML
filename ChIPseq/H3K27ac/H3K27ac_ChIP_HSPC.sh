#!/bin/bash
# by Alexander Fischer,DEC 2021

###############################################################################
###############################################################################
##                                                                           ##
##     Analysis of H3k72ac ChIPseq data in CD34+ HSPCs with Cohesin KDs      ##
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
# general files
CHROMSIZES_HG38="/misc/software/viewer/IGV/IGVTools_2.3.98/genomes/GRCh38.PRI_p10.chrom.sizes"
BLACKLIST_HG38="/misc/data/analysis/generalStuff/annotation/GRCh38/hg38.blacklist.bed"

# general directories
TMPDIR="/loctmp"
## directories with data created by mapChIP.sh pipeline
RAWDIR="${DIR_DATA}/rawData/chromatin/ChIP/CD34/"
FASTQCDIR="${DIR_DATA}/rawData/chromatin/ChIP/CD34/FastQC"
TAGDIR="${DIR_DATA}/processedData/tagDir/chromatin/GRCh38/ChIP/CD34/siRNA_KD"
INPUTDIR="${DIR_DATA}/processedData/tagDir/DNA/GRCh38/Input/CD34"
BWDIR="${DIR_DATA}/processedData/bigWig/chromatin/GRCh38/ChIP/CD34/siRNA_KD"

PROJECTDIR="${DIR_DATA}/analysis/project_cohesin"
WORKDIR_ENH="${PROJECTDIR}/CD34/ChIP_KD_analysis/H3K27ac"
PEAKDIR="${WORKDIR_ENH}/peaks"
SEPEAKDIR="${WORKDIR_ENH}/SEpeaks"
DIFFDIR="${WORKDIR_ENH}/diffPeaks"
FIGURESDIR="${WORKDIR_ENH}/figures"
SCALEDBW="${WORKDIR_ENH}/scaledBigWigs"
MERGEDBW="${WORKDIR_ENH}/mergedBigWigs"

mkdir -p ${WORKDIR_ENH}
mkdir -p ${PEAKDIR}
mkdir -p ${SEPEAKDIR}
mkdir -p ${DIFFDIR}
mkdir -p ${FIGURESDIR}
mkdir -p ${SCALEDBW}
mkdir -p ${MOTIFDIR}
mkdir -p ${MERGEDBW}

                         #####################################
                         #        peak finding for QC        #
                         #####################################
#------------------------------------------------------------------------------------                         

#Declare samples
CTRL_H3K27ac='ChIP_CD34_14_3_siCtrl_H3K27ac ChIP_CD34_14_4_Mock_H3K27ac ChIP_CD34_17_3_siCtrl_H3K27ac ChIP_CD34_18_4_siCtrl_H3K27ac ChIP_CD34_19_2_siCtrl_H3K27ac ChIP_CD34_21_4_siCtrl_H3K27ac ChIP_CD34_22_3_siCtrl_H3K27ac ChIP_CD34_27_4_siCtrl_H3K27ac ChIP_CD34_28_6_siCtrl_H3K27ac ChIP_CD34_24_2_siCtrl_H3K27ac '
SA1KD_H3K27ac='ChIP_CD34_14_1_SA1_2259_4094_H3K27ac ChIP_CD34_17_1_SA1_2259_4094_H3K27ac ChIP_CD34_21_2_SA1_2259_4094_H3K27ac ChIP_CD34_27_3_SA1_2259_4094_H3K27ac ChIP_CD34_28_4_SA1_2259_4094_H3K27ac '
SA2KD_H3K27ac='ChIP_CD34_14_2_SA2_529_1252_H3K27ac ChIP_CD34_17_2_SA2_529_1252_H3K27ac ChIP_CD34_20_5_SA2_529_1252_H3K27ac ChIP_CD34_21_3_SA2_529_1252_H3K27ac ChIP_CD34_22_2_SA2_529_1252_H3K27ac ChIP_CD34_28_5_SA2_529_1252_H3K27ac '
RAD21KD_H3K27ac='ChIP_CD34_18_1_RAD21_467_2031_H3K27ac ChIP_CD34_22_1_RAD21_467_2031_H3K27ac ChIP_CD34_27_1_RAD21_467_2031_H3K27ac ChIP_CD34_28_1_RAD21_467_2031_H3K27ac '
H3K27acTAGDIRS=$CTRL_H3K27ac$SA1KD_H3K27ac$SA2KD_H3K27ac$RAD21KD_H3K27ac
#------------------------------------------------------------------------------------                         

# find peaks in histone style (HOMER) wiht otherwise default settings
cd ${TAGDIR}
for SAMPLE in ${H3K27acTAGDIRS}
do
findPeaks ${SAMPLE} -i ${INPUTDIR}/Input_CD34_merged -style histone -o auto
done
#------------------------------------------------------------------------------------                         

# show values
cd ${TAGDIR}
for SAMPLE in ${H3K27acTAGDIRS}
do
paste <(echo -e "# ${SAMPLE}") <(grep -w "# total peaks =" ${SAMPLE}/regions.txt) <(grep -w "# Approximate IP efficiency" ${SAMPLE}/regions.txt) <(grep -w "# Total tags =" ${SAMPLE}/regions.txt)
done
#save to file
echo -e "Sample_ID\tTotalPeaks\tIPefficiency\tTotalTags" > $PEAKDIR/H3k27peaks.styleHistone.txt
for SAMPLE in ${H3K27acTAGDIRS}; do
    sample_name="${SAMPLE}"
    TotalPeaks=$(grep -w "# total peaks =" $sample_name/regions.txt | awk '{print $5}')
    IPefficiency=$(grep -w "# Approximate IP efficiency =" $sample_name/regions.txt | awk '{print $6}')
    TotalTags=$(grep -w "# Total tags =" $sample_name/regions.txt | awk '{print $5}')
    echo -e "${sample_name}\t${TotalPeaks}\t${IPefficiency}\t${TotalTags}" >> $PEAKDIR/H3k27peaks.styleHistone.txt
done
#------------------------------------------------------------------------------------                         

                         #####################################
                         #       merged tagDirectories       #
                         #####################################

#------------------------------------------------------------------------------------                         
# all replicates combined
cd ${TAGDIR}
makeTagDirectory ChIP_merged_CD34_CTRL_H3K27ac -d ${CTRL_H3K27ac}
makeTagDirectory ChIP_merged_CD34_RAD21KD_H3K27ac -d ${RAD21KD_H3K27ac}
makeTagDirectory ChIP_merged_CD34_SA1KD_H3K27ac -d ${SA1KD_H3K27ac}
makeTagDirectory ChIP_merged_CD34_SA2KD_H3K27ac -d ${SA2KD_H3K27ac}
#------------------------------------------------------------------------------------                         
# merged bigWigs
cd ${TAGDIR}
declare -a conditions=("SA1KD_H3K27ac" "SA2KD_H3K27ac" "RAD21KD_H3K27ac" "CTRL_H3K27ac")
for SAMPLE in "${conditions[@]}";do
makeUCSCfile ChIP_merged_CD34_${SAMPLE} -o ${MERGEDBW}/merged.ChIP_${SAMPLE}.bigWig -fragLength 200 -bigWig ${CHROMSIZES_HG38} -fsize 1e20
bigWigToBedGraph ${MERGEDBW}/merged.ChIP_${SAMPLE}.bigWig ${MERGEDBW}/merged.CD34_${SAMPLE}.bedGraph
done
#------------------------------------------------------------------------------------                         

                     ############################################
                     #      peak finding in merged tagDirs      #
                     ############################################
#------------------------------------------------------------------------------------                         
#peak finding with defined size parameters in replicate merged data
cd ${TAGDIR}
declare -a conditions=("SA1KD_H3K27ac" "SA2KD_H3K27ac" "RAD21KD_H3K27ac" "CTRL_H3K27ac")
for SAMPLE in "${conditions[@]}";do
findPeaks ChIP_merged_CD34_${SAMPLE} -i ${INPUTDIR}/Input_CD34_merged -region -size 250 -L 0 -F 5 -minDist 350 -fdr 0.00001 -ntagThreshold 10 -o ${PEAKDIR}/ChIP_merged_CD34_${SAMPLE}.peaks.txt
done
#------------------------------------------------------------------------------------                         
# peak filtering for blacklisted region and low mappability
declare -a PEAKSETS=("SA1KD_H3K27ac" "SA2KD_H3K27ac" "RAD21KD_H3K27ac" "CTRL_H3K27ac")
cd ${PEAKDIR}
for PEAKSET in "${PEAKSETS[@]}";do
	pos2bed.pl ChIP_merged_CD34_${PEAKSET}.peaks.txt > ${TMPDIR}/tmp.1.bed
	$BEDTOOLS intersect -a ${TMPDIR}/tmp.1.bed -b $BLACKLIST_HG38 -v > ${TMPDIR}/tmp.2.bed
	bed2pos.pl ${TMPDIR}/tmp.2.bed > ${TMPDIR}/tmp.1.txt
	filter4Mappability.sh -p ${TMPDIR}/tmp.1.txt -g hg38 -f 0.8 -s 75
	pos2bed.pl ${TMPDIR}/tmp.1.mapScoreFiltered.txt > ${PEAKSET}.filtered.peaks.bed
	bed2pos.pl ${PEAKSET}.filtered.peaks.bed > ${PEAKSET}.filtered.peaks.txt
done
#------------------------------------------------------------------------------------                         


                         #####################################
                         #        scaling of bigWigs         #
                         #####################################
#to adjust tracks for technical variance HSPC bigwigs are scaled based on the top peaks found in all samples from the merged peakset

# define top ~1500-3000 peaks for scaling
#------------------------------------------------------------------------------------                         
###get ref chromosomes only
awk '$2 ~ /^chr/&& $2 ~ !/^chrM/' ${PEAKDIR}/ChIP_merged_CD34_SA1KD_H3K27ac.peaks.txt > ${PEAKDIR}/ChIP_merged_CD34_SA1KD_H3K27ac.peaks.refChr.txt
awk '$2 ~ /^chr/&& $2 ~ !/^chrM/' ${PEAKDIR}/ChIP_merged_CD34_SA2KD_H3K27ac.peaks.txt > ${PEAKDIR}/ChIP_merged_CD34_SA2KD_H3K27ac.peaks.refChr.txt
awk '$2 ~ /^chr/&& $2 ~ !/^chrM/' ${PEAKDIR}/ChIP_merged_CD34_RAD21KD_H3K27ac.peaks.txt > ${PEAKDIR}/ChIP_merged_CD34_RAD21KD_H3K27ac.peaks.refChr.txt
awk '$2 ~ /^chr/&& $2 ~ !/^chrM/' ${PEAKDIR}/ChIP_merged_CD34_CTRL_H3K27ac.peaks.txt > ${PEAKDIR}/ChIP_merged_CD34_CTRL_H3K27ac.peaks.refChr.txt
###get top only
head -n 3000 ${PEAKDIR}/ChIP_merged_CD34_SA1KD_H3K27ac.peaks.txt | cut -f 1-6 > ${PEAKDIR}/tmpChIP_merged_CD34_SA1KD_H3K27ac.peaks.refChr.top.txt
head -n 3000 ${PEAKDIR}/ChIP_merged_CD34_SA2KD_H3K27ac.peaks.txt | cut -f 1-6 > ${PEAKDIR}/tmpChIP_merged_CD34_SA2KD_H3K27ac.peaks.refChr.top.txt
head -n 3000 ${PEAKDIR}/ChIP_merged_CD34_RAD21KD_H3K27ac.peaks.txt | cut -f 1-6 > ${PEAKDIR}/tmpChIP_merged_CD34_RAD21KD_H3K27ac.peaks.refChr.top.txt
head -n 3000 ${PEAKDIR}/ChIP_merged_CD34_CTRL_H3K27ac.peaks.txt | cut -f 1-6 > ${PEAKDIR}/tmpChIP_merged_CD34_CTRL_H3K27ac.peaks.refChr.top.txt
#merge the top peaks from all sets and get only the ones present in all conditions
mergePeaks ${PEAKDIR}/tmpChIP_merged_CD34_SA1KD_H3K27ac.peaks.refChr.top.txt ${PEAKDIR}/tmpChIP_merged_CD34_SA2KD_H3K27ac.peaks.refChr.top.txt \
${PEAKDIR}/tmpChIP_merged_CD34_RAD21KD_H3K27ac.peaks.refChr.top.txt ${PEAKDIR}/tmpChIP_merged_CD34_CTRL_H3K27ac.peaks.refChr.top.txt \
-code >${TMPDIR}/merge1.txt
awk -v "key=1111" '$7 == key' ${TMPDIR}/merge1.txt >${PEAKDIR}/all_commonCD34.H3K27ac.peaks.refChr.top1600.txt
rm ${PEAKDIR}/tmp*

# calculate median normalization factor
#------------------------------------------------------------------------------------  
cd ${TAGDIR}
annotatePeaks.pl ${PEAKDIR}/all_commonCD34.H3K27ac.peaks.refChr.top1600.txt hg38 -size given -d ${H3K27acTAGDIRS} -noann -nogene -cpu 12 > "${PEAKDIR}/all_commonCD34.H3K27ac.peaks.refChr.top1600.allAnn.txt"
awk 'BEGIN{OFS="\t"} !/chrM/ {sum=0; for(i=8; i<=NF; i++) sum += $i; ave=sum/25; for(j=1;j<=NF;j++){printf "%s\t", $j} print ave }' "${PEAKDIR}/all_commonCD34.H3K27ac.peaks.refChr.top1600.allAnn.txt" | tail -n +2 | cut -f 8-33 > ${TMPDIR}/test2.txt
awk 'BEGIN{OFS="\t"} { for(i=1;i<NF;i++){printf "%s\t", $i/$26}  {print ""} }' ${TMPDIR}/test2.txt >  ${TMPDIR}/test3.txt
allMEDIANS=()
for ((i=1;i<=25;i++));do
allMEDIANS[$i]=$(cut -f ${i} ${TMPDIR}/test3.txt | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }')
done
echo ${allMEDIANS[@]}
rm ${TMPDIR}/test2.txt ${TMPDIR}/test3.txt

# scale bigWigs based on the top peaks
#------------------------------------------------------------------------------------  
cd ${BWDIR}
_DATE=$(date +%s)
COUNT=0
for SAMPLE in ${H3K27acTAGDIRS}; do
    COUNT2=$((COUNT+1))
    FACTOR=$(awk -v "key=${allMEDIANS[$COUNT2]}" 'BEGIN {print 1 / key}')
    echo $FACTOR
	cat >"${TMPDIR}/track.${SAMPLE}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.9/bin:${PATH}
export PATH
cd ${TMPDIR}
scaleBedGraph.pl ${BWDIR}/${SAMPLE}.bigwig -M ${FACTOR} -bigWig hg38 -o ${SCALEDBW}/${SAMPLE}.scaled
EOF
	COUNT=$((COUNT+=1))
	chmod 750 "${TMPDIR}/track.${SAMPLE}.${_DATE}.sh"
	echo "converting track for ${SAMPLE}"
	screen -dm -S track${SAMPLE} bash -c "bash ${TMPDIR}/track.${SAMPLE}.${_DATE}.sh"
done
# loop to check when screen sessions are done
#---------------------------------------------
for SAMPLE in ${H3K27acTAGDIRS}; do 
while [ true ]; do     # Endless loop.
  pid=`screen -S track${NAME} -Q echo '$PID'`  # Get a pid.
  if [[ $pid = *"session"* ]] ; then # If there is none,
    echo -e "\tFinished sample ${NAME}"
    break               # Test next one.
  else
    sleep 10           # Else wait.
  fi
	done
done
#---------------------------------------------



# generate average scaled BigWigs
#------------------------------------------------------------------------------------  
##declare bigwig names by group
CTRL_H3K27ac_bigwigs=()
for i in ${CTRL_H3K27ac}; do
    bGname=${i}.scaled.bigWig
    CTRL_H3K27ac_bigwigs+=("$bGname")
done
echo "${CTRL_H3K27ac_bigwigs[@]}"

SA1KD_H3K27ac_bigwigs=()
for i in ${SA1KD_H3K27ac}; do
    bGname=${i}.scaled.bigWig
    SA1KD_H3K27ac_bigwigs+=("$bGname")
done
echo "${SA1KD_H3K27ac_bigwigs[@]}"

SA2KD_H3K27ac_bigwigs=()
for i in ${SA2KD_H3K27ac}; do
    bGname=${i}.scaled.bigWig
    SA2KD_H3K27ac_bigwigs+=("$bGname")
done
echo "${SA2KD_H3K27ac_bigwigs[@]}"

RAD21KD_H3K27ac_bigwigs=()
for i in ${RAD21KD_H3K27ac}; do
    bGname=${i}.scaled.bigWig
    RAD21KD_H3K27ac_bigwigs+=("$bGname")
done
echo "${RAD21KD_H3K27ac_bigwigs[@]}"
#------------------------------------------------------------------------------------ 
#calculate average scaled bigwigs
cd ${SCALEDBW}
myAverageBigWig.pl -bw ${CTRL_H3K27ac_bigwigs[@]}  -chr ${CHROMSIZES_HG38} -o ${AVEBW}/ave.CD34_CTRL_H3K27ac.scaled.bigWig
myAverageBigWig.pl -bw ${SA1KD_H3K27ac_bigwigs[@]}  -chr ${CHROMSIZES_HG38} -o ${AVEBW}/ave.CD34_SA1KD_H3K27ac.scaled.bigWig
myAverageBigWig.pl -bw ${SA2KD_H3K27ac_bigwigs[@]}  -chr ${CHROMSIZES_HG38} -o ${AVEBW}/ave.CD34_SA2KD_H3K27ac.scaled.bigWig
myAverageBigWig.pl -bw ${RAD21KD_H3K27ac[@]}  -chr ${CHROMSIZES_HG38} -o ${AVEBW}/ave.CD34_RAD21KD_H3K27ac.scaled.bigWig
#------------------------------------------------------------------------------------
#also generate bedgraphs
bigWigToBedGraph ave.CD34_CTRL_H3K27ac.scaled.bigWig ave.CD34_CTRL_H3K27ac.scaled.bedGraph
bigWigToBedGraph ave.CD34_SA1KD_H3K27ac.scaled.bigWig ave.CD34_SA1KD_H3K27ac.scaled.bedGraph
bigWigToBedGraph ave.CD34_SA2KD_H3K27ac.scaled.bigWig ave.CD34_SA2KD_H3K27ac.scaled.bedGraph
bigWigToBedGraph ave.CD34_RAD21KD_H3K27ac.scaled.bigWig ave.CD34_RAD21KD_H3K27ac.scaled.bedGraph
#------------------------------------------------------------------------------------  



            ##############################################################################
            #      Calculate GC and length for cqn correction with refChr peaks          #
            ##############################################################################

#------------------------------------------------------------------------------------  
##check for each condition peakset
declare -a PEAKSETS=("SA1KD_H3K27ac" "SA2KD_H3K27ac" "RAD21KD_H3K27ac" "CTRL_H3K27ac")
for PEAKSET in "${PEAKSETS[@]}";do
	awk '!/^chrM/&&/chr/' ${PEAKDIR}/${PEAKSET}.filtered.peaks.txt > ${PEAKDIR}/${PEAKSET}.filtered.peaks.refChr.txt
	annotatePeaks.pl ${PEAKDIR}/${PEAKSET}.filtered.peaks.refChr.txt hg38 -size given -noann -nogene -CpG -cpu 12 > "${PEAKDIR}/${PEAKSET}.filtered.peaks.tmp.CpGann.txt"
	touch "${PEAKDIR}/${PEAKSET}.filtered.peaks.CpGann.txt"
	echo $'ID\tlength\tgccontent' >> "${PEAKDIR}/${PEAKSET}.filtered.peaks.CpGann.txt"
	awk -v OFS='\t' '{print $1,$4-$3,$9 ; }' <(tail -n+2 "${PEAKDIR}/${PEAKSET}.filtered.peaks.tmp.CpGann.txt") >> "${PEAKDIR}/${PEAKSET}.filtered.peaks.CpGann.txt"
	rm "${PEAKDIR}/${PEAKSET}.filtered.peaks.tmp.CpGann.txt"
done
##run for merged peakset
mergePeaks "${PEAKDIR}/SA1KD_H3K27ac.filtered.peaks.txt" "${PEAKDIR}/SA2KD_H3K27ac.filtered.peaks.txt" "${PEAKDIR}/RAD21KD_H3K27ac.filtered.peaks.txt" "${PEAKDIR}/CTRL_H3K27ac.filtered.peaks.txt"-code > ${PEAKDIR}/mergedCD34_H3K27ac.filtered.peaks.txt
	awk '!/^chrM/&&/chr/' ${PEAKDIR}/mergedCD34_H3K27ac.filtered.peaks.txt > ${PEAKDIR}/mergedCD34_H3K27ac.filtered.peaks.refChr.txt
	annotatePeaks.pl ${PEAKDIR}/mergedCD34_H3K27ac.filtered.peaks.refChr.txt hg38 -size given -noann -nogene -CpG -cpu 12 > "${PEAKDIR}/mergedCD34_H3K27ac.filtered.peaks.tmp.CpGann.txt"
	rm "${PEAKDIR}/mergedCD34_H3K27ac.filtered.peaks.CpGann.txt"
	touch "${PEAKDIR}/mergedCD34_H3K27ac.filtered.peaks.CpGann.txt"
	echo $'ID\tlength\tgccontent' >> "${PEAKDIR}/mergedCD34_H3K27ac.filtered.peaks.CpGann.txt"
	awk -v OFS='\t' '{print $1,$4-$3,$9 ; }' <(tail -n+2 "${PEAKDIR}/mergedCD34_H3K27ac.filtered.peaks.tmp.CpGann.txt") >> "${PEAKDIR}/mergedCD34_H3K27ac.filtered.peaks.CpGann.txt"
	rm "${PEAKDIR}/mergedCD34_H3K27ac.filtered.peaks.tmp.CpGann.txt"

#------------------------------------------------------------------------------------  


                      ###################################################################
                      #       Annnotate individual TAGDIRS to merged peaks set          #
                      ###################################################################

# merged peak table all conditions
#------------------------------------------------------------------------------------  
cd ${TAGDIR}
annotatePeaks.pl "${PEAKDIR}/mergedCD34_H3K27ac.filtered.peaks.txt" hg38 -size given -d ${H3K27acTAGDIRS} -noann -nogene -noadj -cpu 12 > "${PEAKDIR}/mergedCD34_H3K27ac.CD34.ann.txt"
tail -n +2 "${PEAKDIR}/mergedCD34_H3K27ac.CD34.ann.txt" | cut -f1,8-32 > ${TMPDIR}/tmp.1.txt
echo $'ID\tChIP_CD34_14_3_siCtrl_H3K27ac\tChIP_CD34_14_4_Mock_H3K27ac\tChIP_CD34_17_3_siCtrl_H3K27ac\tChIP_CD34_18_4_siCtrl_H3K27ac\tChIP_CD34_19_2_siCtrl_H3K27ac\tChIP_CD34_21_4_siCtrl_H3K27ac\tChIP_CD34_22_3_siCtrl_H3K27ac\tChIP_CD34_27_4_siCtrl_H3K27ac\tChIP_CD34_28_6_siCtrl_H3K27ac\tChIP_CD34_24_2_siCtrl_H3K27ac\tChIP_CD34_14_1_SA1_2259_4094_H3K27ac\tChIP_CD34_17_1_SA1_2259_4094_H3K27ac\tChIP_CD34_21_2_SA1_2259_4094_H3K27ac\tChIP_CD34_27_3_SA1_2259_4094_H3K27ac\tChIP_CD34_28_4_SA1_2259_4094_H3K27ac\tChIP_CD34_14_2_SA2_529_1252_H3K27ac\tChIP_CD34_17_2_SA2_529_1252_H3K27ac\tChIP_CD34_20_5_SA2_529_1252_H3K27ac\tChIP_CD34_21_3_SA2_529_1252_H3K27ac\tChIP_CD34_22_2_SA2_529_1252_H3K27ac\tChIP_CD34_28_5_SA2_529_1252_H3K27ac\tChIP_CD34_18_1_RAD21_467_2031_H3K27ac\tChIP_CD34_22_1_RAD21_467_2031_H3K27ac\tChIP_CD34_27_1_RAD21_467_2031_H3K27ac\tChIP_CD34_28_1_RAD21_467_2031_H3K27ac' | cat - ${TMPDIR}/tmp.1.txt > ${WORKDIR_ENH}/mergedCD34_H3K27ac.CD34.ann.table.txt
#------------------------------------------------------------------------------------  