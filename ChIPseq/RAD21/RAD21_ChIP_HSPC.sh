#!/bin/bash
# by AF, APR 2022


###############################################################################
###############################################################################
##                                                                           ##
##       Analysis of RAD21 ChIPseq in Cohesin KDs in CD34+ HSPCs             ##
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
TAGDIR="${DIR_DATA}/processedData/tagDir/chromatin/GRCh38/ChIP/CD34/siRNA_KD" #tag directories ChIP (HOMER)
BWDIR="${DIR_DATA}/processedData/bigWig/chromatin/GRCh38/ChIP/CD34/siRNA_KD" #bigwigs
INPUTDIR="${DIR_DATA}/processedData/tagDir/DNA/GRCh38/Input/CD34" #tag directories DNA-input (HOMER)


##  directories for analysis
PROJECTDIR="${DIR_DATA}/analysis/project_cohesin"
WORKDIR="${PROJECTDIR}/CD34/ChIP_KD_analysis/Cohesin_CTCF_MED12"
PEAKDIR="${WORKDIR}/peaks"
DIFFPEAKS="${WORKDIR}/diffPeaks"
FIGURESDIR="${WORKDIR}/figures"
RAD21SCALEDBW="${WORKDIR}/RAD21/scaledBigWigs"
MERGEDBW="${WORKDIR}/mergedBigWigs"
MOTIFDIR="${WORKDIR}/motifs"
DTmatrixdir="${WORKDIR}/deeptoolsMatrix"

#AML data for comparisons
AMLRAD21DIR="${PROJECTDIR}/Cohesin_AML/ChIP_analysis"
AMLPEAKS="${AMLRAD21DIR}/peaks"
AMLDIFFPEAKS="${AMLRAD21DIR}/diffPeaks"

# generating new directories
mkdir -p ${WORKDIR}
mkdir -p ${PEAKDIR}
mkdir -p ${DIFFPEAKS}
mkdir -p ${FIGURESDIR}/ghist
mkdir -p ${MERGEDBW}
mkdir -p ${MOTIFDIR}
mkdir -p ${DTmatrixdir}
mkdir -p ${RAD21SCALEDBW}





                         #####################################
                         #        peak finding for QC        #
                         #####################################

###Declare RAD21 ChIP Samples
CTRL_RAD21='ChIP_CD34_14_3_siCtrl_RAD21 ChIP_CD34_14_4_Mock_RAD21 ChIP_CD34_17_3_siCtrl_RAD21 ChIP_CD34_18_4_siCtrl_RAD21 ChIP_CD34_19_2_siCtrl_RAD21 ChIP_CD34_21_4_siCtrl_RAD21 ChIP_CD34_22_3_siCtrl_RAD21 ChIP_CD34_24_2_siCtrl_RAD21 ChIP_CD34_27_4_siCtrl_RAD21 ChIP_CD34_28_6_siCtrl_RAD21 '
SA1KD_RAD21='ChIP_CD34_14_1_SA1_2259_4094_RAD21 ChIP_CD34_17_1_SA1_2259_4094_RAD21 ChIP_CD34_21_2_SA1_2259_4094_RAD21 ChIP_CD34_27_3_SA1_2259_4094_RAD21 ChIP_CD34_28_4_SA1_2259_4094_RAD21 '
SA2KD_RAD21='ChIP_CD34_14_2_SA2_529_1252_RAD21 ChIP_CD34_17_2_SA2_529_1252_RAD21 ChIP_CD34_20_5_SA2_529_1252_RAD21 ChIP_CD34_21_3_SA2_529_1252_RAD21 ChIP_CD34_22_2_SA2_529_1252_RAD21 ChIP_CD34_28_5_SA2_529_1252_RAD21 '
RAD21KD_RAD21='ChIP_CD34_18_1_RAD21_467_2031_RAD21 ChIP_CD34_22_1_RAD21_467_2031_RAD21 ChIP_CD34_27_1_RAD21_467_2031_RAD21 ChIP_CD34_28_1_RAD21_467_2031_RAD21 '
RAD21TAGDIRS=$CTRL_RAD21$SA1KD_RAD21$SA2KD_RAD21$RAD21KD_RAD21

##generate a merged DNA-inputseq tagdir from availabe HSPC DNA-inputs (6 donors combined, all were normal karyotype)
###this will be used as background for peakfinding in all samples
cd ${INPUTDIR}
makeTagDirectory ${INPUTDIR}/Input_CD34_merged/ -d ChIP_IN_CD34_16_3 ChIP_IN_CD34_17_3 IN_CD34_18_4_siCtrl IN_CD34_20_6_siCtrl IN_CD34_21_4_siCtrl IN_CD34_27 IN_CD34_28 -genome hg38 -checkGC


# peak calling with default params
cd ${TAGDIR}
for SAMPLE in ${RAD21TAGDIRS}
do
findPeaks ${SAMPLE} -i ${INPUTDIR}/Input_CD34_merged -style factor -o auto
done

#save peak stats to file
echo -e "Sample_ID\tTotalPeaks\tIPefficiency\tTotalTags" > $PEAKDIR/RAD21peaks.default.txt
for SAMPLE in ${RAD21TAGDIRS}; do
    sample_name="${SAMPLE}"
    TotalPeaks=$(grep -w "# total peaks =" $sample_name/peaks.txt | awk '{print $5}')
    IPefficiency=$(grep -w "# Approximate IP efficiency =" $sample_name/peaks.txt | awk '{print $6}')
    TotalTags=$(grep -w "# Total tags =" $sample_name/peaks.txt | awk '{print $5}')
    echo -e "${sample_name}\t${TotalPeaks}\t${IPefficiency}\t${TotalTags}" >> $PEAKDIR/RAD21peaks.default.txt
done


                         ##################################################
                         #       merged tagDirectories by condition       #
                         ##################################################


# all replicates combined by group
cd ${TAGDIR}
makeTagDirectory ChIP_merged_CD34_CTRL_RAD21 -d ${CTRL_RAD21}
makeTagDirectory ChIP_merged_CD34_SA1KD_RAD21 -d ${SA1KD_RAD21}
makeTagDirectory ChIP_merged_CD34_SA2KD_RAD21 -d ${SA2KD_RAD21}
makeTagDirectory ChIP_merged_CD34_RAD21KD_RAD21 -d ${RAD21KD_RAD21}


#  bigWigs from merged tagdirs (meged bigwigs)
cd ${TAGDIR}
declare -a conditions=("CTRL_RAD21" "SA1KD_RAD21" "SA2KD_RAD21" "RAD21KD_RAD21")
for SAMPLE in "${conditions[@]}";do
makeUCSCfile ChIP_merged_CD34_${SAMPLE} -o ${MERGEDBW}/merged.ChIP_${SAMPLE}.bigWig -fragLength 200 -bigWig ${CHROMSIZES_HG38} -fsize 1e20
done


                     ############################################
                     #      peak finding in merged tagDirs      #
                     ############################################
# peak finding
cd ${TAGDIR}
declare -a conditions=("CTRL_RAD21" "SA1KD_RAD21" "SA2KD_RAD21" "RAD21KD_RAD21")
for SAMPLE in "${conditions[@]}";
do
findPeaks ChIP_merged_CD34_${SAMPLE} -i ${INPUTDIR}/Input_CD34_merged -style factor -tbp 1 -fdr 0.000001 -o ${PEAKDIR}/${SAMPLE}.peaks.txt
done

# peak filtering
declare -a PEAKSETS=("CTRL_RAD21" "SA1KD_RAD21" "SA2KD_RAD21" "RAD21KD_RAD21")
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
                     #             merged RAD21 peak sets             #
                     ##################################################
###combine RAD21 peaks from all conditions
mergePeaks ${PEAKDIR}/CTRL_RAD21.filtered.peaks.txt ${PEAKDIR}/SA1KD_RAD21.filtered.peaks.txt \
${PEAKDIR}/SA2KD_RAD21.filtered.peaks.txt ${PEAKDIR}/RAD21KD_RAD21.filtered.peaks.txt -code > ${PEAKDIR}/CD34_RAD21.filtered.peaks.txt
#look at combined peak number
wc -l ${PEAKDIR}/CD34_RAD21.filtered.peaks.txt



                         ###########################################
                         #        scaling of RAD21 bigWigs         #
                         ###########################################
#to adjust tracks for technical variance HSPC data is scaled based on the top peaks found in all samples from the merged peakset

# define top ~1500 peaks for scaling
declare -a PEAKSETS=("CTRL_RAD21" "SA1KD_RAD21" "SA2KD_RAD21" "RAD21KD_RAD21")

cd ${PEAKDIR}
###get ref chromosomes only
for PEAKSET in "${PEAKSETS[@]}";do
awk '$2 ~ /^chr/&& $2 ~ !/^chrM/' ${PEAKDIR}/${PEAKSET}.peaks.txt > ${PEAKDIR}/${PEAKSET}.peaks.refChr.txt
head -n 3000 ${PEAKDIR}/${PEAKSET}.peaks.refChr.txt | cut -f 1-6 > ${PEAKDIR}/${PEAKSET}.peaks.refChr.top.txt
done
###Create a merged peak set of the top peaks
mergePeaks ${PEAKDIR}/CTRL_RAD21.peaks.refChr.top.txt ${PEAKDIR}/SA1KD_RAD21.peaks.refChr.top.txt \
${PEAKDIR}/SA2KD_RAD21.peaks.refChr.top.txt ${PEAKDIR}/RAD21KD_RAD21.peaks.refChr.top.txt -code >${TMPDIR}/merge1.txt
awk -v "key=1111" '$7 == key' ${TMPDIR}/merge1.txt >${PEAKDIR}/all_commonCD34.RAD21.peaks.refChr.top.txt

# calculate median normalization factor for RAD21
cd ${TAGDIR}
annotatePeaks.pl ${PEAKDIR}/all_commonCD34.RAD21.peaks.refChr.top.txt hg38 -size given -d ${RAD21TAGDIRS} -noann -nogene -cpu 12 > "${PEAKDIR}/all_commonCD34.RAD21.peaks.refChr.top.allAnn.txt"
awk 'BEGIN{OFS="\t"} !/chrM/ {sum=0; for(i=8; i<=NF; i++) sum += $i; ave=sum/25; for(j=1;j<=NF;j++){printf "%s\t", $j} print ave }' "${PEAKDIR}/all_commonCD34.RAD21.peaks.refChr.top.allAnn.txt" | tail -n +2 | cut -f 8-33 > ${TMPDIR}/test2.txt
awk 'BEGIN{OFS="\t"} { for(i=1;i<NF;i++){printf "%s\t", $i/$26}  {print ""} }' ${TMPDIR}/test2.txt >  ${TMPDIR}/test3.txt
allMEDIANSRAD21=()
for ((i=1;i<=26;i++));do
allMEDIANSRAD21[$i]=$(cut -f ${i} ${TMPDIR}/test3.txt | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }')
done
echo ${allMEDIANSRAD21[@]}
rm ${TMPDIR}/test2.txt ${TMPDIR}/test3.txt
#1.05493 0.988721 1.36618 1.04408 0.814035 0.493471 1.3158 0.909887 1.0779 1.09701 1.1394 1.56899 1.55507 0.545412 0.772443 1.31872 1.43662 0.795883 1.39445 1.16105 0.907087 0.384381 0.802768 0.348968 0.511227


# scale bigWigs based on the top peaks
cd ${BWDIR}
_DATE=$(date +%s)
COUNT=0
for SAMPLE in ${RAD21TAGDIRS}; do
	#NAME=${NAMES[$COUNT]}
    COUNT2=$((COUNT+1))
    FACTOR=$(awk -v "key=${allMEDIANSRAD21[$COUNT2]}" 'BEGIN {print 1 / key}')
    echo $FACTOR
	cat >"${TMPDIR}/track.${SAMPLE}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.9/bin:${PATH}
export PATH
cd ${TMPDIR}
scaleBedGraph.pl ${BWDIR}/${SAMPLE}.bigwig -M ${FACTOR} -bigWig hg38 -o ${RAD21SCALEDBW}/${SAMPLE}.scaled
EOF
	COUNT=$((COUNT+=1))
	chmod 750 "${TMPDIR}/track.${SAMPLE}.${_DATE}.sh"
	echo "converting track for ${SAMPLE}"
	screen -dm -S track${SAMPLE} bash -c "bash ${TMPDIR}/track.${SAMPLE}.${_DATE}.sh"
done

# loop to check when screen sessions are done
#---------------------------------------------
for SAMPLE in ${RAD21TAGDIRS}; do 
while [ true ]; do     # Endless loop.
  pid=`screen -S track${SAMPLE} -Q echo '$PID'`  # Get a pid.
  if [[ $pid = *"session"* ]] ; then # If there is none,
    echo -e "\tFinished sample ${SAMPLE}"
    break               # Test next one.
  else
    sleep 10           # Else wait.
  fi
	done
done
#---------------------------------------------
# generate average scaled BigWigs: RAD21 ChIP

cd ${RAD21SCALEDBW}
myAverageBigWig.pl -bw ChIP_CD34_28_6_siCtrl_RAD21.scaled.bigWig ChIP_CD34_14_3_siCtrl_RAD21.scaled.bigWig ChIP_CD34_17_3_siCtrl_RAD21.scaled.bigWig ChIP_CD34_18_4_siCtrl_RAD21.scaled.bigWig ChIP_CD34_19_2_siCtrl_RAD21.scaled.bigWig ChIP_CD34_21_4_siCtrl_RAD21.scaled.bigWig ChIP_CD34_22_3_siCtrl_RAD21.scaled.bigWig ChIP_CD34_24_2_siCtrl_RAD21.scaled.bigWig ChIP_CD34_27_4_siCtrl_RAD21.scaled.bigWig ChIP_CD34_14_4_Mock_RAD21.scaled.bigWig \
-chr ${CHROMSIZES_HG38} -o ave.CTRL.RAD21.scaled.bigWig
bigWigToBedGraph ave.CTRL.RAD21.scaled.bigWig ave.CTRL.RAD21.scaled.bedGraph

myAverageBigWig.pl -bw ChIP_CD34_14_2_SA2_529_1252_RAD21.scaled.bigWig ChIP_CD34_17_2_SA2_529_1252_RAD21.scaled.bigWig ChIP_CD34_20_5_SA2_529_1252_RAD21.scaled.bigWig ChIP_CD34_21_3_SA2_529_1252_RAD21.scaled.bigWig ChIP_CD34_22_2_SA2_529_1252_RAD21.scaled.bigWig ChIP_CD34_28_5_SA2_529_1252_RAD21.scaled.bigWig \
-chr ${CHROMSIZES_HG38} -o ave.SA2KD.RAD21.scaled.bigWig
bigWigToBedGraph ave.SA2KD.RAD21.scaled.bigWig ave.SA2KD.RAD21.scaled.bedGraphSA2KDvsCTRL

myAverageBigWig.pl -bw ChIP_CD34_14_1_SA1_2259_4094_RAD21.scaled.bigWig ChIP_CD34_17_1_SA1_2259_4094_RAD21.scaled.bigWig ChIP_CD34_21_2_SA1_2259_4094_RAD21.scaled.bigWig ChIP_CD34_27_3_SA1_2259_4094_RAD21.scaled.bigWig ChIP_CD34_28_4_SA1_2259_4094_RAD21.scaled.bigWig \
-chr ${CHROMSIZES_HG38} -o ave.SA1KD.RAD21.scaled.bigWig
bigWigToBedGraph ave.SA1KD.RAD21.scaled.bigWig ave.SA1KD.RAD21.scaled.bedGraph

myAverageBigWig.pl -bw ChIP_CD34_18_1_RAD21_467_2031_RAD21.scaled.bigWig ChIP_CD34_22_1_RAD21_467_2031_RAD21.scaled.bigWig ChIP_CD34_27_1_RAD21_467_2031_RAD21.scaled.bigWig ChIP_CD34_28_1_RAD21_467_2031_RAD21.scaled.bigWig \
-chr ${CHROMSIZES_HG38} -o ave.RAD21KD.RAD21.scaled.bigWig
bigWigToBedGraph ave.SA1KD.RAD21.scaled.bigWig ave.RAD21KD.RAD21.scaled.bedGraph


                     #############################################################
                     #                Global RAD21 coverage levels               #
                     #############################################################
##average coverage across all RAD21 peaks comparing KD and CTRLs


##for RAD21 KD use merged bw
annotatePeaks.pl ${PEAKDIR}/CD34_RAD21.filtered.peaks.txt hg38 -size 2000 -hist 25 -ghist -bedGraph "${MERGEDBW}/merged.ChIP_CTRL_RAD21.bedGraph" > "${FIGURESDIR}/ghist/CTRL_RAD21.Peaks.allRAD21pos.ghist.txt"
annotatePeaks.pl ${PEAKDIR}/CD34_RAD21.filtered.peaks.txt hg38 -size 2000 -hist 25 -ghist -bedGraph "${MERGEDBW}/merged.ChIP_RAD21KD_RAD21.bedGraph" > "${FIGURESDIR}/ghist/RAD21KD_RAD21.Peaks.allRAD21pos.ghist.txt"
#for STAG KDs use scaled RAD21 bws
annotatePeaks.pl ${PEAKDIR}/CD34_RAD21.filtered.peaks.txt hg38 -size 2000 -hist 25 -ghist -bedGraph "${RAD21SCALEDBW}/ave.SA1KD.RAD21.scaled.bedGraph" > "${FIGURESDIR}/ghist/SA1KD_RAD21.Peaks.allRAD21pos.scaled.ghist.txt"
annotatePeaks.pl ${PEAKDIR}/CD34_RAD21.filtered.peaks.txt hg38 -size 2000 -hist 25 -ghist -bedGraph "${RAD21SCALEDBW}/ave.SA2KD.RAD21.scaled.bedGraph" > "${FIGURESDIR}/ghist/SA2KD_RAD21.Peaks.allRAD21pos.scaled.ghist.txt"
annotatePeaks.pl ${PEAKDIR}/CD34_RAD21.filtered.peaks.txt hg38 -size 2000 -hist 25 -ghist -bedGraph "${RAD21SCALEDBW}/ave.CTRL.RAD21.scaled.bedGraph" > "${FIGURESDIR}/ghist/CTRL_RAD21.Peaks.allRAD21pos.scaled.ghist.txt"


#generate ghist plots
####separate plots for RAD21 KD vs CTRL
plotHIST.sh -g "${FIGURESDIR}/ghist/SA1KD_RAD21.Peaks.allRAD21pos.ghist.txt ${FIGURESDIR}/ghist/RAD21KD_RAD21.Peaks.allRAD21pos.ghist.txt ${FIGURESDIR}/ghist/CTRL_RAD21.Peaks.allRAD21pos.ghist.txt" \
-s "STAG1KD RAD21KD CTRL-HSPCs" -c "darkgoldenrod mediumvioletred firebrick" -x 1000 -y "0 15" -d ${FIGURESDIR}/ghist -n RAD21.merged.allRAD21pos.SA1KDvsRAD21KDvsCTRL 
####separate plots for SA1KD and SA2KD scaled
plotHIST.sh -g "${FIGURESDIR}/ghist/SA1KD_RAD21.Peaks.allRAD21pos.scaled.ghist.txt ${FIGURESDIR}/ghist/CTRL_RAD21.Peaks.allRAD21pos.scaled.ghist.txt" \
-s "STAG1-KD CTRL-HSPCs" -c "darkgoldenrod firebrick" -x 1000 -y "0 12" -d ${FIGURESDIR}/ghist -n RAD21.merged.allRAD21pos.SA2KDvsCTRL.scaled 

plotHIST.sh -g "${FIGURESDIR}/ghist/SA2KD_RAD21.Peaks.allRAD21pos.scaled.ghist.txt ${FIGURESDIR}/ghist/CTRL_RAD21.Peaks.allRAD21pos.scaled.ghist.txt" \
-s "STAG2-KD CTRL-HSPCs" -c "springgreen firebrick" -x 1000 -y "0 12" -d ${FIGURESDIR}/ghist -n RAD21.merged.allRAD21pos.SA1KDvsCTRL.scaled



##look at average tracks by condition in comparison to AML average tracks at an examplary region:
#STAGmut and STAG KDs
pyGenomeTracks_v3.5.sh --tracks ${FIGURESDIRpyg}/inifiles/Pytracks_SA2mut.avbw.STAGKD.scbw.vsCTRLs.RAD21.average.ini --region chr3:36000000-39000000 --fontSize 28 --dpi 600 -o ${FIGURESDIRpyg}/RAD21tracks/Pytracks_STAGmutSTAGKDvsCTRL.RAD21.average.2.png
#RAD21mut and RAD21 KD
pyGenomeTracks_v3.5.sh --tracks ${FIGURESDIRpyg}/inifiles/Pytracks_RADKD_RADmut.avbw.CTRL_HSPCs.RAD21.average.ini --region chr3:36000000-39000000 --fontSize 28 --dpi 600 -o ${FIGURESDIRpyg}/RAD21tracks/Pytracks_RADmutRADKDvsCTRL.RAD21.average.2.png



                     #############################################################
                     #  Annotation of tagdirs to RAD21 peaks in HSPCs            #
                     #############################################################

# annotation of library-size adjusted data (rlog-transformed)

cd ${TAGDIR}
annotatePeaks.pl "${PEAKDIR}/CD34_RAD21.filtered.peaks.txt" hg38 -size 250 -d ${RAD21TAGDIRS} -raw -cpu 12 > ${WORKDIR}/CD34_RAD21.peaks.ann.txt
#remove non-essential cols and add desired colnames
tail -n +2 ${WORKDIR}/CD34_RAD21.peaks.ann.txt | cut -f1,20-44 > ${TMPDIR}/tmp.1.txt
echo $'ID\t14_CTRL\t14_CTRL\t17_CTRL\t18_CTRL\t19_CTRL\t21_CTRL\t22_CTRL\t24_CTRL\t27_CTRL\t28_CTRL\t14_SA1KD\t17_SA1KD\t21_SA1KD\t27_SA1KD\t28_SA1KD\t14_SA2KD\t17_SA2KD\t20_SA2KD\t21_SA2KD\t22_SA2KD\t28_SA2KD\t18_RAD21KD\t22_RAD21KD\t27_RAD21KD\t28_RAD21KD' | cat - ${TMPDIR}/tmp.1.txt > ${WORKDIR}/CD34_RAD21.peaks.ann.Rinput.txt
cut -f 1-6 ${PEAKDIR}/CD34_RAD21.filtered.${PT}.txt > ${PEAKDIR}/CD34_RAD21.peaks.Rinput.txt



                     #############################################################
                     #  Differential RAD21 peaks in RAD21 KD (norm2total)        #
                     #############################################################
##for the RAD21 KD analysis the data needs to be normalized to toal tags; for this we use the HOMER routine:
getDiffExpression.pl ${WORKDIR}/CD34_RAD21.peaks.ann.txt -peaks \
CTRL CTRL CTRL CTRL CTRL CTRL CTRL CTRL CTRL CTRL SA1KD SA1KD SA1KD SA1KD SA1KD SA2KD SA2KD SA2KD SA2KD SA2KD SA2KD RAD21KD RAD21KD RAD21KD RAD21KD \
-batch 14 14 17 18 19 21 22 24 27 28 14 17 21 27 28 14 17 20 21 22 28 18 22 27 28 -rlog -norm2total > ${DIFFPEAKS}/CD34_RAD21.peaks.DESEQnorm2total.txt

#        Output Stats CTRL vs. RAD21KD:                                                                                                     
#                Total Genes: 78867                                                                                                             
#                Total Up-regulated in RAD21KD vs. CTRL: 5 (0.006%) [log2fold>1, FDR<0.05]                                                  
#                Total Dn-regulated in RAD21KD vs. CTRL: 53162 (67.407%) [log2fold<-1, FDR<0.05] 
                
##for regular diffpeak analyis for SA2/SA1 KD run R script!





            ############################################################
            #  scoring CD34 RAD21 ChIPS based on AML-RAD21-Peaks       #
            ############################################################
##this will allow for direct comparison of peak fold-changes etc. and peak set enrichments

cd ${TAGDIR}
annotatePeaks.pl "${AMLPEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.txt" hg38 -size 250 -d ${RAD21TAGDIRS} -raw -cpu 12 > ${WORKDIR}/CD34_AMLpat_mergePeaks_RAD21.filtered.ann.txt

annotatePeaks.pl "${AMLPEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.stringent.txt" hg38 -size 250 -d ${RAD21TAGDIRS} -raw -cpu 12 > ${WORKDIR}/CD34_AMLpat_mergePeaks_RAD21.stringent.filtered.ann.txt
#rm XY chrom peaks - they were excluded in the patient analysis!
#remove sex chromosomes
grep -v 'chrX' ${WORKDIR}/CD34_AMLpat_mergePeaks_RAD21.filtered.ann.txt > ${TMPDIR}/tmp.filt.txt
grep -v 'chrY' ${TMPDIR}/tmp.filt.txt > ${WORKDIR}/CD34_AMLpat_mergePeaks_RAD21.filtered.XYrm.ann.txt #176736
grep -v 'chrX' ${WORKDIR}/CD34_AMLpat_mergePeaks_RAD21.stringent.filtered.ann.txt > ${TMPDIR}/tmp.filt.txt
grep -v 'chrY' ${TMPDIR}/tmp.filt.txt > ${WORKDIR}/CD34_AMLpat_mergePeaks_RAD21.stringent.filtered.XYrm.ann.txt #
rm ${TMPDIR}/tmp.filt*


#prepare stringent (filtered) AML peakset
tail -n +2 ${WORKDIR}/CD34_AMLpat_mergePeaks_RAD21.stringent.filtered.XYrm.ann.txt | cut -f1,20-44 > ${TMPDIR}/tmp.1.txt
echo $'ID\t14_CTRL\t14_CTRL\t17_CTRL\t18_CTRL\t19_CTRL\t21_CTRL\t22_CTRL\t24_CTRL\t27_CTRL\t28_CTRL\t14_SA1KD\t17_SA1KD\t21_SA1KD\t27_SA1KD\t28_SA1KD\t14_SA2KD\t17_SA2KD\t20_SA2KD\t21_SA2KD\t22_SA2KD\t28_SA2KD\t18_RAD21KD\t22_RAD21KD\t27_RAD21KD\t28_RAD21KD' | cat - ${TMPDIR}/tmp.1.txt > ${WORKDIR}/CD34_AMLpat_mergePeaks_RAD21.stringent.filtered.XYrm.ann.Rinput.txt
cut -f 1-6 ${AMLPEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.stringent.txt > ${PEAKDIR}/tmp.Allpat_mergePeaks_RAD21.stringent.txt

#prepare full (filtered) AML peakset
tail -n +2 ${WORKDIR}/CD34_AMLpat_mergePeaks_RAD21.filtered.XYrm.ann.txt | cut -f1,20-44 > ${TMPDIR}/tmp.1.txt
echo $'ID\t14_CTRL\t14_CTRL\t17_CTRL\t18_CTRL\t19_CTRL\t21_CTRL\t22_CTRL\t24_CTRL\t27_CTRL\t28_CTRL\t14_SA1KD\t17_SA1KD\t21_SA1KD\t27_SA1KD\t28_SA1KD\t14_SA2KD\t17_SA2KD\t20_SA2KD\t21_SA2KD\t22_SA2KD\t28_SA2KD\t18_RAD21KD\t22_RAD21KD\t27_RAD21KD\t28_RAD21KD' | cat - ${TMPDIR}/tmp.1.txt > ${WORKDIR}/CD34_AMLpat_mergePeaks_RAD21.filtered.XYrm.ann.Rinput.txt
cut -f 1-6 ${AMLPEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.txt > ${PEAKDIR}/tmp.Allpat_mergePeaks_RAD21.txt



