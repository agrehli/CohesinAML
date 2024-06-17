#bin/bash
#by Alexander Fischer, Jun 2021
###############################################################################
###############################################################################
##                                                                           ##
##     		 Analysis of ATAC data of siRNA KD in CD34+ HSPCs at d4   		 ##
##                                                                           ##
###############################################################################
###############################################################################

#setting basic path
DIR_SOFT="/misc/software"
DIR_PKG="${DIR_SOFT}/ngs"
DIR_DATA="/misc/data"
OS=$(lsb_release -c |grep "^Codename" | awk -F: '{print $2}' | sed 's/[[:blank:]]//g')

#setting homer environment
PATH_PERL=${DIR_SOFT}/package/perl/perl-5.26.1/bin
PATH_SAMTOOLS=${DIR_PKG}/samtools/samtools-1.6/bin
PATH_HOMER=${DIR_PKG}/homer/v4.11/bin
#PATH_R=${DIR_SOFT}/package/RBioC/3.4.3/bin
export PATH=${PATH_PERL}:${PATH_SAMTOOLS}:${PATH_HOMER}:${PATH}
export PATH

# Defining the program versions
BEDTOOLS="${DIR_PKG}/bedtools/bedtools2-2.27.1/bin/bedtools"

# Defining necessary files
CHROMSIZES_HG38="${DIR_SOFT}/viewer/IGV/IGVTools_2.3.98/genomes/GRCh38.PRI_p10.chrom.sizes"
BLACKLIST_HG38="${DIR_SOFT}/analysis/generalStuff/annotation/GRCh38/hg38.blacklist.bed"

# Defining directories
# general directories
TMPDIR="/loctmp"
## directories with data created by mapATAC.sh pipeline
TMPDIR="/loctmp"
RAWDIR="${DIR_DATA}/rawData/chromatin/ATAC/CD34"
FASTQCDIR="${RAWDIR}/FastQC"
TAGDIR="${DIR_DATA}/processedData/tagDir/chromatin/GRCh38/ATAC/CD34"
BWDIR="${DIR_DATA}/processedData/bigWig/chromatin/GRCh38/ATAC/CD34"
## directories for Analysis
PROJECTDIR="${DIR_DATA}/analysis/project_cohesin"
WORKDIR="${PROJECTDIR}/CD34/ATAC"
FIGURESDIR="${WORKDIR}/figures"
PEAKDIR="${WORKDIR}/peaks"
DIFFDIR="${WORKDIR}/diffPeaks"
MOTIFDIR="${WORKDIR}/motifs"
SCALEDBW="${WORKDIR}/scaledBigWigs"
MERGEDBW="${WORKDIR}/mergedBigWigs"

# directory with ChIP peak positions
CHIPPEAKDIRRAD=${PROJECTDIR}/CD34/ChIP_KD_analysis/Cohesin_CTCF_MED12/peaks

# create the new directories
mkdir ${WORKDIR}
mkdir ${PEAKDIR}
mkdir ${FIGURESDIR}
mkdir ${MOTIFDIR}
mkdir ${MERGEDBW}
mkdir ${DIFFDIR}
mkdir ${SCALEDBW}


                     #########################
                     #     intial QC         #
                     #########################
#------------------------------------------------------------------------------------                         
#Declare samples tagdir names
#Ctrl samples
ctrls='ATAC_CD34_14_3_siCtrl ATAC_CD34_17_3_siCtrl ATAC_CD34_18_4_siCtrl ATAC_CD34_19_2_siCtrl ATAC_CD34_20_6_siCtrl ATAC_CD34_21_4_siCtrl ATAC_CD34_20_7_Mock ATAC_CD34_14_4_Mock ATAC_CD34_22_3_siCtrl ATAC_CD34_27_4_siCtrl ATAC_CD34_28_6_siCtrl ATAC_CD34_24_1_siCtrl '
#SA1KD samples
SA1KDs='ATAC_CD34_14_1_SA1_2259_4094 ATAC_CD34_17_1_SA1_2259_4094 ATAC_CD34_21_2_SA1_2259_4094 ATAC_CD34_20_4_SA1_2259_4094 ATAC_CD34_27_3_SA1 ATAC_CD34_28_4_SA1 '
#SA2KD samples
SA2KDs='ATAC_CD34_14_2_SA2_529_1252 ATAC_CD34_17_2_SA2_529_1252 ATAC_CD34_21_3_SA2_529_1252 ATAC_CD34_22_2_SA2_529_1252 ATAC_CD34_20_5_SA2_529_1252 ATAC_CD34_28_5_SA2 '
#RAD21KD samples
RAD21KDs='ATAC_CD34_18_1_Rad21 ATAC_CD34_20_1_Rad21 ATAC_CD34_22_1_RAD21 ATAC_CD34_27_1_Rad21 ATAC_CD34_28_1_Rad21 '
#SA1_SA2KD samples
SA1SA2KDs='ATAC_CD34_20_2_SA1_SA2 ATAC_CD34_28_2_SA1_SA2 ATAC_CD34_18_2_SA1_SA2 '
ATACTAGDIRS=${ctrls}${SA1KDs}${SA2KDs}${RAD21KDs}${SA1SA2KDs}
#------------------------------------------------------------------------------------                         
#look at peak finding QC (generated per default by mapATAC pipeline)
cd ${TAGDIR}
for SAMPLE in ${ATACTAGDIRS}
do
paste <(echo -e "# ${SAMPLE}") <(grep -w "# total peaks =" ${SAMPLE}/peaks.txt) <(grep -w "# Approximate IP efficiency"  ${SAMPLE}/peaks.txt) <(grep -w "# Total tags ="  ${SAMPLE}/peaks.txt)
done
#------------------------------------------------------------------------------------                         



                     #####################################################
                     #     Peak finding using  findATACpeaks2.0.sh       #
                     #####################################################
# running findATACpeaks2.0.sh script for all samples individually
cd ${TAGDIR}
for SAMPLE in ${ctrls}
do
findATACpeaks2.0.sh -t ${TAGDIR}/${SAMPLE} -o ${PEAKDIR} -g hg38 -s
done

for SAMPLE in ${SA1KDs}
do
findATACpeaks2.0.sh -t ${TAGDIR}/${SAMPLE} -o ${PEAKDIR} -g hg38 -s
done

for SAMPLE in ${SA2KDs}
do
findATACpeaks2.0.sh -t ${TAGDIR}/${SAMPLE} -o ${PEAKDIR} -g hg38 -s
done

for SAMPLE in ${RAD21KDs}
do
findATACpeaks2.0.sh -t ${TAGDIR}/${SAMPLE} -o ${PEAKDIR} -g hg38 -s
done

for SAMPLE in ${SA1SA2KDs}
do
findATACpeaks2.0.sh -t ${TAGDIR}/${SAMPLE} -o ${PEAKDIR} -g hg38 -s
done
#------------------------------------------------------------------------------------                         



                     #####################################################
                     #     Combined TAGDIRECTORIES by condition         #
                     #####################################################

# combined directories for peak finding and bigwig scaling
#------------------------------------------------------------------------------------                         

cd ${TAGDIR}
makeTagDirectory CTRL_combined -d ${ctrls}
findATACpeaks2.0.sh -t CTRL_combined -o ${PEAKDIR} -g hg38 -s

makeTagDirectory SA1KD_combined -d ${SA1KDs}       
findATACpeaks2.0.sh -t SA1KD_combined -o ${PEAKDIR} -g hg38 -s

makeTagDirectory SA2KD_combined -d ${SA2KDs}
findATACpeaks2.0.sh -t SA2KD_combined -o ${PEAKDIR} -g hg38 -s

makeTagDirectory RAD21KD_combined -d ${RAD21KDs}
findATACpeaks2.0.sh -t RAD21KD_combined -o ${PEAKDIR} -g hg38 -s

makeTagDirectory SA1SA2KD_combined -d ${SA1SA2KDs}
findATACpeaks2.0.sh -t SA1SA2KD_combined -o ${PEAKDIR} -g hg38 -s
#------------------------------------------------------------------------------------                         

# combined tag directory merged bigwigs
cd ${TAGDIR}
declare -a conditions=("SA1KD" "SA2KD" "RAD21KD" "CTRL" "SA1SA2KD")
for SAMPLE in "${conditions[@]}";do
makeUCSCfile ${SAMPLE}_combined -o ${MERGEDBW}/merged.ATAC_${SAMPLE}.bigWig -fragLength 200 -bigWig ${CHROMSIZES_HG38} -fsize 1e20
bigWigToBedGraph ${MERGEDBW}/merged.ATAC_${SAMPLE}.bigWig ${MERGEDBW}/merged.ATAC_${SAMPLE}.bedGraph
done
#------------------------------------------------------------------------------------                         



                     #####################################################
                     #  Scaling of ATAC bigwigs and average bigwigs      #
                     #####################################################
# define top ~5000 peaks for scaling
#------------------------------------------------------------------------------------                         

awk '$2 ~ /^chr/&& $2 ~ !/^chrM/' ${PEAKDIR}/CTRL_combined.s150.md250.L2.fdr-5.region.peaks.txt > ${PEAKDIR}/CTRL_combined.s150.md250.L2.fdr-5.region.refChr.txt
awk '$2 ~ /^chr/&& $2 ~ !/^chrM/' ${PEAKDIR}/SA1KD_combined.s150.md250.L2.fdr-5.region.peaks.txt > ${PEAKDIR}/SA1KD_combined.s150.md250.L2.fdr-5.region.refChr.txt
awk '$2 ~ /^chr/&& $2 ~ !/^chrM/' ${PEAKDIR}/SA2KD_combined.s150.md250.L2.fdr-5.region.peaks.txt > ${PEAKDIR}/SA2KD_combined.s150.md250.L2.fdr-5.region.refChr.txt
awk '$2 ~ /^chr/&& $2 ~ !/^chrM/' ${PEAKDIR}/RAD21KD_combined.s150.md250.L2.fdr-5.region.peaks.txt > ${PEAKDIR}/RAD21KD_combined.s150.md250.L2.fdr-5.region.refChr.txt
awk '$2 ~ /^chr/&& $2 ~ !/^chrM/' ${PEAKDIR}/SA1SA2KD_combined.s150.md250.L2.fdr-5.region.peaks.txt > ${PEAKDIR}/SA1SA2KD_combined.s150.md250.L2.fdr-5.region.refChr.txt


head -n 5000 ${PEAKDIR}/CTRL_combined.s150.md250.L2.fdr-5.region.refChr.txt | cut -f 1-6 > ${PEAKDIR}/CTRL_combined.s150.md250.L2.fdr-5.region.top5000.txt
head -n 5000 ${PEAKDIR}/SA1KD_combined.s150.md250.L2.fdr-5.region.refChr.txt | cut -f 1-6 > ${PEAKDIR}/SA1KD_combined.s150.md250.L2.fdr-5.region.top5000.txt
head -n 5000 ${PEAKDIR}/SA2KD_combined.s150.md250.L2.fdr-5.region.refChr.txt | cut -f 1-6 > ${PEAKDIR}/SA2KD_combined.s150.md250.L2.fdr-5.region.top5000.txt
head -n 5000 ${PEAKDIR}/RAD21KD_combined.s150.md250.L2.fdr-5.region.refChr.txt | cut -f 1-6 > ${PEAKDIR}/RAD21KD_combined.s150.md250.L2.fdr-5.region.top5000.txt
head -n 5000 ${PEAKDIR}/SA1SA2KD_combined.s150.md250.L2.fdr-5.region.refChr.txt | cut -f 1-6 > ${PEAKDIR}/SA1SA2KD_combined.s150.md250.L2.fdr-5.region.top5000.txt
#------------------------------------------------------------------------------------                         

# geneerate merged ATAC peak set HSPCs
mergePeaks ${PEAKDIR}/CTRL_combined.s150.md250.L2.fdr-5.region.top5000.txt ${PEAKDIR}/SA1KD_combined.s150.md250.L2.fdr-5.region.top5000.txt ${PEAKDIR}/SA2KD_combined.s150.md250.L2.fdr-5.region.top5000.txt \
${PEAKDIR}/RAD21KD_combined.s150.md250.L2.fdr-5.region.top5000.txt ${PEAKDIR}/SA1SA2KD_combined.s150.md250.L2.fdr-5.region.top5000.txt -code >${PEAKDIR}/merge1.txt

#filter for peaks that are in all conditions of the merged peakset: column $7 contains a "1" for each condition that has the peak --> for 6 conditions : 111111
awk -v "key=11111" '$7 == key' ${PEAKDIR}/merge1.txt >${PEAKDIR}/all_commonCD34UC.s150.md250.L2.fdr-5.region.top3300.txt

##annotate common top peaks
annotatePeaks.pl ${PEAKDIR}/all_commonCD34UC.s150.md250.L2.fdr-5.region.top3300.txt hg38 -size given -d ${ATACTAGDIRS} -noann -nogene -cpu 12 > "${PEAKDIR}/all_commonCD34UC.s150.md250.L2.fdr-5.region.top3300.allAnn.txt"

##calculate average values, adjust for sample number (here: 37)
awk 'BEGIN{OFS="\t"} !/chrM/ {sum=0; for(i=8; i<=NF; i++) sum += $i; ave=sum/37; for(j=1;j<=NF;j++){printf "%s\t", $j} print ave }' "${PEAKDIR}/all_commonCD34UC.s150.md250.L2.fdr-5.region.top3300.allAnn.txt" | cut -f 8-45 | head
#awk 'BEGIN{OFS="\t"} !/chrM/ {sum=0; for(i=8; i<=NF; i++) sum += $i; ave=sum/37; for(j=1;j<=NF;j++){printf "%s\t", $j} print ave }' "${PEAKDIR}/all_commonCD34UC.s150.md250.L2.fdr-5.region.top3300.allAnn.txt" | tail -n +2 | cut -f 8-45 > ${TMPDIR}/test2.txt
awk 'BEGIN{OFS="\t"} { for(i=1;i<NF;i++){printf "%s\t", $i/$38}  {print ""} }' ${TMPDIR}/test2.txt >  ${TMPDIR}/test3.txt
allMEDIANS=()
for ((i=1;i<=37;i++));do
allMEDIANS[$i]=$(cut -f ${i} ${TMPDIR}/test3.txt | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }')
done
echo ${allMEDIANS[@]}
rm ${TMPDIR}/test2.txt ${TMPDIR}/test3.txt

# scale bigWigs based on the top peaks
cd ${BWDIR}
_DATE=$(date +%s)
COUNT=0
for SAMPLE in ${ATACTAGDIRS}; do
    COUNT2=$((COUNT+1))
    FACTOR=$(awk -v "key=${allMEDIANS[$COUNT2]}" 'BEGIN {print 1 / key}')
    echo -e "$SAMPLE $FACTOR\n"
cd ${TMPDIR}
scaleBedGraph.pl ${BWDIR}/${SAMPLE}.bigwig -M ${FACTOR} -bigWig hg38 -o ${SCALEDBW}/${SAMPLE}.scaled
	COUNT=$((COUNT+=1))
done
#------------------------------------------------------------------------------------                         


# generate average scaled BigWigs
cd ${SCALEDBW}
myAverageBigWig.pl -bw ATAC_CD34_14_3_siCtrl.scaled.bigWig ATAC_CD34_17_3_siCtrl.scaled.bigWig ATAC_CD34_18_4_siCtrl.scaled.bigWig ATAC_CD34_19_2_siCtrl.scaled.bigWig ATAC_CD34_20_6_siCtrl.scaled.bigWig ATAC_CD34_21_4_siCtrl.scaled.bigWig ATAC_CD34_20_7_Mock.scaled.bigWig ATAC_CD34_14_4_Mock.scaled.bigWig ATAC_CD34_22_3_siCtrl.scaled.bigWig ATAC_CD34_27_4_siCtrl.scaled.bigWig ATAC_CD34_28_6_siCtrl.scaled.bigWig ATAC_CD34_24_1_siCtrl.scaled.bigWig -chr ${CHROMSIZES_HG38} -o ${SCALEDBW}/ave.ctrl.scaled.bigWig
myAverageBigWig.pl -bw ATAC_CD34_14_1_SA1_2259_4094.scaled.bigWig ATAC_CD34_17_1_SA1_2259_4094.scaled.bigWig ATAC_CD34_21_2_SA1_2259_4094.scaled.bigWig ATAC_CD34_20_4_SA1_2259_4094.scaled.bigWig ATAC_CD34_27_3_SA1.scaled.bigWig ATAC_CD34_28_4_SA1.scaled.bigWig -chr ${CHROMSIZES_HG38} -o ${SCALEDBW}/ave.SA1KD.scaled.bigWig
myAverageBigWig.pl -bw ATAC_CD34_14_2_SA2_529_1252.scaled.bigWig ATAC_CD34_17_2_SA2_529_1252.scaled.bigWig ATAC_CD34_21_3_SA2_529_1252.scaled.bigWig ATAC_CD34_22_2_SA2_529_1252.scaled.bigWig ATAC_CD34_20_5_SA2_529_1252.scaled.bigWig ATAC_CD34_28_5_SA2.scaled.bigWig -chr ${CHROMSIZES_HG38} -o ${SCALEDBW}/ave.SA2KD.scaled.bigWig
myAverageBigWig.pl -bw ATAC_CD34_18_1_Rad21.scaled.bigWig ATAC_CD34_20_1_Rad21.scaled.bigWig ATAC_CD34_22_1_RAD21.scaled.bigWig ATAC_CD34_27_1_Rad21.scaled.bigWig ATAC_CD34_28_1_Rad21.scaled.bigWig -chr ${CHROMSIZES_HG38} -o ${SCALEDBW}/ave.RAD21KD.scaled.bigWig
myAverageBigWig.pl -bw ATAC_CD34_20_2_SA1_SA2.scaled.bigWig ATAC_CD34_28_2_SA1_SA2.scaled.bigWig ATAC_CD34_18_2_SA1_SA2.scaled.bigWig -chr ${CHROMSIZES_HG38} -o ${SCALEDBW}/ave.SA1SA2KD.scaled.bigWig

bigWigToBedGraph ave.ctrl.scaled.bigWig ave.ctrl.scaled.bedGraph
bigWigToBedGraph ave.SA1KD.scaled.bigWig ave.SA1KD.scaled.bedGraph
bigWigToBedGraph ave.SA2KD.scaled.bigWig ave.SA2KD.scaled.bedGraph
bigWigToBedGraph ave.RAD21KD.scaled.bigWig ave.RAD21KD.scaled.bedGraph
bigWigToBedGraph ave.SA1SA2KD.scaled.bigWig ave.SA1SA2KD.scaled.bedGraph
#------------------------------------------------------------------------------------                         





                     ############################################################################
                     #  Generate merged ATAC peakset of all conditions and annotate TAGIDRS     #
                     ############################################################################
cd ${TAGDIR}
#------------------------------------------------------------------------------------      
# merged set for batch correction and MDS etc.
mergePeaks ${PEAKDIR}/CTRL_combined.intersect.peaks.txt ${PEAKDIR}/SA1KD_combined.intersect.peaks.txt ${PEAKDIR}/SA2KD_combined.intersect.peaks.txt ${PEAKDIR}/SA1SA2KD_combined.intersect.peaks.txt ${PEAKDIR}/RAD21KD_combined.intersect.peaks.txt -code > "${PEAKDIR}/CD34_combined.intersect.peaks.txt"
# annotate all TAGDIRS
annotatePeaks.pl "${PEAKDIR}/CD34_combined.intersect.peaks.txt" hg38 -size given -len 65 -d ${ATACTAGDIRS2} -noann -nogene -noadj -cpu 12 > "${PEAKDIR}/CD34_combined.intersect.peaks.CD34.ann.txt"
tail -n +2 "${PEAKDIR}/CD34_combined.intersect.peaks.CD34.ann.txt" | cut -f1,8-40 > ${WORKDIR}/tmp.1.txt
echo $'ID\tATAC_CD34_14_3_siCtrl\tATAC_CD34_17_3_siCtrl\tATAC_CD34_18_4_siCtrl\tATAC_CD34_19_2_siCtrl\tATAC_CD34_20_6_siCtrl\tATAC_CD34_21_4_siCtrl\tATAC_CD34_20_7_Mock\tATAC_CD34_14_4_Mock\tATAC_CD34_22_3_siCtrl\tATAC_CD34_27_4_siCtrl\tATAC_CD34_28_6_siCtrl\tATAC_CD34_24_1_siCtrl\tATAC_CD34_14_1_SA1_2259_4094\tATAC_CD34_17_1_SA1_2259_4094\tATAC_CD34_21_2_SA1_2259_4094\tATAC_CD34_20_4_SA1_2259_4094\tATAC_CD34_27_3_SA1\tATAC_CD34_28_4_SA1\tATAC_CD34_14_2_SA2_529_1252\tATAC_CD34_17_2_SA2_529_1252\tATAC_CD34_21_3_SA2_529_1252\tATAC_CD34_22_2_SA2_529_1252\tATAC_CD34_20_5_SA2_529_1252\tATAC_CD34_28_5_SA2\tATAC_CD34_18_1_Rad21\tATAC_CD34_20_1_Rad21\tATAC_CD34_22_1_RAD21\tATAC_CD34_27_1_Rad21\tATAC_CD34_28_1_Rad21\tATAC_CD34_20_2_SA1_SA2\tATAC_CD34_28_2_SA1_SA2\tATAC_CD34_18_2_SA1_SA2' | cat - ${WORKDIR}/tmp.1.txt > ${WORKDIR}/CD34_combined.intersect.peaks.CD34.ann.table.txt
# get ref chromosomes only
awk '!/^chrM/&&/chr/' ${PEAKDIR}/CD34_combined.intersect.peaks.txt > ${PEAKDIR}/CD34_combined.intersect.peaks.refChr.txt
# annotate GC and length
annotatePeaks.pl ${PEAKDIR}/CD34_combined.intersect.peaks.refChr.txt hg38 -size given -noann -nogene -CpG -cpu 12 > "${PEAKDIR}/CD34_combined.intersect.peaks.tmp.CpGann.txt"
rm "${PEAKDIR}/CD34_combined.intersect.peaks.CpGann.txt"
touch "${PEAKDIR}/CD34_combined.intersect.peaks.CpGann.txt"
echo $'ID\tlength\tgccontent' >> "${PEAKDIR}/CD34_combined.intersect.peaks.CpGann.txt"
awk -v OFS='\t' '{print $1,$4-$3,$9 ; }' <(tail -n+2 "${PEAKDIR}/CD34_combined.intersect.peaks.tmp.CpGann.txt") >> "${PEAKDIR}/CD34_combined.intersect.peaks.CpGann.txt"
rm "${PEAKDIR}/CD34_combined.intersect.peaks.tmp.CpGann.txt"
#------------------------------------------------------------------------------------      


                         #####################################################################
                         #  Annotation a of ATAC data based tp RAD21 peak positions in HSPCs #
                         #####################################################################
#------------------------------------------------------------------------------------      
######using all CD34 RAD21 merged peak
#prepare annotation table using merged peak set
cd ${TAGDIR}
annotatePeaks.pl "${CHIPPEAKDIRRAD}/CD34_RAD21.filtered.peaks.txt" hg38 -size given -len 65 -d ${ATACTAGDIRS2} -noann -nogene -noadj -cpu 12 > "${PEAKDIR}/ATAC_CD34_RAD21.peaks.ann.txt" 
tail -n +2 "${PEAKDIR}/ATAC_CD34_RAD21.peaks.ann.txt"  | cut -f1,8-100 > "${PEAKDIR}/ATAC_CD34_RAD21.peaks.ann.tab.txt" #78867
#------------------------------------------------------------------------------------      







