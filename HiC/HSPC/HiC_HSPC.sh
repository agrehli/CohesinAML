#!/bin/bash
#by Alexander Fischer, AUG 2021

#################################################################################################
#         Analysis of deep sequenced HiC data of Cohesin KD experiments in HSPCs AML samples    #
#################################################################################################

# general paths
TMPDIR="/loctmp"
DIR_PKG="/misc/software/ngs"
DIR_SOFT="/misc/software"
DIR_PKG="${DIR_SOFT}/ngs"
DIR_DATA="/misc/data"
OS=$(lsb_release -c |grep "^Codename" | awk -F: '{print $2}' | sed 's/[[:blank:]]//g')

#setting homer environment
# setting homer environment
DIR_PKG="/misc/software/ngs"
PATH_PERL=${DIR_SOFT}/package/perl/perl-5.26.1/bin
PATH_SAMTOOLS=${DIR_PKG}/samtools/samtools-1.6/bin
PATH_HOMER=${DIR_PKG}/homer/v4.10/bin
PATH_R=${DIR_SOFT}/package/RBioC/3.4.3/bin
PATH_SKEWER=${DIR_PKG}/skewer/skewer-0.2.2
export PATH=${PATH_R}:${PATH_PERL}:${PATH_SAMTOOLS}:${PATH_HOMER}:${PATH_SKEWER}:${PATH}
export PATH
# important tools and files
BEDTOOLS="${DIR_PKG}/bedtools/bedtools2-2.27.1/bin/bedtools"
BOWTIE2="${DIR_PKG}/bowtie/bowtie2-2.3.4-linux-x86_64/bowtie2"
CHROMSIZES_HG38="${DIR_SOFT}/viewer/IGV/IGVTools_2.3.98/genomes/GRCh38.PRI_p10.chrom.sizes"
BAD_HG38="${DIR_DATA}/analysis/generalStuff/annotation/GRCh38/hg38.badRegions.bed"

# Directories for HiC Processing and Analysis
## Processing
RAWHICDATA_CD34="${DIR_DATA}/rawData/3Dchromatin/HiC/CD34"
MAPHICDIR_CD34="${DIR_DATA}/processedData/mapping/3Dchromatin/HiC/GRCh38/CD34"
TAGHICDIR_CD34="${DIR_DATA}/processedData/tagDir/3Dchromatin/HiC/GRCh38/CD34"

## Analysis
WORKDIR_CD34="${DIR_DATA}/analysis/project_cohesin/CD34/HiC_KDs/DeepSeq_Analysis"
FIGURESDIR="${WORKDIR_CD34}/figures"
MATRIXDIR_CD34="${WORKDIR_CD34}/Interactionmatrices"
PCADIR="${WORKDIR_CD34}/PCA"
QCDIR_CD34="${WORKDIR_CD34}/QC"
LOOPDIR_CD34="${WORKDIR_CD34}/loops"

## Data for comparison with AML
AMLLOOPDIR="${DIR_DATA}/analysis/project_cohesin/Cohesin_AML/HiC/loops"

## conditions
KDs="SA1KD SA2KD RAD21KD"

## Define sample names in vectors by condition
HIC_CTRL='HiC_CD34_14_3_siCtrl HiC_CD34_17_3_siCtrl HiC_CD34_18_4_siCtrl HiC_CD34_20_6_siCtrl_Rep1 HiC_CD34_21_4_siCtrl_Rep1 HiC_CD34_22_3_siCtrl HiC_CD34_27_4_siCtrl HiC_CD34_28_6_siCtrl ' 
HIC_SA1KD='HiC_CD34_14_1_SA1_KD HiC_CD34_17_1_SA1_KD HiC_CD34_20_4_SA1_KD HiC_CD34_21_2_SA1_KD HiC_CD34_27_3_SA1_KD HiC_CD34_28_4_SA1_KD '
HIC_SA2KD='HiC_CD34_14_2_SA2_KD HiC_CD34_17_2_SA2_KD HiC_CD34_20_5_SA2_KD HiC_CD34_21_3_SA2_KD HiC_CD34_22_2_SA2_KD HiC_CD34_28_5_SA2_KD '
HIC_RAD21KD='HiC_CD34_18_1_RAD21_KD HiC_CD34_20_1_RAD21_KD HiC_CD34_22_1_RAD21_KD HiC_CD34_27_1_RAD21_KD HiC_CD34_28_1_RAD21_KD '
HIC_CD34=$HIC_CTRL$HIC_SA1KD$HIC_SA2KD$HIC_RAD21KD

##create new directories
mkdir ${WORKDIR_CD34}
mkdir ${FIGURESDIR}
mkdir ${FIGURESDIR}/PC1_Clustering
mkdir ${FIGURESDIR}/Loop_Clustering
mkdir ${FIGURESDIR}/TAD_Clustering
mkdir ${PCADIR}
mkdir ${QCDIR_CD34}
mkdir ${LOOPDIR_CD34}



#------------------------------------------------------------------------------------

                         #####################################
                         #         Sequencing QC             #
                         #####################################

#run fastqc on all files
for sn in ${HIC_CD34};
do
fastqc-0.11.7 ${RAWHICDATA_CD34}/${sn}_R1.fastq.gz -o ${RAWHICDATA_CD34}/FastQC
fastqc-0.11.7 ${RAWHICDATA_CD34}/${sn}_R2.fastq.gz -o ${RAWHICDATA_CD34}/FastQC
done

#multiqc
cd ${RAWHICDATA_CD34}/FastQC
multiqc_1.13 *.fastq.* \
-n HiC_deepseq1_CD34_multiQC -f -o ${RAWHICDATA_CD34}/FastQC/
#------------------------------------------------------------------------------------




                         ##################################################
                         #    Processing and Mapping of the raw data      #
                         ##################################################

#------------------------------------------------------------------------------------ 
mkdir ${QCDIR_CD34}/taginfo_unfiltered/

for sn in ${HIC_CD34};
do
#trimming
homerTools trim -3 GATC -matchStart 20 -min 20 ${RAWHICDATA_CD34}/${sn}_R1.fastq.gz -stats ${RAWHICDATA_CD34}/${sn}_R1.stats.txt
homerTools trim -3 GATC -matchStart 20 -min 20 ${RAWHICDATA_CD34}/${sn}_R2.fastq.gz -stats ${RAWHICDATA_CD34}/${sn}_R2.stats.txt
#mapping using bowtie
${BOWTIE2} -x GRCh38.PRI_p10.refOnly/hg38 -p 8 ${RAWHICDATA_CD34}/${sn}_R1.fastq.gz.trimmed >${MAPHICDIR_CD34}/${sn}_R1.sam 2>${MAPHICDIR_CD34}/${sn}_R1.log
${BOWTIE2} -x GRCh38.PRI_p10.refOnly/hg38 -p 8 ${RAWHICDATA_CD34}/${sn}_R2.fastq.gz.trimmed >${MAPHICDIR_CD34}/${sn}_R2.sam 2>${MAPHICDIR_CD34}/${sn}_R2.log
#tagDir creation with #tbp1: Removes Clonal/PCR duplicates
makeTagDirectory ${TAGHICDIR_CD34}/${sn} ${MAPHICDIR_CD34}/${sn}_R1.sam,${MAPHICDIR_CD34}/${sn}_R2.sam -tbp 1 -genome hg38 -checkGC
cp -r ${TAGHICDIR_CD34}/${sn} ${TAGHICDIR_CD34}/${sn}_filtered
#tagDir filtering: remove self ligations and spikes
makeTagDirectory ${TAGHICDIR_CD34}/${sn}_filtered -update -restrictionSite GATC -both -genome hg38 -removePEbg -removeSelfLigation -removeSpikes 10000 5
echo "cleaning up for ${SAMPLENAME}"
samtools view -S -@12 -b ${MAPHICDIR_CD34}/${sn}_R1.sam > ${MAPHICDIR_CD34}/${sn}_R1.bam
samtools view -S -@12 -b ${MAPHICDIR_CD34}/${sn}_R2.sam > ${MAPHICDIR_CD34}/${sn}_R2.bam
samtools index ${MAPHICDIR_CD34}/${sn}_R1.bam
samtools index ${MAPHICDIR_CD34}/${sn}_R2.bam
rm ${MAPHICDIR_CD34}/${sn}_R1.sam ${MAPHICDIR_CD34}/${sn}_R2.sam
rm ${RAWHICDATA_CD34}/${sn}_R1.fastq.gz.trimmed ${RAWHICDATA_CD34}/${sn}_R2.fastq.gz.trimmed
mv ${TAGHICDIR_CD34}/${sn}/tagInfo.txt ${QCDIR_CD34}/taginfo_unfiltered/${sn}_tagInfo.txt
rm -r ${TAGHICDIR_CD34}/${sn}
echo "finished with ${SAMPLENAME}"
done


                         ##################################################
                         #         QC for HOMER TAG directories           #
                         ##################################################

#------------------------------------------------------------------------------------
### show QC of filtered and unfiltered tagDIRS
cd ${TAGHICDIR_CD34}
for SAMPLE in ${HIC_CD34}
do
paste <(echo -e "# ${SAMPLE}") <(grep -w "genome=hg38" ${QCDIR_CD34}/taginfo_unfiltered/${SAMPLE}_tagInfo.txt) <(grep "localInteractionFraction=" ${QCDIR_CD34}/taginfo_unfiltered/${SAMPLE}_tagInfo.txt) <(grep "interChrInteractionFraction=" ${QCDIR_CD34}/taginfo_unfiltered/${SAMPLE}_tagInfo.txt) <(grep "tagsPerBP=" ${QCDIR_CD34}/taginfo_unfiltered/${SAMPLE}_tagInfo.txt)
paste <(echo -e "# ${SAMPLE}_filtered") <(grep -w "genome=hg38" ${SAMPLE}_filtered/tagInfo.txt) <(grep "localInteractionFraction=" ${SAMPLE}_filtered/tagInfo.txt) <(grep "interChrInteractionFraction=" ${SAMPLE}_filtered/tagInfo.txt) <(grep "tagsPerBP=" ${SAMPLE}_filtered/tagInfo.txt)
done

#------------------------------------------------------------------------------------



                      ##########################################
                      #                                        #
                      #  Chromatin Compartment Analysis (PCA)  #
                      #                                        #
                      ##########################################

cd ${TAGHICDIR_CD34}
## run PCA for 50 kb windows with 25 kb resolution
for sn in ${HIC_CD34};do
runHiCpca.pl auto ${sn}_filtered -res 25000 -window 50000 -genome hg38 -cpu 24
bedGraphToBigWig ${sn}_filtered/${sn}_filtered.25x50kb.PC1.bedGraph ${CHROMSIZES_HG38} ${PCADIR}/${sn}.25x50kb.PC1.bigWig
done
#------------------------------------------------------------------------------------
## generate average tracks  for 25x50 PC1 bedgraphs
cd ${TAGHICDIR_CD34}
bedgraphs25x50_ctrl='HiC_CD34_14_3_siCtrl_filtered/HiC_CD34_14_3_siCtrl_filtered.25x50kb.PC1.bedGraph HiC_CD34_17_3_siCtrl_filtered/HiC_CD34_17_3_siCtrl_filtered.25x50kb.PC1.bedGraph HiC_CD34_18_4_siCtrl_filtered/HiC_CD34_18_4_siCtrl_filtered.25x50kb.PC1.bedGraph HiC_CD34_20_6_siCtrl_Rep1_filtered/HiC_CD34_20_6_siCtrl_Rep1_filtered.25x50kb.PC1.bedGraph HiC_CD34_21_4_siCtrl_Rep1_filtered/HiC_CD34_21_4_siCtrl_Rep1_filtered.25x50kb.PC1.bedGraph HiC_CD34_22_3_siCtrl_filtered/HiC_CD34_22_3_siCtrl_filtered.25x50kb.PC1.bedGraph HiC_CD34_27_4_siCtrl_filtered/HiC_CD34_27_4_siCtrl_filtered.25x50kb.PC1.bedGraph HiC_CD34_28_6_siCtrl_filtered/HiC_CD34_28_6_siCtrl_filtered.25x50kb.PC1.bedGraph'
bedgraphs25x50_RAD21KD='HiC_CD34_18_1_RAD21_KD_filtered/HiC_CD34_18_1_RAD21_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_20_1_RAD21_KD_filtered/HiC_CD34_20_1_RAD21_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_22_1_RAD21_KD_filtered/HiC_CD34_22_1_RAD21_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_27_1_RAD21_KD_filtered/HiC_CD34_27_1_RAD21_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_28_1_RAD21_KD_filtered/HiC_CD34_28_1_RAD21_KD_filtered.25x50kb.PC1.bedGraph'
bedgraphs25x50_SA2KD='HiC_CD34_14_2_SA2_KD_filtered/HiC_CD34_14_2_SA2_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_17_2_SA2_KD_filtered/HiC_CD34_17_2_SA2_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_20_5_SA2_KD_filtered/HiC_CD34_20_5_SA2_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_21_3_SA2_KD_filtered/HiC_CD34_21_3_SA2_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_22_2_SA2_KD_filtered/HiC_CD34_22_2_SA2_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_28_5_SA2_KD_filtered/HiC_CD34_28_5_SA2_KD_filtered.25x50kb.PC1.bedGraph'
bedgraphs25x50_SA1KD='HiC_CD34_14_1_SA1_KD_filtered/HiC_CD34_14_1_SA1_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_17_1_SA1_KD_filtered/HiC_CD34_17_1_SA1_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_20_4_SA1_KD_filtered/HiC_CD34_20_4_SA1_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_21_2_SA1_KD_filtered/HiC_CD34_21_2_SA1_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_27_3_SA1_KD_filtered/HiC_CD34_27_3_SA1_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_28_4_SA1_KD_filtered/HiC_CD34_28_4_SA1_KD_filtered.25x50kb.PC1.bedGraph'

cd ${TAGHICDIR_CD34}
$BEDTOOLS unionbedg -i ${bedgraphs25x50_ctrl} | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' > ${PCADIR}/ave.CD34_siCtrl_50KB.PC1.bedGraph         
bedGraphToBigWig  ${PCADIR}/ave.CD34_siCtrl_50KB.PC1.bedGraph  ${CHROMSIZES_HG38} ${PCADIR}/ave.CD34_siCtrl_50KB.PC1.bigWig

cd ${TAGHICDIR_CD34}
$BEDTOOLS unionbedg -i ${bedgraphs25x50_RAD21KD} | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' > ${PCADIR}/ave.CD34_RAD21KD_50KB.PC1.bedGraph         
bedGraphToBigWig  ${PCADIR}/ave.CD34_RAD21KD_50KB.PC1.bedGraph  ${CHROMSIZES_HG38} ${PCADIR}/ave.CD34_RAD21KD_50KB.PC1.bigWig

cd ${TAGHICDIR_CD34}
$BEDTOOLS unionbedg -i ${bedgraphs25x50_SA1KD} | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' > ${PCADIR}/ave.CD34_SA1KD_50KB.PC1.bedGraph         
bedGraphToBigWig  ${PCADIR}/ave.CD34_SA1KD_50KB.PC1.bedGraph  ${CHROMSIZES_HG38} ${PCADIR}/ave.CD34_SA1KD_50KB.PC1.bigWig

cd ${TAGHICDIR_CD34}
$BEDTOOLS unionbedg -i ${bedgraphs25x50_SA2KD} | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' > ${PCADIR}/ave.CD34_SA2KD_50KB.PC1.bedGraph         
bedGraphToBigWig  ${PCADIR}/ave.CD34_SA2KD_50KB.PC1.bedGraph  ${CHROMSIZES_HG38} ${PCADIR}/ave.CD34_SA2KD_50KB.PC1.bigWig
#------------------------------------------------------------------------------------
# For clustering: Anotate all bedgraph files (#50kb windows) to PC1 regions (any sample from the analysis can be chosen as reference: here 14_1)
bedgraphs25x50='HiC_CD34_14_3_siCtrl_filtered/HiC_CD34_14_3_siCtrl_filtered.25x50kb.PC1.bedGraph HiC_CD34_17_3_siCtrl_filtered/HiC_CD34_17_3_siCtrl_filtered.25x50kb.PC1.bedGraph HiC_CD34_18_4_siCtrl_filtered/HiC_CD34_18_4_siCtrl_filtered.25x50kb.PC1.bedGraph HiC_CD34_20_6_siCtrl_Rep1_filtered/HiC_CD34_20_6_siCtrl_Rep1_filtered.25x50kb.PC1.bedGraph HiC_CD34_21_4_siCtrl_Rep1_filtered/HiC_CD34_21_4_siCtrl_Rep1_filtered.25x50kb.PC1.bedGraph HiC_CD34_22_3_siCtrl_filtered/HiC_CD34_22_3_siCtrl_filtered.25x50kb.PC1.bedGraph HiC_CD34_27_4_siCtrl_filtered/HiC_CD34_27_4_siCtrl_filtered.25x50kb.PC1.bedGraph HiC_CD34_28_6_siCtrl_filtered/HiC_CD34_28_6_siCtrl_filtered.25x50kb.PC1.bedGraph HiC_CD34_14_1_SA1_KD_filtered/HiC_CD34_14_1_SA1_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_17_1_SA1_KD_filtered/HiC_CD34_17_1_SA1_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_20_4_SA1_KD_filtered/HiC_CD34_20_4_SA1_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_21_2_SA1_KD_filtered/HiC_CD34_21_2_SA1_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_27_3_SA1_KD_filtered/HiC_CD34_27_3_SA1_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_28_4_SA1_KD_filtered/HiC_CD34_28_4_SA1_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_14_2_SA2_KD_filtered/HiC_CD34_14_2_SA2_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_17_2_SA2_KD_filtered/HiC_CD34_17_2_SA2_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_20_5_SA2_KD_filtered/HiC_CD34_20_5_SA2_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_21_3_SA2_KD_filtered/HiC_CD34_21_3_SA2_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_22_2_SA2_KD_filtered/HiC_CD34_22_2_SA2_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_28_5_SA2_KD_filtered/HiC_CD34_28_5_SA2_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_18_1_RAD21_KD_filtered/HiC_CD34_18_1_RAD21_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_20_1_RAD21_KD_filtered/HiC_CD34_20_1_RAD21_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_22_1_RAD21_KD_filtered/HiC_CD34_22_1_RAD21_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_27_1_RAD21_KD_filtered/HiC_CD34_27_1_RAD21_KD_filtered.25x50kb.PC1.bedGraph HiC_CD34_28_1_RAD21_KD_filtered/HiC_CD34_28_1_RAD21_KD_filtered.25x50kb.PC1.bedGraph'
###IMPORTANT: use "-noblanks"  option otherwise error for PCA and tSNE plots
cd ${TAGHICDIR_CD34}
annotatePeaks.pl HiC_CD34_14_1_SA1_KD_filtered/HiC_CD34_14_1_SA1_KD_filtered.25x50kb.PC1.txt hg38 -size given -bedGraph ${bedgraphs25x50} -noblanks > ${WORKDIR_CD34}/HiC_CD34_PC1_50KB_annotated.txt #input file for PC1 clustering analysis
tail -n +2 "${WORKDIR_CD34}/HiC_CD34_PC1_50KB_annotated.txt" | cut -f1,20-44 > ${WORKDIR_CD34}/tmp.1.txt
echo $'ID\t14_siCtrl\t17_siCtrl\t18_siCtrl\t20_siCtrl\t21_siCtrl_rep\t22_siCtrl\t27_siCtrl\t28_siCtrl\t14_SA1\t17_SA1\t20_SA1\t21_SA1\t27_SA1\t28_SA1\t14_SA2\t17_SA2\t20_SA2\t21_SA2\t22_SA2\t28_SA2\t18_RAD21\t20_RAD21\t22_RAD21\t27_RAD21\t28_RAD21' | cat - ${WORKDIR_CD34}/tmp.1.txt > ${WORKDIR_CD34}/HiC_CD34_PC1_50KB_annotated.txt
#------------------------------------------------------------------------------------

# Analysis for statistically significant changes of PC1 values
##RAD21 KD vs CTRL
cd ${TAGHICDIR_CD34}
annotatePeaks.pl HiC_CD34_14_1_SA1_KD_filtered/HiC_CD34_14_1_SA1_KD_filtered.25x50kb.PC1.txt hg38 -noblanks -bedGraph \
	${bedgraphs25x50_ctrl} ${bedgraphs25x50KB50_RAD21KD}\
    > ${PCADIR}/CTRLvsRAD21.50KB.PC1.ann.txt

getDiffExpression.pl ${PCADIR}/CTRLvsRAD21.50KB.PC1.ann.txt \
CTRL CTRL CTRL CTRL CTRL CTRL CTRL CTRL RAD21KD RAD21KD RAD21KD RAD21KD RAD21KD \
-pc1 -export RAD21KDvsCTRL.50KB.PC1.regions > RAD21KDvsCTRL.PC1.50KB.diff.txt  

###filter with defiened FDR (Column 6) and logFC (col 5) requirement:
####increased compartmentalization
awk 'BEGIN {FS=OFS="\t"} {print $1, $2, $3, $6, $33, $35}' RAD21KDvsCTRL.PC1.50KB.diff.txt  > tmp.PC1.txt
myFilterFile.pl tmp.PC1.txt -column 6 -upperlimit 0.05 > tmp.PC12.txt
myFilterFile.pl tmp.PC12.txt -column 5 -lowerlimit 0.585 > tmp.PC13.txt 
pos2bed.pl tmp.PC13.txt > tmp.PC1UP.bed
wc -l tmp.PC1UP.bed #479 --> "UP" PC1
####decreased compartmentalization
awk 'BEGIN {FS=OFS="\t"} {print $1, $2, $3, $6, $33, $35}' RAD21KDvsCTRL.PC1.50KB.diff.txt  > tmp.PC1.txt
myFilterFile.pl tmp.PC1.txt -column 6 -upperlimit 0.05 > tmp.PC12.txt
myFilterFile.pl tmp.PC12.txt -column 5 -upperlimit -0.585 > tmp.PC13.txt 
pos2bed.pl tmp.PC13.txt > tmp.PC1Down.bed
wc -l tmp.PC1Down.bed #117 --> "Down" PC1

##SA2 KD vs CTRL
cd ${TAGHICDIR_CD34}
annotatePeaks.pl HiC_CD34_14_1_SA1_KD_filtered/HiC_CD34_14_1_SA1_KD_filtered.25x50kb.PC1.txt hg38 -noblanks -bedGraph \
	${bedgraphs25x50_ctrl} ${bedgraphs25x50_SA2KD}\
    > ${PCADIR}/CTRLvsSA2.50KB.PC1.ann.txt

getDiffExpression.pl ${PCADIR}/CTRLvsSA2.50KB.PC1.ann.txt \
CTRL CTRL CTRL CTRL CTRL CTRL CTRL CTRL SA2KD SA2KD SA2KD SA2KD SA2KD SA2KD \
-pc1 -export SA2KDvsCTRL.50KB.PC1.regions > SA2KDvsCTRL.PC1.50KB.diff.txt  

###filter with defiened FDR (Column 6) and logFC (col 5) requirement:
####increased compartmentalization
awk 'BEGIN {FS=OFS="\t"} {print $1, $2, $3, $6, $34, $36}' SA2KDvsCTRL.PC1.50KB.diff.txt   > tmp.PC1.txt
myFilterFile.pl tmp.PC1.txt -column 6 -upperlimit 0.05 > tmp.PC12.txt
myFilterFile.pl tmp.PC12.txt -column 5 -lowerlimit 0.585 > tmp.PC13.txt 
pos2bed.pl tmp.PC13.txt > tmp.PC1UP.bed
wc -l tmp.PC1UP.bed #0 --> "UP" PC1
####decreased compartmentalization
awk 'BEGIN {FS=OFS="\t"} {print $1, $2, $3, $6, $34, $36}' SA2KDvsCTRL.PC1.50KB.diff.txt   > tmp.PC1.txt
myFilterFile.pl tmp.PC1.txt -column 6 -upperlimit 0.05 > tmp.PC12.txt
myFilterFile.pl tmp.PC12.txt -column 5 -upperlimit -0.585 > tmp.PC13.txt 
pos2bed.pl tmp.PC13.txt > tmp.PC1Down.bed
wc -l tmp.PC1Down.bed #0 --> "Down" PC1

###SA1 KD vs CTRL
cd ${TAGHICDIR_CD34}
annotatePeaks.pl HiC_CD34_14_1_SA1_KD_filtered/HiC_CD34_14_1_SA1_KD_filtered.25x50kb.PC1.txt hg38 -noblanks -bedGraph \
	${bedgraphs25x50_ctrl} ${bedgraphs25x50_SA1KD}\
    > ${PCADIR}/CTRLvsSA1.50KB.PC1.ann.txt
cd ${PCADIR}
getDiffExpression.pl ${PCADIR}/CTRLvsSA1.50KB.PC1.ann.txt \
CTRL CTRL CTRL CTRL CTRL CTRL CTRL CTRL SA1KD SA1KD SA1KD SA1KD SA1KD SA1KD \
-pc1 -export SA1KDvsCTRL.50KB.PC1.regions > SA1KDvsCTRL.PC1.50KB.diff.txt
###filter with defiened FDR (Column 6) and logFC (col 5) requirement:
####increased compartmentalization
awk 'BEGIN {FS=OFS="\t"} {print $1, $2, $3, $6, $34, $36}' SA1KDvsCTRL.PC1.50KB.diff.txt   > tmp.PC1.txt
myFilterFile.pl tmp.PC1.txt -column 6 -upperlimit 0.05 > tmp.PC12.txt
myFilterFile.pl tmp.PC12.txt -column 5 -lowerlimit 0.585 > tmp.PC13.txt 
pos2bed.pl tmp.PC13.txt > tmp.PC1UP.bed
wc -l tmp.PC1UP.bed #0 --> "UP" PC1
####decreased compartmentalization
awk 'BEGIN {FS=OFS="\t"} {print $1, $2, $3, $6, $34, $36}' SA1KDvsCTRL.PC1.50KB.diff.txt   > tmp.PC1.txt
myFilterFile.pl tmp.PC1.txt -column 6 -upperlimit 0.05 > tmp.PC12.txt
myFilterFile.pl tmp.PC12.txt -column 5 -upperlimit -0.585 > tmp.PC13.txt 
pos2bed.pl tmp.PC13.txt > tmp.PC1Down.bed
wc -l tmp.PC1Down.bed #0 --> "Down" PC1
#------------------------------------------------------------------------------------
# generate delta PC1 tracks: KD vs ctrl
## separately - mathced by donor
Donor14='14_1_SA1_KD 14_2_SA2_KD'
#vs: 14_3_siCtrl
Donor17="17_1_SA1_KD 17_2_SA2_KD"
#vs: 17_3_siCtrl
Donor18="18_1_RAD21_KD"
#vs: 18_4_siCtrl
Donor20="20_1_RAD21_KD 20_4_SA1_KD 20_5_SA2_KD"
#vs: 20_6_siCtrl_Rep1
Donor21="21_2_SA1_KD 21_3_SA2_KD"
#vs: 21_4_siCtrl_Rep1
Donor22="22_1_RAD21_KD 22_2_SA2_KD"
#vs 22_3_siCtrl
Donor27="27_1_RAD21_KD 17_2_SA2_KD 27_3_SA1_KD"
#vs: 27_4_siCtrl
Donor28="28_1_RAD21_KD 28_4_SA1_KD 28_5_SA2_KD"
#vs: 28_6_siCtrl

##subtracting Control from KD using the coverage COV.bedGraph generated in the compaction stats analysis 
cd ${TAGHICDIR_CD34}
for sn in ${Donor14}
do                      
	subtractBedGraphs.pl HiC_CD34_${sn}_filtered/HiC_CD34_${sn}_filtered.25x50kb.PC1.bedGraph HiC_CD34_14_3_siCtrl_filtered/HiC_CD34_14_3_siCtrl_filtered.25x50kb.PC1.bedGraph -cov ${COMPACTDIR}/HiC_CD34_14_3_siCtrl.COV.bedGraph -name ${sn}vsCTRL -center > ${PCADIR}/${sn}vsCTRL.deltaPC1.bedGraph  
    bedGraphToBigWig ${PCADIR}/${sn}vsCTRL.deltaPC1.bedGraph ${CHROMSIZES_HG38} ${PCADIR}/${sn}vsCTRL.deltaPC1.bigWig
done
for sn in ${Donor17}
do                      
	subtractBedGraphs.pl HiC_CD34_${sn}_filtered/HiC_CD34_${sn}_filtered.25x50kb.PC1.bedGraph HiC_CD34_17_3_siCtrl_filtered/HiC_CD34_17_3_siCtrl_filtered.25x50kb.PC1.bedGraph -cov ${COMPACTDIR}/HiC_CD34_17_3_siCtrl.COV.bedGraph -name ${sn}vsCTRL -center > ${PCADIR}/${sn}vsCTRL.deltaPC1.bedGraph  
    bedGraphToBigWig ${PCADIR}/${sn}vsCTRL.deltaPC1.bedGraph ${CHROMSIZES_HG38} ${PCADIR}/${sn}vsCTRL.deltaPC1.bigWig
done
for sn in ${Donor18}
do                      
	subtractBedGraphs.pl HiC_CD34_${sn}_filtered/HiC_CD34_${sn}_filtered.25x50kb.PC1.bedGraph HiC_CD34_18_4_siCtrl_filtered/HiC_CD34_18_4_siCtrl_filtered.25x50kb.PC1.bedGraph -cov ${COMPACTDIR}/HiC_CD34_18_4_siCtrl.COV.bedGraph -name ${sn}vsCTRL -center > ${PCADIR}/${sn}vsCTRL.deltaPC1.bedGraph  
    bedGraphToBigWig ${PCADIR}/${sn}vsCTRL.deltaPC1.bedGraph ${CHROMSIZES_HG38} ${PCADIR}/${sn}vsCTRL.deltaPC1.bigWig
done
for sn in ${Donor20}
do                      
	subtractBedGraphs.pl HiC_CD34_${sn}_filtered/HiC_CD34_${sn}_filtered.25x50kb.PC1.bedGraph HiC_CD34_20_6_siCtrl_Rep1_filtered/HiC_CD34_20_6_siCtrl_Rep1_filtered.25x50kb.PC1.bedGraph -cov ${COMPACTDIR}/HiC_CD34_20_6_siCtrl_Rep1.COV.bedGraph -name ${sn}vsCTRL -center > ${PCADIR}/${sn}vsCTRL.deltaPC1.bedGraph  
    bedGraphToBigWig ${PCADIR}/${sn}vsCTRL.deltaPC1.bedGraph ${CHROMSIZES_HG38} ${PCADIR}/${sn}vsCTRL.deltaPC1.bigWig
done
for sn in ${Donor21}
do                      
	subtractBedGraphs.pl HiC_CD34_${sn}_filtered/HiC_CD34_${sn}_filtered.25x50kb.PC1.bedGraph HiC_CD34_21_4_siCtrl_Rep1_filtered/HiC_CD34_21_4_siCtrl_Rep1_filtered.25x50kb.PC1.bedGraph -cov ${COMPACTDIR}/HiC_CD34_21_4_siCtrl_Rep1.COV.bedGraph -name ${sn}vsCTRL -center > ${PCADIR}/${sn}vsCTRL.deltaPC1.bedGraph  
    bedGraphToBigWig ${PCADIR}/${sn}vsCTRL.deltaPC1.bedGraph ${CHROMSIZES_HG38} ${PCADIR}/${sn}vsCTRL.deltaPC1.bigWig
done
for sn in ${Donor22}
do                      
	subtractBedGraphs.pl HiC_CD34_${sn}_filtered/HiC_CD34_${sn}_filtered.25x50kb.PC1.bedGraph HiC_CD34_22_3_siCtrl_filtered/HiC_CD34_22_3_siCtrl_filtered.25x50kb.PC1.bedGraph -cov ${COMPACTDIR}/HiC_CD34_22_3_siCtrl.COV.bedGraph -name ${sn}vsCTRL -center > ${PCADIR}/${sn}vsCTRL.deltaPC1.bedGraph  
    bedGraphToBigWig ${PCADIR}/${sn}vsCTRL.deltaPC1.bedGraph ${CHROMSIZES_HG38} ${PCADIR}/${sn}vsCTRL.deltaPC1.bigWig
done
for sn in ${Donor27}
do                      
	subtractBedGraphs.pl HiC_CD34_${sn}_filtered/HiC_CD34_${sn}_filtered.25x50kb.PC1.bedGraph HiC_CD34_27_4_siCtrl_filtered/HiC_CD34_27_4_siCtrl_filtered.25x50kb.PC1.bedGraph -cov ${COMPACTDIR}/HiC_CD34_27_4_siCtrl.COV.bedGraph -name ${sn}vsCTRL -center > ${PCADIR}/${sn}vsCTRL.deltaPC1.bedGraph  
    bedGraphToBigWig ${PCADIR}/${sn}vsCTRL.deltaPC1.bedGraph ${CHROMSIZES_HG38} ${PCADIR}/${sn}vsCTRL.deltaPC1.bigWig
done
for sn in ${Donor28}
do                      
	subtractBedGraphs.pl HiC_CD34_${sn}_filtered/HiC_CD34_${sn}_filtered.25x50kb.PC1.bedGraph HiC_CD34_28_6_siCtrl_filtered/HiC_CD34_28_6_siCtrl_filtered.25x50kb.PC1.bedGraph -cov ${COMPACTDIR}/HiC_CD34_28_6_siCtrl.COV.bedGraph -name ${sn}vsCTRL -center > ${PCADIR}/${sn}vsCTRL.deltaPC1.bedGraph  
    bedGraphToBigWig ${PCADIR}/${sn}vsCTRL.deltaPC1.bedGraph ${CHROMSIZES_HG38} ${PCADIR}/${sn}vsCTRL.deltaPC1.bigWig
done
#------------------------------------------------------------------------------------
#average delta PC1 tracks by condition: use the donor-matched delta PC1 tracks generated above and average them:
cd ${PCADIR}
##SA1KDvsCTRL
$BEDTOOLS unionbedg -i 14_1_SA1_KDvsCTRL.deltaPC1.bedGraph 17_1_SA1_KDvsCTRL.deltaPC1.bedGraph 20_4_SA1_KDvsCTRL.deltaPC1.bedGraph 21_2_SA1_KDvsCTRL.deltaPC1.bedGraph 27_3_SA1_KDvsCTRL.deltaPC1.bedGraph 28_4_SA1_KDvsCTRL.deltaPC1.bedGraph | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' > ave.SA1KDvsCTRL.deltaPC1.bedGraph      
bedGraphToBigWig  ave.SA1KDvsCTRL.deltaPC1.bedGraph ${CHROMSIZES_HG38} ave.SA1KDvsCTRL.deltaPC1.bigWig
##SA2KDvsCTRL
$BEDTOOLS unionbedg -i 14_2_SA2_KDvsCTRL.deltaPC1.bedGraph 17_2_SA2_KDvsCTRL.deltaPC1.bedGraph 20_5_SA2_KDvsCTRL.deltaPC1.bedGraph 21_3_SA2_KDvsCTRL.deltaPC1.bedGraph 22_2_SA2_KDvsCTRL.deltaPC1.bedGraph 28_5_SA2_KDvsCTRL.deltaPC1.bedGraph | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' > ave.SA2KDvsCTRL.deltaPC1.bedGraph      
bedGraphToBigWig  ave.SA2KDvsCTRL.deltaPC1.bedGraph ${CHROMSIZES_HG38} ave.SA2KDvsCTRL.deltaPC1.bigWig
##RAD21KDvsCTRL
$BEDTOOLS unionbedg -i 18_1_RAD21_KDvsCTRL.deltaPC1.bedGraph 20_1_RAD21_KDvsCTRL.deltaPC1.bedGraph 22_1_RAD21_KDvsCTRL.deltaPC1.bedGraph 27_1_RAD21_KDvsCTRL.deltaPC1.bedGraph 28_1_RAD21_KDvsCTRL.deltaPC1.bedGraph | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' > ave.RAD21KDvsCTRL.deltaPC1.bedGraph      
bedGraphToBigWig  ave.RAD21KDvsCTRL.deltaPC1.bedGraph ${CHROMSIZES_HG38} ave.RAD21KDvsCTRL.deltaPC1.bigWig
#------------------------------------------------------------------------------------




				    ###############################################
			    	#                                             #
                    #      combined filtered tagDirectories       #
					#                                             #
                    ###############################################
#------------------------------------------------------------------------------------
#create merged TAGDIRS by codition, run PCA and find compartments
HIC_CTRL_Filt='HiC_CD34_14_3_siCtrl_filtered HiC_CD34_17_3_siCtrl_filtered HiC_CD34_18_4_siCtrl_filtered HiC_CD34_20_6_siCtrl_Rep1_filtered HiC_CD34_21_4_siCtrl_Rep1_filtered HiC_CD34_22_3_siCtrl_filtered HiC_CD34_27_4_siCtrl_filtered HiC_CD34_28_6_siCtrl_filtered ' 
HIC_SA1KD_Filt='HiC_CD34_14_1_SA1_KD_filtered  HiC_CD34_17_1_SA1_KD_filtered  HiC_CD34_20_4_SA1_KD_filtered  HiC_CD34_21_2_SA1_KD_filtered  HiC_CD34_27_3_SA1_KD_filtered  HiC_CD34_28_4_SA1_KD_filtered  '
HIC_SA2KD_Filt='HiC_CD34_14_2_SA2_KD_filtered  HiC_CD34_17_2_SA2_KD_filtered  HiC_CD34_20_5_SA2_KD_filtered  HiC_CD34_21_3_SA2_KD_filtered  HiC_CD34_22_2_SA2_KD_filtered  HiC_CD34_28_5_SA2_KD_filtered  '
HIC_RAD21KD_Filt='HiC_CD34_18_1_RAD21_KD_filtered  HiC_CD34_20_1_RAD21_KD_filtered  HiC_CD34_22_1_RAD21_KD_filtered  HiC_CD34_27_1_RAD21_KD_filtered  HiC_CD34_28_1_RAD21_KD_filtered  '

cd ${TAGHICDIR_CD34}
makeTagDirectory HiC_CD34_siCtrl_combined -d ${HIC_CTRL_Filt}
runHiCpca.pl auto HiC_CD34_siCtrl_combined -res 25000 -window 100000 -genome hg38 -cpu 12
findHiCCompartments.pl HiC_CD34_siCtrl_combined/HiC_CD34_siCtrl_combined.25x100kb.PC1.txt > HiC_CD34_siCtrl_combined/HiC_CD34_siCtrl_combined_compartments.txt

makeTagDirectory HiC_SA1KD_combined -d ${HIC_SA1KD_Filt}
runHiCpca.pl auto HiC_SA1KD_combined -res 25000 -window 100000 -genome hg38 -cpu 12
findHiCCompartments.pl HiC_SA1KD_combined/HiC_SA1KD_combined.25x100kb.PC1.txt > HiC_SA1KD_combined/HiC_SA1KD_combined_compartments.txt

makeTagDirectory HiC_SA2KD_combined -d ${HIC_SA2KD_Filt}
runHiCpca.pl auto HiC_SA2KD_combined -res 25000 -window 100000 -genome hg38 -cpu 12
findHiCCompartments.pl HiC_SA2KD_combined/HiC_SA2KD_combined.25x100kb.PC1.txt > HiC_SA2KD_combined/HiC_SA2KD_combined_compartments.txt

makeTagDirectory HiC_RAD21KD_combined -d ${HIC_RAD21KD_Filt}
runHiCpca.pl auto HiC_RAD21KD_combined -res 25000 -window 100000 -genome hg38 -cpu 12
findHiCCompartments.pl HiC_RAD21KD_combined/HiC_RAD21KD_combined.25x100kb.PC1.txt > HiC_RAD21KD_combined/HiC_RAD21KD_combined_compartments.txt

#------------------------------------------------------------------------------------
#look at QC-stats of combined tagdirs:
HIC_CD34_combined="HiC_CD34_siCtrl_combined HiC_SA1KD_combined HiC_SA2KD_combined HiC_RAD21KD_combined"

cd ${TAGHICDIR_CD34}
for SAMPLE in ${HIC_CD34_combined}
do
paste <(echo -e "# ${SAMPLE}") <(grep -w "genome" ${SAMPLE}/tagInfo.txt) <(grep "localInteractionFraction=" ${SAMPLE}/tagInfo.txt) <(grep "interChrInteractionFraction=" ${SAMPLE}/tagInfo.txt) <(grep "tagsPerBP=" ${SAMPLE}/tagInfo.txt)
done
#------------------------------------------------------------------------------------



                      ##########################################
                      #                                        #
                      #         Finding TADs and Loops         #
                      #   (merged directories, merged donors)  # 
                      #                                        #
                      ##########################################

# used two resolutions for loop and tad finding and merged 2D files
#------------------------------------------------------------------------------------

# res 2500 window 10000  
declare -a SAMPLES=("HiC_CD34_siCtrl_combined" "HiC_SA1KD_combined" "HiC_SA2KD_combined" "HiC_RAD21KD_combined")                    
cd ${TAGHICDIR_CD34}
for SAMPLE in "${SAMPLES[@]}";do
findTADsAndLoops.pl find ${SAMPLE}/ -cpu 34 -res 2500 -window 10000 -genome hg38 \
-p ${BAD_HG38} -o ${LOOPDIR_CD34}/${SAMPLE}.2.5K.10K
done
# res 3000 window 15000  
declare -a SAMPLES=("HiC_CD34_siCtrl_combined" "HiC_SA1KD_combined" "HiC_SA2KD_combined" "HiC_RAD21KD_combined")                     
cd ${TAGHICDIR_CD34}
for SAMPLE in "${SAMPLES[@]}";do
findTADsAndLoops.pl find ${SAMPLE}/ -cpu 24 -res 3000 -window 15000 -genome hg38 \
-p ${BAD_HG38} -o ${LOOPDIR_CD34}/${SAMPLE}.3K.15K
done
#------------------------------------------------------------------------------------
# merge file and create links file 
declare -a SAMPLES=("HiC_CD34_siCtrl_combined" "HiC_SA1KD_combined" "HiC_SA2KD_combined" "HiC_RAD21KD_combined")
cd ${LOOPDIR_CD34}
for SAMPLE in "${SAMPLES[@]}";do
merge2Dbed.pl ${SAMPLE}.3K.15K.loop.2D.bed ${SAMPLE}.2.5K.10K.loop.2D.bed > ${SAMPLE}.merged.loop.2D.bed
merge2Dbed.pl ${SAMPLE}.3K.15K.tad.2D.bed ${SAMPLE}.2.5K.10K.tad.2D.bed > ${SAMPLE}.merged.tad.2D.bed
cut -f1-6,8 ${SAMPLE}.merged.loop.2D.bed | LC_COLLATE=C sort -k1,1 -k2,2n > ${SAMPLE}.merged.links
done


#count TADs and Loops by conditions
declare -a SAMPLES=("HiC_CD34_siCtrl_combined" "HiC_SA1KD_combined" "HiC_SA2KD_combined" "HiC_RAD21KD_combined")
cd ${LOOPDIR_CD34}
for SAMPLE in "${SAMPLES[@]}";do
wc -l ${SAMPLE}.merged.loop.2D.bed
wc -l ${SAMPLE}.merged.tad.2D.bed
wc -l ${SAMPLE}.5K.20K.tad.2D.bed
wc -l ${SAMPLE}.5K.20K.tad.2D.bed
done
#------------------------------------------------------------------------------------



                     #####################################################
                     # Scoring for differential Loops/TAD analysis       #
                     #####################################################

# all conditions compared - used for clustering
#------------------------------------------------------------------------------------
##create Merged loop and TAD sets of all conditions:
####combined 2.5K + 3K res
cd ${LOOPDIR_CD34}
merge2Dbed.pl HiC_CD34_siCtrl_combined.merged.loop.2D.bed HiC_SA1KD_combined.merged.loop.2D.bed HiC_SA2KD_combined.merged.loop.2D.bed HiC_RAD21KD_combined.merged.loop.2D.bed > CD34.merged.loop.2D.bed 
merge2Dbed.pl HiC_CD34_siCtrl_combined.merged.tad.2D.bed HiC_SA1KD_combined.merged.tad.2D.bed HiC_SA2KD_combined.merged.tad.2D.bed HiC_RAD21KD_combined.merged.tad.2D.bed > CD34.merged.tad.2D.bed 

TAGDIRLIST=${HIC_CTRL_Filt}${HIC_SA1KD_Filt}${HIC_SA2KD_Filt}${HIC_RAD21KD_Filt}

## scoring loops and TADs across all samples using 2.5 +3 K res merged loops: HOMER default HiC normalizaiton
cd ${TAGHICDIR_CD34}
findTADsAndLoops.pl score -tad ${LOOPDIR_CD34}/CD34.merged.tad.2D.bed -loop ${LOOPDIR_CD34}/CD34.merged.loop.2D.bed \
-d ${TAGDIRLIST} -cpu 48 -o ${LOOPDIR_CD34}/CD34.merged
## scoring loops and TADs across all samples using 2.5 +3 K res merged loops: HOMER norm2total normalizaiton (using averegae seq. depth)
findTADsAndLoops.pl score -normTotal 200000000 -loop ${LOOPDIR_CD34}/CD34.merged.loop.2D.bed \
-d ${TAGDIRLIST} -cpu 24 -o ${LOOPDIR_CD34}/CD34.normTotalaverage.merged

## Prepare input tables for clustering
cd ${LOOPDIR_CD34}
cut -f1,11-35 CD34.merged.loop.scores.txt | tail -n +2 > tmp.scores.txt
echo $'Gene\t14_siCtrl\t17_siCtrl\t18_siCtrl\t20_siCtrl\t21_siCtrl_rep\t22_siCtrl\t27_siCtrl\t28_siCtrl\t14_SA1\t17_SA1\t20_SA1\t21_SA1\t27_SA1\t28_SA1\t14_SA2\t17_SA2\t20_SA2\t21_SA2\t22_SA2\t28_SA2\t18_RAD21\t20_RAD21\t22_RAD21\t27_RAD21\t28_RAD21' | cat - tmp.scores.txt > tmp.scores2.txt
cp tmp.scores2.txt CD34.merged.loop.scores.Rinput.txt

cut -f1,11-35 CD34.merged.tad.scores.txt | tail -n +2 > tmp.scores.txt
echo $'Gene\t14_siCtrl\t17_siCtrl\t18_siCtrl\t20_siCtrl\t21_siCtrl_rep\t22_siCtrl\t27_siCtrl\t28_siCtrl\t14_SA1\t17_SA1\t20_SA1\t21_SA1\t27_SA1\t28_SA1\t14_SA2\t17_SA2\t20_SA2\t21_SA2\t22_SA2\t28_SA2\t18_RAD21\t20_RAD21\t22_RAD21\t27_RAD21\t28_RAD21' | cat - tmp.scores.txt > tmp.scores2.txt
cp tmp.scores2.txt CD34.merged.tad.scores.Rinput.txt
#------------------------------------------------------------------------------------


# condition vs CTRL sets - used for differential analysis
#------------------------------------------------------------------------------------
cd ${LOOPDIR_CD34}
##create combined loop/TAD sets of CTRL and KDs for 1v1 comparisons
###for CTRL vs SA1KD (#2.5+3KB res)
merge2Dbed.pl HiC_CD34_siCtrl_combined.merged.loop.2D.bed HiC_SA1KD_combined.merged.loop.2D.bed > SA1KDvsCTRL.merged.loop.2D.bed 
merge2Dbed.pl HiC_CD34_siCtrl_combined.merged.tad.2D.bed HiC_SA1KD_combined.merged.tad.2D.bed > SA1KDvsCTRL.merged.tad.2D.bed 
###for CTRL vs SA2KD (#2.5+3KB res)
merge2Dbed.pl HiC_CD34_siCtrl_combined.merged.loop.2D.bed HiC_SA2KD_combined.merged.loop.2D.bed > SA2KDvsCTRL.merged.loop.2D.bed 
merge2Dbed.pl HiC_CD34_siCtrl_combined.merged.tad.2D.bed HiC_SA2KD_combined.merged.tad.2D.bed > SA2KDvsCTRL.merged.tad.2D.bed 
###for CTRL vs RAD21KD (#2.5+3KB res)
merge2Dbed.pl HiC_CD34_siCtrl_combined.merged.loop.2D.bed HiC_RAD21KD_combined.merged.loop.2D.bed > RAD21KDvsCTRL.merged.loop.2D.bed 
merge2Dbed.pl HiC_CD34_siCtrl_combined.merged.tad.2D.bed HiC_RAD21KD_combined.merged.tad.2D.bed > RAD21KDvsCTRL.merged.tad.2D.bed 


##score tads and loops of KD vs siCTRL comparisons: HOMER default and norm2total (average seq depth) normalizations
#------------------------------------------------------------------------------------
cd ${TAGHICDIR_CD34}
###for CTRL vs SA1KD
findTADsAndLoops.pl score -tad ${LOOPDIR_CD34}/SA1KDvsCTRL.merged.tad.2D.bed -loop ${LOOPDIR_CD34}/SA1KDvsCTRL.merged.loop.2D.bed \
-d ${HIC_CTRL_Filt}${HIC_SA1KD_Filt} -cpu 48 -o ${LOOPDIR_CD34}/SA1KDvsCTRL.merged
findTADsAndLoops.pl score -normTotal 200000000 -loop ${LOOPDIR_CD34}/SA1KDvsCTRL.merged.loop.2D.bed \
-d ${HIC_CTRL_Filt}${HIC_SA1KD_Filt} -cpu 12 -o ${LOOPDIR_CD34}/SA1KDvsCTRL.normTotalaverage.merged
###for CTRL vs SA2KD
findTADsAndLoops.pl score -tad ${LOOPDIR_CD34}/SA2KDvsCTRL.merged.tad.2D.bed -loop ${LOOPDIR_CD34}/SA2KDvsCTRL.merged.loop.2D.bed \
-d ${HIC_CTRL_Filt}${HIC_SA2KD_Filt} -cpu 32 -o ${LOOPDIR_CD34}/SA2KDvsCTRL.merged
findTADsAndLoops.pl score -normTotal 200000000 -loop ${LOOPDIR_CD34}/SA2KDvsCTRL.merged.loop.2D.bed \
-d ${HIC_CTRL_Filt}${HIC_SA2KD_Filt} -cpu 12 -o ${LOOPDIR_CD34}/SA2KDvsCTRL.normTotalaverage.merged
###for CTRL vs RAD21KD
findTADsAndLoops.pl score -tad ${LOOPDIR_CD34}/RAD21KDvsCTRL.merged.tad.2D.bed -loop ${LOOPDIR_CD34}/RAD21KDvsCTRL.merged.loop.2D.bed \
-d ${HIC_CTRL_Filt}${HIC_RAD21KD_Filt} -cpu 48 -o ${LOOPDIR_CD34}/RAD21KDvsCTRL.merged
findTADsAndLoops.pl score -normTotal 200000000 -loop ${LOOPDIR_CD34}/RAD21KDvsCTRL.merged.loop.2D.bed \
-d ${HIC_CTRL_Filt}${HIC_RAD21KD_Filt} -cpu 20 -o ${LOOPDIR_CD34}/RAD21KDvsCTRL.normTotalaverage.merged
#------------------------------------------------------------------------------------

##collect total combined scores: using normTotal scores and default normalized scores in one table per comparison
#------------------------------------------------------------------------------------
cd ${LOOPDIR_CD34}
echo -e "Loop\t14_siCtrl\t17_siCtrl\t18_siCtrl\t20_siCtrl\t21_siCtrl\t22_siCtrl\t27_siCtrl\t28_siCtrl\t14_SA1\t17_SA1\t20_SA1\t21_SA1\t27_SA1\t28_SA1" | \
cat - <(awk 'BEGIN {FS=OFS="\t"} {for(i=11;i<=NF;i++) a[i]+=$i} END{l="sum.n2t";i=11;while(i in a){l=l"\t"a[i];i++};print l}' SA1KDvsCTRL.normTotalaverage.merged.loop.scores.txt) | \
cat - <(awk 'BEGIN {FS=OFS="\t"} {for(i=11;i<=NF;i++) a[i]+=$i} END{l="sum";i=11;while(i in a){l=l"\t"a[i];i++};print l}' SA1KDvsCTRL.merged.loop.scores.txt) > ${LOOPDIR_CD34}/SA1KDvsCTRL.totalCounts.table.txt

echo -e "Loop\t14_siCtrl\t17_siCtrl\t18_siCtrl\t20_siCtrl\t21_siCtrl\t22_siCtrl\t27_siCtrl\t28_siCtrl\t14_SA2\t17_SA2\t20_SA2\t21_SA2\t22_SA2\t28_SA2" | \
cat - <(awk 'BEGIN {FS=OFS="\t"} {for(i=11;i<=NF;i++) a[i]+=$i} END{l="sum.n2t";i=11;while(i in a){l=l"\t"a[i];i++};print l}' SA2KDvsCTRL.normTotalaverage.merged.loop.scores.txt) | \
cat - <(awk 'BEGIN {FS=OFS="\t"} {for(i=11;i<=NF;i++) a[i]+=$i} END{l="sum";i=11;while(i in a){l=l"\t"a[i];i++};print l}' SA2KDvsCTRL.merged.loop.scores.txt) > ${LOOPDIR_CD34}/SA2KDvsCTRL.totalCounts.table.txt

echo -e "Loop\t14_siCtrl\t17_siCtrl\t18_siCtrl\t20_siCtrl\t21_siCtrl\t22_siCtrl\t27_siCtrl\t28_siCtrl\t18_RAD21\t20_RAD21\t22_RAD21\t27_RAD21\t28_RAD21" | \
cat - <(awk 'BEGIN {FS=OFS="\t"} {for(i=11;i<=NF;i++) a[i]+=$i} END{l="sum.n2t";i=11;while(i in a){l=l"\t"a[i];i++};print l}' RAD21KDvsCTRL.normTotalaverage.merged.loop.scores.txt) | \
cat - <(awk 'BEGIN {FS=OFS="\t"} {for(i=11;i<=NF;i++) a[i]+=$i} END{l="sum";i=11;while(i in a){l=l"\t"a[i];i++};print l}' RAD21KDvsCTRL.merged.loop.scores.txt) > ${LOOPDIR_CD34}/RAD21KDvsCTRL.totalCounts.table.txt
#------------------------------------------------------------------------------------


##prepare scoring and loop coordinate tables for R input:
#------------------------------------------------------------------------------------
cd ${LOOPDIR_CD34}
cut -f1,11-23 RAD21KDvsCTRL.merged.loop.scores.txt | tail -n +2 > tmp.scores.txt
echo $'Loop\t14_siCtrl\t17_siCtrl\t18_siCtrl\t20_siCtrl\t21_siCtrl\t22_siCtrl\t27_siCtrl\t28_siCtrl\t18_RAD21\t20_RAD21\t22_RAD21\t27_RAD21\t28_RAD21' | cat - tmp.scores.txt > tmp.scores.RAD21KD.txt
cut -f1,2,3,4,5,6,7 RAD21KDvsCTRL.merged.loop.scores.txt | tail -n +2 > tmp.loopstartstop.txt
echo $'Loop\tchr1\tstart1\tend1\tchr2\tstart2\tend2' | cat - tmp.loopstartstop.txt > tmp.loopstartstop.RAD21KD.txt

cut -f1,11-24 SA2KDvsCTRL.merged.loop.scores.txt | tail -n +2 > tmp.scores.txt
echo $'Loop\t14_siCtrl\t17_siCtrl\t18_siCtrl\t20_siCtrl\t21_siCtrl\t22_siCtrl\t27_siCtrl\t28_siCtrl\t14_SA2\t17_SA2\t20_SA2\t21_SA2\t22_SA2\t28_SA2' | cat - tmp.scores.txt > tmp.scores.SA2KD.txt
cut -f1,2,3,4,5,6,7 SA2KDvsCTRL.merged.loop.scores.txt | tail -n +2 > tmp.loopstartstop.txt
echo $'Loop\tchr1\tstart1\tend1\tchr2\tstart2\tend2' | cat - tmp.loopstartstop.txt > tmp.loopstartstop.SA2KD.txt

cut -f1,11-24 SA1KDvsCTRL.merged.loop.scores.txt | tail -n +2 > tmp.scores.txt
echo $'Loop\t14_siCtrl\t17_siCtrl\t18_siCtrl\t20_siCtrl\t21_siCtrl\t22_siCtrl\t27_siCtrl\t28_siCtrl\t14_SA1\t17_SA1\t20_SA1\t21_SA1\t27_SA1\t28_SA1' | cat - tmp.scores.txt > tmp.scores.SA1KD.txt
cut -f1,2,3,4,5,6,7 SA1KDvsCTRL.merged.loop.scores.txt | tail -n +2 > tmp.loopstartstop.txt
echo $'Loop\tchr1\tstart1\tend1\tchr2\tstart2\tend2' | cat - tmp.loopstartstop.txt > tmp.loopstartstop.SA1KD.txt
#------------------------------------------------------------------------------------



##versions without sex-chromosomes (as done for AMLs) for individual KD vs CTRL loop sets
#------------------------------------------------------------------------------------
LOOPDIR_XYrm="${LOOPDIR_CD34}/withoutXY"
mkdir $LOOPDIR_XYrm
cd ${LOOPDIR_CD34}
for KD in ${KDs};do
grep -v 'chrX' ${LOOPDIR_CD34}/${KD}vsCTRL.merged.loop.scores.txt > tmp.scores3.txt
grep -v 'chrY' tmp.scores3.txt > ${KD}vsCTRL.merged.loop.scores.XYrm.txt
grep -v 'chrX' ${LOOPDIR_CD34}/${KD}vsCTRL.normTotalaverage.merged.loop.scores.txt > tmp.scores3.txt
grep -v 'chrY' tmp.scores3.txt > ${LOOPDIR_XYrm}/${KD}vsCTRL.normTotalaverage.merged.loop.scores.XYrm.txt
done
#collect total scores using normTotala scores (Homeroption for scoring) and default normalized scores in one table excluding XY loops)
cd ${LOOPDIR_XYrm}
echo -e "Loop\t14_siCtrl\t17_siCtrl\t18_siCtrl\t20_siCtrl\t21_siCtrl\t22_siCtrl\t27_siCtrl\t28_siCtrl\t14_SA2\t17_SA2\t20_SA2\t21_SA2\t22_SA2\t28_SA2" | \
cat - <(awk 'BEGIN {FS=OFS="\t"} {for(i=11;i<=NF;i++) a[i]+=$i} END{l="sum.n2t";i=11;while(i in a){l=l"\t"a[i];i++};print l}' SA2KDvsCTRL.normTotalaverage.merged.loop.scores.XYrm.txt) | \
cat - <(awk 'BEGIN {FS=OFS="\t"} {for(i=11;i<=NF;i++) a[i]+=$i} END{l="sum";i=11;while(i in a){l=l"\t"a[i];i++};print l}' SA2KDvsCTRL.merged.loop.scores.XYrm.txt) > SA2KDvsCTRL.totalCounts.table.XYrm.txt

echo -e "Loop\t14_siCtrl\t17_siCtrl\t18_siCtrl\t20_siCtrl\t21_siCtrl\t22_siCtrl\t27_siCtrl\t28_siCtrl\t14_SA1\t17_SA1\t20_SA1\t21_SA1\t27_SA1\t28_SA1" | \
cat - <(awk 'BEGIN {FS=OFS="\t"} {for(i=11;i<=NF;i++) a[i]+=$i} END{l="sum.n2t";i=11;while(i in a){l=l"\t"a[i];i++};print l}' SA1KDvsCTRL.normTotalaverage.merged.loop.scores.XYrm.txt) | \
cat - <(awk 'BEGIN {FS=OFS="\t"} {for(i=11;i<=NF;i++) a[i]+=$i} END{l="sum";i=11;while(i in a){l=l"\t"a[i];i++};print l}' SA1KDvsCTRL.merged.loop.scores.XYrm.txt) > SA1KDvsCTRL.totalCounts.table.XYrm.txt

echo -e "Loop\t14_siCtrl\t17_siCtrl\t18_siCtrl\t20_siCtrl\t21_siCtrl\t22_siCtrl\t27_siCtrl\t28_siCtrl\t18_RAD21\t20_RAD21\t22_RAD21\t27_RAD21\t28_RAD21" | \
cat - <(awk 'BEGIN {FS=OFS="\t"} {for(i=11;i<=NF;i++) a[i]+=$i} END{l="sum.n2t";i=11;while(i in a){l=l"\t"a[i];i++};print l}' RAD21KDvsCTRL.normTotalaverage.merged.loop.scores.XYrm.txt) | \
cat - <(awk 'BEGIN {FS=OFS="\t"} {for(i=11;i<=NF;i++) a[i]+=$i} END{l="sum";i=11;while(i in a){l=l"\t"a[i];i++};print l}' RAD21KDvsCTRL.merged.loop.scores.XYrm.txt) > RAD21KDvsCTRL.totalCounts.table.XYrm.txt
##prepare scoring and loop coordinate tables for R input excludigin XY loops
cd ${LOOPDIR_XYrm}
cut -f1,11-23 RAD21KDvsCTRL.merged.loop.scores.XYrm.txt | tail -n +2 > tmp.scores.txt
echo $'Loop\t14_siCtrl\t17_siCtrl\t18_siCtrl\t20_siCtrl\t21_siCtrl\t22_siCtrl\t27_siCtrl\t28_siCtrl\t18_RAD21\t20_RAD21\t22_RAD21\t27_RAD21\t28_RAD21' | cat - tmp.scores.txt > tmp.scores.RAD21KD.XYrm.txt
cut -f1,2,3,4,5,6,7 RAD21KDvsCTRL.merged.loop.scores.txt | tail -n +2 > tmp.loopstartstop.txt
echo $'Loop\tchr1\tstart1\tend1\tchr2\tstart2\tend2' | cat - tmp.loopstartstop.txt > tmp.loopstartstop.RAD21KD.XYrm.txt

cut -f1,11-24 SA2KDvsCTRL.merged.loop.scores.XYrm.txt | tail -n +2 > tmp.scores.txt
echo $'Loop\t14_siCtrl\t17_siCtrl\t18_siCtrl\t20_siCtrl\t21_siCtrl\t22_siCtrl\t27_siCtrl\t28_siCtrl\t14_SA2\t17_SA2\t20_SA2\t21_SA2\t22_SA2\t28_SA2' | cat - tmp.scores.txt > tmp.scores.SA2KD.XYrm.txt
cut -f1,2,3,4,5,6,7 SA2KDvsCTRL.merged.loop.scores.txt | tail -n +2 > tmp.loopstartstop.txt
echo $'Loop\tchr1\tstart1\tend1\tchr2\tstart2\tend2' | cat - tmp.loopstartstop.txt > tmp.loopstartstop.SA2KD.XYrm.txt

cut -f1,11-24 SA1KDvsCTRL.merged.loop.scores.XYrm.txt | tail -n +2 > tmp.scores.txt
echo $'Loop\t14_siCtrl\t17_siCtrl\t18_siCtrl\t20_siCtrl\t21_siCtrl\t22_siCtrl\t27_siCtrl\t28_siCtrl\t14_SA1\t17_SA1\t20_SA1\t21_SA1\t27_SA1\t28_SA1' | cat - tmp.scores.txt > tmp.scores.SA1KD.XYrm.txt
cut -f1,2,3,4,5,6,7 SA1KDvsCTRL.merged.loop.scores.txt | tail -n +2 > tmp.loopstartstop.txt
echo $'Loop\tchr1\tstart1\tend1\tchr2\tstart2\tend2' | cat - tmp.loopstartstop.txt > tmp.loopstartstop.SA1KD.XYrm.txt
#------------------------------------------------------------------------------------




                     #####################################################
                     #       Differential TAD Analysis using Homer       #
                     #####################################################

#Test differential enrichement of TADs in 1vs1 comparisons (#TADs at 2.5+3KB res)
#-------------------------------------------------------------------------------------------------------------------------------------------------
cd ${LOOPDIR_CD34} 
###SA1 KD vs CTRL
getDiffExpression.pl SA1KDvsCTRL.merged.tad.scores.txt CTRL CTRL CTRL CTRL CTRL CTRL CTRL CTRL SA1KD SA1KD SA1KD SA1KD SA1KD SA1KD \
> SA1KDvsCTRL.merged.diff.tad.txt
###SA2 KD vs CTRL
getDiffExpression.pl SA2KDvsCTRL.merged.tad.scores.txt CTRL CTRL CTRL CTRL CTRL CTRL CTRL CTRL SA2KD SA2KD SA2KD SA2KD SA2KD SA2KD \
> SA2KDvsCTRL.merged.diff.tad.txt
###RAD21 KD vs CTRL
cd ${LOOPDIR_CD34} #no batch correction possible; limma!
getDiffExpression.pl RAD21KDvsCTRL.merged.tad.scores.txt CTRL CTRL CTRL CTRL CTRL CTRL CTRL CTRL RAD21KD RAD21KD RAD21KD RAD21KD RAD21KD \
> RAD21KDvsCTRL.merged.diff.tad.txt
#--------------------------------------------------------------------------------------------------------------------------------------------------





               ###########################################################################
               #         Scoring of loops based on SA2mut/CTRL AML loop set              #
 			   ###########################################################################
####scoring of KD vs CTRL datasets based on the loops detected in STAG2 mutant + CTRL AMLs
####this allows for direct statistical comparisons of differential loops of patients in the KD model
mkdir ${LOOPDIR_CD34}/SA2mutLoopScoring
#--------------------------------------------------------------------------------------------------------------------------------------------------
###########SA2 KDs
cd ${TAGHICDIR_CD34}
#regular normalization
findTADsAndLoops.pl score -loop  ${AMLLOOPDIR}/${SA2mutloops} \
-d ${HIC_CTRL_Filt}${HIC_SA2KD_Filt} -cpu 20 -o ${LOOPDIR_CD34}/SA2mutLoopScoring/SA2KDvsCTRL.SA2mutLoops
#normTotal to average seq depth
findTADsAndLoops.pl score -normTotal 200000000 -loop ${AMLLOOPDIR}/${SA2mutloops} \
-d ${HIC_CTRL_Filt}${HIC_SA2KD_Filt} -cpu 20 -o ${LOOPDIR_CD34}/SA2mutLoopScoring/SA2KDvsCTRL.SA2mutLoops.normTotalaverage
###########SA1 KDs
#regular normalization
findTADsAndLoops.pl score -loop  ${AMLLOOPDIR}/${SA2mutloops} \
-d ${HIC_CTRL_Filt}${HIC_SA1KD_Filt} -cpu 12 -o ${LOOPDIR_CD34}/SA2mutLoopScoring/SA1KDvsCTRL.SA2mutLoops
#normTotal to average seq depth
findTADsAndLoops.pl score -normTotal 200000000 -loop ${AMLLOOPDIR}/${SA2mutloops} \
-d ${HIC_CTRL_Filt}${HIC_SA1KD_Filt} -cpu 12 -o ${LOOPDIR_CD34}/SA2mutLoopScoring/SA1KDvsCTRL.SA2mutLoops.normTotalaverage
###########RAD21 KDs
findTADsAndLoops.pl score -loop  ${AMLLOOPDIR}/${SA2mutloops} \
-d ${HIC_CTRL_Filt}${HIC_RAD21KD_Filt} -cpu 12 -o ${LOOPDIR_CD34}/SA2mutLoopScoring/RAD21KDvsCTRL.SA2mutLoops
#normTotal to average seq depth
findTADsAndLoops.pl score -normTotal 200000000 -loop ${AMLLOOPDIR}/${SA2mutloops} \
-d ${HIC_CTRL_Filt}${HIC_RAD21KD_Filt} -cpu 12 -o ${LOOPDIR_CD34}/SA2mutLoopScoring/RAD21KDvsCTRL.SA2mutLoops.normTotalaverage
#--------------------------------------------------------------------------------------------------------------------------------------------------
#remove XY chromosomes
cd ${LOOPDIR_CD34}/SA2mutLoopScoring/
for KD in ${KDs};do
grep -v 'chrX' ${KD}vsCTRL.SA2mutLoops.loop.scores.txt > tmp.scores3.txt
grep -v 'chrY' tmp.scores3.txt > ${KD}vsCTRL.SA2mutLoops.loop.scores.XYrm.txt
grep -v 'chrX' ${KD}vsCTRL.SA2mutLoops.normTotalaverage.loop.scores.txt > tmp.scores3.txt
grep -v 'chrY' tmp.scores3.txt > ${KD}vsCTRL.SA2mutLoops.normTotalaverage.loop.scores.XYrm.txt
done

#collect total scores using normTotala scores (Homeroption for scoring) and default normalized scores without XY in one table
cd ${LOOPDIR_CD34}/SA2mutLoopScoring/
echo -e "Loop\t14_siCtrl\t17_siCtrl\t18_siCtrl\t20_siCtrl\t21_siCtrl\t22_siCtrl\t27_siCtrl\t28_siCtrl\t14_SA1\t17_SA1\t20_SA1\t21_SA1\t27_SA1\t28_SA1" | \
cat - <(awk 'BEGIN {FS=OFS="\t"} {for(i=11;i<=NF;i++) a[i]+=$i} END{l="sum.n2t";i=11;while(i in a){l=l"\t"a[i];i++};print l}' SA1KDvsCTRL.SA2mutLoops.normTotalaverage.loop.scores.XYrm.txt) | \
cat - <(awk 'BEGIN {FS=OFS="\t"} {for(i=11;i<=NF;i++) a[i]+=$i} END{l="sum";i=11;while(i in a){l=l"\t"a[i];i++};print l}' SA1KDvsCTRL.SA2mutLoops.loop.scores.XYrm.txt) > SA1KDvsCTRL.SA2mutLoops.totalCounts.table.txt

echo -e "Loop\t14_siCtrl\t17_siCtrl\t18_siCtrl\t20_siCtrl\t21_siCtrl\t22_siCtrl\t27_siCtrl\t28_siCtrl\t14_SA2\t17_SA2\t20_SA2\t21_SA2\t22_SA2\t28_SA2" | \
cat - <(awk 'BEGIN {FS=OFS="\t"} {for(i=11;i<=NF;i++) a[i]+=$i} END{l="sum.n2t";i=11;while(i in a){l=l"\t"a[i];i++};print l}' SA2KDvsCTRL.SA2mutLoops.normTotalaverage.loop.scores.XYrm.txt) | \
cat - <(awk 'BEGIN {FS=OFS="\t"} {for(i=11;i<=NF;i++) a[i]+=$i} END{l="sum";i=11;while(i in a){l=l"\t"a[i];i++};print l}' SA2KDvsCTRL.SA2mutLoops.loop.scores.XYrm.txt) > SA2KDvsCTRL.SA2mutLoops.totalCounts.table.txt

echo -e "Loop\t14_siCtrl\t17_siCtrl\t18_siCtrl\t20_siCtrl\t21_siCtrl\t22_siCtrl\t27_siCtrl\t28_siCtrl\t18_RAD21\t20_RAD21\t22_RAD21\t27_RAD21\t28_RAD21" | \
cat - <(awk 'BEGIN {FS=OFS="\t"} {for(i=11;i<=NF;i++) a[i]+=$i} END{l="sum.n2t";i=11;while(i in a){l=l"\t"a[i];i++};print l}' RAD21KDvsCTRL.SA2mutLoops.normTotalaverage.loop.scores.XYrm.txt) | \
cat - <(awk 'BEGIN {FS=OFS="\t"} {for(i=11;i<=NF;i++) a[i]+=$i} END{l="sum";i=11;while(i in a){l=l"\t"a[i];i++};print l}' RAD21KDvsCTRL.SA2mutLoops.loop.scores.XYrm.txt) > RAD21KDvsCTRL.SA2mutLoops.totalCounts.table.txt
#--------------------------------------------------------------------------------------------------------------------------------------------------
##remove cols not required for DESEQ
cd ${LOOPDIR_CD34}/SA2mutLoopScoring/
###SA2KD
cut -f1,11-24 SA2KDvsCTRL.SA2mutLoops.loop.scores.XYrm.txt | tail -n +2 > tmp.scores.txt
echo $'Loop\t14_siCtrl\t17_siCtrl\t18_siCtrl\t20_siCtrl\t21_siCtrl\t22_siCtrl\t27_siCtrl\t28_siCtrl\t14_SA2\t17_SA2\t20_SA2\t21_SA2\t22_SA2\t28_SA2' | cat - tmp.scores.txt > tmp.scores.SA2KD.SA2mut.txt
cut -f1,2,3,4,5,6,7 SA2KDvsCTRL.SA2mutLoops.loop.scores.XYrm.txt | tail -n +2 > tmp.loopstartstop.txt
echo $'Loop\tchr1\tstart1\tend1\tchr2\tstart2\tend2' | cat - tmp.loopstartstop.txt > tmp.loopstartstop.SA2KD.SA2mut.txt
###RAD21KD
cut -f1,11-23 RAD21KDvsCTRL.SA2mutLoops.loop.scores.XYrm.txt | tail -n +2 > tmp.scores.txt
echo $'Loop\t14_siCtrl\t17_siCtrl\t18_siCtrl\t20_siCtrl\t21_siCtrl\t22_siCtrl\t27_siCtrl\t28_siCtrl\t18_RAD21\t20_RAD21\t22_RAD21\t27_RAD21\t28_RAD21' | cat - tmp.scores.txt > tmp.scores.RAD21KD.SA2mut.txt
cut -f1,2,3,4,5,6,7 RAD21KDvsCTRL.SA2mutLoops.loop.scores.XYrm.txt | tail -n +2 > tmp.loopstartstop.txt
echo $'Loop\tchr1\tstart1\tend1\tchr2\tstart2\tend2' | cat - tmp.loopstartstop.txt >tmp.loopstartstop.RAD21KD.SA2mut.txt
###SA1 KD
cut -f1,11-24 SA1KDvsCTRL.SA2mutLoops.loop.scores.XYrm.txt | tail -n +2 > tmp.scores.txt
echo $'Loop\t14_siCtrl\t17_siCtrl\t18_siCtrl\t20_siCtrl\t21_siCtrl\t22_siCtrl\t27_siCtrl\t28_siCtrl\t14_SA1\t17_SA1\t20_SA1\t21_SA1\t27_SA1\t28_SA1' | cat - tmp.scores.txt > tmp.scores.SA1KD.SA2mut.txt
cut -f1,2,3,4,5,6,7 SA1KDvsCTRL.SA2mutLoops.loop.scores.XYrm.txt | tail -n +2 > tmp.loopstartstop.txt
echo $'Loop\tchr1\tstart1\tend1\tchr2\tstart2\tend2' | cat - tmp.loopstartstop.txt > tmp.loopstartstop.SA1KD.SA2mut.txt
#--------------------------------------------------------------------------------------------------------------------------------------------------
