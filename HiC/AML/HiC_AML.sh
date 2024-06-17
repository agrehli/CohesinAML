#!/bin/bash
#by Alexander Fischer, November 2021


#################################################################################
#         Analysis of deep sequenced HiC data of Cohesin mutated AML samples    #
#################################################################################
# general paths
DIR_PKG="/misc/software/ngs"
DIR_SOFT="/misc/software"
DIR_PKG="${DIR_SOFT}/ngs"
DIR_DATA="/misc/data"
OS=$(lsb_release -c |grep "^Codename" | awk -F: '{print $2}' | sed 's/[[:blank:]]//g')

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
RAWHICDATA="${DIR_DATA}/rawData/3Dchromatin/HiC/AML"
MAPHICDIR="${DIR_DATA}/processedData/mapping/3Dchromatin/HiC/GRCh38/AML"
TAGHICDIR="${DIR_DATA}/processedData/tagDir/3Dchromatin/HiC/GRCh38/AML"
## Analysis
WORKDIR="${DIR_DATA}/analysis/project_cohesin/Cohesin_AML/HiC"
FIGURESDIR="${WORKDIR}/figures" 
TMPDIR="/loctmp"
PCADIR="${WORKDIR}/PCA"
QCDIR="${WORKDIR}/QC"
LOOPDIR="${WORKDIR}/loops"

# vectors of all samples named by condition and PatientVialIDs (as used for the file names)
PATIENTS_CTRL='ctr_16911 ctr_18136_Rep1 ctr_18519 ctr_19405_Rep1 ctr_19416 ctr_21047_Rep1 ctr_21290 '
PATIENTS_RAD='RAD21_UKR186_Rep1 RAD21_23039 RAD21_26830 RAD21_38455 '
PATIENTS_SA2='SA2_24743 SA2_27396_Rep1 SA2_29728 SA2_9708_Rep2'
PATIENTS_COAML=$PATIENTS_CTRL$PATIENTS_RAD$PATIENTS_SA2

# create new directories for plots
mkdir ${FIGURESDIR}
mkdir ${FIGURESDIR}/PC1_Clustering
mkdir ${FIGURESDIR}/Loop_Clustering
mkdir ${FIGURESDIR}/TAD_Clustering
mkdir ${PCADIR}
mkdir ${QCDIR}
mkdir ${LOOPDIR}
#------------------------------------------------------------------------------------

                         #####################################
                         #         Sequencing QC             #
                         #####################################

#------------------------------------------------------------------------------------
#run FastQC
for sn in ${PATIENTS_COAML};
do
fastqc-0.11.7 ${RAWHICDATA}/${sn}_R1.fastq.gz -o ${RAWHICDATA}/FastQC
fastqc-0.11.7 ${RAWHICDATA}/${sn}_R2.fastq.gz -o ${RAWHICDATA}/FastQC
done
#------------------------------------------------------------------------------------
#summarize using multiqc
multiqc_1.13 \
${RAWHICDATA}/FastQC/HiC_AML_ctr_16911_R1* \
${RAWHICDATA}/FastQC/HiC_AML_ctr_16911_R2* \
${RAWHICDATA}/FastQC/HiC_AML_ctr_18136_Rep1_R* \
${RAWHICDATA}/FastQC/HiC_AML_ctr_18519_R1* \
${RAWHICDATA}/FastQC/HiC_AML_ctr_18519_R2* \
${RAWHICDATA}/FastQC/HiC_AML_ctr_19405_Rep1_R* \
${RAWHICDATA}/FastQC/HiC_AML_ctr_19416_R1* \
${RAWHICDATA}/FastQC/HiC_AML_ctr_19416_R2* \
${RAWHICDATA}/FastQC/HiC_AML_ctr_21047_Rep1_R* \
${RAWHICDATA}/FastQC/HiC_AML_ctr_21290_R1* \
${RAWHICDATA}/FastQC/HiC_AML_ctr_21290_R2* \
${RAWHICDATA}/FastQC/HiC_RAD21_UKR186_Rep1_R* \
${RAWHICDATA}/FastQC/HiC_AML_RAD21_23039_R1* \
${RAWHICDATA}/FastQC/HiC_AML_RAD21_23039_R2* \
${RAWHICDATA}/FastQC/HiC_AML_RAD21_26830_R1* \
${RAWHICDATA}/FastQC/HiC_AML_RAD21_26830_R2* \
${RAWHICDATA}/FastQC/HiC_AML_RAD21_38455_R1* \
${RAWHICDATA}/FastQC/HiC_AML_RAD21_38455_R2* \
${RAWHICDATA}/FastQC/HiC_AML_SA2_24743_R1* \
${RAWHICDATA}/FastQC/HiC_AML_SA2_24743_R2* \
${RAWHICDATA}/FastQC/HiC_AML_SA2_27396_Rep1_R* \
${RAWHICDATA}/FastQC/HiC_AML_SA2_29728_R1* \
${RAWHICDATA}/FastQC/HiC_AML_SA2_29728_R2* \
${RAWHICDATA}/FastQC/HiC_AML_SA2_9708_Rep2_R* \
-n HiC_deepseq_Cohesin_AML_multiQC -f -o ${RAWHICDATA}/FastQC/
#------------------------------------------------------------------------------------

                         ##################################################
                         #    Processing and Mapping of the raw data      #
                         ##################################################

#------------------------------------------------------------------------------------ 
#pre-processing: trimming, mapping, tagdir, tagdir filtering
#check if all filenames are correct
for sn in ${PATIENTS_COAML};
do
if test -f "${RAWHICDATA}/${sn}_R1.fastq.gz"; then
    echo "${RAWHICDATA}/${sn}_R1.fastq.gz exists."
fi
done
#------------------------------------------------------------------------------------
for sn in ${PATIENTS_COAML};
do
#trimming R1
time skewer --threads 8 -x GATC -m tail  --min 20 ${RAWHICDATA}/${sn}_R1.fastq.gz --compress --output ${RAWHICDATA}/${sn}_R1.fastq.gz.skewer  &>${RAWHICDATA}/skewer_${sn}_R1.log
#trimming R2
time skewer --threads 8 -x GATC -m tail  --min 20 ${RAWHICDATA}/${sn}_R2.fastq.gz --compress --output ${RAWHICDATA}/${sn}_R2.fastq.gz.skewer  &>${RAWHICDATA}/skewer_${sn}_R2.log
#mapping R1 using bowtie
${BOWTIE2} -x GRCh38.PRI_p10.refOnly/hg38 -p 8 ${RAWHICDATA}/${sn}_R1.fastq.gz.skewer-trimmed.fastq.gz >${MAPHICDIR}/${sn}_R1.sam 2>${MAPHICDIR}/${sn}_R1.log
#mapping R2 using bowtie
${BOWTIE2} -x GRCh38.PRI_p10.refOnly/hg38 -p 8 ${RAWHICDATA}/${sn}_R2.fastq.gz.skewer-trimmed.fastq.gz >${MAPHICDIR}/${sn}_R2.sam 2>${MAPHICDIR}/${sn}_R2.log
#tagDir creation
makeTagDirectory ${TAGHICDIR}/${sn} ${MAPHICDIR}/${sn}_R1.sam,${MAPHICDIR}/${sn}_R2.sam -tbp 1 -genome hg38 -checkGC
cp -r ${TAGHICDIR}/${sn} ${TAGHICDIR}/${sn}_filtered
#tagDir filtering
makeTagDirectory ${TAGHICDIR}/${sn}_filtered -update -restrictionSite GATC -both -genome hg38 -removePEbg -removeSelfLigation removeSpikes 10000 5
#cleanup
echo "cleaning up for ${sn}"
samtools view -S -@12 -b ${MAPHICDIR}/${sn}_R1.sam > ${MAPHICDIR}/${sn}_R1.bam
samtools view -S -@12 -b ${MAPHICDIR}/${sn}_R2.sam > ${MAPHICDIR}/${sn}_R2.bam
samtools index ${MAPHICDIR}/${sn}_R1.bam
samtools index ${MAPHICDIR}/${sn}_R2.bam
rm ${MAPHICDIR}/${sn}_R1.sam ${MAPHICDIR}/${sn}_R2.sam
rm ${RAWHICDATA}/${sn}*.fastq.gz.skewer-trimmed.fastq.gz
mv ${TAGHICDIR}/${sn}/tagInfo.txt ${WORKDIR}/taginfo_unfiltered/${sn}_tagInfo.txt
rm -r ${TAGHICDIR}/${sn}
echo "finished with ${sn}"
done
#------------------------------------------------------------------------------------
                         ##################################################
                         #         QC for HOMER TAG directories           #
                         ##################################################

#------------------------------------------------------------------------------------
cd ${TAGHICDIR}
for sn in ${PATIENTS_COAML}
do
paste <(echo -e "# ${sn}") <(grep -w "genome=hg38" ${WORKDIR}/taginfo_unfiltered/${sn}_tagInfo.txt) <(grep "localInteractionFraction=" ${WORKDIR}/taginfo_unfiltered/${sn}_tagInfo.txt) <(grep "interChrInteractionFraction=" ${WORKDIR}/taginfo_unfiltered/${sn}_tagInfo.txt) <(grep "tagsPerBP=" ${WORKDIR}/taginfo_unfiltered/${sn}_tagInfo.txt)
paste <(echo -e "# ${sn}_filtered") <(grep -w "genome=hg38" ${sn}_filtered/tagInfo.txt) <(grep "localInteractionFraction=" ${sn}_filtered/tagInfo.txt) <(grep "interChrInteractionFraction=" ${sn}_filtered/tagInfo.txt) <(grep "tagsPerBP=" ${sn}_filtered/tagInfo.txt)
done
#------------------------------------------------------------------------------------




                      ##########################################
                      #                                        #
                      #  Chromatin Compartment Analysis (PCA)  #
                      #                                        #
                      ##########################################
#------------------------------------------------------------------------------------
cd ${TAGHICDIR}

#50 kb windows with 25 kb resolution
for sn in ${PATIENTS_COAML};do
runHiCpca.pl auto HiC_AML_${sn}_filtered -res 25000 -window 50000 -genome hg38 -cpu 24
bedGraphToBigWig HiC_AML_${sn}_filtered/HiC_AML_${sn}_filtered.25x50kb.PC1.bedGraph ${CHROMSIZES_HG38} ${PCADIR}/HiC_AML_${sn}.25x50kb.PC1.bigWig
done

#define array with names of bedgraphs
bedgraphs25x50=()
for i in ${PATIENTS_COAML}; do
    bGname=HiC_AML_${i}_filtered/HiC_AML_${i}_filtered.25x50kb.PC1.bedGraph
    bedgraphs25x50+=("$bGname")
done
echo "${bedgraphs25x50[@]}"
#------------------------------------------------------------------------------------

#annotate bedGraphs on any PC1.txt file from the set (all have the same "regions" as those are generally specified by window and res)
###IMPORTANT: use "-noblanks"  option otherwise error for R plots
#50kb windows all patients
cd ${TAGHICDIR}
annotatePeaks.pl HiC_AML_ctr_18519_filtered/HiC_AML_ctr_18519_filtered.25x50kb.PC1.txt hg38 -size given -bedGraph ${bedgraphs25x50[@]} -noblanks > ${WORKDIR}/HiC_CohAML_PC1_50KB_annotated.txt
tail -n +2 "${WORKDIR}/HiC_CohAML_PC1_50KB_annotated.txt" | cut -f1,20-44 > ${WORKDIR}/tmp.1.txt
echo $'ID\tctr_16911\tctr_18136\tctr_18519\tctr_19405\tctr_19416\tctr_21047\tctr_21290\tRAD21_UKR186\tRAD21_23039\tRAD21_26830\tRAD21_38455\tSA2_24743\tSA2_27396\tSA2_29728\tSA2_9708' | cat - ${WORKDIR}/tmp.1.txt > ${WORKDIR}/HiC_CohAML_PC1_50KB_annotated.txt
#rm XY chromosomes
grep -v 'chrX' ${WORKDIR}/HiC_CohAML_PC1_50KB_annotated.txt > ${WORKDIR}/tmp.scores3.txt
grep -v 'chrY' ${WORKDIR}/tmp.scores3.txt > ${WORKDIR}/HiC_CohAML_PC1_50KB_annotated.XYrem.txt


#------------------------------------------------------------------------------------
# Analysis for statistically significant changes of PC1 values
## define array with names of bedgraphs by mutation group
bedgraphs25x50_CTRL=()
for i in ${PATIENTS_CTRL}; do
    bGname=HiC_AML_${i}_filtered/HiC_AML_${i}_filtered.25x50kb.PC1.bedGraph
    bedgraphs25x50_CTRL+=("$bGname")
done
echo "${bedgraphs25x50_CTRL[@]}"

bedgraphs25x50_SA2=()
for i in ${PATIENTS_SA2}; do
    bGname=HiC_AML_${i}_filtered/HiC_AML_${i}_filtered.25x50kb.PC1.bedGraph
    bedgraphs25x50_SA2+=("$bGname")
done
echo "${bedgraphs25x50_SA2[@]}"

bedgraphs25x50_RAD21=()
for i in ${PATIENTS_RAD}; do
    bGname=HiC_AML_${i}_filtered/HiC_AML_${i}_filtered.25x50kb.PC1.bedGraph
    bedgraphs25x50_RAD21+=("$bGname")
done
echo "${bedgraphs25x50_RAD21[@]}"

###SA2mut vs CTRL
cd ${TAGHICDIR}
####annotate groups to compare
annotatePeaks.pl HiC_AML_ctr_18519_filtered/HiC_AML_ctr_18519_filtered.25x50kb.PC1.txt hg38 -noblanks -bedGraph \
	${bedgraphs25x50_CTRL[@]} ${bedgraphs25x50_SA2[@]} \
    > ${PCADIR}/CTRLvsSA2.50KB.PC1.ann.txt
####rm XY chrom
grep -v 'chrX'  ${PCADIR}/CTRLvsSA2.50KB.PC1.ann.txt > ${WORKDIR}/tmp.scores3.txt
grep -v 'chrY' ${WORKDIR}/tmp.scores3.txt >  ${PCADIR}/CTRLvsSA2.50KB.PC1.XYrm.ann.txt
cd ${PCADIR}
####run get diff Expression routine of HOMER with the pc1 option which is speicific for HiC PCA data type
getDiffExpression.pl CTRLvsSA2.50KB.PC1.XYrm.ann.txt \
CTRL CTRL CTRL CTRL CTRL CTRL CTRL SA2mut SA2mut SA2mut SA2mut \
-pc1 -export SA2KDvsCTRL.50KB.PC1.XYrm.regions > SA2mutvsCTRL.PC1.50KB.XYrm.diff.txt  

###RAD21mut vs CTRL
cd ${TAGHICDIR}
annotatePeaks.pl HiC_AML_ctr_18519_filtered/HiC_AML_ctr_18519_filtered.25x50kb.PC1.txt hg38 -noblanks -bedGraph \
	${bedgraphs25x50_CTRL[@]} ${bedgraphs25x50_RAD21[@]} \
    > ${PCADIR}/CTRLvsRAD21.50KB.PC1.ann.txt
grep -v 'chrX'  ${PCADIR}/CTRLvsRAD21.50KB.PC1.ann.txt > ${WORKDIR}/tmp.scores3.txt
grep -v 'chrY' ${WORKDIR}/tmp.scores3.txt >  ${PCADIR}/CTRLvsRAD21.50KB.PC1.XYrm.ann.txt
cd ${PCADIR}
getDiffExpression.pl CTRLvsRAD21.50KB.PC1.XYrm.ann.txt \
CTRL CTRL CTRL CTRL CTRL CTRL CTRL RAD21mut RAD21mut RAD21mut RAD21mut \
-pc1 -export RAD21mutvsCTRL.50KB.PC1.XYrm.regions > RAD21mutvsCTRL.PC1.50KB.XYrm.diff.txt  

#------------------------------------------------------------------------------------
# generate average groupwise tracks  for 25x50 PC1 bedgraphs
cd ${TAGHICDIR}
$BEDTOOLS unionbedg -i ${bedgraphs25x50_CTRL[@]} | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' > ${PCADIR}/ave.AML_CTRL_50KB.PC1.bedGraph         
bedGraphToBigWig  ${PCADIR}/ave.AML_CTRL_50KB.PC1.bedGraph ${CHROMSIZES_HG38} ${PCADIR}/ave.AML_CTRL_50KB.PC1.bigWig

$BEDTOOLS unionbedg -i ${bedgraphs25x50_SA2[@]} | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' > ${PCADIR}/ave.AML_SA2mut_50KB.PC1.bedGraph         
bedGraphToBigWig  ${PCADIR}/ave.AML_SA2mut_50KB.PC1.bedGraph ${CHROMSIZES_HG38} ${PCADIR}/ave.AML_SA2mut_50KB.PC1.bigWig

$BEDTOOLS unionbedg -i ${bedgraphs25x50_RAD21[@]} | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' > ${PCADIR}/ave.AML_RAD21mut_50KB.PC1.bedGraph         
bedGraphToBigWig  ${PCADIR}/ave.AML_RAD21mut_50KB.PC1.bedGraph ${CHROMSIZES_HG38} ${PCADIR}/ave.AML_RAD21mut_50KB.PC1.bigWig
#------------------------------------------------------------------------------------
# generate average delta PC1 tracks mut vs ctrl
##requires cov. bedgrapgh for controls
COVbedgraphsCTRL=()
for i in ${PATIENTS_CTRL}; do
    bGname=HiC_AML_${i}.COV.bedGraph
    COVbedgraphsCTRL+=("$bGname")
done
echo "${COVbedgraphsCTRL[@]}"

cd ${COMPACTDIR}
$BEDTOOLS unionbedg -i ${COVbedgraphsCTRL[@]} | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' > ${COMPACTDIR}/ave.AML_CTRL_COV.bedGraph    

##subtract controls from MUT groups
for MUT in ${MUTS}; do
subtractBedGraphs.pl ${PCADIR}/ave.AML_${MUT}_50KB.PC1.bedGraph ${PCADIR}/ave.AML_CTRL_50KB.PC1.bedGraph -cov ${COMPACTDIR}/ave.AML_CTRL_COV.bedGraph -name ${MUT}vsCTRL -center > ${PCADIR}/${MUT}vsCTRL.deltaPC1.bedGraph  
bedGraphToBigWig ${PCADIR}/${MUT}vsCTRL.deltaPC1.bedGraph ${CHROMSIZES_HG38} ${PCADIR}/${MUT}vsCTRL.deltaPC1.bigWig
done

#------------------------------------------------------------------------------------



                      ##########################################
                      #                                        #
                      #         Finding TADs and Loops         #
                      #   (merged directories, merged donors)  # 
                      #                                        #
                      ##########################################
#------------------------------------------------------------------------------------
#create merged TAGDIRS by condition
##define names of filtered tagdir by condition
HIC_CTRL_Filt=()
for i in ${PATIENTS_CTRL}; do
    bGname=HiC_AML_${i}_filtered
    HIC_CTRL_Filt+=("$bGname")
done
echo "${HIC_CTRL_Filt[@]}"

HIC_SA2_Filt=()
for i in ${PATIENTS_SA2}; do
    bGname=HiC_AML_${i}_filtered
    HIC_SA2_Filt+=("$bGname")
done
echo "${HIC_SA2_Filt[@]}"

HIC_RAD_Filt=()
for i in ${PATIENTS_RAD}; do
    bGname=HiC_AML_${i}_filtered
    HIC_RAD_Filt+=("$bGname")
done
echo "${HIC_RAD_Filt[@]}"

##create merged tagdirs
cd ${TAGHICDIR}
makeTagDirectory HiC_AML_CTRL_combined -d ${HIC_CTRL_Filt[@]}

makeTagDirectory HiC_AML_SA2mut_combined -d ${HIC_SA2_Filt[@]}

makeTagDirectory HiC_AML_RAD21mut_combined -d ${HIC_RAD_Filt[@]}

#------------------------------------------------------------------------------------
#Loop and TAD calling
### using two resolutions for loop and tad finding and merged 2D files

##for the 3 patient groups
# res 2500 window 10000 
declare -a SAMPLES=("HiC_AML_CTRL_combined" "HiC_AML_SA2mut_combined" "HiC_AML_RAD21mut_combined" "HiC_adult_CD34_combined")
             cd ${TAGHICDIR}
for SAMPLE in "${SAMPLES[@]}";do
findTADsAndLoops.pl find ${SAMPLE}/ -cpu 12 -res 2500 -window 10000 -genome hg38 \
-p ${BAD_HG38} -o ${LOOPDIR}/${SAMPLE}.2.5K.10K
done

# res 3000 window 15000  
declare -a SAMPLES=("HiC_AML_CTRL_combined" "HiC_AML_SA2mut_combined" "HiC_AML_RAD21mut_combined" "HiC_adult_CD34_combined")    
cd ${TAGHICDIR}
for SAMPLE in "${SAMPLES[@]}";do
findTADsAndLoops.pl find ${SAMPLE}/ -cpu 12 -res 3000 -window 15000 -genome hg38 \
-p ${BAD_HG38} -o ${LOOPDIR}/${SAMPLE}.3K.15K
done

# merge the bed files and create a links file 
declare -a SAMPLES=("HiC_AML_CTRL_combined" "HiC_AML_SA2mut_combined" "HiC_AML_RAD21mut_combined" "HiC_adult_CD34_combined")
cd ${LOOPDIR}
for SAMPLE in "${SAMPLES[@]}";do
merge2Dbed.pl ${SAMPLE}.3K.15K.loop.2D.bed ${SAMPLE}.2.5K.10K.loop.2D.bed > ${SAMPLE}.merged.loop.2D.bed
merge2Dbed.pl ${SAMPLE}.3K.15K.tad.2D.bed ${SAMPLE}.2.5K.10K.tad.2D.bed > ${SAMPLE}.merged.tad.2D.bed
cut -f1-6,8 ${SAMPLE}.merged.loop.2D.bed | LC_COLLATE=C sort -k1,1 -k2,2n > ${SAMPLE}.merged.links
done
#------------------------------------------------------------------------------------

                     #####################################################
                     # Scoring for differential Loops/TAD analysis       #
                     #####################################################

#------------------------------------------------------------------------------------
# Merge the loop and TAD sets of CTRL and comparison group
cd ${LOOPDIR}
## SA2mut vs CTRL #2.5K + 3K res merged sets
merge2Dbed.pl HiC_AML_CTRL_combined.merged.loop.2D.bed HiC_AML_SA2mut_combined.merged.loop.2D.bed > CohAML.SA2vsCTRL.merged.loop.2D.bed 
merge2Dbed.pl HiC_AML_CTRL_combined.merged.tad.2D.bed HiC_AML_SA2mut_combined.merged.tad.2D.bed > CohAML.SA2vsCTRL.merged.tad.2D.bed 

## RAD21mut vs CTRL #2.5K + 3K res merged sets
merge2Dbed.pl HiC_AML_CTRL_combined.merged.loop.2D.bed HiC_AML_RAD21mut_combined.merged.loop.2D.bed > CohAML.RAD21vsCTRL.merged.loop.2D.bed 
merge2Dbed.pl HiC_AML_CTRL_combined.merged.tad.2D.bed HiC_AML_RAD21mut_combined.merged.tad.2D.bed > CohAML.RAD21vsCTRL.merged.tad.2D.bed 

##All patients
merge2Dbed.pl HiC_AML_CTRL_combined.merged.loop.2D.bed HiC_AML_SA2mut_combined.merged.loop.2D.bed HiC_AML_RAD21mut_combined.merged.loop.2D.bed > CohAML.allpat.merged.loop.2D.bed 
merge2Dbed.pl HiC_AML_CTRL_combined.merged.tad.2D.bed HiC_AML_SA2mut_combined.merged.tad.2D.bed HiC_AML_RAD21mut_combined.merged.tad.2D.bed > CohAML.allpat.merged.tad.2D.bed 

#------------------------------------------------------------------------------------

# scoring loops and TADs for the merged sets

## SA2mut vs CTRL
cd ${TAGHICDIR}
###HOMER default normalization
findTADsAndLoops.pl score -tad ${LOOPDIR}/CohAML.SA2vsCTRL.merged.tad.2D.bed -loop ${LOOPDIR}/CohAML.SA2vsCTRL.merged.loop.2D.bed \
-d ${HIC_CTRL_Filt[@]} ${HIC_SA2_Filt[@]} -cpu 28 -o ${LOOPDIR}/CohAML_CTRL_SA2_merged
###normalization to a fixed value (using average seqenecing depth = normTotal)    
cd ${TAGHICDIR}
findTADsAndLoops.pl score -normTotal 500000000 -tad ${LOOPDIR}/CohAML.SA2vsCTRL.merged.tad.2D.bed -loop ${LOOPDIR}/CohAML.SA2vsCTRL.merged.loop.2D.bed \
-d ${HIC_CTRL_Filt[@]} ${HIC_SA2_Filt[@]} -cpu 12 -o ${LOOPDIR}/CohAML_CTRL_SA2_merged.normTotalaverage

     
## RAD21mut vs CTRL
cd ${TAGHICDIR}
###HOMER default normalization
findTADsAndLoops.pl score -tad ${LOOPDIR}/CohAML.RAD21vsCTRL.merged.tad.2D.bed -loop ${LOOPDIR}/CohAML.RAD21vsCTRL.merged.loop.2D.bed \
-d ${HIC_CTRL_Filt[@]} ${HIC_RAD_Filt[@]} -cpu 16 -o ${LOOPDIR}/CohAML_CTRL_RAD21_merged
###normalization to a fixed value (using average seqenecing depth = normTotal)
findTADsAndLoops.pl score -normTotal 500000000 -tad ${LOOPDIR}/CohAML.RAD21vsCTRL.merged.tad.2D.bed -loop ${LOOPDIR}/CohAML.RAD21vsCTRL.merged.loop.2D.bed \
-d ${HIC_CTRL_Filt[@]} ${HIC_RAD_Filt[@]} -cpu 12 -o ${LOOPDIR}/CohAML_CTRL_RAD21_merged.normTotalaverage

## All patients
cd ${TAGHICDIR}
###HOMER default normalization
findTADsAndLoops.pl score -tad ${LOOPDIR}/CohAML.allpat.merged.tad.2D.bed -loop ${LOOPDIR}/CohAML.allpat.merged.loop.2D.bed \
-d ${HIC_CTRL_Filt[@]} ${HIC_SA2_Filt[@]} ${HIC_RAD_Filt[@]} -cpu 16 -o ${LOOPDIR}/CohAML_allpat_merged
###normalization to a fixed value (using average seqenecing depth = normTotal)
findTADsAndLoops.pl score -normTotal 500000000 -tad ${LOOPDIR}/CohAML.allpat.merged.tad.2D.bed -loop ${LOOPDIR}/CohAML.allpat.merged.loop.2D.bed \
-d ${HIC_CTRL_Filt[@]} ${HIC_SA2_Filt[@]} ${HIC_RAD_Filt[@]} -cpu 24 -o ${LOOPDIR}/CohAML_allpat_merged.normTotalaverage
#------------------------------------------------------------------------------------




                     #####################################################
                     #   Prepare Scoring tables for downstream Analyses  #
                     #####################################################

#------------------------------------------------------------------------------------

#Complete patient datasets
cd ${LOOPDIR}
#Loopscores
##keep only numerical columns
cut -f1,11-26 CohAML_allpat_merged.loop.scores.txt | tail -n +2 > tmp.scores.txt
echo $'Gene\tctr_16911\tctr_18136\tctr_18519\tctr_19405\tctr_19416\tctr_21047\tctr_21290\tRAD21_UKR186\tRAD21_23039\tRAD21_26830\tRAD21_38455\tSA2_24743\tSA2_27396\tSA2_29728\tSA2_9708' | cat - tmp.scores.txt > tmp.scores2.txt
##rm X an Y chrom
grep -v 'chrX' tmp.scores2.txt > tmp.scores3.txt
grep -v 'chrY' tmp.scores3.txt > CohAML_allpat_merged.loop.scores.XYrm.Rinput.txt

#TADscores
cd ${LOOPDIR}
cut -f1,11-26 CohAML_allpat_merged.tad.scores.txt | tail -n +2 > tmp.scores.txt
echo $'Gene\tctr_16911\tctr_18136\tctr_18519\tctr_19405\tctr_19416\tctr_21047\tctr_21290\tRAD21_UKR186\tRAD21_23039\tRAD21_26830\tRAD21_38455\tSA2_24743\tSA2_27396\tSA2_29728\tSA2_9708' | cat - tmp.scores.txt > tmp.scores2.txt
##rm X an Y chrom
grep -v 'chrX' tmp.scores2.txt > tmp.scores3.txt
grep -v 'chrY' tmp.scores3.txt > CohAML_allpat_merged.tad.scores.XYrm.Rinput.txt

#------------------------------------------------------------------------------------

# group vs control patient Looppscore datasets
# Total counts for norm2total approach used for Loop differential analysis in DESeq2
cd ${LOOPDIR}
echo -e "Loop\tctr_16911\tctr_18136\tctr_18519\tctr_19405\tctr_19416\tctr_21047\tctr_21290\tSA2_24743\tSA2_27396\tSA2_29728\tSA2_9708" | \
cat - <(awk 'BEGIN {FS=OFS="\t"} {for(i=11;i<=NF;i++) a[i]+=$i} END{l="sum.n2t";i=11;while(i in a){l=l"\t"a[i];i++};print l}' CohAML_CTRL_SA2_merged.normTotalaverage.loop.scores.txt) | \
cat - <(awk 'BEGIN {FS=OFS="\t"} {for(i=11;i<=NF;i++) a[i]+=$i} END{l="sum";i=11;while(i in a){l=l"\t"a[i];i++};print l}' CohAML_CTRL_SA2_merged.loop.scores.XYrm.txt) > ${LOOPDIR}/SA2mutvsCTRL.totalCounts.table.txt

echo -e "Loop\tctr_16911\tctr_18136\tctr_18519\tctr_19405\tctr_19416\tctr_21047\tctr_21290\tRAD21_UKR186_Rep1\tRAD21_23039\tRAD21_26830\tRAD21_38455" | \
cat - <(awk 'BEGIN {FS=OFS="\t"} {for(i=11;i<=NF;i++) a[i]+=$i} END{l="sum.n2t";i=11;while(i in a){l=l"\t"a[i];i++};print l}' CohAML_CTRL_RAD21_merged.normTotalaverage.loop.scores.txt) | \
cat - <(awk 'BEGIN {FS=OFS="\t"} {for(i=11;i<=NF;i++) a[i]+=$i} END{l="sum";i=11;while(i in a){l=l"\t"a[i];i++};print l}' CohAML_CTRL_RAD21_merged.loop.scores.txt) > ${LOOPDIR}/RAD21mutvsCTRL.totalCounts.table.txt

# Homer normalized Loop Score tables
#For SA2mut vs CTRL
cd ${LOOPDIR}
## remove non count columns and add IDs as column names
cut -f1,11-23 CohAML_CTRL_SA2_merged.loop.scores.XYrm.txt | tail -n +2 > tmp.scores.txt
echo $'Loop\tctr_16911\tctr_18136\tctr_18519\tctr_19405\tctr_19416\tctr_21047\tctr_21290\tSA2_24743\tSA2_27396\tSA2_29728\tSA2_9708' | cat - tmp.scores.txt > CohAML_CTRL_SA2_merged.loop.scores.XYrm.Rinput.txt
cut -f1,2,3,7 CohAML_CTRL_SA2_merged.loop.scores.XYrm.txt | tail -n +2 > tmp.loopstartstop.txt
echo $'Loop\tchr\tstart\tend' | cat - tmp.loopstartstop.txt > tmp.loopstartstop2.txt
## information on loop coordinates separately
cut -f1,2,3,4,5,6,7 CohAML_CTRL_SA2_merged.loop.scores.XYrm.txt | tail -n +2 > tmp.loopstartstop.txt
echo $'Loop\tchr1\tstart1\tend1\tchr2\tstart2\tend2' | cat - tmp.loopstartstop.txt > CohAML_CTRL_SA2.loopstartstop.txt

#For RAD21mut vs CTRL
cd ${LOOPDIR}
## remove non count columns and add IDs as column names
cut -f1,11-23 CohAML_CTRL_RAD21_merged.loop.scores.XYrm.txt | tail -n +2 > tmp.scores.txt
echo $'Loop\tctr_16911\tctr_18136\tctr_18519\tctr_19405\tctr_19416\tctr_21047\tctr_21290\tRAD21_UKR186_Rep1\tRAD21_23039\tRAD21_26830\tRAD21_38455' | cat - tmp.scores.txt > CohAML_CTRL_RAD21_merged.loop.scores.XYrm.Rinput.txt
cut -f1,2,3,7 CohAML_CTRL_RAD21_merged.loop.scores.XYrm.txt | tail -n +2 > tmp.loopstartstop.txt
echo $'Loop\tchr\tstart\tend' | cat - tmp.loopstartstop.txt > tmp.loopstartstop2.txt
## information on loop coordinates separately
cut -f1,2,3,4,5,6,7 CohAML_CTRL_RAD21_merged.loop.scores.XYrm.txt | tail -n +2 > tmp.loopstartstop.txt
echo $'Loop\tchr1\tstart1\tend1\tchr2\tstart2\tend2' | cat - tmp.loopstartstop.txt > CohAML_CTRL_RAD21.loopstartstop.txt
#------------------------------------------------------------------------------------





                     #####################################################
                     #       Differential TAD Analysis using Homer       #
                     #####################################################

#------------------------------------------------------------------------------------
#test differential enrichement of TADs in groupwise comparisons #using TADs at 2.5+3KB res merged dataset
###SA2mut vs CTRL
cd ${LOOPDIR}
grep -v 'chrX' CohAML_CTRL_SA2_merged.tad.scores.txt > tmp.scores3.txt
grep -v 'chrY' tmp.scores3.txt > CohAML_CTRL_SA2_merged.tad.scores.XYrm.txt
getDiffExpression.pl CohAML_CTRL_SA2_merged.tad.scores.XYrm.txt  CTRL CTRL CTRL CTRL CTRL CTRL CTRL SA2mut SA2mut SA2mut SA2mut\
> SA2mutvsCTRL.merged.diff.XYrm.tad.txt
####RAD21mut vs CTRL
grep -v 'chrX' CohAML_CTRL_RAD21_merged.tad.scores.txt > tmp.scores3.txt
grep -v 'chrY' tmp.scores3.txt > CohAML_CTRL_RAD21_merged.tad.scores.XYrm.txt
cd ${LOOPDIR}
getDiffExpression.pl CohAML_CTRL_RAD21_merged.tad.scores.XYrm.txt CTRL CTRL CTRL CTRL CTRL CTRL CTRL RAD21mut RAD21mut RAD21mut RAD21mut\
> RAD21mutvsCTRL.merged.diff.tad.txt
#------------------------------------------------------------------------------------
