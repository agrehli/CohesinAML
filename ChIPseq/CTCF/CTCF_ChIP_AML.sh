#!/bin/bash
# by Alexander Fischer, July 2021


###############################################################################
###############################################################################
##                                                                           ##
##       Analysis of CTCF ChIPseq in AML patients                            ##
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


# general files
CHROMSIZES_HG38="${DIR_SOFT}/viewer/IGV/IGVTools_2.3.98/genomes/GRCh38.PRI_p10.chrom.sizes"
CHAINFILE="${DIR_PKG}/genome/chainFiles/hg19ToHg38.over.chain"
BLACKLIST_HG38="${DIR_DATA}/analysis/generalStuff/annotation/GRCh38/hg38.blacklist.bed"

# general directories
TMPDIR="/loctmp"
## directories with data created by mapChIP.sh pipeline
TAGDIR="${DIR_DATA}/processedData/tagDir/chromatin/GRCh38/ChIP/Cohesin_AML" #Tag directories (HOMER)
BWDIR="${DIR_DATA}/processedData/bigWig/chromatin/GRCh38/ChIP/Cohesin_AML" #bigwigs
INPUTDIR="${DIR_DATA}/processedData/tagDir/DNA/GRCh38/Input/AML" # Tag directorie (HOMER) of input DNA-sequencing data
CNVDIR="${DIR_DATA}/processedData/mapping/DNA/GRCh38/Input/AML/CNVdata" # processed input DNA-sequencing data CNV FreeC output

## directories for analysis
PROJECTDIR="${DIR_DATA}/analysis/project_cohesin/Cohesin_AML"
WORKDIR="${PROJECTDIR}/ChIP_analysis"
PEAKDIR="${WORKDIR}/peaks"
ANNDIR="${WORKDIR}/annotation_tables"
DIFFPEAKS="${WORKDIR}/diffPeaks"
FIGURESDIR="${WORKDIR}/figures"
MERGEDBW="${WORKDIR}/mergedBigWigs"
MOTIFDIR="${WORKDIR}/motifs"


                         #####################################
                         #        peak finding for QC        #
                         #####################################
#------------------------------------------------------------------------------------                         
                        
#Declare Samples
PATIENTS_RAD2='Rad21_27650 Rad21_28420 Rad21_29047 Rad21_29992 Rad21_38455 Rad21_UKR186 Rad21_23039 Rad21_26830 '
PATIENTS_RAD4='27650 28420 29047 29992 38455 UKR186 23039 26830 '
PATIENTS_SA2='SA2_BB007 SA2_15365 SA2_15640 SA2_22464 SA2_25458 SA2_31501 SA2_9708 SA2_24743 SA2_27396 SA2_29728 '
PATIENTS_CTRL='ctr_11374 ctr_14094 ctr_18367 ctr_19049 ctr_19284 ctr_20180 ctr_20458 ctr_20899 ctr_22023 ctr_9873 ctr_16911 ctr_18136 ctr_18519 ctr_19405 ctr_19416 ctr_21047 ctr_21290 '
ALLPAT_CTCFchips=$PATIENTS_CTRL$PATIENTS_RAD2$PATIENTS_SA2
#------------------------------------------------------------------------------------                         

cd ${TAGDIR}
#loop through patients and call peaks with basic parameters (HOMER)
#CTRL AML
for patID in ${PATIENTS_CTRL};do
findPeaks ${TAGDIR}/ChIP_AML_${patID}_CTCF -i ${INPUTDIR}/ChIP_IN_AML_${patID} -style factor -o auto
done
#SA2 AML
for patID in ${PATIENTS_SA2};do
findPeaks ${TAGDIR}/ChIP_AML_${patID}_CTCF -i ${INPUTDIR}/ChIP_IN_AML_${patID} -style factor -o auto
done
#RAD21 AML
for patID in ${PATIENTS_RAD4};do
findPeaks ${TAGDIR}/ChIP_AML_Rad21_${patID}_CTCF -i ${INPUTDIR}/ChIP_IN_AML_RAD21_${patID} -style factor -o auto
done

#------------------------------------------------------------------------------------                         

# show values
cd ${TAGDIR}
for SAMPLE in ${ALLPAT_CTCFchips};do
paste <(echo -e "# ChIP_AML_${SAMPLE}_CTCF") <(grep -w "# total peaks =" ChIP_AML_${SAMPLE}_CTCF/peaks.txt) <(grep -w "# Approximate IP efficiency" ChIP_AML_${SAMPLE}_CTCF/peaks.txt) <(grep -w "# Total tags =" ChIP_AML_${SAMPLE}_CTCF/peaks.txt)
done
#save to file
echo -e "Sample_ID\tTotalPeaks\tIPefficiency\tTotalTags" > $PEAKDIR/CTCFpeaks.default.txt
for SAMPLE in ${ALLPAT_CTCFchips}; do
    sample_name="ChIP_AML_${SAMPLE}_CTCF"
    TotalPeaks=$(grep -w "# total peaks =" $sample_name/peaks.txt | awk '{print $5}')
    IPefficiency=$(grep -w "# Approximate IP efficiency =" $sample_name/peaks.txt | awk '{print $6}')
    TotalTags=$(grep -w "# Total tags =" $sample_name/peaks.txt | awk '{print $5}')
    echo -e "${sample_name}\t${TotalPeaks}\t${IPefficiency}\t${TotalTags}" >> $PEAKDIR/CTCFpeaks.default.txt
done
#------------------------------------------------------------------------------------                         



                         #####################################
                         #        normalize TAGDIRS by CNV   #
                         #####################################

#------------------------------------------------------------------------------------                         
#CTRL AML
for patID in ${PATIENTS_CTRL};do
normalizeTagDirByCopyNumber.pl ${TAGDIR}/ChIP_AML_${patID}_CTCF -cnv ${CNVDIR}/ChIP_IN_AML_${patID}.sam_CNVs
done
#STAG2 AML
for patID in ${PATIENTS_SA2};do
normalizeTagDirByCopyNumber.pl ${TAGDIR}/ChIP_AML_${patID}_CTCF -cnv ${CNVDIR}/ChIP_IN_AML_${patID}.sam_CNVs
done
#RAD21 AML
for patID in ${PATIENTS_RAD4};do
normalizeTagDirByCopyNumber.pl ${TAGDIR}/ChIP_AML_Rad21_${patID}_CTCF -cnv ${CNVDIR}/ChIP_IN_AML_RAD21_${patID}.sam_CNVs
done
#------------------------------------------------------------------------------------                         

                    #############################################
                    #  peak finding for QC in CNV norm tagdirs  #
                    #############################################
#------------------------------------------------------------------------------------                         

#loop through  patients and find peaks with specified stringency settings
#CTRL AML
for patID in ${PATIENTS_CTRL};do
findPeaks ${TAGDIR}/ChIP_AML_${patID}_CTCF_CNVnorm -i ${INPUTDIR}/ChIP_IN_AML_${patID} -style factor -o auto
findPeaks ${TAGDIR}/ChIP_AML_${patID}_CTCF_CNVnorm -i ${INPUTDIR}/ChIP_IN_AML_${patID} -style factor -tbp 1 -fdr 0.000001 -ntagThreshold 10 -o ${TAGDIR}/ChIP_AML_${patID}_CTCF_CNVnorm/stringentPeaks.txt
done
#SA2 AML
for patID in ${PATIENTS_SA2};do
findPeaks ${TAGDIR}/ChIP_AML_${patID}_CTCF_CNVnorm -i ${INPUTDIR}/ChIP_IN_AML_${patID} -style factor -o auto
findPeaks ${TAGDIR}/ChIP_AML_${patID}_CTCF_CNVnorm -i ${INPUTDIR}/ChIP_IN_AML_${patID} -style factor -tbp 1 -fdr 0.000001 -ntagThreshold 10 -o ${TAGDIR}/ChIP_AML_${patID}_CTCF_CNVnorm/stringentPeaks.txt
done
#RAD21 AML
for patID in ${PATIENTS_RAD4};do
findPeaks ${TAGDIR}/ChIP_AML_Rad21_${patID}_CTCF_CNVnorm -i ${INPUTDIR}/ChIP_IN_AML_RAD21_${patID} -style factor -o auto
findPeaks ${TAGDIR}/ChIP_AML_Rad21_${patID}_CTCF_CNVnorm -i ${INPUTDIR}/ChIP_IN_AML_RAD21_${patID} -style factor -tbp 1 -fdr 0.000001 -ntagThreshold 10 -o ${TAGDIR}/ChIP_AML_Rad21_${patID}_CTCF_CNVnorm/stringentPeaks.txt
done
#------------------------------------------------------------------------------------                         
#show values for all CTCF chips
cd ${TAGDIR}
for SAMPLE in ${ALLPAT_CTCFchips};do
paste <(echo -e "# ChIP_AML_${SAMPLE}_CTCF_CNVnorm") <(grep -w "# total peaks =" ChIP_AML_${SAMPLE}_CTCF_CNVnorm/peaks.txt) <(grep -w "# Approximate IP efficiency" ChIP_AML_${SAMPLE}_CTCF_CNVnorm/peaks.txt) <(grep -w "# Total tags =" ChIP_AML_${SAMPLE}_CTCF_CNVnorm/peaks.txt)
paste <(echo -e "# ChIP_AML_${SAMPLE}_CTCF_CNVnorm") <(grep -w "# total peaks =" ChIP_AML_${SAMPLE}_CTCF_CNVnorm/stringentPeaks.txt) <(grep -w "# Approximate IP efficiency" ChIP_AML_${SAMPLE}_CTCF_CNVnorm/stringentPeaks.txt) <(grep -w "# Total tags =" ChIP_AML_${SAMPLE}_CTCF_CNVnorm/stringentPeaks.txt)
done

#save to files
echo -e "Sample_ID\tTotalPeaks\tIPefficiency\tTotalTags" > $PEAKDIR/CTCFpeaks.CNVnorm.txt
for SAMPLE in ${ALLPAT_CTCFchips}; do
    sample_name="ChIP_AML_${SAMPLE}_CTCF_CNVnorm"
    TotalPeaks=$(grep -w "# total peaks =" $sample_name/peaks.txt | awk '{print $5}')
    IPefficiency=$(grep -w "# Approximate IP efficiency =" $sample_name/peaks.txt | awk '{print $6}')
    TotalTags=$(grep -w "# Total tags =" $sample_name/peaks.txt | awk '{print $5}')
    echo -e "${sample_name}\t${TotalPeaks}\t${IPefficiency}\t${TotalTags}" >> $PEAKDIR/CTCFpeaks.CNVnorm.txt
done

echo -e "Sample_ID\tTotalPeaks\tIPefficiency\tTotalTags" > $PEAKDIR/CTCFpeaks.strigent.CNVnorm.txt
for SAMPLE in ${ALLPAT_CTCFchips}; do
    sample_name="ChIP_AML_${SAMPLE}_CTCF_CNVnorm"
    TotalPeaks=$(grep -w "# total peaks =" $sample_name/stringentPeaks.txt | awk '{print $5}')
    IPefficiency=$(grep -w "# Approximate IP efficiency =" $sample_name/stringentPeaks.txt | awk '{print $6}')
    TotalTags=$(grep -w "# Total tags =" $sample_name/stringentPeaks.txt | awk '{print $5}')
    echo -e "${sample_name}\t${TotalPeaks}\t${IPefficiency}\t${TotalTags}" >> $PEAKDIR/CTCFpeaks.strigent.CNVnorm.txt
done
#------------------------------------------------------------------------------------                         

                     ######################################
                     #CNV norm BIGwigs and average BIGwigs#
                     ######################################
#------------------------------------------------------------------------------------                         
# CNV norm. bw file: create bigwig from corrected tagdir
for SAMPLE in ${ALLPAT_CTCFchips}
do
makeUCSCfile ${TAGDIR}/ChIP_AML_${SAMPLE}_CTCF*_CNVnorm -o ${BWDIR}/ChIP_AML_${SAMPLE}_CTCF_CNVnorm.bigWig -fragLength 200 -bigWig ${CHROMSIZES_HG38} -fsize 1e20
done
#------------------------------------------------------------------------------------                         
# average BW file by condition (create from CNVnorm BWs)
##define vectors with bigwig names by group
normRad21CTCFbw=()
for i in ${PATIENTS_RAD2}; do
    bGname=ChIP_AML_${i}_CTCF_CNVnorm.bigWig
    normRad21CTCFbw+=("$bGname")
done
echo "${normRad21CTCFbw[@]}"

normSA2CTCFbw=()
for i in ${PATIENTS_SA2}; do
    bGname=ChIP_AML_${i}_CTCF_CNVnorm.bigWig
    normSA2CTCFbw+=("$bGname")
done
echo "${normSA2CTCFbw[@]}"

normCTRLCTCFbw=()
for i in ${PATIENTS_CTRL}; do
    bGname=ChIP_AML_${i}_CTCF_CNVnorm.bigWig
    normCTRLCTCFbw+=("$bGname")
done
echo "${normCTRLCTCFbw[@]}"
#------------------------------------------------------------------------------------                         
##run AverageBigwig.pl script and create bedgraphs in addition
cd ${BWDIR}
#CTRL AML
myAverageBigWig.pl -bw ${normCTRLCTCFbw[@]} -chr ${CHROMSIZES_HG38} -o ${MERGEDBW}/ave.AML_CTRL_CTCF_CNVnorm.bigwig
bigWigToBedGraph ${MERGEDBW}/ave.AML_CTRL_CTCF_CNVnorm.bigwig ${MERGEDBW}/ave.AML_CTRL_CTCF.CNVnorm.bedGraph
#STAG2mut AML
myAverageBigWig.pl -bw ${normSA2CTCFbw[@]} -chr ${CHROMSIZES_HG38} -o ${MERGEDBW}/ave.AML_SA2mut_CTCF_CNVnorm.bigwig
bigWigToBedGraph ${MERGEDBW}/ave.AML_SA2mut_CTCF_CNVnorm.bigwig ${MERGEDBW}/ave.AML_SA2mut_CTCF.CNVnorm.bedGraph
#CTRL AML
myAverageBigWig.pl -bw ${normCTRLCTCFbw[@]} -chr ${CHROMSIZES_HG38} -o ${MERGEDBW}/ave.AML_CTRL_CTCF_CNVnorm.bigwig
bigWigToBedGraph ${MERGEDBW}/ave.AML_CTRL_CTCF_CNVnorm.bigwig ${MERGEDBW}/ave.AML_CTRL_CTCF.CNVnorm.bedGraph
#RAD21mut AML
myAverageBigWig.pl -bw ${normRad21CTCFbw[@]} -chr ${CHROMSIZES_HG38} -o ${MERGEDBW}/ave.AML_Rad21mut_CTCF_CNVnorm.bigwig
bigWigToBedGraph ${MERGEDBW}/ave.AML_Rad21mut_CTCF_CNVnorm.bigwig ${MERGEDBW}/ave.AML_Rad21mut_CTCF.CNVnorm.bedGraph
#------------------------------------------------------------------------------------                         


                     ################################################
                     #merged Peaksets based on individual peakfiles #
                     ################################################
#------------------------------------------------------------------------------------                         
#regular peaks
CTCFpeaks_norm=()
for i in ${ALLPAT_CTCFchips}; do
    bGname=ChIP_AML_${i}_CTCF*CNVnorm/peaks.txt
    CTCFpeaks_norm+=("$bGname")
done
echo "${CTCFpeaks_norm[@]}"

cd ${TAGDIR}
mergePeaks ${CTCFpeaks_norm[@]} -code > ${PEAKDIR}/Allpat_mergePeaks_CTCF.txt
#filter merged peakset
pos2bed.pl ${PEAKDIR}/Allpat_mergePeaks_CTCF.txt > ${TMPDIR}/tmp.1.bed
	$BEDTOOLS intersect -a ${TMPDIR}/tmp.1.bed -b $BLACKLIST_HG38 -v > ${TMPDIR}/tmp.2.bed
	bed2pos.pl ${TMPDIR}/tmp.2.bed > ${TMPDIR}/tmp.1.txt
	filter4Mappability.sh -p ${TMPDIR}/tmp.1.txt -g hg38 -f 0.8 -s 75
	pos2bed.pl ${TMPDIR}/tmp.1.mapScoreFiltered.txt > ${PEAKDIR}/Allpat_mergePeaks_CTCF.filtered.peaks.bed
	bed2pos.pl ${PEAKDIR}/Allpat_mergePeaks_CTCF.filtered.peaks.bed > ${PEAKDIR}/Allpat_mergePeaks_CTCF.filtered.peaks.txt

#stringent peaks
CTCFstrpeaks_norm=()
for i in ${ALLPAT_CTCFchips}; do
    bGname=ChIP_AML_${i}_CTCF*CNVnorm/stringentPeaks.txt
    CTCFstrpeaks_norm+=("$bGname")
done
echo "${CTCFstrpeaks_norm[@]}"

mergePeaks ${CTCFstrpeaks_norm[@]} -code > ${PEAKDIR}/Allpat_mergePeaks_CTCF_stringent.txt
#filter merged peaksets
pos2bed.pl ${PEAKDIR}/Allpat_mergePeaks_CTCF_stringent.txt > ${TMPDIR}/tmp.1.bed
	$BEDTOOLS intersect -a ${TMPDIR}/tmp.1.bed -b $BLACKLIST_HG38 -v > ${TMPDIR}/tmp.2.bed
	bed2pos.pl ${TMPDIR}/tmp.2.bed > ${TMPDIR}/tmp.1.txt
	filter4Mappability.sh -p ${TMPDIR}/tmp.1.txt -g hg38 -f 0.8 -s 75
	pos2bed.pl ${TMPDIR}/tmp.1.mapScoreFiltered.txt > ${PEAKDIR}/Allpat_mergePeaks_CTCF_stringent.filtered.peaks.bed
	bed2pos.pl ${PEAKDIR}/Allpat_mergePeaks_CTCF_stringent.filtered.peaks.bed > ${PEAKDIR}/Allpat_mergePeaks_CTCF_stringent.filtered.peaks.txt
#------------------------------------------------------------------------------------                         



                     ####################################################
                     #     Annotation of tagdirs to peak sets           #
                     ####################################################
#------------------------------------------------------------------------------------ 
##define vectors with tagdir names by group
normRad21CTCF=()
for i in ${PATIENTS_RAD}; do
    bGname=ChIP_AML_${i}_CTCF_CNVnorm
    normRad21CTCF+=("$bGname")
done
echo "${normRad21CTCF[@]}"

normSA2CTCF=()
for i in ${PATIENTS_SA2}; do
    bGname=ChIP_AML_${i}_CTCF_CNVnorm
    normSA2CTCF+=("$bGname")
done
echo "${normSA2CTCF[@]}"

normCTRLCTCF=()
for i in ${PATIENTS_CTRL}; do
    bGname=ChIP_AML_${i}_CTCF_CNVnorm
    normCTRLCTCF+=("$bGname")
done
echo "${normCTRLCTCF[@]}"
#------------------------------------------------------------------------------------                         
###annotate all CNVnorm tagdirs on the filtered merged peak set, filter out X and Y!
#use regular peak set
annotatePeaks.pl "${PEAKDIR}/Allpat_mergePeaks_CTCF.filtered.peaks.txt" hg38 -size 250 -d \
${normCTRLCTCF[@]} ${normRad21CTCF[@]} ${normSA2CTCF[@]} \
-raw -cpu 24 > ${ANNDIR}/Allpat_CTCF.peaks.ann.txt
#remove sex chromosomes
grep -v 'chrX' ${ANNDIR}/Allpat_CTCF.peaks.ann.txt > ${TMPDIR}/tmp.filt.txt
grep -v 'chrY' ${TMPDIR}/tmp.filt.txt > ${ANNDIR}/Allpat_CTCF.peaks.ann.filt.txt
#remove non-essential cols
#cut -f2-7 --complement ${TMPDIR}/CTCF.filt2.txt > ${ANNDIR}/Allpat_CTCF.peaks.ann.filt.txt
rm ${TMPDIR}/tmp.filt*
tail -n +2 ${ANNDIR}/Allpat_CTCF.peaks.ann.filt.txt | cut -f1,20-100 > ${ANNDIR}/Allpat_CTCF.peaks.ann.filt.Rinput.txt
cut -f 1-6 ${PEAKDIR}/Allpat_mergePeaks_CTCF.filtered.peaks.txt > ${ANNDIR}/CTCF.filtered.peaks.Rinput.txt
#------------------------------------------------------------------------------------                         


