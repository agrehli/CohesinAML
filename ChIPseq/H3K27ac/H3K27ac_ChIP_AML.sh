#!/bin/bash
# by Alexander Fischer, DEC 2021

###############################################################################
###############################################################################
##                                                                           ##
# Analysis of H3k72ac ChIPseq data  in Cohesin-mutant AML patients           ##
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

## directories with data created by mapChIP.sh pipeline
TAGDIR="${DIR_DATA}/processedData/tagDir/chromatin/GRCh38/ChIP/Cohesin_AML" #Tag directories (HOMER)
BWDIR="${DIR_DATA}/processedData/bigWig/chromatin/GRCh38/ChIP/Cohesin_AML" #bigwigs
INPUTDIR="${DIR_DATA}/processedData/tagDir/DNA/GRCh38/Input/AML" # Tag directorie (HOMER) of input DNA-sequencing data
CNVDIR="${DIR_DATA}/processedData/mapping/DNA/GRCh38/Input/AML/CNVdata" # processed input DNA-sequencing data CNV FreeC output

## directories for analysis
PROJECTDIR="${DIR_DATA}/analysis/project_cohesin/Cohesin_AML"
WORKDIR="${PROJECTDIR}/ChIP_analysis/H3K27ac"
PEAKDIR="${WORKDIR}/peaks"
ANNDIR="${WORKDIR}/annotation_tables"
DIFFDIR="${WORKDIR}/diffPeaks"
FIGURESDIR="${WORKDIR}/figures"
MERGEDBW="${WORKDIR}/mergedBigWigs"
MOTIFDIR="${WORKDIR}/motifs"

#create the new directories
mkdir $WORKDIR
mkdir $PEAKDIR
mkdir $ANNDIR
mkdir $DIFFDIR
mkdir $FIGURESDIR
mkdir $MERGEDBW
mkdir $MOTIFDIR



                         #####################################
                         #        peak finding for QC        #
                        #####################################
#------------------------------------------------------------------------------------                         

#Declare samples
#all Cohesin AML + CTRL
PATIENTS_RAD='Rad21_23039 Rad21_26830 Rad21_27650 Rad21_28420 Rad21_29047 Rad21_29992 Rad21_38455 Rad21_UKR186 '
PATIENTS_RAD4='23039 26830 27650 28420 29047 29992 38455 UKR186'
PATIENTS_SA2='SA2_15365 SA2_15640 SA2_19142 SA2_22464 SA2_22525 SA2_24743 SA2_25458 SA2_27396 SA2_27493 SA2_28041 SA2_29728 SA2_31501 SA2_32635 SA2_9708 SA2_BB007 '
PATIENTS_CTRL='ctr_11374 ctr_14094 ctr_16911 ctr_18136 ctr_18367 ctr_18519 ctr_19049 ctr_19284 ctr_19405 ctr_19416 ctr_20180 ctr_20458 ctr_20899 ctr_21047 ctr_21290 ctr_22023 ctr_9873 '
ALLPAT=$PATIENTS_RAD$PATIENTS_SA2$PATIENTS_CTRL

#------------------------------------------------------------------------------------                         

#loop through patients and find peaks in histone style (HOMER) wiht otherwise default settings
#CTRL AML
for patID in ${PATIENTS_CTRL};do
findPeaks ${TAGDIR}/ChIP_AML_${patID}_H3K27ac -i ${INPUTDIR}/ChIP_IN_AML_${patID} -style histone -o auto
done
for patID in ${PATIENTS_SA2};do
findPeaks ${TAGDIR}/ChIP_AML_${patID}_H3K27ac -i ${INPUTDIR}/ChIP_IN_AML_${patID} -style histone -o auto
done
for patID in ${PATIENTS_RAD4};do
findPeaks ${TAGDIR}/ChIP_AML_Rad21_${patID}_H3K27ac -i ${INPUTDIR}/ChIP_IN_AML_RAD21_${patID} -style histone -o auto
done
#------------------------------------------------------------------------------------                         
# show values
cd ${TAGDIR}
for SAMPLE in ${ALLPAT}
do
paste <(echo -e "# ChIP_AML_${SAMPLE}_H3K27ac") <(grep -w "# total peaks =" ChIP_AML_${SAMPLE}_H3K27ac/regions.txt) <(grep -w "# Approximate IP efficiency" ChIP_AML_${SAMPLE}_H3K27ac/regions.txt) <(grep -w "# Total tags =" ChIP_AML_${SAMPLE}_H3K27ac/regions.txt)
done
#save to file
echo -e "Sample_ID\tTotalPeaks\tIPefficiency\tTotalTags" > $PEAKDIR/H3k27peaks.styleHistone.txt
for SAMPLE in ${ALLPAT}; do
    sample_name="ChIP_AML_${SAMPLE}_H3K27ac"
    TotalPeaks=$(grep -w "# total peaks =" $sample_name/regions.txt | awk '{print $5}')
    IPefficiency=$(grep -w "# Approximate IP efficiency =" $sample_name/regions.txt | awk '{print $6}')
    TotalTags=$(grep -w "# Total tags =" $sample_name/regions.txt | awk '{print $5}')
    echo -e "${sample_name}\t${TotalPeaks}\t${IPefficiency}\t${TotalTags}" >> $PEAKDIR/H3k27peaks.styleHistone.txt
done
#------------------------------------------------------------------------------------                         


                         ###########################################################
                         #   normalize TAGDIRS by CNV and find regions/peaks       #
                         ###########################################################
#loop through  patients and find peaks with specified stringency and size settings
#------------------------------------------------------------------------------------ 
cd ${TAGDIR}
#CTRL AML
for patID in ${PATIENTS_CTRL};do
normalizeTagDirByCopyNumber.pl ${TAGDIR}/ChIP_AML_${patID}_H3K27ac -cnv ${CNVDIR}/ChIP_IN_AML_${patID}.sam_CNVs
findPeaks ${TAGDIR}/ChIP_AML_${patID}_H3K27ac_CNVnorm -i ${INPUTDIR}/ChIP_IN_AML_${patID} -region -size 250 -L 0 -F 5 -minDist 350 -fdr 0.00001 -ntagThreshold 10 -o auto
done
#STAG2 AML
for patID in ${PATIENTS_SA2};do
normalizeTagDirByCopyNumber.pl ${TAGDIR}/ChIP_AML_${patID}_H3K27ac -cnv ${CNVDIR}/ChIP_IN_AML_${patID}.sam_CNVs
findPeaks ${TAGDIR}/ChIP_AML_${patID}_H3K27ac_CNVnorm -i ${INPUTDIR}/ChIP_IN_AML_${patID} -region -size 250 -L 0 -F 5 -minDist 350 -fdr 0.00001 -ntagThreshold 10 -o auto
done
#RAD21 AML
for patID in ${PATIENTS_RAD4};do
normalizeTagDirByCopyNumber.pl ${TAGDIR}/ChIP_AML_Rad21_${patID}_H3K27ac -cnv ${CNVDIR}/ChIP_IN_AML_RAD21_${patID}.sam_CNVs
findPeaks ${TAGDIR}/ChIP_AML_Rad21_${patID}_H3K27ac_CNVnorm -i ${INPUTDIR}/ChIP_IN_AML_RAD21_${patID} -region -size 250 -L 0 -F 5 -minDist 350 -fdr 0.00001 -ntagThreshold 10 -o auto
done
#------------------------------------------------------------------------------------ 
#save peakstats to file
cd ${TAGDIR}
echo -e "Sample_ID\tTotalPeaks\tIPefficiency\tTotalTags" > $PEAKDIR/H3k27peaks.CNVnorm.summary.txt
for SAMPLE in ${ALLPAT}; do
    sample_name="ChIP_AML_${SAMPLE}_H3K27ac_CNVnorm"
    TotalPeaks=$(grep -w "# total peaks =" $sample_name/peaks.txt | awk '{print $5}')
    IPefficiency=$(grep -w "# Approximate IP efficiency =" $sample_name/peaks.txt | awk '{print $6}')
    TotalTags=$(grep -w "# Total tags =" $sample_name/peaks.txt | awk '{print $5}')
    echo -e "${sample_name}\t${TotalPeaks}\t${IPefficiency}\t${TotalTags}" >> $PEAKDIR/H3k27peaks.CNVnorm.summary.txt
done
#------------------------------------------------------------------------------------ 



                     ######################################
                     #CNV norm BIGwigs and average BIGwigs#
                     ######################################
#------------------------------------------------------------------------------------ 
#CNV norm bw file
for SAMPLE in ${ALLPAT};do
makeUCSCfile ${TAGDIR}/ChIP_AML_${SAMPLE}_H3K27ac_CNVnorm -o ${BWDIR}/ChIP_AML_${SAMPLE}_H3K27ac_CNVnorm.bigWig -fragLength 200 -bigWig ${CHROMSIZES_HG38} -fsize 1e20
done

#------------------------------------------------------------------------------------ 
#define list of names of individual BW files by condition
normRad21bw=()
for i in ${PATIENTS_RAD}; do
    bGname=ChIP_AML_${i}_H3K27ac_CNVnorm.bigWig
    normRad21bw+=("$bGname")
done
echo "${normRad21bw[@]}"

normSA2bw=()
for i in ${PATIENTS_SA2}; do
    bGname=ChIP_AML_${i}_H3K27ac_CNVnorm.bigWig
    normSA2bw+=("$bGname")
done
echo "${normSA2bw[@]}"

normCTRLbw=()
for i in ${PATIENTS_CTRL}; do
    bGname=ChIP_AML_${i}_H3K27ac_CNVnorm.bigWig
    normCTRLbw+=("$bGname")
done
echo "${normCTRLbw[@]}"
#------------------------------------------------------------------------------------ 

#create average bigwigs (and bedgraphs) by patient group
cd ${BWDIR}
#CTRL AML
myAverageBigWig.pl -bw ${normCTRLbw[@]} -chr ${CHROMSIZES_HG38} -o ${MERGEDBW}/ave.AML_CTRL_H3K27ac_CNVnorm.bigwig
bigWigToBedGraph ${MERGEDBW}/ave.AML_CTRL_H3K27ac_CNVnorm.bigwig ${MERGEDBW}/ave.AML_CTRL_H3K27ac.CNVnorm.bedGraph

#STAG2mut AML
myAverageBigWig.pl -bw ${normSA2bw[@]} -chr ${CHROMSIZES_HG38} -o ${MERGEDBW}/ave.AML_SA2mut_H3K27ac_CNVnorm.bigwig
bigWigToBedGraph ${MERGEDBW}/ave.AML_SA2mut_H3K27ac_CNVnorm.bigwig ${MERGEDBW}/ave.AML_SA2mut_H3K27ac.CNVnorm.bedGraph

#RAD21mut AML
myAverageBigWig.pl -bw ${normRad21bw[@]} -chr ${CHROMSIZES_HG38} -o ${MERGEDBW}/ave.AML_Rad21mut_H3K27ac_CNVnorm.bigwig
bigWigToBedGraph ${MERGEDBW}/ave.AML_Rad21mut_H3K27ac_CNVnorm.bigwig ${MERGEDBW}/ave.AML_Rad21mut_H3K27ac.CNVnorm.bedGraph
#------------------------------------------------------------------------------------ 



                     ###################################################
                     # merged Peaksets based on individual peakfiles   #
                     ###################################################
#------------------------------------------------------------------------------------ 
#define names of Peakfiles
#peaks called from CNVnorm tagdirs
H3K27ac_peaks_norm=()
for i in ${ALLPAT}; do
    bGname=ChIP_AML_${i}_H3K27ac_CNVnorm/peaks.txt
    H3K27ac_peaks_norm+=("$bGname")
done
echo "${H3K27ac_peaks_norm[@]}"
#merge peaks of individual patients
cd ${TAGDIR}
mergePeaks ${H3K27ac_peaks_norm[@]} -code > ${PEAKDIR}/Allpat_mergePeaks_H3K27ac.CNVnorm.txt
#82942 peaks 
#filter merged peakset
pos2bed.pl ${PEAKDIR}/Allpat_mergePeaks_H3K27ac.CNVnorm.txt > ${TMPDIR}/tmp.1.bed
	$BEDTOOLS intersect -a ${TMPDIR}/tmp.1.bed -b $BLACKLIST_HG38 -v > ${TMPDIR}/tmp.2.bed
	bed2pos.pl ${TMPDIR}/tmp.2.bed > ${TMPDIR}/tmp.1.txt
	filter4Mappability.sh -p ${TMPDIR}/tmp.1.txt -g hg38 -f 0.8 -s 75
	pos2bed.pl ${TMPDIR}/tmp.1.mapScoreFiltered.txt > ${PEAKDIR}/Allpat_mergePeaks_H3K27ac.CNVnorm.filtered.peaks.bed
	bed2pos.pl ${PEAKDIR}/Allpat_mergePeaks_H3K27ac.CNVnorm.filtered.peaks.bed > ${PEAKDIR}/Allpat_mergePeaks_H3K27ac.CNVnorm.filtered.peaks.txt
#filter out X and Y chromosomes
grep -v 'chrX' ${PEAKDIR}/Allpat_mergePeaks_H3K27ac.CNVnorm.filtered.peaks.txt > ${TMPDIR}/H3K27ac.filt.txt
grep -v 'chrY' ${TMPDIR}/H3K27ac.filt.txt > ${PEAKDIR}/Allpat_mergePeaks_H3K27ac.CNVnorm.filtered.XYrem.peaks.txt
wc -l ${PEAKDIR}/Allpat_mergePeaks_H3K27ac.CNVnorm.filtered.XYrem.peaks.txt
##leaves: 77297  peaks
#------------------------------------------------------------------------------------ 





                         ################################################
                         #    annotation of Tagdirs to merged peakset   #
                         ################################################
#------------------------------------------------------------------------------------ 
##define vectors with tagdir names by group
normRad21=()
for i in ${PATIENTS_RAD}; do
    bGname=ChIP_AML_${i}_H3K27ac_CNVnorm
    normRad21+=("$bGname")
done
echo "${normRad21[@]}"

normSA2=()
for i in ${PATIENTS_SA2}; do
    bGname=ChIP_AML_${i}_H3K27ac_CNVnorm
    normSA2+=("$bGname")
done
echo "${normSA2[@]}"

normCTRL=()
for i in ${PATIENTS_CTRL}; do
    bGname=ChIP_AML_${i}_H3K27ac_CNVnorm
    normCTRL+=("$bGname")
done
echo "${normCTRL[@]}"

#------------------------------------------------------------------------------------ 
###annotate all CNVnorm tagdirs on the filtered merged peak set without X and Y
cd ${TAGDIR}
annotatePeaks.pl "${PEAKDIR}/Allpat_mergePeaks_H3K27ac.CNVnorm.filtered.XYrem.peaks.txt" hg38 -size given -d \
${normCTRL[@]} ${normRad21[@]} ${normSA2[@]} \
 -raw -cpu 24 > ${ANNDIR}/Allpat_H3K27ac.XYrem.CNVnorm.peaks.ann.txt
#dimensions
wc -l ${ANNDIR}/Allpat_H3K27ac.XYrem.CNVnorm.peaks.ann.txt #77298

##prepare filtered annotation table for R by removing the first non-count  cols (19)
#check number of samples #40
echo ${normCTRL[@]} ${normRad21[@]} ${normSA2[@]} | wc -w
ncols=$(expr 19 + $(echo ${normCTRL[@]} ${normRad21[@]} ${normSA2[@]} | wc -w))
##removal the first 19 cols
cut -f1,20-${ncols} "${ANNDIR}/Allpat_H3K27ac.XYrem.CNVnorm.peaks.ann.txt" > ${ANNDIR}/Allpat_H3K27ac.XYrem.CNVnorm.peaks.ann.Rinput.txt
#------------------------------------------------------------------------------------ 





            #############################################################################
            #  determine GC and length for cqn normalization with refChr peaks          #
            #############################################################################
#------------------------------------------------------------------------------------ 

	awk '!/^chrM/&&/chr/' ${PEAKDIR}/Allpat_mergePeaks_H3K27ac.CNVnorm.filtered.XYrem.peaks.txt > ${PEAKDIR}/Allpat_mergePeaks_H3K27ac.filtered.XYrem.peaks.refChr.txt #  77291 peak regions
	annotatePeaks.pl ${PEAKDIR}/Allpat_mergePeaks_H3K27ac.filtered.XYrem.peaks.refChr.txt hg38 -size given -noann -nogene -CpG -cpu 12 > "${PEAKDIR}/Allpat_mergePeaks_H3K27ac.filtered.peaks.tmp.CpGann.txt"
	rm "${PEAKDIR}/Allpat_mergePeaks_H3K27ac.filtered.XYrem.peaks.CpGann.txt"
	touch "${PEAKDIR}/Allpat_mergePeaks_H3K27ac.filtered.XYrem.peaks.CpGann.txt"
	echo $'ID\tlength\tgccontent' >> "${PEAKDIR}/Allpat_mergePeaks_H3K27ac.filtered.XYrem.peaks.CpGann.txt"
	awk -v OFS='\t' '{print $1,$4-$3,$9 ; }' <(tail -n+2 "${PEAKDIR}/Allpat_mergePeaks_H3K27ac.filtered.peaks.tmp.CpGann.txt") >> "${PEAKDIR}/Allpat_mergePeaks_H3K27ac.filtered.XYrem.peaks.CpGann.txt"  #length: 77292 
	rm "${PEAKDIR}/Allpat_mergePeaks_H3K27ac.filtered.peaks.tmp.CpGann.txt"

#------------------------------------------------------------------------------------ 

