#!/bin/bash
# by AF, APR 2022

###############################################################################
###############################################################################
##                                                                           ##
##               Analysis of RAD21 ChIPseq in AML patients                   ##
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

## new directories for analysis
PROJECTDIR="${DIR_DATA}/analysis/project_cohesin"
WORKDIR="${PROJECTDIR}/Cohesin_AML/ChIP_analysis"
PEAKDIR="${WORKDIR}/peaks"
ANNDIR="${WORKDIR}/annotation_tables"
DIFFPEAKS="${WORKDIR}/diffPeaks"
FIGURESDIR="${WORKDIR}/figures"
FIGURESDIRpyg="${WORKDIR}/figures/pygenomtracks"
MERGEDBW="${WORKDIR}/mergedBigWigs"
MOTIFDIR="${WORKDIR}/motifs"

mkdir $WORKDIR
mkdir $ANNDIR
mkdir $PEAKDIR
mkdir $DIFFPEAKS
mkdir $FIGURESDIR
mkdir $FIGURESDIRpyg
mkdir $MERGEDBW
mkdir $MOTIFDIR


                         #####################################
                         #        peak finding for QC        #
                         #####################################
#------------------------------------------------------------------------------------                         
# Declare Samples
TARGETS2='Rad21'
PATIENTS_RAD2='Rad21_23039 Rad21_26830 Rad21_27650 Rad21_28420 Rad21_29047 Rad21_29992 Rad21_38455 Rad21_UKR186 Rad21_AML0117 '
PATIENTS_RAD4='23039 26830 27650 28420 29047 29992 38455 UKR186 AML0117'
PATIENTS_SA2_2='SA2_9708 SA2_15365 SA2_15640 SA2_19142 SA2_22464 SA2_22525 SA2_24743 SA2_25458 SA2_27396 SA2_27493 SA2_28041 SA2_29728 SA2_31501 SA2_32635 SA2_BB007 '
PATIENTS_CTRL2='ctr_9873 ctr_10156 ctr_11374 ctr_14094 ctr_16911 ctr_18136 ctr_18367 ctr_18519 ctr_19049 ctr_19284 ctr_19405 ctr_19416 ctr_19504 ctr_20180 ctr_20458 ctr_20899 ctr_21047 ctr_21290 ctr_22023 ctr_22759 ctr_22874 '
ALLPAT_RAD21chips=$PATIENTS_CTRL2$PATIENTS_RAD2$PATIENTS_SA2_2
#------------------------------------------------------------------------------------

# Loop through patients and find peaks with default parameters using the DNA input TAGdir as background
cd ${TAGDIR}
#CTRL AML
for patID in ${PATIENTS_CTRL2}; do
findPeaks ${TAGDIR}/ChIP_AML_${patID}_Rad21 -i ${INPUTDIR}/ChIP_IN_AML_${patID} -style factor -o auto
done
#SA2 AML
for patID in ${PATIENTS_SA2_2};do
findPeaks ${TAGDIR}/ChIP_AML_${patID}_Rad21 -i ${INPUTDIR}/ChIP_IN_AML_${patID} -style factor -o auto
done
#RAD21 AML
for patID in ${PATIENTS_RAD4};do
findPeaks ${TAGDIR}/ChIP_AML_Rad21_${patID}_Rad21 -i ${INPUTDIR}/ChIP_IN_AML_RAD21_${patID} -style factor -o auto
done
#------------------------------------------------------------------------------------
# show values
cd ${TAGDIR}
for SAMPLE in ${ALLPAT_RAD21chips};do
paste <(echo -e "# ChIP_AML_${SAMPLE}_Rad21") <(grep -w "# total peaks =" ChIP_AML_${SAMPLE}_Rad21/peaks.txt) <(grep -w "# Approximate IP efficiency" ChIP_AML_${SAMPLE}_Rad21/peaks.txt) <(grep -w "# Total tags =" ChIP_AML_${SAMPLE}_Rad21/peaks.txt)
done
#save to file
echo -e "Sample_ID\tTotalPeaks\tIPefficiency\tTotalTags" > $PEAKDIR/RAD21peaks.default.txt
for SAMPLE in ${ALLPAT_RAD21chips}; do
    sample_name="ChIP_AML_${SAMPLE}_Rad21"
    TotalPeaks=$(grep -w "# total peaks =" $sample_name/peaks.txt | awk '{print $5}')
    IPefficiency=$(grep -w "# Approximate IP efficiency =" $sample_name/peaks.txt | awk '{print $6}')
    TotalTags=$(grep -w "# Total tags =" $sample_name/peaks.txt | awk '{print $5}')
    echo -e "${sample_name}\t${TotalPeaks}\t${IPefficiency}\t${TotalTags}" >> $PEAKDIR/RAD21peaks.default.txt
done
#------------------------------------------------------------------------------------



                         #####################################
                         #        normalize TAGDIRS by CNV   #
                         #####################################

#------------------------------------------------------------------------------------
#loop through targets and patients
#CTRL AML
for patID in ${PATIENTS_CTRL2};
do
normalizeTagDirByCopyNumber.pl ${TAGDIR}/ChIP_AML_${patID}_Rad21 -cnv ${CNVDIR}/ChIP_IN_AML_${patID}.sam_CNVs
done
#STAG2 AML
for patID in ${PATIENTS_SA2_2};do
normalizeTagDirByCopyNumber.pl ${TAGDIR}/ChIP_AML_${patID}_Rad21 -cnv ${CNVDIR}/ChIP_IN_AML_${patID}.sam_CNVs
done
#RAD21 AML
for patID in ${PATIENTS_RAD4};do
normalizeTagDirByCopyNumber.pl ${TAGDIR}/ChIP_AML_Rad21_${patID}_Rad21 -cnv ${CNVDIR}/ChIP_IN_AML_RAD21_${patID}.sam_CNVs
done
#------------------------------------------------------------------------------------


                    #############################################
                    #       peak finding CNV norm tagdirs       #
                    #############################################
#------------------------------------------------------------------------------------
#loop through  patients and find peaks with specified stringency settings
#CTRL AML
for patID in ${PATIENTS_CTRL2};do
findPeaks ${TAGDIR}/ChIP_AML_${patID}_Rad21_CNVnorm -i ${INPUTDIR}/ChIP_IN_AML_${patID} -style factor -o auto
findPeaks ${TAGDIR}/ChIP_AML_${patID}_Rad21_CNVnorm -i ${INPUTDIR}/ChIP_IN_AML_${patID} -style factor -tbp 1 -fdr 0.000001 -ntagThreshold 10 -o ${TAGDIR}/ChIP_AML_${patID}_Rad21_CNVnorm/stringentPeaks.txt
done
#SA2 AML
for patID in ${PATIENTS_SA2_2};do
findPeaks ${TAGDIR}/ChIP_AML_${patID}_Rad21_CNVnorm -i ${INPUTDIR}/ChIP_IN_AML_${patID} -style factor -o auto
findPeaks ${TAGDIR}/ChIP_AML_${patID}_Rad21_CNVnorm -i ${INPUTDIR}/ChIP_IN_AML_${patID} -style factor -tbp 1 -fdr 0.000001 -ntagThreshold 10 -o ${TAGDIR}/ChIP_AML_${patID}_Rad21_CNVnorm/stringentPeaks.txt
done
#RAD21 AML
for patID in ${PATIENTS_RAD4};do
findPeaks ${TAGDIR}/ChIP_AML_Rad21_${patID}_Rad21_CNVnorm -i ${INPUTDIR}/ChIP_IN_AML_RAD21_${patID} -style factor -o auto
findPeaks ${TAGDIR}/ChIP_AML_Rad21_${patID}_Rad21_CNVnorm -i ${INPUTDIR}/ChIP_IN_AML_RAD21_${patID} -style factor -tbp 1 -fdr 0.000001 -ntagThreshold 10 -o ${TAGDIR}/ChIP_AML_Rad21_${patID}_Rad21_CNVnorm/stringentPeaks.txt
done
#------------------------------------------------------------------------------------
# show peak finding values for all RAD21 chips after CNVnorm
cd ${TAGDIR}
##standard calling parmeters
for SAMPLE in ${ALLPAT_RAD21chips};do
paste <(echo -e "# ChIP_AML_${SAMPLE}_Rad21_CNVnorm") <(grep -w "# total peaks =" ChIP_AML_${SAMPLE}_Rad21_CNVnorm/peaks.txt) <(grep -w "# Approximate IP efficiency" ChIP_AML_${SAMPLE}_Rad21_CNVnorm/peaks.txt) <(grep -w "# Total tags =" ChIP_AML_${SAMPLE}_Rad21_CNVnorm/peaks.txt)
done
##stringent calling parmeters
for SAMPLE in ${ALLPAT_RAD21chips};do
paste <(echo -e "# ChIP_AML_${SAMPLE}_Rad21_CNVnorm") <(grep -w "# total peaks =" ChIP_AML_${SAMPLE}_Rad21_CNVnorm/stringentPeaks.txt) <(grep -w "# Approximate IP efficiency" ChIP_AML_${SAMPLE}_Rad21_CNVnorm/stringentPeaks.txt) <(grep -w "# Total tags =" ChIP_AML_${SAMPLE}_Rad21_CNVnorm/stringentPeaks.txt)
done

#save to files
echo -e "Sample_ID\tTotalPeaks\tIPefficiency\tTotalTags" > $PEAKDIR/RAD21peaks.CNVnorm.txt
for SAMPLE in ${ALLPAT_RAD21chips}; do
    sample_name="ChIP_AML_${SAMPLE}_Rad21_CNVnorm"
    TotalPeaks=$(grep -w "# total peaks =" $sample_name/peaks.txt | awk '{print $5}')
    IPefficiency=$(grep -w "# Approximate IP efficiency =" $sample_name/peaks.txt | awk '{print $6}')
    TotalTags=$(grep -w "# Total tags =" $sample_name/peaks.txt | awk '{print $5}')
    echo -e "${sample_name}\t${TotalPeaks}\t${IPefficiency}\t${TotalTags}" >> $PEAKDIR/RAD21peaks.CNVnorm.txt
done

echo -e "Sample_ID\tTotalPeaks\tIPefficiency\tTotalTags" > $PEAKDIR/RAD21peaks.strigent.CNVnorm.txt
for SAMPLE in ${ALLPAT_RAD21chips}; do
    sample_name="ChIP_AML_${SAMPLE}_Rad21_CNVnorm"
    TotalPeaks=$(grep -w "# total peaks =" $sample_name/stringentPeaks.txt | awk '{print $5}')
    IPefficiency=$(grep -w "# Approximate IP efficiency =" $sample_name/stringentPeaks.txt | awk '{print $6}')
    TotalTags=$(grep -w "# Total tags =" $sample_name/stringentPeaks.txt | awk '{print $5}')
    echo -e "${sample_name}\t${TotalPeaks}\t${IPefficiency}\t${TotalTags}" >> $PEAKDIR/RAD21peaks.strigent.CNVnorm.txt
done
#------------------------------------------------------------------------------------

                     ######################################
                     #CNV norm BIGwigs and average BIGwigs#
                     ######################################
#------------------------------------------------------------------------------------
# CNV norm. bw file: create bigwig from corrected tagdir
for SAMPLE in ${ALLPAT_RAD21chips};do
makeUCSCfile ${TAGDIR}/ChIP_AML_${SAMPLE}_Rad21*_CNVnorm -o ${BWDIR}/ChIP_AML_${SAMPLE}_Rad21_CNVnorm.bigWig -fragLength 200 -bigWig ${CHROMSIZES_HG38} -fsize 1e20
done
#------------------------------------------------------------------------------------
# average BW file by condition (create from CNVnorm BWs)
##define vectors with bigwig names by group
normRad21Rad21bw=()
for i in ${PATIENTS_RAD2}; do
    bGname=ChIP_AML_${i}_Rad21_CNVnorm.bigWig
    normRad21Rad21bw+=("$bGname")
done
echo "${normRad21Rad21bw[@]}"

normSA2Rad21bw=()
for i in ${PATIENTS_SA2_2}; do
    bGname=ChIP_AML_${i}_Rad21_CNVnorm.bigWig
    normSA2Rad21bw+=("$bGname")
done
echo "${normSA2Rad21bw[@]}"

normCTRLRad21bw=()
for i in ${PATIENTS_CTRL2}; do
    bGname=ChIP_AML_${i}_Rad21_CNVnorm.bigWig
    normCTRLRad21bw+=("$bGname")
done
echo "${normCTRLRad21bw[@]}"
#------------------------------------------------------------------------------------
##run AverageBigwig.pl script and create bedgraphs in addition
cd ${BWDIR}
#CTRL AML
myAverageBigWig.pl -bw ${normCTRLRad21bw[@]} -chr ${CHROMSIZES_HG38} -o ${MERGEDBW}/ave.AML_CTRL_RAD21_CNVnorm.bigwig
bigWigToBedGraph ${MERGEDBW}/ave.AML_CTRL_RAD21_CNVnorm.bigwig ${MERGEDBW}/ave.AML_CTRL_RAD21.CNVnorm.bedGraph
#STAG2mut AML
myAverageBigWig.pl -bw ${normSA2Rad21bw[@]} -chr ${CHROMSIZES_HG38} -o ${MERGEDBW}/ave.AML_SA2mut_RAD21_CNVnorm.bigwig
bigWigToBedGraph ${MERGEDBW}/ave.AML_SA2mut_RAD21_CNVnorm.bigwig ${MERGEDBW}/ave.AML_SA2mut_RAD21.CNVnorm.bedGraph
#CTRL AML
myAverageBigWig.pl -bw ${normCTRLRad21bw[@]} -chr ${CHROMSIZES_HG38} -o ${MERGEDBW}/ave.AML_CTRL_RAD21_CNVnorm.bigwig
bigWigToBedGraph ${MERGEDBW}/ave.AML_CTRL_RAD21_CNVnorm.bigwig ${MERGEDBW}/ave.AML_CTRL_RAD21.CNVnorm.bedGraph
#RAD21mut AML
myAverageBigWig.pl -bw ${normRad21Rad21bw[@]} -chr ${CHROMSIZES_HG38} -o ${MERGEDBW}/ave.AML_Rad21mut_RAD21_CNVnorm.bigwig
bigWigToBedGraph ${MERGEDBW}/ave.AML_Rad21mut_RAD21_CNVnorm.bigwig ${MERGEDBW}/ave.AML_Rad21mut_RAD21.CNVnorm.bedGraph
#------------------------------------------------------------------------------------
                       
##look at individual patient tracks at an exemplary locus:
##pygenometracks using ini file with all bigwigs listed ##SA2mut vs RAD21mut vs CTRL all pat
pyGenomeTracks_v3.5.sh --tracks ${FIGURESDIRpyg}/inifiles/Pytracks_SA2_RAD21mutvsCTRL.RAD21.allPat.ini --region chr3:36750000-38750000 --dpi 600 --fontSize 28  -o ${FIGURESDIRpyg}/Pytracks_SA2_RAD21mutvsCTRL.RAD21.allPat.ITGA9.png
#------------------------------------------------------------------------------------  

                     ################################################
                     #merged Peaksets based on individual peakfiles #
                     ################################################
#create merged peak set based on peaks called in individual patients for standard and strigent peak sets
#------------------------------------------------------------------------------------
#regular peaks
##define vector with peak files
RAD21peaks_norm=()
for i in ${ALLPAT_RAD21chips}; do
    bGname=ChIP_AML_${i}_Rad21*CNVnorm/peaks.txt
    RAD21peaks_norm+=("$bGname")
done
echo "${RAD21peaks_norm[@]}"
##merge and filter peaks
cd ${TAGDIR}
mergePeaks ${RAD21peaks_norm[@]} -code > ${PEAKDIR}/Allpat_mergePeaks_RAD21.txt #merging
###filtering for blacklisted regions and low mappability
pos2bed.pl ${PEAKDIR}/Allpat_mergePeaks_RAD21.txt > ${TMPDIR}/tmp.1.bed
	$BEDTOOLS intersect -a ${TMPDIR}/tmp.1.bed -b $BLACKLIST_HG38 -v > ${TMPDIR}/tmp.2.bed
	bed2pos.pl ${TMPDIR}/tmp.2.bed > ${TMPDIR}/tmp.1.txt
	filter4Mappability.sh -p ${TMPDIR}/tmp.1.txt -g hg38 -f 0.8 -s 75
	pos2bed.pl ${TMPDIR}/tmp.1.mapScoreFiltered.txt > ${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.bed
	bed2pos.pl ${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.bed > ${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.txt
#------------------------------------------------------------------------------------
#stringent peaks
RAD21strpeaks_norm=()
for i in ${ALLPAT_RAD21chips}; do
    bGname=ChIP_AML_${i}_Rad21*CNVnorm/stringentPeaks.txt
    RAD21strpeaks_norm+=("$bGname")
done
echo "${RAD21strpeaks_norm[@]}"
##merge and filter peaks
mergePeaks ${RAD21strpeaks_norm[@]} -code > ${PEAKDIR}/Allpat_mergePeaks_RAD21_stringent.txt
###filtering for blacklisted regions and low mappability
pos2bed.pl ${PEAKDIR}/Allpat_mergePeaks_RAD21_stringent.txt > ${TMPDIR}/tmp.1.bed
	$BEDTOOLS intersect -a ${TMPDIR}/tmp.1.bed -b $BLACKLIST_HG38 -v > ${TMPDIR}/tmp.2.bed
	bed2pos.pl ${TMPDIR}/tmp.2.bed > ${TMPDIR}/tmp.1.txt
	filter4Mappability.sh -p ${TMPDIR}/tmp.1.txt -g hg38 -f 0.8 -s 75
	pos2bed.pl ${TMPDIR}/tmp.1.mapScoreFiltered.txt > ${PEAKDIR}/Allpat_mergePeaks_RAD21_stringent.filtered.peaks.bed
	bed2pos.pl ${PEAKDIR}/Allpat_mergePeaks_RAD21_stringent.filtered.peaks.bed > ${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.stringent.txt
#------------------------------------------------------------------------------------



                     ####################################################
                     #     Global RAD21 coverage ghist plots            #
                     ####################################################
##visualize average RAD21 coverage across all RAD21 peaks detectecd
#------------------------------------------------------------------------------------
mkdir ${FIGURESDIR}/ghist/
mkdir ${FIGURESDIR}/ghist/RAD21

patgroups="CTRL SA2mut Rad21mut"
_DATE=$(date +%s)
for group in ${patgroups};do
	cat >"${TMPDIR}/ghist.${group}.${_DATE}.sh" <<EOF
#!/bin/bash
#setting homer environment
export PATH=/misc/software/package/RBioC/3.4.3/bin:/misc/software/package/perl/perl-5.26.1/bin:/misc/software/ngs/samtools/samtools-1.6/bin:/misc/software/ngs/homer/v4.11/bin:${PATH}
export PATH
cd /loctmp
#use CNVnormed averaged bedgraphs
annotatePeaks.pl ${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_${group}_RAD21.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/RAD21/RAD21_${group}.RAD21pos.all.ghist.txt"
annotatePeaks.pl ${PEAKDIR}/Allpat_mergePeaks_RAD21_stringent.filtered.peaks.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_${group}_RAD21.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/RAD21/RAD21_${group}.RAD21pos.stringent.all.ghist.txt"
EOF
	chmod 750 "${TMPDIR}/ghist.${group}.${_DATE}.sh"
	echo "plotting ghists for ${group} "
	screen -dm -S ${group}${_DATE} bash -c "bash ${TMPDIR}/ghist.${group}.${_DATE}.sh"
done

#------------------------------------------------------------------------------------
#plot using plotHIST.sh script (R-based plot)
plotHIST.sh -g "RAD21_CTRL.RAD21pos.all.ghist.txt RAD21_SA2mut.RAD21pos.all.ghist.txt" \
-s "CTRL-AML STAG2-mut" -c "firebrick3 seagreen3" -x 1000 -y "0 15" -d ${FIGURESDIR}/ghist/RAD21 -n RAD21.all.CNVnorm.SA2mut.vs.CTRL
plotHIST.sh -g "RAD21_CTRL.RAD21pos.all.ghist.txt RAD21_RAD21mut.RAD21pos.all.ghist.txt" \
-s "CTRL-AML RAD21-mut" -c "firebrick3 mediumvioletred" -x 1000 -y "0 15" -d ${FIGURESDIR}/ghist/RAD21 -n RAD21.all.RAD21mut.vs.CTRL

plotHIST.sh -g "RAD21_CTRL.RAD21pos.stringent.all.ghist.txt RAD21_SA2mut.RAD21pos.stringent.all.ghist.txt" \
-s "CTRL-AML STAG2-mut" -c "firebrick3 seagreen3" -x 1000 -y "0 15" -d ${FIGURESDIR}/ghist/RAD21 -n RAD21.all.stringent.CNVnorm.SA2mut.vs.CTRL
plotHIST.sh -g "RAD21_CTRL.RAD21pos.all.stringent.ghist.txt RAD21_RAD21mut.RAD21pos.all.stringent.ghist.txt" \
-s "CTRL-AML RAD21-mut" -c "firebrick3 mediumvioletred" -x 1000 -y "0 15" -d ${FIGURESDIR}/ghist/RAD21 -n RAD21.all.stringent.RAD21mut.vs.CTRL

#------------------------------------------------------------------------------------




                     ####################################################
                     #     Annotation of tagdirs to peak sets           #
                     ####################################################
##define vectors with tagdir names by group
normRad21Rad21=()
for i in ${PATIENTS_RAD2}; do
    bGname=ChIP_AML_${i}_Rad21_CNVnorm
    normRad21Rad21+=("$bGname")
done
echo "${normRad21Rad21[@]}"

normSA2Rad21=()
for i in ${PATIENTS_SA2_2}; do
    bGname=ChIP_AML_${i}_Rad21_CNVnorm
    normSA2Rad21+=("$bGname")
done
echo "${normSA2Rad21[@]}"

normCTRLRad21=()
for i in ${PATIENTS_CTRL2}; do
    bGname=ChIP_AML_${i}_Rad21_CNVnorm
    normCTRLRad21+=("$bGname")
done
echo "${normCTRLRad21[@]}"

###annotate all CNVnorm tagdirs on the filtered merged peak set, filter out X and Y! These Tables will be used for differential peak analysis
#------------------------------------------------------------------------------------
cd ${TAGDIR}
#regular peak set
annotatePeaks.pl "${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.txt" hg38 -size 250 -d \
${normCTRLRad21[@]} ${normRad21Rad21[@]} ${normSA2Rad21[@]} \
-raw -cpu 24 > ${ANNDIR}/Allpat_RAD21.peaks.ann.txt
#remove sex chromosomes
grep -v 'chrX' ${ANNDIR}/Allpat_RAD21.peaks.ann.txt > ${TMPDIR}/tmp.filt.txt
grep -v 'chrY' ${TMPDIR}/tmp.filt.txt > ${ANNDIR}/Allpat_RAD21.peaks.ann.filt.txt
rm ${TMPDIR}/tmp.filt*
#prepare for R input for downstream analyses
tail -n +2 ${ANNDIR}/Allpat_RAD21.peaks.ann.filt.txt | cut -f1,20-100 > ${ANNDIR}/Allpat_RAD21.peaks.ann.filt.Rinput.txt
cut -f 1-6 ${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.txt > ${ANNDIR}/RAD21.filtered.peaks.Rinput.txt

#------------------------------------------------------------------------------------
#stringent peak set
annotatePeaks.pl "${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.stringent.txt" hg38 -size 250 -d \
${normCTRLRad21[@]} ${normRad21Rad21[@]} ${normSA2Rad21[@]} \
-raw -cpu 24 > ${ANNDIR}/Allpat_RAD21.peaks.stringent.ann.txt
#remove sex chromosomes
grep -v 'chrX' ${ANNDIR}/Allpat_RAD21.peaks.stringent.ann.txt > ${TMPDIR}/tmp.filt.txt
grep -v 'chrY' ${TMPDIR}/tmp.filt.txt > ${ANNDIR}/Allpat_RAD21.peaks.stringent.ann.filt.txt
rm ${TMPDIR}/tmp.filt*
#prepare for R input for downstream analyses
tail -n +2 ${ANNDIR}/Allpat_RAD21.peaks.stringent.ann.filt.txt | cut -f1,20-100 > ${ANNDIR}/Allpat_RAD21.peaks.stringent.ann.filt.Rinput.txt
cut -f 1-6 ${PEAKDIR}/Allpat_mergePeaks_RAD21_stringent.filtered.peaks.txt > ${ANNDIR}/RAD21_stringent.filtered.peaks.Rinput.txt
#------------------------------------------------------------------------------------

