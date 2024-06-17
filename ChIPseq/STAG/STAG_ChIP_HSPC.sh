
#!/bin/bash
# by Alexander Fischer, APR 2022

###############################################################################
###############################################################################
##                                                                           ##
## Analysis of STAG1 and STAG2 ChIPseq  in CD34 HSPCs with Cohesin KD        ##
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
DIFFPEAKS="${WORKDIR}/diffPeaks/SA1vsSA2"
FIGURESDIR="${WORKDIR}/figures"
FIGURESDIRpyg="${FIGURESDIR}/pygenometracks" #output of pygenometracks
RAD21SCALEDBW="${WORKDIR}/RAD21/scaledBigWigs"
MERGEDBW="${WORKDIR}/mergedBigWigs"
MOTIFDIR="${WORKDIR}/motifs"
DTmatrixdir="${WORKDIR}/deeptoolsMatrix"

mkdir $DIFFPEAKS

##Declare samples
###SA1
CTRL_SA1='ChIP_CD34_14_3_siCtrl_SA1_DQI ChIP_CD34_14_4_Mock_SA1_DQI ChIP_CD34_17_3_siCtrl_SA1_DQI ChIP_CD34_21_4_siCtrl_SA1 ChIP_CD34_21_4_siCtrl_SA1-DQI_Rep1 ChIP_CD34_22_3_siCtrl_SA1 ChIP_CD34_28_6_siCtrl_SA1-DQI_Rep1 '
SA1KD_SA1='ChIP_CD34_14_1_SA1_2259_4094_SA1_DQI ChIP_CD34_17_1_SA1_2259_4094_SA1_DQI ChIP_CD34_21_2_SA1_2259_4094_SA1 ChIP_CD34_27_3_SA1_2259_4094_SA1-DQI ChIP_CD34_28_4_SA1_2259_4094_SA1-DQI '
SA2KD_SA1='ChIP_CD34_14_2_SA2_529_1252_SA1_DQI ChIP_CD34_17_2_SA2_529_1252_SA1_DQI ChIP_CD34_22_2_SA2_529_1252_SA1 ChIP_CD34_21_3_SA2_529_1252_SA1-DQI '
SA1TAGDIRS=$CTRL_SA1$SA1KD_SA1$SA2KD_SA1
###SA2
CTRL_SA2='ChIP_CD34_14_3_siCtrl_SA2 ChIP_CD34_14_4_Mock_SA2 ChIP_CD34_17_3_siCtrl_SA2 ChIP_CD34_18_4_siCtrl_SA2 ChIP_CD34_21_4_siCtrl_SA2 ChIP_CD34_22_3_siCtrl_SA2 ChIP_CD34_28_6_siCtrl_SA2 '
SA1KD_SA2='ChIP_CD34_14_1_SA1_2259_4094_SA2 ChIP_CD34_17_1_SA1_2259_4094_SA2 ChIP_CD34_21_2_SA1_2259_4094_SA2 ChIP_CD34_27_3_SA1_2259_4094_SA2 ChIP_CD34_28_4_SA1_2259_4094_SA2 '
SA2KD_SA2='ChIP_CD34_14_2_SA2_529_1252_SA2 ChIP_CD34_17_2_SA2_529_1252_SA2 ChIP_CD34_20_5_SA2_529_1252_SA2 ChIP_CD34_21_3_SA2_529_1252_SA2 ChIP_CD34_22_2_SA2_529_1252_SA2 ChIP_CD34_28_5_SA2_529_1252_SA2 '
SA2TAGDIRS=$CTRL_SA2$SA1KD_SA2$SA2KD_SA2



                         #####################################
                         #        peak finding for QC        #
                         #####################################
# SA1 PEAKS
for SAMPLE in ${SA1TAGDIRS};do
findPeaks ${SAMPLE} -i ${INPUTDIR}/Input_CD34_merged -style factor -o auto
done
# SA2 PEAKS
for SAMPLE in ${SA2TAGDIRS};do
findPeaks ${SAMPLE} -i ${INPUTDIR}/Input_CD34_merged -style factor -o auto
done
#save peak stats to file
STAGTAGDIRS=$SA1TAGDIRS$SA2TAGDIRS
echo -e "Sample_ID\tTotalPeaks\tIPefficiency\tTotalTags" > $PEAKDIR/STAG1_2peaks.default.txt
for SAMPLE in ${STAGTAGDIRS}; do
    sample_name="${SAMPLE}"
    TotalPeaks=$(grep -w "# total peaks =" $sample_name/peaks.txt | awk '{print $5}')
    IPefficiency=$(grep -w "# Approximate IP efficiency =" $sample_name/peaks.txt | awk '{print $6}')
    TotalTags=$(grep -w "# Total tags =" $sample_name/peaks.txt | awk '{print $5}')
    echo -e "${sample_name}\t${TotalPeaks}\t${IPefficiency}\t${TotalTags}" >> $PEAKDIR/STAG1_2peaks.default.txt
done


                         #####################################
                         #       merged tagDirectories       #
                         #####################################

# all replicates combined by group and antibody
cd ${TAGDIR}
makeTagDirectory ChIP_merged_CD34_CTRL_SA1 -d ${CTRL_SA1}
makeTagDirectory ChIP_merged_CD34_SA1KD_SA1 -d ${SA1KD_SA1}
makeTagDirectory ChIP_merged_CD34_SA2KD_SA1 -d ${SA2KD_SA1}
cd ${TAGDIR}
makeTagDirectory ChIP_merged_CD34_CTRL_SA2 -d ${CTRL_SA2}
makeTagDirectory ChIP_merged_CD34_SA1KD_SA2 -d ${SA1KD_SA2}
makeTagDirectory ChIP_merged_CD34_SA2KD_SA2 -d ${SA2KD_SA2}


# merged bigWigs
cd ${TAGDIR}
declare -a conditions=("CTRL_SA1" "SA1KD_SA1" "SA2KD_SA1" "CTRL_SA2" "SA1KD_SA2" "SA2KD_SA2")
for SAMPLE in "${conditions[@]}";do
makeUCSCfile ChIP_merged_CD34_${SAMPLE} -o ${MERGEDBW}/merged.ChIP_${SAMPLE}.bigWig -fragLength 200 -bigWig ${CHROMSIZES_HG38} -fsize 1e20
done

#------------------------------------------------------------------------------------                         
##look at averge tracks at an exemplary locus (small zoomed in version):
pyGenomeTracks_v3.5.sh --tracks ${FIGURESDIRpyg}/inifiles/Pytracks_SA2KD_CTRL_HSPCs.SA2.average.ini --region chr3:37000000-38000000 --fontSize 40 --dpi 600 -o ${FIGURESDIRpyg}/Pytracks_SA2KD_CTRL_HSPCs.SA2.average.1v1.png
pyGenomeTracks_v3.5.sh --tracks ${FIGURESDIRpyg}/inifiles/Pytracks_SA2KD_CTRL_HSPCs.SA1.average.ini --region chr3:37000000-38000000 --fontSize 40 --dpi 600 -o ${FIGURESDIRpyg}/Pytracks_SA2KD_CTRL_HSPCs.SA1.average.1v1.png

pyGenomeTracks_v3.5.sh --tracks ${FIGURESDIRpyg}/inifiles/Pytracks_SA1KD_CTRL_HSPCs.SA2.average.ini --region chr3:37000000-38000000 --fontSize 40 --dpi 600 -o ${FIGURESDIRpyg}/Pytracks_SA1KD_CTRL_HSPCs.SA2.average.1v1.png
pyGenomeTracks_v3.5.sh --tracks ${FIGURESDIRpyg}/inifiles/Pytracks_SA1KD_CTRL_HSPCs.SA1.average.ini --region chr3:37000000-38000000 --fontSize 40 --dpi 600 -o ${FIGURESDIRpyg}/Pytracks_SA1KD_CTRL_HSPCs.SA1.average.1v1.png
#------------------------------------------------------------------------------------                         


                     ############################################
                     #      peak finding in merged tagDirs      #
                     ############################################
# peak finding
cd ${TAGDIR}
for SAMPLE in "${conditions[@]}";do
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

                     ############################################
                     #             merged peak sets             #
                     ############################################
###SA1
mergePeaks ${PEAKDIR}/CTRL_SA1.filtered.stringentPeaks.txt ${PEAKDIR}/SA1KD_SA1.filtered.stringentPeaks.txt \
${PEAKDIR}/SA2KD_SA1.filtered.stringentPeaks.txt -code > ${PEAKDIR}/CD34_SA1.filtered.stringentPeaks.txt
###SA2
mergePeaks ${PEAKDIR}/CTRL_SA2.filtered.stringentPeaks.txt ${PEAKDIR}/SA1KD_SA2.filtered.stringentPeaks.txt \
${PEAKDIR}/SA2KD_SA2.filtered.stringentPeaks.txt -code > ${PEAKDIR}/CD34_SA2.filtered.stringentPeaks.txt
wc -l ${PEAKDIR}/CD34_SA1.filtered.stringentPeaks.txt
wc -l ${PEAKDIR}/CD34_SA2.filtered.stringentPeaks.txt





                     ########################################################
                     #            annotation to RAD21 postitions            #
                     ########################################################
# annotation of CTRL-HSPC STAG1 and STAG2 ChIP TAGDIRS to all RAD21 positions (merged, filtered peakset from all CD34 conditions)
##these annotation tables are required for the STAG dominance analyses 
cd ${TAGDIR}
annotatePeaks.pl "${PEAKDIR}/CD34_RAD21.filtered.peaks.txt" hg38 -size 250 -d ${CTRL_SA1}${CTRL_SA2} -raw -cpu 12 > ${WORKDIR}/CTRL_SA1S2A2_allCD34_RAD21pos.peaks.ann.txt
tail -n +2 ${WORKDIR}/CTRL_SA1S2A2_allCD34_RAD21pos.peaks.ann.txt | cut -f1,20-33 > ${TMPDIR}/tmp.1.txt
echo $'ID\t14_SA1\t14_SA1\t17_SA1\t21_SA1\t21_SA1\t22_SA1\t28_SA1\t14_SA2\t14_SA2\t17_SA2\t18_SA2\t21_SA2\t22_SA2\t28_SA2' | cat - ${TMPDIR}/tmp.1.txt > ${WORKDIR}/CTRL_SA1S2A2_allCD34_RAD21pos.peaks.ann.Rinput.txt
cut -f 1-6 ${PEAKDIR}/CD34_RAD21.filtered.peaks.txt > ${TMPDIR}/tmp.peaks.txt






                     ########################################################
                     #      annotation to STAG1/STAG2 postitions            #
                     ########################################################

annotatePeaks.pl "${PEAKDIR}/CD34_SA1.filtered.peaks.txt" hg38 -size 250 -d ${SA1TAGDIRS} -raw -cpu 12 > ${WORKDIR}/CD34_SA1.peaks.ann.txt
annotatePeaks.pl "${PEAKDIR}/CD34_SA2.filtered.peaks.txt" hg38 -size 250 -d ${SA2TAGDIRS} -raw -cpu 12 > ${WORKDIR}/CD34_SA2.peaks.ann.txt



                     ########################################################
                     #   Differential Analysis with norm to total tags      #
                     ########################################################
#for global KD target factors: use DESEQ2 with normalization to total tag count -->use HOMER routine

#STAG2 peaks
getDiffExpression.pl ${WORKDIR}/CD34_SA2.peaks.ann.txt -peaks \
CTRL CTRL CTRL CTRL CTRL CTRL CTRL SA1KD SA1KD SA1KD SA1KD SA1KD SA2KD SA2KD SA2KD SA2KD SA2KD SA2KD \
-batch 14 14 17 18 21 22 28 14 17 21 27 28 14 17 20 21 22 28 -rlog -norm2total > ${DIFFPEAKS}/CD34_SA2.peaks.DESEQnorm2total.txt

#STAG1 peaks
getDiffExpression.pl ${WORKDIR}/CD34_SA1.peaks.ann.txt -peaks \
CTRL CTRL CTRL CTRL CTRL CTRL CTRL SA1KD SA1KD SA1KD SA1KD SA1KD SA2KD SA2KD SA2KD SA2KD \
-batch 14 14 17 21 21 22 28 14 17 21 27 28 14 17 22 21 -rlog -norm2total > ${DIFFPEAKS}/CD34_SA1.peaks.DESEQnorm2total.txt
