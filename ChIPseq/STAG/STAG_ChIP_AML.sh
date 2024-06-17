#!/bin/bash
# by Alexander Fischer, July 2021

###############################################################################
###############################################################################
##                                                                           ##
##       Analysis of STAG2 and STAG1 ChIPseq in AML patients                 ##
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

##  directories for analysis
PROJECTDIR="${DIR_DATA}/analysis/project_cohesin"
WORKDIR="${PROJECTDIR}/Cohesin_AML/ChIP_analysis"
PEAKDIR="${WORKDIR}/peaks"
ANNDIR="${WORKDIR}/annotation_tables"
DIFFPEAKS="${WORKDIR}/diffPeaks"
FIGURESDIR="${WORKDIR}/figures"
FIGURESDIRpyg="${FIGURESDIR}/pygenometracks" #output for pygenometracks figures
MERGEDBW="${WORKDIR}/mergedBigWigs"
MOTIFDIR="${WORKDIR}/motifs"

mkdir $FIGURESDIRpyg 
mkdir ${FIGURESDIRpyg}/inifiles/ #copy ini files for pygenometracks here!

## Directory for HSPC peak comparisons:
PEAKDIRCD34="${PROJECTDIR}/CD34/ChIP_KD_analysis/Cohesin_CTCF_MED12/peaks"


                        #######################################
                        #######################################
                        ##    STAG2 Chipseq Analysis         ##
                        #######################################
                        #######################################   


                         ###########################################
                         #  STAG2 (SA2) peak finding for QC        #
                         ###########################################
#------------------------------------------------------------------------------------                         
#Declare Samples
PATIENTS_SA2_SA2chips='SA2_9708 SA2_24743  SA2_29728 SA2_27396 SA2_15365 SA2_15640 SA2_19142 SA2_22464 SA2_22525 SA2_25458 SA2_27493 SA2_28041 SA2_31501 SA2_32635 SA2_BB007 '
PATIENTS_RAD21_SA2chips='Rad21_38455 Rad21_26830 Rad21_23039 Rad21_UKR186 RAD21_29047 RAD21_27650 RAD21_28420 RAD21_29992 '
PATIENTS_CTRL_SA2chips='ctr_18136 ctr_16911 ctr_19416 ctr_21047 ctr_18519 ctr_19405 ctr_21290 ctr_9873 ctr_11374 ctr_14094 ctr_18367 ctr_19049 ctr_19284 ctr_20180 ctr_20458 ctr_20899 ctr_22023 '
ALLPAT_SA2chips=$PATIENTS_SA2_SA1chips$PATIENTS_RAD21_SA2chips$PATIENTS_CTRL_SA1chips

#loop through targets and patients
for patID in ${ALLPAT_SA2chips};do
findPeaks ${TAGDIR}/ChIP_AML_${patID}_SA2 -i ${INPUTDIR}/ChIP_IN_AML_${patID} -style factor -o auto
done

#save peak finding stats to file
echo -e "Sample_ID\tTotalPeaks\tIPefficiency\tTotalTags" > $PEAKDIR/STAG2peaks.default.txt
for SAMPLE in ${ALLPAT_SA2chips}; do
    sample_name="ChIP_AML_${SAMPLE}_SA2"
    TotalPeaks=$(grep -w "# total peaks =" $sample_name/peaks.txt | awk '{print $5}')
    IPefficiency=$(grep -w "# Approximate IP efficiency =" $sample_name/peaks.txt | awk '{print $6}')
    TotalTags=$(grep -w "# Total tags =" $sample_name/peaks.txt | awk '{print $5}')
    echo -e "${sample_name}\t${TotalPeaks}\t${IPefficiency}\t${TotalTags}" >> $PEAKDIR/SA2peaks.default.txt
done

                         #####################################
                         #    normalize SA2 TAGDIRS by CNV   #
                         #####################################

#Usage: normalizeTagDirByCopyNumber.pl <tagDir> -cnv <CNV file>
for patID in ${ALLPAT_SA2chips};do
normalizeTagDirByCopyNumber.pl ${TAGDIR}/ChIP_AML_${patID}_SA2 -cnv ${CNVDIR}/ChIP_IN_AML_${patID}.sam_CNVs
done

#------------------------------------------------------------------------------------                         

                    #############################################
                    #  SA2   peak finding  in CNV norm tagdirs  #
                    #############################################

#loop through patients
for patID in ${ALLPAT_SA1chips};do
findPeaks ${TAGDIR}/ChIP_AML_${patID}_SA2_CNVnorm -i ${INPUTDIR}/ChIP_IN_AML_${patID} -style factor -o auto
done
done

#save to files
echo -e "Sample_ID\tTotalPeaks\tIPefficiency\tTotalTags" > $PEAKDIR/STAG2peaks.CNVnorm.txt
for SAMPLE in ${ALLPAT_SA1chips}; do
    sample_name="ChIP_AML_${SAMPLE}_SA2_CNVnorm"
    TotalPeaks=$(grep -w "# total peaks =" $sample_name/peaks.txt | awk '{print $5}')
    IPefficiency=$(grep -w "# Approximate IP efficiency =" $sample_name/peaks.txt | awk '{print $6}')
    TotalTags=$(grep -w "# Total tags =" $sample_name/peaks.txt | awk '{print $5}')
    echo -e "${sample_name}\t${TotalPeaks}\t${IPefficiency}\t${TotalTags}" >> $PEAKDIR/STAG2peaks.CNVnorm.txt
done
#------------------------------------------------------------------------------------                         

                     ###########################################
                     #SA2 CNV norm BIGwigs and average BIGwigs #
                     ###########################################

#generate CNV norm bw file
for SAMPLE in ${ALLPAT_SA2chips}
do
makeUCSCfile ${TAGDIR}/ChIP_AML_${SAMPLE}_SA2*_CNVnorm -o ${BWDIR}/ChIP_AML_${SAMPLE}_SA2_CNVnorm.bigWig -fragLength 200 -bigWig ${CHROMSIZES_HG38} -fsize 1e20
done
#------------------------------------------------------------------------------------                         
#assemble average BW file name by condition
normSA2SA2bw=()
for i in ${PATIENTS_SA2_SA2chips}; do
    bGname=ChIP_AML_${i}_SA2_CNVnorm.bigWig
    normSA2SA2bw+=("$bGname")
done
echo "${normSA2SA2bw[@]}"

normRAD21SA2bw=()
for i in ${PATIENTS_RAD21_SA2chips}; do
    bGname=ChIP_AML_${i}_SA2_CNVnorm.bigWig
    normRAD21SA2bw+=("$bGname")
done
echo "${normRAD21SA2bw[@]}"

normCTRLSA2bw=()
for i in ${PATIENTS_CTRL_SA2chips}; do
    bGname=ChIP_AML_${i}_SA2_CNVnorm.bigWig
    normCTRLSA2bw+=("$bGname")
done
echo "${normCTRLSA2bw[@]}"

#------------------------------------------------------------------------------------                         
##run AverageBigwig.pl script and create bedgraphs in addition
cd ${BWDIR}
#CTRL AML
myAverageBigWig.pl -bw ${normCTRLSA2bw[@]} -chr ${CHROMSIZES_HG38} -o ${MERGEDBW}/ave.AML_CTRL_SA2_CNVnorm.bigwig
bigWigToBedGraph ${MERGEDBW}/ave.AML_CTRL_SA2_CNVnorm.bigwig ${MERGEDBW}/ave.AML_CTRL_SA2.CNVnorm.bedGraph

#STAG2mut AML
myAverageBigWig.pl -bw ${normSA2SA2bw[@]} -chr ${CHROMSIZES_HG38} -o ${MERGEDBW}/ave.AML_SA2mut_SA2_CNVnorm.bigwig
bigWigToBedGraph ${MERGEDBW}/ave.AML_SA2mut_SA2_CNVnorm.bigwig ${MERGEDBW}/ave.AML_SA2mut_SA2.CNVnorm.bedGraph

#RAD21mut AML
myAverageBigWig.pl -bw ${normRAD21SA2bw[@]} -chr ${CHROMSIZES_HG38} -o ${MERGEDBW}/ave.AML_RAD21mut_SA2_CNVnorm.bigwig
bigWigToBedGraph ${MERGEDBW}/ave.AML_RAD21mut_SA2_CNVnorm.bigwig ${MERGEDBW}/ave.AML_RAD21mut_SA2.CNVnorm.bedGraph


#------------------------------------------------------------------------------------                         


                     #####################################################
                     #merged STAG2 Peakset based on individual peakfiles #
                     #####################################################
#------------------------------------------------------------------------------------
##define vector with peak files
SA2peaks_norm=()
for i in ${ALLPAT_SA2chips}; do
    bGname=ChIP_AML_${i}_SA1*CNVnorm/peaks.txt
    SA2peaks_norm+=("$bGname")
done
echo "${SA2peaks_norm[@]}"
##merge and filter peaks
cd ${TAGDIR}
mergePeaks ${SA2peaks_norm[@]} -code > ${PEAKDIR}/Allpat_mergePeaks_SA2.txt
###filtering for blacklisted regions and low mappability
pos2bed.pl ${PEAKDIR}/Allpat_mergePeaks_SA2.txt > ${TMPDIR}/tmp.1.bed
	$BEDTOOLS intersect -a ${TMPDIR}/tmp.1.bed -b $BLACKLIST_HG38 -v > ${TMPDIR}/tmp.2.bed
	bed2pos.pl ${TMPDIR}/tmp.2.bed > ${TMPDIR}/tmp.1.txt
	filter4Mappability.sh -p ${TMPDIR}/tmp.1.txt -g hg38 -f 0.8 -s 75
	pos2bed.pl ${TMPDIR}/tmp.1.mapScoreFiltered.txt > ${PEAKDIR}/Allpat_merAllpat_mergePeaks_SA2gePeaks_SA1.filtered.peaks.bed
	bed2pos.pl ${PEAKDIR}/Allpat_mergePeaks_SA2.filtered.peaks.bed > ${PEAKDIR}/Allpat_mergePeaks_SA2.filtered.peaks.txt    
#------------------------------------------------------------------------------------  
##look at individual patient tracks at an exemplary locus:
mkdir ${FIGURESDIRpyg}/SA2tracks/
##pygenometracks using ini file with all bigwigs listed
pyGenomeTracks_v3.5.sh --tracks ${FIGURESDIRpyg}/inifiles/Pytracks_SA2_RAD21mutvsCTRL.SA2.allPat.ini --region chr3:36550000-38550000 --dpi 600 --fontSize 28 -o ${FIGURESDIRpyg}/SA2tracks/Pytracks_SA2_RAD21mutvsCTRL.SA2.allPat.ITGA9.png
#------------------------------------------------------------------------------------  




                         ###############################################################
                         #   Annotation of STAG2 Tagdirectories to STAG2 peak set      #
                         ###############################################################
#------------------------------------------------------------------------------------ 
##define vectors with tagdir names by group
normSA2=()
for i in ${PATIENTS_SA2}; do
    bGname=ChIP_AML_${i}_SA2_CNVnorm
    normSA2+=("$bGname")
done
echo "${normSA2[@]}"

normCTRL=()
for i in ${PATIENTS_CTRL}; do
    bGname=ChIP_AML_${i}_SA2_CNVnorm
    normCTRL+=("$bGname")
done
echo "${normCTRL[@]}"

###annotate all CNVnorm tagdirs on the filtered merged peak set, filter out X and Y! These Tables will be used for differential peak analysis
#------------------------------------------------------------------------------------
annotatePeaks.pl "${PEAKDIR}/Allpat_mergePeaks_SA2.filtered.peaks.txt" hg38 -size 250 -d \
${normCTRL[@]} ${normSA2[@]} \
-raw -cpu 24 > ${WORKDIR}/Allpat_mergePeaks_SA2.filtered.peaks.ann.txt
#remove sex chromosomes
grep -v 'chrX' ${WORKDIR}/Allpat_mergePeaks_SA2.filtered.peaks.ann.txt > ${TMPDIR}/SA2_2.filt.txt
grep -v 'chrY' ${TMPDIR}/SA2_2.filt.txt > ${WORKDIR}/Allpat_mergePeaks_SA2.filtered.peaks.ann.filt.txt
rm ${TMPDIR}/SA2_2.filt*
#------------------------------------------------------------------------------------  

                         ###############################################################
                         #       STAG2 diffPeak analysis DESEQ: SA2mut vs CTRL         #
                         ###############################################################
#------------------------------------------------------------------------------------  

#Run with HOMER getDiffExpression wiht norm to total tag count
getDiffExpression.pl ${WORKDIR}/Allpat_mergePeaks_SA2.filtered.peaks.ann.filt.txt -peaks \
ctr ctr ctr ctr ctr ctr ctr ctr ctr ctr ctr ctr ctr ctr ctr ctr ctr SA2 SA2 SA2 SA2 SA2 SA2 SA2 SA2 SA2 SA2 SA2 SA2 SA2 SA2 SA2 \
-rlog -norm2total > ${DIFFPEAKS}/Allpat_SA2.peaks.DESEQnorm2total.txt
#------------------------------------------------------------------------------------  
#------------------------------------------------------------------------------------                         



                        #######################################
                        #######################################
                        ##    STAG1 Chipseq Analysis         ##
                        #######################################
                        #######################################   


                         ###########################################
                         #  STAG1 (SA1) peak finding for QC        #
                         ###########################################
#------------------------------------------------------------------------------------                         

#Declare Samples
PATIENTS_SA2_SA1chips='SA2_9708 SA2_15365 SA2_24743 SA2_27396 SA2_29728 '
PATIENTS_RAD21_SA1chips='Rad21_23039 Rad21_38455 RAD21_27650 RAD21_28420'
PATIENTS_CTRL_SA1chips='ctr_16911 ctr_18136 ctr_18519 ctr_19405 ctr_19416 ctr_21047 '
ALLPAT_SA1chips=$PATIENTS_SA2_SA1chips$PATIENTS_RAD21_SA1chips$PATIENTS_CTRL_SA1chips

#loop through targets and patients
for patID in ${ALLPAT_SA1chips};do
findPeaks ${TAGDIR}/ChIP_AML_${patID}_SA1 -i ${INPUTDIR}/ChIP_IN_AML_${patID} -style factor -o auto
done


#save peak finding stats to file
echo -e "Sample_ID\tTotalPeaks\tIPefficiency\tTotalTags" > $PEAKDIR/STAG1peaks.default.txt
for SAMPLE in ${ALLPAT_SA1chips}; do
    sample_name="ChIP_AML_${SAMPLE}_SA1"
    TotalPeaks=$(grep -w "# total peaks =" $sample_name/peaks.txt | awk '{print $5}')
    IPefficiency=$(grep -w "# Approximate IP efficiency =" $sample_name/peaks.txt | awk '{print $6}')
    TotalTags=$(grep -w "# Total tags =" $sample_name/peaks.txt | awk '{print $5}')
    echo -e "${sample_name}\t${TotalPeaks}\t${IPefficiency}\t${TotalTags}" >> $PEAKDIR/SA1peaks.default.txt
done
#------------------------------------------------------------------------------------                         

                         #####################################
                         #    normalize SA1 TAGDIRS by CNV   #
                         #####################################

#Usage: normalizeTagDirByCopyNumber.pl <tagDir> -cnv <CNV file>
for patID in ${ALLPAT_SA1chips};do
normalizeTagDirByCopyNumber.pl ${TAGDIR}/ChIP_AML_${patID}_SA1 -cnv ${CNVDIR}/ChIP_IN_AML_${patID}.sam_CNVs
done

#------------------------------------------------------------------------------------                         

                    #############################################
                    #  SA1   peak finding  in CNV norm tagdirs  #
                    #############################################

#loop through patients
for patID in ${ALLPAT_SA1chips};do
findPeaks ${TAGDIR}/ChIP_AML_${patID}_SA1_CNVnorm -i ${INPUTDIR}/ChIP_IN_AML_${patID} -style factor -o auto
done


#save to files
echo -e "Sample_ID\tTotalPeaks\tIPefficiency\tTotalTags" > $PEAKDIR/STAG1peaks.CNVnorm.txt
for SAMPLE in ${ALLPAT_SA1chips}; do
    sample_name="ChIP_AML_${SAMPLE}_SA1_CNVnorm"
    TotalPeaks=$(grep -w "# total peaks =" $sample_name/peaks.txt | awk '{print $5}')
    IPefficiency=$(grep -w "# Approximate IP efficiency =" $sample_name/peaks.txt | awk '{print $6}')
    TotalTags=$(grep -w "# Total tags =" $sample_name/peaks.txt | awk '{print $5}')
    echo -e "${sample_name}\t${TotalPeaks}\t${IPefficiency}\t${TotalTags}" >> $PEAKDIR/STAG1peaks.CNVnorm.txt
done
#------------------------------------------------------------------------------------                         


                     ############################################
                     #SA1: CNV norm BIGwigs and average BIGwigs #
                     ############################################

#generate CNV norm bw file
for SAMPLE in ${ALLPAT_SA1chips}
do
makeUCSCfile ${TAGDIR}/ChIP_AML_${SAMPLE}_SA1*_CNVnorm -o ${BWDIR}/ChIP_AML_${SAMPLE}_SA1_CNVnorm.bigWig -fragLength 200 -bigWig ${CHROMSIZES_HG38} -fsize 1e20
done
#------------------------------------------------------------------------------------                         
#assemble average BW file name by condition
normSA2SA1bw=()
for i in ${PATIENTS_SA2_SA1chips}; do
    bGname=ChIP_AML_${i}_SA1_CNVnorm.bigWig
    normSA2SA1bw+=("$bGname")
done
echo "${normSA2SA1bw[@]}"

normRAD21SA1bw=()
for i in ${PATIENTS_SA2_SA1chips}; do
    bGname=ChIP_AML_${i}_SA1_CNVnorm.bigWig
    normRAD21SA1bw+=("$bGname")
done
echo "${normRAD21SA1bw[@]}"

normCTRLSA1bw=()
for i in ${PATIENTS_CTRL_SA1chips}; do
    bGname=ChIP_AML_${i}_SA1_CNVnorm.bigWig
    normCTRLSA1bw+=("$bGname")
done
echo "${normCTRLSA1bw[@]}"
#------------------------------------------------------------------------------------                         
##run AverageBigwig.pl script and create bedgraphs in addition
cd ${BWDIR}
#CTRL AML
myAverageBigWig.pl -bw ${normCTRLSA1bw[@]} -chr ${CHROMSIZES_HG38} -o ${MERGEDBW}/ave.AML_CTRL_SA1_CNVnorm.bigwig
bigWigToBedGraph ${MERGEDBW}/ave.AML_CTRL_SA1_CNVnorm.bigwig ${MERGEDBW}/ave.AML_CTRL_SA1.CNVnorm.bedGraph

#RAD21mut AML
myAverageBigWig.pl -bw ${normSA2SA1bw[@]} -chr ${CHROMSIZES_HG38} -o ${MERGEDBW}/ave.AML_RAD21mut_SA1_CNVnorm.bigwig
bigWigToBedGraph ${MERGEDBW}/ave.AML_RAD21mut_SA1_CNVnorm.bigwig ${MERGEDBW}/ave.AML_RAD21mut_SA1.CNVnorm.bedGraph

#STAG2mut AML
myAverageBigWig.pl -bw ${normSA2SA1bw[@]} -chr ${CHROMSIZES_HG38} -o ${MERGEDBW}/ave.AML_SA2mut_SA1_CNVnorm.bigwig
bigWigToBedGraph ${MERGEDBW}/ave.AML_SA2mut_SA1_CNVnorm.bigwig ${MERGEDBW}/ave.AML_SA2mut_SA1.CNVnorm.bedGraph
#------------------------------------------------------------------------------------                         
##look at individual patient tracks at an exemplary locus:
mkdir ${FIGURESDIRpyg}/SA1tracks/
##pygenometracks using ini file with all bigwigs listed
pyGenomeTracks_v3.5.sh --tracks ${FIGURESDIRpyg}/inifiles/Pytracks_SA2_RAD21mutvsCTRL.SA1.allPat.ini --region chr3:36550000-38550000 --dpi 600 --fontSize 28 -o ${FIGURESDIRpyg}/SA1tracks/Pytracks_SA2_RAD21mutvsCTRL.SA1.allPat.ITGA9.png
#------------------------------------------------------------------------------------  


                     #####################################################
                     #merged STAG1 Peakset based on individual peakfiles #
                     #####################################################
#------------------------------------------------------------------------------------
##define vector with peak files
SA1peaks_norm=()
for i in ${ALLPAT_SA1chips}; do
    bGname=ChIP_AML_${i}_SA1*CNVnorm/peaks.txt
    SA1peaks_norm+=("$bGname")
done
echo "${SA1peaks_norm[@]}"
##merge and filter peaks
cd ${TAGDIR}
mergePeaks ${SA1peaks_norm[@]} -code > ${PEAKDIR}/Allpat_mergePeaks_SA1.txt
###filtering for blacklisted regions and low mappability
pos2bed.pl ${PEAKDIR}/Allpat_mergePeaks_SA1.txt > ${TMPDIR}/tmp.1.bed
	$BEDTOOLS intersect -a ${TMPDIR}/tmp.1.bed -b $BLACKLIST_HG38 -v > ${TMPDIR}/tmp.2.bed
	bed2pos.pl ${TMPDIR}/tmp.2.bed > ${TMPDIR}/tmp.1.txt
	filter4Mappability.sh -p ${TMPDIR}/tmp.1.txt -g hg38 -f 0.8 -s 75
	pos2bed.pl ${TMPDIR}/tmp.1.mapScoreFiltered.txt > ${PEAKDIR}/Allpat_mergePeaks_SA1.filtered.peaks.bed
	bed2pos.pl ${PEAKDIR}/Allpat_mergePeaks_SA1.filtered.peaks.bed > ${PEAKDIR}/Allpat_mergePeaks_SA1.filtered.peaks.txt
#------------------------------------------------------------------------------------                         






                        #######################################
                        #######################################
                        ##  STAG2 vs STAG1 Chipseq compared  ##
                        #######################################
                        #######################################   

#------------------------------------------------------------------------------------                         
                         #############################################
                         # Average Histogram plots  "global" levels  #
                         #############################################
#annotate to total stringent filtered RAD21peak merged positions of all AML pat SA1 and SA2
mkdir ${FIGURESDIR}/ghist/SA1vsSA2/
#using the stringent filtered RAD21 peakset (based on all patients) and the averaged bedgraphs
annotatePeaks.pl ${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.stringent.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_CTRL_SA2.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/SA1vsSA2/AML_CTRL.str.RAD21pos.CNVnorm.SA2.ghist.txt"
annotatePeaks.pl ${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.stringent.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_SA2mut_SA2.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/SA1vsSA2/AML_SA2mut.str.RAD21pos.CNVnorm.SA2.ghist.txt"
annotatePeaks.pl ${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.stringent.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_RAD21mut_SA2.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/SA1vsSA2/AML_RAD21mut.str.RAD21pos.CNVnorm.SA2.ghist.txt"

annotatePeaks.pl ${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.stringent.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_CTRL_SA1.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/SA1vsSA2/AML_CTRL.str.RAD21pos.CNVnorm.SA1.ghist.txt"
annotatePeaks.pl ${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.stringent.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_SA2mut_SA1.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/SA1vsSA2/AML_SA2mut.str.RAD21pos.CNVnorm.SA1.ghist.txt"
annotatePeaks.pl ${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.stringent.txt hg38 -size 2000 -hist 25 -ghist -bedGraph ${MERGEDBW}/ave.AML_RAD21mut_SA1.CNVnorm.bedGraph > "${FIGURESDIR}/ghist/SA1vsSA2/AML_RAD21mut.str.RAD21pos.CNVnorm.SA1.ghist.txt"


# generate STAG2 mut vs CTRL AML histogram plots
cd ${FIGURESDIR}/ghist/SA1vsSA2

plotHIST.sh -g "AML_CTRL.str.RAD21pos.CNVnorm.SA1.ghist.txt AML_SA2mut.str.RAD21pos.CNVnorm.SA1.ghist.txt" \
-s "AML_CTRL_SA1 AML_SA2mut_SA1" -c "firebrick seagreen3" -x 1000 -y "0 9" -d ${FIGURESDIR}/ghist/SA1vsSA2 -n SA1.str.RAD21pos_CNVnorm.SA2mut.vs.CTRL

plotHIST.sh -g "AML_CTRL.str.RAD21pos.CNVnorm.SA1.ghist.txt AML_RAD21mut.str.RAD21pos.CNVnorm.SA1.ghist.txt" \
-s "AML_CTRL_SA1 AML_RAD21mut_SA1" -c "firebrick plum1" -x 1000 -y "0 9" -d ${FIGURESDIR}/ghist/SA1vsSA2 -n SA1.str.RAD21pos_CNVnorm.RAD21mut.vs.CTRL


plotHIST.sh -g "AML_CTRL.str.RAD21pos.CNVnorm.SA2.ghist.txt AML_SA2mut.str.RAD21pos.CNVnorm.SA2.ghist.txt" \
-s "AML_CTRL_SA2 AML_SA2mut_SA2" -c "firebrick seagreen3" -x 1000 -y "0 9" -d ${FIGURESDIR}/ghist/SA1vsSA2 -n SA2.str.RAD21pos_CNVnorm.SA2mut.vs.CTRL

plotHIST.sh -g "AML_CTRL.str.RAD21pos.CNVnorm.SA1.ghist.txt AML_RAD21mut.str.RAD21pos.CNVnorm.SA1.ghist.txt" \
-s "AML_CTRL_SA2 AML_RAD21mut_SA2" -c "firebrick plum1" -x 1000 -y "0 9" -d ${FIGURESDIR}/ghist/SA1vsSA2 -n SA2.str.RAD21pos_CNVnorm.RAD21mut.vs.CTRL

#------------------------------------------------------------------------------------                         
##look at averge tracks at an exemplary locus (small zoomed in version):
pyGenomeTracks_v3.5.sh --tracks ${FIGURESDIRpyg}/inifiles/Pytracks_SA2mut_CTRL_AML.SA2.average.ini --region chr3:37000000-38000000 --fontSize 40 --dpi 600 -o ${FIGURESDIRpyg}/Pytracks_SA2mutvsCTRL.SA2.average.1v1.png
pyGenomeTracks_v3.5.sh --tracks ${FIGURESDIRpyg}/inifiles/Pytracks_SA2mut_CTRL_AML.SA1.average.ini --region chr3:37000000-38000000 --fontSize 40 --dpi 600 -o ${FIGURESDIRpyg}/Pytracks_SA2mutvsCTRL.SA1.average.1v1.png
#------------------------------------------------------------------------------------                         





                         #############################################################################
                         # Annotate CTRL-AMLS with matched STAG1/STAG2 datasets to all cohesin sites #
                         #############################################################################
##these annotation tables are required for the STAG dominance analyses in cohesin wt AMLs

mkdir ${DIFFPEAKS}/SA1vsSA2

#annotate matched controls to various RAD21 peaksets
#------------------------------------------------------------------------------------                         
##Declare samples
PATIENTS_CTRL_matched='ctr_16911 ctr_18136 ctr_18519 ctr_19405 ctr_19416 ctr_21047 '

CTRLSA1top=()
for i in ${PATIENTS_CTRL_matched}; do
    bGname=ChIP_AML_${i}_SA1-DQI_CNVnorm
    CTRLSA1top+=("$bGname")
done
echo "${CTRLSA1top[@]}"

CTRLSA2top=()
for i in ${PATIENTS_CTRL_matched}; do
    bGname=ChIP_AML_${i}_SA2_CNVnorm
    CTRLSA2top+=("$bGname")
done
echo "${CTRLSA2top[@]}"

#------------------------------------------------------------------------------------                         
#annotate using the AML RAD21 peaksets
###annotate SA1 and SA2 data to stringent RAD21 peak set
cd ${TAGDIR}
annotatePeaks.pl "${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.stringent.txt" hg38 -size 250 -d \
${CTRLSA1top[@]} ${CTRLSA2top[@]} \
-raw -cpu 24 > ${ANNDIR}/CTRL.AML.matched.SA1.SA2.strRAD21.peaks.ann.txt
####remove sex chromosomes and columns not required
grep -v 'chrX' ${ANNDIR}/CTRL.AML.matched.SA1.SA2.strRAD21.peaks.ann.txt > ${TMPDIR}/tmp.filt.txt
grep -v 'chrY' ${TMPDIR}/tmp.filt.txt > ${ANNDIR}/CTRL.AML.matched.SA1.SA2.strRAD21.peaks.ann.filt.txt
tail -n +2 "${ANNDIR}/CTRL.AML.matched.SA1.SA2.strRAD21.peaks.ann.filt.txt" | cut -f1,20-100 > "${ANNDIR}/CTRL.AML.matched.SA1.SA2.strRAD21.peaks.ann.tab.filt.txt"

###annotate SA1 and SA2 data to regular RAD21 peak set
annotatePeaks.pl "${PEAKDIR}/Allpat_mergePeaks_RAD21.filtered.peaks.txt" hg38 -size 250 -d \
${CTRLSA1top[@]} ${CTRLSA2top[@]} \
-raw -cpu 24 > ${ANNDIR}/CTRL.AML.matched.SA1.SA2.RAD21.peaks.ann.txt
####remove sex chromosomes
grep -v 'chrX' ${ANNDIR}/CTRL.AML.matched.SA1.SA2.RAD21.peaks.ann.txt > ${TMPDIR}/tmp.filt.txt
grep -v 'chrY' ${TMPDIR}/tmp.filt.txt > ${ANNDIR}/CTRL.AML.matched.SA1.SA2.RAD21.peaks.ann.filt.txt
tail -n +2 "${ANNDIR}/CTRL.AML.matched.SA1.SA2.RAD21.peaks.ann.filt.txt" | cut -f1,20-100 > "${ANNDIR}/CTRL.AML.matched.SA1.SA2.RAD21.peaks.ann.tab.filt.txt"
#------------------------------------------------------------------------------------ 
#annotate using the HSPC RAD21 peakset
####this analyis will be requred to check for STAG1vsSTAG2 enrichments in the CTRL-AML dataset based on HSPC peak coordinates 
cd ${TAGDIR}
annotatePeaks.pl "${PEAKDIRCD34}/CD34_RAD21.filtered.peaks.txt" hg38 -size 250 -d \
${CTRLSA1top[@]} ${CTRLSA2top[@]} \
-raw -cpu 24 > ${ANNDIR}/CTRL.AML.matched.SA1.SA2.CD34RAD21.peaks.ann.txt
###remove sex chromosomes #might keep them for this analysis though
grep -v 'chrX' ${ANNDIR}/CTRL.AML.matched.SA1.SA2.CD34RAD21.peaks.ann.txt > ${TMPDIR}/tmp.filt.txt
grep -v 'chrY' ${TMPDIR}/tmp.filt.txt > ${ANNDIR}/CTRL.AML.matched.SA1.SA2.CD34RAD21.peaks.ann.filt.txt
#------------------------------------------------------------------------------------   



                         #############################################################################
                         # Annotate all AML RAD21 ChIP tagdirs to the HSPC RAD21 peakset             #
                         #############################################################################
#this analyis will be requred to check for RAD21 enrichments in the complete AML dataset based on HSPC peak coordinates
#------------------------------------------------------------------------------------   
##define vectors with RAD21 tagdirs names by group
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
#------------------------------------------------------------------------------------   
#CD34 peak set --> can be used for peak set enrichments etc at positions defined in HSPCs
cd ${TAGDIR}
annotatePeaks.pl "${PEAKDIRCD34}/CD34_RAD21.filtered.peaks.txt" hg38 -size 250 -d \
${normCTRLRad21[@]} ${normRad21Rad21[@]} ${normSA2Rad21[@]} \
-raw -cpu 24 > ${ANNDIR}/Allpat_CD34.RAD21.peaks.ann.txt
tail -n +2 ${ANNDIR}/Allpat_CD34.RAD21.peaks.ann.txt | cut -f1,20-100 > ${TMPDIR}/tmp.1.txt
cp ${TMPDIR}/tmp.1.txt ${ANNDIR}/Allpat_CD34.RAD21.peaks.ann.Rinput.txt
#remove sex chromosomes
grep -v 'chrX' ${ANNDIR}/Allpat_CD34.RAD21.peaks.ann.Rinput.txt > ${TMPDIR}/tmp.filt.txt
grep -v 'chrY' ${TMPDIR}/tmp.filt.txt > ${ANNDIR}/Allpat_CD34.RAD21.peaks.ann.filt.Rinput.txt
rm ${TMPDIR}/tmp.filt*
#------------------------------------------------------------------------------------   


