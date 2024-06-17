
#bin/bash
#by Alexander Fischer, Oct 2021
###############################################################################
###############################################################################
##                                                                           ##
##     		 Analysis of ATAC data of Cohesin mut AML patients       		       ##
##                                                                           ##
###############################################################################
###############################################################################

#setting basic paths

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


RAWDIR="${DIR_DATA}/rawData/chromatin/ATAC/Cohesin_AML"
FASTQCDIR="${RAWDIR}/FastQC"
TAGDIR="${DIR_DATA}/processedData/tagDir/chromatin/GRCh38/ATAC/AML/Cohesin_AML"
BWDIR="${DIR_DATA}/processedData/bigWig/chromatin/GRCh38/ATAC/AML/Cohesin_AML"

# general directories
TMPDIR="/loctmp"
## directories with data created by mapChIP.sh pipeline
INPUTDIR="${DIR_DATA}/processedData/tagDir/DNA/GRCh38/Input/AML"
CNVDIR="${DIR_DATA}/processedData/mapping/DNA/GRCh38/Input/AML/CNVdata"

##  directories for analysis
PROJECTDIR="${DIR_DATA}/analysis/project_cohesin"
WORKDIR="${PROJECTDIR}/Cohesin_AML/ATAC"
FIGURESDIR="${WORKDIR}/figures"
PEAKDIR_CNVnorm="${WORKDIR}/peaks_CNVnorm"
DIFFDIR="${WORKDIR}/diffPeaks"
MERGEDBW="${WORKDIR}/mergedBigWigs"

##RAD21 peaks
CHIPPEAKDIRRAD="${PROJECTDIR}/Cohesin_AML/ChIP_analysis/peaks"

# Defining necessary files
CHROMSIZES_HG38="${DIR_SOFT}/viewer/IGV/IGVTools_2.3.98/genomes/GRCh38.PRI_p10.chrom.sizes"
BLACKLIST_HG38="${DIR_SOFT}/analysis/generalStuff/annotation/GRCh38/hg38.blacklist.bed"

#create new dirs
mkdir ${WORKDIR}
mkdir ${PEAKDIR_CNVnorm}
mkdir ${MERGEDBW}
mkdir ${DIFFDIR}


                     #########################
                     #     intial QC         #
                     #########################
#------------------------------------------------------------------------------------                         

#Declare samples
CTRLpat='ATAC_AML_ctr_10156 ATAC_AML_ctr_16911 ATAC_AML_ctr_18136 ATAC_AML_ctr_19416 ATAC_AML_ctr_21047 '
SA2mut='ATAC_AML_SA2_22525 ATAC_AML_SA2_24743 ATAC_AML_SA2_27493 ATAC_AML_SA2_29728 '
RAD21mut='ATAC_AML_RAD21_27650 ATAC_AML_RAD21_AML0117 ATAC_AML_RAD21_UKR186'
ATACTADIRSpats=$CTRLpat$SA2mut$RAD21mut
#------------------------------------------------------------------------------------                         

#look at peak finding QC (generated per default by mapATAC pipeline)
cd ${TAGDIR}
for SAMPLE in ${ATACTADIRSpats};do
paste <(echo -e "# ${SAMPLE}") <(grep -w "# total peaks =" ${SAMPLE}/peaks.txt) <(grep -w "# Approximate IP efficiency"  ${SAMPLE}/peaks.txt) <(grep -w "# Total tags ="  ${SAMPLE}/peaks.txt)
done
#------------------------------------------------------------------------------------                         



                     ###########################################
                     #     CNV normalization of TAGDIRS        #
                     ###########################################
#------------------------------------------------------------------------------------                         
##CNV normalization using DNAinputs
cd ${TAGDIR}
for patID in ${ATACTADIRSpats};
do
normalizeTagDirByCopyNumber.pl ${TAGDIR}/ATAC_AML_${patID} -cnv ${CNVDIR}/ChIP_IN_AML_${patID}.sam_CNVs
done
#------------------------------------------------------------------------------------                         


                     #####################################################
                     #     Peak finding using  findATACpeaks2.0.sh       #
                     #####################################################
# running findATACpeaks2.0.sh for CNV norm tagdirs in different screen sessions in parallel
#------------------------------------------------------------------------------------                         
for SAMPLE in ${ATACTADIRSpats}
do
findATACpeaks2.0.sh -t ${TAGDIR}/${SAMPLE}_CNVnorm -o ${PEAKDIR_CNVnorm} -g hg38 -s
done
#------------------------------------------------------------------------------------                         


                     #####################################################
                     #     Average BigWigs by patient group              #
                     #####################################################
# vectors with bigwig names
#------------------------------------------------------------------------------------                         

bigwigsSA2=()
for i in ${SA2mut}; do
    bGname=${i}.bigwig
    bigwigsSA2+=("$bGname")
done
echo "${bigwigsSA2[@]}"

bigwigsCTRL=()
for i in ${CTRLpat}; do
    bGname=${i}.bigwig
    bigwigsCTRL+=("$bGname")
done
echo "${bigwigsCTRL[@]}"

bigwigsRAD=()
for i in ${RAD21mut}; do
    bGname=${i}.bigwig
    bigwigsRAD+=("$bGname")
done
echo "${bigwigsRAD[@]}"
#------------------------------------------------------------------------------------                         

# generate average bigwigs
cd ${BWDIR}
myAverageBigWig.pl -bw ${bigwigsSA2[@]} -chr ${CHROMSIZES_HG38} -o ${MERGEDBW}/ave.SA2mut.ATAC.bigWig
myAverageBigWig.pl -bw ${bigwigsCTRL[@]} -chr ${CHROMSIZES_HG38} -o ${MERGEDBW}/ave.CTRL.ATAC.bigWig
myAverageBigWig.pl -bw ${bigwigsRAD[@]} -chr ${CHROMSIZES_HG38} -o ${MERGEDBW}/ave.RAD21mut.ATAC.bigWig

#------------------------------------------------------------------------------------                         

bigWigToBedGraph ${MERGEDBW}/ave.SA2mut.ATAC.bigWig ${MERGEDBW}/ave.SA2mut.ATAC.bedGraph
bigWigToBedGraph ${MERGEDBW}/ave.CTRL.ATAC.bigWig ${MERGEDBW}/ave.CTRL.ATAC.bedGraph
bigWigToBedGraph ${MERGEDBW}/ave.RAD21mut.ATAC.bigWig ${MERGEDBW}/ave.RAD21mut.ATAC.bedGraph
#------------------------------------------------------------------------------------                         



                     ###############################################################
                     #     Annotation of ATAC Tagdirs to RAD21 postions            #
                     ###############################################################

#get coverage /FC etc for each RAD21 peak -> this will be useful for correlation analyses as there will be 1value/peak
cd ${TAGDIR}
#------------------------------------------------------------------------------------                         

#using all RAD21 peaks
annotatePeaks.pl "${CHIPPEAKDIRRAD}/Allpat_mergePeaks_RAD21.filtered.peaks.txt" hg38 -size given -len 65 -d ${goodqual} -noann -nogene -noadj -cpu 12 > "${PEAKDIR_CNVnorm}/ATAC_Cohesin_AML.Allpat_mergePeaks_RAD21.filtered.peaks.ann.txt"   ##74425 
tail -n +2 "${PEAKDIR_CNVnorm}/ATAC_Cohesin_AML.Allpat_mergePeaks_RAD21.filtered.peaks.ann.txt" | cut -f1,8-52 > "${PEAKDIR_CNVnorm}/ATAC_Cohesin_AML.Allpat_mergePeaks_RAD21.filtered.peaks.ann.tab.txt" #181528

#using stringent RAD21 peaks
annotatePeaks.pl "${CHIPPEAKDIRRAD}/Allpat_mergePeaks_RAD21_stringent.filtered.peaks.txt" hg38 -size given -len 65 -d ${goodqual} -noann -nogene -noadj -cpu 12 > "${PEAKDIR_CNVnorm}/ATAC_Cohesin_AML.Allpat_mergePeaks_strRAD21.filtered.peaks.ann.txt" 
tail -n +2 "${PEAKDIR_CNVnorm}/ATAC_Cohesin_AML.Allpat_mergePeaks_strRAD21.filtered.peaks.ann.txt" | cut -f1,8-52 > "${PEAKDIR_CNVnorm}/ATAC_Cohesin_AML.Allpat_mergePeaks_strRAD21.filtered.peaks.ann.tab.txt" #128342
#------------------------------------------------------------------------------------                         




