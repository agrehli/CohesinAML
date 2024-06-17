#bin/bash
#by Alex F, Jan 2024


################################################################################
#               RNAseq - average coverage tracks for CD34 HSPC groups          #
################################################################################

DIR_SOFT=/misc/software
DIR_DATA=/misc/data
PROJDIR=$DIR_DATA/analysis/project_cohesin/CD34
BWDIR=$DIR_DATA/processedData/bigWig/RNA/GRCh38/RNAseq/CD34 #contains individual bigWigs 
AVBWDIR=$PROJDIR/RNAseq/AverageRNAcoverageTracks
CHROMSIZES_HG38=$DIR_SOFT/viewer/IGV/IGVTools_2.3.98/genomes/GRCh38.PRI_p10.chrom.sizes
TMPDIR=/loctmp/

mkdir $AVBWDIR


KDMETA=$PROJDIR/RNAseq/Metadata/Metadata_RNAseq_HSPCs_cohesin_KD.txt #contains only IDs used in study


# Define an array containing the RNA sample ids including suffix
ALLsamps=()
while IFS= read -r sn; do
    name=${sn}.Unique.for.bigwig
    ALLsamps+=("$name")
done < <(tail -n +2 ${KDMETA} | awk -F '\t' '{print $1}')
echo ${ALLsamps[@]}

# Print the array elements to check if everything is there
for sn in "${ALLsamps[@]}"; do
    echo "$sn"
done
#look at complete array
echo ${ALLsamps[@]}


#The forward bigwigs were the wrong way around ... need to be inverted 
##NOTE: (wrong option in mapping pipeline used, invert strand option needs to be set at 0 for Illumina TrueSeq mRNA )
cd $BWDIR
for sn in "${ALLsamps[@]}"; do
bigWigToBedGraph $sn ${TMPDIR}/$sn.bg
awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4*(-1)}' ${TMPDIR}/$sn.bg > ${TMPDIR}/$sn.inv.bg
bedGraphToBigWig ${TMPDIR}/$sn.inv.bg ${CHROMSIZES_HG38} $BWDIR/$sn
rm ${TMPDIR}/$sn.bg ${TMPDIR}/$sn.inv.bg
done


# Get only SA2 KD files
STAG2KDBW=()
for element in "${ALLsamps[@]}"; do # Loop through the elements of ALLTAGDIRS
    # Check if the element contains the string SA2 but not SA1 (to avoid the double KD samples)
    if [[ $element == *_SA2_* ]] && [[ $element != *_SA1_* ]]; then
        # If it does, add it to the new array
        STAG2KDBW+=("$element")
    fi
done
echo ${STAG2KDBW[@]}

# Get only SA1 KD files
STAG1KDBW=()
for element in "${ALLsamps[@]}"; do # Loop through the elements of ALLTAGDIRS
    # Check if the element contains the string SA1 but not SA2 (to avoid the double KD samples)
    if [[ $element == *_SA1_* ]] && [[ $element != *_SA2_* ]]; then
        # If it does, add it to the new array
        STAG1KDBW+=("$element")
    fi
done
echo ${STAG1KDBW[@]}

# Get only SA1+SA2 KD files
STAG12KDBW=()
for element in "${ALLsamps[@]}"; do # Loop through the elements of ALLTAGDIRS
    # Check if the element contains the string SA1 but not SA2 (to avoid the double KD samples)
    if [[ $element == *_SA1_* ]] && [[ $element == *_SA2_* ]]; then
        # If it does, add it to the new array
        STAG12KDBW+=("$element")
    fi
done
echo ${STAG12KDBW[@]}

# Get only RAD21 KD files
RAD21KDBW=()
for element in "${ALLsamps[@]}"; do # Loop through the elements of ALLTAGDIRS
    # Check if the element contains the string
    if [[ $element == *RAD21* ]]; then
        # If it does, add it to the new array
        RAD21KDBW+=("$element")
    fi
done
echo ${RAD21KDBW[@]}



# Get only CTRL files (all types of ctrls combined!)
ctrlBW=()
for element in "${ALLsamps[@]}"; do # Loop through the elements of ALLTAGDIRS
    # Check if the element contains the string
    if [[ $element == *siCtrl* ]] || [[ $element == *Mock* ]] || [[ $element == *untreated* ]]; then
        # If it does, add it to the new array
        ctrlBW+=("$element")
    fi
done

echo ${ctrlBW[@]}

for sn in ${ctrlBW[@]};
do
if test -f "${sn}"; then
    echo "${sn} exists."
fi
done




# Calculate average bigwigs
cd $BWDIR
myAverageBigWig.pl -bw ${ctrlBW[@]} -chr ${CHROMSIZES_HG38} -o ${AVBWDIR}/ave.CD34_CTRL_RNA.Unique.for.bigwig
myAverageBigWig.pl -bw ${STAG2KDBW[@]} -chr ${CHROMSIZES_HG38} -o ${AVBWDIR}/ave.CD34_STAG2KD_RNA.Unique.for.bigwig
myAverageBigWig.pl -bw ${STAG1KDBW[@]} -chr ${CHROMSIZES_HG38} -o ${AVBWDIR}/ave.CD34_STAG1KD_RNA.Unique.for.bigwig
myAverageBigWig.pl -bw ${STAG12KDBW[@]} -chr ${CHROMSIZES_HG38} -o ${AVBWDIR}/ave.CD34_STAG1_2KD_RNA.Unique.for.bigwig
myAverageBigWig.pl -bw ${RAD21KDBW[@]} -chr ${CHROMSIZES_HG38} -o ${AVBWDIR}/ave.CD34_RAD21KD_RNA.Unique.for.bigwig







#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

##for reverse bigwigs
# Define an array containing the RNA sample ids including suffix
ALLsamps1=()
while IFS= read -r sn; do
    name=${sn}.Unique.rev.bigwig
    ALLsamps1+=("$name")
done < <(tail -n +2 ${KDMETA} | awk -F '\t' '{print $1}')
echo ${ALLsamps1[@]}
cd $BWDIR
for sn in "${ALLsamps1[@]}"; do
bigWigToBedGraph $sn ${TMPDIR}/$sn.bg
awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4*(-1)}' ${TMPDIR}/$sn.bg > ${TMPDIR}/$sn.inv.bg
bedGraphToBigWig ${TMPDIR}/$sn.inv.bg ${CHROMSIZES_HG38} $BWDIR/$sn
rm ${TMPDIR}/$sn.bg ${TMPDIR}/$sn.inv.bg
done


# Get only SA2 KD files
STAG2KDBW=()
for element in "${ALLsamps1[@]}"; do # Loop through the elements of ALLTAGDIRS
    # Check if the element contains the string SA2 but not SA1 (to avoid the double KD samples)
    if [[ $element == *_SA2_* ]] && [[ $element != *_SA1_* ]]; then
        # If it does, add it to the new array
        STAG2KDBW+=("$element")
    fi
done
echo ${STAG2KDBW[@]}

# Get only SA1 KD files
STAG1KDBW=()
for element in "${ALLsamps1[@]}"; do # Loop through the elements of ALLTAGDIRS
    # Check if the element contains the string SA1 but not SA2 (to avoid the double KD samples)
    if [[ $element == *_SA1_* ]] && [[ $element != *_SA2_* ]]; then
        # If it does, add it to the new array
        STAG1KDBW+=("$element")
    fi
done
echo ${STAG1KDBW[@]}

# Get only SA1+SA2 KD files
STAG12KDBW=()
for element in "${ALLsamps1[@]}"; do # Loop through the elements of ALLTAGDIRS
    # Check if the element contains the string SA1 but not SA2 (to avoid the double KD samples)
    if [[ $element == *_SA1_* ]] && [[ $element == *_SA2_* ]]; then
        # If it does, add it to the new array
        STAG12KDBW+=("$element")
    fi
done
echo ${STAG12KDBW[@]}

# Get only RAD21 KD files
RAD21KDBW=()
for element in "${ALLsamps1[@]}"; do # Loop through the elements of ALLTAGDIRS
    # Check if the element contains the string
    if [[ $element == *RAD21* ]]; then
        # If it does, add it to the new array
        RAD21KDBW+=("$element")
    fi
done
echo ${RAD21KDBW[@]}



# Get only CTRL files (all types of ctrls combined!)
ctrlBW=()
for element in "${ALLsamps1[@]}"; do # Loop through the elements of ALLTAGDIRS
    # Check if the element contains the string
    if [[ $element == *siCtrl* ]] || [[ $element == *Mock* ]] || [[ $element == *untreated* ]]; then
        # If it does, add it to the new array
        ctrlBW+=("$element")
    fi
done

# Calculate average bigwigs
cd $BWDIR
myAverageBigWig.pl -bw ${ctrlBW[@]} -chr ${CHROMSIZES_HG38} -minus -o ${AVBWDIR}/ave.CD34_CTRL_RNA.Unique.rev.bigwig
myAverageBigWig.pl -bw ${STAG2KDBW[@]} -chr ${CHROMSIZES_HG38} -minus -o ${AVBWDIR}/ave.CD34_STAG2KD_RNA.Unique.rev.bigwig
myAverageBigWig.pl -bw ${STAG1KDBW[@]} -chr ${CHROMSIZES_HG38} -minus -o ${AVBWDIR}/ave.CD34_STAG1KD_RNA.Unique.rev.bigwig
myAverageBigWig.pl -bw ${STAG12KDBW[@]} -chr ${CHROMSIZES_HG38} -minus -o ${AVBWDIR}/ave.CD34_STAG1_2KD_RNA.Unique.rev.bigwig
myAverageBigWig.pl -bw ${RAD21KDBW[@]} -chr ${CHROMSIZES_HG38} -minus -o ${AVBWDIR}/ave.CD34_RAD21KD_RNA.Unique.rev.bigwig



