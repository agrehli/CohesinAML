#bin/bash
#by Alex F, Jan 2024


################################################################################
#               RNAseq - average coverage tracks for AML groups                #
################################################################################

DIR_SOFT=/misc/software
DIR_DATA=/misc/data
PROJDIR=$DIR_DATA/analysis/project_cohesin/Cohesin_AML
BWDIR=$DIR_DATA/processedData/bigWig/RNA/GRCh38/RNAseq/Cohesin_AML #contains individual bigWigs 
AVBWDIR=$PROJDIR/RNAseq/AverageRNAcoverageTracks
CHROMSIZES_HG38=$DIR_SOFT/viewer/IGV/IGVTools_2.3.98/genomes/GRCh38.PRI_p10.chrom.sizes

mkdir $AVBWDIR


AMLMETA=$PROJDIR/RNAseq/Metadata/RNAseq_Metadata_AML_STAG2_RAD21.txt #contains only IDs used in study


# Define an array containing the RNA sample ids including suffix
ALLsamps=()
while IFS= read -r sn; do
    name=${sn}.Unique.for.bigwig
    ALLsamps+=("$name")
done < <(tail -n +2 ${AMLMETA} | awk -F '\t' '{print $1}')
echo ${ALLsamps[@]}

# Print the array elements to check if everything is there
for sn in "${ALLsamps[@]}"; do
    echo "$sn"
done
#look at complete array
echo ${ALLsamps[@]}

# Get only STAG2mut files
STAG2mutBW=()
for element in "${ALLsamps[@]}"; do # Loop through the elements of ALLTAGDIRS
    # Check if the element contains the string
    if [[ $element == *AML_STAG2* ]]; then
        # If it does, add it to the new array
        STAG2mutBW+=("$element")
    fi
done

echo ${STAG2mutBW[@]}

# Get only RAD21mut files
RAD21mutBW=()
for element in "${ALLsamps[@]}"; do # Loop through the elements of ALLTAGDIRS
    # Check if the element contains the string
    if [[ $element == *AML_RAD21* ]]; then
        # If it does, add it to the new array
        RAD21mutBW+=("$element")
    fi
done

echo ${RAD21mutBW[@]}

# Get only CTRL files
ctrlBW=()
for element in "${ALLsamps[@]}"; do # Loop through the elements of ALLTAGDIRS
    # Check if the element contains the string
    if [[ $element == *AML_CTRL* ]]; then
        # If it does, add it to the new array
        ctrlBW+=("$element")
    fi
done

echo ${ctrlBW[@]}


# Calculate average bigwigs
cd $BWDIR
myAverageBigWig.pl -bw ${ctrlBW[@]} -chr ${CHROMSIZES_HG38} -o ${AVBWDIR}/ave.AML_CTRL_RNA.Unique.for.bigwig
myAverageBigWig.pl -bw ${STAG2mutBW[@]} -chr ${CHROMSIZES_HG38} -o ${AVBWDIR}/ave.AML_STAG2mut_RNA.Unique.for.bigwig
myAverageBigWig.pl -bw ${RAD21mutBW[@]} -chr ${CHROMSIZES_HG38} -o ${AVBWDIR}/ave.AML_RAD21mut_RNA.Unique.for.bigwig





#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

##same for reverse strand  bigwigs

# Define an array containing the RNA sample ids including suffix
ALLsampsrev=()
while IFS= read -r sn; do
    name=${sn}.Unique.rev.bigwig
    ALLsampsrev+=("$name")
done < <(tail -n +2 ${AMLMETA} | awk -F '\t' '{print $1}')
echo ${ALLsampsrev[@]}

# Print the array elements to check if everything is there
for sn in "${ALLsampsrev[@]}"; do
    echo "$sn"
done
#look at complete array
echo ${ALLsampsrev[@]}

# Get only STAG2mut files
STAG2mutBW=()
for element in "${ALLsampsrev[@]}"; do # Loop through the elements of ALLTAGDIRS
    # Check if the element contains the string
    if [[ $element == *AML_STAG2* ]]; then
        # If it does, add it to the new array
        STAG2mutBW+=("$element")
    fi
done

echo ${STAG2mutBW[@]}

# Get only RAD21mut files
RAD21mutBW=()
for element in "${ALLsampsrev[@]}"; do # Loop through the elements of ALLTAGDIRS
    # Check if the element contains the string
    if [[ $element == *AML_RAD21* ]]; then
        # If it does, add it to the new array
        RAD21mutBW+=("$element")
    fi
done

echo ${RAD21mutBW[@]}

# Get only CTRL files
ctrlBW=()
for element in "${ALLsampsrev[@]}"; do # Loop through the elements of ALLTAGDIRS
    # Check if the element contains the string
    if [[ $element == *AML_CTRL* ]]; then
        # If it does, add it to the new array
        ctrlBW+=("$element")
    fi
done

echo ${ctrlBW[@]}


# Calculate average bigwigs for reverse strands (use the -minus option!)
cd $BWDIR
myAverageBigWig.pl -bw ${ctrlBW[@]} -chr ${CHROMSIZES_HG38} -minus -o ${AVBWDIR}/ave.AML_CTRL_RNA.Unique.rev.bigwig
myAverageBigWig.pl -bw ${STAG2mutBW[@]} -chr ${CHROMSIZES_HG38} -minus -o ${AVBWDIR}/ave.AML_STAG2mut_RNA.Unique.rev.bigwig
myAverageBigWig.pl -bw ${RAD21mutBW[@]} -chr ${CHROMSIZES_HG38} -minus -o ${AVBWDIR}/ave.AML_RAD21mut_RNA.Unique.rev.bigwig
