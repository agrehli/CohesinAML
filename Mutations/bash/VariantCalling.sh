#!/bin/bash
set -e

#--------------------------------------------------------------------------------------------------
# Name        : VariantAnalysis.sh
# Author      : Roger Mulet
# Version     :
# Description :
# Usage       : VariantAnalysis.sh -r Reference Genome -t tumor.bam -n normal.bam 
#--------------------------------------------------------------------------------------------------

start=`date +%s`
LOGN=1

[ ! $(module is-loaded strelka) ] && module load strelka
[ ! $(module is-loaded samtools) ] && module load samtools
[ ! $(module is-loaded bcftools) ] && module load bcftools
[ ! $(module is-loaded GATK) ] && module load GATK

JAVA7='/tools/java7/bin/java'
JAVA8='/tools/java8/bin/java'

STRELKA_GERMLINE="configureStrelkaGermlineWorkflow.py"
STRELKA_SOMATIC="configureStrelkaSomaticWorkflow.py"
MUTECT='/tools/MuTect/mutect/target/mutect-1.1.7.jar'
PICARD_TOOLS='/tools/picard-tools/picard-2.17.10.jar'
PINDEL='/tools/pindel/pindel'
PINDEL2VCF='/tools/pindel/pindel2vcf'
UNION_ANNOVAR='/tools/set_operators/union_annovar'
VARSCAN='/tools/varscan/VarScan.v2.4.2.jar'

#--------------------------------------------------------------------------------------------------
#Get input arguments
#--------------------------------------------------------------------------------------------------

display_usage()  {

cat << HEREDOC
Usage: VariantAnalysis.sh -r referenceSequence -t tumor.bam [-n normal.bam -b regions.bed -s SNPDB -w] [-h]

Optional arguments:

  -w|--wgs        WGS mode. Pindel will not be run (too slow) and --pileup-regions will be used in AnnotateBamStatistics [FALSE]
  -b|--bed        BED file containing regions for variant calling, e.g. exome
  -s|--snp	  DBSNP for AnnotateBamStatistics. If not provided, it will be detected from the reference genome
  -n|--normal     BAM file of the normal sample used as control
  -o|--output     Output prefix for the variant files. By default, use the input minus the extension
  -T|--threads    Number of threads for the annotateBamStatistics step
    
HEREDOC
}

if [[ $# == 0 ]]
then
	echo -e 'You have not provided any arguments (see -h for help)\n'
	exit -1
fi

ALL_TOOLS='mutect,mutect2,haplotype_caller,bcftools,strelka,varscan,pindel'
WGS=FALSE
THREADS=4

while [[ $# -gt 0 ]]
do
        case "$1" in
			-r|--reference)
			REFERENCE=$2
			shift
			;;
			-t|--tumor)
			TUMOR=$2
			shift
			;;
			-n|--normal)
			NORMAL=$2
			shift
			;;
    		-w|--wgs)
	        WGS=TRUE
	        ;;
	        -b|--bed)
	        BED=$2
	        shift
	        ;;
            -o|--output)
            PREFIX=$2
            shift
            ;;
            -u|--use)
            USE=$2
            shift
            ;;
			-T|--threads)
			THREADS=$2
            shift
            ;;
            -s|--snp)
            SNPDB=$2
            shift
			;;
			-h|--help)
			display_usage
			exit 0
			;;
			*) # No more options
			echo -e "Error: Invalid option -$2\n"
		        display_usage
			exit 1
			;;
        esac
shift
done  

echo ''

#--------------------------------------------------------------------------------------------------
#Check input arguments
#--------------------------------------------------------------------------------------------------

if [[ -z $REFERENCE ]]; then
	echo 'Error: No reference genome specified.'
    exit -1
   
elif [[ ! -f $REFERENCE ]]; then
	echo -e "Error: Specified reference file $REFERENCE is not a valid file or does not exist\n"
	exit -1

fi

if [[ -z $TUMOR ]]; then
	echo 'Error: Please specify a tumor sample file (-h for help)'
	exit -1
fi

if [[ ! -f $TUMOR ]]; then
	echo "Error: Specified tumor sample file $TUMOR is not a valid file or does not exist"
	exit -1
fi

if [[ -z $PREFIX ]]; then

    PREFIX=`echo $TUMOR | sed 's:\.bam$::'` # Base name of the bam file. Also basename $TUMOR

fi

if [[ -z $NORMAL ]]; then
	echo -e 'Warning: Please note that a normal sample is missing. The analysis will continue without it\n'	

elif [[ ! -f $NORMAL ]]; then
	echo "Error: Specified normal sample file $NORMAL is not a valid file or does not exist\n"
	exit -1

fi

if [[ ! -z $BED && ! -f $BED ]]; then
	echo 'Error: The BED file specified does not exist (-h for help)'
	exit -1
fi

if [[ -z $USE ]]; then

    echo 'Warning: No USE parameter specified. All tools will be used: $USE'
    USE=$ALL_TOOLS

fi    

# Redirect stdout ( > ) into a named pipe ( >() ) running "tee".
rm -f ${PREFIX}_VariantCalling.log
exec > >(tee -a ${PREFIX}_VariantCalling.log) 2>&1

TOOLS_USED=${USE//,/ }

CHECK_TOOLS=$(comm -23 <(sort <(echo -e ${USE//,/\\n})) <(sort <(echo -e ${ALL_TOOLS//,/\\n})) )

if [[ ! -z $CHECK_TOOLS ]]; then

    echo -e "ERROR: Tool $CHECK_TOOLS is not available for analysis! Chose only from $ALL_TOOLS" && exit -1

fi

#-----------------------------------------------------------------
#Tools, folders and helper functions
#-----------------------------------------------------------------

# Determine what annotation script will be used based on the reference
if [[ $REFERENCE =~ hg19 || $REFERENCE =~ GRCh37 || $REFERENCE =~ v37 ]]; then
	BUILD=hg19
    if [[ -z $SNPDB ]]; then SNPDB=/genomes/homo_sapiens/hg19/snp_db/common_snp151_20180423_chr_clean.vcf.gz; fi
elif [[ $REFERENCE =~ hg38 || $REFERENCE =~ GRCh38 ]]; then
	BUILD=hg38
    if [[ -z $SNPDB ]]; then SNPDB=/genomes/homo_sapiens/hg38/snp_db/common_snp151_20180418_chr_clean.vcf.gz; fi
fi

if [[ $WGS == TRUE ]]; then # --pileup-regions enabled for annotatebamstatistics

    ANNOTATE="$ANNOTATE -w"

fi

# Display options:
echo -e "> Alignment files: $TUMOR, $NORMAL"
echo -e "> Reference genome: $REFERENCE"
echo -e "> Selected regions: $BED"
echo -e "> Annotation for build: $BUILD"
echo -e "> Tools to be used: $USE"
sleep 5

#--------------------------------------------------------------------------------------------------
#Declare annotation functions
#--------------------------------------------------------------------------------------------------

function convert_to_annovar() { 

    if [[ ! -f $2 ]]; then

	    local OUTPUT_PREFIX=`echo $1 | sed  -r 's:\.vcf(.gz)?$::'`

        if [[ $(bcftools view -h $1 | grep bcftools_norm | wc -l) == 0 ]]; then

            echo " - Info: left-normalize indels and split multiallelic variants $1"
            # http://annovar.openbioinformatics.org/en/latest/articles/VCF

            bcftools norm -Ov -m-both -f $REFERENCE $1 > ${OUTPUT_PREFIX}_norm.vcf

	        if [[ $? -ne 0 ]]
	        then

		        echo "	Error: left-align and normalize indels failed: $1"
                rm ${OUTPUT_PREFIX}_norm.vcf
                exit -1

	        fi

            mv ${OUTPUT_PREFIX}_norm.vcf ${OUTPUT_PREFIX}.vcf
            bgzip -f ${OUTPUT_PREFIX}.vcf
            bcftools index ${OUTPUT_PREFIX}.vcf.gz & echo ""

        fi

	    echo -e " - Info: Convert to annovar variant results $1\n"

	    convert2annovar.pl --format vcf4old $1 | cut -f1-5 | awk 'length($4)< 5000 && length($5) < 5000' > $2
	    if [[ $? -ne 0 ]]
	    then

		    echo "Error: Convert to annovar variant results failed $1"
		    exit -1

	    fi

    else

        echo -e " - Info: $2 already found. Skipping..."

    fi

}

function combine_and_annotate()
{
    local ANNOVAR_FILES=$1
    local ANNOVAR_NAMES=$2
    local BAM=$(echo $3 | sed "s/,$//") # remove trailing commas
    local OUTPUT_PREFIX=$4

    echo $ANNOVAR_FILES ${ANNOVAR_NAMES}

    ## Combine ANNOVAR files AND remove non primary assembly variants
	echo " - Info: Combine annovar files ${OUTPUT_PREFIX}"
	cat $(echo $ANNOVAR_FILES | tr ',' '\n') | sort | uniq | awk '$1 ~ /chr[0-9XYM]+$/ || $1 ~ /^[0-9XYMT]+$/' > variant_calling/${OUTPUT_PREFIX}.annovar

    ## Annotate combined files
	echo " - Info: Annotate combined annovar file ${OUTPUT_PREFIX}"

    if [[ $BUILD == 'hg19' ]]; then

	    table_annovar.pl variant_calling/${OUTPUT_PREFIX}.annovar /genomes/homo_sapiens/hg19/annovar/ --remove --buildver hg19 --protocol refGene,ensGene,band,phastConsElements46way,tfbs,wgEncodeRegTfbsClustered,segdup,wgEncodeBroadHmmK562HMM,avsnp150,ucsc_var151,\
popfreq_all,gnomad_exome,cosmic94,cosmic94_FLAGGED,clinvar_20220320,ljb26_all --operation g,g,r,r,r,r,r,r,f,f,f,f,f,f,f,f --nastring " " # dbnsfp35a, gerp++gt2, esp6500siv2_all

    elif [[ $BUILD == 'hg38' ]]; then

        table_annovar.pl variant_calling/${OUTPUT_PREFIX}.annovar /genomes/homo_sapiens/hg38/annovar/ --remove --buildver hg38 --protocol refGene,ensGene,band,phastConsElements100way,gwascatalog,\
avsnp150,ucsc_var151,gnomad_exome,cosmic98,clinvar_20221231,ljb26_all --operation g,g,r,r,r,f,f,f,f,f,f --nastring " " # ljb26_all, dbnsfp35a -- no way to select specific columns, we have to cut them afterwards, wgEncodeBroadHistoneNhekH3k27acStdSig -- empty

    else

        echo -e "This build is not supported. You will have to modify this script"
        exit -1

    fi

	if [[ $? -ne 0 ]]
	then

		echo "Error: Annotate combined annovar file failed ${OUTPUT_PREFIX}"
		return 1

	fi

	rm variant_calling/${OUTPUT_PREFIX}.annovar
    ## Remove very long indels and move to .txt 
    awk 'length($4)< 20000 && length($5) < 20000' variant_calling/${OUTPUT_PREFIX}.annovar.${BUILD}_multianno.txt > variant_calling/${OUTPUT_PREFIX}.txt
    rm variant_calling/${OUTPUT_PREFIX}.annovar.${BUILD}_multianno.txt

    ## Fix tier annotation
	sed -i 's:bed:tier1:;s:bed2:tier2:;s:bed3:tier3:;s:bed4:tier4:;s:Name=NA:true:g' variant_calling/${OUTPUT_PREFIX}.txt
	
	## Annotate bam statistics
	echo " - Info: Annotate bam statistics ${OUTPUT_PREFIX} with $BAM"
 
      /tools/annotateBamStatistics/annotatebamstatistics_roger -a variant_calling/${OUTPUT_PREFIX}.txt -b ${BAM} -f $REFERENCE -v $SNPDB --max-mismatches 1 --min-diff-suboptimal 5 --min-alignment-score 36 --min-base-score 30 --min-mismatch-base-score 10 -t $THREADS -12345678 > variant_calling/${OUTPUT_PREFIX}.temp # --min-mismatch-base-score 10 -> detect more mismatches   

	if [[ $(tail -1 ${OUTPUT_PREFIX}_VariantCalling.log) != "Info: Done" && $? -ne 0 ]]
	then

		echo "Error: Annotate bam statistics failed ${OUTPUT_PREFIX}"
		return 1

	fi

    mv variant_calling/${OUTPUT_PREFIX}.temp variant_calling/${OUTPUT_PREFIX}.txt

    ## Indicate presence of every variant in input files	
	echo -e " - Info: Annotate individual annovar files ${OUTPUT_PREFIX}"

    /tools/assign_annovar_files/assign_annovar_files --table-file variant_calling/${OUTPUT_PREFIX}.txt --annovar-files $ANNOVAR_FILES --annovar-names $ANNOVAR_NAMES > variant_calling/${OUTPUT_PREFIX}.temp

	if [[ $? -ne 0 ]]
	then

		echo "Error: Annotate individual annovar files failed ${OUTPUT_PREFIX}"
		return 1

	fi
	
    mv variant_calling/${OUTPUT_PREFIX}.temp variant_calling/${OUTPUT_PREFIX}.txt

	return 0
}

#==================================================================================================
# Single point mutations
#==================================================================================================

#--------------------------------------------------------------------------------------------------
#run Mutect
#--------------------------------------------------------------------------------------------------

function run_mutect() {

    if [[ ! -e variant_calling/mutect/${PREFIX}_mutect.vcf && ! -e variant_calling/mutect/${PREFIX}_mutect.vcf.gz ]]; then

        echo -e "\n- Information $LOGN: variant calling with MuTecT" && LOGN=$((LOGN+1))

        if [[ -z $NORMAL ]]; then

            if [[ -z $BED ]]; then

                time $JAVA7 -Xmx6G -jar $MUTECT -U ALLOW_N_CIGAR_READS -T MuTect -R $REFERENCE --gap_events_threshold 100 --strand_artifact_power_threshold 0.9 --pir_mad_threshold 5 --pir_median_threshold 1 --normal_lod 0.6 --max_alt_allele_in_normal_fraction 0.10 --max_alt_alleles_in_normal_count 6 -rf BadCigar --input_file:tumor $TUMOR -o ${PREFIX}_mutect.out -vcf ${PREFIX}_mutect.vcf        

            else

                time $JAVA7 -Xmx6G -jar $MUTECT -U ALLOW_N_CIGAR_READS -T MuTect -R $REFERENCE --gap_events_threshold 100 --strand_artifact_power_threshold 0.9 --pir_mad_threshold 5 --pir_median_threshold 1 --normal_lod 0.6 --max_alt_allele_in_normal_fraction 0.10 --max_alt_alleles_in_normal_count 6 -rf BadCigar --input_file:tumor $TUMOR -L $BED -o ${PREFIX}_mutect.out -vcf ${PREFIX}_mutect.vcf

            fi

        else

            if [[ -z $BED ]]; then

                time $JAVA7 -Xmx6G -jar $MUTECT -U ALLOW_N_CIGAR_READS -T MuTect -R $REFERENCE --gap_events_threshold 100 --strand_artifact_power_threshold 0.9 --pir_mad_threshold 5 --pir_median_threshold 1 --normal_lod 0.6 --max_alt_allele_in_normal_fraction 0.10 --max_alt_alleles_in_normal_count 6 -rf BadCigar --input_file:tumor $TUMOR --input_file:normal $NORMAL -o ${PREFIX}_mutect.out -vcf ${PREFIX}_mutect.vcf

            else

                time $JAVA7 -Xmx6G -jar $MUTECT -U ALLOW_N_CIGAR_READS -T MuTect -R $REFERENCE --gap_events_threshold 100 --strand_artifact_power_threshold 0.9 --pir_mad_threshold 5 --pir_median_threshold 1 --normal_lod 0.6 --max_alt_allele_in_normal_fraction 0.10 --max_alt_alleles_in_normal_count 6 -rf BadCigar --input_file:tumor $TUMOR --input_file:normal $NORMAL -L $BED -o ${PREFIX}_mutect.out -vcf ${PREFIX}_mutect.vcf

            fi            

        fi


        if [[ $? -ne 0 ]]; then
            echo 'Error: Mutect failed!'
            exit -1

        fi

        awk '$0~/^#/ || $7=="PASS"' ${PREFIX}_mutect.vcf > ${PREFIX}_mutect_pass.vcf
        mv ${PREFIX}_mutect_pass.vcf ${PREFIX}_mutect.vcf
        rm ${PREFIX}_mutect.out

        mv ${PREFIX}_mutect* variant_calling/mutect/

        bgzip variant_calling/mutect/${PREFIX}_mutect.vcf
        bcftools index variant_calling/mutect/${PREFIX}_mutect.vcf.gz

    else 

        echo -e "\n- Information 1: MuTecT file already found. Skipping..."

    fi

}

#--------------------------------------------------------------------------------------------------
#run Mutect2
#--------------------------------------------------------------------------------------------------

function run_mutect2() {

    TUMOR_SM=`samtools view -H $TUMOR | grep '^@RG' | head -1 | grep -o 'SM:[^ 	]*' | sed 's;SM:;;'`
    if [[ ! -z $NORMAL ]]; then
        NORMAL_SM=`samtools view -H $NORMAL | grep '^@RG' | head -1 | grep -o 'SM:[^ 	]*' | sed 's;SM:;;'`
    fi

    # WARNING: For our somatic analysis that uses alt-aware and post-alt processed alignments to GRCh38, we disable a specific read filter with --disable-read-filter
    # MateOnSameContigOrNoMappedMateReadFilter. This filter removes from analysis paired reads whose mate maps to a different contig

    if [[ ! -e variant_calling/mutect2/${PREFIX}_mutect2.vcf && ! -e variant_calling/mutect2/${PREFIX}_mutect2.vcf.gz ]]; then

        echo -e "\n- Information $LOGN: variant calling with MuTect2\n" && LOGN=$((LOGN+1))

        if [[ -z $NORMAL ]]; then

            if [[ -z $BED ]]; then

                gatk Mutect2 -R $REFERENCE -I $TUMOR -tumor $TUMOR_SM -O ${PREFIX}_mutect2.vcf --read-filter GoodCigarReadFilter

            else

                gatk Mutect2 -R $REFERENCE -I $TUMOR -tumor $TUMOR_SM -O ${PREFIX}_mutect2.vcf --read-filter GoodCigarReadFilter -L $BED 
        
            fi

        else

            if [[ -z $BED ]]; then

                gatk Mutect2 -R $REFERENCE -I $TUMOR -I $NORMAL -tumor $TUMOR_SM -normal $NORMAL_SM -O ${PREFIX}_mutect2.vcf --read-filter GoodCigarReadFilter

            else

                gatk Mutect2 -R $REFERENCE -I $TUMOR -I $NORMAL -tumor $TUMOR_SM -normal $NORMAL_SM -O ${PREFIX}_mutect2.vcf --read-filter GoodCigarReadFilter -L $BED            

            fi

        fi

        if [[ $? -ne 0 ]]; then
            echo "Error: Mutect2 failed - $PREFIX"
            exit -1
        fi

        #awk '$0~/^#/ || $7=="PASS"' ${PREFIX}_mutect2.vcf > ${PREFIX}_mutect2_pass.vcf
        #mv ${PREFIX}_mutect2_pass.vcf ${PREFIX}_mutect2.vcf

        mv ${PREFIX}_mutect2* variant_calling/mutect2/

        bgzip variant_calling/mutect2/${PREFIX}_mutect2.vcf
        bcftools index variant_calling/mutect2/${PREFIX}_mutect2.vcf.gz

    else

        echo -e "\n- Information 2: MuTecT2 file already found. Skipping..."

    fi

}

#--------------------------------------------------------------------------------------------------
#run HaplotypeCaller
#--------------------------------------------------------------------------------------------------

function run_haplotype_caller() {

    if [[ ! -e variant_calling/haplotype_caller/${PREFIX}_hcaller.vcf && ! -e variant_calling/haplotype_caller/${PREFIX}_hcaller.vcf.gz ]]; then

        echo -e "\n- Information $LOGN: variant calling with Haplotype Caller\n" && LOGN=$((LOGN=1))

        if [[ -z $BED ]]; then

            gatk HaplotypeCaller -R $REFERENCE -I $TUMOR -O variant_calling/haplotype_caller/${PREFIX}_hcaller.vcf

        else

            gatk HaplotypeCaller -R $REFERENCE -I $TUMOR -O variant_calling/haplotype_caller/${PREFIX}_hcaller.vcf -L $BED
    
        fi            

        if [[ $? -ne 0 ]]
        then

            echo "Error: Haplotype caller failed $NAME"
            return 1

        fi

        bgzip variant_calling/haplotype_caller/${PREFIX}_hcaller.vcf
        bcftools index variant_calling/haplotype_caller/${PREFIX}_hcaller.vcf.gz

    else

        echo -e "\n- Information $LOGN: Haplotype caller file already found. Skipping..." && LOGN=$((LOGN=1))

    fi

}

#--------------------------------------------------------------------------------------------------
#run strelka
#--------------------------------------------------------------------------------------------------

function run_strelka() {

    if [[ ! -e variant_calling/strelka/${PREFIX}_strelka.vcf && ! -e variant_calling/strelka/${PREFIX}_strelka.vcf.gz ]]; then

        echo -e "\n- Information $LOGN: variant calling with Strelka\n" && LOGN=$((LOGN+1))
        module load python/2.7.18

        # Check that contig names match FASTA (necessary for strelka)
        # for i in $(grep ">" $REFERENCE | sed "s/^>//" | awk '{print $1}'); do
        #     CHR=$(samtools view -H $TUMOR | grep -P "@SQ.*$i\t")
        #     echo $CHR
        #     if [[ -z $CHR ]]; then
        #         echo -e "Contig $i from reference genome $REFERENCE not found in $TUMOR. Contigs must be the same for Stelka\n" && exit -1
        #     fi
        # done
        
        ## Configure strelka ##	

        # --exome / --targeted: turns off high-depth filters, designed to exclude pericentromeric reference in WGS. Recommended when --cal-Regions is provided.
        
        echo "Info: Configure strelka $NAME"

        if [[ ! -z $BED ]]; then

            if [[ ! -f ${BED}.gz.tbi ]]; then
                bgzip -f -c $BED > ${BED}.gz ; tabix -p bed ${BED}.gz
            fi

            if [[ -z $NORMAL ]]; then

                $STRELKA_GERMLINE --bam=$TUMOR --referenceFasta $REFERENCE --targeted --callRegions=${BED}.gz --runDir=variant_calling/strelka/${PREFIX}

            else

                $STRELKA_SOMATIC --tumorBam $TUMOR --normalBam $NORMAL --referenceFasta $REFERENCE --targeted --callRegions=${BED}.gz --runDir=variant_calling/strelka/${PREFIX}           
            
            fi

        else

            if [[ -z $NORMAL ]]; then

                $STRELKA_GERMLINE --bam=$TUMOR --referenceFasta $REFERENCE --runDir=variant_calling/strelka/${PREFIX}


            else

                $STRELKA_SOMATIC --tumorBam $TUMOR --normalBam $NORMAL --referenceFasta $REFERENCE --runDir=variant_calling/strelka/${PREFIX}


            fi

        fi

        if [[ $? -ne 0 ]]
        then 

            echo "Error: Configure strelka failed $NAME"
            exit -1

        fi

        ## Run Strelka ##
        
        echo "Info: Run strelka $NAME"

        variant_calling/strelka/${PREFIX}/runWorkflow.py --mode=local --jobs=$THREADS --memGb=8

        if [[ $? -ne 0 ]]
        then 

            echo "Error: Run strelka failed $NAME"
            exit -1

        fi

        if [[ -z $NORMAL ]]; then

            zcat variant_calling/strelka/${PREFIX}/results/variants/variants.vcf.gz | awk -F'\t' '$1~/^#/ || $7=="PASS"' > variant_calling/strelka/${PREFIX}_strelka.vcf

        else

            bcftools concat -a variant_calling/strelka/${PREFIX}/results/variants/somatic*vcf.gz | awk -F'\t' '$1~/^#/ || $7=="PASS"' > variant_calling/strelka/${PREFIX}_strelka.vcf 
            # output: somatic.snvs.vcf.gz , somatic.indels.vcf.gz
    
        fi

        #rm -r variant_calling/strelka/${PREFIX}
        bgzip variant_calling/strelka/${PREFIX}_strelka.vcf 
        bcftools index variant_calling/strelka/${PREFIX}_strelka.vcf.gz
        
    else

        echo -e "\n- Information $LOGN: Strelka file already found. Skipping" && LOGN=$((LOGN+1))

    fi

}

#--------------------------------------------------------------------------------------------------
#run bcftools
#--------------------------------------------------------------------------------------------------

function run_bcftools() {

    # -q: minimum MAPQ, -Q minimum base quality, -A: do not skip anomalous read pairs, -d max depth, -L max indel depth, -m multiallelic caller
    # NOTE: -p only affects the probability in classical variant caller

    if [[ ! -e variant_calling/bcftools/${PREFIX}_bcftools.vcf && ! -e variant_calling/bcftools/${PREFIX}_bcftools.vcf.gz ]]; then

        echo -e "\n- Information $LOGN: variant calling with Bcftools\n" && LOGN=$((LOGN+1))

        if [[ -z $BED ]]; then

            bcftools mpileup -q 5 -Q 20 -A -d 100000 -L 100000 -a AD,DP,SP -f $REFERENCE $TUMOR | bcftools call -mv -Ov - | awk -F '\t' -v OFS='\t' '$4 !~ /[RYSWKMBDHV]/ && $5 !~ /[RYSWKMBDHV]/ ||  $1 ~ /^#/{print $0}' | bcftools filter -Oz -i 'INFO/DP > 2' 1> ${PREFIX}_bcftools.vcf.gz

        else

            bcftools mpileup -q 5 -Q 20 -A -d 100000 -L 100000 -a AD,DP,SP -f $REFERENCE -R $BED $TUMOR | bcftools call -mv -Ov - | awk -F '\t' -v OFS='\t' '$4 !~ /[RYSWKMBDHV]/ && $5 !~ /[RYSWKMBDHV]/ ||  $1 ~ /^#/{print $0}' | bcftools filter -Oz -i 'INFO/DP > 2' 1> ${PREFIX}_bcftools.vcf.gz

        fi

        if [[ $? -ne 0 ]]
        then

            echo 'Error: bcftools failed'
            exit -1

        else

            mv ${PREFIX}_bcftools.vcf.gz variant_calling/bcftools/

        fi

        bcftools index variant_calling/bcftools/${PREFIX}_bcftools.vcf.gz

    else

        echo -e "\n- Information 5: Bcftools file already found. Skipping"

    fi

}

#--------------------------------------------------------------------------------------------------
#run varscan2
#--------------------------------------------------------------------------------------------------

function run_varscan() {

    if [[ ! -e variant_calling/varscan/${PREFIX}_varscan.vcf && ! -e variant_calling/varscan/${PREFIX}_varscan.vcf.gz ]]; then

        echo -e "\n- Information $LOGN: variant calling with Varscan\n" && LOGN=$((LOGN+1))

        # NOTE: VCF does not accept ambiguous nucleotides, only N, but Varscan does not convert those to N. Therefore, we add an additional step with awk to fix this

        # Variant AND copy number calling (for separate use)
        # -B = disables SAMtools BAQ computation, which is turned on by default in SAMtools mpileup, but occasionally too stringent for variant calling.     
        # For a single sample mpileup and pileup operate similarly and the output is identical. When you give it multiple BAM files, however, SAMtools mpileup 
        # generates a multi-sample pileup format that must be processed with the mpileup2* commands in VarScan.  

        # time $SAMTOOLS mpileup -B -d 100000 -L 100000 -q 1 -f $REFERENCE $TUMOR  | tee >(java -jar $VARSCAN somatic --mpileup ${PREFIX}_varscan --output-vcf) | awk -F"\t" '$4 > 0 && $7 > 0' | java -jar $VARSCAN copynumber varScan --mpileup 1 ${PREFIX}_varscan

        if [[ -z $NORMAL ]]; then

            if [[ -z $BED ]]; then
                time samtools mpileup -B -d 100000 -L 100000 -q 1 -f $REFERENCE $TUMOR | java -jar $VARSCAN mpileup2cns --output-vcf 1 --variants | awk -F '\t' -v OFS='\t' '$4 !~ /[RYSWKMBDHV]/ && $5 !~ /[RYSWKMBDHV]/ ||  $1 ~ /^#/{print $0}' > variant_calling/varscan/${PREFIX}_varscan.vcf

            else

                time samtools mpileup -B -d 100000 -L 100000 -q 1 -l $BED -f $REFERENCE $TUMOR | java -jar $VARSCAN mpileup2cns --output-vcf 1 --variants | awk -F '\t' -v OFS='\t' '$4 !~ /[RYSWKMBDHVN]/ && $5 !~ /[RYSWKMBDHVN]/ || $1 ~ /^#/{print $0}' > variant_calling/varscan/${PREFIX}_varscan.vcf

            fi

        elif [[ ! -z $NORMAL ]]; then

            if [[ -z $BED ]]; then

                time samtools mpileup -B -d 100000 -L 100000 -q 1 -f $REFERENCE $NORMAL $TUMOR | java -jar $VARSCAN somatic --mpileup variant_calling/varscan/${PREFIX}_varscan --output-vcf

            else

                time samtools mpileup -B -d 100000 -L 100000 -q 1 -l $BED -f $REFERENCE $NORMAL $TUMOR | java -jar $VARSCAN somatic --mpileup variant_calling/varscan/${PREFIX}_varscan --output-vcf 

            fi

        fi

        if [[ $? -ne 0 ]]
        then

            echo 'Error: Varscan failed!'
            exit -1

        fi

        if [[ -z $NORMAL ]]; then

            bgzip variant_calling/varscan/${PREFIX}_varscan.vcf # Compress

        elif [[ ! -z $NORMAL ]]; then    
        # Join lists of SNPs and INDELS
            ls variant_calling/varscan/${PREFIX}_varscan*vcf | xargs -P 10 -I {} bgzip {}
            ls variant_calling/varscan/${PREFIX}_varscan*vcf.gz | xargs -P 10 -I {} tabix -p vcf {}
            
            # We remove variants were ambiguous nucleotides are called with RYSWKMBDHVN
            bcftools concat -a variant_calling/varscan/${PREFIX}*vcf.gz | awk -F '\t' -v OFS='\t' '$4 !~ /[RYSWKMBDHVN]/ && $5 !~ /[RYSWKMBDHVN]/ || $1 ~ /^#/{print $0}' | bgzip > variant_calling/varscan/${PREFIX}_varscan.vcf.gz # bcftools concat -a to allow overlaps (coordinates can overlap each other)    

        fi   

        bcftools index variant_calling/varscan/${PREFIX}_varscan.vcf.gz # index for bcftools norm
        #rm variant_calling/varscan/${PREFIX}_varscan*vcf.gz*

    else

        echo -e "\n- Information $LOGN: Varscan file already found. Skipping" && LOGN=$((LOGN+1))

    fi

}

#==================================================================================================
# SMALL INDELS
#==================================================================================================

#--------------------------------------------------------------------------------------------------
#run pindel (large indels)
#--------------------------------------------------------------------------------------------------

function run_pindel() {

    if [[ ! -e variant_calling/pindel/${PREFIX}_pindel.vcf && ! -e variant_calling/pindel/${PREFIX}_pindel.vcf.gz ]]; then

        echo -e "\n- Information $LOGN: variant calling with Pindel\n" && LOGN=$((LOGN+1))

        READ_SIZE=$(samtools view $TUMOR | head -100000 | awk '{if(length($10) > L)L=length($10)}END{print(L)}')

        ## PINDEL CONFIGURATION: INSERT SIZE ##
        INSERT_METRICS=$(find .. -wholename "*/$PREFIX*insert_size_metrics")

        if [[ -z $INSERT_METRICS ]]; then

            gatk CollectInsertSizeMetrics --INPUT $TUMOR --OUTPUT $PREFIX.insert_size_metrics --Histogram_FILE /dev/null --VALIDATION_STRINGENCY LENIENT

            if [[ $? -ne 0 ]]; then

                echo -e "CollectInsertSizeMetrics failed!" && exit -1

            fi

            INSERT_METRICS=$(find . -name "*$PREFIX*insert_size_metrics")

        fi

        INSERT_SIZE=$(awk 'NR==8{print int($6+0.5)}' $INSERT_METRICS) # Retrieves mean insert size and rounds it up

        if [[ -z $INSERT_METRICS || $INSERT_SIZE < $READ_SIZE ]]; then

            echo 'Warning: Could not generate insert size metrics. Insert size will be retrieved from BAM file'
            INSERT_SIZE=$(samtools view -F 4 $TUMOR | head -100000 | awk -v RSIZE=$READ_SIZE 'length($9)<6 && $9 > RSIZE{sum = sum + sqrt($9*$9); n = n+1}END{print(sum/n)}') # sqrt to have only positive values
        fi	

        echo "${TUMOR} $INSERT_SIZE ${PREFIX}"	> ${PREFIX}_pindel_config.txt # The conf file requires an insert size

        ## RUN PINDEL ##
        if [[ ! -f ${TUMOR%.bam}.bam.bai ]]; then mv ${TUMOR%.bam}.bai ${TUMOR%.bam}.bam.bai; fi

        if [[ -z $BED ]]; then

            time $PINDEL -i ${PREFIX}_pindel_config.txt -f $REFERENCE -o ${PREFIX}_pindel -T $THREADS -L ${PREFIX}_pindel.log ; RET=$? # -T = threads

        else

            time $PINDEL -i ${PREFIX}_pindel_config.txt -f $REFERENCE -o ${PREFIX}_pindel -j $BED -T $THREADS -L ${PREFIX}_pindel.log ; RET=$? # -T = threads, -j = include region

        fi    
            
        if [[ ${RET} -ne 0 ]]
        then

            echo "Error: Pindel failed $TUMOR!"
            exit -1

        fi

        # Move big structural variants to the corresponding folder
        mkdir -p structural_variants/pindel
        mv -t structural_variants/pindel ${PREFIX}_pindel_LI* ${PREFIX}_pindel_INV* ${PREFIX}_pindel_INT* ${PREFIX}_pindel_CloseEndMapped* ${PREFIX}_pindel_RP* ${PREFIX}_pindel_BP*

        # Generate VCF files for small structural variants and compress them
        ls ${PREFIX}_pindel* | grep -vP "log|config" | grep -v vcf | xargs -P 10 -I {} $PINDEL2VCF -p {} -r $REFERENCE -R RefGene -d 20140901
        ls ${PREFIX}_pindel*vcf | xargs -P 10 -I {} bgzip {}
        ls ${PREFIX}_pindel*vcf.gz | xargs -P 10 -I {} tabix -f -p vcf {}

        # Concatenate all pindel files in one and annotate (use wildcards in case some files do not exist)
        # Filter out variants
        bcftools concat -a ${PREFIX}_pindel*vcf.gz | bcftools view -i 'FMT/AD[0:1] > 2' | bcftools sort -Oz -o ${PREFIX}_pindel.vcf.gz
        bcftools index ${PREFIX}_pindel.vcf.gz

        mv ${PREFIX}_pindel* variant_calling/pindel 

    else

        echo -e "\n- Information $LOGN: Pindel file already found. Skipping" && LOGN=$((LOGN+1))

    fi

}

#--------------------------------------------------------------------------------------------------
# Create output folders
#--------------------------------------------------------------------------------------------------

function create_folder() {

    mkdir -p $1

    if [[ $? -ne 0 ]]; then
	    echo 'Error: Could not create output folder'
	exit 1
    fi

}

create_folder variant_calling

for tool in $TOOLS_USED; do

    create_folder variant_calling/$tool

done

#--------------------------------------------------------------------------------------------------
# Run all the tools
#--------------------------------------------------------------------------------------------------

for tool in ${TOOLS_USED}; do

    run_${tool}

done

#==================================================================================================
# Annotate with Annovar and BAM statistics
#==================================================================================================

echo -e "\n- Information $LOGN: Convert VCF files to Annovar format\n" && LOGN=$((LOGN+1))

for tool in ${TOOLS_USED}; do

    convert_to_annovar variant_calling/$tool/${PREFIX}_${tool}.vcf.gz variant_calling/$tool/${PREFIX}_${tool}.annovar

done

echo -e "\n- Information $LOGN: Combine variants from multiple callers and annotate with AnnotateBamStatistics\n" && LOGN=$((LOGN+1))

# Arguments: 1) Annovar files; 2) Annovar names; 3) BAM; 4) OUTPUT_PREFIX
combine_and_annotate $(for tool in $TOOLS_USED; do echo -n "variant_calling/$tool/${PREFIX}_${tool}.annovar,"; done | sed "s/,$//") $USE ${TUMOR},${NORMAL} ${PREFIX}

if [[ $? -ne 0 ]]
then

	echo 'Error: combine_and_annotate failed'
	exit -1

fi

#--------------------------------------------------------------------------------------------------
#Add pindel data
#--------------------------------------------------------------------------------------------------

# replace with cgpvaf!!!
if [[ $TOOLS_USED =~ "pindel" ]]; then

    /data/roger/scripts/add_vaf_pindel.py -v variant_calling/${PREFIX}.txt -p variant_calling/pindel/${PREFIX}_pindel.vcf.gz -r $REFERENCE > variant_calling/${PREFIX}_pindel.txt

    if [[ $? -ne 0 ]]
    then

        echo -e "Error: Pindel VAF annotation failed"
        echo -e "Command: /data/roger/scripts/add_vaf_pindel.py -v variant_calling/${PREFIX}.txt -p variant_calling/pindel/${PREFIX}_pindel.vcf.gz -r $REFERENCE"
        exit -1

    fi

fi

#--------------------------------------------------------------------------------------------------
#done
#--------------------------------------------------------------------------------------------------

## Compress all the input files and remove invalid_input from Annovar

find -type f -name ${PREFIX}*invalid_input | xargs -I {} rm {}

end=`date +%s`
echo -e "Variant calling was performed in $((end-start)) seconds"

#----------------------------------------------------------------
#Done
#----------------------------------------------------------------

exit
