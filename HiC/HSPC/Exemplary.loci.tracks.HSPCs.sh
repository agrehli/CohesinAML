#!/bin/bash
#by Alex F
#by Alex F JUL 2022

#########################################################################################################
#       Average ChIP/ATAC/RNA-seq tracks for HiC interaction matrices of HSPC samples                   #
#                                       for loci of interest                                            #
#########################################################################################################

#general paths
DIR_PKG="/misc/software/ngs"
DIR_SOFT="/misc/software"
DIR_PKG="${DIR_SOFT}/ngs"
DIR_DATA="/misc/data"
OS=$(lsb_release -c |grep "^Codename" | awk -F: '{print $2}' | sed 's/[[:blank:]]//g')

#files
HG38GTF="${DIR_PKG}/genome/annotation/GRCh38.PRI_p10/gencode.v27.primary_assembly.annotation.gtf"
GENEGTF="${DIR_DATA}/analysis/generalStuff/annotation/GRCh38/gencode27.genes.labeled.txt"

#directories
WORKDIR_CD34="${DIR_DATA}/analysis/project_cohesin/CD34/HiC_KDs/DeepSeq_Analysis"
MATRIXDIR_CD34="${WORKDIR_CD34}/Interactionmatrices"
FIGURESDIR="${WORKDIR_CD34}/figures"
TMPDIR="/loctmp"
LOOPDIR_CD34="${WORKDIR_CD34}/loops"
INIDIR="${FIGURESDIR}/pygenometracks"


#Variables
KDs="SA1KD SA2KD RAD21KD"
STAGKDs="SA2KD SA1KD"

##new dirs for figures output
mkdir ${FIGURESDIR}/interactionmaps/paper_loci/
mkdir ${FIGURESDIR}/interactionmaps/paper_loci/inclRNA

#####loci considered for paper figures - including ATAC/CTCF/RNAseq tracks
gene_chr_region_size=(
"SOCS2" "chr12" "92855000-94355000" "1.5MB"
"KIF17" "chr1" "19968017-21468016" "1.5MB"
"ADGRA2" "chr8" "37034191-38534190" "1.5MB"
"SLC1A3" "chr5" "35846588-37346587" "1.5MB"
"ILDR2" "chr1" "166225540-167725539" "1.5MB"
"DACT1" "chr14" "57883967-59383966" "1.5MB"
"PAWR" "chr12" "78940964-80440963" "1.5MB"
"ITGA9" "chr3" "36702115-38202114" "1.5MB"
"MYCT1" "chr6" "151947897-153447896" "1.5MB"
"SH2D2A" "chr1" "156066848-157566847" "1.5MB"
)
for KD in ${STAGKDs};do
for (( idx=0 ; idx<${#gene_chr_region_size[@]} ; idx+=4 )) ; do
    gene=${gene_chr_region_size[idx]}
    chr=${gene_chr_region_size[idx+1]}
    region=${gene_chr_region_size[idx+2]}
    size=${gene_chr_region_size[idx+3]}
pyGenomeTracks_v3.5.sh --tracks ${INIDIR}/Pytracks_CTRLvs${KD}.STAG.RAD21.max5.Loops.H3K.CTCF.ATAC.RNA.ini --region ${chr}:${region} --dpi 900 -o ${FIGURESDIR}/interactionmaps/paper_loci/inclRNA/Pytracks_CTRLvs${KD}.STAG.RAD21.H3K.CTCF.ATAC.RNA.loops.${gene}.${size}.png
done
done

#####loci considered for paper figures - versions with other fontsize scale etc
gene_chr_region_size=(
"PAWR" "chr12" "78940964-80440963" "1.5MB"
"SH2D2A" "chr1" "156066848-157566847" "1.5MB"
)
for KD in ${STAGKDs};do
KD="SA1KD"
for (( idx=0 ; idx<${#gene_chr_region_size[@]} ; idx+=4 )) ; do
    gene=${gene_chr_region_size[idx]}
    chr=${gene_chr_region_size[idx+1]}
    region=${gene_chr_region_size[idx+2]}
    size=${gene_chr_region_size[idx+3]}
pyGenomeTracks_v3.5.sh --tracks ${INIDIR}/Pytracks_CTRLvs${KD}.STAG.RAD21.Loops.H3K.CTCF.ATAC.max10.RNA.ini  --region ${chr}:${region} --dpi 900 -o ${FIGURESDIR}/interactionmaps/paper_loci/inclRNA/Pytracks_CTRL.vs${KD}.STAG.RAD21.H3K.CTCF.ATAC.RNA.loops.${gene}.${size}.max10.png
done
done

