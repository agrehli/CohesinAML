# CohesinAML
Pipelines and scripts to reproduce the analysis in Fischer et al. (2024)  

# Description
The repository is structured based on data packages that were analysed partially dependent from each other. This includes the following analyses in the order listed:

- Scripts to analyse mutation screening data (???)
- Scripts and input csv files to analyse and visualize HSPC culture expansion rates, genome editing rates and colony forming untit (CFU) results ([HSPC_cultures folder](HSPC_cultures/))
- Pipelines for mapping ChIP/ATAC/RNA-seq and routine scripts used for several downstream analyses ([Pipelines folder](Pipelines/))
- Scripts to analyse and visualize RNA-seq results after mapping ([RNAseq folder](RNAseq/))
- Scripts to analyse and visualize single cell RNA-seq data ([scRNAseq folder](scRNAseq/))
- Scripts and input csv files to analyse and visualize actin-normalized western blots ([WB_actin_norm folder](WB_actin_norm/))
- Scripts to analyse and visualize ATAC-seq results after mapping ([ATACseq folder](ATACseq/))
- Scripts to analyse and visualize ChIP-seq results after mapping, contains subdirectories by antibody ([ChIPseq folder](ChIPseq/))
- Scripts for mapping, analyzing and visualizing Hi-C data, inlcluding matrix processing via the FAN-C suite ([HiC folder](HiC/)) 
- Scripts to integrate and correlate transcriptomic results with epigentic assays (ATAC/ChIP/HiC)([Epigenetic_GEX_correlations folder](Epigenetic_GEX_correlations/)) 
