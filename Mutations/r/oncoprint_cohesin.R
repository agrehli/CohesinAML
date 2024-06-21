
## SOURCES OF ONCOPRINT / WATERFALL PLOT

# MAF FORMAT : https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/

# https://www.biostars.org/p/181159/
# http://bioconductor.org/packages/release/bioc/vignettes/GenVisR/inst/doc/Intro.html
# https://www.bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html#2_generating_maf_files
# https://bioconductor.org/packages/release/bioc/vignettes/ComplexHeatmap/inst/doc/s8.oncoprint.html

# Filtering mutations with MutSigCV: http://software.broadinstitute.org/cancer/software/genepattern/modules/docs/MutSigCV

#--------------------------------------------------------------------------------------------------
# Compare with Francois' analysis
#--------------------------------------------------------------------------------------------------

library(data.table)

## Load data


roger$Review[is.na(roger$Review)] <- 1
roger <- roger[grep("T|AML|UKR|BB",Tumor_Sample_Barcode,invert = T)]

myeloid <- fread('/genomes/hg19/myeloid_panel.txt',col.names='Gene')
roger <- roger[Gene.refGene  %in% myeloid$Gene]

roger[Gene.refGene == 'FLT3' & Ref=='-']$Gene.refGene <- 'FLT3-ITD'
roger[Gene.refGene == 'FLT3' & grepl("exon20",AAChange.refGene)]$Gene.refGene <- 'FLT3 TKD'

fran <- openxlsx::read.xlsx('/data/bas/1_enhancer_analysis/NGS_variant_analyse_final_FK01092018.xlsx')
fran <- as.data.table(fran)
fran2 <- unique(fran[,.(samplenr,Gene,Hotspot_Filter,AAChange,mutatie_final_FK01092018,Filter_new,Var.Perc=as.double(`Var-Perc1`))])

## Combine Roger + Francois 

# Precision or positive predictive value (TP / FP + TP)

merged <- merge(roger,fran2,by.x=c('Tumor_Sample_Barcode','Gene.refGene'),by.y=c('samplenr','Gene'),all.x=T)
# View(merged[,.(Tumor_Sample_Barcode,Gene.refGene,Review,mutatie_final_FK01092018,AAChange)])

same_aa <- unlist(sapply(1:nrow(merged),function(x) { grepl(merged[x]$AAChange,merged[x]$AAChange.refGene) }))

merged2 <- unique(merged[is.na(mutatie_final_FK01092018) | same_aa == 1 ])

true_positive <- merged2[Tumor_Sample_Barcode %in% fran2$samplenr]
true_positive[Gene.refGene == 'FLT3-ITD']$mutatie_final_FK01092018 <- 1 # Francois missed many true FLT3-ITD

true_positive[Review == 1 & mutatie_final_FK01092018 == 1,.N]/true_positive[Review == 1,.N] # 85% precision

true_positive[Review %in% c(0,0.5) & mutatie_final_FK01092018 %in% c(2,3,4),.N]/
  true_positive[mutatie_final_FK01092018 %in% c(2,3,4),.N] # 50% specificity (TN/TN+FP) of my review, NOT the calling

fwrite(merged2,'variants_francois2.txt',sep='\t')

# Sensitivity or true positive rate (TP/TP+FN)

fran.true <- fran2[mutatie_final_FK01092018 == 1]
fran.true <- fran.true[samplenr %in% roger$Tumor_Sample_Barcode]

sensitivity <-  merge(roger,fran.true,by.x=c('Tumor_Sample_Barcode','Gene.refGene'),by.y=c('samplenr','Gene'),all.y=T,allow.cartesian = T)
same_aa2 <- unlist(sapply(1:nrow(sensitivity),function(x) { grepl(sensitivity[x]$AAChange,sensitivity[x]$AAChange.refGene) }))

# fwrite(sensitivity,'sensitivity_francois.txt',sep='\t')

sensitivity <- sensitivity[ same_aa2 == 1 | is.na(AAChange.refGene),] 
sensitivity[Review == 1,.N]/sensitivity[,.N] # 79% (only reviewed as 1)
sensitivity[!is.na(AAChange.refGene),.N]/sensitivity[,.N] # 84.3% (total) --> Most due to low VAF
sensitivity[!is.na(AAChange.refGene) & Var.Perc > 0.15,.N]/sensitivity[Var.Perc > 0.15,.N]  # 97.5% if only variants with high VAF

#--------------------------------------------------------------------------------------------------
# Generate oncoplot
#--------------------------------------------------------------------------------------------------

## IMPORT DATA

setwd('/data/bas/1_enhancer_analysis/exome/variant_analysis/myeloid')
library(maftools)
library(data.table)

metadata <- as.data.table(openxlsx::read.xlsx('/data/bas/1_enhancer_analysis/chip_h3k27ac/metadata/mut_data_20181004.xlsx'))
metadata$label <- gsub("_.*","",metadata$label)
metadata$fab_class <- gsub("(M[0-9]).*","\\1",metadata$fab_class)
metadata$fab_class <- ifelse(metadata$fab_class %in% c("RAEB-t","M0", "M1", "M2", "M3", "M4", "M5","M6"),metadata$fab_class,'Other')
clinicaldata <- metadata[,.(Tumor_Sample_Barcode=label,Sex=sex,FAB_classification=fab_class)]

maf <- annovarToMaf('all_variants_reviewed_diagnostics.txt')
# Remove inconvenient genes
maf <- maf[Review_Diagnostics != 0 & Review_Diagnostics != 0.5]

cohesin = maf[grep("^(SMC1[AB]|SMC3|RAD21|SMC5|STAG1|STAG2|STAG3)$",maf$Hugo_Symbol)][,Hugo_Symbol:='COHESIN']
table(maf[grep("^(SMC1[AB]|SMC3|RAD21|SMC5|STAG1|STAG2|STAG3)$",maf$Hugo_Symbol)]$Hugo_Symbol)
maf = rbind(maf,cohesin)

aml = read.maf(maf=maf,removeDuplicatedVariants=F,clinicalData = clinicaldata)

## GENERATE PLOTS

#Color coding for FAB classification; try getAnnotations(x = laml) to see available annotations.
fabcolors = RColorBrewer::brewer.pal(n = 8,name = 'Spectral')
names(fabcolors) = c("RAEB-t","M0", "M1", "M2", "M3", "M4", "M5","Other")
fabcolors = list(FAB_classification = fabcolors)

png(file='Oncoplot_Exome.png',width=1400,900)
oncoplot(maf = aml, top=25, clinicalFeatures = 'FAB_classification', annotationColor = fabcolors,fontSize = 14)
dev.off()

png(file='Oncostrip_Cohesin.png',width=1400,300)
oncostrip(maf = aml, genes = unique(grep("^(SMC1[AB]|SMC3|RAD21|SMC5|STAG1|STAG2|STAG3|COHESIN)$",maf$Hugo_Symbol,value=T))) 
dev.off()

# Additional plots
lollipopPlot(maf = aml, gene = 'DNMT3A', AACol = 'Protein_Change', refSeqID = 'NM_175629', collapsePosLabel = TRUE, cBioPortal = TRUE)
plotVaf(maf = aml, vafCol = 'sec HQ ratio')

# Advice to prevent errors in repeat regions https://www.biostars.org/p/301059

# Mutations in MDS: SRSF2 - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3837510/
# FLT3:  Although always leading to an in-frame transcript, FLT3-ITDs vary in length (between three to over 400 nucleotides)

## Somatic interactions

png(file='SomaticInteractions_Exome.png',width=1000,1000)
somaticInteractions(maf = subsetMaf(maf=aml,query='Hugo_Symbol!="COHESIN"',mafObj=T), top = 25, pvalue = c(0.05, 0.1),fontSize=1)
dev.off()

## Detecting cancer driver genes based on positional clustering
aml.sig = oncodrive(maf = aml, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = aml.sig, fdrCutOff = 0.1, useFraction = TRUE)

## Pfam domains

aml.pfam = pfamDomains(maf = aml, AACol = 'AAChange', top = 20)

#Run main function with maximum 6 signatures. 
# library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
# aml.tnm = trinucleotideMatrix(maf = aml, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
# library('NMF')
# aml.sign = extractSignatures(mat = aml.tnm, nTry = 6, plotBestFitRes = FALSE)

#--------------------------------------------------------------------------------------------------
# Oncoplot Diagnostics
#--------------------------------------------------------------------------------------------------

fran.filtered <- fran[mutatie_FK22072018 == 1,.(samplenr,Chr,Start,End,Ref,Alt,Func.refGene=Func,Gene.refGene=Gene,
                                                GeneDetail.refGene=Gene,ExonicFunc.refGene=ExonicFunc,AAChange.refGene=AAChange)]
fran.filtered <- fran.filtered[!is.na(Chr)]

fwrite(fran.filtered,file='caca.txt',sep='\t')

maf.diag <- annovarToMaf(annovar = 'caca.txt',tsbCol='samplenr')
cohesin = maf.diag[grep("^(SMC1[AB]|SMC3|RAD21|SMC5|STAG1|STAG2|STAG3)$",maf$Hugo_Symbol)][,Hugo_Symbol:='COHESIN']
table(maf.diag[grep("^(SMC1[AB]|SMC3|RAD21|SMC5|STAG1|STAG2|STAG3)$",maf$Hugo_Symbol)]$Hugo_Symbol)
maf.diag = rbind(maf.diag,cohesin)

aml.diag = read.maf(maf=maf.diag,removeDuplicatedVariants=F)

png(file='Oncoplot_Diagnostics.png',width=1400,900)
oncoplot(maf = aml.diag, top=25, annotationColor = fabcolors,fontSize = 14)
dev.off()

png(file='Oncostrip_Diagnostics_Cohesin.png',width=1400,300)
oncostrip(maf = aml, genes = unique(grep("^(SMC1[AB]|SMC3|RAD21|SMC5|STAG1|STAG2|STAG3|COHESIN)$",maf$Hugo_Symbol,value=T))) 
dev.off()

## Somatic interactions

png(file='SomaticInteractions_Diagnostics.png',width=1000,1000)
somaticInteractions(maf = subsetMaf(maf=aml.diag,query='Hugo_Symbol!="COHESIN"',mafObj=T), top = 25, pvalue = c(0.05, 0.1),fontSize=1)
dev.off()

