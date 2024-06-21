library(openxlsx)
library(data.table)
library(parallel)
library(ComplexHeatmap)
library(trackViewer) # For lollipop plot

source("~/scripts/get_protein_domains.R")

#===============================================================================
# Custom screen analysis
#===============================================================================

madex <- fread("/data/bas/1_enhancer_analysis/metadata/madex_data_all.csv")
patients.panel <- fread("/data/stanley/cohesin_panel/tested_patients.txt",header=F)

## Miseq data
miseq <- data.table(read.xlsx("/data/bas/1_enhancer_analysis/metadata/TSM_mutaties_compleet_diagnose_v27012020.xlsx"))
miseq <- miseq[,.(label=samplenr,Gene=genetarget,Func,ExonicFunc,AAChange,VAF=vaf,source='Illumina')]

## Cohesin panel

structural.genes <- fread("/data/claudia_gebhard/cohesins/mutations/structural_genes.txt",header = F)[V1 != "PHF6"]
cohesin.panel <- data.table(read.xlsx("/data/claudia_gebhard/cohesins/mutations/summary_variants_filtered.xlsx"))

cohesin.panel <- cohesin.panel[Gene.refGene %in% c(structural.genes$V1,unique(miseq$Gene)) & Artifact == "No"]
cohesin.panel[Gene.refGene == 'FLT3' & Ref=='-']$Gene.refGene <- 'FLT3-ITD'
cohesin.panel[Gene.refGene == 'FLT3' & grepl("exon20",AAChange.refGene)]$Gene.refGene <- 'FLT3 TKD'
cohesin.panel <- cohesin.panel[,.(label=Tumor_Sample_Barcode,Gene=Gene.refGene,Func=Func.refGene,
                                  ExonicFunc=ExonicFunc.refGene,AAChange=AAChange.refGene,VAF=alt.freq,source='Cohesin')]

## Combine data

merged <- rbind(cohesin.panel,miseq[label %in% patients.panel$V1])
merged <- merged[,if( all(c("Cohesin","Illumina") %in% .SD$source) ){.SD[source == "Illumina"]}else{.SD},by=c("label","Gene")]

merged$class <- ifelse(grepl("splicing",merged$Func),
                              "splicing",merged$ExonicFunc)
merged$class <- gsub(".*frameshift.*","indel",merged$class)
merged$class <- gsub("nonsynonymous SNV","missense_mutation",merged$class)
merged$class <- gsub("stopgain|startloss|stoploss","nonsense_mutation",merged$class)
merged <- merged[class != "unknown"] # Combine panel and Miseq data

## Add cytogenetics from MADEX

molten.madex <- melt(madex[,.(label,inv16,`7q`,`11q23`,t8_21)],id.vars='label',variable.name = 'Gene')
molten.madex <- molten.madex[value == "pos" & label %in% patients.panel$V1]
molten.madex$class <- "structural_variant"

merged.madex <- rbind(merged[,.(label,Gene,class)],molten.madex[,.(label,Gene,class)])

#-------------------------------------------------------------------------------
# Color functions for Oncoprint
#-------------------------------------------------------------------------------

pal.oncoplot <- MetBrewer::met.brewer("Lakota",5)
col = c(indel = pal.oncoplot[1], missense_mutation = pal.oncoplot[2],nonsense_mutation=pal.oncoplot[3],
        splicing = pal.oncoplot[4],unknown="black",structural_variant=pal.oncoplot[5],mutation="firebrick")
# We can change the thickness of the bars by altering the w-unit parameter
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.65, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big red
  indel = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.65, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = pal.oncoplot[1], col = NA))
  },
  mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.65, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "firebrick", col = NA))
  },
  # small blue
  exonic_mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.65, "mm"), h*0.33, 
              gp = gpar(fill = "steelblue", col = NA))
  },
  # small blue
  missense_mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.65, "mm"), h*0.33, 
              gp = gpar(fill = pal.oncoplot[2], col = NA))
  },
  # small blue
  nonsense_mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.65, "mm"), h*0.33, 
              gp = gpar(fill = pal.oncoplot[3], col = NA))
  },
  # small green
  splicing = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.65, "mm"), h*0.33, 
              gp = gpar(fill = pal.oncoplot[4], col = NA))
  },
  noncoding = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.65, "mm"), h*0.33, 
              gp = gpar(fill = "black", col = NA))
  },
  structural_variant = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.65, "mm"), h*0.33, 
              gp = gpar(fill = pal.oncoplot[5], col = NA))
  }
)

#-------------------------------------------------------------------------------
# Generate oncoprint
#-------------------------------------------------------------------------------

muts <- dcast(merged.madex, Gene ~ label,value.var="class",
              fun.aggregate=function(x) head(x, 1),fill="")
MAT <- as.matrix(muts,rownames='Gene')
mutcounts <- apply(MAT,1,function(x) sum(x != ""))
MAT <- MAT[(mutcounts / ncol(MAT))*100 >=2.5,] # Remove muts below 2.5%
mutcounts <- apply(MAT,1,function(x) sum(x != ""))
MAT <- MAT[order(mutcounts,decreasing = T),]
# Add a row without any STAG1 mutations
MAT <- rbind(MAT,matrix("",nrow=1,ncol=ncol(MAT),dimnames=list("STAG1")))
# Add columns with patients that have 0 mutations
MAT <- cbind(MAT,matrix("",ncol=sum(!patients.panel$V1 %in% colnames(MAT)),nrow=nrow(MAT),
       dim=list(NULL,patients.panel$V1[!patients.panel$V1 %in% colnames(MAT)])))

split.cohesins <- factor(ifelse(rownames(MAT) %in% structural.genes$V1,"Cohesin","Other"),
                         levels=c("Cohesin","Other"))
order.cohesins <- c(which(split.cohesins == "Cohesin"),which(split.cohesins == "Other"))

pdf("/data/claudia_gebhard/cohesins/oncoplot_unbiased.pdf",width=12,height=8)
draw(oncoPrint(MAT, alter_fun = alter_fun, col = col,row_order=1:nrow(MAT),
               column_title = sprintf("OncoPrint (n = %s)",ncol(MAT)),
               pct_digits=2,row_names_gp = gpar(col= ifelse(rownames(MAT) == "STAG2","red","black"),fontface='bold'))
     )
dev.off()

pdf("/data/claudia_gebhard/cohesins/oncoplot_split_all.pdf",width=15/2,height=16/2)
draw(oncoPrint(MAT, alter_fun = alter_fun, col = col,row_order=order.cohesins,
               column_title = sprintf("Mutation screening (n = %s)",ncol(MAT)),name="Alterations",
               column_title_gp=gpar(fontface='bold'),
               row_split=split.cohesins,pct_digits=2,border=T,remove_empty_columns = F,
               row_names_gp = gpar(col= ifelse(rownames(MAT) == "STAG2","red","black"),fontface='bold',fontsize=10),
               heatmap_legend_param=list(nrow=2),row_title_gp=gpar(fontface='bold')),
     heatmap_legend_side = "bottom")
dev.off()

fwrite(MAT,"/data/claudia_gebhard/cohesins/analyses/mutation_matrix.txt",row.names=F,col.names=T,sep='\t')

#-------------------------------------------------------------------------------
# Show association between genes - small cohort
#-------------------------------------------------------------------------------

# Good example here: https://www.nature.com/articles/nature10351

muts2 <- dcast(merged.madex,label ~ Gene,value.var='class',fun.aggregate=function(x)ifelse(length(x)>0,1,0))
sample_id <- unlist(muts2[,1])
intMat <- as.matrix(muts2 [,-1])

interaction_pairs <- function(skMat,sample_id,LOG=T) {
  
  s <- mcmapply(function(i) {
    message(i)
    interactions <- mapply(function(j) {
      # If we provide the two vectors separately to the F-test, all four categories need to be present (00,01,10,11). However,
      # due to NA it may happen that there are no patients with all these categories for a given combination of genes
      f <- try(fisher.test(table(skMat[,i], skMat[,j])), silent=TRUE);
      if(LOG == T){
        if(class(f)=="try-error") NA else ifelse(f$estimate>1, -log10(f$p.val),log10(f$p.val)) 
        pval <- setNames(round(if(class(f)=="try-error") NA else ifelse(f$estimate>1, -log10(f$p.val),log10(f$p.val)),4),'p.value')  
      } else {
        if(class(f)=="try-error") NA else ifelse(f$estimate>1, -f$p.val,f$p.val) 
        pval <- setNames(round(if(class(f)=="try-error") NA else ifelse(f$estimate>1,f$p.val,f$p.val),4),'p.value')
      }
      pairs <- factor(as.data.table(skMat)[,c(i,j),with=FALSE][,paste0(get(colnames(skMat)[i]),get(colnames(skMat)[j]))],
                      levels=c('00','01','10','11'))
      samples <- paste0(na.omit(sample_id[pairs == "11"]),collapse=",")
      pairs <- table(pairs)
      oddsRatio <-  unname(if(class(f)=="try-error") f=NA else f$estimate)
      c(gene1=colnames(skMat)[i],gene2=colnames(skMat)[j],pval,oddsRatio=oddsRatio,pairs,samples=samples)
    },1:ncol(skMat),SIMPLIFY = FALSE,USE.NAMES = TRUE)

    m <- do.call(rbind,interactions)
    as.data.table(m)
    
  },1:ncol(skMat),SIMPLIFY = F,USE.NAMES = T,mc.cores= parallel::detectCores()-2) 
  return(s)  
}

lolo <- interaction_pairs(intMat,sample_id)

# int.pairs.new <- rbindlist(int.pairs.new)

jeje <- do.call(rbind,lolo)
jeje$p.value <- as.numeric(jeje$p.value)

heatmap.stag2 <- ggplot(jeje[gene1=="STAG2" & gene2!="STAG2"],aes(x=gene2,y=1,fill=p.value)) + 
  geom_tile(color = "white",lwd = 1,linetype = 1) + coord_fixed() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_gradientn(colors=c("blue","grey","red"),breaks=c(-3,0,3),
                       limits=c(-2,2),oob = scales::squish) +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  geom_point(aes(size=abs(p.value) > -log10(0.05)),shape=1) + 
  scale_size_manual(values=c(-3,3)) +
  labs(x="Genes mutated in AML",fill="log10(p.value)",size="Significant")
  
  
  # geom_point(data=jeje[gene1=="STAG2" & abs(p.value) > -log10(0.05)],
  #            aes(x=gene2,y=1),size=0.5,show.legend = T,shape=4) +

ggsave(file="/data/claudia_gebhard/cohesins/heatmap_stag2.pdf",heatmap.stag2,width=12,height=4)

#-------------------------------------------------------------------------------
# Male/female bias of mutations
#-------------------------------------------------------------------------------

pooled <- data.table(read.xlsx("/data/bas/1_enhancer_analysis/metadata/Pooled_AML_database_mol_data_v11112020_FK23022021.xlsx"))
pooled$sex <- c("F"="female","M"="male")[pooled$sex]

genders <- fread("/data/stanley/cohesin_panel/gender_samples.txt",col.names=c("label","reads","gender"))
genders$label <- sapply(genders$label,function(x)strsplit(x,'_')[[1]][1])

miseq <- data.table(read.xlsx("/data/bas/1_enhancer_analysis/metadata/TSM_mutaties_compleet_diagnose_v27012020.xlsx"))
cohesin.panel2 <- merge(cohesin.panel,unique(miseq[,.(label=as.numeric(samplenr), pid=`trial-pid`)]),by="label",all.x=T)
cohesin.panel2 <- merge(cohesin.panel2,pooled[,.(pid=`trialnr-pid`,sex)],by='pid',all.x=T)

cohesin.panel2 <- merge(cohesin.panel2,unique(madex[,.(label=as.numeric(label),sex.madex=sex)]),by="label",all.x=T)
cohesin.panel2 <- merge(cohesin.panel2,genders[,.(label=as.numeric(label),sex.ngs=gender)],by="label")
cohesin.panel2$sex.madex[cohesin.panel2$sex.madex == "NULL"] <- NA

cohesin.panel2[,final_sex := as.factor(ifelse(is.na(sex),ifelse(is.na(sex.madex),sex.ngs,sex.madex),sex))]

p.sex.small <- ggplot(cohesin.panel2[Gene %in% c("STAG2","RAD21")],aes(x=final_sex,y=VAF)) + 
  geom_violin(aes(color=final_sex)) + stat_summary(aes(color=final_sex),fun.data = "mean_sdl",  
          fun.args = list(mult = 1),  geom = "pointrange") +
  geom_jitter(size=0.5) + facet_wrap(~Gene) + theme_cowplot() +
  ggsci::scale_color_npg(guide=F) + labs(x="") + 
  geom_signif(comparisons = list(c("male","female")),test='t.test')

ggsave(p.sex.small,filename="/data/claudia_gebhard/cohesins/analyses/VAFbias_plot_small.pdf",width=5,height=5)

# Calculate significance
casted.small <- dcast(cohesin.panel2[!is.na(sex)],formula=label+sex~Gene,fun.aggregate = function(x)ifelse(length(x)>0,"mut","wt"))

test.sex <- function(casted) {
    res <- mapply(FUN=function(i){
    contingency <- table(casted[,c(i,"sex"),with=F])
    f <- try(fisher.test(contingency), silent=TRUE)
    pval <- if(class(f)=="try-error") NA else f$p.value
    oddsRatio <-  unname(if(class(f)=="try-error") f=NA else f$estimate)
    pairs <- factor(casted[,c(i,"sex"),with=F][,paste0(get(i),'_',sex)],
                    levels=c('wt_male','wt_female','mut_male','mut_female'))
    c(pval=pval,oddsRatio=oddsRatio,table(pairs))
  },colnames(casted[,-c("label","sex")]))
    return(data.table(t(res),keep.rownames = 'Gene'))
}
  
# Sanity check

sig.bias.small <- test.sex(casted.small)
unique(cohesin.panel2[!is.na(sex)][!label %in% cohesin.panel2[Gene=="STAG2"]$label][,.(final_sex,label)])[,.N,by=final_sex]

# Alternatives would include table() or dt[,,by=]. The latter does not count empty observations.
# To make it work in dplyr we need factors in the 'group_by'
# https://stackoverflow.com/questions/25956178/proper-idiom-for-adding-zero-count-rows-in-tidyr-dplyr
sex.counts.panel <- cohesin.panel2 %>% group_by(Gene,final_sex,.drop=F) %>% 
  summarise(cnt = n(),.groups="drop") %>% group_by(Gene) %>% mutate(perc = cnt*100/sum(cnt)) %>% 
  merge(sig.bias.small,by='Gene',all.x=T) %>% arrange(final_sex,perc)

p.sexbias.small <- sex.counts.panel %>% mutate(Gene = factor(Gene,levels=unique(sex.counts.panel$Gene))) %>%
 ggplot(aes(x=Gene,y=perc,fill=final_sex)) + geom_bar(position='stack',stat='identity') +
         theme_cowplot() + ggsci::scale_fill_npg(name="sex") +
  geom_text(aes(y=101,label = ifelse(pval < 0.05, "*", ""),color=ifelse(oddsRatio > 1,"female","male")),size=5,show.legend=F) +
  labs(x="",y="Percentage of cases") + theme(axis.text.x = element_text(angle = 45,hjust=1,size=8)) +
  geom_hline(yintercept=50,linetype=2)  

ggsave(p.sexbias.small,filename="/data/claudia_gebhard/cohesins/analyses/Sexbias_plot_small.pdf",width=6,height=5)

#===============================================================================
# Show association between genes - large cohort
#===============================================================================

muts.large <- fread("/data/bas/1_enhancer_analysis/metadata/TSM_mutaties_diagnose_pos_neg_v18122020.csv")
muts.large[,52] <- NULL # Remove second BRAF column that contains wrong results
for (col in names(muts.large)[2:ncol(muts.large)]) set (x = muts.large,j = col,value = ifelse(muts.large[[col]] == "pos",1,0))

sample_id.large <- unlist(muts.large[,1])
intMat.large <- as.matrix(muts.large[,-1])

int.pairs.large <- interaction_pairs(intMat.large,sample_id.large)

int.pairs.large2 <- do.call(rbind,int.pairs.large )
for (j in c("00","01","10","11","p.value",'oddsRatio')) {set(int.pairs.large2,j=j,value=as.numeric(int.pairs.new[[j]]))}

col_fun = colorRamp2(c(-10, 0, 10), c("blue", "grey", "red"))

heatmap_pairs_gene <- function(gene,int.pairs.large2) {
  
  mat.pairs <- t(data.frame(na.omit(int.pairs.large2[gene1==gene & gene2!=gene,.(gene2,p.value)]),row.names='gene2'))
  counts.samples <- na.omit(int.pairs.large2[gene1==gene & gene2!=gene,.(gene2,p.value,`11`)])
  heat.pairs <- Heatmap(mat.pairs,show_column_dend = F,col=col_fun,show_row_names=F,show_heatmap_legend = F,
                        rect_gp = gpar(col = "white", lwd = 2),column_names_gp = gpar(fontsize=8),column_names_rot = 45,
                        top_annotation = HeatmapAnnotation(`# samples` = anno_barplot(as.numeric(counts.samples$`11`)),
                                                           annotation_name_gp = gpar(fontsize=8)),
                        column_title=sprintf("%s co-mutation analysis (n = %s)",gene,nrow(intMat.large)),
                        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                        #       heatmap_legend_param = list( direction = "horizontal"),height=unit(0.5,"cm"),width=unit(20,"cm"),
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          if (abs(mat.pairs[i,j]) > 2) grid.circle(x = x, y = y, r = width*9, gp = gpar(fill = "black"))
                          else if (abs(mat.pairs[i,j]) > 1.3) grid.circle(x = x, y = y, r = width*9, gp = gpar(fill = NA))}
  )
  return(heat.pairs)
}

lgd = Legend(col_fun = col_fun, title = "-log(p-value)",direction = "horizontal")
lgd.dots = Legend(labels = c("p-value < 0.05","p-value < 0.01"),direction='vertical',
             graphics = list(
               function(x, y, w, h) grid.circle(x, y, r=w/4, gp = gpar(col = "black",fill=NA)),
               function(x, y, w, h) grid.circle(x, y, r=w/4,gp = gpar(col = "black",fill="black"))
             ))

heat.pairs.STAG2 <- heatmap_pairs_gene("STAG2",int.pairs.large2)
heat.pairs.RAD21 <- heatmap_pairs_gene("RAD21",int.pairs.large2)

# We cannot draw it as a list of heatmaps because then columns are aligned to be in the same 
# order as in the first heatmap
p1 <- grid::grid.grabExpr(draw(heat.pairs.STAG2)); p2 <- grid::grid.grabExpr(draw(heat.pairs.RAD21))
p3 <- grid::grid.grabExpr(draw(packLegend(lgd,lgd.dots,direction='horizontal')))
#p2 <- grid::grid.grabExpr(draw(heat.pairs.RAD21,annotation_legend_list = packLegend(lgd,lgd.dots,direction='horizontal'),annotation_legend_side='bottom'))

pdf(file="/data/claudia_gebhard/cohesins/heatmap_STAG2_RAD21.pdf",width=12,height=5)
plot_grid(p1,p2,p3,ncol=1)
dev.off()

heat.pairs.large2 <- ggplot(na.omit(int.pairs.large2[gene1=="STAG2" & gene2!="STAG2"]),aes(x=gene2,y=1,fill=p.value)) + 
  geom_tile(color = "white",lwd = 1,linetype = 1) + coord_fixed() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_gradientn(colors=c("blue","grey","red"),breaks=c(-5,0,5),
                       limits=c(-2,2),oob = scales::squish) +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  geom_point(aes(size=abs(p.value) > -log10(0.05)),shape=1) + 
  scale_size_manual(values=c(-3,3)) +
  labs(x="Genes mutated in AML",fill="log10(p.value)",size="Significant")

## Save results

int.pairs.new <- int.pairs.large2

int.pairs.output <- int.pairs.new[gene2 %in% c("STAG2","RAD21"),.(gene1,gene2,p.value,mut_wt=`01`,mut_wt_freq=`01`*100/(`01`+`00`),
                                                     mut_mut=`11`,mut_mut=`11`*100/(`11`+`10`) )]

write.xlsx(int.pairs.output,"/data/claudia_gebhard/cohesins/table_comutation.xls")

# Frequency of mutations in AML
mapply(FUN=function(gene){length(unique(miseq[Gene == gene]$samplenr))*100/length(unique(miseq$samplenr)) },
       c("STAG2","RAD21","SMC1A","SMC3"),USE.NAMES = T)

#-------------------------------------------------------------------------------
# Association ASXL1 / SRSF2 in STAG2-mutant AML
#-------------------------------------------------------------------------------

# ASXL1 and SRSF2 are associated in STAG2-mutant AML
asxl1.srsf2 <- muts.large[STAG2 == 1][,.(ASXL1,SRSF2)]

asxl1.srsf2.tab <- table(as.matrix(asxl1.srsf2[,1]), as.matrix(asxl1.srsf2[,2]))
dimnames(asxl1.srsf2.tab) <- list(c("ASXL1wt","ASXL1mut"),c("SRSF2wt","SRSF2mut"))

f.test <- fisher.test(asxl1.srsf2.tab)

slices <- as.vector(asxl1.srsf2.tab)
lbls <- c("ASXL1wt/SRSF2wt", "ASXL1mut/SRSF2wt", "ASXL1wt/SRSF2mut", "ASXL1mut/SRSF2mut")
pct <- round(slices/sum(slices)*100)
pct <- paste(pct,"%",sep="")
#lbls <- paste(lbls, " (",pct,"%)",sep="" ) # add percents to labels
pdf(file="/data/claudia_gebhard/cohesins/ASXL1_SRSF2_pie.pdf",width=7,height=7)
pie(slices,labels = pct, col=ggsci::pal_nejm()(4),cex=0.8,
    main="ASXL1 and SRSF2 mutations in STAG2-mutated AML")
legend("topright", lbls, cex = 0.65,
       fill = ggsci::pal_nejm()(4))
dev.off()

## Oncoplot

MAT.asxl1 <- t(muts.large[STAG2 == 1][,.(ASXL1,SRSF2,STAG2)])
MAT.asxl1[MAT.asxl1 == 1] <- "mutation"
MAT.asxl1[MAT.asxl1 == 0] <- ""

pdf(file="/data/claudia_gebhard/cohesins/ASXL1_SRSF2_oncoplot.pdf",width=10,height=1)
draw(oncoPrint(MAT.asxl1,alter_fun = alter_fun,col=col,top_annotation = NULL,
               column_title="Co-occurrence of SRSF2 and ASXL1 mutations in STAG2-mutated patients",
               row_names_gp = gpar(col= "black",fontface='bold'))
)
dev.off()

## Waffle plot

library(waffle)

pdf(file="/data/claudia_gebhard/cohesins/ASXL1_SRSF2_waffle.pdf",width=7,height=7)
waffle(setNames(slices,lbls), rows = 8,colors= ggsci::pal_nejm()(4),
       title="ASXL1 and SRSF2 mutations in STAG2-mutated AML",
       xlab=sprintf("1 square = %s patient",1))
dev.off()

## Barplot

# https://rstudio-pubs-static.s3.amazonaws.com/291083_2b0374fddc464ed08b4eb16c95d84075.html

pdf(file="/data/claudia_gebhard/cohesins/ASXL1_SRSF2_barplot.pdf",width=7,height=7)
ggplot(data=asxl1.srsf2[,.(merge=paste0(ASXL1,SRSF2))],aes(x=merge)) +
  geom_bar(aes(y = ..prop.., fill = factor(..x..),group=1), stat="count") +
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ,group=1), stat= "count", vjust = -.5) + 
  scale_x_discrete(limits=c("00","11","01","10"),labels=c("ASXL1wt/SRSF2wt","ASXL1mut/SRSF2mut","ASXL1wt/SRSF2mut","ASXL1mut/SRSF2wt")) +
  labs(y="Percentage of patients",x="") + scale_fill_wsj(guide="none") + theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  scale_y_continuous(labels = scales::percent_format()) +  ylim(c(0,0.5)) +
  ggtitle("ASXL1 and SRSF2 mutations in STAG2-mutated AML")
dev.off()

#-------------------------------------------------------------------------------
# Association ASXL1 / SRSF2 in other AMLs
#-------------------------------------------------------------------------------

# We can see that a small percentage of STAG2-WT have ASXL1/SRSF2, but that only
# indicates that STAG2-mut are indeed enriched
slices.wt <- as.vector(asxl1.srsf2.wt.tab)
lbls <- c("ASXL1wt/SRSF2wt", "ASXL1mut/SRSF2wt", "ASXL1wt/SRSF2mut", "ASXL1mut/SRSF2mut")
pct.wt <- round(slices.wt/sum(slices.wt)*100)
pct.wt <- paste(pct.wt,"%",sep="")
#lbls <- paste(lbls, " (",pct,"%)",sep="" ) # add percents to labels
pdf(file="/data/claudia_gebhard/cohesins/mutations/ASXL1_SRSF2_STAG2wt_pie.pdf",width=7,height=7)
pie(slices.wt,labels = pct.wt, col=ggsci::pal_nejm()(4),cex=0.8,
    main="ASXL1 and SRSF2 mutations in STAG2-wild type AML")
legend("topright", lbls, cex = 0.65,
       fill = ggsci::pal_nejm()(4))
dev.off()

# Thus, better to depict what percentage of ASXL1-mut also have SRSF2 either in
# STAG2-mut or STAG2-wt AML

prop.SRSF2 <- muts.large[,.(ASXL1,STAG2,SRSF2)][,.(prop_SRSF2=sum(SRSF2)*100/.N),by=.(STAG2,ASXL1)]
prop.SRSF2$ASXL1 <- ifelse(prop.SRSF2$ASXL1 == 1,"Mutated","WT")
prop.SRSF2$STAG2 <- ifelse(prop.SRSF2$STAG2 == 1,"STAG2-mutant","STAG2-WT")

ggplot(prop.SRSF2,aes(x=ASXL1,y=prop_SRSF2)) + geom_bar(stat='identity') + 
  facet_wrap(~STAG2) + theme_cowplot() + labs(y="% of cases with SRSF2 mutations")

## To find if STAG2 is more associated with combination of mutations than WT:
muts.reduced <- muts.large[,.(STAG2,ASXL1,SRSF2)]
for (j in colnames(muts.reduced)) set(muts.reduced, j = j, value = factor(muts.reduced[[j]]))
for(i in colnames(muts.reduced)) { levels(muts.reduced[[i]]) <- c("WT","MUT") }

tab.STAG2.comutated <- muts.reduced[,.(STAG2,ASXL1_SRSF2=paste0(ASXL1,'_',SRSF2))][
  ,.N,by=.(STAG2,ASXL1_SRSF2)][,perc:=N*100/sum(N),by='STAG2']
tab.STAG2.comutated$STAG2 <- ifelse(tab.STAG2.comutated$STAG2 =="MUT","STAG2-mutant","Other AML")
# tab.STAG2.comutated <- tab.STAG2.comutated %>% group_by(STAG2) %>% mutate(percent = N*100/sum(N)) %>% arrange(STAG2)

pdf("/data/claudia_gebhard/cohesins/mutations/comparison_ASXL1SRSF2_WT_others.pdf",width=8,height=5)
ggplot(tab.STAG2.comutated,aes(x=ASXL1_SRSF2,y=perc,fill=ASXL1_SRSF2)) + geom_bar(stat='identity') + 
  facet_wrap(~STAG2) + theme_cowplot() + ggsci::scale_fill_d3() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + labs(y="Percentage of cases") + 
  geom_text(aes( label = round(perc,1), y=perc), vjust = -.5) + ylim(c(0,100))
dev.off()

# Ideally we would like to test whether there are significant differences between associations
# and then do a post-hoc test like Tukey HSD with ANOVA. But since that does not exist, we 
# compute chi.square and look at the residuals, measuring how much a variable contributes to pvalue
# https://stats.stackexchange.com/questions/70107/fishers-exact-test-in-r-2x4-table

# Chi-square
df.tab <- data.frame(dcast(tab.STAG2.comutated[,.(STAG2,ASXL1_SRSF2,N)],formula=STAG2~ASXL1_SRSF2),row.names='STAG2')
chisq.test(df.tab)$residuals

# Fisher's exact test
fisher.test(matrix(unlist(df.tab[1,]),nrow=2,ncol=2)) # Other
fisher.test(matrix(unlist(df.tab[2,]),nrow=2,ncol=2)) # STAG2-mutant

#===============================================================================
# Analysis Miseq data
#===============================================================================

miseq <- data.table(read.xlsx("/data/bas/1_enhancer_analysis/metadata/TSM_mutaties_compleet_diagnose_v27012020.xlsx"))

miseq$class <- ifelse(grepl("splicing",miseq$Func),
                       "splicing",miseq$ExonicFunc)
miseq$class <- gsub(".*frameshift.*","indel",miseq$class)
miseq$class <- gsub("nonsynonymous SNV","missense_mutation",miseq$class)
miseq$class <- gsub("stopgain|startloss|stoploss","nonsense_mutation",miseq$class)
miseq <- miseq[class != "unknown"] # Combine panel and Miseq data

muts.mi <- dcast(miseq[Gene %in% c("STAG2","RAD21","SMC3","SMC1A")], Gene ~ samplenr ,value.var="class",
              fun.aggregate=function(x) head(x, 1),fill="")
MAT.mi <- as.matrix(muts.mi,rownames='Gene')
mutcounts <- apply(MAT.mi,1,function(x) sum(x != ""))
MAT.mi <- MAT.mi[(mutcounts / ncol(MAT.mi))*100 >=2.5,] # Remove muts below 2.5%
mutcounts <- apply(MAT.mi,1,function(x) sum(x != "")) # Again to reorder
MAT.mi <- MAT.mi[order(mutcounts,decreasing = T),]
# Add columns with patients that have 0 mutations
MAT.mi <- cbind(MAT.mi,matrix("",ncol=sum(!unique(miseq$samplenr) %in% colnames(MAT.mi)),nrow=nrow(MAT.mi),
                        dim=list(NULL,unique(miseq$samplenr)[!unique(miseq$samplenr) %in% colnames(MAT.mi)])))

pdf("/data/claudia_gebhard/cohesins/oncoplot_bigcohort.pdf",width=14,height=2)
draw(oncoPrint(MAT.mi, alter_fun = alter_fun, col = col,
               column_title = sprintf("Mutation screening (n = 2765)"),name="Alterations",
               column_title_gp=gpar(fontface='bold'),top_annotation=NULL,
               pct_digits=1,border=T,remove_empty_columns = F,
               row_names_gp = gpar(col= ifelse(rownames(MAT.mi) == "STAG2","red","black"),fontface='bold',fontsize=10),
               heatmap_legend_param=list(nrow=2),row_title_gp=gpar(fontface='bold')),
     heatmap_legend_side = "bottom")
dev.off()


#-------------------------------------------------------------------------------
# Distribution of VAFs for STAG2, ASXL1, NPM1 and SRSF2 in patients
#-------------------------------------------------------------------------------

# The muts.large file does not contain the VAFs, but this is largely the same data
# contained in the miseq file. All cases with STAG mutations are in the miseq dataset:

sum(!muts.large[STAG2 == 1]$`hovon-pid` %in% miseq$`trial-pid`)

miseq.stag2 <- miseq[samplenr %in% miseq[genetarget == "STAG2"]$samplenr][
  genetarget %in% c("STAG2","SRSF2","ASXL1","BCOR","CEBPA","RUNX1","IDH2","IDH1","TET2",
                    "RAD21","GATA2","FLT3-ITD","KIT","DNMT3A","NPM1","FLT3-ITD",
                    "SMC1A","EZH2","SMC3")]
miseq.stag2 <- miseq.stag2[,.(label=samplenr,`trial-pid`,Gene=genetarget,aachange,Func,vaf)]

# Separate columns for cases with multiple mutations in one gene
miseq.stag2[,Gene2 := make.unique(.SD[['Gene']]),by='label']
casted <- dcast(miseq.stag2[,.(label,`trial-pid`,Gene2,vaf)],formula = label+`trial-pid` ~ Gene2,value.var="vaf",
                fun.aggregate=mean,fill=as.numeric(NA))
# Take only the max VAF, which is what we want to determine clonality
casted2 <- dcast(miseq.stag2[,.(label,`trial-pid`,Gene,vaf)],formula = label+`trial-pid` ~ Gene,value.var="vaf",
                fun.aggregate=max,fill=as.numeric(NA))

table(miseq.stag2[Gene == 'STAG2']$label)
# casted[label == 24881] 
miseq.stag2$in_plot <- miseq.stag2$`trial-pid` %in% muts.large$`hovon-pid`

openxlsx::write.xlsx(list(miseq.stag2,casted,casted2),"/data/claudia_gebhard/cohesins/STAG2_comutations.xlsx",
                     na.string="",overwrite = T,sheetName=c("AllVariants","PerSample","MaxVAF"),
                     firstRow = TRUE,headerStyle=createStyle(textDecoration='bold'))

#fwrite(casted2,file="/data/claudia_gebhard/cohesins/STAG2mutAML_comutations.tsv",sep='\t')

#-------------------------------------------------------------------------------
# Splicing mutations in STAG2
#-------------------------------------------------------------------------------

openxlsx::write.xlsx(miseq[Gene == "STAG2" & Func == 'splicing',.(samplenr,Gene,Begin,End,Ref,Var,Func,`Var-Perc1`)],
                     "/data/claudia_gebhard/cohesins/STAG2_splicing.xlsx", firstRow = TRUE,headerStyle=createStyle(textDecoration='bold'))

#-------------------------------------------------------------------------------
# Lollipop plot mutations in STAG2
#-------------------------------------------------------------------------------

## MAFtools ##
 # Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Variant_Classification, Variant_Type and Tumor_Sample_Barcode.

## Trackviewer ##

make_lolliplot <- function(gene,isoform,mutations,domains,len_protein,...) {
  
  muts <- data.table(aachange=gsub(sprintf(".*%s.*?p.([A-Z]?[0-9_]+[A-Za-z]+).*",isoform),"\\1",
                                           mutations[Gene == gene ]$AAChange),
                     position=as.numeric(gsub(sprintf(".*%s.*?p.[A-Z]?([0-9]+).*",isoform),"\\1",
                                              mutations[Gene == gene ]$AAChange,perl=T)),
                     type=mutations[Gene == gene ]$ExonicFunc)
  # NOTE: Some mutations of different types share the same position, so they have been counted separately
  # For example, Y1044fs results from either deletion or insertion
  muts[!is.na(position),number := .N,by=.(aachange,type)]        
  muts.unique <- unique(muts[!is.na(position)][order(position)])   
  
  pal.lolliplot <- MetBrewer::met.brewer("Signac",length(unique(na.omit(mutations$ExonicFunc))))
  colors.lolli <- setNames(c(pal.lolliplot),sort(unique(na.omit(mutations$ExonicFunc))))
  # With parameter.draw, only the lollipops (if jitter=node) will be jittered. The rest will behave
  # as with the jitter=label parameter. It is not necessary to specify which names should be plotted
   
  mut.pos <- GRanges("chr1", IRanges(muts.unique$position, width=1, 
                                           name=muts.unique$aachange,
                                           dashline.col=ifelse(muts.unique$number > 2,
                                                               "black",NA),
                                          label.parameter.draw=ifelse(muts.unique$number > 2,
                                                         TRUE,FALSE),
                                           score=muts.unique$number,
                                           color=colors.lolli[muts.unique$type]))
  mut.pos$label.parameter.rot <- 45
  yaxis <- seq(0,max(mut.pos$score),2) # Custom y axis to avoid decimals
  lolliplot(mut.pos,domains,cex=0.8,ylab="# mutations",
            ranges=GRanges("chr1",IRanges(start=1,end=len_protein)), # Length protein
            xaxis=c(1,seq(200,len_protein-200,by=200),len_protein), # Custom x axis
            yaxis=setNames(yaxis,yaxis),
            legend=list(labels=names(colors.lolli),fill=colors.lolli),...)

}

# STAG2

features_STAG2 <- get_features_interpro(geneName_to_uniProt("STAG2"))

domains <- as(features_STAG2[type == "domain",.(chr="chr1",start,end,name,
              fill=c("#FF8833","#51C6E6"),height=c(0.03,0.06))],"GRanges")
names(domains) <- c("STAG","SCD")
#GenomicRanges::setdiff(GRanges("chr1",IRanges(start=1,end=features_STAG2$protein_length[1])),domains)

# We remove stoploss to prevent it from being included in the legend
# The legend is made from all mutations to make sure STAG2 and RAD21 are the same
pdf("/data/claudia_gebhard/cohesins/STAG2_Lolliplot_JitterLabel.pdf",width=20.4,height=4.2)
make_lolliplot("STAG2","NM_006603",miseq,
               domains,len_protein=features_STAG2$protein_length[1],
               label_on_feature=F,jitter='node')
grid.text("STAG2 (NM_006603)", x=.5, y=.97, just="top",gp=gpar(cex=1.5, fontface="bold"))
dev.off()

length(unique(miseq[Gene == "STAG2"]$`trial-pid`)) # 205 -- note that the muts.large file only contains 194...

# RAD21
features_RAD21 <- get_features_interpro(geneName_to_uniProt("RAD21"))
# Domains obtained from PMID: 32687945
domains.RAD21 <- GRanges("chr1",IRanges(start=c(1,362,558),end=c(103,403,628)),
                  fill=MetBrewer::met.brewer(name = "Hokusai1",3),
                   height=c(0.04,0.05,0.04))
names(domains.RAD21) <- c("SMC3","STAG1/2","SMC1A")

pdf("/data/claudia_gebhard/cohesins/RAD21_Lolliplot_JitterLabel.pdf",width=20,height=3.3)
make_lolliplot("RAD21","NM_006265",miseq,
               domains.RAD21,len_protein=features_RAD21$protein_length[1],
               label_on_feature=F,jitter='label')
grid.text("RAD21 (NM_006265)", x=.5, y=.9, just="top",gp=gpar(cex=1.5, fontface="bold"))
dev.off()

length(unique(miseq[Gene == "RAD21"]$`trial-pid`)) # 112 

write.xlsx(list(miseq[Gene == "STAG2"],miseq[Gene == "RAD21"]),"/data/claudia_gebhard/cohesins/lolliplot_datasets.xlsx")
# https://jianhong.github.io/trackViewer/articles/lollipopPlot.html

#-------------------------------------------------------------------------------
# Female/male bias: large cohort
#-------------------------------------------------------------------------------

miseq.sex <- merge(miseq[,.(trialpid=`trial-pid`,Gene,VAF=vaf,label=samplenr)],unique(pooled[,.(trialpid=`trialnr-pid`,sex)]),by="trialpid",all.x=T)
miseq.sex <- miseq.sex[!is.na(sex)]
miseq.sex$Gene <- gsub("[:;].*","",miseq.sex$Gene)

p.sex.large <- ggplot(miseq.sex[Gene %in% c("STAG2","RAD21")],aes(x=sex,y=VAF)) + 
  geom_violin(aes(color=sex)) + stat_summary(aes(color=sex),fun.data = "mean_sdl",  
                                                   fun.args = list(mult = 1),  geom = "pointrange") +
  geom_jitter(size=0.5) + facet_wrap(~Gene) + theme_cowplot() +
  ggsci::scale_color_npg(guide=F) + labs(x="") + 
  geom_signif(comparisons = list(c("male","female")),test='t.test')

ggsave(p.sex.large,filename="/data/claudia_gebhard/cohesins/analyses/VAFbias_plot_large.pdf",width=5,height=5)

casted.large <- dcast(miseq.sex,label+sex~Gene,fun.aggregate=function(x)ifelse(length(x) > 0,"mut","wt"))
sig.bias.large <- test.sex(casted.large)
unique(miseq.sex[!is.na(sex)][!label %in% miseq.sex[Gene=="STAG2"]$label][,.(sex,label)])[,.N,by=sex] # Unmutated cases

# Alternatives would include table() or dt[,,by=]. The latter does not count empty observations.
# To make it work in dplyr we need factors in the 'group_by'
# https://stackoverflow.com/questions/25956178/proper-idiom-for-adding-zero-count-rows-in-tidyr-dplyr
sex.counts.miseq <- miseq.sex %>% mutate(sex = as.factor(sex)) %>% group_by(Gene,sex,.drop=F) %>% 
  summarise(cnt = n(),.groups="drop") %>% group_by(Gene) %>% mutate(perc = cnt*100/sum(cnt)) %>%
  merge(sig.bias.large,by='Gene',all.x=T) %>% arrange(sex,perc)

p.sexbias.large <- sex.counts.miseq %>% mutate(Gene = factor(Gene,levels=unique(sex.counts.miseq$Gene))) %>%
  ggplot(aes(x=Gene,y=perc,fill=sex)) + geom_bar(position='stack',stat='identity') +
  theme_cowplot() + ggsci::scale_fill_npg(name="sex") +
  geom_text(aes(y=101,label=ifelse(pval < 0.05,"*",""),color=ifelse(oddsRatio > 1,"female","male")),size=5,show.legend=F) +
  labs(x="",y="Percentage of cases") + theme(axis.text.x = element_text(angle = 45,hjust=1,size=8)) + 
  geom_hline(yintercept = 50,linetype=2)

ggsave(p.sexbias.large,filename="/data/claudia_gebhard/cohesins/analyses/Sexbias_plot_large.pdf",width=8,height=5)