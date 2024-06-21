#!/usr/bin/env Rscript
# 7B_cohesin_e_p_diff_loops_size_split.R
# 22_12_05 L 190-333

### Load annotated table of loops and select differential loops (padj <= 0.05 & log2FoldChange >0.585)
### between patients with Cohesin (SA2) mutations compared to control samples
### Besides being differential loops, their anchors or 'edges' have to overlap Cohesin-associated enhancers
### then, split into size groups and save coordinates as an input for Aggregate Peak Analysis (APA)


# R version: 4.2.2
# Differential loops were defined using DEseq2: padj <= 0.05 & log2FoldChange >0.585 in HiC_AML.diffloops.DESeq2.rmd
# in addition, differential loops were classefied based on the overlap to H3K27ac and Cohesin in HiC_AML.diff.LoopAnchors.sh
# Here, Loops whose anchors had overlap to either of these enhancer-related features are subset and classified by size


# Load packages and plotting functions for accompanying MA-plot
	pacman::p_load("ggplot2","wesanderson","RColorBrewer","gridExtra","tidyverse","ggrepel","rcartocolor")
	# Custom ggplot theme for Fischer et al
	theme_custom <- function(psize=10, legend.position.p="bottom", aspect.ratio.p=1){
	  theme_bw(base_family='sans') +
	    theme(
	      text = element_text(size = psize, colour='#333333'),
	      plot.margin = unit(c(2,2,2,2), "pt"),
	      aspect.ratio = aspect.ratio.p,

	      panel.grid = element_blank(),
	      panel.grid.major = element_line(colour="grey", size = (0.5)),
	      panel.background = element_rect(fill = "white"), 
	      panel.border = element_rect(fill=NA, colour = "#333333", size=1),
	      panel.spacing = unit(0.5, "lines"),

	      plot.title=element_text(size=psize, colour='#333333'), 
	      axis.title.y=element_text(size=psize, colour='#333333'),
	      axis.title.x=element_text(size=psize, colour='#333333'), 
	      axis.text.x=element_text(size=psize, margin=margin(1,2,2,2,"pt"), colour='#333333'),
	      axis.text.y=element_text(size=psize, margin=margin(2,1,2,2,"pt"), colour='#333333'),
	      
	      legend.text=element_text(size=psize, colour='#333333'),
	      legend.title=element_text(size=psize, colour='#333333'),
	      legend.key.size = unit(0.5, 'lines'),
	      legend.margin=margin(0,0,0,0),
	      legend.box.margin=margin(0,0,0,0, "pt"),
	      legend.position=legend.position.p, 
	      strip.text = element_text(size=psize, colour='#333333'),
	      strip.text.x = element_text(size=psize, colour='#333333'),
	      strip.background = element_blank()
	    )
	}	

# Set input and output paths
	wdir <- '/private/'
	setwd(wdir)
	date <- '22_12_05'
	indir2 <- paste0(wdir, '/LoopCoordinates/DifferentialAnchorsEnhancer/')
	printdir <- paste0(indir2,date,'.')
	outdir2 <- paste0(indir2,date,'.')

	# Anchors:
		samplesAnchors <- c('SA2mutvsCTRL.loopAnchors.strenghtened.CohesinAssEnhancersOverlap.bed','SA2mutvsCTRL.loopAnchors.strenghtened.H3K27acOverlap.bed','SA2mutvsCTRL.loopAnchors.weakened.CohesinAssEnhancersOverlap.bed','SA2mutvsCTRL.loopAnchors.weakened.H3K27acOverlap.bed')
	# Loops
		indir <- paste0(wdir, '/LoopCoordinates/DifferentialLoopAnalysis/')
		samples <- 'SA2mutvsCTRL'
		end <- '.loops.txt'


	# For this analysis, we focused on differential loops defined using DEseq2 with significance thresholds (padj < 0.05) between ctrl Patients and Patients with STAG2 mutations

	# Go through each set of loop anchors overlapping Cohesin-associated loops and split by
	# cohesin anchor presence (both,minus,plus,either), and distance (0.1-1, 1-10, 10-25, 25-100 Mb)
		for (i in samplesAnchors){
			print(i)
			# load anchors
			diffAnchors <- as.data.frame(read.table(paste0(indir2,i), header=F, sep=''))
			# print(head(diffAnchors))
			
			# retrieve original loop IDs
			diffAnchorsIDsPlusEdge <- unique(unlist(lapply(diffAnchors$V4, strsplit,',')))
			diffAnchorsIDs <- as.data.frame(diffAnchorsIDsPlusEdge) %>%
				separate(diffAnchorsIDsPlusEdge, c('loopID','edge'), sep="\\_")

			# separate into both, minus and right edge
			bothAnchors <- diffAnchorsIDs$loopID[duplicated(diffAnchorsIDs$loopID)]
			minusAnchors <- diffAnchorsIDs[!duplicated(diffAnchorsIDs$loopID),] %>% filter(edge=='minus')
			plusAnchors <- diffAnchorsIDs[!duplicated(diffAnchorsIDs$loopID),] %>% filter(edge=='plus')
			eitherAnchors <- unique(c(bothAnchors,minusAnchors$loopID,plusAnchors$loopID))
			print(length(bothAnchors)) #[1] 161
			print(nrow(minusAnchors)) #[1] 151
			print(nrow(plusAnchors)) #[1] 302
			print(length(eitherAnchors))

			# load source loops and annotate SizeGroups
			diffLoops <- as.data.frame(read.table(paste0(indir,samples,end), header=T, sep=''))
			SizeGroups <- c('<50kb','50kb-100kb','100kb-500kb',
				'500kb-1Mb',
				'1Mb-2Mb','2Mb-5Mb','>5Mb')
			diffLoops2 <- diffLoops %>% 
							mutate(extstart1=start1-15000,extstart2=start2-15000,
									extend1=end1+15000,extend2=end2+15000) %>%
								mutate(Size=end2-start1) %>%
								mutate(SizeGroup=ifelse(Size < 50000, '<50kb',
													ifelse(Size < 100000, '50kb-100kb',
														ifelse(Size < 500000, '100kb-500kb',
															ifelse(Size < 1000000, '500kb-1Mb',
																ifelse(Size < 2000000, '1Mb-2Mb',
																	ifelse(Size < 5000000, '2Mb-5Mb','>5Mb'))))))) %>%
								mutate(SizeGroup=factor(SizeGroup, levels=SizeGroups, ordered=TRUE))
			print(with(diffLoops2 %>% filter(Row.names%in%eitherAnchors),table(SizeGroup)))

			# subset loops and save coordinates as bedpe files
			subcols <- c('chr1','start1','end1','chr2','start2','end2')
			write.table(diffLoops2 %>% filter(Row.names%in%bothAnchors) %>% select(subcols),
				paste0(outdir2,i,'.source_loops.',samples,'.loops.AssociatedAnchorIs.both.SizeGroup.All.bedpe'),
				quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
			write.table(diffLoops2 %>% filter(Row.names%in%minusAnchors$loopID) %>% select(subcols),
				paste0(outdir2,i,'.source_loops.',samples,'.loops.AssociatedAnchorIs.minus.SizeGroup.All.bedpe'),
				quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
			write.table(diffLoops2 %>% filter(Row.names%in%plusAnchors$loopID) %>% select(subcols),
				paste0(outdir2,i,'.source_loops.',samples,'.loops.AssociatedAnchorIs.plus.SizeGroup.All.bedpe'),
				quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
			write.table(diffLoops2 %>% filter(Row.names%in%eitherAnchors) %>% select(subcols),
				paste0(outdir2,i,'.source_loops.',samples,'.loops.AssociatedAnchorIs.either.SizeGroup.All.bedpe'),
				quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")


			# save as extended bedpe files for visualisation:
			subcols2 <- c('chr1','extstart1','extend1','chr2','extstart2','extend2')
			write.table(diffLoops2 %>% filter(Row.names%in%bothAnchors) %>% select(subcols2),
				paste0(outdir2,i,'.source_loops.',samples,'.loops.AssociatedAnchorIs.both.SizeGroup.All.extended.bedpe'),
				quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
			write.table(diffLoops2 %>% filter(Row.names%in%minusAnchors$loopID) %>% select(subcols2),
				paste0(outdir2,i,'.source_loops.',samples,'.loops.AssociatedAnchorIs.minus.SizeGroup.All.extended.bedpe'),
				quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
			write.table(diffLoops2 %>% filter(Row.names%in%plusAnchors$loopID) %>% select(subcols2),
				paste0(outdir2,i,'.source_loops.',samples,'.loops.AssociatedAnchorIs.plus.SizeGroup.All.extended.bedpe'),
				quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
			write.table(diffLoops2 %>% filter(Row.names%in%eitherAnchors) %>% select(subcols2),
				paste0(outdir2,i,'.source_loops.',samples,'.loops.AssociatedAnchorIs.either.SizeGroup.All.extended.bedpe'),
				quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

			#Print individual size groups separately
				for (j in seq_along(SizeGroups)){
					# Save as bedpe files
					subcols <- c('chr1','start1','end1','chr2','start2','end2')
					write.table(diffLoops2 %>% filter(SizeGroup==SizeGroups[j], Row.names%in%bothAnchors) %>% select(subcols),
						paste0(outdir2,i,'.source_loops.',samples,'.loops.AssociatedAnchorIs.both.SizeGroup.',j,'.bedpe'),
						quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
					write.table(diffLoops2 %>% filter(SizeGroup==SizeGroups[j], Row.names%in%minusAnchors$loopID) %>% select(subcols),
						paste0(outdir2,i,'.source_loops.',samples,'.loops.AssociatedAnchorIs.minus.SizeGroup.',j,'.bedpe'),
						quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
					write.table(diffLoops2 %>% filter(SizeGroup==SizeGroups[j], Row.names%in%plusAnchors$loopID) %>% select(subcols),
						paste0(outdir2,i,'.source_loops.',samples,'.loops.AssociatedAnchorIs.plus.SizeGroup.',j,'.bedpe'),
						quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
					write.table(diffLoops2 %>% filter(SizeGroup==SizeGroups[j], Row.names%in%eitherAnchors) %>% select(subcols),
						paste0(outdir2,i,'.source_loops.',samples,'.loops.AssociatedAnchorIs.either.SizeGroup.',j,'.bedpe'),
						quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

					# save as extended bedpe files for visualisation:
					subcols2 <- c('chr1','extstart1','extend1','chr2','extstart2','extend2')
					write.table(diffLoops2 %>% filter(SizeGroup==SizeGroups[j], Row.names%in%bothAnchors) %>% select(subcols2),
						paste0(outdir2,i,'.source_loops.',samples,'.loops.AssociatedAnchorIs.both.SizeGroup.',j,'.extended.bedpe'),
						quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
					write.table(diffLoops2 %>% filter(SizeGroup==SizeGroups[j], Row.names%in%minusAnchors$loopID) %>% select(subcols2),
						paste0(outdir2,i,'.source_loops.',samples,'.loops.AssociatedAnchorIs.minus.SizeGroup.',j,'.extended.bedpe'),
						quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
					write.table(diffLoops2 %>% filter(SizeGroup==SizeGroups[j], Row.names%in%plusAnchors$loopID) %>% select(subcols2),
						paste0(outdir2,i,'.source_loops.',samples,'.loops.AssociatedAnchorIs.plus.SizeGroup.',j,'.extended.bedpe'),
						quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
					write.table(diffLoops2 %>% filter(SizeGroup==SizeGroups[j], Row.names%in%eitherAnchors) %>% select(subcols2),
						paste0(outdir2,i,'.source_loops.',samples,'.loops.AssociatedAnchorIs.either.SizeGroup.',j,'.extended.bedpe'),
						quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
				}
		}