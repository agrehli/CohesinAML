#!/usr/bin/env Rscript
# 6B_diff_loops_size_split.R

### Load annotated table of loops and select differential loops (padj <= 0.05 & log2FoldChange >0.585)
### between patients/mutant CD34 cells with Cohesin knockdown/mutations compared to control samples
### then, split into size groups and save coordinates as an input for Aggregate Peak Analysis (APA)

# R version: 4.2.2
# Differential loops were defined using DEseq2: padj <= 0.05 & log2FoldChange >0.585 in HiC_AML.diffloops.DESeq2.rmd


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
	indir <- paste0(wdir,'/LoopCoordinates/DifferentialLoopAnalysis/')
	printdir <- paste0(indir,date,'.')
	outdir1 <- paste0(indir,date,'.')
	# outdir2 <- paste0(wdir, '/APA/DifferentialLoopAnalysis/',date,'.')

# Go through each set of differential loops called by Alex using DEseq2 and save them as bedpe files split by 
# significance thresholds (padj < 0.001, 0.01,0.05), change in direction (strengthened/weakened), and distance (0.1-1, 1-10, 10-25, 25-100 Mb)
	samples2 <- c('RAD21KDvsCTRL','SA1KDvsCTRL','SA2KDvsCTRL','RAD21mutvsCTRL','SA2mutvsCTRL')
	for (i in samples2){
		print(i)

		#load differential loops
		samples <- i#'SA2mutvsCTRL'
		end <- '.loops.txt'
		diffLoops <- as.data.frame(read.table(paste0(indir,samples,end), header=T, sep=''))
		# summary(diffLoops)

		#add extended coordinates for visualisation in HiGlass
		diffLoops <-  diffLoops %>% 
				mutate(extstart1=start1-15000,extstart2=start2-15000,
						extend1=end1+15000,extend2=end2+15000)

		# table(diffLoops$chr1,diffLoops$chr2)
		#save the complete set together as bedpe files
		subcols <- c('chr1','start1','end1','chr2','start2','end2')
		subcols2 <- c('chr1','extstart1','extend1','chr2','extstart2','extend2')
		write.table(diffLoops[,subcols],
					paste0(outdir1,samples,'.loops.all.bed.bedpe'),
					quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
		write.table(diffLoops[,subcols2],
					paste0(outdir1,samples,'.loops.all.extended.bedpe'),
					quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

		# Select differential loops using 3 distinct significance thresholds as well as loop sizes and make bedpe files with each
		# Since the results didn't change between thresholds, we used: padj<0.05 and log2FoldChange >0.585
		y <- lapply(c(0.001, 0.01, 0.05),
			function(threshold){
				# Differential loop strength plot coloured by FDR
				# threshold <- 0.05
				print(threshold)

				#add differential, change direction and distance information
									# (0.1-1, 1-10, 10-25, 25-100 Mb) #too optimistic
				# 0-50000 50000-100000 100000-500000 500000-1000000 1000000-2000000 2000000-5000000
				SizeGroups <- c('<50kb','50kb-100kb','100kb-500kb','500kb-1Mb','1Mb-2Mb','2Mb-5Mb','>5Mb')
				diffLoops2 <- diffLoops %>% 
								mutate(DifferentialLoop=( padj<threshold & abs(log2FoldChange)>0.585)) %>%
								arrange(DifferentialLoop) %>%
								mutate(Category=ifelse(DifferentialLoop & log2FoldChange>0, 'Strengthened',
													ifelse(DifferentialLoop & log2FoldChange<0, 'Weakened', 'NS'))) %>% 
								mutate(Category=factor(Category,levels=c('Strengthened','NS','Weakened'),ordered=TRUE)) %>%
								mutate(Size=end2-start1) %>%
								mutate(SizeGroup=ifelse(Size < 50000, '<50kb',
													ifelse(Size < 100000, '50kb-100kb',
														ifelse(Size < 500000, '100kb-500kb',
															ifelse(Size < 1000000, '500kb-1Mb',
																ifelse(Size < 2000000, '1Mb-2Mb',
																	ifelse(Size < 5000000, '2Mb-5Mb','>5Mb'))))))) %>%
								mutate(SizeGroup=factor(SizeGroup, levels=SizeGroups, ordered=TRUE))
				
				# How many differential loops are there?
				newLevels <- paste(c(levels(diffLoops2$Category)), c(table(diffLoops2$Category)), sep='\n') %>%
						setNames(c('Strengthened','NS','Weakened'))
				print(newLevels)
				diffLoops2 <- diffLoops2 %>% 
								mutate(Category=recode(Category, !!!newLevels))
				# What about split by distance
				print(with(diffLoops2,table(Category,SizeGroup)))

				# MA-plot
					DiffCols <- c("#B31B21", "darkgray", "#1465AC") %>% setNames(newLevels)
					maxL2FC <- max(abs(diffLoops2$log2FoldChange))
					q <- qplot(data=diffLoops2,
								x=baseMean,y=log2FoldChange,col=Category) +
						theme_custom(aspect.ratio=1,legend.position.p="top") +
						ylim((-1*maxL2FC),maxL2FC) +
						scale_colour_manual(values=DiffCols) +
						labs(title=paste0('Differential loops comparing ',samples,'\npadj <=', threshold, ' & log2FoldChange >0.585'),
							y=paste0('Loop strength log2FC(',samples,')'))
					ggsave(paste0(outdir1,samples,'.loops.FDR.',threshold,'.A.png'), 
						q, width = 12, height = 12, dpi=600, units='cm')
					q <- qplot(data=diffLoops2,
								x=baseMean,y=log2FoldChange,col=SizeGroup) +
						theme_custom(aspect.ratio=1,legend.position.p="top") +
						ylim((-1*maxL2FC),maxL2FC) +
						labs(title=paste0('Differential loops comparing ',samples,'\npadj <=', threshold, ' & log2FoldChange >0.585'),
							y=paste0('Loop strength log2FC(',samples,')'))
					ggsave(paste0(outdir1,samples,'.loops.FDR.',threshold,'.B.png'), 
						q, width = 12, height = 12, dpi=600, units='cm')

				#Print all sizes together

					# Save as bedpe files
					write.table(diffLoops2 %>% filter(DifferentialLoop,log2FoldChange>0) %>%  select(subcols),
						paste0(outdir1,samples,'.loops.Strengthened.FDR.',threshold,'.SizeGroup.All.bedpe'),
						quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
					write.table(diffLoops2 %>% filter(DifferentialLoop,log2FoldChange<0) %>%  select(subcols),
						paste0(outdir1,samples,'.loops.Weakened.FDR.',threshold,'.SizeGroup.All.bedpe'),
						quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
					# save as extended bedpe files for visualisation:
					write.table(diffLoops2 %>% filter(DifferentialLoop,log2FoldChange>0) %>%  select(subcols2),
						paste0(outdir1,samples,'.loops.Strengthened.FDR.',threshold,'.SizeGroup.All.extended.bedpe'),
						quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
					write.table(diffLoops2 %>% filter(DifferentialLoop,log2FoldChange<0) %>%  select(subcols2),
						paste0(outdir1,samples,'.loops.Weakened.FDR.',threshold,'.SizeGroup.All.extended.bedpe'),
						quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")


				#Print individual size groups separately
				for (j in seq_along(SizeGroups)){
					# Save as bedpe files
					write.table(diffLoops2 %>% filter(DifferentialLoop, SizeGroup==SizeGroups[j],log2FoldChange>0) %>%  select(subcols),
						paste0(outdir1,samples,'.loops.Strengthened.FDR.',threshold,'.SizeGroup.',j,'.bedpe'),
						quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
					write.table(diffLoops2 %>% filter(DifferentialLoop, SizeGroup==SizeGroups[j],log2FoldChange<0) %>%  select(subcols),
						paste0(outdir1,samples,'.loops.Weakened.FDR.',threshold,'.SizeGroup.',j,'.bedpe'),
						quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
					# save as extended bedpe files for visualisation:
					write.table(diffLoops2 %>% filter(DifferentialLoop, SizeGroup==SizeGroups[j],log2FoldChange>0) %>%  select(subcols2),
						paste0(outdir1,samples,'.loops.Strengthened.FDR.',threshold,'.SizeGroup.',j,'.extended.bedpe'),
						quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
					write.table(diffLoops2 %>% filter(DifferentialLoop, SizeGroup==SizeGroups[j],log2FoldChange<0) %>%  select(subcols2),
						paste0(outdir1,samples,'.loops.Weakened.FDR.',threshold,'.SizeGroup.',j,'.extended.bedpe'),
						quote = FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
				}
				
			})


	}