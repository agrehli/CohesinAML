#!/usr/bin/env Rscript
# 6B_6Ga_ggplot.R

### Plot Aggregate Peak Analysis (APA) results from coolpup.py over merged Hi-C Patient matrices
### Hi-C signal was aggregated over loop anchors overlapping 'weakened' or 'strengthened' loops 
### between patients with STAG2 mutations compared to control patients

### loops are divided by size and APAs performed independently for each class for 6B
### for 6G all loops are kept together

# R version: 4.2.2
# colorscale: logarithmic, centred on 1, log Fold Change range [0.25-4]
# number of loops considered: Weakened: 1519 Strengthened: 808
# Strengthened: 0,0,131,311,362,4,0 in size groups '<50kb','50kb-100kb','100kb-500kb','500kb-1Mb','1Mb-2Mb','2Mb-5Mb'
# Weakened: 5,169,900,220,222,3,0 in size groups '<50kb','50kb-100kb','100kb-500kb','500kb-1Mb','1Mb-2Mb','2Mb-5Mb'
# only size groups with >10 loops for both Strengthened and Weakened loops were plotted ('100kb-500kb','500kb-1Mb','1Mb-2Mb')


# Load required packages and functions
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
	date <- '2038_01_19'#'23_03_22'
	indir <- paste0(wdir,'/APA/DifferentialLoopAnalysis/')
	printdir <- paste0(indir,date,'.')

# Read all APA plots together
	#find matching coolpup.py output
	paths <- list.files(path = indir, recursive = TRUE,
			pattern = "bedpe\\.txt$", full.names = TRUE)

	# Start by plotting merged Cohesin patient samples together:
	paths <- paths[grepl('merged',paths)] #only merged samples
	paths <- paths[!grepl('merged_CD34',paths)] #only merged Cohesin patient samples

	# load coolpup.py results as matrices and flatten them out.
	APAs.l <- lapply(paths, function(x) {
			out <- as.data.frame(read.table(x, header=F, sep=''))
			colnames(out) <- rownames(out) <- paste0('ID',seq(-200,200,length.out=41))
			out <- reshape2::melt(as.matrix(out))
			out <- out %>%
				mutate(file = basename(x)) %>%
				separate(file, c('matrix','loops'), sep="\\.10kb") %>%
				separate(matrix, c('AML','sample'), sep="\\.")
			return(out)
			})
	dh(APAs.l[[1]]) #[1] 1681    6
	APAs.df <- do.call(rbind,APAs.l)
	dh(APAs.df) #680805
	APAs.df$Var2 <- factor(APAs.df$Var2,levels=rev(unique(APAs.df$Var2)),ordered = TRUE)
	APAs.df$value2 <- format(round(APAs.df$value, digits=2), nsmall = 2)
	table(APAs.df$loops)
	table(APAs.df$sample)

# Select loop size groups  to plot
	SizeGroups <- c('<50kb','50kb-100kb','100kb-500kb','500kb-1Mb','1Mb-2Mb','2Mb-5Mb','>5Mb','All') %>%
		setNames(c('1','2','3','4','5','6','7','All'))
	# subset for only relevant groups since the rest have less than 10 loops!
	# these are: '100kb-500kb','500kb-1Mb','1Mb-2Mb'
	# removed the 'All' category for the final plots, see below
	SizeGroups <- SizeGroups[c(3:5,8)] 
	APAs.b.df <- APAs.df %>%
		mutate(loopsToSplit=loops) %>%
		separate(loopsToSplit, c('nothing','over','date','comparison','loops','category','FDR','Zero','padj','SG','SizeGroup_b','bed','text'), sep="\\.") %>%
		select(-c('nothing','over','loops','FDR','Zero','SG','bed','text')) %>%
		mutate(loops3= paste0(sample,comparison,category,SizeGroup_b,padj)) %>%
		mutate(category=factor(category,levels=c('Strengthened','Weakened'),ordered=TRUE)) %>%
		filter(SizeGroup_b %in% names(SizeGroups)) %>%
		mutate(SizeGroup=factor(SizeGroup_b, levels=names(SizeGroups),ordered=T)) %>%
		mutate(SizeGroup=recode(SizeGroup, !!!SizeGroups))
	table(APAs.b.df$category)
	table(APAs.b.df$SizeGroup)

# Define names: Cohesin project - Cohesin patients - merged matrices
	model <- c('ctrl','stag2mut','rad21mu') %>% 
		setNames(c('merged_aml_ctrl','merged_sa2_mut','merged_rad21_mut'))
	conditions <- c('CTRL-AML','STAG2 mut','RAD21 mut') %>%
		setNames(c('ctrl','stag2mut','rad21mu'))
	AML_cohesin.c1 <- c('#df1413','#58812f','#623a87') %>%
		setNames(conditions)

# Get loop score:
# trimmed mean of the enrichment scores at the center 10 kb of the aggregated matrix is displayed in the bottom left corner as a quantification of enrichment.
	MeanLoopstrength <- APAs.b.df %>%
		filter(is.finite(value)) %>%
		filter(Var1 %in% c('ID0')) %>%
		filter(Var2 %in% c('ID0')) %>%
		group_by(loops3) %>%
		summarise(mean=round(mean(value,na.rm=TRUE),digits=2),
				median=median(value,na.rm=TRUE))
	APAs.df2 <- APAs.b.df %>%
		left_join(MeanLoopstrength,by=c('loops3'='loops3'))

# [2023-03-16]	Add number of loops:
	LoopsPath <- paste0(wdir,'/LoopCoordinates/DifferentialLoopAnalysis/')
	library("fs")
	NumberofLoops <- APAs.b.df %>%
		filter(is.finite(value)) %>%
		filter(Var1 %in% c('ID0')) %>%
		filter(Var2 %in% c('ID0')) %>%
		mutate(command = paste0('wc -l < ', LoopsPath, '22_12_05.',comparison,'.loops.',category,'.FDR.0.',padj,'.SizeGroup.',SizeGroup_b,'.bedpe'),
			n=1)
	NumberofLoops$n <- sapply(NumberofLoops$command, function(x) system(x,intern = TRUE))

	APAs.df3 <- APAs.df2  %>%
		left_join(NumberofLoops %>% select(loops3,n),by=c('loops3'='loops3'))
	
# Rename variables for consistency
	cd34.df <- APAs.df3 %>%
		filter(sample %in% names(model)) %>%
		mutate(sample=factor(sample, levels=names(model),ordered=T)) %>%
		mutate(condition=factor(sample, levels=names(model),ordered=T)) %>%
		mutate(condition=recode(condition, !!!model)) %>%
		mutate(condition=factor(condition, levels=names(conditions),ordered=T)) %>%
		mutate(condition=recode(condition, !!!conditions))  %>%
		arrange(desc(sample))
	head(cd34.df)
	table(cd34.df$condition)
	table(cd34.df$sample)
	table(cd34.df$category)
	# table(cd34.df$loops3)

# colour palettes to try:
	palettes2 <- c('RdBu')

# 6B - Subpanel A - Ctrl patient samples - split by distance
	cd34.b.df <- cd34.df %>%
				filter(sample=='merged_aml_ctrl',
					comparison=='SA2mutvsCTRL',
					padj=='05',
					SizeGroup!='All')
	# plot parameters
	breaks <- c(0.25,1,4) # add labels at symmetric distance from 1 in logarithmic scale
	lims <- breaks[c(1,3)] # limits for plotting range
	labels <- as.character(breaks)

	# ensure all pixels pixels are within colour plotting range. 
	# mean value (text) was calculated before this
	cd34.b.df$value[cd34.b.df$value>lims[2]] <- lims[2]
	cd34.b.df$value[cd34.b.df$value<lims[1]] <- lims[1]

	#make plot skeleton
	q2 <- qplot(data=cd34.b.df, x=Var1, y=Var2,
			geom='blank') +
		geom_tile(size=1,aes(fill=value)) +
		geom_text(aes(label= round(mean,2)), colour='#333333', size=3,
			vjust = "inward", hjust = "inward",
			data=subset(cd34.b.df,Var1=='ID180' & Var2=='ID-180')) +
		geom_text(aes(label= paste0('n=',n)), colour='#333333', size=3,
			vjust = "inward", hjust = "inward",
			data=subset(cd34.b.df,Var1=='ID180' & Var2=='ID180')) +
		theme_custom() +
		facet_grid(category~SizeGroup) +
		labs(title=paste0('Differential loops comparing SA2mutvsCTRL\npadj <= 0.05 & log2FoldChange >0.585')) +
		theme(	axis.title.x=element_blank(), axis.title.y=element_blank(), 
				axis.text.x=element_blank(), axis.ticks.x=element_blank(),
				axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.key.size = unit(1, 'cm')) +
		guides( fill = guide_colourbar(title.position='top', title.hjust = 0.5, label.position='bottom',nrow=1)) + 
		scale_fill_distiller('enrichment', trans = 'log10',
								palette = "RdBu", direction = -1, 
								limits=lims, breaks=breaks, labels=labels)

		ggsave(paste0(printdir,'6B_Ctrl.png'),
			q2, width = 24, height = 10, dpi=600, units='cm')
		ggsave(paste0(printdir,'6B_Ctrl.pdf'),
			q2, width = 24, height = 10, dpi=600, units='cm', useDingbats=FALSE)

# 6B - Subpanel B - STAG2 patient samples - split by distance
	cd34.b.df <- cd34.df %>%
				filter(sample=='merged_sa2_mut',
					comparison=='SA2mutvsCTRL',
					padj=='05',
					SizeGroup!='All')
	# plot
	breaks <- c(0.25,1,4)
	lims <- breaks[c(1,3)]
	labels <- as.character(breaks)

	# ensure all pixels pixels are within colour plotting range. 
	# mean value (text) was calculated before this
	cd34.b.df$value[cd34.b.df$value>lims[2]] <- lims[2]
	cd34.b.df$value[cd34.b.df$value<lims[1]] <- lims[1]

	#make plot skeleton
	q2 <- qplot(data=cd34.b.df, x=Var1, y=Var2,
			geom='blank') +
		geom_tile(size=1,aes(fill=value)) +
		geom_text(aes(label= round(mean,2)), colour='#333333', size=3,
			vjust = "inward", hjust = "inward",
			data=subset(cd34.b.df,Var1=='ID180' & Var2=='ID-180')) +
		geom_text(aes(label= paste0('n=',n)), colour='#333333', size=3,
			vjust = "inward", hjust = "inward",
			data=subset(cd34.b.df,Var1=='ID180' & Var2=='ID180')) +
		theme_custom() +
		facet_grid(category~SizeGroup) +
		labs(title=paste0('Differential loops comparing SA2mutvsCTRL\npadj <= 0.05 & log2FoldChange >0.585')) +
		theme(	axis.title.x=element_blank(), axis.title.y=element_blank(), 
				axis.text.x=element_blank(), axis.ticks.x=element_blank(),
				axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.key.size = unit(1, 'cm')) +
		guides( fill = guide_colourbar(title.position='top', title.hjust = 0.5, label.position='bottom',nrow=1)) + 
		scale_fill_distiller('enrichment', trans = 'log10',
								palette = "RdBu", direction = -1, 
								limits=lims, breaks=breaks, labels=labels)

		ggsave(paste0(printdir,'6B_STAG2mut.png'),
			q2, width = 24, height = 10, dpi=600, units='cm')
		ggsave(paste0(printdir,'6B_STAG2mut.pdf'),
			q2, width = 24, height = 10, dpi=600, units='cm', useDingbats=FALSE)

# 6G - Subpanel A - Every patient sample - not split by distance
	cd34.b.df <- cd34.df %>%
				filter(
					comparison=='SA2mutvsCTRL',
					padj=='05',
					SizeGroup=='All')
	# plot
	breaks <- c(0.25,1,4)
	lims <- breaks[c(1,3)]
	labels <- as.character(breaks)

	# ensure all pixels pixels are within colour plotting range. 
	# mean value (text) was calculated before this
	cd34.b.df$value[cd34.b.df$value>lims[2]] <- lims[2]
	cd34.b.df$value[cd34.b.df$value<lims[1]] <- lims[1]

	#make plot skeleton
	q2 <- qplot(data=cd34.b.df, x=Var1, y=Var2,
			geom='blank') +
		geom_tile(size=1,aes(fill=value)) +
		geom_text(aes(label= round(mean,2)), colour='#333333', size=3,
			vjust = "inward", hjust = "inward",
			data=subset(cd34.b.df,Var1=='ID180' & Var2=='ID-180')) +
		geom_text(aes(label= paste0('n=',n)), colour='#333333', size=3,
			vjust = "inward", hjust = "inward",
			data=subset(cd34.b.df,Var1=='ID180' & Var2=='ID180')) +
		theme_custom() +
		facet_grid(category~condition) +
		labs(title=paste0('Differential loops comparing SA2mutvsCTRL\npadj <= 0.05 & log2FoldChange >0.585')) +
		theme(	axis.title.x=element_blank(), axis.title.y=element_blank(), 
				axis.text.x=element_blank(), axis.ticks.x=element_blank(),
				axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.key.size = unit(1, 'cm')) +
		guides( fill = guide_colourbar(title.position='top', title.hjust = 0.5, label.position='bottom',nrow=1)) + 
		scale_fill_distiller('enrichment', trans = 'log10',
								palette = "RdBu", direction = -1, 
								limits=lims, breaks=breaks, labels=labels)

		ggsave(paste0(printdir,'6G_Patients.png'),
			q2, width = 24, height = 10, dpi=600, units='cm')
		ggsave(paste0(printdir,'6G_Patients.pdf'),
			q2, width = 24, height = 10, dpi=600, units='cm', useDingbats=FALSE)