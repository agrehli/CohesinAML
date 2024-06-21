#!/usr/bin/env Rscript
# 7B_ggplot.R

### Plot Aggregate Peak Analysis (APA) results from coolpup.py over merged Hi-C Patient matrices
### Hi-C signal was aggregated over loop anchors overlapping 'weakened' or 'strengthened' loops
### between patients with STAG2 mutations compared to control patients
### Besides being differential loops, their anchors or 'edges' have to overlap Cohesin-associated enhancers
### loops are not divided by size, but analyses taking this into account were also performed with similar results
### The APA analyses were performed using coolpup.py in 7B_coolpuppy.sh

# R version: 4.2.2
# colorscale: logarithmic, centred on 1, log Fold Change range [0.25-4]


# 1. Re-plot APAs over differential loops using merged CD34 matrices

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
		date <- '23_03_22' #'23_01_04' #'23_01_03'
		indir <- paste0(wdir,'/APA/DifferentialAnchorsEnhancer/')
		printdir <- paste0(indir,date,'.')

	# Read all APA plots together
		#find matching coolpup.py output
		paths <- list.files(path = indir, recursive = TRUE,
				pattern = "\\.bedpe$", full.names = TRUE)
		length(paths) #[1] 3008

	# Start by plotting merged CD34 samples together:
		paths <- paths[grepl('merged_CD34',paths)] #only merged CD34 samples
		length(paths) #[1] 256
		# paths <- paths[!grepl('Norm.expected',paths)] #only merged CD34 samples
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
		dh(APAs.df) #430336
		APAs.df$Var2 <- factor(APAs.df$Var2,levels=rev(unique(APAs.df$Var2)),ordered = TRUE)
		APAs.df$value2 <- format(round(APAs.df$value, digits=2), nsmall = 2)
		table(APAs.df$loops)
		table(APAs.df$sample)

	# Select loop size groups  to plot
		SizeGroups <- c('<50kb','50kb-100kb','100kb-500kb','500kb-1Mb','1Mb-2Mb','2Mb-5Mb','>5Mb','All') %>%
			setNames(c('1','2','3','4','5','6','7','All'))
		SizeGroups <- SizeGroups[c(3:5,8)] #subset for only relevant groups since the rest have less than 10 loops!
		APAs.b.df <- APAs.df %>%
			mutate(loopsToSplit=loops) %>%
			separate(loopsToSplit, c('nothing','over','date','comparison','loopsAnchors','category','type','bed','source_loops','comparisonb','loops','AssociatedAnchorIs','edge','SG','SizeGroup','bed2'), sep="\\.") %>%
			select(-c('nothing','over','loopsAnchors','AssociatedAnchorIs','bed','source_loops','loops','SG','bed2')) %>%
			mutate(category_b=category) %>%
			mutate(category=ifelse(category=='weakened','Weakened','Strengthened')) %>%
			mutate(loops3= paste0(sample,comparison,category,type,edge,SizeGroup)) %>%
			mutate(category=factor(category,levels=c('Strengthened','Weakened'),ordered=TRUE)) %>%
			filter(SizeGroup %in% names(SizeGroups)) %>%
			mutate(SizeGroup=factor(SizeGroup, levels=names(SizeGroups),ordered=T)) %>%
			mutate(SizeGroup=recode(SizeGroup, !!!SizeGroups))
		table(APAs.b.df$sample)
		table(APAs.b.df$comparison)
		table(APAs.b.df$category)
		table(APAs.b.df$category_b)
		table(APAs.b.df$type)
		table(APAs.b.df$edge)
		table(APAs.b.df$SizeGroup)

	# Define names: CD34 model - merged matrices
		model <- c('siCtrl','stag1','stag2','rad21') %>%
				setNames(c('merged_CD34_siCtrl','merged_CD34_stag1','merged_CD34_stag2','merged_CD34_rad21'))
		conditions <- c('CTRL-HSPCs','STAG1 KD','STAG2 KD','RAD21 KD') %>%
			setNames(c('siCtrl','stag1','stag2','rad21'))
		AML_cohesin.c1 <- c('#df1413','#c58f07','#58812f','#623a87') %>%
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

	# [2023-03-22] Add number of loops:
		library("fs")
		LoopsPath <- indir
		NumberofLoops <- APAs.b.df %>%
			filter(is.finite(value)) %>%
			filter(Var1 %in% c('ID0')) %>%
			filter(Var2 %in% c('ID0')) %>%
			mutate(command = paste0('wc -l < ', LoopsPath, '22_11_22.',comparison,'.loopAnchors.',category_b,'.',type,'.bed.source_loops.SA2mutvsCTRL.loops.AssociatedAnchorIs.',edge,'.bed'),
				n=1)
			# 22_11_22.SA2mutvsCTRL.loopAnchors.strenghtened.CohesinAssEnhancersOverlap.bed.source_loops.SA2mutvsCTRL.loops.AssociatedAnchorIs.either.bed
			# 22_11_22.SA2mutvsCTRL.loopAnchors.weakened.CohesinAssEnhancersOverlap.bed.source_loops.SA2mutvsCTRL.loops.AssociatedAnchorIs.either.bed

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
		table(cd34.df$condition)
		table(cd34.df$sample)
		table(cd34.df$category)

	# colour palettes to try:
		palettes2 <- c('RdBu')

# 7B - Subpanel A - CD34 model - merged matrices
	
	# Select comparison/FDR combination
		cd34.b.df <- cd34.df %>%
				filter(comparison=='SA2mutvsCTRL',
						type=='CohesinAssEnhancersOverlap',
						edge=='either',
						SizeGroup=="All")
		print(table(cd34.b.df$sample))

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
			facet_grid(category~condition) +
			labs(title=paste0('Differential loops comparing SA2mutvsCTRL\npadj <= 0.05 & log2FoldChange >0.585')) +
			theme(	axis.title.x=element_blank(), axis.title.y=element_blank(), 
					axis.text.x=element_blank(), axis.ticks.x=element_blank(),
					axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.key.size = unit(1, 'cm')) +
			guides( fill = guide_colourbar(title.position='top', title.hjust = 0.5, label.position='bottom',nrow=1)) + 
			scale_fill_distiller('enrichment', trans = 'log10',
									palette = "RdBu", direction = -1, 
									limits=lims, breaks=breaks, labels=labels)	
		ggsave(paste0(printdir,'7B_CD34model.png'),
			q2, width = 24, height = 10, dpi=600, units='cm')
		ggsave(paste0(printdir,'7B_CD34model.pdf'),
			q2, width = 24, height = 10, dpi=600, units='cm', useDingbats=FALSE)

# 2. Re-plot APAs over differential loops using merged Cohesin patient matrices
	
	# Read all APA plots together again, but now focusing on merged patient matrices
		
		#find matching coolpup.py output
		paths <- list.files(path = indir, recursive = TRUE,
				pattern = "\\.bedpe$", full.names = TRUE)
		length(paths) #[1] 3008

	# Now merged Cohesin patient samples together
		paths <- paths[grepl('merged',paths)] #only merged Cohesin patient samples
		paths <- paths[!grepl('merged_CD34',paths)] #only merged Cohesin patient samples
		length(paths) #[1] 192

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
		dh(APAs.df) #322752
		APAs.df$Var2 <- factor(APAs.df$Var2,levels=rev(unique(APAs.df$Var2)),ordered = TRUE)
		APAs.df$value2 <- format(round(APAs.df$value, digits=2), nsmall = 2)
		table(APAs.df$loops)
		table(APAs.df$sample)

	# Select loop size groups  to plot
		SizeGroups <- c('<50kb','50kb-100kb','100kb-500kb','500kb-1Mb','1Mb-2Mb','2Mb-5Mb','>5Mb','All') %>%
			setNames(c('1','2','3','4','5','6','7','All'))
		SizeGroups <- SizeGroups[c(3:5,8)] #subset for only relevant groups since the rest have less than 10 loops!
		APAs.b.df <- APAs.df %>%
			mutate(loopsToSplit=loops) %>%
			separate(loopsToSplit, c('nothing','over','date','comparison','loopsAnchors','category','type','bed','source_loops','comparisonb','loops','AssociatedAnchorIs','edge','SG','SizeGroup','bed2'), sep="\\.") %>%
			select(-c('nothing','over','loopsAnchors','AssociatedAnchorIs','bed','source_loops','loops','SG','bed2')) %>%
			mutate(category_b=category) %>%
			mutate(category=ifelse(category=='weakened','Weakened','Strengthened')) %>%
			mutate(loops3= paste0(sample,comparison,category,type,edge,SizeGroup)) %>%
			mutate(category=factor(category,levels=c('Strengthened','Weakened'),ordered=TRUE)) %>%
			filter(SizeGroup %in% names(SizeGroups)) %>%
			mutate(SizeGroup=factor(SizeGroup, levels=names(SizeGroups),ordered=T)) %>%
			mutate(SizeGroup=recode(SizeGroup, !!!SizeGroups))
		table(APAs.b.df$sample)
		table(APAs.b.df$comparison)
		table(APAs.b.df$category)
		table(APAs.b.df$type)
		table(APAs.b.df$edge)
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

	# [2023-03-22] Add number of loops:
		library("fs")
		LoopsPath <- indir
		NumberofLoops <- APAs.b.df %>%
			filter(is.finite(value)) %>%
			filter(Var1 %in% c('ID0')) %>%
			filter(Var2 %in% c('ID0')) %>%
			mutate(command = paste0('wc -l < ', LoopsPath, '22_11_22.',comparison,'.loopAnchors.',category_b,'.',type,'.bed.source_loops.SA2mutvsCTRL.loops.AssociatedAnchorIs.',edge,'.bed'),
				n=1)
			# 22_11_22.SA2mutvsCTRL.loopAnchors.strenghtened.CohesinAssEnhancersOverlap.bed.source_loops.SA2mutvsCTRL.loops.AssociatedAnchorIs.either.bed
			# 22_11_22.SA2mutvsCTRL.loopAnchors.weakened.CohesinAssEnhancersOverlap.bed.source_loops.SA2mutvsCTRL.loops.AssociatedAnchorIs.either.bed

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
		table(cd34.df$condition)
		table(cd34.df$sample)
		table(cd34.df$category)
		# table(cd34.df$loops3)

	# colour palettes to try:
		palettes2 <- c('RdBu')

# 7B - Subpanel B - Cohesin mutation patient samples - merged matrices
	
	# Select comparison/FDR combination
		cd34.b.df <- cd34.df %>%
				filter(comparison=='SA2mutvsCTRL',
						type=='CohesinAssEnhancersOverlap',
						edge=='either',
						SizeGroup=="All")
		print(table(cd34.b.df$sample))

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
			facet_grid(category~condition) +
			labs(title=paste0('Differential loops comparing SA2mutvsCTRL\npadj <= 0.05 & log2FoldChange >0.585')) +
			theme(	axis.title.x=element_blank(), axis.title.y=element_blank(), 
					axis.text.x=element_blank(), axis.ticks.x=element_blank(),
					axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.key.size = unit(1, 'cm')) +
			guides( fill = guide_colourbar(title.position='top', title.hjust = 0.5, label.position='bottom',nrow=1)) + 
			scale_fill_distiller('enrichment', trans = 'log10',
									palette = "RdBu", direction = -1, 
									limits=lims, breaks=breaks, labels=labels)	
		ggsave(paste0(printdir,'7B_Patients.png'),
			q2, width = 24, height = 10, dpi=600, units='cm')
		ggsave(paste0(printdir,'7B_Patients.pdf'),
			q2, width = 24, height = 10, dpi=600, units='cm', useDingbats=FALSE)

