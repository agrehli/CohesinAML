#!/usr/bin/env Rscript
# S6E_S6F_ggplot.R

### Plot Aggregate Peak Analysis (APA) results from coolpup.py over indivdual patient or replicate Hi-C matrices of Ctrl and Cohesin deficient cells/patients.
### Hi-C signal was aggregated over loop anchors overlapping 'weakened' or 'strengthened' loops 
### between patients with STAG2 mutations compared to control patients
### loops are not divided by size, but analyses taking this into account were also performed with similar results
### The APA analyses were performed using coolpup.py in S6E_coolpuppy.sh

# R version: 4.2.2
# colorscale: logarithmic, centred on 1, log Fold Change range [0.25-4]


# 1. Re-plot APAs over differential loops using individual replicate CD34 matrices

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
		date <- '23_03_23' #'22_12_21'#'22_12_07'#'22_12_05'#'22_11_28'#'22_11_22'
		indir <- paste0(wdir,'/APA/DifferentialLoopAnalysis/')
		printdir <- paste0(indir,date,'.')

	# Read all APA plots together		
		#find matching coolpup.py output
		paths <- list.files(path = indir, recursive = TRUE,
			pattern = "bedpe\\.txt$", full.names = TRUE)

		# Now plotting individual replicate CD34 samples
		paths <- paths[grepl('\\.CD34_',paths)] #only individual replicate CD34 samples

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
		dh(APAs.df) #1933150 when together
		APAs.df$Var2 <- factor(APAs.df$Var2,levels=rev(unique(APAs.df$Var2)),ordered = TRUE)
		APAs.df$value2 <- format(round(APAs.df$value, digits=2), nsmall = 2)
		table(APAs.df$loops)
		table(APAs.df$sample)

	# Select loop size groups to plot
		SizeGroups <- c('<50kb','50kb-100kb','100kb-500kb','500kb-1Mb','1Mb-2Mb','2Mb-5Mb','>5Mb','All') %>%
			setNames(c('1','2','3','4','5','6','7','All'))
		SizeGroups <- SizeGroups[c(3:5,8)] #subset for only relevant groups since the rest have less than 10 loops!
		APAs.b.df <- APAs.df %>%
			mutate(loopsToSplit=loops) %>%
			separate(loopsToSplit, c('nothing','over','date','comparison','loops','category','FDR','Zero','padj','SG','SizeGroup','bed','text'), sep="\\.") %>%
			select(-c('nothing','over','loops','FDR','Zero','SG','bed','text')) %>%
			mutate(loops3= paste0(sample,comparison,category,SizeGroup,padj)) %>%
			mutate(category=factor(category,levels=c('Strengthened','Weakened'),ordered=TRUE)) %>%
			filter(SizeGroup %in% names(SizeGroups)) %>%
			mutate(SizeGroup=factor(SizeGroup, levels=names(SizeGroups),ordered=T)) %>%
			mutate(SizeGroup=recode(SizeGroup, !!!SizeGroups))
		table(APAs.b.df$category)
		table(APAs.b.df$SizeGroup)

	# Define names: # Cohesin project - CD34 model - replicates
		model <- c('14_siCtrl','17_siCtrl','18_siCtrl','20_siCtrl','21_siCtrl','22_siCtrl','27_siCtrl','28_siCtrl',
					'14_stag1','17_stag1','20_stag1','21_stag1','27_stag1','28_stag1',
					'14_stag2','17_stag2','20_stag2','21_stag2','22_stag2','28_stag2',
					'18_rad21','20_rad21','22_rad21','27_rad21','28_rad21') %>%
				setNames(c('CD34_14_3_siCtrl_R1','CD34_17_3_siCtrl_R1','CD34_18_4_siCtrl_R1','CD34_20_6_siCtrl_Rep1_R1','CD34_21_4_siCtrl_Rep1_R1','CD34_22_3_siCtrl_R1','CD34_27_4_siCtrl_R1','CD34_28_6_siCtrl_R1',
					'CD34_14_1_SA1_KD_R1','CD34_17_1_SA1_KD_R1','CD34_20_4_SA1_KD_R1','CD34_21_2_SA1_KD_R1','CD34_27_3_SA1_KD_R1','CD34_28_4_SA1_KD_R1',
					'CD34_14_2_SA2_KD_R1','CD34_17_2_SA2_KD_R1','CD34_20_5_SA2_KD_R1','CD34_21_3_SA2_KD_R1','CD34_22_2_SA2_KD_R1','CD34_28_5_SA2_KD_R1',
					'CD34_18_1_RAD21_KD_R1','CD34_20_1_RAD21_KD_R1','CD34_22_1_RAD21_KD_R1','CD34_27_1_RAD21_KD_R1','CD34_28_1_RAD21_KD_R1'))

		replicates <- c(paste0('Ctrl_',1:8),paste0('SA1_',1:6),paste0('SA2_',1:6),paste0('RAD21_',1:5)) %>%
				setNames(c('CD34_14_3_siCtrl_R1','CD34_17_3_siCtrl_R1','CD34_18_4_siCtrl_R1','CD34_20_6_siCtrl_Rep1_R1','CD34_21_4_siCtrl_Rep1_R1','CD34_22_3_siCtrl_R1','CD34_27_4_siCtrl_R1','CD34_28_6_siCtrl_R1',
					'CD34_14_1_SA1_KD_R1','CD34_17_1_SA1_KD_R1','CD34_20_4_SA1_KD_R1','CD34_21_2_SA1_KD_R1','CD34_27_3_SA1_KD_R1','CD34_28_4_SA1_KD_R1',
					'CD34_14_2_SA2_KD_R1','CD34_17_2_SA2_KD_R1','CD34_20_5_SA2_KD_R1','CD34_21_3_SA2_KD_R1','CD34_22_2_SA2_KD_R1','CD34_28_5_SA2_KD_R1',
					'CD34_18_1_RAD21_KD_R1','CD34_20_1_RAD21_KD_R1','CD34_22_1_RAD21_KD_R1','CD34_27_1_RAD21_KD_R1','CD34_28_1_RAD21_KD_R1'))

		conditions <- c('CTRL\nHSPCs','STAG1\nKD','STAG2\nKD','RAD21\nKD') %>%
				setNames(c('siCtrl','stag1','stag2','rad21'))
		
		AML_cohesin.c1 <- c('#B22222','#B8860B','#2E8B57','#C71585') %>%
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
	
	# Rename variables for consistency
		cd34.df <- APAs.df2 %>%
			filter(sample %in% names(model)) %>%
			mutate(sample=factor(sample, levels=names(model),ordered=T)) %>%
			mutate(condition=factor(sample, levels=names(model),ordered=T)) %>%
			mutate(condition=recode(condition, !!!model)) %>%
			separate(condition,c('replicate','condition')) %>%
			mutate(condition=factor(condition, levels=names(conditions),ordered=T)) %>%
			mutate(condition=recode(condition, !!!conditions)) %>%
			mutate(replicate2=factor(sample, levels=names(replicates),ordered=T)) %>%
			mutate(replicate2=recode(replicate2, !!!replicates)) %>%
			separate(replicate2, c('c3','replicate3'), sep="\\_", remove=FALSE) %>%
			mutate(replicate3=factor(replicate3, levels=1:8,ordered=T)) %>%
			arrange(desc(sample))
		table(cd34.df$replicate3)
		table(cd34.df$condition)
		table(cd34.df$sample)
		table(cd34.df$category)
		# table(cd34.df$loops3)

	# colour palettes to try:
		palettes2 <- c('RdBu')

# Data not shown - APA - individual CD34 replicates - not split by distance

	# # Re-plot 4x2 'Replicate vs Condition' once per comparison, FDR combination and change direction
		comparisons <- unique(cd34.df$comparison)
		fdrs <- unique(cd34.df$padj)
		categories <- unique(cd34.df$category)
		
		for (i in comparisons){
			for (j in fdrs){
				for (l in categories){

					#Select comparison/FDR combination
					print(paste0(i,'_',j,'_',l))
					cd34.b.df <- cd34.df %>%
						filter(comparison==i,padj==j,category==l,SizeGroup=="All")
					print(table(cd34.b.df$sample))
					print(table(cd34.b.df$condition))
					print(table(cd34.b.df$replicate3))

					# plot parameters
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
						geom_text(aes(label= round(mean,2)), colour='#333333', size=2,
							vjust = "inward", hjust = "inward",
							data=subset(cd34.b.df,Var1=='ID180' & Var2=='ID-180')) +
						theme_custom() +
						facet_grid(condition~replicate3) +
						labs(title=paste0(l, ' loops comparing ',i,'\npadj <= 0.', j, ' & log2FoldChange >0.585')) +
						theme(	axis.title.x=element_blank(), axis.title.y=element_blank(), 
								axis.text.x=element_blank(), axis.ticks.x=element_blank(),
								axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.key.size = unit(1, 'cm')) +
						guides( fill = guide_colourbar(title.position='top', title.hjust = 0.5, label.position='bottom',nrow=1))

					##fill distiller:
					for (k in c(palettes2)){
						q3 <- q2 + scale_fill_distiller('enrichment', trans = 'log10',
												palette = k, direction = -1, 
												limits=lims, breaks=breaks, labels=labels)	
						ggsave(paste0(printdir,'AML2021.APA.DifferentialLoops.',i,'.FDR.',j,'.samples.individual_CD34_matrices.SizeGroup.All.category.',l,'.Colour.',k,'.png'),
							q3, width = 24, height = 16, dpi=600, units='cm')
						ggsave(paste0(printdir,'AML2021.APA.DifferentialLoops.',i,'.FDR.',j,'.samples.individual_CD34_matrices.SizeGroup.All.category.',l,'.Colour.',k,'.pdf'),
							q3, width = 24, height = 16, dpi=600, units='cm', useDingbats=FALSE)
					}
				}
			}
		}

# S6E - APA - individual patients with cohesin mutations - not split by distance

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
		date <- '23_03_23' #'23_03_16'#'23_01_03'#'22_12_07'#'22_12_05'#'22_11_28'#'22_11_22'
		indir <- paste0(wdir,'/APA/DifferentialLoopAnalysis/')
		printdir <- paste0(indir,date,'.')

	# Read all APA plots together		
		#find matching coolpup.py output
		paths <- list.files(path = indir, recursive = TRUE,
			pattern = "bedpe\\.txt$", full.names = TRUE)

		# Now plotting individual patients with cohesin mutations
		paths <- paths[!grepl('\\.CD34_',paths)]
		paths <- paths[!grepl('merged_',paths)] #only individual patients with cohesin mutations

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
		dh(APAs.df) #1159890
		APAs.df$Var2 <- factor(APAs.df$Var2,levels=rev(unique(APAs.df$Var2)),ordered = TRUE)
		APAs.df$value2 <- format(round(APAs.df$value, digits=2), nsmall = 2)
		table(APAs.df$loops)
		table(APAs.df$sample)

	# Select loop size groups to plot
		SizeGroups <- c('<50kb','50kb-100kb','100kb-500kb','500kb-1Mb','1Mb-2Mb','2Mb-5Mb','>5Mb','All') %>%
			setNames(c('1','2','3','4','5','6','7','All'))
		SizeGroups <- SizeGroups[c(3:5,8)] #subset for only relevant groups since the rest have less than 10 loops!
		APAs.b.df <- APAs.df %>%
			mutate(loopsToSplit=loops) %>%
			separate(loopsToSplit, c('nothing','over','date','comparison','loops','category','FDR','Zero','padj','SG','SizeGroup','bed','text'), sep="\\.") %>%
			select(-c('nothing','over','loops','FDR','Zero','SG','bed','text')) %>%
			mutate(loops3= paste0(sample,comparison,category,SizeGroup,padj)) %>%
			mutate(category=factor(category,levels=c('Strengthened','Weakened'),ordered=TRUE)) %>%
			filter(SizeGroup %in% names(SizeGroups)) %>%
			mutate(SizeGroup=factor(SizeGroup, levels=names(SizeGroups),ordered=T)) %>%
			mutate(SizeGroup=recode(SizeGroup, !!!SizeGroups))
		table(APAs.b.df$category)
		table(APAs.b.df$SizeGroup)

	# Define names: # Cohesin project - CD34 model - replicates
			patients <- c( 'AMLctr_16911','AMLctr_18136','AMLctr_18519','AMLctr_19405','AMLctr_19416','AMLctr_21047','AMLctr_21290',
							'AMLRAD21_UKR186','AMLRAD21_23039','AMLRAD21_26830','AMLRAD21_38455',
							'AMLSA2_24743','AMLSA2_27396','AMLSA2_29728','AMLSA2_9708') %>%
					setNames(c('AML_ctr_16911_R1','AML_ctr_18136_Rep1_R1','AML_ctr_18519_R1','AML_ctr_19405_Rep1_R1','AML_ctr_19416_R1','AML_ctr_21047_Rep1_R1','AML_ctr_21290_R1',
					'RAD21_UKR186_Rep1_R1','AML_RAD21_23039_R1','AML_RAD21_26830_R1','AML_RAD21_38455_R1',
					'AML_SA2_24743_R1','AML_SA2_27396_Rep1_R1','AML_SA2_29728_R1','AML_SA2_9708_Rep2_S56_L002_R1'))

			patients2 <- c(paste0('AMLctr_',1:7),paste0('AMLRAD21_',1:4),paste0('AMLSA2_',1:4)) %>%
					setNames(c('AML_ctr_16911_R1','AML_ctr_18136_Rep1_R1','AML_ctr_18519_R1','AML_ctr_19405_Rep1_R1','AML_ctr_19416_R1','AML_ctr_21047_Rep1_R1','AML_ctr_21290_R1',
					'RAD21_UKR186_Rep1_R1','AML_RAD21_23039_R1','AML_RAD21_26830_R1','AML_RAD21_38455_R1',
					'AML_SA2_24743_R1','AML_SA2_27396_Rep1_R1','AML_SA2_29728_R1','AML_SA2_9708_Rep2_S56_L002_R1'))

			patients4 <- c('2236','1551','1747','3323','3488','6246','5285','UKR186','7314','12557','41580a','12514','12567','24603','2193') %>%
					setNames(c('AML_ctr_16911_R1','AML_ctr_18136_Rep1_R1','AML_ctr_18519_R1','AML_ctr_19405_Rep1_R1','AML_ctr_19416_R1','AML_ctr_21047_Rep1_R1','AML_ctr_21290_R1',
					'RAD21_UKR186_Rep1_R1','AML_RAD21_23039_R1','AML_RAD21_26830_R1','AML_RAD21_38455_R1',
					'AML_SA2_24743_R1','AML_SA2_27396_Rep1_R1','AML_SA2_29728_R1','AML_SA2_9708_Rep2_S56_L002_R1'))

			conditions <- c('CTRL\nAML','STAG2\nmut','RAD21\nmut') %>%
				setNames(c('AMLctr','AMLSA2','AMLRAD21'))

			AML_cohesin_patients.c1 <- c('#B22222','#2E8B57','#C71585') %>%
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
	
	# Rename variables for consistency
		patients.df <- APAs.df2 %>%
				filter(sample %in% names(patients)) %>%
				mutate(sample=factor(sample, levels=names(patients),ordered=T)) %>%
				mutate(patient=factor(sample, levels=names(patients),ordered=T)) %>%
				mutate(patient=recode(patient, !!!patients)) %>%
				separate(patient,c('condition','patient')) %>%
				mutate(condition=factor(condition, levels=names(conditions),ordered=T)) %>%
				mutate(condition=recode(condition, !!!conditions)) %>%
				mutate(patient2=factor(sample, levels=names(patients2),ordered=T)) %>%
				mutate(patient2=recode(patient2, !!!patients2)) %>%
				separate(patient2, c('c3','patient3'), sep="\\_", remove=FALSE) %>%
				mutate(patient3=factor(patient3, levels=1:7,ordered=T)) %>%
				mutate(patient4=factor(sample, levels=names(patients4),ordered=T)) %>%
				mutate(patient4=recode(patient4, !!!patients4)) %>%
				arrange(desc(sample))
		head(patients.df)
		table(patients.df$patient3)
		table(patients.df$patient4)
		table(patients.df$sample)
		table(patients.df$patient)
		table(patients.df$condition)
		table(patients.df$category)

	# colour palettes to try:
		palettes2 <- c('RdBu')

	# S6E - Subpanel A - Strengthened loops - Patients
		patients.b.df <- patients.df %>%
					filter(category=='Strengthened',
						comparison=='SA2mutvsCTRL',
						padj=='05',
						SizeGroup=='All')
		# plot
		breaks <- c(0.25,1,4)
		lims <- breaks[c(1,3)]
		labels <- as.character(breaks)

		# ensure all pixels pixels are within colour plotting range. 
		# mean value (text) was calculated before this
		patients.b.df$value[patients.b.df$value>lims[2]] <- lims[2]
		patients.b.df$value[patients.b.df$value<lims[1]] <- lims[1]

		#make plot skeleton
		q2 <- qplot(data=patients.b.df, x=Var1, y=Var2,
				geom='blank') +
			geom_tile(size=1,aes(fill=value)) +
			geom_text(aes(label= round(mean,2)), colour='#333333', size=3,
				vjust = "inward", hjust = "inward",
				data=subset(patients.b.df,Var1=='ID180' & Var2=='ID-180')) +
			geom_text(aes(label= patient4), colour='#333333', size=3,
				vjust = "inward", hjust = "inward",
				data=subset(patients.b.df,Var1=='ID180' & Var2=='ID180')) +
			theme_custom() +
			facet_grid(condition~patient3) +
			# facet_grid(category~SizeGroup) +
			labs(title=paste0('Strengthened loops comparing SA2mutvsCTRL\npadj <= 0.05 & log2FoldChange >0.585')) +
			theme(	axis.title.x=element_blank(), axis.title.y=element_blank(), 
					axis.text.x=element_blank(), axis.ticks.x=element_blank(),
					axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.key.size = unit(1, 'cm')) +
			guides( fill = guide_colourbar(title.position='top', title.hjust = 0.5, label.position='bottom',nrow=1)) + 
			scale_fill_distiller('enrichment', trans = 'log10',
									palette = "RdBu", direction = -1, 
									limits=lims, breaks=breaks, labels=labels)

			ggsave(paste0(printdir,'S6E_Strengthened.png'),
				q2, width = 24, height = 12, dpi=600, units='cm')
			ggsave(paste0(printdir,'S6E_Strengthened.pdf'),
				q2, width = 24, height = 12, dpi=600, units='cm', useDingbats=FALSE)

	# S6E - Subpanel B - Weakened loops - Patients
		patients.b.df <- patients.df %>%
					filter(category=='Weakened',
						comparison=='SA2mutvsCTRL',
						padj=='05',
						SizeGroup=='All')
		# plot
		breaks <- c(0.25,1,4)
		lims <- breaks[c(1,3)]
		labels <- as.character(breaks)

		# ensure all pixels pixels are within colour plotting range. 
		# mean value (text) was calculated before this
		patients.b.df$value[patients.b.df$value>lims[2]] <- lims[2]
		patients.b.df$value[patients.b.df$value<lims[1]] <- lims[1]

		#make plot skeleton
		q2 <- qplot(data=patients.b.df, x=Var1, y=Var2,
				geom='blank') +
			geom_tile(size=1,aes(fill=value)) +
			geom_text(aes(label= round(mean,2)), colour='#333333', size=3,
				vjust = "inward", hjust = "inward",
				data=subset(patients.b.df,Var1=='ID180' & Var2=='ID-180')) +
			geom_text(aes(label= patient4), colour='#333333', size=3,
				vjust = "inward", hjust = "inward",
				data=subset(patients.b.df,Var1=='ID180' & Var2=='ID180')) +
			theme_custom() +
			facet_grid(condition~patient3) +
			# facet_grid(category~SizeGroup) +
			labs(title=paste0('Weakened loops comparing SA2mutvsCTRL\npadj <= 0.05 & log2FoldChange >0.585')) +
			theme(	axis.title.x=element_blank(), axis.title.y=element_blank(), 
					axis.text.x=element_blank(), axis.ticks.x=element_blank(),
					axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.key.size = unit(1, 'cm')) +
			guides( fill = guide_colourbar(title.position='top', title.hjust = 0.5, label.position='bottom',nrow=1)) + 
			scale_fill_distiller('enrichment', trans = 'log10',
									palette = "RdBu", direction = -1, 
									limits=lims, breaks=breaks, labels=labels)

			ggsave(paste0(printdir,'S6E_Weakened.png'),
				q2, width = 24, height = 12, dpi=600, units='cm')
			ggsave(paste0(printdir,'S6E_Weakened.pdf'),
				q2, width = 24, height = 12, dpi=600, units='cm', useDingbats=FALSE)

# S6F - Barplot quantifying S6E - individual patients with cohesin mutations - not split by distance

	# [2023-03-23 19:16] Replot barplot with correct colours, font-, dot- and line- sizes

	# S6F - barplot quantifying S6E

		i = "SA2mutvsCTRL"
		print(i)

			patients.b.df <- patients.df %>%
				filter(comparison==i,padj=="05",Var1=='ID0',Var2=='ID0',SizeGroup=="All") %>%
				mutate(ID=paste0(condition, category))
			print(table(patients.b.df$sample))
			print(table(patients.b.df$condition))
			print(table(patients.b.df$category))

		# error bar:
			error_bar <- patients.b.df %>%
			group_by(ID) %>%
			summarise(mean_mean=mean(mean,na.rm=TRUE),
					  sd_mean=sd(mean,na.rm=TRUE),
					  ymin_mean=mean_mean-sd_mean,
					  ymax_mean=mean_mean+sd_mean)
			patients.b.df <- patients.b.df %>%
				left_join(error_bar,by='ID')

		# barplot
			set.seed(10) ; q2 <- qplot(data=patients.b.df, 
				x=condition, y=mean, 
				col=condition, geom='blank') +
			geom_bar(data=subset(patients.b.df,patient3==1),
				aes(y = mean_mean),
				col='black',fill='darkgrey',stat="identity", 
				position=position_dodge(), width = 0.75) +
			geom_point(data=subset(patients.b.df,patient3==1),
				aes(y = mean_mean),size=1,col='black') +
			geom_errorbar(data=subset(patients.b.df,patient3==1),
				aes(ymin = ymin_mean,
					ymax = ymax_mean),
				width = 0.15, size=1,col='black') +
			geom_jitter(size=1.5,width = 0.15) +
			facet_grid(rows = vars(category)) +
			theme_custom(psize=10) + 
			theme(axis.title.x=element_blank(),legend.position = "none") +
			labs(y=paste0('Loop strength')) +
			ylim(c(0,6)) +
			scale_colour_manual('',values = AML_cohesin_patients.c1) +
			scale_fill_manual('',values = AML_cohesin_patients.c1)

			set.seed(10)
			ggsave(paste0(printdir,'S6F.png'),
					q2, width = 5.5, height = 11, dpi=600, units='cm')
			set.seed(10)
			ggsave(paste0(printdir,'S6F.pdf'),
				q2, width = 5.5, height = 11, dpi=600, units='cm', useDingbats=FALSE)
			write.table(patients.b.df,
				 file = paste0(printdir,'S6F.txt'),
				 quote = FALSE,sep = "\t",
				 row.names = FALSE,col.names = TRUE)
			saveRDS(patients.b.df,
				 file = paste0(printdir,'S6F.rds'))

# Data not shown - Barplot - individual patients and CD34 KD replicates with cohesin mutations/depletion - not split by distance

		# I'll make a version in which I plot both together: Replicates from CD34 Cohesin KDs and Patients with Cohesin mutations

		# merge dataframes contating both pieces of information:
			patients.df$patient4 <- NULL
			patients.df$patient <- NULL
			patients.df$patient2 <- NULL
			patients.df$model <- 'AML patients'

			cd34.c.df <- cd34.b.df
			cd34.c.df$replicate <- NULL
			cd34.c.df$replicate2 <- NULL
			cd34.c.df$model <- 'HSPCs'
			head(patients.df)
			head(cd34.c.df)
			
			colnames(patients.df) <- colnames(cd34.c.df)
			both.df <- rbind(patients.df,cd34.c.df)
			dim(both.df)

			summary(both.df$condition)
			summary(both.df$condition)

			write.csv(both.df,
					 file = paste0(printdir,'S6F_both.txt'),
					 quote = TRUE,
					 row.names = FALSE)
			saveRDS(both.df,
					 file = paste0(printdir,'S6F_both.rds'))

		# keep axes free within facet using facet_grid2			
			devtools::install_github("teunbrand/ggh4x")
			date <- '23_03_23'			
			both.df <- readRDS('23_03_23.S6F_both.rds')

		# define combined colour scale
			conditions <- c('CTRL\nHSPCs','STAG1\nKD','STAG2\nKD','RAD21\nKD') %>%
					setNames(c('siCtrl','stag1','stag2','rad21'))	
			AML_cohesin.c1 <- c('#B22222','#B8860B','#2E8B57','#C71585') %>%
					setNames(conditions)
			conditions2 <- c('CTRL\nAML','STAG2\nmut','RAD21\nmut') %>%
					setNames(c('AMLctr','AMLSA2','AMLRAD21'))
			AML_cohesin_patients.c1 <- c('#B22222','#2E8B57','#C71585') %>%
					setNames(conditions2)
			AML_cohesin.c2 <- c(AML_cohesin.c1,AML_cohesin_patients.c1)

		# ggplot
			set.seed(10) ; q2 <- qplot(data=both.df, 
				x=condition, y=mean, 
				col=condition, geom='blank') +
			geom_bar(data=subset(both.df,replicate3==1),
				aes(y = mean_mean),
				col='black',fill='darkgrey',stat="identity", 
				position=position_dodge(), width = 0.75) +
			geom_point(data=subset(both.df,replicate3==1),
				aes(y = mean_mean),size=1,col='black') +
			geom_errorbar(data=subset(both.df,replicate3==1),
				aes(ymin = ymin_mean,
					ymax = ymax_mean),
				width = 0.15, size=1,col='black') +
			geom_jitter(size=1.5,width = 0.15) +			
			ggh4x::facet_grid2(category~model, 
				scales = "free", independent = "y") +
			theme_custom(psize=8) + 
			theme(axis.title.x=element_blank(),legend.position = "none") +
			labs(title='Differential loops in STAG2 mut patients', y=paste0('Loop strength')) +
			# ylim(c(0,4)) +
			scale_colour_manual('',values = AML_cohesin.c2) +
			scale_fill_manual('',values = AML_cohesin.c2)

			set.seed(10)
			ggsave(paste0(printdir,'S6F_both.b.png'),
					q2, width = 10, height = 10, dpi=600, units='cm')
			set.seed(10)
			ggsave(paste0(printdir,'S6F_both.b.pdf'),
				q2, width = 10, height = 10, dpi=600, units='cm', useDingbats=FALSE)

		# alternatively: with order of categories reversed
			both.df$category <- factor(both.df$category,levels=c('Weakened','Strengthened'),ordered=TRUE)

			set.seed(10) ; q2 <- qplot(data=both.df, 
				x=condition, y=mean, 
				col=condition, geom='blank') +
			geom_bar(data=subset(both.df,replicate3==1),
				aes(y = mean_mean),
				col='black',fill='darkgrey',stat="identity", 
				position=position_dodge(), width = 0.75) +
			geom_point(data=subset(both.df,replicate3==1),
				aes(y = mean_mean),size=1,col='black') +
			geom_errorbar(data=subset(both.df,replicate3==1),
				aes(ymin = ymin_mean,
					ymax = ymax_mean),
				width = 0.15, size=1,col='black') +
			geom_jitter(size=1.5,width = 0.15) +			
			ggh4x::facet_grid2(category~model, 
				scales = "free", independent = "y") +
			theme_custom(psize=8) + 
			theme(axis.title.x=element_blank(),legend.position = "none") +
			labs(title='Differential loops in STAG2 mut patients', y=paste0('Loop strength')) +
			# ylim(c(0,4)) +
			scale_colour_manual('',values = AML_cohesin.c2) +
			scale_fill_manual('',values = AML_cohesin.c2)

			set.seed(10)
			ggsave(paste0(printdir,'S6F_both.c.png'),
					q2, width = 10, height = 10, dpi=600, units='cm')
			set.seed(10)
			ggsave(paste0(printdir,'S6F_both.c.pdf'),
				q2, width = 10, height = 10, dpi=600, units='cm', useDingbats=FALSE)
