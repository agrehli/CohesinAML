#!/usr/bin/env Rscript
# S7B_barplots.R
# 23_01_03 L 733-794 and 527-649, as well as 1594-1652 and 1369-1510

### Plot Aggregate Peak Analysis (APA) results from coolpup.py over individual patient or replicate Hi-C matrices of Ctrl and Cohesin deficient cells/patients.
### Hi-C signal was aggregated over loop anchors overlapping 'weakened' or 'strengthened' loops
### between patients with STAG2 mutations compared to control patients
### Besides being differential loops, their anchors or 'edges' have to overlap Cohesin-associated enhancers
### loops are not divided by size, but analyses taking this into account were also performed with similar results
### The APA analyses were performed using coolpup.py in S7B_coolpuppy.sh

# R version: 4.2.2
# colorscale: logarithmic, centred on 1, log Fold Change range [0.25-4]


# 1. Re-plot APAs over differential loops overlaping cohesin-associated enhancers using individual replicate CD34 matrices
	# [2022-12-08 16:19] folowing up from code in 22_11_27

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
		date <- '23_01_12' #'23_01_04' #'23_01_03'
		indir <- paste0(wdir,'/APA/DifferentialAnchorsEnhancer/')
		printdir <- paste0(indir,date,'.')

	# Read all APA plots together
		#find matching coolpup.py output
		paths <- list.files(path = indir, recursive = TRUE,
				pattern = "\\.bedpe$", full.names = TRUE)
		length(paths) #[1] 3008

	# Start by plotting individual CD34 samples:
		paths <- paths[grepl('\\.CD34_',paths)] #only individual replicate CD34 samples
		length(paths) #[1] 1600

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
		dh(APAs.df) #2689600 when together
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

	# Define names: CD34 model - replicates
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
		
		conditions <- c('si Ctrl','SA1 KD','SA2 KD','RAD21 KD') %>%
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

	# colour palettes to try:
		palettes2 <- c('RdBu')

# S7B - Subpanel A - Barplots - CD34 model - individual matrices
	
	# # Re-plot as 'Replicate vs Condition' barplot
		comparisons <- unique(cd34.df$comparison) # [1] "SA2mutvsCTRL" #only one 
		types <- unique(cd34.df$type) #[1] "CohesinAssEnhancersOverlap" "H3K27acOverlap"  
		edges <- unique(cd34.df$edge) #[1] "both"   "either" "minus"  "plus" #results are the same with either...
		categories <- unique(cd34.df$category)
		
		for (i in comparisons){
			for (j in types){
				for (l in edges){
					for (m in categories){
						# subset
							print(paste0(i,'_',j,'_',l))

							cd34.b.df <- cd34.df %>%
								filter(comparison==i,type==j,edge==l,category==m,Var1=='ID0',Var2=='ID0',SizeGroup=="All")
							print(table(cd34.b.df$sample))
							print(table(cd34.b.df$condition))
							print(table(cd34.b.df$category))

						# error bar:
							error_bar <- cd34.b.df %>%
							group_by(condition) %>%
							summarise(mean_mean=mean(mean,na.rm=TRUE),
									  sd_mean=sd(mean,na.rm=TRUE),
									  ymin_mean=mean_mean-sd_mean,
									  ymax_mean=mean_mean+sd_mean)
							cd34.b.df <- cd34.b.df %>%
								left_join(error_bar,by='condition')

						# barplot
							q2 <- qplot(data=cd34.b.df, 
								x=condition, y=mean, 
								col=condition, geom='blank') +
							geom_bar(data=subset(cd34.b.df,replicate3==1),
								aes(y = mean_mean),
								col='black',fill='darkgrey',stat="identity", 
								position=position_dodge(), width = 0.75) +
							geom_point(data=subset(cd34.b.df,replicate3==1),
								aes(y = mean_mean),size=1,col='black') +
							geom_errorbar(data=subset(cd34.b.df,replicate3==1),
								aes(ymin = ymin_mean,
									ymax = ymax_mean),
								width = 0.15, size=1,col='black') +
							geom_jitter(size=1,width = 0.15) +
							theme_custom(psize=8) + 
							theme(axis.title.x=element_blank(),legend.position = "none") +
							labs(title=paste0(m, ' loops comparing ',i,'\npadj <= 0.05 & log2FoldChange >0.585'),
								y=paste0('Loop strength')) +
							scale_colour_manual('',values = AML_cohesin.c1) +
							scale_fill_manual('',values = AML_cohesin.c1)

							ggsave(paste0(printdir,'AML2021.MeanLoopStrength.DifferentialAnchorsEnhancer.',i,'.overlappingType.',j,'.edgeType.',l,'.samples.individual_CD34_matrices.SizeGroup.All.category.',m,'.png'), 
									q2, width = 8, height = 8, dpi=600, units='cm')
							ggsave(paste0(printdir,'AML2021.MeanLoopStrength.DifferentialAnchorsEnhancer.',i,'.overlappingType.',j,'.edgeType.',l,'.samples.individual_CD34_matrices.SizeGroup.All.category.',m,'.pdf'), 
								q2, width = 8, height = 8, dpi=600, units='cm', useDingbats=FALSE)
					}
				}
			}
		}
		

# 2. Re-plot APAs over differential loops using individual Cohesin patient matrices

	# [2022-12-20 08:26] following up from code in 22_11_27
	
	# Read all APA plots together again, but now focusing on individual patient matrices
		
		#find matching coolpup.py output
		paths <- list.files(path = indir, recursive = TRUE,
				pattern = "\\.bedpe$", full.names = TRUE)
		length(paths) #[1] 3008

	# Now individual Cohesin patient samples together
		paths <- paths[!grepl('\\.CD34_',paths)]
		paths <- paths[!grepl('merged_',paths)]
		length(paths) #[1] 960

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
		dh(APAs.df) #1613760
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

	# Define names: Cohesin project - Cohesin patients - individual matrices
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
			
			conditions <- c('AML\nCtrl','AML\nSA2 \u394','AML\nRAD21 \u394') %>%
				setNames(c('AMLctr','AMLSA2','AMLRAD21'))
			
			AML_cohesin_patients.c1 <- c('#df1413','#623a87','#58812f') %>%
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
				arrange(desc(sample))
		table(patients.df$patient3)
		table(patients.df$sample)
		table(patients.df$patient)
		table(patients.df$condition)
		# table(patients.df$loops2)
		table(patients.df$category)

	# colour palettes to try:
		palettes2 <- c('RdBu')

# S7B - Subpanel B - Cohesin mutation patient samples - individual matrices
	
# # Re-plot as 'Patient vs Condition' barplot
		comparisons <- unique(patients.df$comparison)
		types <- unique(patients.df$type) #[1] "CohesinAssEnhancersOverlap" "H3K27acOverlap"  
		edges <- unique(patients.df$edge) #[1] "both"   "either" "minus"  "plus" #results are the same with either...
		categories <- unique(patients.df$category)
		
		for (i in comparisons){
			for (j in types){
				for (l in edges){
					for (m in categories){
						# subset
							print(paste0(i,'_',j,'_',l,'_',m))

							patients.b.df <- patients.df %>%
								filter(comparison==i,type==j,edge==l,category==m,Var1=='ID0',Var2=='ID0',SizeGroup=="All")
							print(table(patients.b.df$sample))
							print(table(patients.b.df$condition))
							print(table(patients.b.df$category))

						# error bar:
							error_bar <- patients.b.df %>%
							group_by(condition) %>%
							summarise(mean_mean=mean(mean,na.rm=TRUE),
									  sd_mean=sd(mean,na.rm=TRUE),
									  ymin_mean=mean_mean-sd_mean,
									  ymax_mean=mean_mean+sd_mean)
							patients.b.df <- patients.b.df %>%
								left_join(error_bar,by='condition')

						# barplot
							q2 <- qplot(data=patients.b.df, 
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
							geom_jitter(size=1,width = 0.15) +
							theme_custom(psize=8) + 
							theme(axis.title.x=element_blank(),legend.position = "none") +
							labs(title=paste0(m, ' loops comparing ',i,'\npadj <= 0.05 & log2FoldChange >0.585'),
								y=paste0('Loop strength')) +
							scale_colour_manual('',values = AML_cohesin_patients.c1) +
							scale_fill_manual('',values = AML_cohesin_patients.c1)

							ggsave(paste0(printdir,'AML2021.MeanLoopStrength.DifferentialAnchorsEnhancer.',i,'.overlappingType.',j,'.edgeType.',l,'.samples.individual_Cohesin_patients.SizeGroup.All.category.',m,'.png'),
									q2, width = 8, height = 8, dpi=600, units='cm')
							ggsave(paste0(printdir,'AML2021.MeanLoopStrength.DifferentialAnchorsEnhancer.',i,'.overlappingType.',j,'.edgeType.',l,'.samples.individual_Cohesin_patients.SizeGroup.All.category.',m,'.pdf'),
								q2, width = 8, height = 8, dpi=600, units='cm', useDingbats=FALSE)
					}
				}
			}
		}
			