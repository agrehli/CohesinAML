#!/tools/conda/bin/python
# Modified so that read_feather works [/usr/bin/env python2]
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 11:15:36 2018

@Name: heuristic_filter.py
@Author: Roger Mulet
@Version: 0.2
@Description: Heuristic filtering
@Usage: 

"""
import sys
import pandas as pd
import warnings
import time
import os
import pybedtools
from natsort import index_natsorted,order_by_index
warnings.filterwarnings("ignore",'This pattern has match groups')
warnings.filterwarnings("ignore",'read_table')
pd.options.mode.chained_assignment = None  # default='warn' --> Remove pandas warning
# pandas.set_option('max_columns',10)

# annovar = '/data/bas/1_enhancer_analysis/exome/bam/variant_calling/2190_180720_A00383_0011_AHFL2NDMXX_pindel.txt'

#--------------------------------------------------------------------------------------------------
#Testing parameters
#--------------------------------------------------------------------------------------------------

# MIN_DEPTH_TUMOR = 8
# MIN_ALT_FREQ_TUMOR = 0.1
# MIN_ALT_DEPTH_TUMOR = 4
# POPMAX = 0.1
# 
# MIN_DEPTH_CONTROL=10
# MAX_ALT_DEPTH_CONTROL=5
# MAX_ALT_FREQ_CONTROL=0.05
# DIF_ALT_FREQ_CONTROL=0.08
# 
# MAX_AA_RATIO = 0.1
# MAX_EXCEED_RATIO = 0.1
# MAX_SMALL_RATIO = 0.2
# MAX_NON_SNP = 0.15
# MIN_HQ_RATIO = 0.15
# MAX_CLIPPED_ALT = 0.80

#normal = 'no'
#ftype = 'exonic_splicing'

#annovar = '/home/roger/Documents/Projects/Aniko/ALL_44_combined.txt'
#hotspots = '/data/bas/1_enhancer_analysis/Hotspots_TS_Myeloid_FK05072018.txt'
#
#bias_capture = '/data/bas/1_enhancer_analysis/exome/capture/MedExome_hg19_capture_targets.bed'
#blacklist = '/data/bas/1_enhancer_analysis/exome/variant_analysis/thijs1/raw/blacklist_genes.txt'
#whitelist = '/data/bas/1_enhancer_analysis/exome/variant_analysis/thijs1/raw/whitelist_claudia.txt'
#blacklist_pos = '/data/bas/1_enhancer_analysis/exome/variant_analysis/thijs1/raw/blacklist_claudia.txt'
#whitelist_pos = '/data/bas/1_enhancer_analysis/exome/variant_analysis/thijs1/raw/whitelist_pos_claudia.txt'
#repeats = '/data/bas/1_enhancer_analysis/exome/variant_analysis/thijs1/raw/repeats_simple.bed'

# https://stackoverflow.com/questions/43871626/pandas-read-hdf-very-slow-for-non-numeric-data
#control = '/data/bas/1_enhancer_analysis/exome/control/variant_analysis/T00839_germline_combined.feather'
# %timeit pfeather = pd.read_feather('/data/bas/1_enhancer_analysis/exome/control/variant_analysis/T00839_germline_combined.feather')

#--------------------------------------------------------------------------------------------------
#Helper functions
#--------------------------------------------------------------------------------------------------

def check_files(input_file):

	if isinstance(input_file,list): # Recursively applied if the file itself is a list
		[check_files(item) for item in file]
		
	else:	
		if input_file is not None and os.path.exists(input_file) == False:	
			raise ValueError("File {} not found!".format(input_file))

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')		
			
#--------------------------------------------------------------------------------------------------
#Get input arguments
#--------------------------------------------------------------------------------------------------

import argparse

parser = argparse.ArgumentParser(description='Script to filter variants annotated with Annovar and AnnotateBamStatistics. The program uses high quality reads (HQ),\
                                 where HQ is defined',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(dest='annovar',type=str,help="Annovar-annotated TXT file containing variants of interest")
parser.add_argument('-t','--type',required=False,type=str,default='relevant',help="What variants should be kept? You can provide a comma separated list [relevant/exonic/exonic_splicing/regulatory]")
parser.add_argument('-n','--normal',required=False,default=False,action='store_true',help="Indicate whether we have a normal sample")
parser.add_argument('--normal_suffix',required=False,default='.1',help="Suffix to identify normal samples")
parser.add_argument('-c','--control',required=False,type=str,help="Annovar annotated file used as a control. You can provide an HDF5 file for speed")
parser.add_argument('-b','--blacklist',required=False,type=str,help="Blacklist containing genes to be discarded")
parser.add_argument('-B','--blacklist_positions',required=False,type=str,help="Blacklist containing genomic coordinates of REGIONS to be discarded")
parser.add_argument('-R','--repeats',required=False,type=str,help="Repeat mask for filtering ONLY indels which fall entirely in any of these regions. Low complexity (simple repeats) are recommended")
parser.add_argument('-w','--whitelist',required=False,type=str,help="Whitelist containing genes to be considered. However, variants in those genes might be removed by other filters")
parser.add_argument('-W','--whitelist_positions',required=False,type=str,help="Whitelist containing regions to be considered. However, variants in those regions might be removed by other filters")
parser.add_argument('-K','--keep_positions',required=False,type=str,help="Keep-list containing genomic coordinates to be kept (after basic depth/frequency filtering).\
                    Please note that this overrides any other filtering")
parser.add_argument('--hotspots',required=False,type=str,help="File containing known hotspots of myeloid mutations. A new column will be added and these variants will be kept")
parser.add_argument('--flt3',default='yes',type=str2bool,help="Activate the FLT3 special mode: FLT3 variants will be kept regardless of BAM statistics")
parser.add_argument('-s','--somatic',default=False,action='store_true',help="Attempt to do somatic filtering: variants with benign prediction")
parser.add_argument('-S','--self_chain',type=str,default=None,help="Path to SelfChain file with genome self alignment")
parser.add_argument('-p','--plus',action='store_true',help="Filter with metrics from AnnotateBamStatistics")
parser.add_argument('-o','--output',required=False,type=str,help="Basename of the output file. If missing, output to stdout")
parser.add_argument('-m','--majority',required=False,type=int,default=0,help="Select only samples detected at least by N callers")

parser.add_argument('--min_depth_tumor',type=int,default=8,help="At least X reads in the tumor")
parser.add_argument('--min_alt_freq_tumor',type=float,default=0.095,help="At least a variant allele frequency of X in the tumor (divided by 2 since there are 2 alleles)")
parser.add_argument('--min_alt_freq_pindel',type=float,default=0.15,help="At least a variant allele frequency of X in the tumor, as calculated by Pindel for indels")
parser.add_argument('--min_alt_depth_tumor',type=int,default=4,help="At least X reads showing the variant (Get rid of most of the bogus loci)")
parser.add_argument('--min_hq_ratio',type=float,default=0.15,help="At least X fraction of the total reads must be of high quality")
parser.add_argument('--popmax',type=float,default=0.01,help="Maximum population frequency lower than X. If unavailable, it is taken anyway. Set to 1 to include all")
parser.add_argument('--min_depth_control',type=int,default=8,help="We would like to have at least X reads in the control")
parser.add_argument('--max_alt_depth_control',type=int,default=2,help="At most X variant read (or a maximum VAF of X, see below)")
parser.add_argument('--max_alt_freq_control',type=float,default=0.05,help="At most a variant allele frequency of X in the control")
parser.add_argument('--dif_alt_freq',type=int,default=0.20,help="The alt freq of the tumor should be at least X higher than that of control")

parser.add_argument('--max_aa_ratio',type=float,default=0.10,help="Maximum proportion of reads with an alternative alignment")
parser.add_argument('--max_exceed_ratio',type=float,default=0.10,help="Maximum proportion of reads with superior alternative alignment score (XS)")
parser.add_argument('--max_small_ratio',type=float,default=0.20,help="Maximum proportion of reads with slightly smaller alternative alignment score (normally 5)")
parser.add_argument('--max_non_snp',type=float,default=0.25,help="Maximum ratio of reads with N (normally 2) mismatches not reported in SNPdb")
parser.add_argument('--max_clipped_alt',type=float,default=0.60,help="Maximum ratio of clipped reads supporting ALT allele")
parser.add_argument('--max_proximity',type=float,default=None,help="Maximum distance allowed between SNPs; those that are close will be discarded (recommended 5)")
parser.add_argument('--max_cluster',type=float,default=3,help="Maximum size of SNP clusters; larger clusters are going to be removed")

group = parser.add_mutually_exclusive_group(required=False)
group.add_argument('--bias',required=False,type=str2bool,default='yes',help="Enable or disable bias filter, i.e. based on whether all reads align to same strand")
group.add_argument('--bias_capture',required=False,type=str,help="Enable the bias filter with capture correction -- i.e. regions outside the capture will not be filtered\
                    because the outer borders of the capture tend to exhibit only reads from one strand")

args = parser.parse_args()

# Files
annovar = args.annovar
control = args.control
whitelist = args.whitelist
whitelist_pos = args.whitelist_positions
blacklists = args.blacklist.split(",") if args.blacklist is not None else None
blacklists_pos = args.blacklist_positions.split(",") if args.blacklist_positions is not None else None
repeats = args.repeats
bias_capture = args.bias_capture
hotspots = args.hotspots
self_chain = args.self_chain
keep_pos = args.keep_positions

# Check files
for file in [annovar,control,whitelist,whitelist_pos,blacklists,blacklists_pos,\
repeats,bias_capture,hotspots,self_chain,keep_pos]:

	check_files(file)

# Options
ftype = args.type
normal = args.normal
output = args.output
bias = args.bias
stats_plus = args.plus
flt3_mode = args.flt3
somatic = args.somatic

NORMAL_SUFFIX = args.normal_suffix
MAJORITY_VOTING = args.majority
MIN_DEPTH_TUMOR = args.min_depth_tumor
MIN_ALT_FREQ_TUMOR = args.min_alt_freq_tumor
MIN_ALT_DEPTH_TUMOR = args.min_alt_depth_tumor
MIN_HQ_RATIO = args.min_hq_ratio
POPMAX = args.popmax

MIN_DEPTH_CONTROL = args.min_depth_control
MAX_ALT_DEPTH_CONTROL = args.max_alt_depth_control
MAX_ALT_FREQ_CONTROL = args.max_alt_freq_control
DIF_ALT_FREQ_CONTROL = args.dif_alt_freq

MAX_AA_RATIO = args.max_aa_ratio
MAX_EXCEED_RATIO = args.max_exceed_ratio
MAX_SMALL_RATIO = args.max_small_ratio
MAX_NON_SNP = args.max_non_snp
MAX_CLIPPED_ALT = args.max_clipped_alt
MAX_PROXIMITY = args.max_proximity
MAX_CLUSTER = args.max_cluster
MAX_GENOTYPE_ERROR = 0.4
MIN_ALT_DEPTH_TUMOR_PINDEL = 6
MIN_ALT_FREQ_TUMOR_PINDEL = args.min_alt_freq_pindel
MIN_DEPTH_TUMOR_PINDEL = 12

start_time = time.process_time()

#==================================================================================================
#Import file and filter entries
#==================================================================================================

# Identify columns and add new ones: Keep, ntools, depth
variants = pd.read_csv(annovar,sep='\t')
TOOLS = variants.columns.isin(['bcftools','mutect','mutect2','strelka','hcaller','varscan','pindel','vardict','pisces','caveman'])
DTYPES = variants.dtypes.to_dict()
variants['ntools'] = variants.iloc[:,TOOLS].sum(axis=1)
variants['Keep'] = [False]*len(variants)
column_cosmic = variants.columns[variants.columns.str.contains(pat='cosmic',case=False)][0]
column_cosmic_flagged = variants.columns[variants.columns.str.contains(pat='cosmic.*flagged',case=False)]

try:
	DEPTH_COL = variants.columns[variants.columns.str.contains('\)depth')].to_list()[0]
except IndexError:
	DEPTH_COL = variants.columns[variants.columns.isin(['depth'])].to_list()[0]
variants = variants.rename(columns={DEPTH_COL:'depth'})

# If left DF is empty at any point, the order of the columns will change
# https://github.com/pandas-dev/pandas/issues/9937#event-2070332719


if 'alt depth' + NORMAL_SUFFIX in variants.columns and normal == False:

	sys.stderr.write("WARNING: Columns ending with {} are detected. Did you forget to indicate that a normal is present?\n".format(NORMAL_SUFFIX))
	sys.stderr.flush()     

#--------------------------------------------------------------------------------------------------
#Functional filter: 1st filter that removes most of the variants
#--------------------------------------------------------------------------------------------------

if ftype == 'exonic':
    
    variants_rel = variants[(variants['Func.refGene'].str.contains('exonic')) & (~variants['Func.refGene'].str.contains('ncRNA')) &\
                                   (~variants['ExonicFunc.refGene'].str.contains('^synonymous',na=False))]

elif ftype == 'exonic_splicing':
    
    variants_rel = variants[(variants['Func.refGene'].str.contains('exonic') | variants['Func.refGene'].str.contains('splicing')) &\
    (~variants['Func.refGene'].str.contains('ncRNA',na=False)) & (~variants['ExonicFunc.refGene'].str.contains('^synonymous',na=False)) ]
    
elif ftype == 'regulatory':
    
    variants_rel = variants[(variants['Func.refGene'].str.contains('UTR') | variants['Func.refGene'].str.contains('upstream') |\
                                       variants['Func.refGene'].str.contains('splicing',na=False)) & (~variants['Func.refGene'].str.contains('ncRNA',na=False))]

else:
    
    variants_rel = variants
    
sys.stderr.write("- {}/{} variants have been removed by functional filter {}. {} remain.\n".format(
	len(variants)-len(variants_rel),len(variants),ftype,len(variants_rel))); sys.stderr.flush()     

#--------------------------------------------------------------------------------------------------
# Filter by majority voting (if enabled)
#--------------------------------------------------------------------------------------------------

if MAJORITY_VOTING > 0:

    LESS_THAN_MAJORITY = (variants_rel['ntools'] < MAJORITY_VOTING)
    variants_rel = variants_rel[~LESS_THAN_MAJORITY]
    	
    sys.stderr.write("- {}/{} variants have been removed by majority voting filtering. {} remain.\n".format(
	sum(LESS_THAN_MAJORITY),len(variants_rel)+sum(LESS_THAN_MAJORITY),len(variants_rel) ) ); sys.stderr.flush()     

#--------------------------------------------------------------------------------------------------
#Depth and population filters
#--------------------------------------------------------------------------------------------------

variants_rel['indel'] = ( (variants_rel['Ref'] == "-") | (variants_rel['Alt'] == "-") |\
            (variants_rel['Ref'].map(len) != variants_rel['Alt'].map(len) ) )

## Recalculate PopFreqMax based on various large scale genome projects
sys.stderr.write("INFO: PopFreqMax has been recalculated using other population databases\n");  sys.stderr.flush() 
pop_columns = variants_rel.columns.str.contains("PopFreqMax|ExAC|gnomAD") # or isin(['PopFreqMax','1000G_ALL'...
variants_rel.loc[:,pop_columns] = variants_rel.loc[:,pop_columns].apply(pd.to_numeric,errors='coerce')
# Insert PopFreqMax before the first population column
variants_rel.insert(min([id for id,pop in enumerate(pop_columns) if pop]),\
'PopFreqMax',variants_rel.iloc[:,pop_columns].max(axis=1,skipna=True)) 
# Convert PopFreqMax to numeric (it contains empty strings)
variants_rel['PopFreqMax'] = pd.to_numeric(variants_rel['PopFreqMax'],errors='coerce')

if flt3_mode:

    # Select only FLT3 variants_rel with nonframeshift insertion
    keep_variants = variants_rel[(variants_rel['Gene.refGene'] == 'FLT3') & (variants_rel['ExonicFunc.refGene'] == 'nonframeshift insertion')]
    # Take those that affect exons 14/15, multiple of 3 and select the largest one
    try:
        keep_variants = keep_variants[ (keep_variants['AAChange.refGene'].str.contains(pat="exon14|exon15")) &\
                                      (keep_variants['Alt'].map(len) % 3 == 0 )]

        if len(keep_variants) > 0 and 'pindel_alt' in variants_rel.columns and sum(keep_variants['pindel_alt'].isna()) == 0:

            variants_rel.loc[keep_variants['pindel_alt_freq'].idxmax(),'Keep'] = True

        elif len(keep_variants) > 0:
            # Prioritize pindel results, but higher score if we also see it with other tools
            TOP_ROW = keep_variants.sort_values(by=['HQ alt depth','pindel','ntools','HQ depth'],ascending=False, axis=0).index[0]          
            variants_rel.loc[TOP_ROW,'Keep'] = True # keep_variants['HQ depth'].idxmax()

    except ValueError:
        pass   

# print(variants_rel[variants_rel['Gene.refGene']=='RAD21'].loc[:,['Ref','Alt','depth','alt freq','pindel_alt_freq','pindel_depth']] );sys.exit()

## TUMOR-ONLY FILTER ##

LEN_START = len(variants_rel)        
NOT_BOGUS = ( (variants_rel['Keep'] == True) | (variants_rel['Ref'] != variants_rel['Alt']) & (variants_rel['HQ depth'] >= MIN_DEPTH_TUMOR) &\
                             (variants_rel['HQ alt freq'] >= MIN_ALT_FREQ_TUMOR) & (variants_rel['HQ alt depth'] >= MIN_ALT_DEPTH_TUMOR))                           
NOT_BOGUS_INDEL = ( (variants_rel['Keep'] == True) | ( variants_rel['indel'] == True ) & (variants_rel['depth'] >= MIN_DEPTH_TUMOR*1.2) &\
                   (variants_rel['alt freq'] >= MIN_ALT_FREQ_TUMOR*1.2) & (variants_rel['alt depth'] >= MIN_ALT_DEPTH_TUMOR*1.2))                           
NOT_BOGUS = (NOT_BOGUS | NOT_BOGUS_INDEL)

LOW_HQ_RATIO = (variants_rel['HQ ratio'] < MIN_HQ_RATIO)

if 'pindel_alt' in variants_rel.columns: # We set more stringent thresholds because it's reads instead of fragments
    
    NOT_BOGUS_PINDEL = ( (variants_rel['pindel_depth'] != 0) & ( (variants_rel['Ref'].map(len) > 2) | (variants_rel['Alt'].map(len) > 2) ) &\
                        ( (variants_rel['alt depth'] == 0) | (variants_rel['depth'] == variants_rel['alt depth']) ) &\
   (variants_rel['pindel_depth'] >= MIN_DEPTH_TUMOR_PINDEL) & (variants_rel['pindel_alt_freq'] >= MIN_ALT_FREQ_TUMOR_PINDEL) & (variants_rel['pindel_alt'] >= MIN_ALT_DEPTH_TUMOR_PINDEL) )
                                 
    # Filter out identical variants_rel, normally introduced by Pindel
    PINDEL_DUPLICATES = ( (variants_rel['End'].duplicated(keep=False)) & (variants_rel['pindel'] == 1) &\
                         (variants_rel['ExonicFunc.refGene'].str.contains('substitution')) )                                    
    
    NOT_BOGUS = ( ( NOT_BOGUS | NOT_BOGUS_PINDEL) & ~ PINDEL_DUPLICATES )

variants_rel = variants_rel[ NOT_BOGUS & ~ LOW_HQ_RATIO ]

# Remove non primary chromosomes for Mathijs program
variants_rel = variants_rel[variants_rel['Chr'].str.contains('chr[0-9XYMm]+$')]
# Remove extremely long indels
variants_rel = variants_rel[variants_rel['Alt'].map(len) < 500]
variants_rel = variants_rel[variants_rel['Ref'].map(len) < 500]

sys.stderr.write(f"- {LEN_START-len(variants_rel)}/{LEN_START} variants have been removed by depth and frequency filters. {len(variants_rel)} remain.\n")
sys.stderr.write(f"[MIN_DEPTH_TUMOR: {MIN_DEPTH_TUMOR} (Pindel: {MIN_DEPTH_TUMOR_PINDEL}) | MIN_ALT_DEPTH_TUMOR: {MIN_ALT_DEPTH_TUMOR} (Pindel: {MIN_ALT_DEPTH_TUMOR_PINDEL}) | \
MIN_ALT_FREQ_TUMOR: {MIN_ALT_FREQ_TUMOR} (Pindel: {MIN_ALT_FREQ_TUMOR_PINDEL})]\n")
sys.stderr.flush() 

if hotspots is not None and len(variants_rel) > 0: # Otherwise we have issues with merging...

    HOTSPOTS = pd.read_csv(hotspots,sep='\t')
    aachange = variants_rel['AAChange.refGene'].str.contains(pat='|'.join(HOTSPOTS['KeyScript'].tolist()).replace('.','.*'))
    variants_rel['Hotspot'] = aachange
    variants_rel['Keep'] = (variants_rel['Keep'] == True) | (variants_rel['Hotspot'] == True)   

cancer_hematology = variants_rel[column_cosmic].str.extract(
        pat='(\d+)\(haematopoietic_and_lymphoid_tissue\)',expand=False).astype(float).fillna(0)

if 'CLNSIG' in variants_rel.columns:
    variants_rel['Keep'] = (variants_rel['Keep'] == True) | ( (variants_rel['CLNSIG'].str.contains('Pathogenic')) &
    (cancer_hematology >= 20) )

if 'PopFreqMax' in variants_rel.columns:
	# We remove common variants not related to hematological malignancies. Here we do not take into account alt freq
	# because even if they are not germline variants, if they are common they are unlikely to be cancer drivers
	# We add a small exception for CHIP genes because some DNMT3A variants are quite frequent

    NOT_SNP = ( (variants_rel['PopFreqMax'] <= POPMAX) | (variants_rel['PopFreqMax'].isnull()) | 
            ( (cancer_hematology >= 5) & (variants_rel['PopFreqMax'] <= 0.1) ) |  
			( (cancer_hematology >= 2) & (variants_rel['PopFreqMax'] <= 0.1) & (variants_rel['Gene.refGene'].isin(['DNMT3A','TET2','ASXL1'])) ) )     
    variants_rel = variants_rel[ (variants_rel['Keep'] == True) | NOT_SNP ]

    sys.stderr.write("- {}/{} variants have been removed by population frequency (freq > {}). {} remain.\n".format(
    sum(~NOT_SNP),len(NOT_SNP),POPMAX,len(variants_rel))); sys.stderr.flush()  

## BIAS FILTER ##

# Filter on bias, i.e. exclude variants with all reads mapping to the same strand  
# To reduce the chances of removing true variants with few reads, we add a condition of minimum depth

if bias or bias_capture is not None:
    
    BIAS_SITES = (variants_rel['HQ alt depth'] > 10) & ( (variants_rel['HQ alt bias'] == 0) | (variants_rel['HQ alt bias'] == 1) ) 
    
if bias_capture is not None and len(variants_rel) > 0:
    # To further limit false negatives, we only filter by bias variants in the exome capture. Regions close to the borders
    # of the capture may display biased reads because only one strand is captured   
       
    CAPTURE = pd.read_csv(bias_capture,sep='\t')
    CAPTURE.columns = ['chr','start','end','gene'] + ["col" + str(x) for x in range(4,CAPTURE.shape[1])]
    BORDER_START = pd.concat([CAPTURE['chr'],(CAPTURE['start']-200),CAPTURE['start'].rename('end')],axis=1)
    BORDER_END = pd.concat([CAPTURE['chr'],CAPTURE['end'].rename('start'),(CAPTURE['end']+200)],axis=1)
    BORDERS = pd.concat([BORDER_START,BORDER_END],ignore_index=True)    
    BORDERS_BED = pybedtools.BedTool.from_dataframe(df=BORDERS)
    
    variants_bedtool = pybedtools.BedTool.from_dataframe(df=variants_rel.iloc[:,0:6])   
    variants_in_borders = variants_bedtool.intersect(BORDERS_BED,c=True,f=1) # Overlap with borders
    IN_BORDERS = variants_in_borders.to_dataframe().set_index(variants_rel.index)
    IN_BORDERS = (IN_BORDERS.iloc[:,-1] != 0)
    BIAS_SITES = ( BIAS_SITES & ~IN_BORDERS)        
   
    variants_rel = variants_rel [ ~ BIAS_SITES | variants_rel['Keep'] == True ]    
    sys.stderr.write("- {}/{} variants have been removed by bias filter. {} remain.\n".format(
            sum(BIAS_SITES),len(BIAS_SITES),len(variants_rel))); sys.stderr.flush()      

#--------------------------------------------------------------------------------------------------
#Tumor vs Normal samples if available
#--------------------------------------------------------------------------------------------------

if control is not None: # If an external control is provided, we add it to the list    
    
    if 'HQ depth.1' in variants_rel.columns.tolist():
        
        sys.stderr.write('Warning: Columns from control sample already detected.\n'); sys.stderr.flush()
    
    if control.lower().endswith('.feather'):        
        CONTROL = pd.read_feather(control)        
    elif control.lower().endswith('.hdf5'):        
        CONTROL = pd.read_hdf(control)        
    else:
        CONTROL = pd.read_csv(control,sep='\t')
    
    # Keep only relevant columns: variant position and HQ bam statistics
    DEPTH_COL_CTRL = CONTROL.columns[CONTROL.columns.str.contains('\)depth')].to_list()
    if len(DEPTH_COL_CTRL) > 0:
        CONTROL = CONTROL.rename(columns={DEPTH_COL_CTRL[0]:'depth'})
    variants_control = CONTROL.loc[:,CONTROL.columns.str.contains('^(HQ|depth|alt|Chr|Start|End|Ref|Alt|pindel_)')]      
    variants_rel = variants_rel.merge(variants_control,how='left',sort=False,on=['Chr','Start','End','Ref','Alt'],suffixes=['',NORMAL_SUFFIX])  

if normal == True or control is not None:

	# Replace the ()depth column name with the same pattern used in other columns
    DEPTH_COL_CTRL = variants_rel.columns[variants_rel.columns.str.contains('\)depth')].to_list()    
    if len(DEPTH_COL_CTRL) > 0:
        variants_rel = variants_rel.rename(columns={DEPTH_COL_CTRL[0]:'depth'+NORMAL_SUFFIX})

    cancer_hematology = variants_rel[column_cosmic].str.extract(
	pat='(\d+)\(haematopoietic_and_lymphoid_tissue\)',expand=False).astype(float).fillna(0)

    # We remove variants for which 1) we have data; 2) we certainly know they are in the control    
    PRESENT_IN_CONTROL = (variants_rel['depth.1'] >= MIN_DEPTH_CONTROL) &\
    ((variants_rel['alt depth.1'] >= MAX_ALT_DEPTH_CONTROL) | (variants_rel['alt freq.1'] >= MAX_ALT_FREQ_CONTROL)) &\
    ((variants_rel['alt freq'] - variants_rel['alt freq.1']) <= DIF_ALT_FREQ_CONTROL) 

    # If called by pindel, we remove those with a) relatively high VAF, b) not in COSMIC
    if 'pindel_alt.1' in variants_rel:
        
        PRESENT_PINDEL = ( (variants_rel['pindel_depth.1'] != 0) & ( (variants_rel['pindel_depth.1'] >= MIN_DEPTH_CONTROL) & \
                          (variants_rel['pindel_alt_freq.1'] >= MAX_ALT_FREQ_CONTROL) & (variants_rel['pindel_alt.1'] >= MAX_ALT_DEPTH_CONTROL) &\
                                 ( (variants_rel['pindel_alt_freq'] - variants_rel['pindel_alt_freq.1']) <= DIF_ALT_FREQ_CONTROL*2) ) |\
                                 ( (variants_rel['pindel_alt.1'] > 0) & (cancer_hematology < 5) ) )
        PRESENT_IN_CONTROL = (PRESENT_IN_CONTROL | PRESENT_PINDEL)

    variants_rel = variants_rel[ ~PRESENT_IN_CONTROL ]          
    
    sys.stderr.write("- {}/{} variants have been removed because they were detected in control sample. {} remain.\n".format(
    sum(PRESENT_IN_CONTROL),len(variants_rel)+sum(PRESENT_IN_CONTROL),len(variants_rel)))
    sys.stderr.write(f"[MIN_DEPTH_CONTROL: {MIN_DEPTH_CONTROL} | MAX_ALT_DEPTH_CONTROL: {MAX_ALT_DEPTH_CONTROL} | MAX_ALT_FREQ_CONTROL: {MAX_ALT_FREQ_CONTROL}]\n")
    sys.stderr.flush()

#--------------------------------------------------------------------------------------------------
# Whitelist and blacklist based filters
#--------------------------------------------------------------------------------------------------

# WHITELIST / BLACKLIST    
if whitelist is not None:

    WHITELIST = pd.read_csv(whitelist,header=None)
    PATTERN = "|".join(WHITELIST[0].tolist())
    IN_WHITELIST = variants_rel['Gene.refGene'].str.contains(pat=PATTERN)
    variants_rel = variants_rel[ IN_WHITELIST ]
    
    sys.stderr.write("- {}/{} variants have been removed because genes were not in the supplied whitelist. {} remain.\n".format(
    sum(~IN_WHITELIST),len(IN_WHITELIST),len(variants_rel))) ; sys.stderr.flush() 
    
if whitelist_pos is not None: # It used to be based on single locations, now whole regions

    WHITELIST_BED = pybedtools.BedTool(whitelist_pos)    
    LEN_START = len(variants_rel)
    variants_rel.BedTool.from_dataframe(variants_rel)    
    
    variants_in_regions = variants_bedtool.intersect(WHITELIST_BED,u=True,f=1)
    variants_rel = variants_in_regions.to_dataframe(names=variants_rel.columns.tolist(),na_values='.',dtype=DTYPES)
    os.remove(annovar+'.temp')
    
    sys.stderr.write("- {}/{} variants have been removed because positions were not in the supplied whitelist. {} remain.\n".format(
    LEN_START+len(variants_rel),LEN_START,len(variants_rel))) ; sys.stderr.flush()
    
if keep_pos is not None and len(variants_rel) > 0: # Otherwise we have issues with merging...
    KEEP_POS = pd.read_csv(keep_pos,sep='\t',names=['Chr','Start','End','Other'])
    KEEP_POS = KEEP_POS[~KEEP_POS.duplicated()]
    merged_keep = variants_rel.merge(KEEP_POS,how='left',indicator='_KEEP',on=['Chr','Start','End'],sort=False)    
    merged_keep = merged_keep.set_index(variants_rel.index)
    
    variants_rel['Keep'] = (variants_rel['Keep'] == True) | (merged_keep['_KEEP'] == 'both' )  

if blacklists is not None: 
	for blacklist in blacklists:
		# We DO NOT use patterns to prevent: 1) the removal of similarly named genes; 2) intergenic regions with that name
		BLACKLIST = pd.read_csv(blacklist,header=None,sep="\t")
		IN_BLACKLIST = (variants_rel['Gene.refGene'].isin(BLACKLIST[0])) & (variants_rel['Keep'] == False)
		variants_rel = variants_rel[ ~IN_BLACKLIST ]    
		
		sys.stderr.write("- {}/{} variants have been removed because genes were in blacklist {}. {} remain.\n".format(
		sum(IN_BLACKLIST),len(~IN_BLACKLIST),os.path.basename(blacklist),len(variants_rel))); sys.stderr.flush()

if blacklists_pos is not None: # It used to be based on single locations, now whole regions
	"""
	The pybedtools function .from_dataframe does not work if the last field contains NaN (it sees one less column)
	Instead, we save it with NA represented by ' ', but this turns colums with such values into strings...
	"""
	for blacklist_pos in blacklists_pos:
		BLACKLIST_BED = pybedtools.BedTool(blacklist_pos)    
		LEN_START = len(variants_rel)
	variants_bedtool = pybedtools.BedTool.from_dataframe(variants_rel)
	# Make sure we do not remove variants labeled as "keep"
	variants_keep_bedtool = pybedtools.BedTool.from_dataframe(variants_rel[variants_rel['Keep']==True])
	blacklist_minus_keep = BLACKLIST_BED.intersect(variants_keep_bedtool,v=True)

	variants_in_regions = variants_bedtool.intersect(blacklist_minus_keep,v=True)
	variants_rel = variants_in_regions.to_dataframe(names=variants_rel.columns.tolist(),na_values='.',dtype=DTYPES)   

	sys.stderr.write("- {}/{} variants have been removed because positions were in blacklist {}. {} remain.\n".format(
	LEN_START-len(variants_rel),LEN_START,os.path.basename(blacklist_pos),len(variants_rel))); sys.stderr.flush()     

#--------------------------------------------------------------------------------------------------
# Additional filters to remove false positives
#---------------------------------------------------------------------------------------------
    
if repeats is not None: # Remove only repeats for indels. IMPORTANT: We could use len > 1 to include long substitutions
    
    REPEATS_BED = pybedtools.BedTool(repeats)    
    LEN_START = len(variants_rel)
    indels_bedtool = pybedtools.BedTool.from_dataframe(
            variants_rel[(variants_rel['indel'] == True) & (variants_rel['Keep'] == False)])    
    
    indels_not_repeats = indels_bedtool.intersect(REPEATS_BED,v=True)
    indels_not_repeats = indels_not_repeats.to_dataframe(names=variants_rel.columns.tolist(),na_values='.',dtype=DTYPES)

    variants_rel = pd.concat([variants_rel[(variants_rel['indel'] == False) | (variants_rel['Keep'] == True)],
    indels_not_repeats], ignore_index=True)
    
    sys.stderr.write("- {}/{} variants have been removed because they were indels falling in repeat regions. {} remain.\n".format(
    LEN_START-len(variants_rel),LEN_START,len(variants_rel))); sys.stderr.flush()  

# Filter based on SNP clusters
if MAX_PROXIMITY is not None:
    # Sort by coordinate        
    variants_rel = variants_rel.reindex(index=order_by_index(variants_rel.index,index_natsorted(
            zip(variants_rel['Chr'],variants_rel['Start'])))).reset_index(drop=True)
    # Calculate difference between positions (one row and the previous)
    variants_rel['diff'] = variants_rel['Start'].diff()
    # If absolute value (!) is smaller than max proximity, we take that row and the previous one
    CLOSE_VARIANTS = ( (variants_rel['diff'].abs() < MAX_PROXIMITY) | (variants_rel['diff'].abs() < MAX_PROXIMITY).shift(-1) )
    
    # https://stackoverflow.com/questions/45964740/python-pandas-cumsum-with-reset-everytime-there-is-a-0
    C = CLOSE_VARIANTS.cumsum()-CLOSE_VARIANTS.cumsum().where(~CLOSE_VARIANTS).ffill()
    C = C.sort_index(ascending=False).tolist()
    for i in range(1,len(C)):
        if C[i] < C[i-1] and C[i] != 0 and i != 0:
            C[i] = C[i-1]
    variants_rel['diff2'] = C[::-1]

    variants_rel = variants_rel[ (variants_rel['Keep'] == True) | (variants_rel['diff2'] <= MAX_CLUSTER) | (variants_rel['indel'] == True) ]
    variants_rel.drop(columns=['diff','diff2'],inplace=True)
    sys.stderr.write("- {}/{} variants have been removed based on proximity (distance <= {} and cluster >= {}). {} remain.\n".format(
         len(CLOSE_VARIANTS)-len(variants_rel),len(CLOSE_VARIANTS),MAX_PROXIMITY,MAX_CLUSTER,len(variants_rel))); sys.stderr.flush() 

# Filter SNVs in highly repetitive genomic regions - sequence identity of 95% to another region based on selfChain Link
if self_chain is not None:

    SELF_CHAIN = pd.read_csv(self_chain,sep='\t',names=['bin','score','tName','tSize','tStart',
                                                        'tEnd','qName','qSize','qStrand','qStart','qEnd','id','normScore'])
    SELF_CHAIN = SELF_CHAIN[~(SELF_CHAIN['tStart'] == SELF_CHAIN['qStart']) & (SELF_CHAIN['tName'] == SELF_CHAIN['qName'])]
    SELF_CHAIN = SELF_CHAIN[SELF_CHAIN['normScore'] >= 95]        
    SELF_CHAIN = pd.concat([SELF_CHAIN[['tName','tStart','tEnd']].rename(columns={'tName':'qName','tStart':'qStart','tEnd':'qEnd'}),
				 SELF_CHAIN[['qName','qStart','qEnd']]]).drop_duplicates()

    SELF_CHAIN_BED = pybedtools.BedTool.from_dataframe(SELF_CHAIN)    
    LEN_START = len(variants_rel)
    variants_bedtool = pybedtools.BedTool.from_dataframe(variants_rel)
    
    variants_in_regions = variants_bedtool.intersect(SELF_CHAIN_BED,v=True)
    variants_rel = variants_in_regions.to_dataframe(names=variants_rel.columns.tolist(),na_values='.',dtype=DTYPES)
    
    sys.stderr.write("- {}/{} variants have been removed because positions were in highly repetitive regions (selfChain). {} remain.\n".format(
    LEN_START-len(variants_rel),LEN_START,len(variants_rel))); sys.stderr.flush()     

#--------------------------------------------------------------------------------------------------
#Additional statistics and kept variants addition
#--------------------------------------------------------------------------------------------------

## STATS PLUS: Consider also applying to indels? Maybe NOT clipped though...
    
if stats_plus:
    
    if 'suboptimalAlt' in variants_rel.columns and 'mismatchesTooHighAlt' in variants_rel.columns:  
         # Check rna_vaf_pipeline for suboptimal condition  
        BOGUS_SNV = ( (variants_rel['indel'] == False) & ( ( (variants_rel['alternativesAlt'] >= MAX_AA_RATIO) &\
                        (variants_rel['suboptimalAlt'] >= MAX_SMALL_RATIO) ) | ( (variants_rel['alternatives'] >= MAX_AA_RATIO) &\
                        (variants_rel['suboptimal'] >= MAX_SMALL_RATIO) ) | (variants_rel['mismatchesTooHighAlt'] >= MAX_NON_SNP) |\
                        (variants_rel['clippedAlt'] >= MAX_CLIPPED_ALT) | (variants_rel['genotypingError'] >= MAX_GENOTYPE_ERROR) ) )
        BOGUS_INDEL = ( (variants_rel['indel'] == True) & ( ( (variants_rel['alternativesAlt'] >= MAX_AA_RATIO*2) &\
                        (variants_rel['suboptimalAlt'] >= MAX_SMALL_RATIO) ) | ( (variants_rel['alternatives'] >= MAX_AA_RATIO*2) &\
                        (variants_rel['suboptimal'] >= MAX_SMALL_RATIO) ) | (variants_rel['mismatchesTooHigh'] >= MAX_NON_SNP*2) ) )
        BOGUS_LOCI = (BOGUS_SNV | BOGUS_INDEL)
        variants_rel = variants_rel[ ~BOGUS_LOCI | (variants_rel['Keep'] == True) ]     
        
        sys.stderr.write("- {}/{} variants have been removed by additional statistics. {} remain.\n".format(
        sum(BOGUS_LOCI),len(variants_rel)+sum(BOGUS_LOCI),len(variants_rel))); sys.stderr.flush()     
        
    else:

        sys.stderr.write('WARNING: Additional BAM statistics not detected.\n') ; sys.stderr.flush()

elif 'suboptimalAlt' in variants_rel.columns and 'mismatchesTooHighAlt' in variants_rel.columns:

    sys.stderr.write('WARNING: Additional BAM statistics detected, but stats_plus not specified\n'); sys.stderr.flush() 

elif 'Alternative_alignment_var' in variants_rel.columns and 'NonSNP_variants_var' in variants_rel.columns:

    sys.stderr.write('WARNING: Mathijs additional BAM statistics detected. Please run another version of this script that can handle them') 

## FILTER BASED ON PREDICTION

def discretize(series,bins,labels):

     blank_rows = ( series.isin(['','.',' ']) | series.isnull() )

     dseries = series.mask(~blank_rows,pd.cut(series[~blank_rows].astype(float),\
bins=bins,include_lowest=True,labels=labels))

     return(dseries)

if somatic:

     # Discretize as suggested in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6386394
     variants_rel.insert(variants_rel.columns.get_loc('VEST3_score')+1,'VEST3_pred',\
	discretize(variants_rel['VEST3_score'],bins=[0,0.5,1],labels=['B','P']))
     variants_rel.insert(variants_rel.columns.get_loc('CADD_phred')+1,'CADD_pred',\
	discretize(variants_rel['CADD_phred'],bins=[0,20,99],labels=['B','P']))

     predictions = { 'LJB2_PP2_HDIV_Pred': ['D','P'], 'LJB2_PolyPhen2_HVAR_Pred': ['D','P'], 'LJB2_LRT_Pred': ['D'], \
      'LJB2_MutationTaster_Pred':['D','A'],'LJB_MutationAssessor_Pred':['medium','high'],'LRT_pred':['D'],'SIFT_pred':['D'],
      'Polyphen2_HDIV_pred':['D','P'],'Polyphen2_HVAR_pred':['D','P'],'MutationTaster_pred':['D','A'],'MutationAssessor_pred':['M','H'],
      'FATHMM_pred':['D'],'RadialSVM_pred':['D'],'LR_pred':['D'],'CADD_pred':['P'],'VEST3_pred':['P']} # For MutatinAssessor, high and medium FUNCTIONAL 

     deleterious = pd.Series([0]*len(variants_rel),index=variants_rel.index)
     known = pd.Series([0]*len(variants_rel),index=variants_rel.index)

     for column_prediction in variants_rel.columns[variants_rel.columns.str.contains('_Pred|_pred')]:

         deleterious = deleterious + variants_rel[column_prediction].isin(predictions[column_prediction])
         # Known have NO columns (predictions) with . / ' ' / NaN
         known = known + ~(variants_rel[column_prediction].isin(['','.',' ']) | variants_rel[column_prediction].isnull() )

     cancer_hematology = variants_rel[column_cosmic].str.extract(
             pat='(\d+)\(haematopoietic_and_lymphoid_tissue\)',expand=False).astype(float).fillna(0)
     # cancer = variants_rel[column_cosmic[0]].str.contains(pat='COS') -- any cancer

     # Majority voting: 50% of the prediction tools say it's deleterious
     # We also keep unknown (known -== 0) because we cannot exclude the possibily they are harmful
     # We remove FLAGGED cosmic variants as long as they are germline and only rare in COSMIC
     SOMATIC_LOCI = ( (deleterious >= 0.5*known) | (cancer_hematology >= 2) | (known == 0) )
     if len(column_cosmic_flagged) > 0:
         SOMATIC_LOCI = ( SOMATIC_LOCI & ~( (variants_rel[column_cosmic_flagged[0]] == 'FLAGGED') &
                         (variants_rel['alt freq'].between(0.44,0.56)) & (cancer_hematology < 15) ) )

     SOMATIC_LOCI = ( SOMATIC_LOCI | (variants_rel['Keep'] == True) |\
                     (variants_rel['ExonicFunc.refGene'].str.contains('stopgain|frameshift|startloss')) )

     variants_rel = variants_rel[ SOMATIC_LOCI ]
     
     sys.stderr.write("- {}/{} variants have been removed by deleterious prediction. {} remain.\n".format(
        sum(~SOMATIC_LOCI),len(variants_rel)+sum(~SOMATIC_LOCI),sum(SOMATIC_LOCI))); sys.stderr.flush()

if 'CLNSIG' in variants_rel.columns:
     BENIGN = variants_rel['CLNSIG'].str.contains('Benign|benign')
	
     variants_rel = variants_rel[ ~BENIGN ]
	
     sys.stderr.write("- {}/{} variants have been removed because they are reported as benign in ClinVar. {} remain.\n".format(
	sum(BENIGN),len(BENIGN),len(variants_rel))); sys.stderr.flush()

else:

     sys.stderr.write("WARNING: ClinVar annotation is missing. Annotation will proceed"); sys.stderr.flush()

if 'wgEncodeBroadHmmK562HMM' in variants_rel.columns:

    LEN_START = len(variants_rel)    
    cancer_hematology = variants_rel[column_cosmic].str.extract(
             pat='(\d+)\(haematopoietic_and_lymphoid_tissue\)',expand=False).astype(float).fillna(0)
    HMM_REPETITIVE = (variants_rel['wgEncodeBroadHmmK562HMM'].str.contains('14_Repetitive|15_Repetitive'))
    variants_rel = variants_rel [ ~HMM_REPETITIVE | (cancer_hematology >= 5) ]

    sys.stderr.write("- {}/{} variants have been removed because they are in repetitive regions of HMM K562. {} remain.\n".format(
        LEN_START-len(variants_rel),LEN_START,len(variants_rel))); sys.stderr.flush()

#--------------------------------------------------------------------------------------------------
#Adjust name and output to file
#--------------------------------------------------------------------------------------------------

# variants_rel.columns = variants_rel.columns.str.replace('.1','') # We can keep them to allow refiltering

# Sort by coordinate        
variants_rel = variants_rel.reindex(index=order_by_index(variants_rel.index,index_natsorted(
        zip(variants_rel['Chr'],variants_rel['Start'])))).reset_index(drop=True)

if output is None:

    variants_rel.to_csv(path_or_buf=sys.stdout,sep='\t',index=False,float_format='%.6g') 
    
else:
    
    variants_rel.to_csv(path_or_buf=output,sep='\t',index=False,float_format='%.6g')    

sys.stderr.write('\nINFO: Variant filtering finished in {} seconds\n'.format(time.process_time() - start_time)); sys.stderr.flush() 
