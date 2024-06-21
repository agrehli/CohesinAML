#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 09:26:49 2019

@author: roger
"""

variants = "/data/bas/1_enhancer_analysis/exome/bam/variant_calling/1316_180720_A00383_0011_AHFL2NDMXX.txt"
pindel =  "/data/bas/1_enhancer_analysis/exome/bam/variant_calling/pindel/1316_180720_A00383_0011_AHFL2NDMXX_pindel.vcf.gz"

#------------------------------------------------------------------------------
# LIBRARIES AND CUSTOM FUNCTIONS
#------------------------------------------------------------------------------

import sys
import os
import argparse
import pandas as pd
import subprocess

class ProcessException(Exception):
    """
    https://gis.stackexchange.com/questions/260853/pyqgis-merge-raster-layers
    """
    pass

def execute(command):
    """ 
    Function to constantly print print subprocess while it is running AND get the error code. 
    NOTE: To remove stderr, redirect it to the terminal: FNULL = open(os.devnull, 'w'); Popen(stderr=FNULL)
    https://stackoverflow.com/questions/11269575/how-to-hide-output-of-subprocess-in-python-2-7
    Use bash instead of /bin/sh: executable='/bin/bash'
    """     
    
    try:
        from StringIO import StringIO
    except ImportError: # For python 3
        from io import StringIO
   
    #FNULL = open(os.devnull, 'w')
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE,executable='/bin/bash')   
   
    output2 = process.communicate()[0]
    exitCode = process.returncode

    if (exitCode == 0):
        return StringIO(output2.decode('utf-8'))
    else:
        raise ProcessException(command, exitCode, output2)        
        
#------------------------------------------------------------------------------
# PARSE ARGUMENTS
#------------------------------------------------------------------------------
        
# variants = "/data/leonie/3q_capture/bam/variant_calling/atypical/09H019074_150430_SN583_0139_AC5WGDACXX.txt"
# pindel = "/data/leonie/3q_capture/bam/variant_calling/pindel/09H019074_150430_SN583_0139_AC5WGDACXX_pindel.vcf.gz"

parser = argparse.ArgumentParser(description='Script to add pindel supporting reads to text files processed with Annovar')
parser.add_argument('-v','--variants',type=str,help="Annovar file to be annotated (txt)")
parser.add_argument('-p','--pindel',type=str,help="Pindel merged VCF file")
parser.add_argument('-r','--reference',type=str,help="Reference FASTA file (for bcftools norm)")
parser.add_argument('-o','--output',default='stdout',type=str,help="Name of the output file. By default, print to the stdout (redirect with >)")
parser.add_argument('-m','--mode',default='add',type=str,help="Should the Pindel counts replace depth/alt or be added separately? (replace,add)")

args = parser.parse_args()
variants = args.variants
pindel = args.pindel
output = args.output
reference = args.reference
mode = args.mode

if not os.path.isfile(reference):    
    sys.exit("File {} does not exist".format(reference))

#------------------------------------------------------------------------------
# LIBRARIES AND CUSTOM FUNCTIONS
#------------------------------------------------------------------------------

sys.stderr.write('Info 1. Read Annovar variant file\n')
df = pd.read_csv(variants,sep='\t')

# Convert to Annovar format to ensure we have the same positions 
sys.stderr.write('Info 2. Normalize and process Pindel VCF with Annovar\n')

convert_annovar = "/tools/annovar/convert2annovar.pl --format vcf4old --includeinfo <(bcftools norm -m-both -f {} {})".format(reference,pindel)
out = execute(convert_annovar) 
vcf = pd.read_csv(out,sep='\t',names=['Chr','Start','End','Ref','Alt','chr2','pos2','id','ref2','alt2','qual','filter','info','format','sample'])

sys.stderr.write('Info 3. Merge variant data with Pindel results\n')

# Split the counts column and keep duplicates with largest pindel alt
vcf[['genotype','counts']] = vcf['sample'].str.split(":",expand=True)
vcf[['pindel_ref','pindel_alt']] = vcf['counts'].str.split(',',expand=True).apply(pd.to_numeric,downcast='integer')
vcf['pindel_depth'] = vcf['pindel_ref'] + vcf['pindel_alt']
vcf['pindel_alt_freq'] = vcf['pindel_alt'] / vcf['pindel_depth']

vcf2 = vcf[['Chr','Start','End','Ref','Alt','pindel_depth','pindel_alt','pindel_alt_freq']]
# https://stackoverflow.com/questions/15705630/python-getting-the-row-which-has-the-max-value-in-groups-using-groupby
vcf2 = vcf2.sort_values('pindel_alt',ascending=False).drop_duplicates(['Chr','Start','End','Ref','Alt'])

# Merge Annovar variants with Pindel counts
merged = df.merge(vcf2,on=['Chr','Start','End','Ref','Alt'],how='left')

if mode == 'replace':
    
    merged.loc[(merged['alt depth'] == 0) & (merged['pindel_alt'] > 0),['depth','alt depth','alt freq']] =\
    merged.loc[(merged['alt depth'] == 0) & (merged['pindel_alt'] > 0),['pindel_depth','pindel_alt','pindel_alt_freq']]

if output == 'stdout':
    try: # Prevent broken pipe error with pipelines in bash
        merged.to_csv(sys.stdout,sep='\t',index=False,float_format='%.6g')  
    except IOError:  # stdout is closed, no point in continuing   
        try:
            sys.stdout.close()  # Attempt to close them explicitly
        except IOError:
            pass
else:
    merged.to_csv(output,sep='\t',index=False,float_format='%.6g')  
