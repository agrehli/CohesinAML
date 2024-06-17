#/usr/local/bin/pyEnv_sccoda-0.1.8.sh"
#/misc/software/ngs/sccoda/scCODA-0.1.8/pyEnv39/bin/python

#scCODA analysis of cluster disribution in STAG2KO HSPC experiment
#by Alexander Fischer
# July 2022

#import libraries
import sys
# Setup
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import pickle as pkl
import anndata as ad
import arviz as az

#import sccoda
from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz
#to use other model also import
from sccoda.model import other_models as om
import numpy as np
import patsy as pt
import matplotlib.pyplot as plt

# Data preparation
cell_counts = pd.read_csv("/misc/data/analysis/project_cohesin/CD34/scRNAseq/custerdistribution.summary.csv")
cell_counts.rename(columns={'Unnamed: 0':'sample'}, inplace=True)
#add condition and donor as extra columns
cell_counts[["Condition","donor"]]=cell_counts["sample"].str.split("_", expand=True)
print(cell_counts)


# Convert data to anndata object
data_all = dat.from_pandas(cell_counts, covariate_columns=["sample","Condition","donor"])
print(data_all.obs)

#grouped boxplot split by condition
viz.boxplots(data_all, feature_name="Condition",add_dots=True,cmap="Reds")
figure = plt.gcf() ####gets the current figure
figure.set_size_inches(10, 8)
plt.savefig('/misc/data/analysis/project_cohesin/CD34/scRNAseq/scCODA_analysis/grade_boxplots.png', dpi = 300)

# Stacked barplot for each sample
viz.stacked_barplot(data_all, feature_name="sample")
figure = plt.gcf()
figure.set_size_inches(25, 10)
plt.savefig('/misc/data/analysis/project_cohesin/CD34/scRNAseq/scCODA_analysis/stacked_barplot.png', dpi = 200)

# Stacked barplot by condition
viz.stacked_barplot(data_all, feature_name="Condition")
figure = plt.gcf()
figure.set_size_inches(10, 10)
plt.savefig('/misc/data/analysis/project_cohesin/CD34/scRNAseq/scCODA_analysis/stacked_barplot.condition.png', dpi = 200)

#Finding a reference cell type: presence (share of non-zero samples) over all samples for each cell type versus its dispersion in relative abundance
viz.rel_abundance_dispersion_plot(
    data=data_all,
    abundant_threshold=0.9
)
figure = plt.gcf()
figure.set_size_inches(10, 10)
plt.savefig('/misc/data/analysis/project_cohesin/CD34/scRNAseq/scCODA_analysis/rel_abundance_dispersion_plot.png', dpi = 200)

####BCELL has lowest dispersion could be used as reference

#Statistical model with reference cell type (BCELL is also selected if set to "automatic")
model_all = mod.CompositionalAnalysis(data_all, formula="Condition + donor", reference_cell_type="BCELL")
#calculate results (HMC sampling)
sim_results = model_all.sample_hmc()
#result summary
sim_results.summary()

#show only credible effects (default FDR 0.05)
sim_results.set_fdr(est_fdr=0.05) ##HSC only bocomes sig if fdr = 0.25 !
sim_results.summary()
print(sim_results.credible_effects()) #only NEUT and MONO


#trace plot
az.plot_trace(
    sim_results,
    divergences=False,
    var_names=["alpha", "beta"],
    coords={"cell_type": sim_results.posterior.coords["cell_type_nb"]},
)
plt.savefig('/misc/data/analysis/project_cohesin/CD34/scRNAseq/scCODA_analysis/trace.png')


# saving
##dataframes
sim_results.effect_df.to_csv("/misc/data/analysis/project_cohesin/CD34/scRNAseq/scCODA_analysis/sccoda_general_effect.csv", sep='\t')
sim_results.to_dataframe().to_csv("/misc/data/analysis/project_cohesin/CD34/scRNAseq/scCODA_analysis/sccoda_general_all.csv", sep='\t')
##object 
path = "/misc/data/analysis/project_cohesin/CD34/scRNAseq/scCODA_analysis/sim_results"
sim_results.save(path)