import os
import scanpy as sc
import pandas as pd
import numpy as np
import glob
import re
import anndata as ad
from joblib import parallel_backend
import warnings
import logging

sc.settings.n_jobs = 32
warnings.filterwarnings("ignore")

pwd = os.getcwd()

def build_effect_size_anndata(
    results_dir,
    glob_pattern,
    file_pattern,
    test,
    adata,
    subclass,
    celltype,
    blacklisted_genes=[],
    filter_genes=True,
    normalize=True,
    effect_size_cutoff=10,
):
    results_files = sorted(glob.glob(os.path.join(results_dir, glob_pattern)))
    subclasses = [re.sub("([^_]+)_(.*)(_[0-9]{1,2})?(-SEAAD)?$", "\\1", os.path.basename(i).replace(file_pattern, "").replace("Lamp5_Lhx6", "Lamp5 Lhx6")) for i in results_files]
    supertypes = [re.sub("([^_]+)_(.*)$", "\\2", os.path.basename(i).replace(file_pattern, "").replace("Lamp5_Lhx6", "Lamp5 Lhx6")) for i in results_files]
    supertypes = [i.replace("Lamp5 Lhx6", "Lamp5_Lhx6").replace("L2 3 IT", "L2/3 IT").replace("L5 6 NP", "L5/6 NP") for i in supertypes]
    effect_sizes = pd.DataFrame(np.zeros((len(adata.obs[celltype].cat.categories), adata.shape[1])), columns=adata.var_names, index=adata.obs[celltype].cat.categories)
    pvalues = pd.DataFrame(np.ones((len(adata.obs[celltype].cat.categories), adata.shape[1])), columns=adata.var_names, index=adata.obs[celltype].cat.categories)
    std_errors = pd.DataFrame(np.ones((len(adata.obs[celltype].cat.categories), adata.shape[1])), columns=adata.var_names, index=adata.obs[celltype].cat.categories)
    
    for i,j in enumerate(results_files):
        print(os.path.basename(j))
        results = pd.read_csv(j, index_col=0)
        effect_sizes.loc[supertypes[i], results.index] = results.loc[:, "logFC_" + test]
        pvalues.loc[supertypes[i], results.index] = results.loc[:, "p_" + test]
        std_errors.loc[supertypes[i], results.index] = results.loc[:, "se_" + test]
    
    effect_sizes = ad.AnnData(X=effect_sizes)
    subclasses = adata.obs.loc[:, [subclass, celltype]].drop_duplicates().sort_values(by=celltype).loc[:, subclass].to_list()
    effect_sizes.obs[subclass] = subclasses
    effect_sizes.obs[subclass] = effect_sizes.obs[subclass].astype("category")
    
    pvalues = ad.AnnData(X=pvalues)
    pvalues.obs[subclass] = subclasses
    std_errors = ad.AnnData(X=std_errors)
    std_errors.obs[subclass] = subclasses

    effect_sizes = effect_sizes.T
    pvalues = pvalues.T
    std_errors = std_errors.T

    effect_sizes.var["Class"] = "Non-neuronal and non-neural"
    effect_sizes.var.loc[[i in ["Sst", "Sst Chodl", "Pvalb", "Chandelier", "Vip", "Sncg", "Pax6", "Lamp5", "Lamp5 Lhx6"] for i in effect_sizes.var["Subclass"]], "Class"] = "Neuronal: GABAergic"
    effect_sizes.var.loc[[i in ["L2/3 IT", "L4 IT", "L5 IT", "L6 IT", "L6 IT Car3", "L5 ET", "L5/6 NP", "L6 CT", "L6b"] for i in effect_sizes.var["Subclass"]], "Class"] = "Neuronal: Glutamatergic"

    pvalues.var["Class"] = "Non-neuronal and non-neural"
    pvalues.var.loc[[i in ["Sst", "Sst Chodl", "Pvalb", "Chandelier", "Vip", "Sncg", "Pax6", "Lamp5", "Lamp5 Lhx6"] for i in pvalues.var["Subclass"]], "Class"] = "Neuronal: GABAergic"
    pvalues.var.loc[[i in ["L2/3 IT", "L4 IT", "L5 IT", "L6 IT", "L6 IT Car3", "L5 ET", "L5/6 NP", "L6 CT", "L6b"] for i in pvalues.var["Subclass"]], "Class"] = "Neuronal: Glutamatergic"

    std_errors.var["Class"] = "Non-neuronal and non-neural"
    std_errors.var.loc[[i in ["Sst", "Sst Chodl", "Pvalb", "Chandelier", "Vip", "Sncg", "Pax6", "Lamp5", "Lamp5 Lhx6"] for i in std_errors.var["Subclass"]], "Class"] = "Neuronal: GABAergic"
    std_errors.var.loc[[i in ["L2/3 IT", "L4 IT", "L5 IT", "L6 IT", "L6 IT Car3", "L5 ET", "L5/6 NP", "L6 CT", "L6b"] for i in std_errors.var["Subclass"]], "Class"] = "Neuronal: Glutamatergic"


    effect_sizes.var["Subclass"] = effect_sizes.var["Subclass"].cat.reorder_categories(
        adata.obs["Subclass"].cat.categories,
    )
    effect_sizes = effect_sizes[:, adata.obs["Supertype"].cat.categories].copy()
    effect_sizes.var["Supertype"] = effect_sizes.var.index.copy()

    if normalize == True:
        pvalues.layers["pvalues"] = pvalues.X.copy()
        pvalues.X = -np.log10(pvalues.X) * np.sign(effect_sizes.X)
        
        effect_sizes.X = np.nan_to_num(effect_sizes.X)
        effect_sizes.X[effect_sizes.X > effect_size_cutoff] = 0
        effect_sizes.X[effect_sizes.X < -1 * effect_size_cutoff] = 0
        effect_sizes.layers["effect_sizes"] = effect_sizes.X.copy()
        effect_sizes.X = effect_sizes.X / std_errors.X
        effect_sizes.X = np.nan_to_num(effect_sizes.X)

    if len(blacklisted_genes) > 0:
        effect_sizes = effect_sizes[[i not in blacklisted_genes for i in effect_sizes.obs_names], :].copy()
        pvalues = pvalues[[i not in blacklisted_genes for i in pvalues.obs_names], :].copy()
        std_errors = std_errors[[i not in blacklisted_genes for i in std_errors.obs_names], :].copy()
        
    return effect_sizes, pvalues, std_errors

def construct_gene_graph(
    mean_expression,
    fraction_expressed,
    effect_sizes_early,
    effect_sizes_late,
    prefix,
    scale_expression=True,
    scale_effect_sizes=True,
    max_exp_z=2,
    min_mean_exp=0.035,
    min_frac_exp=0.04,
    n_neighbors=15,
    aggregate_metrics=False,
    vectors=None
):

    genes_to_keep = mean_expression.obs_names[
        (mean_expression.X.max(axis=1) > min_mean_exp) &
        (fraction_expressed.X.max(axis=1) > min_frac_exp)
    ]
    
    if scale_expression == True:
        mean_expression = mean_expression.T
        sc.pp.scale(mean_expression)
        mean_expression = mean_expression.T

    mean_expression.X[mean_expression.X > max_exp_z] = max_exp_z
    mean_expression.X[mean_expression.X < -1 * max_exp_z] = -1 * max_exp_z
    mean_expression = mean_expression[effect_sizes_early.obs_names, :].copy()
    
    if aggregate_metrics == True:
        for i,j in vectors.items():
            effect_sizes_early.obs["effect_sizes_early_" + str(i)] = effect_sizes_early[:, j == 1].X.mean(axis=1).tolist()
            effect_sizes_late.obs["effect_sizes_late_" + str(i)] = effect_sizes_late[:, j == 1].X.mean(axis=1).tolist()
            mean_expression.obs["mean_expression_" + str(i)] = mean_expression[:, j == 1].X.mean(axis=1).tolist()

    effect_sizes_all = ad.concat([effect_sizes_early, effect_sizes_late, mean_expression], axis=1, keys=["early", "late", "mean"], index_unique="_", merge="first")
    effect_sizes_all = effect_sizes_all[genes_to_keep].copy()
    
    with parallel_backend('threading', n_jobs=32):
        sc.pp.neighbors(effect_sizes_all, use_rep="X", n_neighbors=n_neighbors)
        sc.tl.umap(effect_sizes_all, min_dist=0.3)
                    
    return effect_sizes_all