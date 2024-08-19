import os
import numpy as np
import anndata
import re
import pandas as pd
import matplotlib.pyplot as plt
import csv
import scanpy as sc
import json
import hashlib
import scvi
import glob
import seaborn as sns
import copy
import random
import scipy.sparse as sp_sparse
import scipy.stats as sp_stats
import transcriptomic_clustering as tc
from iterative_scANVI import *
from datetime import datetime
from joblib import parallel_backend
from joblib import Parallel, delayed
from igraph import *
import warnings

warnings.filterwarnings("ignore")

'''
Integrates and predicts labels for a query dataset iteratively using scVI and scANVI (Xu et al 2021 Mol Syst Biol)

Input arguements:
adata_query: (AnnData) Object with unknown cells

adata_ref: (AnnData) Object with reference labels

labels_keys: (list) Ordered list of labels to iteratively predict and split (e.g. ["class", "subclass", "cluster"])

output_dir: (str) Location to store trained models and final results CSV

**kwargs: (dict) Passed to several functions, details below:

    layer: (None or str, default None) None if unnormalized counts are in AnnData.X, else a str where they are stored in AnnData.layers
    
    batch_key: (None or str, default None) Name of the batch variable pass to scVI and scANVI (e.g. "donor_name")
    
    categorical_covariate_keys: (list) List of categorical covariates to pass to scVI and scANVI (e.g. ["donor_name", "sex"])
    
    continuous_covariate_keys: (list) List of continuous covariates to pass to scVI and scANVI (e.g. ["n_genes"])
    
    use_hvg: (bool, default True) Whether to calculate and include highly variable genes from the reference dataset to scVI and scANVI
    
    use_de: (bool, default True) Whether to calculate and include differentially genes from the reference dataset to scVI and scANVI
    
    n_top_genes: (int or list, default 2000) Number of highly variable genes to pass in each iteration. Not used if use_hvg is False
    
    n_downsample_ref: (int or list, default 1000) Number of cells to downsample reference groups too in each iteration when 
    calculcating differentially expressed genes
    
    n_ref_genes: (int or list, default 500): Number of differentially expressed genes per group to pass in each iteration from the 
    reference dataset to scVI and scANVI
    
    max_epochs_scVI: (int or list, default 200) Number of epochs to train scVI for in each iteration
    
    max_epochs_scANVI: (int or list, default 20) Number of epochs to train scANVI for in each iteration 
    
    min_accuracy: (float, default 0.85) Minimum accuracy for label self-projection needed to avoid an accuracy flag in cleanup (See output).
    
    plot_confusion: (bool, default True) Whether to plot the confusion matrix after scANVI label prediction
    
    scVI_model_args: (dict, default {"n_layer": 2}) kwargs passed to scvi.model.SCVI
    
    scANVI_model_args: (dict, default {}) kwargs passed to scvi.model.SCANVI.from_scvi_model

Outputs:
scVI models in output_dir/scVI_models

scANVI models in output_dir/scANVI_models

CSV named "iterative_scANVI_results.<date>.csv" in output_dir with cells as rows and the following columns:

    <Label>_stash (str) Original labels from reference and query datasets
    
    <Label> (str) Original labels from reference and query datasets, coerced to "Unknown" if its predicted label disagreed with ground truth
    
    <Label>_scANVI: (str) Predicted type for each level of the hierarchy (if [labels_keys] was ["class", "subclass"] class_scANVI
    and subclass_scANVI would exist)
    
    <Label>_conf_scANVI: (float) Probability of the prediction
    
    <All Possible Label Values>: (float) Probability a cell is that label (e.g. glia, excitatory, inhibitory or Astro, Endo, VLMC, etc)
    
    cleanup: ("Unknown" or "Reference" or "Accuracy" or "Confidence") Column indicated whether a cell was flagged as reference, 
    having a label with a low prediction accuracy, or having a confidence signficantly different from others in its group
    
Example Usage:
iterative_scANVI_kwargs = {
    "categorical_covariate_keys": ["donor_name"]
    "continuous_covariate_keys": ["n_genes"]
    "n_top_genes": [5000, 2000, 2000]
}

iterative_scANVI(
    adata_query, 
    adata_ref,
    labels_keys=["class", "subclass", "cluster"]
    output_dir=os.path.join("scANVI_output"), 
    **iterative_scANVI_kwargs
)

'''

def iterative_scANVI(adata_query, adata_ref, labels_keys, output_dir, **kwargs):
    
    default_kwargs = {
        "layer": "UMIs",
        "batch_key": None,
        "categorical_covariate_keys": None,
        "continuous_covariate_keys": None,
        "use_hvg": True,
        "use_de": True,
        "n_top_genes": 2000,
        "n_downsample_ref": 1000,
        "n_ref_genes": 500,
        "user_genes": None,
        "max_epochs_scVI": 200,
        "max_epochs_scANVI": 20,
        "min_accuracy": 0.85,
        "plot_confusion": True,
        "scVI_model_args": {"n_layers": 2},
        "scANVI_model_args": {},
    }
    
    kwargs = {**default_kwargs, **kwargs}
    
    layer = kwargs["layer"]
    batch_key = kwargs["batch_key"]
    categorical_covariate_keys = kwargs["categorical_covariate_keys"]
    continuous_covariate_keys = kwargs["continuous_covariate_keys"]
    use_hvg = kwargs["use_hvg"]
    use_de = kwargs["use_de"]
    n_top_genes = kwargs["n_top_genes"]
    n_downsample_ref = kwargs["n_downsample_ref"]
    n_ref_genes = kwargs["n_ref_genes"]
    user_genes = kwargs["user_genes"]
    max_epochs_scVI = kwargs["max_epochs_scVI"]
    max_epochs_scANVI = kwargs["max_epochs_scANVI"]
    min_accuracy = kwargs["min_accuracy"]
    plot_confusion = kwargs["plot_confusion"]
    scVI_model_args = kwargs["scVI_model_args"]
    scANVI_model_args = kwargs["scANVI_model_args"]
    
    if isinstance(adata_query, anndata.AnnData) is False or isinstance(adata_ref, anndata.AnnData) is False:
        raise TypeError("One or more of the AnnData objects have an incorrect type,")

    if adata_query.shape[1] != adata_ref.shape[1]:
        common_labels = np.intersect1d(adata_query.var_names, adata_ref.var_names)
        adata_ref = adata_ref[:, common_labels].copy()
        adata_query = adata_query[:, common_labels].copy()
        warnings.warn("Reference and query AnnData objects have different shapes, using " + str(len(common_labels)) + " common labels. This may have a deterimental effect on model performance.")
    
    try:
        iter(labels_keys)
    except TypeError:
        raise TypeError("labels_keys must be an iterable.")
        
    if all([i in adata_ref.obs.columns for i in labels_keys]) is False:
        raise KeyError("One or more labels_keys do not exist in the reference AnnData object.")        
    
    if layer != None:
        if layer not in adata_query.layers.keys() or layer not in adata_ref.layers.keys():
            raise KeyError("Layer " + layer + "does not exist in both AnnData objects.")
    
    else:
        raise ValueError("You must provide a layer key that contains raw UMIs or counts. AnnData.X is assumed to be log-normalized expression data.")
    
    if batch_key is not None:
        if batch_key not in adata_query.obs.columns or batch_key not in adata_ref.obs.columns:
            raise KeyError("Batch key is not in both AnnData objects.")
    try:
        if all([i in adata_query.obs.columns for i in categorical_covariate_keys]) is False or all([i in adata_ref.obs.columns for i in categorical_covariate_keys]) is False:
            raise KeyError("One or more categorical covariates do not exist in both AnnData objects.")
    except:
        pass
    
    try:    
        if all([i in adata_query.obs.columns for i in continuous_covariate_keys]) is False or all([i in adata_ref.obs.columns for i in continuous_covariate_keys]) is False:
            warnings.warn("One or more continuous covariates do not exist in both AnnData objects. This can be safely ignored if they are calculated later.")
    except:
        pass
                          
    for i in [use_hvg, use_de, n_downsample_ref, n_ref_genes, user_genes, max_epochs_scVI, max_epochs_scANVI]:
        if isinstance(i, bool) or isinstance(i, str) or isinstance(i, int) or i is None:
            i = [i]
        if len(i) != 1 and len(i) != len(labels_keys):
            raise ValueError(str(i) + " should be either a single value used in each iteration or a list of values of equal length to labels_keys.")
               
    adata = adata_query.concatenate(adata_ref, index_unique=None)
    del adata_query
    
    
    for i,j in enumerate(labels_keys):
        get_model_genes_kwargs = {
            "groupby": j,
            "use_hvg": use_hvg[i] if isinstance(use_hvg, list) else use_hvg,
            "use_de": use_de[i] if isinstance(use_de, list) else use_de,
            "n_top_genes": n_top_genes[i] if isinstance(n_top_genes, list) else n_top_genes,
            "n_downsample_ref": n_downsample_ref[i] if isinstance(n_downsample_ref, list) else n_downsample_ref,
            "n_ref_genes": n_ref_genes[i] if isinstance(n_ref_genes, list) else n_ref_genes,
            "user_genes": user_genes[i] if isinstance(user_genes, list) else None
        }
        run_scVI_kwargs = {
            "layer": layer,
            "max_epochs_scVI": max_epochs_scVI[i] if isinstance(max_epochs_scVI, list) else max_epochs_scVI,
            "batch_key": batch_key,
            "categorical_covariate_keys": categorical_covariate_keys,
            "continuous_covariate_keys": continuous_covariate_keys,
            "scVI_model_args": scVI_model_args
        }
        run_scANVI_kwargs = {
            "layer": layer,
            "max_epochs_scANVI": max_epochs_scANVI[i] if isinstance(max_epochs_scANVI, list) else max_epochs_scANVI,
            "batch_key": batch_key,
            "categorical_covariate_keys": categorical_covariate_keys,
            "continuous_covariate_keys": continuous_covariate_keys,
            "labels_key": j,
            "scANVI_model_args": scANVI_model_args
        }
                
        adata_ref.obs[j] = adata_ref.obs[j].astype("category")
        adata.obs[j] = adata.obs[j].astype('object')
        adata.obs.loc[adata.obs[j].isna(), j] = "Unknown"
        adata.obs[j] = adata.obs[j].astype('category')
        adata.obs[j + "_stash"] = adata.obs[j].copy()

        # To do: Unify the if i == 0 and else: blocks.  
        if i == 0:
            print(str(datetime.now()) + " --- Training model predicting " + j + " on entire dataset")
            model_name = hashlib.md5(str(json.dumps({**get_model_genes_kwargs, **run_scVI_kwargs})).replace("/", " ").encode()).hexdigest()
            label_model_name = hashlib.md5(str(json.dumps({**get_model_genes_kwargs, **run_scANVI_kwargs})).replace("/", " ").encode()).hexdigest()

            if os.path.exists(os.path.join(output_dir, "scANVI_models", label_model_name)) is False:
                
                if os.path.exists(os.path.join(output_dir, "scVI_models", model_name)) is False:
                    if user_genes is None:
                        markers = get_model_genes(adata_ref, **get_model_genes_kwargs)
                    else:
                        markers = user_genes[i]
                    model = run_scVI(adata[:, markers], **run_scVI_kwargs)
                    model.save(os.path.join(output_dir, "scVI_models", model_name))

                else:
                    markers = pd.read_csv(os.path.join(output_dir, "scVI_models", model_name, "var_names.csv"), header=None)
                    markers = markers[0].to_list()
                    model = scvi.model.SCVI.load(os.path.join(output_dir, "scVI_models", model_name), adata[:, markers])
                
                label_model, probabilities = run_scANVI(adata[:, markers], model=model, **run_scANVI_kwargs)            
                label_model.save(os.path.join(output_dir, "scANVI_models", label_model_name))
                probabilities.to_csv(os.path.join(output_dir, "scANVI_models", label_model_name, "probabilities.csv"))
                    
            else:
                probabilities = pd.read_csv(os.path.join(output_dir, "scANVI_models", label_model_name, "probabilities.csv"), index_col=0)
            probabilities.index = [str(l) for l in probabilities.index]
            probabilities.dropna(axis=1, how='all', inplace=True)
            probabilities.drop([l for l in probabilities.columns if l.startswith("_")], axis=1, inplace=True)
            tmp = pd.merge(adata.obs, probabilities, how="left", left_index=True, right_index=True)
            for l in [m for m in tmp.columns if m.endswith("_y")]:
                l = l.replace("_y", "")
                tmp[l + "_x"] = tmp[l + "_x"].astype("object")
                tmp[l + "_y"] = tmp[l + "_y"].astype("object")
                tmp[l] = tmp[l + "_y"].fillna(tmp[l + "_x"])
                tmp[l] = tmp[l].astype("category")
                tmp.drop([l + "_y", l + "_x"], axis=1, inplace=True)
            adata.obs = tmp.copy()

            conf_mat = adata.obs.groupby([j, j + "_scANVI"]).size().unstack(fill_value=0)

            if plot_confusion is True:
                display(conf_mat)

            conf_mat = conf_mat.div(conf_mat.sum(axis=1), axis=0)
            conf_mat[conf_mat.isna()] = 0
            conf_mat = conf_mat.reindex(sorted(conf_mat.columns), axis=1)
            conf_mat = conf_mat.reindex(sorted(conf_mat.index), axis=0)

            for l in conf_mat.index:
                if l == "Unknown":
                    continue
                elif l not in conf_mat.columns:
                    print("WARNING: Label " + l + " fell below accruacy threshold " + str(min_accuracy) + " on reference cells. Label was not used.")
                elif conf_mat.loc[l,l] < min_accuracy:
                    print("WARNING: Label " + l + " fell below accruacy threshold " + str(min_accuracy) + " on reference cells. Accuracy=" + str(conf_mat.loc[l,l]))

            if plot_confusion is True:
                plt.figure(figsize=(10, 10))
                ax = plt.pcolor(conf_mat)
                ax = plt.xticks(np.arange(0.5, len(conf_mat.columns), 1), conf_mat.columns, rotation=90)
                ax = plt.yticks(np.arange(0.5, len(conf_mat.index), 1), conf_mat.index)
                plt.xlabel("Predicted")
                plt.ylabel("Observed")
                plt.colorbar()
                plt.show()
                
        else:
            for z,k in enumerate(adata_ref.obs[labels_keys[i - 1]].cat.categories):
                if k == "Unknown":
                    continue
                    
                print(str(datetime.now()) + " --- Training model predicting " + j + " on " + labels_keys[i - 1] + "_scANVI=" + k + " subsetted dataset")
                model_name = hashlib.md5(str(json.dumps({**{labels_keys[i - 1]: k}, **get_model_genes_kwargs, **run_scVI_kwargs})).replace("/", " ").encode()).hexdigest()
                label_model_name = hashlib.md5(str(json.dumps({**{labels_keys[i - 1]: k}, **get_model_genes_kwargs, **run_scANVI_kwargs})).replace("/", " ").encode()).hexdigest()
                
                cells = adata.obs[labels_keys[i - 1] + "_scANVI"] == k
                if any(cells) is False:
                    continue
                    
                adata.obs[labels_keys[i - 1]] = adata.obs[labels_keys[i - 1]].astype("object")
                adata.obs[labels_keys[i - 1] + "_scANVI"] = adata.obs[labels_keys[i - 1] + "_scANVI"].astype("object")
                adata.obs.loc[(adata.obs[labels_keys[i - 1]] != adata.obs[labels_keys[i - 1] + "_scANVI"]) & (cells), j] = "Unknown"
                adata.obs[labels_keys[i - 1]] = adata.obs[labels_keys[i - 1]].astype("category")
                adata.obs[labels_keys[i - 1] + "_scANVI"] = adata.obs[labels_keys[i - 1] + "_scANVI"].astype("category")
                
                if os.path.exists(os.path.join(output_dir, "scANVI_models", label_model_name)) is False:

                    if os.path.exists(os.path.join(output_dir, "scVI_models", model_name)) is False:
                        ref_cells = adata_ref.obs[labels_keys[i - 1]] == k

                        if user_genes is None:
                            markers = get_model_genes(adata_ref[ref_cells], **get_model_genes_kwargs)
                        else:
                            if isinstance(user_genes[i], dict):
                                markers = user_genes[i][k]
                            else:
                                markers = user_genes[i]

                        model = run_scVI(adata[cells, markers], **run_scVI_kwargs)
                        model.save(os.path.join(output_dir, "scVI_models", model_name))

                    else:
                        markers = pd.read_csv(os.path.join(output_dir, "scVI_models", model_name, "var_names.csv"), header=None)
                        markers = markers[0].to_list()
                        model = scvi.model.SCVI.load(os.path.join(output_dir, "scVI_models", model_name), adata[cells, markers])
                    
                    try:
                        label_model, probabilities = run_scANVI(adata[cells, markers], model=model, **run_scANVI_kwargs)
                        label_model.save(os.path.join(output_dir, "scANVI_models", label_model_name))
                        probabilities.to_csv(os.path.join(output_dir, "scANVI_models", label_model_name, "probabilities.csv"))
                    except IndexError:
                        continue

                else:
                    probabilities = pd.read_csv(os.path.join(output_dir, "scANVI_models", label_model_name, "probabilities.csv"), index_col=0)
                probabilities.index = [str(l) for l in probabilities.index]
                probabilities.dropna(axis=1, how='all', inplace=True)
                probabilities.drop([l for l in probabilities.columns if l.startswith("_")], axis=1, inplace=True)
                tmp = pd.merge(adata.obs, probabilities, how="left", left_index=True, right_index=True)
                for l in [m for m in tmp.columns if m.endswith("_y")]:
                    l = l.replace("_y", "")
                    tmp[l + "_x"] = tmp[l + "_x"].astype("object")
                    tmp[l + "_y"] = tmp[l + "_y"].astype("object")
                    tmp[l] = tmp[l + "_y"].fillna(tmp[l + "_x"])
                    tmp[l] = tmp[l].astype("category")
                    tmp.drop([l + "_y", l + "_x"], axis=1, inplace=True)
                adata.obs = tmp.copy()
                
                confs = [l for l in adata.obs.columns if l.endswith("_conf_scANVI")]
                
                for l in confs:
                    adata.obs[l] = adata.obs[l].astype('float')

                conf_mat = adata[cells].obs.groupby([j, j + "_scANVI"]).size().unstack(fill_value=0)

                if plot_confusion is True:
                    display(conf_mat)

                conf_mat = conf_mat.div(conf_mat.sum(axis=1), axis=0)
                conf_mat[conf_mat.isna()] = 0
                conf_mat = conf_mat.reindex(sorted(conf_mat.columns), axis=1)
                conf_mat = conf_mat.reindex(sorted(conf_mat.index), axis=0)
                
                for l in conf_mat.index:
                    if l == "Unknown":
                        continue
                    elif l not in conf_mat.columns:
                        print("WARNING: Label " + l + " fell below accruacy threshold " + str(min_accuracy) + " on reference cells. Label was not used.")
                    elif conf_mat.loc[l,l] < min_accuracy:
                        print("WARNING: Label " + l + " fell below accruacy threshold " + str(min_accuracy) + " on reference cells. Accuracy=" + str(conf_mat.loc[l,l]))

                if plot_confusion is True:
                    plt.figure(figsize=(10, 10))
                    ax = plt.pcolor(conf_mat)
                    ax = plt.xticks(np.arange(0.5, len(conf_mat.columns), 1), conf_mat.columns, rotation=90)
                    ax = plt.yticks(np.arange(0.5, len(conf_mat.index), 1), conf_mat.index)
                    plt.xlabel("Predicted")
                    plt.ylabel("Observed")
                    plt.colorbar()
                    plt.show()
                
    tmp = pd.concat([adata.obs.iloc[:, int(np.where(adata.obs.columns == labels_keys[0] + "_scANVI")[0]):], adata.obs.iloc[:, np.where([i in labels_keys for i in adata.obs.columns])[0]]], axis=1)
    tmp.to_csv(os.path.join(output_dir, "iterative_scANVI_results." + str(datetime.date(datetime.now())) + ".csv"))

'''
Creates a list of genes that are either highly variable or differentially expressed in the given labels
Called from iterative_scANVI.

Input arguements:
adata_ref: (AnnData) object with reference labels
**kwargs: (dict) Passed to several functions, details below:
    
    use_hvg: (bool, default True) Whether to calculate and include highly variable genes from the reference dataset to scVI and scANVI
    
    use_de: (bool, default True) Whether to calculate and include differentially genes from the reference dataset to scVI and scANVI
    
    n_top_genes: (int or list, default 2000) Number of highly variable genes to pass in each iteration. Not used if use_hvg is False
    
    n_downsample_ref: (int or list, default 1000) Number of cells to downsample reference groups too in each iteration when 
    calculcating differentially expressed genes
    
    n_ref_genes: (int or list, default 500): Number of differentially expressed genes per group to pass in each iteration from the 
    reference dataset to scVI and scANVI

Outputs:
Returns list of unique markers from the procedure above
'''
                             
def get_model_genes(adata_ref, **kwargs):
    
    groupby = kwargs["groupby"]
    use_hvg = kwargs["use_hvg"]
    use_de = kwargs["use_de"]
    n_top_genes = kwargs["n_top_genes"]
    n_downsample_ref = kwargs["n_downsample_ref"]
    n_ref_genes = kwargs["n_ref_genes"]
    user_genes = kwargs["user_genes"]
    
    markers = []
    
    if "log1p" not in adata_ref.uns_keys():    
        sc.pp.normalize_total(adata_ref, target_sum=1e4)
        sc.pp.log1p(adata_ref)
    
    if use_hvg is True:
        try:
            sc.pp.highly_variable_genes(adata_ref, flavor="seurat_v3", n_top_genes=n_top_genes, layer="UMIs")
        except:
            sc.pp.highly_variable_genes(adata_ref, min_mean=1, min_disp=0.5)    
        markers = adata_ref.var[adata_ref.var.highly_variable == True].index.to_list()
        
    if use_de is True:
        ref_counts = adata_ref.obs[groupby].value_counts()
        adata_ref = adata_ref[~(adata_ref.obs[groupby].isin(ref_counts[ref_counts < 15].index))]
        
        if np.setdiff1d(adata_ref.obs[groupby].cat.categories, "Unknown").shape[0] > 1:
            
            cells = []
            for i in adata_ref.obs[groupby].cat.categories:
                tmp_cells = adata_ref[adata_ref.obs[groupby] == i].obs_names.to_list()
                
                if len(tmp_cells) > n_downsample_ref:
                    cells = cells + random.sample(tmp_cells, k=n_downsample_ref)

                else:
                    cells.extend(tmp_cells)

            adata_ref = adata_ref[cells]
            
            sc.tl.rank_genes_groups(adata_ref, method="wilcoxon", tie_correct=True, groupby=groupby, pts=True)

            result = adata_ref.uns['rank_genes_groups']
            groups = result['names'].dtype.names
            marker_genes = {group: pd.DataFrame({key: result[key][group] for key in ['names', 'pvals_adj', 'logfoldchanges']}) for group in groups}

            for group in groups:
                marker_genes[group]['pts'] = result['pts'][group][result['names'][group]].to_list()
                marker_genes[group]['pts_rest'] = result['pts_rest'][group][result['names'][group]].to_list()
                marker_genes[group].index = marker_genes[group].names
                marker_genes[group].drop(columns=['names'], inplace=True)
                tmp_genes = marker_genes[group].copy()
                tmp_genes = tmp_genes[tmp_genes.pvals_adj < 0.05]
                tmp_genes.sort_values(by="logfoldchanges", axis=0, inplace=True, ascending=False)
                markers.extend(tmp_genes.head(n_ref_genes).index.to_list())
        else:
            warnings.warn(groupby + " contains only one label. Differentially expressed genes were NOT included in the model")
    
    return np.unique(markers)

'''
Wrapper for scVI
Called from iterative_scANVI.

Input arguements:
adata: (AnnData) Merged AnnData object
**kwargs: (dict) Passed to several functions, details below:

    layer: (None or str, default None) None if unnormalized counts are in AnnData.X, else a str where they are stored in AnnData.layers
    
    categorical_covariate_keys: (list) List of categorical covariates to pass to scVI and scANVI (e.g. ["donor_name"])
    
    continuous_covariate_keys: (list) List of continuous covariates to pass to scVI and scANVI (e.g. ["n_genes"])
    
    max_epochs_scVI: (int, default 200) Number of epochs to train scVI
    
    scVI_model_args: (None or dict) kwargs passed to scvi.model.SCVI


Outputs:
Returns trained scVI model
'''

def run_scVI(adata, **kwargs):
    layer = kwargs["layer"]
    max_epochs_scVI = kwargs["max_epochs_scVI"]
    batch_key = kwargs["batch_key"]
    categorical_covariate_keys = kwargs["categorical_covariate_keys"]
    continuous_covariate_keys = kwargs["continuous_covariate_keys"]
    scVI_model_args = kwargs["scVI_model_args"]
        
    adata = adata.copy()
        
    scvi.model.SCVI.setup_anndata(
        adata,
        layer=layer,
        batch_key=batch_key,
        categorical_covariate_keys=categorical_covariate_keys,
        continuous_covariate_keys=continuous_covariate_keys
    )
    model = scvi.model.SCVI(adata, **scVI_model_args)
    model.train(max_epochs=max_epochs_scVI, early_stopping=True)
    
    return model

'''
Wrapper for scANVI
Called from iterative_scANVI.

Input arguements:
adata: (AnnData) Merged AnnData object
**kwargs: (dict) Passed to several functions, details below:

    layer: (None or str, default None) None if unnormalized counts are in AnnData.X, else a str where they are stored in AnnData.layers
    
    categorical_covariate_keys: (list) List of categorical covariates to pass to scVI and scANVI (e.g. ["donor_name"])
    
    continuous_covariate_keys: (list) List of continuous covariates to pass to scVI and scANVI (e.g. ["n_genes"])
    
    max_epochs_scANVI: (int, default 20) Number of epochs to train scANVI
    
    scANVI_model_args: (None or dict) kwargs passed to scvi.model.SCANVI.from_scvi_model


Outputs:
Returns tupple with trained scANVI model and label predictions/probabilities (pd.DataFrame)
'''
                              
def run_scANVI(adata, model, **kwargs):
    layer = kwargs["layer"]
    max_epochs_scANVI = kwargs["max_epochs_scANVI"]
    batch_key = kwargs["batch_key"]
    categorical_covariate_keys = kwargs["categorical_covariate_keys"]
    continuous_covariate_keys = kwargs["continuous_covariate_keys"]
    labels_key = kwargs["labels_key"]
    scANVI_model_args = kwargs["scANVI_model_args"]
        
    adata = adata.copy()
    
    scvi.model.SCANVI.setup_anndata(
        adata,
        layer=layer,
        batch_key=batch_key,
        categorical_covariate_keys=categorical_covariate_keys,
        continuous_covariate_keys=continuous_covariate_keys,
        labels_key=labels_key
    )
    label_model = scvi.model.SCANVI.from_scvi_model(model, "Unknown", adata=adata, **scANVI_model_args)
    label_model.train(max_epochs=max_epochs_scANVI, early_stopping=True)
    
    adata.obs[labels_key + "_scANVI"] = label_model.predict()
    adata.obs[labels_key + "_scANVI"] = adata.obs[labels_key + "_scANVI"].astype("category")

    probabilities = label_model.predict(soft=True)

    tmp = pd.merge(adata.obs, probabilities, how="left", left_index=True, right_index=True)
    for l in [m for m in tmp.columns if m.endswith("_y")]:
        l = l.replace("_y", "")
        tmp[l + "_x"] = tmp[l + "_x"].astype("object")
        tmp[l + "_y"] = tmp[l + "_y"].astype("object")
        tmp[l] = tmp[l + "_y"].fillna(tmp[l + "_x"])
        tmp[l] = tmp[l].astype("category")
        tmp.drop([l + "_y", l + "_x"], axis=1, inplace=True)
    adata.obs = tmp.copy()

    adata.obs[labels_key + "_conf_scANVI"] = 0
    adata.obs[labels_key + "_conf_scANVI"] = adata.obs[labels_key + "_conf_scANVI"].astype("float")
    adata.obs = adata.obs.copy()

    for i in adata.obs[labels_key + "_scANVI"].cat.categories:
        adata.obs[i] = adata.obs[i].astype("float")
        adata.obs.loc[adata.obs[labels_key + "_scANVI"] == i, labels_key + "_conf_scANVI"] = adata.obs[adata.obs[labels_key + "_scANVI"] == i][i]
    
    to_pass = [labels_key + "_scANVI", labels_key + "_conf_scANVI"]
    to_pass.extend(adata.obs[labels_key + "_scANVI"].cat.categories)
    probabilities = adata.obs.loc[:, to_pass]
    
    return (label_model, probabilities)

'''
Writes AnnData objects to disk that contain scANVI results and UMAP projections based on the latent representation.

Input arguements:
adata_query: (AnnData) Object with unknown cells

adata_ref: (AnnData) Object with reference labels

split_key: (str) scANVI metadata value to iteratively subset and split on (e.g. subclass_scANVI)

groupby: (str) Label predicted within the split_key (e.g. cluster if split_key is subclass_scANVI)

output_dir: (str) Location to write AnnData object

date: (str) Datestamp on the iterative_scANVI results file in the output_dir

model_args: (dict): Changes to made to scVI_model_args during training (e.g. {"n_top_genes": 5000})

**kwargs: (dict) Passed to several functions, details below:

    n_cores: (int, default 1) Number of CPU cores to use when constructing the nearest neighbor graph
    
    normalize_data: (bool, default False) Whether to log-normalize AnnData.X
    
    calculate_umap: (bool, default True) Whether to project cells into 2D with UMAP


Outputs:
Writes AnnData objects to disk at output_dir/<split_key value>_scANVI.<date>.h5ad

Example Usage:
save_anndata_kwargs = {
    **{"n_cores": 32}
}

save_anndata(
    adata=adata,
    adata_ref=adata_ref,
    split_key="subclass_scANVI",
    groupby="cluster",
    output_dir=os.path.join("output_scANVI"),
    results_file="iterative_scANVI_results.2022-02-14.csv",
    model_args={"n_top_genes": 5000},
    **save_anndata_kwargs
)
'''

def save_anndata(adata_query, adata_ref, split_key, groupby, output_dir, date, diagnostic_plots=None, model_args={}, **kwargs):
        
    default_kwargs = {
        "n_cores": 1,
        "normalize_data": False,
        "calculate_umap": True,
        "cluster_cells": False
    }
    
    kwargs = {**default_kwargs, **kwargs}
    
    n_cores = kwargs["n_cores"]
    normalize_data = kwargs["normalize_data"]
    calculate_umap = kwargs["calculate_umap"]
    cluster_cells = kwargs["cluster_cells"]
    
    results_file = "iterative_scANVI_results." + date + ".csv"
    
    if isinstance(adata_query, anndata.AnnData) is False or isinstance(adata_ref, anndata.AnnData) is False:
        raise TypeError("One or more of the AnnData objects have an incorrect type,")
        
    if adata_query.shape[1] != adata_ref.shape[1]:
        common_labels = np.intersect1d(adata_query.var_names, adata_ref.var_names)
        adata_ref = adata_ref[:, common_labels].copy()
        adata_query = adata_query[:, common_labels].copy()
        print("WARNING: Reference and query AnnData objects have different shapes, using " + str(len(common_labels)) + " common labels. This may have a deterimental effect on model performance.")
    
    if split_key is not None:
        if split_key in adata_ref.obs.columns is False or split_key in adata_query.obs.columns is False:
            raise KeyError("One or more labels_keys do not exist in the reference AnnData object.")
        
    if os.path.exists(os.path.join(output_dir, results_file)) is False:
        raise ValueError("Output directory lacks an iterative scANVI results file.")
    
    adata = adata_query.concatenate(adata_ref, index_unique=None)
    del adata_query
    
    try:
        scANVI_results = pd.read_csv(os.path.join(output_dir, results_file), index_col=0)
        scANVI_results.index = [str(l) for l in scANVI_results.index]

        if scANVI_results.shape[0] != adata.shape[0]:
            common_cells = np.intersect1d(adata.obs_names, scANVI_results.index)
            adata = adata[common_cells].copy()
            print("WARNING: Mismatch between cells in scANVI results and merged AnnData object, using " + str(len(common_cells)) + " common cells. Was this expected?") 
            
        if groupby in adata.obs.columns:
            adata.obs.drop([groupby], axis=1, inplace=True)
        
        scANVI_results = scANVI_results.loc[:, np.setdiff1d(scANVI_results.columns, adata.obs.columns)]
        
        adata.obs = pd.concat([adata.obs, scANVI_results.loc[adata.obs_names, :]], axis=1)
        
    except:
        print("WARNING: Error merging scANVI results, saving AnnData without them.")
        pass
        
    if os.path.exists(os.path.join(output_dir, "objects")) is False:
        os.makedirs(os.path.join(output_dir, "objects"))
        
    for j in adata.obs.columns:
        if adata.obs[j].dtype == bool:
            print("Correcting bool for " + j)
            adata.obs[j] = adata.obs[j].astype("object")
            adata.obs[j] = adata.obs[j].replace({True: "True", False: "False"})

    for i in adata.obs.columns[adata.obs.isna().sum(axis=0) > 0]:
        if any(adata.obs[i].notna()) is False:
            print("Dropping no-value column " + i)
            adata.obs.drop([i], axis=1, inplace=True)
                
        else:            
            replace_with = ""
            
            if isinstance(adata.obs.loc[adata.obs[i].notna(), i][0], np.float64) is True or isinstance(adata.obs.loc[adata.obs[i].notna(), i][0], np.float32) is True:
                replace_with = 0.0

            if isinstance(adata.obs[i].dtype, pd.core.dtypes.dtypes.CategoricalDtype) is True:
                adata.obs[i] = adata.obs[i].astype("object")
            
            print("Replacing NaNs with " + str(replace_with) + " for " + i + " with dtype " + str(type(adata.obs.loc[adata.obs[i].notna(), i][0])))
            adata.obs.loc[adata.obs[i].isna(), i] = replace_with
            
            if isinstance(adata.obs.loc[(adata.obs[i].notna()) & (adata.obs[i] != ""), i][0], bool) is True:
                adata.obs[i] = [str(l) for l in adata.obs[i]]
            
    if split_key != None:
            
        if isinstance(adata.obs[split_key].dtype, pd.core.dtypes.dtypes.CategoricalDtype) is False:
            adata.obs[split_key] = adata.obs[split_key].astype("category")

        if len(adata.obs[split_key].cat.categories) == 1:
            print("WARNING: Chosen split key has only 1 category.")
            
        splits = adata.obs[split_key].cat.categories
        model_split_key = split_key.replace("_scANVI", "")
        
    else:
        splits = ["All"]
        model_split_key = None

    for i in splits:            
        if i == "" or os.path.exists(os.path.join(output_dir, "objects", i.replace("/", " ") + "_scANVI." + date + ".h5ad")) == True:
            continue
            
        elif i != "All":
            cells = adata.obs[split_key] == i
            sub = adata[cells].copy()
        
        else:
            cells = adata.obs_names
            sub = adata
            
        print(i)
        
        with parallel_backend('threading', n_jobs=n_cores):
            if normalize_data is True:
                sc.pp.normalize_total(sub, 1e4)
                sc.pp.log1p(sub)

            if calculate_umap is True:
                model_name, label_model_name = get_model_names(model_split_key, i, groupby, **model_args)
                
                if os.path.exists(os.path.join(output_dir, "scANVI_models", label_model_name)) is False:
                    # To do, implement option to perform 1-off training of an scVI model
                    print("WARNING: Cannot find the scVI model, did you run iterative_scANVI?")
                
                else:
                    markers = pd.read_csv(os.path.join(output_dir, "scANVI_models", label_model_name, "var_names.csv"), header=None)
                    markers = markers[0].to_list()
                    sub_markers_only = sub[:, markers].copy()
                    
                    try:
                        sub_markers_only.obs[groupby] = sub_markers_only.obs[groupby].astype("object")
                        sub_markers_only.obs[split_key] = sub_markers_only.obs[split_key].astype("object")
                        sub_markers_only.obs[split_key.replace("_scANVI", "")] = sub_markers_only.obs[split_key.replace("_scANVI", "")].astype("object")
                        sub_markers_only.obs.loc[(sub_markers_only.obs[split_key.replace("_scANVI", "")] != sub_markers_only.obs[split_key]), groupby] = "Unknown"
                        sub_markers_only.obs.loc[sub_markers_only.obs[groupby] == "", groupby] = "Unknown"
                        sub_markers_only.obs[groupby] = sub_markers_only.obs[groupby].astype("category")
                    except KeyError:
                        pass
                    
                    scvi.model.SCANVI.setup_anndata(
                        sub_markers_only,
                        layer=model_args["layer"],
                        batch_key=model_args["batch_key"],
                        categorical_covariate_keys=model_args["categorical_covariate_keys"],
                        continuous_covariate_keys=model_args["continuous_covariate_keys"],
                        labels_key=groupby
                    )

                    label_model = scvi.model.SCANVI.load(os.path.join(output_dir, "scANVI_models", label_model_name), sub_markers_only)
                    sub.obsm["X_scVI"] = label_model.get_latent_representation()
                    sc.pp.neighbors(sub, use_rep="X_scVI")
                    sc.tl.umap(sub)
                    
        if diagnostic_plots != None:
            plt.rcParams["figure.figsize"] = (10,10)
            
            sc.pp.subsample(sub, fraction=1)
                        
            for k in diagnostic_plots:
                sub.obs[k] = sub.obs[k].replace("", np.nan)
            
            sc.pl.umap(sub,
                       color=diagnostic_plots,
                       size=50,
                       ncols=2,
                       cmap="YlGnBu",
                       na_color="lightgrey",
                       sort_order=False)
            
            sub.obs[groupby] = sub.obs[groupby].replace("Unknown", np.nan)
            
            sc.pl.umap(sub,
                       color=[groupby, groupby + "_scANVI", groupby + "_conf_scANVI"],
                       size=50,
                       ncols=2,
                       na_color="lightgrey",
                       legend_loc="on data",
                       sort_order=False)
        
        if cluster_cells is True:
            
            sc.tl.leiden(sub, resolution=5)
            sub.obs["leiden_original"] = sub.obs["leiden"].copy()
            sub.obs["metacell"] = 0
            
            while sub.obs["metacell"].equals(sub.obs["leiden"]) == False:
                projected_sub = sub[:, 0:10].copy()
                projected_sub.X = sub.obsm["X_scVI"]
                projected_sub.var_names = [str(j) + "_scVI" for j in range(10)]
                projected_sub.var.drop(["gene_ids"], axis=1, inplace=True)
                del projected_sub.layers["UMIs"]

                obs_by_cluster = dict()
                tmp = sub.obs["leiden"].reset_index().drop("index_name", axis=1).reset_index()
                for j in tmp["leiden"].cat.categories:
                    obs_by_cluster[j] = tmp.loc[tmp["leiden"] == j, "index"].to_list()

                thresholds = {
                    'q1_thresh': 0.5,
                    'q2_thresh': None,
                    'cluster_size_thresh': 15,
                    'qdiff_thresh': 0.7,
                    'padj_thresh': 0.05,
                    'lfc_thresh': 1.0,
                    'score_thresh': 200,
                    'low_thresh': 1,
                    'min_genes': 1,
                }
                
                cluster_assignments_after_merging = tc.merge_clusters(
                    adata_norm=sub,
                    adata_reduced=projected_sub,
                    cluster_assignments=obs_by_cluster,
                    cluster_by_obs=sub.obs["leiden"].to_list(),
                    thresholds=thresholds,
                    de_method='ebayes'
                )
 
                if len(cluster_assignments_after_merging) > 1:
                    cluster_assignments_after_merging = cluster_assignments_after_merging[0]

                cluster_by_obs_after_merging = np.zeros(sub.shape[0], dtype=int)
                for cluster, obs in cluster_assignments_after_merging.items():
                    cluster_by_obs_after_merging[obs] = cluster

                sub.obs["metacell"] = sub.obs["leiden"].copy()
                sub.obs["metacell"] = sub.obs["metacell"].astype("category")
                sub.obs["leiden"] = cluster_by_obs_after_merging.astype("str")
                sub.obs["leiden"] = sub.obs["leiden"].astype("category")

            sub.obs["leiden_merged"] = sub.obs["leiden"].copy()
            sub.obs["leiden"] = sub.obs["leiden_original"].copy()
            sub.obs.drop(["metacell"], axis=1, inplace=True)
            
            tmp = pd.concat(
                [
                    sub.obs[["fraction_mito", "leiden_merged"]].groupby("leiden_merged").mean(),
                    sub.obs[["doublet_score", "leiden_merged"]].groupby("leiden_merged").mean()
                ],
                axis=1
            )
                        
            tmp.columns = ["cluster_fraction_mito", "cluster_doublet_score"]
            
            conf_mat= sub.obs.groupby(["leiden_merged", "donor_name"]).size().unstack(fill_value=0)
            conf_mat = conf_mat.div(conf_mat.sum(axis=1), axis=0)
            conf_mat[conf_mat.isna()] = 0
            tmp["cluster_donor_split"] = conf_mat.max(axis=1)
            
            sub.obs = sub.obs.merge(tmp, how="left", left_on="leiden_merged", right_index=True)
            
        sub.write(os.path.join(output_dir, "objects", i.replace("/", " ") + "_scANVI." + date + ".h5ad"))

'''
Gets the scVI and scANVI model name based on args passed.
Called by save_anndata.

Input arguements:
split_key: (str) scANVI metadata value to iteratively subset and split on (e.g. subclass_scANVI)

split_value: (str) Specific value for the split_key (e.g. Astro)

groupby: (str) Label predicted within the split_key (e.g. cluster if split_key is subclass_scANVI)

**kwargs: (dict) Passed to construct model_name and label_model_name

    layer: (None or str, default None) None if unnormalized counts are in AnnData.X, else a str where they are stored in AnnData.layers
    
    categorical_covariate_keys: (list) List of categorical covariates to pass to scVI and scANVI (e.g. ["donor_name"])
    
    continuous_covariate_keys: (list) List of continuous covariates to pass to scVI and scANVI (e.g. ["n_genes"])
    
    use_hvg: (bool, default True) Whether to calculate and include highly variable genes from the reference dataset to scVI and scANVI
    
    use_de: (bool, default True) Whether to calculate and include differentially genes from the reference dataset to scVI and scANVI
    
    n_top_genes: (int or list, default 2000) Number of highly variable genes to pass in each iteration. Not used if use_hvg is False
    
    n_downsample_ref: (int or list, default 1000) Number of cells to downsample reference groups too in each iteration when 
    calculcating differentially expressed genes
    
    n_ref_genes: (int or list, default 500): Number of differentially expressed genes per group to pass in each iteration from the 
    reference dataset to scVI and scANVI
    
    max_epochs_scVI: (int or list, default 200) Number of epochs to train scVI for in each iteration
    
    max_epochs_scANVI: (int or list, default 20) Number of epochs to train scANVI for in each iteration 
    
    scVI_model_args: (None or dict, default {"n_layer": 2}) kwargs passed to scvi.model.SCVI
    
    scANVI_model_args: (None or dict) kwargs passed to scvi.model.SCANVI.from_scvi_model


Outputs:
Returns tupple with scVI and scANVI model names (str)
'''

def get_model_names(split_key, split_value, groupby, **kwargs):
    
    default_kwargs = {
        "layer": "UMIs",
        "batch_key": None,
        "categorical_covariate_keys": None,
        "continuous_covariate_keys": None,
        "use_hvg": True,
        "use_de": True,
        "n_top_genes": 2000,
        "n_downsample_ref": 1000,
        "n_ref_genes": 500,
        "max_epochs_scVI": 200,
        "max_epochs_scANVI": 20,
        "scVI_model_args": {"n_layers": 2},
        "scANVI_model_args": {},
        "user_genes": None
    }
        
    kwargs = {**default_kwargs, **kwargs}
    
    layer = kwargs["layer"]
    batch_key = kwargs["batch_key"]
    categorical_covariate_keys = kwargs["categorical_covariate_keys"]
    continuous_covariate_keys = kwargs["continuous_covariate_keys"]
    use_hvg = kwargs["use_hvg"]
    use_de = kwargs["use_de"]
    n_top_genes = kwargs["n_top_genes"]
    n_downsample_ref = kwargs["n_downsample_ref"]
    n_ref_genes = kwargs["n_ref_genes"]
    max_epochs_scVI = kwargs["max_epochs_scVI"]
    max_epochs_scANVI = kwargs["max_epochs_scANVI"]
    scVI_model_args = kwargs["scVI_model_args"]
    scANVI_model_args = kwargs["scANVI_model_args"]
    user_genes = kwargs["user_genes"]
    
    get_model_genes_kwargs = {
        "groupby": groupby,
        "use_hvg": use_hvg,
        "use_de": use_de,
        "n_top_genes": n_top_genes,
        "n_downsample_ref": n_downsample_ref,
        "n_ref_genes": n_ref_genes,
        "user_genes": user_genes
    }
    run_scVI_kwargs = {
        "layer": layer,
        "max_epochs_scVI": max_epochs_scVI,
        "batch_key": batch_key,
        "categorical_covariate_keys": categorical_covariate_keys,
        "continuous_covariate_keys": continuous_covariate_keys,
        "scVI_model_args": scVI_model_args,
    }
    run_scANVI_kwargs = {
        "layer": layer,
        "max_epochs_scANVI": max_epochs_scANVI,
        "batch_key": batch_key,
        "categorical_covariate_keys": categorical_covariate_keys,
        "continuous_covariate_keys": continuous_covariate_keys,
        "labels_key": groupby,
        "scANVI_model_args": scANVI_model_args,
    }
    
    if split_key == None:
        model_name = hashlib.md5(str(json.dumps({**get_model_genes_kwargs, **run_scVI_kwargs})).replace("/", " ").encode()).hexdigest()
        label_model_name = hashlib.md5(str(json.dumps({**get_model_genes_kwargs, **run_scANVI_kwargs})).replace("/", " ").encode()).hexdigest()
        
    else:
        model_name = hashlib.md5(str(json.dumps({**{split_key: split_value}, **get_model_genes_kwargs, **run_scVI_kwargs})).replace("/", " ").encode()).hexdigest()
        label_model_name = hashlib.md5(str(json.dumps({**{split_key: split_value}, **get_model_genes_kwargs, **run_scANVI_kwargs})).replace("/", " ").encode()).hexdigest()
        
    return (model_name, label_model_name)

def clean_taxonomies(groups, splitby, reference_key, object_dir, model_args={}, **kwargs):
    
    default_kwargs = {
       "n_cores": 1,
        "layer": "UMIs",
        "categorical_covariate_keys": None,
        "continuous_covariate_keys": None,
        "diagnostic_plots": None,
        "use_markers": False,
        "rerun_scVI": True,
        "refine_supertypes": True,
        "save_results": True,
    }
        
    kwargs = {**default_kwargs, **kwargs}
        
    n_cores = kwargs["n_cores"]
    layer = kwargs["layer"]
    categorical_covariate_keys = kwargs["categorical_covariate_keys"]
    continuous_covariate_keys = kwargs["continuous_covariate_keys"]
    diagnostic_plots = kwargs["diagnostic_plots"] 
    use_markers = kwargs["use_markers"]
    rerun_scVI = kwargs["rerun_scVI"]
    refine_supertypes = kwargs["refine_supertypes"]
    save_results = kwargs["save_results"]

    
    if isinstance(groups, dict) is False:
        raise TypeError("Groups must be a dictionary with format 'name': 'type'.")
    
    for i,j in groups.items():
        print(i)
        
        if j["type"] == "glia":
            default_cutoffs = {
                "cluster_donor_split": (0.5, "lt"),
                "nFeature_RNA": (1500, "gt"),
                "fraction_mito": (0.05, "lt"),
                "doublet_score": (0.1, "lt")
            }
        elif j["type"] == "inhibitory":
            default_cutoffs = {
                "cluster_donor_split": (0.5, "lt"),
                "nFeature_RNA": (1500, "gt"),
                "fraction_mito": (0.05, "lt"),
                "doublet_score": (0.4, "lt")
            }
        else:
            default_cutoffs = {
                "cluster_donor_split": (0.5, "lt"),
                "nFeature_RNA": (2000, "gt"),
                "fraction_mito": (0.05, "lt"),
                "doublet_score": (0.4, "lt")
            }
            
        if "cutoffs" in j.keys():
            cutoffs = {**default_cutoffs, **j["cutoffs"]}
        else:
            cutoffs = default_cutoffs
        
        h5ad = glob.glob(os.path.join(object_dir, i.replace("/", " ") + "_*.h5ad"))
        updates = glob.glob(os.path.join(object_dir, i.replace("/", " "), "*-cell-labels-*.csv"))
        
        if len(h5ad) == 0:
            raise ValueError("AnnData object with prefix " + i + " was not found.")
        
        if len(h5ad) > 1 or len(updates) > 1:
            h5ad.sort(reverse=True)
            updates.sort(reverse=True)
            print("WARNING: There are multiple AnnData objects or cell label updates in your search, picking the first one.")
            
        h5ad = h5ad[0]
            
        adata = sc.read_h5ad(h5ad)
        to_keep = [True] * adata.shape[0]

        output = pd.DataFrame(adata.obs[splitby])
        output[splitby] = output[splitby].astype("object")

        if len(updates) != 0:
            updates = updates[0]
            updates = pd.read_csv(updates, index_col=0, skiprows=2)

            flags = updates.columns[updates.columns.str.endswith("_flag")]
            for k in flags:
                to_keep = (to_keep) & (updates[k] == "unassigned")
                output.loc[updates[k] != "unassigned", splitby] = k

            if "cluster_changes" in updates.columns:
                adata.obs[splitby] = adata.obs[splitby].astype("object")
                tmp = updates.loc[updates["cluster_changes"] != "unassigned", "cluster_changes"]
                tmp = tmp.str.replace("L2 3", "L2/3")
                adata.obs.loc[updates["cluster_changes"] != "unassigned", splitby] = tmp
                output.loc[updates["cluster_changes"] != "unassigned", splitby] = tmp
                adata.obs[splitby] = adata.obs[splitby].astype("category")

            if "cluster_reassign" in updates.columns:
                to_keep = (to_keep) & (updates["cluster_reassign"] == "unassigned")
                tmp = updates.loc[updates["cluster_reassign"] != "unassigned", "cluster_reassign"]
                tmp = [k + "_reassign" for k in tmp]
                output.loc[updates["cluster_reassign"] != "unassigned", splitby] = tmp

        for k,l in cutoffs.items():
            if k == "percent_mito":
                k = "fraction_mito"
            if l[1] == "lt":
                to_keep = (to_keep) & (adata.obs[k] < l[0])
                output.loc[~(adata.obs[k] < l[0]), splitby] = k
            if l[1] == "gt":
                to_keep = (to_keep) & (adata.obs[k] > l[0])
                output.loc[~(adata.obs[k] > l[0]), splitby] = k

        adata = adata[to_keep].copy()
        markers = adata.var_names
        if rerun_scVI is True or "X_umap" not in adata.obsm.keys():
            if os.path.exists(os.path.join(object_dir, i.replace("/", " "), "scVI_model")) is False:
                if use_markers is True:
                    get_model_genes_kwargs = {
                        "groupby": splitby,
                        "use_hvg": True,
                        "use_de": True,
                        "n_top_genes": 2000,
                        "n_downsample_ref": 1000,
                        "n_ref_genes": 500,
                        "user_genes": None
                    }
                    markers = get_model_genes(adata[adata.obs[reference_key] == 1], **get_model_genes_kwargs)
            else:
                if use_markers is True:
                    markers = pd.read_csv(os.path.join(object_dir, i.replace("/", " "), "scVI_model", "var_names.csv"), header=None)
                    markers = markers[0].to_list()

            default_model_args = {
                "n_layers": 2,
                "dispersion": "gene-label",
                "n_latent": 10
            }

            sub = adata[:, markers].copy()

            model_args = {**default_model_args, **model_args}

            scvi.model.SCVI.setup_anndata(
                sub,
                layer=layer,
                categorical_covariate_keys=categorical_covariate_keys,
                continuous_covariate_keys=continuous_covariate_keys,
                labels_key=splitby
            )

            if os.path.exists(os.path.join(object_dir, i.replace("/", " "), "scVI_model")) is False:
                model = scvi.model.SCVI(sub, **model_args)
                model.train(max_epochs=200, early_stopping=True)
                model.save(os.path.join(object_dir, i.replace("/", " "), "scVI_model"))

            else:
                model = scvi.model.SCVI.load(os.path.join(object_dir, i.replace("/", " "), "scVI_model"), sub)

            with parallel_backend('threading', n_jobs=n_cores):
                adata.obsm["X_scVI"] = model.get_latent_representation()
                sc.pp.neighbors(adata, use_rep="X_scVI")
                sc.tl.umap(adata)
            
        if diagnostic_plots != None:
            plt.rcParams["figure.figsize"] = (10,10)
            
            sc.pp.subsample(adata, fraction=1)
                        
            for k in diagnostic_plots:
                adata.obs[k] = adata.obs[k].replace("", np.nan)
            
            sc.pl.umap(adata,
                       color=diagnostic_plots,
                       size=50,
                       ncols=2,
                       cmap="YlGnBu",
                       na_color="lightgrey",
                       sort_order=False)
            
            adata.obs[splitby.replace("_scANVI", "")] = adata.obs[splitby.replace("_scANVI", "")].replace("Unknown", np.nan)
            
            sc.pl.umap(adata,
                       color=[splitby, splitby.replace("_scANVI", ""), splitby.replace("_scANVI", "_conf_scANVI")],
                       size=50,
                       ncols=2,
                       na_color="lightgrey",
                       legend_loc="on data",
                       sort_order=False)
            
        if refine_supertypes == True:
            
            sc.tl.leiden(adata, resolution=5, key_added="leiden2")
            adata.obs["leiden_original2"] = adata.obs["leiden2"].copy()
            adata.obs["metacell"] = 0
            
            while adata.obs["metacell"].equals(adata.obs["leiden2"]) == False:
                
                output[splitby + "_leiden"] = output[splitby].copy()
                projected_adata = adata[:, 0:model_args["n_latent"]].copy()
                projected_adata.X = adata.obsm["X_scVI"]
                projected_adata.var_names = [str(k) + "_scVI" for k in range(model_args["n_latent"])]
                projected_adata.var.drop(["gene_ids"], axis=1, inplace=True)
                del projected_adata.layers[layer]

                obs_by_cluster = dict()
                tmp = adata.obs["leiden2"].reset_index().drop("index_name", axis=1).reset_index()
                for k in tmp["leiden2"].cat.categories:
                    obs_by_cluster[k] = tmp.loc[tmp["leiden2"] == k, "index"].to_list()
                    
                if j["type"] != "glia":
                    merge_thresholds = {
                        'q1_thresh': 0.75,
                        'q2_thresh': None,
                        'cluster_size_thresh': 100,
                        'qdiff_thresh': 0.60,
                        'padj_thresh': 1e-5,
                        'lfc_thresh': 1.3,
                        'score_thresh': 200,
                        'low_thresh': 1,
                        'min_genes': 1
                    }
                else:
                    merge_thresholds = {
                        'q1_thresh': 0.50,
                        'q2_thresh': None,
                        'cluster_size_thresh': 100,
                        'qdiff_thresh': 0.60,
                        'padj_thresh': 1e-5,
                        'lfc_thresh': 0.9,
                        'score_thresh': 200,
                        'low_thresh': 1,
                        'min_genes': 1
                    }

                cluster_assignments_after_merging, marker_genes = tc.merge_clusters(
                    adata_norm=adata,
                    adata_reduced=projected_adata,
                    cluster_assignments=obs_by_cluster,
                    cluster_by_obs=adata.obs["leiden2"].to_list(),
                    thresholds=merge_thresholds,
                    de_method='ebayes'
                )

                cluster_by_obs_after_merging = np.zeros(sub.shape[0], dtype=int)
                for cluster, obs in cluster_assignments_after_merging.items():
                    cluster_by_obs_after_merging[obs] = cluster

                adata.obs["metacell"] = adata.obs["leiden2"].copy()
                adata.obs["metacell"] = adata.obs["metacell"].astype("category")
                adata.obs["leiden2"] = cluster_by_obs_after_merging.astype("str")
                adata.obs["leiden2"] = adata.obs["leiden2"].astype("category")
        
            adata.obs["leiden_merged2"] = adata.obs["leiden2"].copy()
            adata.obs["leiden2"] = adata.obs["leiden_original2"].copy()
            adata.obs.drop(["metacell"], axis=1, inplace=True)
        
            sc.pl.umap(adata,
                       color=["leiden_merged2"],
                       ncols=1,
                       legend_loc="on data")
            
            plt.figure(figsize=(10,10))
            conf_mat = adata.obs.groupby(["leiden_merged2", splitby]).size().unstack(fill_value=0)
            conf_mat = conf_mat.div(conf_mat.sum(axis=1), axis=0)
            conf_mat[conf_mat.isna()] = 0
            if 1 not in conf_mat.shape:
                ax = sns.clustermap(conf_mat, method="ward", yticklabels=1, xticklabels=1);
                plt.setp(ax.ax_heatmap.yaxis.get_majorticklabels(), rotation=0);

            plt.figure(figsize=(10,10))
            conf_mat = adata.obs.groupby(["leiden_merged2", "donor_name"]).size().unstack(fill_value=0)
            conf_mat = conf_mat.div(conf_mat.sum(axis=1), axis=0)
            conf_mat[conf_mat.isna()] = 0
            if 1 not in conf_mat.shape:
                ax = sns.clustermap(conf_mat, method="ward", yticklabels=1, xticklabels=1);
                plt.setp(ax.ax_heatmap.yaxis.get_majorticklabels(), rotation=0);


            leiden_counts = adata.obs.groupby([splitby, "leiden_merged2"]).size()
            for k in adata.obs[splitby].cat.categories:
                l = sp_stats.entropy(leiden_counts[k] / leiden_counts[k].sum())
                print(str(k) + " has entropy " + str(l))  

            adata.obs[splitby + "_leiden"] = adata.obs[splitby].copy()
            supertype_counts = adata.obs.groupby(["leiden_merged2", splitby.replace("_scANVI", "")]).size()
            reference_counts = adata.obs.loc[:, ["leiden_merged2", reference_key]].groupby(["leiden_merged2"]).sum() / adata.obs[reference_key].value_counts().loc[1]
            for k in adata.obs["leiden_merged2"].cat.categories:
                l = np.nanmax(supertype_counts[k] / adata.obs.groupby([splitby.replace("_scANVI", "")]).size())
                m = reference_counts.loc[k, :].max()
                print(str(k) + " has supertype max " + str(l) + " and reference percent " + str(m))

                if l < 0.1 and m < 0.1:
                    print("--Keeping for analysis")
                    adata.obs[splitby + "_leiden"] = adata.obs[splitby + "_leiden"].astype("object")
                    adata.obs.loc[adata.obs["leiden_merged2"] == k, splitby + "_leiden"] = i + "_Unknown_" + str(k)
                    output.loc[adata.obs_names[adata.obs["leiden_merged2"] == k], splitby + "_leiden"] = i + "_Unknown_" + str(k)
                    adata.obs[splitby + "_leiden"] = adata.obs[splitby + "_leiden"].astype("category")
                    
            sc.pl.umap(adata, color=splitby + "_leiden", legend_loc="on data", size=10)
        if save_results is True:
            output.to_csv(os.path.join(object_dir, i.replace("/", " "), "corrections.csv"))