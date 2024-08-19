import numpy as np
import scanpy as sc

def compute_cell_quality(adata_mvi, cell_idx):
    idx = np.where(adata_mvi.uns['neighbors']['connectivities'][cell_idx].todense()>0)[1]
    df = adata_mvi.obs[["Doublet_or_LowQuality", "modality"]].iloc[idx]
    df = df.loc[df["modality"].isin(["paired", "expression"])]
    ratio = np.sum(df["Doublet_or_LowQuality"] == "RNA doublet or LQ cells") / df["Doublet_or_LowQuality"].shape[0]
    return ratio

def compute_cell_mixing(adata_mvi, cell_idx):
    idx = np.where(adata_mvi.uns['neighbors']['connectivities'][cell_idx].todense()>0)[1]
    df = adata_mvi.obs[["modality"]].iloc[idx]
    ratio = np.sum(df["modality"] == adata_mvi.obs["modality"][cell_idx]) / df["modality"].shape[0]
    return ratio / theoretic_score

def compute_label_purity(adata_mvi, cell_idx):
    idx = np.where(adata_mvi.uns['neighbors']['connectivities'][cell_idx].todense()>0)[1]
    df = adata_mvi.obs[["subclass_scANVI"]].iloc[idx]
    df = df.loc[~df["subclass_scANVI"].isnull()]
    u, c = np.unique(df, return_counts=True)
    if np.size(c) == 0:
        ratio = 0
        label = 'NA'
    else:
        ratio = c[np.argmax(c)] / c.sum()
        label = u[np.argmax(c)]
    
    return ratio, label