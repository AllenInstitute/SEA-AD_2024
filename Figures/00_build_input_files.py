import os
import scanpy as sc
import numpy as np
import pandas as pd
from scipy import sparse as sp_sparse
from datetime import datetime

pwd = os.getcwd()

# For Figure 1 and Extended Data Figure 1 and 2
print(str(datetime.now()) + " -- Starting Figure 1 and Extended Data Figure 1 and 2")

# From https://sea-ad-single-cell-profiling.s3.amazonaws.com/index.html#MTG/RNAseq/
adata = sc.read_h5ad(os.path.join(pwd, "input", "SEAAD_MTG_RNAseq_all-nuclei.2024-02-13.h5ad"))
nuclear_cytosolic = sc.get.obs_df(adata, ["MT-CO1", "MT-ND3", "MALAT1", "MEG3", "Severely Affected Donor"], layer="UMIs")
nuclear_cytosolic.to_csv(os.path.join(pwd, "input", "Figure 1 and Extended Data Figure 1 and 2", "nuclear_cytosolic.csv"))
del adata
del nuclear_cytosolic


# For Figure 2 and Extended Data Figure 3 and 4
print(str(datetime.now()) + " -- Starting Figure 2 and Extended Data Figure 3 and 4")


# For Figure 3 and Extended Data 8
print(str(datetime.now()) + " -- Starting Figure 3 and Extended Data Figure 8")

# From https://sea-ad-single-cell-profiling.s3.amazonaws.com/index.html#MTG/RNAseq/
adata = sc.read_h5ad(os.path.join(pwd, "input", "SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"))
adata.X = sp_sparse.csr_matrix(adata.X.shape, dtype=np.int32)
del adata.layers["UMIs"]
adata.write(os.path.join(pwd, "input", "Figure 3 and Extended Data 8", "SEAAD_MTG_RNAseq_final-nuclei_no_data.2024-02-13.h5ad"))
del adata


# For Figure 5 and Extended Data 10
print(str(datetime.now()) + " -- Starting Figure 5 and Extended Data Figure 10")


# For Figure 6 and Extended Data 11
print(str(datetime.now()) + " -- Starting Figure 6 and Extended Data Figure 11")

# From https://sea-ad-single-cell-profiling.s3.amazonaws.com/index.html#MTG/RNAseq/
adata = sc.read_h5ad(os.path.join(pwd, "input", "SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"))
adata = adata[:, ["MAPK8", "NGF", "MME", "LNX2", "ATP5MPL", "RPL4"]].copy()
del adata.layers["UMIs"]
adata.write(os.path.join(pwd, "input", "Figure 6 and Extended Data Figure 11", "SEAAD_MTG_RNAseq_final-nuclei_limited.2024-02-13.h5ad"))
del adata

# Figure 7
print(str(datetime.now()) + " -- Starting Figure 7")

# From https://sea-ad-single-cell-profiling.s3.amazonaws.com/index.html#MTG/RNAseq/
adata = sc.read_h5ad(os.path.join(pwd, "input", "SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"))
adata = adata[:, ["RUNX1", "CTSC", "JAK3", "APOE"]].copy()
del adata.layers["UMIs"]
adata.write(os.path.join(pwd, "input", "Figure 7", "SEAAD_MTG_RNAseq_final-nuclei_limited.2024-02-13.h5ad"))
del adata


# Figure 8
print(str(datetime.now()) + " -- Starting Figure 8")

# From https://sea-ad-single-cell-profiling.s3.amazonaws.com/index.html#MTG/RNAseq/
adata = sc.read_h5ad(os.path.join(pwd, "input", "SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"))
adata = adata[:, ["PSEN1", "PSEN2", "BACE1", "BACE2", "PSENEN", "NCSTN", "APH1A", "APP", "OMG", "DHCR24", "OLIG2", "IGF1", "IGF2", "PDGFA", "PDGFB", "PDGFC"]].copy()
del adata.layers["UMIs"]
adata.write(os.path.join(pwd, "input", "Figure 8", "SEAAD_MTG_RNAseq_final-nuclei_limited.2024-02-13.h5ad"))
del adata


# Extended Data Figure 5
print(str(datetime.now()) + " -- Starting Extended Data Figure 5")

adata = sc.read_h5ad(os.path.join(pwd, "input", "SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"))
adata_ref = adata[adata.obs["Neurotypical reference"] == "True", :].copy()
marker_genes = {}
cells = []
# Downsample reference dataset for differential expression (DE)
for i in adata_ref.obs["Supertype (non-expanded)"].cat.categories:
    
    if i == "Lamp5_Lhx6_1":
        print("Skipping " + i)
        continue
    
    tmp = adata_ref[adata_ref.obs["Supertype (non-expanded)"] == i].obs_names.to_list()
    
    if len(tmp) > 1000:
        cells = cells + random.sample(tmp, k=100)
    
    else:
        cells = cells + tmp
adata_ds = adata_ref[cells].copy()
# Compute DE
for i in adata_ds.obs["Subclass"].cat.categories:
    
    adata_tmp = adata_ds[adata_ds.obs["Subclass"] == i]
    adata_tmp.obs["Supertype (non-expanded)"] = adata_tmp.obs["Supertype (non-expanded)"].cat.remove_unused_categories()
    
    sc.tl.rank_genes_groups(adata_tmp,
                            groupby="Supertype (non-expanded)",
                            method="wilcoxon",
                            tie_correct=True)
    
    result = adata_tmp.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    tmp = {group: pd.DataFrame({key: result[key][group] for key in ['names', 'pvals_adj', 'logfoldchanges']}) for group in groups}
    
    marker_genes = {**marker_genes, **tmp}

# Compute supertype signature score
splitby = "Supertype (non-expanded)"
adata.obs["Supertype Signature Score"] = 0
for i in adata.obs[splitby].cat.categories:
    if i == "Lamp5_Lhx6_1":
        continue
    print(str(datetime.now()) + " -- " + i)
    for j in adata.obs["Continuous Pseudo-progression Score"].dropna().unique():
        if adata[(adata.obs["Used in analysis"] == True) & (adata.obs["Continuous Pseudo-progression Score"] == j) & (adata.obs[splitby] == i), :].shape[0] < 5 or adata[(adata.obs["Used in analysis"] == True) & (adata.obs["Neurotypical reference"] == "True") & (adata.obs[splitby] == i), :].shape[0] < 5:
            continue
        loss = 0
        for k in marker_genes[i].iloc[:50].loc[:, "names"]:
            reference = adata[(adata.obs["Used in analysis"] == True) & (adata.obs["Neurotypical reference"] == "True") & (adata.obs[splitby] == i), k].X.toarray()
            seaad = adata[(adata.obs["Used in analysis"] == True) & (adata.obs["Continuous Pseudo-progression Score"] == j) & (adata.obs[splitby] == i), k].X.toarray()
            loss += np.abs(sp_stats.ttest_ind(reference, seaad)[0].max())

        adata.obs.loc[(adata.obs["Continuous Pseudo-progression Score"] == j) & (adata.obs[splitby] == i), "Supertype Signature Score"] = -loss

adata.X = sp_sparse.csr_matrix(adata.X.shape, dtype=np.int32)
del adata.layers["UMIs"]
adata.write(os.path.join(pwd, "input", "Extended Data Figure 5", "SEAAD_MTG_RNAseq_final-nuclei_no_data_signature_score.2024-02-13.h5ad"))
del adata


# Extended Data Figure 6
print(str(datetime.now()) + " -- Starting Extended Data Figure 6")

adata = sc.read_h5ad(os.path.join(pwd, "input", "SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"))
del adata.layers["UMIs"]
adata = adata[:, []].copy()
adata.write(os.path.join(pwd, "input", "Extended Data Figure 6", "SEAAD_MTG_RNAseq_final-nuclei_no_data.2024-02-13.h5ad"))
del adata

adata = sc.read_h5ad(os.path.join(pwd, "input", "SEAAD_MTG_ATACseq_all-nuclei.2024-02-13.h5ad"))
adata = adata[:, []].copy()
adata.write(os.path.join(pwd, "input", "Extended Data Figure 6", "SEAAD_MTG_ATACseq_all-nuclei_no_data.2024-02-13.h5ad"))
del adata
