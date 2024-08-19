import os
import re
import glob
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import pyreadstat
from sklearn.neighbors import KNeighborsTransformer
from datetime import datetime
from igraph import *
import warnings

sc.settings.n_jobs = 32
warnings.filterwarnings("ignore")

pwd = os.getcwd()

def build_anndata(datadir, outfile, to_correct=[], to_str=[], to_drop=[], metadata_file=None, metadata_name=None, adata_key=None, external_key=None, cell_ids=None, remove_arid=False):
    adata = []

    for n, i in enumerate(os.listdir(datadir)):
        if os.path.isdir(os.path.join(datadir, i)) is False or i.startswith("."):
            continue
            
        print(str(datetime.now()) + " -- " + i)

        tmp = sc.read_10x_mtx(os.path.join(datadir, i))
        tmp.obs = pd.read_csv(os.path.join(datadir, i, "samp.dat.csv"), index_col=0)

        for j in to_correct:
            if tmp.obs[j].dtype == 'int64':
                continue

            tmp.obs[j] = tmp.obs[j].str.replace(",|%", "").astype('float64')
            
        for j in to_str:
            tmp.obs[j] = [str(l) for l in tmp.obs[j].fillna("")]

        adata.append(tmp)

    first = adata.pop()
    adata = first.concatenate(adata, index_unique=None)
    if metadata_file != None:
        if metadata_file.endswith(".txt") or metadata_file.endswith(".tsv"):
            metadata = pd.read_csv(os.path.join(pwd, metadata_file), index_col=0, delimiter="\t")
        elif metadata_file.endswith(".csv"):
            metadata = pd.read_csv(os.path.join(pwd, metadata_file), index_col=0, delimiter=",")
        elif metadata_file.endswith(".sav"):
            metadata = pd.read_spss(os.path.join(pwd, metadata_file))
        else:
            raise ValueError("The metadata file has an unrecognized extension, possible options are tsv, txt, csv, and sav")
            
        if adata_key in adata.obs.columns and external_key in metadata.columns:
            if remove_arid == True:
                adata.obs[cell_ids] = [re.sub("([^_]+)_([^-]+)-[0-9]+$", "\\1_\\2", i) for i in adata.obs[cell_ids]]
                adata.obs.index = adata.obs[cell_ids].copy()
                adata.obs.index.name = "index_name"
            tmp = adata.obs.merge(metadata, left_on=adata_key, right_on=external_key)
            tmp.index = tmp[cell_ids].copy()
            tmp.index.name = "index_name"
            adata = adata[tmp.index].copy()
            adata.obs = tmp.copy()
            adata.uns[metadata_name] = metadata_file
        else:
            raise KeyError("Either the AnnData object or external metadata lacks the supplied merge_on column.") 
        
    if to_drop != []:
        adata.obs.drop(to_drop, axis=1, inplace=True)

    adata.layers["UMIs"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
        
    adata.write(filename=outfile)
    
def build_external_anndata(doi, external_datasets):
    print(" ================================= " + doi + " ================================= ")
    i = external_datasets[doi]
    
    if os.path.exists(os.path.join(pwd, "output", "PFC_" + i["name"] + "_external_singleomeCR6." + str(datetime.date(datetime.now())) + ".h5ad")) == True:
        print("--Skipping")
        return

    results = glob.glob(os.path.join(pwd, "input", "External_AD_singleomeCR6", doi, i["dir"], i["glob_pattern"]))
    adata = []
    
    libraries = pd.read_csv(os.path.join(pwd, "input", "External_AD_singleomeCR6", "metadata", i["name"], "libraries.csv"), index_col=0)
    donors = pd.read_csv(os.path.join(pwd, "input", "External_AD_singleomeCR6", "metadata", i["name"], "donors.csv"), index_col=0)
    
    for j in results:
        if os.path.isdir(j) != True:
            continue
        print(os.path.basename(j))
        
        try:
            tmp = sc.read_10x_h5(os.path.join(j, "outs", "filtered_feature_bc_matrix.h5"))
        except:
            print("WARNING: CANNOT READ LIBRARY")
            continue
        
        if os.path.basename(j) not in libraries["library_prep"].to_list():
            print("WARNING: NO LIBRARY METADATA")
        
        tmp.var = tmp.var.reset_index()
        tmp.var.loc[tmp.var["index"].duplicated(), "index"] = (tmp.var.loc[tmp.var["index"].duplicated(), "index"] + " " + tmp.var.loc[tmp.var["index"].duplicated(), "gene_ids"]).to_list()
        tmp.var.index = tmp.var["index"].copy()
        tmp.var = tmp.var.drop(["index"], axis=1)

        tmp.obs["Number of UMIs"] = tmp.X.sum(axis=1)
        tmp.obs["Genes detected"] = (tmp.X > 0).sum(axis=1)
        tmp.obs["Fraction mitochondrial UMIs"] = tmp[:, tmp.var_names.str.startswith('MT-')].X.sum(axis=1) / tmp.X.sum(axis=1)
        
        obs_names = tmp.obs_names.copy()
        obs_names = [k.replace("-1", "") + "-" + os.path.basename(j) for k in obs_names]
        
        # Merge in library metadata, alignment statitstics, and donor metadata
        if "barcode" not in libraries.columns.to_list():
            tmp.obs["library_prep"] = os.path.basename(j)
            tmp_obs = tmp.obs.merge(libraries, left_on="library_prep", right_on="library_prep", how="left")
        else:
            tmp.obs_names = obs_names.copy()
            libraries["library_prep_barcode"] = libraries["barcode"] + "-" + libraries["library_prep"]
            tmp_obs = tmp.obs.merge(libraries, left_index=True, right_on="library_prep_barcode", how="left")
            tmp_obs.index = tmp_obs["library_prep_barcode"].copy()
            tmp_obs = tmp_obs.drop(["library_prep_barcode"], axis=1)
        
        metrics_summary = pd.read_csv(os.path.join(j, "outs", "metrics_summary.csv"))
        metrics_summary = fix_metrics_summary(metrics_summary)
        metrics_summary["library_prep"] = os.path.basename(j)
        tmp_obs = tmp_obs.merge(metrics_summary, left_on="library_prep", right_on="library_prep", how="left")
        tmp_obs.index = obs_names.copy()
        tmp.obs = tmp_obs.copy()
        tmp = tmp[~tmp.obs["Donor ID"].isna(), :].copy()
        
        if tmp.shape[0] == 0:
            print("WARNING: NO CELLS WERE PRESENT FOR THIS LIBRARY")
            continue

        tmp.obs = tmp.obs.merge(donors, left_on="Donor ID", right_index=True, how="left")

        # Compute doublet scores
        tmp.obs["Doublet score"] = compute_doublet_scores(tmp)
        
        with pd.option_context("display.max_columns", None):
            display(tmp.obs.head())
        
        # Compute log normalized UP10K
        tmp.layers["UMIs"] = tmp.X.copy()
        sc.pp.normalize_total(tmp, target_sum=1e4)
        sc.pp.log1p(tmp)
        
        adata.append(tmp)
    
    # Merge AnnData objects together
    first = adata.pop()
    adata = first.concatenate(adata, index_unique=None)
    
    # Write out AnnData objects for each dataset
    adata.write(os.path.join(pwd, "output", "PFC_" + i["name"] + "_external_singleomeCR6." + str(datetime.date(datetime.now())) + ".h5ad"))
    
    
def compute_doublet_scores(adata, proportion_artificial=0.2):
    k = np.int64(np.round(np.min([100, adata.shape[0] * 0.01])))
    
    n_doublets = np.int64(np.round(adata.shape[0] / (1 - proportion_artificial) - adata.shape[0]))
    real_cells_1 = np.random.choice(adata.obs_names, size=n_doublets, replace=True)
    real_cells_2 = np.random.choice(adata.obs_names, size=n_doublets, replace=True)
    doublet_X = adata[real_cells_1, :].X + adata[real_cells_2, :].X
    doublet_obs_names = ["X" + str(i) for i in range(n_doublets)]
    doublet_adata = ad.AnnData(X=doublet_X, obs=pd.DataFrame(index=doublet_obs_names), var=pd.DataFrame(index=adata.var_names))
    
    adata = adata.concatenate(doublet_adata, index_unique=None)
    adata.obs["doublet_cell"] = adata.obs_names.isin(doublet_obs_names)
    adata.obs["doublet_cell"] = adata.obs["doublet_cell"].astype("category")
    
    adata.layers["UMIs"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=5000, flavor="seurat_v3", layer="UMIs")
        adata.layers["UMIs"]
    except:
        sc.pp.highly_variable_genes(adata, min_mean=1, min_disp=0.5)
        del adata.layers["UMIs"]

    adata_sub = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_sub)
    sc.pp.pca(adata_sub)
    
    v = adata_sub.uns['pca']['variance']
    n_pcs = np.max(np.where(((v - np.mean(v)) / np.std(v)) > 0)[0])
    
    knn = KNeighborsTransformer(
        n_neighbors=k,
        algorithm="kd_tree",
        n_jobs=-1,
    )
    knn = knn.fit(adata_sub.obsm["X_pca"][:, :n_pcs])
    dist, idx = knn.kneighbors()
    
    knn_mapper = KNeighborsTransformer(
        n_neighbors=10,
        algorithm="kd_tree",
        n_jobs=-1,
    )
    knn_mapper = knn_mapper.fit(adata_sub[adata_sub.obs["doublet_cell"] == False, :].obsm["X_pca"][:, :n_pcs])
    dist1, _ = knn_mapper.kneighbors(adata_sub[adata_sub.obs["doublet_cell"] == True, :].obsm["X_pca"][:, :n_pcs])
    
    dist_th = np.mean(dist1) + (1.64 * np.std(dist1))
    freq = (dist < dist_th) & (idx > adata[adata.obs["doublet_cell"] == False, :].shape[0])
    score1 = freq.mean(axis=1)
    score2 = freq[:, :np.int64(np.ceil(k/2))].mean(axis=1)
    adata.obs["Doublet score"] = [np.max([score1[i], score2[i]]) for i in range(adata.shape[0])]
    
    return adata.obs.loc[~adata.obs_names.isin(doublet_obs_names), "Doublet score"]

def fix_metrics_summary(metrics_summary):
    to_remove_commas = [
        "Estimated Number of Cells",
        "Mean Reads per Cell",
        "Median Genes per Cell",
        "Number of Reads",
        "Total Genes Detected",
        "Median UMI Counts per Cell"
    ]
    to_remove_percent_and_divide = [
        "Valid Barcodes",
        "Sequencing Saturation",
        "Q30 Bases in Barcode",
        "Q30 Bases in RNA Read",
        "Q30 Bases in UMI",
        "Reads Mapped to Genome",
        "Reads Mapped Confidently to Genome",
        "Reads Mapped Confidently to Intergenic Regions",
        "Reads Mapped Confidently to Intronic Regions",
        "Reads Mapped Confidently to Exonic Regions",
        "Reads Mapped Confidently to Transcriptome",
        "Reads Mapped Antisense to Gene",
        "Fraction Reads in Cells"
    ]
    to_rename = {
        "Estimated Number of Cells": "Estimated_number_of_cells",
        "Fraction Reads in Cells": "GEX_Fraction_of_transcriptomic_reads_in_cells",
        "Mean Reads per Cell": "GEX_Mean_raw_reads_per_cell",
        "Median Genes per Cell": "GEX_Median_genes_per_cell",
        "Number of Reads": "GEX_number_of_reads",
        "Valid Barcodes": "GEX_Valid Barcodes",
        "Sequencing Saturation": "GEX_sequencing_saturation",
        "median_umi_counts_per_cell": "GEX_Median_UMI_counts_per_cell",
        "Q30 Bases in RNA Read": "GEX_Q30_bases_in_read_2",
        "Q30 Bases in Barcode": "GEX_Q30_bases_in_barcode",
        "Q30 Bases in UMI": "GEX_Q30_bases_in_UMI",
        "Reads Mapped Antisense to Gene": "GEX_Reads_mapped_antisense_to_gene",
        "Reads Mapped Confidently to Exonic Regions": "GEX_Reads_mapped_confidently_to_exonic_regions",
        "Reads Mapped Confidently to Genome": "GEX_Reads_mapped_confidently_to_genome",
        "Reads Mapped Confidently to Intergenic Regions": "GEX_Reads_mapped_confidently_to_intergenic_regions",
        "Reads Mapped Confidently to Intronic Regions": "GEX_Reads_mapped_confidently_to_intronic_regions",
        "Reads Mapped Confidently to Transcriptome": "GEX_Reads_mapped_confidently_to_transcriptome",
        "Reads Mapped to Genome": "GEX_Reads_mapped_to_genome",
        "Total Genes Detected": "GEX_Total_genes_detected",
        "Median UMI Counts per Cell": "GEX_Median_UMI_counts_per_cell"
    }
    for k in to_remove_commas:
        metrics_summary[k] = metrics_summary[k].astype("str")
        metrics_summary[k] = metrics_summary[k].str.replace(",", "")
        metrics_summary[k] = metrics_summary[k].astype("int")
    for k in to_remove_percent_and_divide:
        metrics_summary[k] = metrics_summary[k].str.replace("%", "")
        metrics_summary[k] = metrics_summary[k].astype("float")
        metrics_summary[k] = metrics_summary[k] / 100
    metrics_summary = metrics_summary.rename(
        to_rename,
        axis=1
    )
    return metrics_summary
    
