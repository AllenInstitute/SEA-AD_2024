import os
import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sc
import copy
import pickle
from matplotlib import pyplot as plt
import warnings
from datetime import datetime
from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from scipy import sparse as sp_sparse

warnings.filterwarnings("ignore")
pwd = os.getcwd()

def run_scCODA(cell_count, random_effect, split_key, split_value, labels_keys, tests, region, covariates, formula):
    for j in labels_keys:
        print(str(datetime.now()) + " -- Starting labels_keys=" + j)
        cell_count = cell_count.loc[[i in split_value for i in cell_count[split_key]], :]
        counts = sp_sparse.csr_matrix(np.zeros((cell_count.shape[0], 0)), dtype=np.float32)
        adata = ad.AnnData(counts)
        adata.obs_names = cell_count.index
        adata.obs = cell_count

        abundances = dat.from_scanpy(
            adata,
            cell_type_identifier=j,
            sample_identifier=random_effect,
        )
        tmp = adata.obs.loc[:, covariates]
        tmp.drop_duplicates(inplace=True)
        tmp.index = tmp[random_effect].copy()
        abundances.obs = tmp.loc[abundances.obs.index, :]
        abundances.obs.index.name = "random_effect"
            
        outfile = " ".join(split_value) + "_" + j + "_results.csv"
        for i in abundances.obs.columns:
            if abundances.obs[i].dtype == "category":
                abundances.obs[i] = abundances.obs[i].astype("object")
                abundances.obs[i] = [str(n) for n in abundances.obs[i]]
                
        if os.path.exists(os.path.join(pwd, "output", region, "objects")) == False:
            os.makedirs(os.path.join(pwd, "output", region, "objects"))
            
        abundances.write(os.path.join(pwd, "output", region, "objects", " ".join(split_value) + "_" + j + "_abundances.h5ad"))
              
        for k in tests:
            print(str(datetime.now()) + " -- testing across " + k)
            i = 0
            cell_types = abundances.var.index
            try:
                del results_table
            except:
                pass
            
            if os.path.exists(os.path.join(pwd, "output", region, k, outfile)) == True:
                continue
            
            while i < cell_types.shape[0]:    

                ct = cell_types[i]
                print("Reference: " + ct)

                pickle_file = ct.replace("/", " ").replace(":", "") + ".pkl"
                if os.path.exists(os.path.join(pwd, "output", region, k, pickle_file)):
                    with open(os.path.join(pwd, "output", region, k, pickle_file), 'rb') as pickle_load:
                        results = pickle.load(pickle_load)

                else:
                    models = mod.CompositionalAnalysis(abundances, formula=formula + k, reference_cell_type=ct)
                    results = models.sample_hmc()
                    if os.path.exists(os.path.join(pwd, "output", region, k)) == False:
                        os.makedirs(os.path.join(pwd, "output", region, k))
                        
                    results.save(pickle_file)

                accepted = results.sample_stats["is_accepted"].to_numpy()

                if accepted.sum() / accepted.shape[1] < 0.6:
                    del results
                    os.remove(pickle_file)
                    print("Did not achieve acceptance threshold, trying again")
                    continue

                else:
                    print("Converged!")
                    i = i + 1

                current_results = results.summary_prepare()[1]
                current_results = current_results.reset_index()
                current_results["Reference Cell Type"] = ct

                try:
                    results_table = pd.concat([results_table, current_results], axis=0)
                except:
                    results_table = current_results.copy()

        if os.path.exists(os.path.join(pwd, "output", region, k)) == False:
            os.makedirs(os.path.join(pwd, "output", region, k))
                        
        results_table.to_csv(os.path.join(pwd, "output", region, k, outfile))
        
