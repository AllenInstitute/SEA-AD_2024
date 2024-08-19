import os
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import copy
import re
from joblib import Parallel, delayed
import warnings
from datetime import datetime
from igraph import *
import scvi
import rpy2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import argparse

sc.settings.n_jobs = 32

parser = argparse.ArgumentParser()
parser.add_argument("--threads", help="threads")
parser.add_argument("--target", help="target")
parser.add_argument("--region", help="region")
parser.add_argument("--split_key", help="split_key")
parser.add_argument("--covariates", help="covariates")
parser.add_argument("--random_effect", help="random_effect")
parser.add_argument("--covariate_formula", help="covariate_formula")
parser.add_argument("--tests", help="tests")
parser.add_argument("--layer", help="layer")
parser.add_argument("--offset_variable", help="offset_variable")


args = parser.parse_args()
pwd = os.getcwd()
warnings.filterwarnings("ignore")

def run_de_tests(sub, region, target, split_key, outgroup, covariates, random_effect, covariate_formula, tests, layer, offset_variable):
    r_Matrix = importr("Matrix")
    ro.r("library(stats)")
    ro.r("library(nebula)")
    ro.r("library(here)")
    
    if os.path.exists(os.path.join(pwd, "output", region, "versus_all", outgroup.replace("/", " ") + "_vs_all_DE.csv")) is False:
    
        print(str(datetime.now()) + " -- Starting " + outgroup.replace("/", " ") + " versus all", flush=True)
        sub.obs["comparison"] = "0"
        sub.obs.loc[sub.obs[split_key] == outgroup, "comparison"] = "1"

        if layer != "X":
            tmp = sub.layers[layer].T.tocoo()
        else:
            tmp = sub.X.T.tocoo()

        
        counts = r_Matrix.sparseMatrix(
            i=ro.IntVector(tmp.row + 1),
            j=ro.IntVector(tmp.col + 1),
            x=ro.FloatVector(tmp.data),
            dims=ro.IntVector(tmp.shape)
        )
        
        with localconverter(ro.default_converter + pandas2ri.converter):
            obs = ro.conversion.py2rpy(sub.obs)
        
        ro.r.assign("obs", obs)
        ro.r.assign("counts", counts)
        ro.r("covariates <- c('comparison', " + str(covariates).replace("[", "").replace("]", "") + ")")
        ro.r("df <- model.matrix(~" + covariate_formula + "comparison, data=obs[,covariates])")

        with localconverter(ro.default_converter + pandas2ri.converter):
            df = ro.r("as.data.frame(df)")
                
        ro.r("data_g <- group_cell(count=counts, id=obs$" + random_effect + ", offset=obs$" + offset_variable + ", pred=df)")
        ordered = ro.r("data_g")
        if isinstance(ordered, rpy2.rinterface_lib.sexp.NULLType) is True:
            ro.r("re <- nebula(counts, obs$" + random_effect + ", offset=obs$" + offset_variable + ", pred=df, covariance=TRUE)")
            ro.r("saveRDS(re, '" + os.path.join(pwd, "output", region, "versus_all", outgroup.replace("/", " ") + "_vs_all_DE.rds") + "')")
        else:
            ro.r("re <- nebula(data_g$count, data_g$id, offset=data_g$offset, pred=data_g$pred, covariance=TRUE)")
            ro.r("saveRDS(re, '" + os.path.join(pwd, "output", region, "versus_all", outgroup.replace("/", " ") + "_vs_all_DE.rds") + "')")

        with localconverter(ro.default_converter + pandas2ri.converter):
            results = ro.r("re$summary")
            covariance = ro.r("re$covariance")
            overdispersion = ro.r("re$overdispersion")
            convergence = ro.r("re$convergence")
            random_effect = ro.r("re$random_effect")

        results = results.loc[convergence == 1, :]
        results = results.loc[overdispersion["Subject"] < 1, :]
        gene_index = [j - 1 for j in results["gene_id"].to_list()]
        results.index = adata.var_names[gene_index]

        results.to_csv(os.path.join(pwd, "output", region, "versus_all", outgroup.replace("/", " ") + "_versus_all_DE.csv"))
        print(str(datetime.now()) + " -- " + outgroup.replace("/", " ") + " versus all was written to disk", flush=True)

    else:
        print(str(datetime.now()) + " -- Skipping " + outgroup.replace("/", " ") + " versus all (already exists)", flush=True)
    
    sub = sub[sub.obs[split_key] == outgroup]
    starting_features = sub.shape[1]
    
    if layer != "X":
        counts_per_cell = sub.layers[layer].sum(axis=0) / sub.shape[0]
    else:
        counts_per_cell = sub.X.sum(axis=0) / sub.shape[0]
    
    sub = sub[:, counts_per_cell > 0.005].copy()
    ending_features = sub.shape[1]
    if starting_features != ending_features:
        print(str(datetime.now()) + " -- Removing " + str(starting_features - ending_features) + " features from " + outgroup.replace("/", " ") + " for low numbers of counts per cell", flush=True)
        
    if layer != "X":
        global_nonzero_values = (sub.layers[layer] > 0).sum(axis=0)
    else:
        global_nonzero_values = (sub.X > 0).sum(axis=0)
    to_remove = []
    for z in covariates:
        try:
            sub.obs[z] = (sub.obs[z] - sub.obs[z].min()) / sub.obs[z].max()
            print(str(datetime.now()) + " -- Detected " + z + " as an integer or float, applying a min-max normalization", flush=True)
        except:
            pass
            
        if sub.obs[z].unique().shape[0] < 2:
            to_remove.append(z)
            covariate_formula = covariate_formula.replace(z + " + ", "")
            print(str(datetime.now()) + " -- Removing " + z + " from the covariate formula", flush=True)
            continue
        
        if sub.obs[z].unique().shape[0] > 10:
            continue
        
        for y in sub.obs[z].unique():
            if layer != "X":
                spec_nonzero_values = (sub[sub.obs[z] == y, :].layers[layer] > 0).sum(axis=0)
            else:
                spec_nonzero_values = (sub[sub.obs[z] == y, :].X > 0).sum(axis=0)
            if len(sub[:, (spec_nonzero_values == 0) & (global_nonzero_values != 0)].var_names) > 0:
                    print(str(datetime.now()) + " -- Adding 3 pseudocounts to " + z + "=" + str(y) + " for " + str(len(sub[:, (spec_nonzero_values == 0) & (global_nonzero_values != 0)].var_names)) + " genes in " + outgroup.replace("/", " ") + " because all values are 0 for that covariate", flush=True)
                    for x in sub[:, (spec_nonzero_values == 0) & (global_nonzero_values != 0)].var_names:
                        random_three = np.random.choice(sub[sub.obs[z] == y].obs_names, 3)
                        if layer != "X":
                            sub.layers[layer][[i in random_three for i in sub.obs_names], np.where(sub.var_names == x)[0][0]] = 1
                        else:
                            sub.X[[i in random_three for i in sub.obs_names], np.where(sub.var_names == x)[0][0]] = 1
    
    covariates = np.setdiff1d(covariates, to_remove).tolist()

    if layer != "X":
        tmp = sub.layers[layer].T.tocoo()
    else:
        tmp = sub.X.T.tocoo()
        
    counts = r_Matrix.sparseMatrix(
        i=ro.IntVector(tmp.row + 1),
        j=ro.IntVector(tmp.col + 1),
        x=ro.FloatVector(tmp.data),
        dims=ro.IntVector(tmp.shape)
    )
    with localconverter(ro.default_converter + pandas2ri.converter):
        obs = ro.conversion.py2rpy(sub.obs)
    
    ro.r.assign("obs", obs)
    ro.r.assign("counts", counts)
    
    if tests == [""]:
        return
    
    for i in tests:
        if os.path.exists(os.path.join(pwd, "output", region, i)) == False:
            os.makedirs(os.path.join(pwd, "output", region, i))
            
        if os.path.exists(os.path.join(pwd, "output", region, i, target.replace("/", " ") + "_" + outgroup.replace("/", " ") + "_across_" + i + "_DE.csv")) == False:
            print(str(datetime.now()) + " -- Starting " + outgroup.replace("/", " ") + " across " + i)
            ro.r("covariates <- c('" + i + "', " + str(covariates).replace("[", "").replace("]", "") + ")")
            ro.r("df <- model.matrix(~" + covariate_formula + i + ", data=obs[,covariates])")
            ro.r("data_g <- group_cell(count=counts, id=obs$" + random_effect + ", offset=obs$" + offset_variable + ", pred=df)")

            with localconverter(ro.default_converter + pandas2ri.converter):
                df = ro.r("as.data.frame(df)")
            
            ordered = ro.r("data_g")
            if isinstance(ordered, rpy2.rinterface_lib.sexp.NULLType) == True:
                ro.r("re <- nebula(counts, obs$" + random_effect + ", pred=df, offset=obs$" + offset_variable + ", covariance=TRUE)")
                ro.r("saveRDS(re, '" + os.path.join(pwd, "output", region, i, target.replace("/", " ") + "_" + outgroup.replace("/", " ") + "_across_" + i + "_DE.rds") + "')")

            else:
                ro.r("re <- nebula(data_g$count, data_g$id, pred=data_g$pred, offset=data_g$offset, covariance=TRUE)")
                ro.r("saveRDS(re, '" + os.path.join(pwd, "output", region, i, target.replace("/", " ") + "_" + outgroup.replace("/", " ") + "_across_" + i + "_DE.rds") + "')")

            with localconverter(ro.default_converter + pandas2ri.converter):
                results = ro.r("re$summary")
                overdispersion = ro.r("re$overdispersion")
                convergence = ro.r("re$convergence")
                random_effect = ro.r("re$random_effect")

            results = results.loc[convergence == 1, :]
            gene_index = [j - 1 for j in results["gene_id"].to_list()]
            results.index = sub.var_names[gene_index]
            results.to_csv(os.path.join(pwd, "output", region, i, target.replace("/", " ") + "_" + outgroup.replace("/", " ") + "_across_" + i + "_DE.csv"))
            print(str(datetime.now()) + " -- " + outgroup.replace("/", " ") + " along " + i + " was written to disk", flush=True)
        
        else:
            print(str(datetime.now()) + " -- Skipping " + outgroup.replace("/", " ") + " along " + i + " (already exists)", flush=True)
    return

threads = args.threads
threads = int(threads)
target = args.target
region = args.region
split_key = args.split_key
covariates = args.covariates
covariates = [a for a in covariates.split(",")]
covariate_formula = args.covariate_formula
tests = args.tests
tests = [a for a in tests.split(",")]
layer = args.layer
adata = sc.read_h5ad(os.path.join(pwd, "tmp", region, target.replace("/", " ") + ".h5ad"))
adata.obs[split_key] = adata.obs[split_key].astype("category")
offset_variable = args.offset_variable

Parallel(n_jobs=threads)(
    delayed(run_de_tests)(
        adata, region, target, split_key, b, covariates, covariate_formula, tests, layer, offset_variable) for b in adata.obs[split_key].cat.categories
)