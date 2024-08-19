import numpy as np
import pandas as pd
import scanpy as sc

def create_atac_anndata(base_path, library_prep_ar_id, mtx_ext, peak_ext, barcode_ext):
    adata = sc.read_mtx(os.path.join(base_path, library_prep_ar_id + mtx_ext))
    coords = pd.read_csv(
        os.path.join(base_path, library_prep_ar_id + peak_ext),
        sep="\t",
        header=None,
        index_col=None,
    )
    coords.rename({0: "chr", 1: "start", 2: "end"}, axis="columns", inplace=True)
    coords.set_index(
        coords.chr.astype(str)
        + ":"
        + coords.start.astype(str)
        + "-"
        + coords.end.astype(str),
        inplace=True,
    )
    coords.index = coords.index.astype(str)
    
    cell_annot = pd.read_csv(
        os.path.join(base_path, library_prep_ar_id + barcode_ext), 
        sep="-", 
        header=None, 
        index_col=None
    )
    cell_annot.rename({0: "barcode", 1: "batch_id"}, axis="columns", inplace=True)
    cell_annot["library_prep"] = library_prep_ar_id.split("-")[0]
    cell_annot["sample_id"] = cell_annot["barcode"] + "-" + library_prep_ar_id
    cell_annot["barcode"] = cell_annot["barcode"] + "-" + cell_annot["library_prep"]
    cell_annot.set_index("barcode", inplace=True)
    cell_annot.index = cell_annot.index.astype(str)
    
    adata.obs = cell_annot
    adata.var = coords
    adata.var["modality"] = "Peaks"
    return adata.copy()
