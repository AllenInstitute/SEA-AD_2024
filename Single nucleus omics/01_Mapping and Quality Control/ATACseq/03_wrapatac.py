import sys
sys.path.insert(0, '/allen/programs/celltypes/workgroups/hct/SEA-AD/Integration/python_supp_scripts/mvi_atac_processor')
from wlib import *


ad = sc.read_h5ad('pre_pipeline_atac_Astro.h5ad')

required_fields = {
    'barcodes': 'bc_pipeline_atac',
    'library_id': 'libraryprep_arid',
    'modality': 'modality_pipeline', 
    'path_to_fragment': 'path_to_fragments',
    'bed': 'path_to_bed',
    'group_by': 'modality_pipeline',
    'split_by': 'label_transfer',
}
selected_modality = ['paired', 'atac']

job_id = 'Astro'
temp_dir = '/allen/programs/celltypes/workgroups/hct/SEA-AD/Integration/'
temp_dir += 'multivi_subclasses/Astro/'


## checkpoint path
## working dir
## temp dir

#exps = ['L8AT_220209_01_B05-1160415286','L8AT_220316_01_H10-1167718971', 'L8XR_220210_02_D05-1160243206']

#tmp = ad[ad.obs.exp_id.isin(exps)].copy()

ad_atac = wrapper_build_atac(ad, required_fields, temp_dir, job_id, selected_modality)

ad_atac = sc.read_h5ad(ad_atac)


ad_rna = sc.read_h5ad('pre_pipeline_rna_Astro.h5ad')
ad_atac.var['modality'] = 'Peaks'

ad_rna.obs = ad_rna.obs.set_index('sample_id')

merged_path = aggregate_2_modalities(ad_rna, ad_atac, f'multi_{job_id}.h5ad', force=True)
merged = sc.read_h5ad(merged_path)

