# notebooks for processing and analysis of spatial data for the SEA_AD project

# Initial data processing
	- 01_map_filtering.R: filters the segmented MERSCOPE results to remove poor quality cells before mapping process begins. Additionally converts the results to a scrattch-mapping readable format.
	- 02_scattch_mapping_premade_taxonomy.R: Runs scrattch-mapping script on the filtered segmented MERSCOPE data using the SEA-AD taxonomy, and returns a csv of mapped results for each file.
	- 03_create_complete_h5ad: combines segmentation and mapping results to create initial anndata objects for each processed section, then combines said anndata objects and metadata information about each into a final combined anndata object.
