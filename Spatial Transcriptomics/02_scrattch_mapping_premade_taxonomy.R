## Load scrattch.mapping - be sure to run on 0.16 as loading anndata breaks on 0.41
library(scrattch.mapping)

#grab all files with a properly filtered and transposed cell by gene table.
mpath <- "/allen/programs/celltypes/workgroups/hct/SEA-AD/MERSCOPE/resegmentation_and_corr_mapping_data/mtg_manuscript/"
path_list <- list.files(path=mpath ,pattern="cbg_cpum_T.csv", recursive=TRUE) 

for (cpath in path_list) {
  #skip any sections that have already been mapped
  full_path <- file.path(mpath, cpath)
  raw_path <- sub('cbg_filtered_T.csv','', full_path)
  map_list <- list.files(path=raw_path, pattern='*_mapped.csv', recursive = TRUE) 
  map_path <- file.path(raw_path, map_list[1])
  if (!file.exists(map_path)){
	print(raw_path)
	to_map_data = read.csv(full_path, sep=",", row.names=1)
	query.data <- as.matrix(to_map_data) #counts matrix of data to be annotated 
	
	#read in taxonomy anndata
	AIT.anndata = read_h5ad('//allen/programs/celltypes/workgroups/hct/SEA-AD/MERSCOPE/resegmentation_and_corr_mapping_data/mtg_taxonomy_spatial_gene_subset.h5ad')
	
	# Map! Returns a data.frame with mapping results.
	mapping = taxonomy_mapping(AIT.anndata=AIT.anndata,
							   query.data=query.data,
							   label.cols="supertype_scANVI_leiden", "subclass_label", "class_label"), 
	 							corr.map=TRUE,
								tree.map=FALSE,
								seurat.map=FALSE)
	#save out results
	res_path <- sub('_T.csv','_mapped.csv', full_path)
	print('writing')
	write.csv(mapping, res_path)
  }
}
  