#tracker = read.csv('/allen/programs/celltypes/workgroups/hct/emilyg/resegmented_pfc_data/new_mtg_tracker.csv', sep=',')
#path_list <- as.list(tracker$isilon.location) #grab all segmented file paths
#path_list = path_list[nzchar(path_list)] #remove empty char paths (failed files, not-yet-run, etc)
base_path <- '/allen/programs/celltypes/workgroups/hct/SEA-AD/MERSCOPE/resegmentation_and_corr_mapping_data/mtg_manuscript/'
path_list <- list.files(base_path)

min_genes <- 3
min_total_reads <- 30
min_vol <- 100

upper_bound_reads <- 4000
upper_genes_read <- 130
## Compute log CPM
for (cpath in path_list) {
  cpath <- paste0(base_path, cpath)
  test_path <- paste0(cpath, "/cbg_cpum.csv")
  print(test_path)
  if (!file.exists(test_path)) {
    cbg_ops <- list.files(path = cpath, pattern = "cellpose-cell-by-gene.csv", recursive = TRUE)
    cbg <- data.table::fread(paste0(cpath, '/', cbg_ops[1]), 
                             header=TRUE,
                             stringsAsFactors=FALSE,
                             check.names = FALSE
    )
    # turn the first row (containing cell id's) to rownames
    cbg <- tibble::column_to_rownames(cbg,var="cell")
    
    # load metadata for associated cell by gene table
    md_ops <- list.files(path = cpath, pattern = "cellpose_metadata.csv", recursive = TRUE)
    metadata <- data.table::fread(paste0(cpath, '/', md_ops[1]),
                                  header=TRUE,
                                  stringsAsFactors=FALSE,
                                  check.names = FALSE
    )
    # turn the first row (containing cell id's) to rownames
    metadata <- tibble::column_to_rownames(metadata,var="id")
    
    metadata <- metadata[match(rownames(cbg), rownames(metadata)),]
    cbg <- cbg[match(rownames(metadata),rownames(cbg)),]
    
    blanks <- dplyr::select(cbg, contains("Blank"))
    cbg <- dplyr::select(cbg,-contains("Blank"))
    
    metadata$genes_detected <- rowSums(cbg != 0)
    metadata$total_reads <- rowSums(cbg)
    
    upper_bound_area <- 3 * (median(metadata$volume[metadata$volume > 100]))
    
    metadata <- metadata %>%
      mutate(
        cell_qc = if_else(
          genes_detected < min_genes |
            total_reads < min_total_reads |
            volume < min_vol |
            volume > upper_bound_area |
            genes_detected > upper_genes_read |
            total_reads > upper_bound_reads,
          "Low",
          "High"
        )
      )
    
    metadata$min_genes <- min_genes
    metadata$min_total_reads <- min_total_reads
    metadata$min_vol <- min_vol
    metadata$upper_bound_area <- upper_bound_area
    metadata$upper_bound_reads <- upper_bound_reads
    metadata$upper_genes_read <- upper_genes_read
    
    # combine blanks and cbg back together
    cbg <- merge(cbg,blanks,by=0)
    # turn the first row (containing cell id's) to rownames
    cbg <- tibble::column_to_rownames(cbg,var="Row.names")
    # order rows in metadata file to match the one in cell by gene table
    cbg <- cbg[match(rownames(metadata),rownames(cbg)),]
    
    # filter metadata to only contain high quality cells
    metadata_subset <- metadata %>% dplyr::filter(cell_qc == "High") 
    # filter cell by gene to only contain high quality cells
    to_keep_org <- intersect(rownames(cbg),rownames(metadata_subset))
    cbg_filtered <- cbg[to_keep_org,]
    cbg_filtered <- as.matrix(cbg_filtered)
    # normalize cells by volume (counts per ??m)
    cbg_cpum <- (cbg_filtered / metadata_subset$volume*1000)
    # log normalize data
    cbg_norm <- log2(cbg_cpum+1)
    # order rows between cell by genes tables and metadata to match each other
    cbg_norm <- cbg_norm[match(rownames(cbg_cpum),rownames(cbg_norm)),]
    cbg_filtered <- cbg_filtered[match(rownames(cbg_cpum),rownames(cbg_filtered)),]

    #write each version that will go into final anndata
    write.csv(cbg_norm, paste0(cpath,"/cbg_cpum.csv"), row.names=TRUE)
	write.csv(t(cbg_filtered), paste0(cpath,"/cbg_filtered_T.csv"), row.names=TRUE) #just this file used in mapping
	write.csv(cbg_filtered, paste0(cpath,"/cbg_filtered.csv"), row.names=TRUE) 
    write.csv(metadata, paste0(cpath, "/metadata_updated.csv"), row.names=TRUE)
    print(paste0(cpath,"/cbg_filtered_T.csv"))
  }
}