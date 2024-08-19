suppressPackageStartupMessages(
  {
    library(ArchR)
    library(stringr)
    library(dplyr)
  }
)

addArchRGenome("hg38")

output_directory <- "MTG_AD_Non-neuronal"

dir <- "input/"
outdir <- "output/"

inputf_fls <- paste0(dir, "ATAC_AD_Center_Grant_update.csv")
inputf <- read.csv(inputf_fls, header = T)

inputf$arrow_path <- paste0(dir, "ArrowFiles/",inputf$arrow_file)

barcode.files <- paste0(dir, "Non-neuronal_Supertype_Annotation_with_ArchR_index.csv")
barcodes_anno <- read.csv(barcode.files, header = TRUE)
cell_barcodes <- barcodes_anno$archr_index

ArrowFiles <- inputf$arrow_path

# Build ArchR file from pre-made Arrow files
proj1 <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = output_directory,
  copyArrows = FALSE
)

# Subset ArchR file on Non-neuronal cells
barcode_proj1 <- proj1$cellNames
barcode.shared <- intersect(barcode_proj1, cell_barcodes)

proj1 <- subsetArchRProject(
  ArchRProj = proj1,
  cells = barcode.shared,
  outputDirectory = output_directory,
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)

## Add barcodes_anno metadata into the ArchR object.
meta <- as.data.frame(getCellColData(proj1))

meta <- meta %>%
  rownames_to_column('sample_plus_barcode') %>%
  mutate(archr_index = sample_plus_barcode) %>%
  left_join(barcodes_anno, by = c("archr_index" = "archr_index")) %>%
  column_to_rownames('sample_plus_barcode')

meta <- as(meta, "DataFrame")
proj1@cellColData <- meta

saveArchRProject(ArchRProj = proj1, outputDirectory = output_directory, load = FALSE)


# Call peaks with MACS2
proj1 <- addGroupCoverages(ArchRProj = proj1, 
                           groupBy = "Supertype",
                           force = TRUE
                           )

pathToMacs2 <- findMacs2()
proj1 <- addReproduciblePeakSet(
  ArchRProj = proj1, 
  groupBy = "Supertype", 
  pathToMacs2 = pathToMacs2,
  force = TRUE
)

proj1 <- addPeakMatrix(ArchRProj = proj1, force = TRUE)

saveArchRProject(ArchRProj = proj1, outputDirectory = output_directory, load = FALSE)

## Test pycistopic object
peak_matrix <- getMatrixFromProject(ArchRProj = proj1, useMatrix = "PeakMatrix")

# Access assay data
assay_data <- assays(peak_matrix)$PeakMatrix

# Save sparse matrix in Matrix Market format
writeMM(assay_data, paste0(outdir, "sparse_matrix.mtx"))

# Assuming you have rowRanges(peak_matrix)
row_ranges <- rowRanges(peak_matrix)

# Extract seqnames and ranges
seqnames <- seqnames(row_ranges)
ranges <- ranges(row_ranges)

# Create a list with "seqnames:ranges" format
seqnames_ranges_list <- paste(seqnames, ":", ranges, sep = "")

# Print the list
print(seqnames_ranges_list)

# Convert the list to a data frame with a single column
df <- data.frame(Seqnames_Ranges = unlist(seqnames_ranges_list))

# Save the data frame to a CSV file
write.csv(df, paste0(outdir, "seqnames_ranges.csv"), row.names = FALSE)

cell_names <- colnames(peak_matrix)
df <- data.frame(cells = unlist(cell_names))
write.csv(df, paste0(outdir, "cell_names.csv"), row.names = FALSE)

# Save cell metadata to a CSV file.
meta <- as.data.frame(getCellColData(proj1))
write.csv(meta, paste0(outdir, "cell_metadata.csv"))

