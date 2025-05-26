# Install required libraries if not already installed
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

# Load libraries
library(readr)
library(Seurat)
library(dplyr)

Sys.setenv("VROOM_CONNECTION_SIZE" = 500000)  # Adjust the size as needed
# Read the gzipped CSV file
# Specify the directory containing the files

file="C:/Users/silas/Downloads/GSE156632_RAW/GSM4735375_RCC7t.csv.gz"
count_matrix <- read_csv(file)

# Optionally, you might want to process the count_matrix (e.g., ensuring consistent column names or row names)
count_matrix <- as.data.frame(count_matrix)
count_matrix <- count_matrix[!duplicated(count_matrix$Symbol), ]
# Set rownames to the 'Gene_ID' column
rownames(count_matrix) <- count_matrix$Symbol

count_matrix <- count_matrix[, -c(1, 2)] # Remove the gene column
# Add the count_matrix to the list




# Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = count_matrix)
seurat_obj <- NormalizeData(seurat_obj)

# Define marker genes for each cell type
markers <- list(
  Tumor = c("CRYAB", "CD24", "NDUFA4L2", "NNMT", "RARRES2"),  # Example markers
  Tcell = c("NKG7", "CCL5", "GNLY", "CD52", "KLRB1"),
  Epithelial = c("MT1G", "ALDOB", "MT1H", "GPX3", "MIOX"),
  Endothelial = c("PLVAP", "SPRY1", "SPARCL1", "IGFBP3", "VWF"),
  Myeloid = c("APOC1", "C1QB", "CCL3","CCL3L3", "HLA-DRA"),
  Fibroblast = c("RGS5", "TAGLN", "MGP", "MYL9"),
  Bcell = c("IGLC3", "IGHG3", "IGKC", "IGLC2", "IGHA1")
)

for (cell_type in names(markers)) {
  missing_genes <- markers[[cell_type]][!markers[[cell_type]] %in% rownames(seurat_obj)]
  if (length(missing_genes) > 0) {
    warning(paste("The following features are not present for", cell_type, ":", paste(missing_genes, collapse = ", ")))
  }
}

for (cell_type in names(markers)) {
  seurat_obj <- AddModuleScore(seurat_obj, features = list(markers[[cell_type]]), name = cell_type)
}


# Extract scores
scores <- seurat_obj@meta.data[, grep("^[a-zA-Z]+1$", colnames(seurat_obj@meta.data))]

# Assign the cell type with the highest score
seurat_obj$cell_type <- apply(scores, 1, function(row) names(row)[which.max(row)])


# Export the annotated cell metadata
annotated_metadata <- seurat_obj@meta.data
write.csv(annotated_metadata, 
          paste0(file, "cell_annotation.csv"), 
          row.names = TRUE)

