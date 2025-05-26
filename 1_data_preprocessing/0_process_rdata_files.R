# Load Seurat library
.libPaths(c("~/R/libs", .libPaths()))
library(Seurat)
obj <- readRDS("/data/sr933/scRCC/Tcell_CD8_other_ha.1.RDS")
print(class(obj))

# Assume `seurat_obj` is your SeuratObject
# Extract the expression matrix (e.g., raw counts)
expression_matrix <- GetAssayData(object = obj, slot = "counts")


# Convert to a data.frame
expression_df <- as.data.frame(as.matrix(expression_matrix))

# Save as a CSV file
saveRDS(expression_df, file = "/data/sr933/scRCC/CD8_cells_expression_matrix.rds")