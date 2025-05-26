# Specify the library path
.libPaths(c("~/R/libs", .libPaths()))
lib_path <- "~/R/libs"
Sys.setenv("VROOM_CONNECTION_SIZE" = 500000)
# Load necessary library
library(readr)
dir_path="/data/sr933/scRCC validation/GSE202813_RAW"
# List all files in the directory with the .csv.gz extension and excluding "annotation" in the filename
files <- list.files(path = dir_path, pattern = "\\.csv\\.gz$", full.names = TRUE)
files <- files[!grepl("annotation", files)]  # Exclude files with 'annotation' in the name

# Loop over the filtered files
for (file in files) {

    # Read the compressed file and save as a CSV
data <- read_csv(file)
output_file <- sub("\\.gz$", "", file)
# Write the data to a CSV file
write_csv(data, output_file)

# Print message confirming the operation
cat("File has been successfully saved as:", output_file, "\n")
}




