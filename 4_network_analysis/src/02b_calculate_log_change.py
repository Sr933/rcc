import numpy as np
import pandas as pd
import pickle

# Define the path to the .pkl file
input_path = "/data/sr933/scRCC/combined_data/RCC_data_dict.pkl"

# Load the dictionary from the .pkl file
with open(input_path, "rb") as f:
    data_dict = pickle.load(f)

# Access the contents of the dictionary
X_combined = data_dict["X"].T
y_labels = data_dict["y"]
gene_list = data_dict["Genes"]
# Separate the expression data and labels
X_combined = data_dict["X"].T  # Genes as rows, samples as columns
y_labels = data_dict["y"]      # Labels
gene_list = data_dict["Genes"] # Gene names

# Convert to DataFrame for easier manipulation
X_df = pd.DataFrame(X_combined.T, index=gene_list)
y_series = pd.Series(y_labels)

# Filter samples based on y=0 and y≠0
X_y0 = X_df.loc[:, y_series == 0]  # Samples where y=0
X_other = X_df.loc[:, y_series != 0]  # Samples where y≠0

# Calculate mean expression for each gene in each group
mean_y0 = X_y0.mean(axis=1)
mean_other = X_other.mean(axis=1)

# Calculate log2 fold change
log2_fold_change = np.log2((mean_y0 + 1e-8) / (mean_other + 1e-8))  # Add a small value to avoid division by zero

# Store results in a DataFrame
log2_fc_df = pd.DataFrame({
    "Gene": gene_list,
    "log2FoldChange": log2_fold_change
})

# Sort by absolute log2 fold change (optional)
#log2_fc_df["abs_log2FoldChange"] = log2_fc_df["log2FoldChange"].abs()
#log2_fc_df = log2_fc_df.sort_values(by="abs_log2FoldChange", ascending=False)
log2_fc_df.to_csv("/home/sr933/rcc/4_network_analysis/data/log2_fold_change.csv", index=False)
