import pandas as pd

# File paths
gene_input = "/home/sr933/rcc/data/Target_genes.csv"
node_output = "../data/key_proteins.txt"

# Read the data
gene_df = pd.read_csv(gene_input, header=None)
node_df = pd.read_csv(node_output)

# Find common items in the first column
common_items = gene_df.iloc[:, 0][gene_df.iloc[:, 0].isin(node_df.iloc[:, 0])]

# Print the common items
print("Common items:")
print(len(common_items.tolist()))

for gene in common_items.tolist():
    print(gene)
