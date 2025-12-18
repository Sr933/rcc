import pandas as pd
import ast

# Function to check if genes match the specified criteria
def check_genes_in_targets(row, gene_list_df):
    # Convert string representations of lists to actual lists
    degree_0_genes = ast.literal_eval(row['0degree.target.keyprotein'])
    degree_1_genes = ast.literal_eval(row['1degree.target.keyprotein'])

    # Flatten the lists in case they are nested
    degree_0_genes_flat = [gene for sublist in degree_0_genes for gene in (sublist if isinstance(sublist, list) else [sublist])]
    degree_1_genes_flat = [gene for sublist in degree_1_genes for gene in (sublist if isinstance(sublist, list) else [sublist])]

    if len(degree_0_genes_flat) >= 10:
        return False

    # Check against the gene list with "Druggable" criteria
    for _, gene_row in gene_list_df.iterrows():
        gene = gene_row['Gene']
        #druggable = gene_row['Druggable']

        if gene in degree_0_genes_flat or gene in degree_1_genes_flat:
            return True
        
        

    return False

# Load the gene list CSV
# Assuming the file has columns "Gene" and "Druggable" (with values "Yes" or "No")
gene_list_df = pd.read_csv('/home/sr933/rcc/4_network_analysis/data/key_proteins_crispr.txt')

# List of CSVs to load
csv_files = [
    #'/home/sr933/rcc/4_network_analysis/data/drug_network_proximity_results_all_drugs_reduced.csv',
    #'/home/sr933/rcc/4_network_analysis/data/drug_network_proximity_results_up_down_fda.csv',
    #'/home/sr933/rcc/4_network_analysis/data/drug_network_proximity_results_up_fda_no_common_essentials.csv'
    #'/home/sr933/rcc/4_network_analysis/data/drug_network_proximity_results_up_down_all.csv'
    #'/home/sr933/rcc/4_network_analysis/data/drug_network_proximity_results_manually_filtered.csv'
    '/home/sr933/rcc/4_network_analysis/data/drug_network_proximity_results_manually_filtered_fda.csv'
]

# Initialize an empty DataFrame to store the concatenated results
final_df = pd.DataFrame()

# Loop through each CSV file
for file in csv_files:
    df = pd.read_csv(file)

    # Filter out rows where 'z' >= -2 before applying the gene check function
    df_filtered = df[df['z'] <= -2]

    # Apply the function to check for genes in degree 0 or degree 1 targets
    df_filtered = df_filtered[df_filtered.apply(lambda row: check_genes_in_targets(row, gene_list_df), axis=1)]

    # Concatenate the filtered rows to the final dataframe
    final_df = pd.concat([final_df, df_filtered])

# Reset index for the final dataframe
final_df.reset_index(drop=True, inplace=True)
final_df_sorted = final_df.sort_values(by='z', ascending=True)
final_df_sorted = final_df_sorted.drop_duplicates(subset=['0degree.target.keyprotein', '1degree.target.keyprotein',  '2degree.target.keyprotein'], keep='first')

# Export the final DataFrame to a CSV
final_df_sorted.to_csv('manual_results_all.csv', index=False)
