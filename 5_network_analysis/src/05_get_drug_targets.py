import pandas as pd
import ast
# Load the two CSV files

def get_targets(file_path, out_path):
    df1 = pd.read_csv(file_path)
    #df2 = pd.read_csv(file_path2)

    # Filter rows where z <= -2 in both DataFrames
    filtered_df1 = df1[df1['z'] <= -2.5]
    #filtered_df2 = df2[df2['z'] <= -2]
    filtered_df1=filtered_df1.sort_values(by='z', ascending=True)
    #filtered_df1 = filtered_df1.drop_duplicates(subset=['0degree.target.keyprotein', '1degree.target.keyprotein',  '2degree.target.keyprotein'], keep='first')
    #final_df_sorted = filtered_df1.drop_duplicates(subset=['0degree.target.keyprotein', '1degree.target.keyprotein',  '2degree.target.keyprotein'], keep='first')
    filtered_df1.to_csv(out_path, index=False)

file_path_all = "/home/sr933/rcc/4_network_analysis/data/drug_network_proximity_results_up_down_all_no_common_essentials.csv"  # Replace with the path to the first CSV file
out_all= "filtered_drugs_all.csv"

get_targets(file_path_all, out_all)

file_path_fda = "/home/sr933/rcc/4_network_analysis/data/drug_network_proximity_results_up_down_fda_no_common_essentials.csv"  # Replace with the path to the second CSV file
out_all= "filtered_drugs_fda.csv"

get_targets(file_path_fda, out_all)
