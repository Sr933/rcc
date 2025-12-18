import pandas as pd

# Load the original CSV file
file_path = '/home/sr933/rcc/4_network_analysis/src/PubChem_compound_id_list.csv'  # Replace with your file path
data = pd.read_csv(file_path)

# Extract the desired columns
columns_to_extract = ['cmpdname']
filtered_data = data[columns_to_extract]
filtered_data = filtered_data.drop_duplicates(subset=['cmpdname'], keep='first')
# Save the filtered data to a new CSV file
output_path = 'all_drug_names_only.csv'  # Replace with your desired output file path
filtered_data.to_csv(output_path, index=False)

print(f"Extracted columns saved to {output_path}")
