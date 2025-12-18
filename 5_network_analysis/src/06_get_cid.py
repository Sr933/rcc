import pandas as pd

# Load the CSV file
file_path = '/home/sr933/rcc/4_network_analysis/src/manual_results_all.csv'  # Replace with the path to your CSV file
data = pd.read_csv(file_path)

# Extract the 'Drug' column
drug_ids = data['Drug']

# Filter numeric drug IDs only
numeric_drug_ids = [drug_id[4:] for drug_id in drug_ids if drug_id.startswith("CID") and drug_id[4:].isdigit()]
numeric_drug_ids = list(dict.fromkeys(numeric_drug_ids)) # Remove duplicates
print(len(numeric_drug_ids))
# Print the result as a comma-separated list
if numeric_drug_ids:
    print(','.join(numeric_drug_ids))
else:
    print("No numeric drug IDs found.")
