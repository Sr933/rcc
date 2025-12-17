import pickle
from collections import Counter

import numpy as np
import pandas as pd
from imblearn.over_sampling import SMOTE
from lightgbm import LGBMClassifier
from sklearn.base import clone
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

# Define the path to the .pkl file
input_path = "/home/sr933/rcc_data/scRCC/combined_data/RCC_data_dict.pkl"

# Load the dictionary from the .pkl file
with open(input_path, "rb") as f:
    data_dict = pickle.load(f)

# Access the contents of the dictionary
X_combined = data_dict["X"].T
y_labels = data_dict["y"]
gene_list = data_dict["Genes"]

# Set all elements >= 1 to 1
y_labels[y_labels >= 1] = 1


##Do SMOTE

def data_smote(x,y):
    oversample = SMOTE()
    x_new, y_new = oversample.fit_resample(x, y)

    counter_old = Counter(y)
    print(counter_old)
    counter = Counter(y_new)
    print(counter)

    #Make features binary
    x_new=np.rint(x_new)
    return x_new, y_new

X_train, y_train= data_smote(X_combined, y_labels)

# Convert arrays to pandas structures for easier column handling
X_train_df = pd.DataFrame(X_train, columns=gene_list)
y_train_series = pd.Series(y_train)

print("Starting feature selection")

base_model = LGBMClassifier(
    n_estimators=1000,       # 1000 trees
    learning_rate=0.05,      # Small learning rate
    max_depth=10,            # Limit tree depth to 10
    num_leaves=31,           # Number of leaves per tree
)

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
fold_splits = list(cv.split(X_train_df, y_train_series))
num_folds = len(fold_splits)

high_auc_genes_by_fold = {fold_idx: set() for fold_idx in range(num_folds)}
results = []

for gene in gene_list:
    fold_aucs = []
    for fold_idx, (train_idx, test_idx) in enumerate(fold_splits):
        model = clone(base_model)
        X_train_fold = X_train_df.iloc[train_idx][[gene]]
        y_train_fold = y_train_series.iloc[train_idx]
        X_test_fold = X_train_df.iloc[test_idx][[gene]]
        y_test_fold = y_train_series.iloc[test_idx]

        model.fit(X_train_fold, y_train_fold)
        y_pred_proba = model.predict_proba(X_test_fold)[:, 1]
        auc = roc_auc_score(y_test_fold, y_pred_proba)
        fold_aucs.append(auc)

        if auc > 0.8:
            high_auc_genes_by_fold[fold_idx].add(gene)

    mean_auc = float(np.mean(fold_aucs))
    std_auc = float(np.std(fold_aucs, ddof=1)) if len(fold_aucs) > 1 else 0.0

    gene_result = {
        "Gene": gene,
        "Mean_AUC": mean_auc,
        "Std_AUC": std_auc,
    }
    for fold_idx, auc in enumerate(fold_aucs, start=1):
        gene_result[f"Fold{fold_idx}_AUC"] = auc

    results.append(gene_result)

results_df = pd.DataFrame(results).set_index("Gene")
results_df.to_csv("/home/sr933/rcc/data/LightGBM_Features.csv")

print("Feature selection results with mean and std AUC saved to LightGBM_Features.csv.")

# Prepare high AUC (>0.8) summaries
high_auc_records = []
for fold_idx, genes in high_auc_genes_by_fold.items():
    for gene in sorted(genes):
        high_auc_records.append({"Fold": fold_idx + 1, "Gene": gene})

high_auc_df = pd.DataFrame(high_auc_records)
high_auc_df.to_csv("/home/sr933/rcc/data/LightGBM_high_auc_genes_per_fold.csv", index=False)

if high_auc_genes_by_fold:
    overlap_genes = set.intersection(*high_auc_genes_by_fold.values()) if all(high_auc_genes_by_fold.values()) else set()
else:
    overlap_genes = set()

overlap_df = pd.DataFrame({"Gene": sorted(overlap_genes)})
overlap_df.to_csv("/home/sr933/rcc/data/LightGBM_high_auc_fold_overlap.csv", index=False)

summary_rows = [
    {"Fold": fold_idx + 1, "Genes_over_0.8": len(genes)}
    for fold_idx, genes in high_auc_genes_by_fold.items()
]
summary_rows.append({"Fold": "Intersection_all_folds", "Genes_over_0.8": len(overlap_genes)})
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv("/home/sr933/rcc/data/LightGBM_high_auc_summary.csv", index=False)

print("High AUC gene summaries saved to data directory.")
