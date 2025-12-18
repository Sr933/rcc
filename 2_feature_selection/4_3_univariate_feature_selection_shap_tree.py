import pickle
from feature_engine.selection import SelectBySingleFeaturePerformance
from lightgbm import LGBMClassifier
import pandas as pd
from imblearn.over_sampling import SMOTE
from collections import Counter
import numpy as np

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

print("Starting feature selection")

model = LGBMClassifier(
    n_estimators=1000,       # 1000 trees
    learning_rate=0.05,      # Small learning rate
    max_depth=10,            # Limit tree depth to 10
    num_leaves=31,           # Number of leaves per tree
)
# ------------------ New: 5-fold CV LightGBM + SelectFromModel + SHAP ------------------
from sklearn.model_selection import StratifiedKFold
from lightgbm import LGBMClassifier
from sklearn.feature_selection import SelectFromModel
import shap

print("Starting 5-fold CV training with LightGBM (no pre-selection)...")

skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

topk = 100
per_fold_selected = {}
per_fold_shap = {}

for fold, (train_idx, valid_idx) in enumerate(skf.split(X_train, y_train), start=1):
    print(f"Fold {fold}: training LightGBM...")
    X_tr, X_val = X_train[train_idx], X_train[valid_idx]
    y_tr, y_val = y_train[train_idx], y_train[valid_idx]

    model = LGBMClassifier(
    n_estimators=1000,       # 1000 trees
    learning_rate=0.05,      # Small learning rate
    max_depth=10,            # Limit tree depth to 10
    num_leaves=31,           # Number of leaves per tree
)
    model.fit(X_tr, y_tr, eval_set=[(X_val, y_val)], eval_metric='auc')

    # SelectFromModel using feature importances
    sfm = SelectFromModel(model, prefit=True, threshold='mean')
    mask = sfm.get_support()
    importances = model.feature_importances_

    # Rank by model importance and take topk
    idx_sorted = np.argsort(importances)[::-1]
    topk_idx_importance = idx_sorted[:topk]
    topk_features_importance = [gene_list[i] for i in topk_idx_importance]
    per_fold_selected[f'fold_{fold}'] = topk_features_importance

    # Compute SHAP values (TreeExplainer)
    try:
        explainer = shap.TreeExplainer(model)
        shap_vals = explainer.shap_values(X_tr)
        if isinstance(shap_vals, list):
            shap_vals = shap_vals[1]
    except Exception as e:
        print(f"SHAP TreeExplainer failed on fold {fold}: {e}. Trying KernelExplainer (slower)")
        explainer = shap.KernelExplainer(model.predict_proba, X_tr[:100])
        shap_vals = explainer.shap_values(X_tr, nsamples=100)
        if isinstance(shap_vals, list):
            shap_vals = shap_vals[1]

    mean_abs_shap = np.mean(np.abs(shap_vals), axis=0)
    topk_idx_shap = np.argsort(mean_abs_shap)[::-1][:topk]
    topk_features_shap = [gene_list[i] for i in topk_idx_shap]
    per_fold_shap[f'fold_{fold}'] = topk_features_shap

    # Save per-fold CSVs
    pd.DataFrame({'feature_importance_rank': [gene_list[i] for i in idx_sorted],
                  'importance': importances[idx_sorted]}).head(topk).to_csv(f'/home/sr933/rcc_data/data/fold_{fold}_lgbm_importance_top{topk}.csv', index=False)

    pd.DataFrame({'feature': topk_features_shap, 'mean_abs_shap': mean_abs_shap[topk_idx_shap]}).to_csv(f'/home/sr933/rcc_data/data/fold_{fold}_shap_top{topk}.csv', index=False)

    print(f"Fold {fold} done: saved top-{topk} by importance and SHAP")

# Aggregate across folds: simple frequency-based ranking
from collections import Counter

def aggregate_topk(dict_of_lists, topk=100):
    all_items = []
    for v in dict_of_lists.values():
        all_items.extend(v)
    counter = Counter(all_items)
    most_common = [f for f, _ in counter.most_common(topk)]
    return most_common, counter

agg_imp, counter_imp = aggregate_topk(per_fold_selected, topk=topk)
agg_shap, counter_shap = aggregate_topk(per_fold_shap, topk=topk)

pd.DataFrame({'feature': agg_imp}).to_csv(f'/home/sr933/rcc_data/data/aggregated_lgbm_top{topk}_by_freq.csv', index=False)
pd.DataFrame({'feature': agg_shap}).to_csv(f'/home/sr933/rcc_data/data/aggregated_shap_top{topk}_by_freq.csv', index=False)

print(f"Saved aggregated top-{topk} feature lists to /home/sr933/rcc_data/data/")
