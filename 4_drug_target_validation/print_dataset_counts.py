#!/usr/bin/env python3
"""
Scan dataset files under this folder and print a small table with counts of tumour vs non-tumour samples.
Heuristics used:
 - For CSV/TSV: look for columns named 'label', 'y', 'tumor', 'tumour', 'class'.
   Accept values like 0/1, 'tumor'/'normal', 'tumour'/'normal', 'non-tumour', 'non_tumour'.
 - For PKL: attempt to load (pickle) and look for keys 'X','y' or dataframe-like objects.
 - For NPZ/NPY: skip unless a companion label file exists.

Outputs per-file counts and an aggregated table.
"""
import os
import pandas as pd
import pickle
from collections import defaultdict

ROOT = "/home/sr933/rcc_data"
DATA_EXTS = ['.csv', '.tsv', '.pkl', '.pickle']

LABEL_COL_CANDIDATES = ['label', 'y', 'target', 'tumor', 'tumour', 'class']


def infer_counts_from_df(df):
    df_cols = [c.lower() for c in df.columns]
    label_col = None
    for cand in LABEL_COL_CANDIDATES:
        if cand in df_cols:
            label_col = df.columns[df_cols.index(cand)]
            break
    if label_col is None:
        return None
    labels = df[label_col]
    # normalize labels to 0/1
    mapping = {}
    unique = pd.Series(labels.dropna().unique()).astype(str).str.lower()
    # Try mapping common values
    for val in unique:
        if val in ['0', 'false', 'neg', 'normal', 'non-tumour', 'nontumour', 'non-tumour', 'non_tumour', 'non tumour', 'non']:
            mapping[val] = 0
        elif val in ['1', 'true', 'pos', 'tumour', 'tumor', 'cancer', 'tumour', 'tumorous']:
            mapping[val] = 1
    if not mapping:
        return None
    norm = labels.astype(str).str.lower().map(mapping).dropna()
    counts = norm.value_counts().to_dict()
    return {'tumour': int(counts.get(1, 0)), 'non_tumour': int(counts.get(0, 0))}


def scan():
    results = []
    agg = defaultdict(int)
    for root, dirs, files in os.walk(ROOT):
        for fn in files:
            path = os.path.join(root, fn)
            _, ext = os.path.splitext(fn.lower())
            if ext not in DATA_EXTS:
                continue
            counts = None
            try:
                if ext in ['.csv', '.tsv']:
                    sep = ',' if ext == '.csv' else '\t'
                    df = pd.read_csv(path, sep=sep)
                    counts = infer_counts_from_df(df)
                elif ext in ['.pkl', '.pickle']:
                    with open(path, 'rb') as f:
                        obj = pickle.load(f)
                    # dict with X,y
                    if isinstance(obj, dict):
                        if 'y' in obj:
                            y = obj['y']
                            # assume binary 0/1 or labels
                            try:
                                s = pd.Series(y)
                                val_counts = s.value_counts().to_dict()
                                tumour = int(val_counts.get(1, 0) + val_counts.get('tumour', 0) + val_counts.get('tumor', 0))
                                non_tumour = int(val_counts.get(0, 0) + val_counts.get('normal', 0))
                                counts = {'tumour': tumour, 'non_tumour': non_tumour}
                            except Exception:
                                counts = None
                        elif hasattr(obj, 'shape') and len(obj.shape) == 2:
                            counts = None
                    elif hasattr(obj, 'columns'):
                        counts = infer_counts_from_df(obj)
            except Exception:
                counts = None

            results.append((path, counts))
            if counts:
                agg['tumour'] += counts['tumour']
                agg['non_tumour'] += counts['non_tumour']
#!/usr/bin/env python3
"""
Print a table with counts of tumour vs non-tumour samples for all datasets
as they are loaded in the notebooks in this folder. This mirrors the exact
loading and label rules used in the notebooks, without performing any
feature selection or model training.
"""
import os
import sys
import pickle
from typing import Dict, Tuple

import pandas as pd


def count_american() -> Tuple[str, Dict[str, int]]:
    name = "American"
    pkl_path = "/home/sr933/rcc_data/scRCC validation/processed data/scRCC_validation_data.pkl"
    if not os.path.exists(pkl_path):
        return name, {"tumour": -1, "non_tumour": -1}
    with open(pkl_path, "rb") as f:
        data = pickle.load(f)
    # X_data is a DataFrame with genes as index and samples as columns
    X_data = data.get("X_data")
    if not isinstance(X_data, pd.DataFrame):
        return name, {"tumour": -1, "non_tumour": -1}
    cols = list(X_data.columns)
    # Notebook rule: y_test = [0 if 'Tumor' in col else 1]
    tumour = sum(1 for c in cols if "Tumor" in c)
    non_tumour = len(cols) - tumour
    return name, {"tumour": tumour, "non_tumour": non_tumour}


def count_chinese() -> Tuple[str, Dict[str, int]]:
    name = "Chinese"
    pkl_path = "/home/sr933/rcc_data/scRCC validation/processed data/scRCC_validation_dataset_china.pkl"
    if not os.path.exists(pkl_path):
        return name, {"tumour": -1, "non_tumour": -1}
    with open(pkl_path, "rb") as f:
        data = pickle.load(f)
    y = data.get("y_test")
    if y is None:
        return name, {"tumour": -1, "non_tumour": -1}
    tumour = sum(1 for v in y if v == "Tumor1")
    non_tumour = len(y) - tumour
    return name, {"tumour": tumour, "non_tumour": non_tumour}


def count_metastasis() -> Tuple[str, Dict[str, int]]:
    name = "Metastasis"
    pkl_path = "/home/sr933/rcc_data/scRCC validation/processed data/scRCC_validation_dataset_metastatic.pkl"
    if not os.path.exists(pkl_path):
        return name, {"tumour": -1, "non_tumour": -1}
    with open(pkl_path, "rb") as f:
        data = pickle.load(f)
    y = data.get("y_test")
    if y is None:
        return name, {"tumour": -1, "non_tumour": -1}
    tumour = sum(1 for v in y if v == "Tumor1")
    non_tumour = len(y) - tumour
    return name, {"tumour": tumour, "non_tumour": non_tumour}

def count_train() -> Tuple[str, Dict[str, int]]:
    name = "Train"
    pkl_path = "/home/sr933/rcc_data/scRCC/combined_data/RCC_data_dict.pkl"
    if not os.path.exists(pkl_path):
        return name, {"tumour": -1, "non_tumour": -1}
    with open(pkl_path, "rb") as f:
        data = pickle.load(f)
    y = data.get("y")
    if y is None:
        return name, {"tumour": -1, "non_tumour": -1}
    tumour = sum(1 for v in y if v == 1)
    non_tumour = len(y) - tumour
    return name, {"tumour": tumour, "non_tumour": non_tumour}


def count_lithuanian() -> Tuple[str, Dict[str, int]]:
    name = "Lithuanian"
    pkl_path = "/home/sr933/rcc_data/scRCC validation/processed data/scRCC_validation_dataset_lithuanian.pkl"
    if not os.path.exists(pkl_path):
        return name, {"tumour": -1, "non_tumour": -1}
    with open(pkl_path, "rb") as f:
        data = pickle.load(f)
    y = data.get("y_test")
    if y is None:
        return name, {"tumour": -1, "non_tumour": -1}
    tumour = sum(1 for v in y if v == "Tumor cells")
    non_tumour = len(y) - tumour
    return name, {"tumour": tumour, "non_tumour": non_tumour}


def count_tcga() -> Tuple[str, Dict[str, int]]:
    name = "TCGA_KIRC"
    path = "/home/sr933/rcc_data/Bulk RNA/TCGA Bulk RNA/TCGA.KIRC.sampleMap_HiSeqV2.gz"
    if not os.path.exists(path):
        return name, {"tumour": -1, "non_tumour": -1}
    # Read only header to get sample column names
    df = pd.read_csv(path, sep="\t", nrows=0)
    cols = list(df.columns)[1:]  # first column is gene/sample identifier
    # Notebook rule: Tumor if '-01' in sample barcode
    tumour = sum(1 for c in cols if "-01" in c)
    non_tumour = len(cols) - tumour
    return name, {"tumour": tumour, "non_tumour": non_tumour}


def count_czech() -> Tuple[str, Dict[str, int]]:
    name = "Czech_GSE167093"
    folder = "/home/sr933/rcc_data/Bulk RNA/GSE167093"
    if not os.path.isdir(folder):
        return name, {"tumour": -1, "non_tumour": -1}
    all_cols = []
    for fn in os.listdir(folder):
        if not fn.endswith((".txt", ".txt.gz")):
            continue
        path = os.path.join(folder, fn)
        # Mixed compression handling, mirror notebook special-case
        try:
            if fn == "GSE167093_non_norm_matrix13.txt":
                df = pd.read_csv(path, sep='\t', nrows=0)
            else:
                df = pd.read_csv(path, sep='\t', compression='gzip', nrows=0)
        except Exception:
            continue
        # Only consider columns that contain 'AVG' as in the notebook
        cols = [c for c in df.columns if 'AVG' in c]
        all_cols.extend(cols)
    # Labels: Tumor if 'TU' in column name else Benign
    tumour = sum(1 for c in all_cols if 'TU' in c)
    non_tumour = len(all_cols) - tumour
    return name, {"tumour": tumour, "non_tumour": non_tumour}


def main():
    datasets = [
        count_american,
        count_chinese,
        count_metastasis,
        count_lithuanian,
        count_tcga,
        count_czech,
        count_train,
    ]
    rows = []
    for fn in datasets:
        name, counts = fn()
        rows.append({"dataset": name, **counts})

    df = pd.DataFrame(rows)
    # Sort for readability
    df = df[["dataset", "tumour", "non_tumour"]]
    print(df.to_string(index=False))

    # Also save to CSV next to this script
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dataset_counts.csv")
    df.to_csv(out_path, index=False)
    print(f"\nSaved counts to {out_path}")


if __name__ == '__main__':
    main()
