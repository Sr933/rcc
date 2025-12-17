#!/usr/bin/env python3
"""Compute overlaps between aggregated feature lists (LGBM, SHAP, Logistic) and target genes."""

from pathlib import Path
from typing import Dict, Iterable, Tuple

import pandas as pd

AGG_LGBM_PATH = Path("/home/sr933/rcc_data/data/aggregated_lgbm_top100_by_freq.csv")
AGG_SHAP_PATH = Path("/home/sr933/rcc_data/data/aggregated_shap_top100_by_freq.csv")
TARGET_PATH = Path("/home/sr933/rcc/4_network_analysis/data/key_proteins_muanually_filtered.txt")
TARGET_LOGISTIC_PATH = Path("/home/sr933/rcc/data/Target_genes_Logistic.csv")

OUTPUT_DIR = Path("/home/sr933/rcc/data")
SUMMARY_CSV = OUTPUT_DIR / "gene_overlap_summary.csv"
DETAILS_CSV = OUTPUT_DIR / "gene_overlap_details.csv"


def load_gene_set(csv_path: Path) -> Tuple[str, set]:
    """Load first column of CSV as a set of gene names."""
    if not csv_path.exists():
        raise FileNotFoundError(f"Cannot find file: {csv_path}")
    df = pd.read_csv(csv_path)
    first_col = df.columns[0]
    genes = set(df[first_col].astype(str).str.strip())
    return first_col, genes


def summarize_overlap(source: str, source_genes: set, reference: str, reference_genes: set) -> Dict[str, object]:
    overlap = sorted(source_genes & reference_genes)
    source_only = sorted(source_genes - reference_genes)
    reference_only = sorted(reference_genes - source_genes)
    return {
        "source": source,
        "reference": reference,
        "source_count": len(source_genes),
        "reference_count": len(reference_genes),
        "overlap_count": len(overlap),
        "source_only_count": len(source_only),
        "reference_only_count": len(reference_only),
        "overlap_genes": overlap,
        "source_only_genes": source_only,
        "reference_only_genes": reference_only,
    }


def expand(name: str, genes: Iterable[str]) -> pd.DataFrame:
    return pd.DataFrame({"category": name, "gene": list(genes)})


def main() -> None:
    # Load gene sets
    _, agg_lgbm = load_gene_set(AGG_LGBM_PATH)
    _, agg_shap = load_gene_set(AGG_SHAP_PATH)
    _, target = load_gene_set(TARGET_PATH)
    _, target_log = load_gene_set(TARGET_LOGISTIC_PATH)

    comparisons = [
        summarize_overlap("aggregated_lgbm", agg_lgbm, "target", target),
        summarize_overlap("aggregated_shap", agg_shap, "target", target),
        summarize_overlap("target_logistic", target_log, "target", target),
        summarize_overlap("aggregated_lgbm", agg_lgbm, "target_logistic", target_log),
        summarize_overlap("aggregated_shap", agg_shap, "target_logistic", target_log),
    ]

    summary_df = pd.DataFrame(
        [{k: v for k, v in comp.items() if not k.endswith("_genes")} for comp in comparisons]
    )
    summary_df.to_csv(SUMMARY_CSV, index=False)

    detail_frames = []
    for comp in comparisons:
        prefix = f"{comp['source']}__vs__{comp['reference']}"
        detail_frames.append(expand(f"{prefix}__overlap", comp["overlap_genes"]))
        detail_frames.append(expand(f"{prefix}__source_only", comp["source_only_genes"]))
        detail_frames.append(expand(f"{prefix}__reference_only", comp["reference_only_genes"]))

    details_df = pd.concat(detail_frames, ignore_index=True)
    details_df.to_csv(DETAILS_CSV, index=False)

    print("Gene overlap summary (counts):")
    print(summary_df.to_string(index=False))
    print(f"\nSummary saved to {SUMMARY_CSV}")
    print(f"Detailed gene lists saved to {DETAILS_CSV}")


if __name__ == "__main__":
    main()