"""Compares the tables written by the scripts with the published Supplementary Data files.

Usage: python validate_against_supplementary_data.py <folder holding 'Supplementary Data 1.xlsx' ... 'Supplementary Data 7.xlsx'>
"""
import os
import sys
import numpy as np
import pandas as pd
import common as C

SHEETS = {   # result file: (Supplementary Data number, sheet, key columns)
    "S1_strip_QC_zones.csv": (1, "S1_strip_QC_zones", ["strip"]),
    "S1b_medullary_index.csv": (1, "S1b_medullary_index", ["specimen", "zone"]),
    "S2_region_stats.csv": (2, "S2_region_stats", ["gene", "zone"]),
    "S3_genomewide.csv": (3, "S3_genomewide", ["gene", "zone"]),
    "S4a_zone_sensitivity.csv": (4, "S4a_zone_sensitivity", ["variant", "gene", "zone"]),
    "S5_macaque_all_pairings.csv": (5, "S5_macaque_all_pairings", ["gene"]),
    "S6a_GTEx_adjustment.csv": (6, "S6a_GTEx_adjustment", ["gene"]),
    "S6b_GTEx_sex.csv": (6, "S6b_GTEx_sex", ["gene"]),
    "S6c_GTEx_strata.csv": (6, "S6c_GTEx_strata", ["gene"]),
    "S7_H295R_Tukey.csv": (7, "S7_H295R_Tukey", ["readout", "group1", "group2"]),
}
folder = sys.argv[1]
for f, (n, sheet, keys) in SHEETS.items():
    path = os.path.join(C.OUT, f)
    if not os.path.exists(path):
        print(f"{f:32s} not generated"); continue
    mine = pd.read_csv(path).set_index(keys)
    ref = pd.read_excel(os.path.join(folder, f"Supplementary Data {n}.xlsx"), sheet_name=sheet).set_index(keys)
    rows = mine.index.intersection(ref.index)
    cols = [c for c in ref.columns if c in mine.columns]
    a, b = mine.loc[rows, cols], ref.loc[rows, cols]
    worst, mismatched = 0.0, 0
    for c in cols:
        if pd.api.types.is_numeric_dtype(b[c]) and b[c].dtype != bool:
            d = np.abs(a[c].astype(float) - b[c].astype(float))
            worst = max(worst, float(np.nanmax(d)) if d.notna().any() else 0.0)
        else:
            mismatched += int((a[c].astype(str) != b[c].astype(str)).sum())
    print(f"{f:32s} rows {len(rows)}/{len(ref)}  columns {len(cols)}/{ref.shape[1]}  "
          f"largest numeric difference {worst:.2e}  non-numeric mismatches {mismatched}")
