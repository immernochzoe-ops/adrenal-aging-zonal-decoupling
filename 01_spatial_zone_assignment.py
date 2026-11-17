"""Strip-level quality control, marker-panel scores, zone assignment and medullary admixture
(Supplementary Data 1; Supplementary Figures 1a-c and 3a,b).

Input : per-strip protein matrices of the young and aged human adrenal (see README, 'Data').
Output: results/S1_strip_QC_zones.csv, results/S1b_medullary_index.csv
"""
import numpy as np
import pandas as pd
import common as C

rows, admix = [], []
for specimen in ("young", "old"):
    raw = C.load_strip_matrix(specimen)
    log2, n_quant = C.log2_tissue_strips(raw)
    scores = C.panel_scores(log2)
    zone = C.assign_zones(scores)
    for strip in raw.columns:
        kept = strip in log2.columns
        rows.append(dict(specimen=specimen, strip=strip, axis=strip.split("_")[1][0], n_quantified=int(n_quant[strip]),
                         total_intensity=raw[strip].sum(), kept=kept,
                         zone=zone[strip] if kept else "excluded (no tissue)",
                         **({k: scores.loc[strip, k] for k in ("ZG", "ZF", "ZR", "Medulla")} if kept else {})))
    mi = C.medullary_index(log2)
    print(f"{specimen}: {log2.shape[1]} tissue-bearing strips of {raw.shape[1]}; "
          f"medullary index vs CHGA r = {np.corrcoef(mi, log2.loc['CHGA'])[0, 1]:.2f}")
    for z in ("ZG", "ZF", "ZR", "Medulla"):
        s = zone[zone == z].index
        admix.append(dict(specimen=specimen, zone=z, n_strips=len(s), medullary_index=round(mi[s].mean(), 5),
                          CHGA=round(log2.loc["CHGA", s].mean(), 5)))
pd.DataFrame(rows).to_csv(C.out_path("S1_strip_QC_zones.csv"), index=False)
pd.DataFrame(admix).to_csv(C.out_path("S1b_medullary_index.csv"), index=False)
print(pd.DataFrame(admix).to_string(index=False))
