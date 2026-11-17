"""Robustness of the spatial effect sizes to the zone-assignment rule (Supplementary Data 4, sheet S4a;
Supplementary Figures 1d,e and 3d).

Seven variants are run on the same per-strip log2 values as the main analysis: the rule as used, a strict argmax
rule without the medulla threshold, a stricter coverage filter (strips quantifying >= 4,000 protein groups), and
four alternative marker panels.

Output: results/S4a_zone_sensitivity.csv
"""
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import common as C

GENES = ["LDLR", "SCARB1", "SOAT1", "DHCR24", "SULT2A1", "CYB5A", "AKR1C3", "CYP17A1", "PAPSS2", "HSD3B2", "FDX1",
         "FDXR", "POR", "CYB5B", "CYB5R3"]


def panels(**changes):
    p = {k: list(v) for k, v in C.PANELS.items()}
    p.update(changes)
    return p


VARIANTS = {   # name: (marker panels, assignment rule, minimum protein groups per strip)
    "zone rule as used": (C.PANELS, "as_used", 1),
    "strict argmax rule": (C.PANELS, "argmax", 1),
    "coverage filter >= 4,000": (C.PANELS, "as_used", 4000),
    "ZF without FDX1": (panels(ZF=["SCARB1", "HSD3B2"]), "as_used", 1),
    "ZF = HSD3B2 only": (panels(ZF=["HSD3B2"]), "as_used", 1),
    "ZG + CALB2": (panels(ZG=["DACH1", "VSNL1", "CALB2"]), "as_used", 1),
    "ZR expanded": (panels(ZR=["SULT2A1", "CYB5A", "AKR1C3", "PAPSS2"]), "as_used", 1),
}

raw = {s: C.load_strip_matrix(s) for s in ("young", "old")}
base = {}
rows = []
for name, (pan, rule, min_q) in VARIANTS.items():
    L, zone = {}, {}
    for s in raw:
        L[s], _ = C.log2_tissue_strips(raw[s], min_quantified=min_q)
        zone[s] = C.assign_zones(C.panel_scores(L[s], pan), rule)
    if name == "zone rule as used":
        base = zone
    agree = {s: (zone[s] == base[s][zone[s].index]).mean() if pan is not C.PANELS else np.nan for s in raw}
    for z in ("ZR", "ZF", "ZG"):
        y, o = zone["young"][zone["young"] == z].index, zone["old"][zone["old"] == z].index
        a, b = L["young"].loc[GENES, y].values, L["old"].loc[GENES, o].values
        p = stats.mannwhitneyu(b, a, axis=1, alternative="two-sided").pvalue
        q = multipletests(p, method="fdr_bh")[1]
        for i, g in enumerate(GENES):
            rows.append(dict(variant=name, gene=g, zone=z, nY=len(y), nO=len(o),
                             log2_OvY=round(b[i].mean() - a[i].mean(), 4), p=round(p[i], 4), q=round(q[i], 4),
                             agree_young=round(agree["young"], 3), agree_old=round(agree["old"], 3)))
t = pd.DataFrame(rows)
t.to_csv(C.out_path("S4a_zone_sensitivity.csv"), index=False)
chol = t[t.gene.isin(["LDLR", "SCARB1", "SOAT1", "DHCR24"]) & t.zone.isin(["ZR", "ZF"])]
print("cholesterol-acquisition proteins negative in ZR and ZF under every variant:", bool((chol.log2_OvY < 0).all()))
print(t[t.zone == "ZR"].pivot(index="gene", columns="variant", values="log2_OvY").loc[GENES, list(VARIANTS)].round(2).to_string())
