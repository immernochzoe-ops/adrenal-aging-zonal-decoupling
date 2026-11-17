"""Rhesus macaque adrenal proteomics: all four aged-versus-young comparisons (Supplementary Data 5; Fig. 4;
Supplementary Figure 4d).

Input : MacaqueAdrenal_DIA_protein_groups.tsv (OMIX023586; Spectronaut protein-group table, four animals).
Output: results/S5_macaque_all_pairings.csv

Protein groups quantified in at least two animals are kept (9,026). Each is mapped to the first gene symbol of the
group that is not a locus placeholder; where several groups share a symbol the one with most peptides is kept.
Missing quantities are set to a floor (0.1, log2 = -3.32) and flagged; flagged proteins are excluded from every
count below, which guards against fold changes that only reflect the floor.
"""
import re
import numpy as np
import pandas as pd
import common as C
import os

THRESH = 0.3                      # |log2 fold change| required in each comparison
FLOOR = 0.1
PLACEHOLDER = re.compile(r"^(EGK_|LOC\d|MAMU-)")
ANIMALS = ["Y1", "Y2", "O1", "O2"]
PAIRS = {"O1_Y1": ("O1", "Y1"), "O2_Y2": ("O2", "Y2"), "O1_Y2": ("O1", "Y2"), "O2_Y1": ("O2", "Y1")}
PANEL = ["LDLR", "SCARB1", "SOAT1", "DHCR24", "SULT2A1", "PAPSS2", "POR", "CYB5A", "CYB5B", "CYP17A1", "CYP11A1",
         "CYP11B1", "CYP21A2", "HSD3B2", "FDX1", "FDXR", "STAR", "AKR1C3", "PTGES", "PTGES3"]


def first_symbol(genes):
    for tok in str(genes).split(";"):
        if tok and tok != "nan" and not PLACEHOLDER.match(tok):
            return tok
    return None


pg = pd.read_csv(os.path.join(C.DATA, "MacaqueAdrenal_DIA_protein_groups.tsv"), sep="\t")
pg.columns = ["group", "genes", "description", "npep", "qvalue", "cscore", "mw"] + ANIMALS
pg = pg[pg[ANIMALS].notna().sum(axis=1) >= 2].copy()
print(len(pg), "protein groups quantified in at least two animals")
pg["gene"] = pg["genes"].map(first_symbol)
pg["Accession"] = pg["group"].str.split(";").str[0]
t = (pg[pg.gene.notna()].sort_values("npep", kind="stable").drop_duplicates("gene", keep="last")
     .sort_values("npep", ascending=False, kind="stable").reset_index(drop=True))
print(len(t), "genes")

log2 = np.log2(t[ANIMALS].fillna(0).clip(lower=FLOOR))
t["anyfloor"] = t[ANIMALS].isna().any(axis=1)
for name, (o, y) in PAIRS.items():
    t[name] = log2[o] - log2[y]
for a in ANIMALS:
    t["log2_ADR_" + a] = log2[a]
fc = t[list(PAIRS)]
t["mean4"] = fc.mean(axis=1)
t["all4_down"], t["all4_up"] = (fc < -THRESH).all(axis=1), (fc > THRESH).all(axis=1)
t["orig_down"], t["orig_up"] = (fc[["O1_Y1", "O2_Y2"]] < -THRESH).all(axis=1), (fc[["O1_Y1", "O2_Y2"]] > THRESH).all(axis=1)
t["alt_down"], t["alt_up"] = (fc[["O1_Y2", "O2_Y1"]] < -THRESH).all(axis=1), (fc[["O1_Y2", "O2_Y1"]] > THRESH).all(axis=1)
cols = (["gene", "Accession", "group", "npep", "anyfloor"] + list(PAIRS) + ["log2_ADR_" + a for a in ANIMALS] +
        ["mean4", "all4_down", "all4_up", "orig_down", "alt_down", "orig_up", "alt_up"])
t[cols].to_csv(C.out_path("S5_macaque_all_pairings.csv"), index=False)

ok = ~t.anyfloor
print(ok.sum(), "proteins without floor values")
print("concordant in all four aged-versus-young comparisons:", int(((t.all4_down | t.all4_up) & ok).sum()))
# calibration: the two partitions of the four animals that ignore age
for a, b in [(("Y1", "O1"), ("Y2", "O2")), (("Y1", "O2"), ("Y2", "O1"))]:
    f = pd.DataFrame({f"{x}_{y}": log2[x] - log2[y] for x in a for y in b})
    print(f"  age-blind partition {a} vs {b}:", int((((f < -THRESH).all(axis=1) | (f > THRESH).all(axis=1)) & ok).sum()))
panel = t[t.gene.isin(PANEL)].set_index("gene")
print("\nsteroidogenic panel, lower in all four comparisons:", sorted(panel.index[panel.all4_down & ~panel.anyfloor]))
print(panel.loc[[g for g in PANEL if g in panel.index], list(PAIRS) + ["mean4", "anyfloor", "all4_down"]].round(2).to_string())
