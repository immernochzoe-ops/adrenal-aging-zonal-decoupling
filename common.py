"""Shared helpers for the spatial-proteomics analysis."""
import os
import numpy as np
import pandas as pd

DATA = os.environ.get("ADRENAL_DATA", "data")       # folder holding the input files (see README)
OUT = os.environ.get("ADRENAL_OUT", "results")      # folder the scripts write to

PANELS = {                                           # zone-marker panels (Methods: Cortical zone assignment)
    "Medulla": ["CHGA", "TH", "PNMT", "DBH"],
    "ZR": ["SULT2A1", "CYB5A"],
    "ZF": ["FDX1", "SCARB1", "HSD3B2"],               # "non-ZR cortex" panel
    "ZG": ["DACH1", "VSNL1"],
}
MEDULLA_Z = 0.75                                     # a strip is medulla if its medullary panel z exceeds this
MEDULLA_33 = ("ABAT AP3B2 CAMK2B CHD5 CRYZ DBH DDC DPP6 DPYSL5 ELAVL4 EML5 F5 GNAO1 HMGCLL1 KIF1A KIF5C LSAMP "
              "MAP6 MAP7 MAPK10 NRCAM PCSK1 PCSK2 PHYHIPL PPFIA2 SEPTIN5 SEZ6L2 SHROOM3 SNAP91 SYP TH VAT1L VWA5B2").split()
ZONES = ["ZR", "ZF", "ZG", "Medulla"]
SPECIMENS = {"young": "HumanAdrenal_Young16y", "old": "HumanAdrenal_Aged67y"}


def out_path(name):
    os.makedirs(OUT, exist_ok=True)
    return os.path.join(OUT, name)


def load_strip_matrix(specimen):
    """Per-strip protein matrix of one specimen (rows = gene symbols, columns = strips, linear scale).

    Accepts either one file <prefix>_strip_protein_matrix.csv or the pair
    <prefix>_horizontal_strips_protein_matrix.csv / <prefix>_vertical_strips_protein_matrix.csv."""
    prefix = os.path.join(DATA, SPECIMENS[specimen])
    single = prefix + "_strip_protein_matrix.csv"
    if os.path.exists(single):
        d = pd.read_csv(single)
    else:
        h = pd.read_csv(prefix + "_horizontal_strips_protein_matrix.csv")
        v = pd.read_csv(prefix + "_vertical_strips_protein_matrix.csv")
        assert (h["Gene_symbol"] == v["Gene_symbol"]).all()
        d = pd.concat([h, v.drop(columns=["Gene_ID", "Gene_symbol", "Description"])], axis=1)
    d = d.set_index("Gene_symbol").drop(columns=["Gene_ID", "Description"])
    return d.astype(float)


def log2_tissue_strips(raw, min_quantified=1):
    """log2(x + 1) of the strips that carry tissue (at least `min_quantified` protein groups quantified)."""
    n_quant = (raw > 0).sum()
    kept = n_quant[n_quant >= min_quantified].index
    return np.log2(raw[kept] + 1), n_quant


def panel_scores(log2, panels=PANELS):
    """Mean z-score of each marker panel per strip; z-scores are taken across strips within the specimen."""
    z = log2.sub(log2.mean(axis=1), axis=0).div(log2.std(axis=1, ddof=1), axis=0)
    return pd.DataFrame({k: z.loc[[g for g in v if g in z.index]].mean() for k, v in panels.items()})


def assign_zones(scores, rule="as_used"):
    """Zone of each strip. 'as_used': medulla if the medullary panel exceeds MEDULLA_Z, otherwise the cortical
    panel with the highest score. 'argmax': the panel with the highest score among all four."""
    if rule == "argmax":
        return scores[["ZR", "ZF", "ZG", "Medulla"]].idxmax(axis=1)
    cortical = scores[["ZR", "ZF", "ZG"]].idxmax(axis=1)
    return cortical.where(scores["Medulla"] <= MEDULLA_Z, "Medulla")


def medullary_index(log2):
    """Mean log2 abundance of the 33 medulla-restricted proteins, per strip."""
    return log2.loc[[g for g in MEDULLA_33 if g in log2.index]].mean()
