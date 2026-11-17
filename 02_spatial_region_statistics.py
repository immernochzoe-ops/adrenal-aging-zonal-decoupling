"""Region-level aged-versus-young statistics for every protein quantified in both specimens, by zone
(Supplementary Data 2 and 3; Fig. 2a,c; Supplementary Figure 3c,e).

For each protein and zone: log2(aged/young) as the difference of mean log2(x+1) across strips, difference of medians,
two-sided Mann-Whitney U test on per-strip values, Benjamini-Hochberg q across all proteins of the zone, the
percentile rank of the effect within the zone, and the age coefficient after entering the medullary admixture
index as a covariate (reported, but not used to support any claim).

Output: results/S3_genomewide.csv, results/S2_region_stats.csv
"""
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import common as C

PANEL = ["SCARB1", "CHGA", "CYB5A", "CYP17A1", "DBH", "DHCR24", "CYB5R3", "FDX1", "FDXR", "HSD3B2", "LDLR", "PNMT",
         "POR", "SULT2A1", "TH", "VSNL1", "AKR1C3", "PAPSS2", "PTGES", "PTGES3", "CYB5B", "DACH1", "SOAT1"]

L, zone, mi = {}, {}, {}
for specimen in ("young", "old"):
    L[specimen], _ = C.log2_tissue_strips(C.load_strip_matrix(specimen))
    zone[specimen] = C.assign_zones(C.panel_scores(L[specimen]))
    mi[specimen] = C.medullary_index(L[specimen])
genes = [g for g in L["young"].index if g in L["old"].index]          # quantified in both specimens


def region_stats(y_strips, o_strips):
    a = L["young"].loc[genes, y_strips].values
    b = L["old"].loc[genes, o_strips].values
    eff = b.mean(axis=1) - a.mean(axis=1)
    p = stats.mannwhitneyu(b, a, axis=1, alternative="two-sided").pvalue
    # age coefficient with the medullary admixture index as covariate (ordinary least squares, all genes at once)
    X = np.c_[np.ones(len(y_strips) + len(o_strips)), np.r_[np.zeros(len(y_strips)), np.ones(len(o_strips))],
              np.r_[mi["young"][y_strips].values, mi["old"][o_strips].values]]
    Yv = np.c_[a, b].T
    beta, *_ = np.linalg.lstsq(X, Yv, rcond=None)
    resid = Yv - X @ beta
    dof = X.shape[0] - X.shape[1]
    se = np.sqrt((resid ** 2).sum(axis=0) / dof * np.linalg.inv(X.T @ X)[1, 1])
    with np.errstate(invalid="ignore", divide="ignore"):
        p_adj = 2 * stats.t.sf(np.abs(beta[1] / se), dof)
    p_adj = np.where(np.isfinite(p_adj), p_adj, 1.0)     # protein not quantified in any strip of the zone: no test, P = 1
    n = len(genes)
    return pd.DataFrame(dict(gene=genes, nY=len(y_strips), nO=len(o_strips), log2_OvY=eff,
                             dmedian=np.median(b, axis=1) - np.median(a, axis=1), p_MWU=p,
                             beta_age_medadj=beta[1], p_medadj=p_adj,
                             q_MWU=multipletests(p, method="fdr_bh")[1], q_medadj=multipletests(p_adj, method="fdr_bh")[1],
                             pct_raw=stats.rankdata(eff) / n * 100, pct_medadj=stats.rankdata(beta[1]) / n * 100))


if __name__ == "__main__":
    tabs = []
    for z in C.ZONES:
        t = region_stats(zone["young"][zone["young"] == z].index, zone["old"][zone["old"] == z].index)
        t.insert(1, "zone", z)
        tabs.append(t)
    full = pd.concat(tabs, ignore_index=True)
    full.to_csv(C.out_path("S3_genomewide.csv"), index=False)
    full[full.gene.isin(PANEL)].to_csv(C.out_path("S2_region_stats.csv"), index=False)
    zr = full[full.zone == "ZR"].set_index("gene")
    print(f"{len(genes)} proteins quantified in both specimens")
    print(zr.loc[["LDLR", "SCARB1", "DHCR24", "SOAT1", "FDXR", "CYB5B", "POR", "SULT2A1", "HSD3B2", "AKR1C3", "CYB5A",
                  "CYP17A1", "PAPSS2"], ["log2_OvY", "q_MWU", "pct_raw"]].round(3).to_string())
