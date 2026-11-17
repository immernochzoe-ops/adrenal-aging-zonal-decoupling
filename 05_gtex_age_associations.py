"""GTEx v10 adrenal gland: association of expression with donor age under covariate, composition and
sex-stratified analysis (Table 1; Fig. 6c; Supplementary Figure 6d,e; Supplementary Data 6).

Inputs (see README):
  gene_tpm_v10_adrenal_gland.gct                              GTEx v10 adrenal TPM matrix
  GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt       sex, age bracket, Hardy scale
  GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt        RNA integrity number, ischaemic time
  DECONV_MuSiC_coarse_cortex.csv, DECONV_Bisque_coarse_cortex.csv   deconvolved fractions (06_gtex_deconvolution.Rmd)
  GTEx_ZR_marker_proxy.csv                                    data-derived ZR marker proxy per sample
Outputs: results/S6a_GTEx_adjustment.csv, S6b_GTEx_sex.csv, S6c_GTEx_strata.csv

Associations are Spearman correlations with ordinal age bracket. Adjusted values are partial Spearman
correlations: expression and age are rank-transformed, both are regressed on the covariates, and the
residuals are correlated (P from the t distribution with n - 2 - k degrees of freedom).
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import common as C

GENES_A = ["CYB5A", "SULT2A1", "AKR1C3", "OSBPL6", "PAPSS2", "LDLR", "SCARB1", "SOAT1", "HMGCR", "INSIG1", "CYP11B1",
           "CYP11B2", "CYP17A1", "STAR", "CYP11A1", "NR5A1", "FDX1", "FDXR", "HSD3B2", "DACH1", "VSNL1"]
GENES_B = ["CYB5A", "SULT2A1", "AKR1C3", "PAPSS2", "LDLR", "SCARB1", "SOAT1", "HMGCR", "INSIG1", "CYP11B1", "CYP11B2",
           "DACH1", "VSNL1"]                                   # Table 1
GENES_C = ["CYB5A", "SULT2A1", "AKR1C3", "OSBPL6", "CYP17A1", "PAPSS2", "LDLR", "SCARB1", "SOAT1", "HMGCR", "INSIG1",
           "DHCR7", "CYP11B1", "CYP11B2", "CYP11A1", "STAR", "NR5A1", "FDX1", "FDXR", "HSD3B2", "MC2R", "CYB5B", "POR",
           "PTGES", "PTGES3", "CHGA", "TH", "IGFBP2", "PSME4", "DHCR24", "VSNL1", "DACH1", "CDKN1A"]
p = lambda f: os.path.join(C.DATA, f)

tpm = pd.read_csv(p("gene_tpm_v10_adrenal_gland.gct"), sep="\t", skiprows=2)
expr = tpm.drop_duplicates("Description").set_index("Description").drop(columns=["Name"]).astype(float)
meta = pd.DataFrame(index=expr.columns)
meta["SUBJID"] = ["-".join(s.split("-")[:2]) for s in meta.index]
meta = meta.join(pd.read_csv(p("GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt"), sep="\t").set_index("SUBJID"), on="SUBJID")
sa = pd.read_csv(p("GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"), sep="\t", usecols=["SAMPID", "SMRIN", "SMTSISCH"])
meta = meta.join(sa.set_index("SAMPID"))
meta["age"] = meta.AGE.map({a: i + 1.0 for i, a in enumerate(sorted(meta.AGE.unique()))})   # ordinal age bracket
meta["male"] = (meta.SEX == 1).astype(float)
meta["vent"] = (meta.DTHHRDY == 0).astype(float)                                            # Hardy 0 = ventilator case
hardy = pd.get_dummies(meta.DTHHRDY.astype(int), prefix="hardy", drop_first=True).astype(float)
meta["music"] = pd.read_csv(p("DECONV_MuSiC_coarse_cortex.csv")).set_index("SampleID")["Cortex"]
meta["bisque"] = pd.read_csv(p("DECONV_Bisque_coarse_cortex.csv")).set_index("SampleID")["Cortex"]
meta["zr_proxy"] = pd.read_csv(p("GTEx_ZR_marker_proxy.csv"), index_col=0)["ZR"]
print(len(meta), "donors;", int(meta.male.sum()), "men,", int((1 - meta.male).sum()), "women;",
      int(meta.vent.sum()), "ventilator cases")
print("ventilator vs other, mean age-bracket midpoint:",
      meta.groupby("vent").AGE.apply(lambda s: s.str[:2].astype(int).add(4.5).mean()).round(1).to_dict())


def partial_spearman(x, y, covariates=None):
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    k = 0
    if covariates is not None and covariates.shape[1]:
        Z = np.c_[np.ones(len(rx)), covariates]
        rx = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]
        ry = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
        k = np.linalg.matrix_rank(Z) - 1
    r = np.corrcoef(rx, ry)[0, 1]
    n = len(rx)
    return r, 2 * stats.t.sf(abs(r * np.sqrt((n - 2 - k) / (1 - r ** 2))), n - 2 - k)


def assoc(gene, cov_cols=(), with_hardy=False, subset=None):
    idx = meta.index if subset is None else meta.index[subset]
    cov = meta.loc[idx, list(cov_cols)].values
    if with_hardy:
        h = hardy.loc[idx]
        cov = np.c_[cov, h.loc[:, h.var() > 0].values]
    return partial_spearman(expr.loc[gene, idx].values, meta.loc[idx, "age"].values, cov)


DONOR = ("male", "SMRIN", "SMTSISCH")          # donor covariates (Hardy scale is added as indicator variables)


def table(genes, columns):
    rows = []
    for g in genes:
        row = {"gene": g}
        for name, kw in columns.items():
            row["rho_" + name], row["p_" + name] = assoc(g, **kw)
        rows.append(row)
    return pd.DataFrame(rows)


# S6a: progressive adjustment for deconvolved composition and donor covariates
a = table(GENES_A, {"raw": {}, "MuSiC": dict(cov_cols=("music",)), "Bisque": dict(cov_cols=("bisque",)),
                    "both": dict(cov_cols=("music", "bisque")), "covar": dict(cov_cols=DONOR, with_hardy=True),
                    "covar_both": dict(cov_cols=DONOR + ("music", "bisque"), with_hardy=True),
                    "covar_ZRproxy": dict(cov_cols=DONOR + ("zr_proxy",), with_hardy=True)})
a.to_csv(C.out_path("S6a_GTEx_adjustment.csv"), index=False)

# S6b: sex-stratified (Table 1)
men, women = (meta.male == 1).values, (meta.male == 0).values
rows = []
for g in GENES_B:
    r = {"gene": g}
    r["rho_both"], r["p_both"] = assoc(g)
    r["rhoadj_both"], r["padj_both"] = assoc(g, DONOR, True)
    for name, s in (("male", men), ("female", women)):
        r["rho_" + name], r["p_" + name] = assoc(g, subset=s)
        r["rhoadj_" + name], r["padj_" + name] = assoc(g, ("SMRIN", "SMTSISCH"), True, subset=s)
    r["median_tpm"] = expr.loc[g].median()
    rows.append(r)
b = pd.DataFrame(rows)
b.insert(3, "q_both", multipletests(b.p_both, method="fdr_bh")[1])     # BH across the genes of Table 1
b.to_csv(C.out_path("S6b_GTEx_sex.csv"), index=False)

# S6c: each covariate set, and within the ventilator and non-ventilator strata
vent, nonvent = (meta.vent == 1).values, (meta.vent == 0).values
c = table(GENES_C, {"raw": {}, "sex": dict(cov_cols=("male",)), "sex_vent": dict(cov_cols=("male", "vent")),
                    "sex_hardy": dict(cov_cols=("male",), with_hardy=True), "full": dict(cov_cols=DONOR, with_hardy=True),
                    "in_vent": dict(subset=vent), "in_nonvent": dict(subset=nonvent)})   # unadjusted within stratum
c.insert(c.columns.get_loc("p_in_vent") + 1, "n_vent", int(vent.sum()))
c["n_nonvent"] = int(nonvent.sum())
c["median_tpm"] = expr.loc[GENES_C].median(axis=1).values
c.to_csv(C.out_path("S6c_GTEx_strata.csv"), index=False)
print(b[["gene", "rho_both", "q_both", "rhoadj_both", "padj_both", "rho_male", "p_male", "rho_female", "p_female"]].round(4).to_string(index=False))
