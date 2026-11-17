"""NCI-H295R experiment: one-way ANOVA with Tukey's test on secreted DHEA-S and PGE2 (Fig. 7b,c;
Supplementary Data 7).

Input : H295R_ELISA_supernatant.csv  (columns: readout, group, well, concentration; DHEA-S in ng/mL, PGE2 in ng/L;
        three wells per group, values as measured in the supernatant)
Output: results/S7_H295R_Tukey.csv
"""
import os
import pandas as pd
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import common as C

d = pd.read_csv(os.path.join(C.DATA, "H295R_ELISA_supernatant.csv"))
rows = []
for readout, g in d.groupby("readout", sort=False):
    f, p = stats.f_oneway(*[x.concentration.values for _, x in g.groupby("group")])
    print(f"{readout}: one-way ANOVA F = {f:.1f}, P = {p:.2e}; group means",
          g.groupby("group").concentration.mean().round(1).to_dict())
    res = pairwise_tukeyhsd(g.concentration, g.group)
    for (g1, g2, diff, padj, *_ ) in res.summary().data[1:]:
        rows.append(dict(readout=readout, group1=g1, group2=g2, meandiff=diff, p_adj=padj))
t = pd.DataFrame(rows)
t.to_csv(C.out_path("S7_H295R_Tukey.csv"), index=False)
print(t.to_string(index=False))
