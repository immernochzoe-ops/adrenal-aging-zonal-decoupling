# adrenal-aging-zonal-decoupling

Analysis code for the article

> Ding G., Xu Y., Guo T., Dai Y. & Feng C. *Spatially resolved proteomics reveals zonally structured
> transcriptome–proteome decoupling in the aging human adrenal cortex.*

Each script regenerates one block of the quantitative results and writes the corresponding table of the
Supplementary Data to `results/`. `validate_against_supplementary_data.py` compares those tables with the published
files.

## Scripts

| Script | Analysis | Reproduces |
|---|---|---|
| `01_spatial_zone_assignment.py` | strip-level quality control, marker-panel z-scores, zone assignment, medullary admixture index | Supplementary Data 1 |
| `02_spatial_region_statistics.py` | aged-versus-young statistics per zone for all 10,350 proteins: effect size, Mann–Whitney test, BH q, percentile rank, admixture-adjusted coefficient | Supplementary Data 2 and 3 |
| `03_spatial_zone_sensitivity.py` | seven zone-assignment variants (alternative marker panels, strict argmax rule, coverage filter) | Supplementary Data 4 (sheet S4a) |
| `04_macaque_proteomics_pairings.py` | rhesus macaque proteomics under all four aged-versus-young pairings, with the age-blind calibration | Supplementary Data 5 |
| `05_gtex_age_associations.py` | GTEx v10 adrenal: Spearman and partial Spearman age associations with donor covariates, deconvolved composition, sex stratification and ventilator strata | Table 1, Supplementary Data 6 |
| `06_singlecell_meta_and_gtex_deconvolution.Rmd` | human adrenal single-cell meta-analysis (consensus zonal signature, per-donor module scores) and GTEx deconvolution with MuSiC and BisqueRNA (R notebook, as run) | inputs of `05`, single-cell results |
| `07_h295r_elisa_statistics.py` | one-way ANOVA and Tukey's test for the NCI-H295R ELISA readouts | Supplementary Data 7 |
| `common.py` | marker panels, loaders and shared functions | |

Not in this repository: the Flow2spatial per-pixel reconstruction (the published Flow2spatial method, ref. 9 of the
article), the haematoxylin–eosin morphometry and held-out signature analysis (Supplementary Data 9), processing of the
laser-microdissection arrays (Supplementary Data 4, sheet S4b), Western blot densitometry and figure rendering.

## Agreement with the published tables

Run on the study data, the scripts reproduce Supplementary Data 1–4 (S4a) and 6 to floating-point precision
(largest difference < 1e-11) and Supplementary Data 7 to the rounding of the input concentrations. For Supplementary
Data 5, 7,344 of 7,347 genes are identical in every value; three genes (MXRA7, RWDD3, STK11) are each represented by
two protein groups with equal peptide counts and the published table kept the other group (none of the three is
reported in the article), and for 53 groups the published `Accession` label names a different member of the same
protein group. The `d_vent` column of sheet S6c is not regenerated.

## Data

Put the input files in `data/` (or point `ADRENAL_DATA` to another folder).

| File | Source |
|---|---|
| `MacaqueAdrenal_DIA_protein_groups.tsv` | OMIX023586 (https://ngdc.cncb.ac.cn/omix/release/OMIX023586) |
| `HumanAdrenal_Young16y_…` and `HumanAdrenal_Aged67y_…strips_protein_matrix.csv` | human spatial proteomics; available from the corresponding authors on reasonable request (see the article's Data availability statement) |
| `gene_tpm_v10_adrenal_gland.gct`, `GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt`, `GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt` | GTEx portal, v10 open-access files (https://gtexportal.org) |
| `DECONV_MuSiC_coarse_cortex.csv`, `DECONV_Bisque_coarse_cortex.csv`, `GTEx_ZR_marker_proxy.csv` | provided here in `data/`: per-sample deconvolved fractions (output of script 06) and the expression-derived ZR marker proxy |
| `H295R_ELISA_supernatant.csv` | provided here in `data/`: measured supernatant concentrations, three wells per group |

Rhesus macaque single-nucleus RNA-seq reads: GSA CRA047097 (BioProject PRJCA069995).

## Running

```
pip install -r requirements.txt
python 01_spatial_zone_assignment.py
python 02_spatial_region_statistics.py
python 03_spatial_zone_sensitivity.py
python 04_macaque_proteomics_pairings.py
python 05_gtex_age_associations.py
python 07_h295r_elisa_statistics.py
python validate_against_supplementary_data.py "<folder with Supplementary Data 1.xlsx ... 7.xlsx>"
```

Python 3.10, numpy 2.2, pandas 2.3, scipy 1.15 and statsmodels 0.14 were used. Script 06 needs R with Seurat, MuSiC
and BisqueRNA.

## Licence

MIT (see `LICENSE`).
