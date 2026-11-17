# data

Small inputs provided with the code:

- `DECONV_MuSiC_coarse_cortex.csv`, `DECONV_Bisque_coarse_cortex.csv` — deconvolved cell-type fractions of the 295 GTEx v10
  adrenal samples (MuSiC and BisqueRNA; output of `06_singlecell_meta_and_gtex_deconvolution.Rmd`).
- `GTEx_ZR_marker_proxy.csv` — expression-derived zone marker proxies per GTEx sample; column `ZR` is used by
  `05_gtex_age_associations.py`.
- `H295R_ELISA_supernatant.csv` — secreted DHEA-S (ng/mL) and PGE2 (ng/L), three wells per group, as measured.

Add the larger inputs listed in the top-level README (macaque protein-group table from OMIX023586, GTEx v10 files,
human strip matrices) to this folder before running the scripts.
