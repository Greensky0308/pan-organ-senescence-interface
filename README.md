# pan-organ-senescence-interface

Analysis code for a pan-organ **senescence–disease interface**: nine human organ
case-control transcriptome cohorts are intersected with a cross-context
senescence secretory consensus (SenCat) to ask which genes recur across organs,
how the shared interface behaves in normal tissue with age, and whether the
shared core is causally related to disease.

This repository contains the analysis scripts only. It does not ship any data.

## Layout

```
01_assemble.py            assemble per-organ expression + group tables (out/cohorts/)
01b_assemble_more.py      assemble the remaining organs
01c_liver_cohort.py       build the liver cohort (cirrhosis, GSE6764)
02_de.R                   per-organ disease signature (limma/limma-voom)
03_interface.py           disease-senescence overlap + permutation enrichment
04_validate.py            interface module score + pathway enrichment (g:Profiler)
06_multisig.py            repeat with five additional senescence signatures
07_gtex_age.py            GTEx age correlation / aging breadth
09_gtex_tau_all.py        GTEx tissue specificity (tau)
10_subtype_consistency.py cross-subtype consistency per interface gene
24_crossorgan_null.py     cross-organ sharing null (margin-preserving permutation)
11_mr_eqtl.py             cis-eQTL two-sample MR
12_mr_pqtl.py             cis-pQTL two-sample MR
13_coloc_input.py         pull cis summary statistics for colocalisation
14_coloc.R                colocalisation, full cis region
15_screen_mr_coloc.py     lead-SNP Wald screen across interface genes
19_coloc_narrow.R         colocalisation, narrow window around the lead variant
20_validate_lyz.py        independent-cohort colocalisation check
21_liver_cirrhosis.py     liver (cirrhosis) eQTL / pQTL MR
22_liver_cirrhosis_coloc.R liver (cirrhosis) colocalisation
23_liver_tissue_eqtl_cirrhosis.py liver (cirrhosis) tissue-eQTL MR
25_liver_screen_mr_cirrhosis.py   liver (cirrhosis) lead-SNP screen
26_tissue_eqtl_mr.py      tissue-specific cis-eQTL MR across organs
27_validate_build.py      build the remaining validation cohorts
28_validate_kidney.py     kidney validation cohort
29_validation_cohort.py   build the liver (cirrhosis) validation cohort
30_validate_enrich.py     enrichment in the independent validation cohorts
34_figures_5fig.py        main figures 1-5
39_fig4_v2.py             validation figure
37_make_tables.py         main tables
38_make_supp_tables.py    supplementary tables
```

## Environment

```
pip install -r requirements.txt
```

R packages: `limma`, `edgeR`, `coloc`, `data.table`.

Fonts: figures use Arial.

## Configuration

No credentials or machine paths are hard-coded. Set as environment variables
where needed:

- `HTTPS_PROXY` — outbound HTTP proxy (optional).
- `OPENGWAS_JWT` — OpenGWAS access string, for the pQTL queries in
  `12_mr_pqtl.py`, `13_coloc_input.py`, `15_screen_mr_coloc.py`,
  `20_validate_lyz.py`, `21_liver_cirrhosis.py`.
- `KID_DIR` — path to the kidney prior expression data (optional).
- `FINNGEN_DIR` — directory holding the FinnGen R12 summary statistics.

Intermediate downloads are cached under `cache/`.

## Data

No data are redistributed here. Source accessions and provenance are listed in
`data/README.md`. Place downloaded files under `data/` (or the paths you set via
the environment variables above).

## Running

Run the scripts in the order listed in **Layout**: cohort assembly → `02_de.R`
→ `03_interface.py` → the downstream analyses → colocalisation / MR → figures →
tables. All intermediate and output files are written under `out/` (figures in
`out/figures/`, tables in `out/tables/`).

## License

MIT — see `LICENSE`.
