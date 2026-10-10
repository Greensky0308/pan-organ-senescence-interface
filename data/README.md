# Data sources

No data are redistributed with this repository. Download each resource from its
public source and place it under `data/` (or the directory you set through the
environment variables in the top-level README). Accessions below identify the
exact series used.

## Discovery cohorts (one per organ, case vs control)

| Organ | Disease | GEO accession | Platform |
|---|---|---|---|
| Kidney | Diabetic kidney disease | GSE104954 | microarray |
| Liver | Cirrhosis | GSE6764 | microarray (GPL570) |
| Heart | Dilated cardiomyopathy | GSE141910 | RNA-seq |
| Lung | Idiopathic pulmonary fibrosis | GSE150910 | RNA-seq |
| Colon | Inflammatory bowel disease | GSE75214 | microarray |
| Skin | Psoriasis | GSE13355 | microarray |
| Artery | Atherosclerosis | GSE43292 | microarray |
| Esophagus | Esophageal adenocarcinoma | GSE13898 | microarray |
| Brain | Alzheimer disease | GSE5281 | microarray |

## Independent validation cohorts

| Organ | GEO accession | Platform |
|---|---|---|
| Kidney | GSE30529 | microarray |
| Liver | GSE25097 | microarray (GPL10687) |
| Heart | GSE79962 | microarray |
| Lung | GSE24206 | microarray |
| Colon | GSE16879 | microarray |
| Skin | GSE14905 | microarray |
| Artery | GSE28829 | microarray |
| Esophagus | GSE26886 | microarray |
| Brain | GSE28146 | microarray |
| Brain (additional) | GSE48350 | microarray |

## Senescence signatures

Gene lists are read from `raw/` in the analysis workspace:

| Signature | Type | Source |
|---|---|---|
| SenCat | secretory (SASP consensus) | assembled from experimentally supported secretory factors reported across senescence models |
| SenMayo | secretory (SASP) | SenMayo (Saul et al.) |
| Fridman | senescence-up | Fridman & Tainsky |
| Reactome SASP | secretory (SASP) | Reactome |
| Reactome cellular senescence | cell-cycle | Reactome |
| GO:0090398 | cell-cycle | Gene Ontology |

## Age / tissue specificity

- **GTEx v8** (expression, age, eQTL). `07_gtex_age.py`, `09_gtex_tau_all.py`,
  and `23_liver_tissue_eqtl_cirrhosis.py` read GTEx-derived tables.

## Causal analysis inputs

- **FinnGen release 12** summary statistics (per-endpoint `.gz` tabix files).
  Public. Endpoints used are the disease matching each organ.
- **OpenGWAS** plasma pQTL panels (Folkersen et al., Suhre et al., GWAS Catalog)
  and cis-eQTL summary statistics. Queried through the OpenGWAS API with a
  personal access string (`OPENGWAS_JWT`).

## Access notes

All resources above are publicly available. None is controlled-access. No
individual-level genotype or phenotype data are used.
