# P4HA1 in breast cancer: spatial transcriptomics analysis

Analysis code for a study of **P4HA1** (collagen prolyl 4-hydroxylase alpha-1) in
breast cancer (BRCA), using Visium spatial transcriptomics from the HTAN WUSTL
atlas.

This is the breast cancer port of the pancreatic cancer pipeline in
[p4ha1-pdac-spatial](https://github.com/Mahditajbakhsh/p4ha1-pdac-spatial). The
question, the gene modules, the models and the section structure are the same, so
the two cancers can be compared directly. Only the cohort and the pooling unit
change.

Izzi Lab, ECM-Hypoxia unit, Faculty of Biochemistry and Molecular Medicine,
University of Oulu.

> **Status: code only, no results yet.** The pipeline has not been run on the BRCA
> data. It was checked by parsing it in R and by running the new pooling, file
> matching and sections 21 and 22 on synthetic numbers. Unpublished work in
> progress; do not redistribute.

---

## The cohort

Source: Mo et al., *Nature* 2024, doi:[10.1038/s41586-024-08087-4](https://doi.org/10.1038/s41586-024-08087-4).
The list of sections comes from `Sample_ID_Lookup_table_v1.xlsx`, sheet Visium, in
[ding-lab/ST_subclone_publication](https://github.com/ding-lab/ST_subclone_publication)
(folder `Data_access`), and is stored here as `docs/brca_visium_samples.csv`.

| | Sections | Patients |
|---|---|---|
| All BRCA Visium | 54 | 22 |
| Primary tumours (main analysis) | 49 | 21 |
| Metastasis (sensitivity only) | 5 | 1 (HT268B1) |

As a check, the same table gives 23 PDAC sections from 16 patients, which matches
the PDAC pipeline.

All 54 sections have a Level 4 Seurat object on the HTAN portal: checked against
the atlas Level 4 manifest `syn61518654`, where every BRCA file name matches its
`sample_id` exactly (38 OCT, 16 FFPE).

## What differs from the PDAC pipeline

1. **The patient is the pooling unit.** Several BRCA patients contribute 4 to 6
   serial sections cut 5 to 100 um apart. Pooling over sections, as in PDAC, would
   let one tumour count up to six times. `meta_pool()` first combines the sections
   of each patient (fixed effect, inverse variance, the same collapse PDAC section
   21 used) and then pools patients with DerSimonian-Laird. It stops with an error
   if it is called without section IDs, so no step can fall back to section level
   silently. Section 21 reports the section level result as the sensitivity.
2. **Primary tumours only** in the main run. `INCLUDE_METASTASIS <- TRUE` adds the
   5 metastatic sections and writes to a separate output folder; section 22 then
   reports the result with and without HT268B1.
3. **Leave one patient out** (section 22) instead of leave one section out.
4. **HT339B1 and HT339B2** have different HTAN participant IDs but share a prefix.
   Section 21 repeats the pooling with them merged, as PDAC did for HT270P1 and
   HT270P2.
5. **One project root.** All paths come from the `P4HA1_BRCA_ROOT` environment
   variable.
6. **PDAC numbers removed from the logic.** Comments that quote PDAC results are
   marked as PDAC, and each section header from 7 onward carries a note saying so.
   Every code change is tagged `BRCA:` in the script.

## Data issues found in the BRCA objects

* **AddModuleScore cannot form 20 bins in three sections** (HT262B1-S1H2Fc2U1Z1Bs1,
  HT486B1-S1H1Fs1U1Bp1, HT397B1-S1H3A1U41Bp1). In these sparse sections 700 to
  1,800 genes have exactly one UMI in one spot, so they share one identical mean
  expression and Seurat's binning fails. The pipeline now tries 20 bins first
  (unchanged for every other section, same seed, same score) and only on that
  error retries with 15 and then 10 bins. Each fallback is logged and written to
  `nbin_fallback_log.csv` in the output folder. In the first run, before this
  fix, the PDAC code dropped these sections silently from sections 3 to 11.
* **HT397B1-S1H3A1U61Bp1 has no P4HA1 in the object** (7,498 features only) and is
  excluded by the existing check in section 3.
* Some sections lack individual module genes (CA9, HAS1, CXCL1, CXCL2). Modules are
  scored from the genes present, as in PDAC.

## Decisions kept from PDAC, and what they mean here

* **Same gene sets.** `iCAF` and `myCAF` modules come from pancreatic cancer work.
  They were kept unchanged so the two cancers are directly comparable; a breast
  specific CAF definition would be a separate sensitivity analysis.
* The `Epithelial` module includes CDH1 (E-cadherin), which is typically lost in
  invasive lobular carcinoma. The histological subtype of each section is not in
  the sample table, so this was not checked.
* `RND_MIN_K <- 15` (section 16) was kept. It now counts patients: 15 of 21.
* The circularity rule is unchanged: any module containing P4HA1, P4HB or a direct
  enzymatic partner is never a predictor. Sections 11.0, 15.1, 18.1 and 24.2
  still re-audit this and halt on failure.

## Repository layout

```
pipeline/p4ha1_ecm_pipeline_brca.R   the complete pipeline, sections 0 to 26
docs/brca_visium_samples.csv         the 54 BRCA sections, patient, primary/metastasis
docs/section_map.md                  what each section does, with line numbers
docs/data_setup.md                   where to put the data and the reused inputs
```

## Running it

1. Install the packages below.
2. Put the data in place as described in `docs/data_setup.md`.
3. From the repository root, in a fresh R session:

```r
Sys.setenv(P4HA1_BRCA_ROOT = "/path/to/p4ha1_brca")
source("pipeline/p4ha1_ecm_pipeline_brca.R")
```

The script runs top to bottom, as in PDAC. Later sections read objects and run
state flags defined earlier, and `assert_run_state()` halts if a section is
re-entered below its own settings block.

### Required packages

`Seurat`, `SeuratObject`, `ggplot2`, `dplyr`, `tidyr`, `patchwork`

### Optional packages

`viridis`, `spdep`, `ppcor`, `metafor` (three level model in section 21),
`synapser` (only for manifest mode).

## Known limitations

* **Spatial autocorrelation is not modelled in the regressions.** Model confidence
  intervals are too narrow. Sections 16 (random genes) and 26 (spatial permutation)
  measure how much; prefer their empirical percentiles over model p values.
* **The within patient collapse treats serial sections as independent**, so each
  patient's variance is too small, although each patient still counts once for k.
  The three level model in section 21 does not make that assumption.
* **mRNA, not protein.** All results are transcript level associations.
* **Section 20** depends on per spot ESTIMATE output inside each Seurat object.
  This was verified for PDAC only; sections without it are skipped and logged.

## Reproducibility

Unchanged from PDAC: `SEED <- 1234` is re-applied before every `AddModuleScore`,
outputs from section 13 onward carry a run id, and `ADJUST_DEPTH <- TRUE` writes
`_depthadj` files next to the unadjusted ones.
