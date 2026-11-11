# Data setup

Everything the pipeline reads hangs off one folder, set with the environment
variable `P4HA1_BRCA_ROOT` (default `~/p4ha1_brca`). Outputs are written under the
same folder.

```
<root>/
  data/
    seurat/                    the BRCA *-SeuratObj.rds files      (required)
    manifest.csv               Synapse manifest                    (optional)
    estimate_audit/
      ESTIMATE_signatures_used.csv                                 (section 20)
    upr_sets/                  MSigDB UPR gene set TSVs            (section 24)
  cohort_primary/              output of the main run (created)
  cohort_with_metastasis/      output when INCLUDE_METASTASIS <- TRUE (created)
```

## 1. The Seurat objects

The PDAC run used the Level 4 Visium Seurat objects of the HTAN WUSTL atlas. The
BRCA objects come from the same place: the HTAN Data Portal
(https://humantumoratlas.org), atlas HTAN WUSTL, Visium, Level 4, cancer type
breast. Access goes through Synapse, as for PDAC.

Put the `*-SeuratObj.rds` files in `data/seurat/`. The file name without
`-SeuratObj.rds` is matched to `sample_id` in `docs/brca_visium_samples.csv`:
first exactly, then by prefix in either direction (HTAN file names can carry a
suffix). A file that matches nothing halts the run, so no unlisted section enters
the cohort unnoticed; set `ALLOW_UNLISTED_FILES <- TRUE` to skip such files instead.
A file that matches more than one section also halts the run.

Checked against the Level 4 manifest `syn61518654` (132 files, all cancer types):
all 54 BRCA sections have a Level 4 Seurat object (`File Format` RData, `Visium
Workflow Type` Seurat4), matched on `HTAN Parent Biospecimen ID`, and every file
name stem equals the `sample_id` in the table exactly. 38 are OCT and 16 FFPE. The
run log still lists any expected section that is not in `data/seurat/`.

Portal filters that give exactly these 54 files: Atlas HTAN WUSTL, Assay 10X
Visium, Level 4, Organ Breast.

## 2. Inputs reused from the PDAC project

These two are gene lists, not patient data, so the files from the PDAC project can
be copied over unchanged.

* `ESTIMATE_signatures_used.csv` (section 20). In PDAC it was produced by
  `estimate_overlap_audit.R`, which is not in either repository. Copy the file from
  `~/Downloads/estimate_audit/` of the PDAC setup. Without it, section 20 stops.
* `upr_sets/` (section 24): the MSigDB UPR gene set TSVs used in PDAC. Copy the
  folder from `~/Downloads/upr_sets/` of the PDAC setup.

## 3. Optional manifest

`syn61518654.csv` is the Level 4 manifest of the whole WUSTL atlas (132 files,
all six cancer types), the same file the PDAC run used. Save only its 54 BRCA rows
as `data/manifest.csv`. In folder mode (the default) it is only used to label each
section OCT or FFPE; if it is absent those columns stay empty. In manifest mode the
full 132 row file would halt the run, because the other cancers are not in the
BRCA sample table.

## 4. The metastasis sensitivity run

Set `INCLUDE_METASTASIS <- TRUE` in section 0 and run again. Output goes to
`cohort_with_metastasis/`, so the main run is not touched. Then run section 22
with `INCLUDE_METASTASIS <- TRUE` to compare the pooled result with and without
patient HT268B1.
