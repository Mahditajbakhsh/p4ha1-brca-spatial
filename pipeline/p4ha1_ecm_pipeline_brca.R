################################################################################
# P4HA1 / ECM / hypoxia in spatial transcriptomics  --  BREAST CANCER (BRCA)
# cohort pipeline v6
#
# COHORT AND DESIGN
# (B1) 54 Visium sections from 22 patients, HTAN WUSTL atlas (Mo et al.,
#      Nature 2024, doi:10.1038/s41586-024-08087-4). Sections, patients and
#      primary / metastasis status come from docs/brca_visium_samples.csv
#      (Sample_ID_Lookup_table_v1.xlsx, sheet Visium,
#      github.com/ding-lab/ST_subclone_publication).
# (B2) PRIMARY TUMOURS ONLY in the main run. The 5 metastatic sections (all
#      from HT268B1) enter only with INCLUDE_METASTASIS <- TRUE, which writes
#      to a separate output folder.
# (B3) THE PATIENT IS THE POOLING UNIT. Serial sections of one patient are
#      collapsed (fixed effect, inverse variance), then patients are pooled
#      with DerSimonian-Laird. POOL_LEVEL <- "section" pools sections instead;
#      section 21 reports it as a sensitivity analysis.
# (B4) ONE PROJECT ROOT. All paths are derived from P4HA1_BRCA_ROOT.
#
# CELL-TYPE AND HYPOXIA MODULES (breast-specific; section 1)
#   Hypoxia     Buffa hypoxia metagene, MSigDB BUFFA_HYPOXIA_METAGENE (M34030,
#               50 genes) minus P4HA1, which is the outcome: 49 genes.
#               Buffa et al., Br J Cancer 2010, doi:10.1038/sj.bjc.6605450;
#               derived in breast and head and neck cancer. Among the best
#               performing signatures in tumours in Di Giovannantonio et al.,
#               Cell Genomics 2025, doi:10.1016/j.xgen.2025.100764.
#   iCAF, mCAF  Breast CAF classification, Cords et al., Nat Commun 2023,
#               doi:10.1038/s41467-023-39762-1, Table 1 (named genes only;
#               "collagens" in the mCAF entry is not expanded).
#   Epithelial  EPCAM, luminal keratins KRT8 / KRT18 / KRT19, and the basal /
#               myoepithelial keratins KRT5 / KRT14 (Reed et al., Nat Genet
#               2024, doi:10.1038/s41588-024-01688-9). CDH1 is not used: it is
#               lost in invasive lobular carcinoma.
#   The collagen, ECM, ER stress, ERAD and metabolic modules describe P4HA1
#   biology, not a cancer type; their sources are given where they are defined.
#
# ---- CHANGES IN v6 -----------------------------------------------------------
#   Breast-specific Hypoxia, iCAF, mCAF (replaces mCAF) and Epithelial
#   modules as above. Section 9 tests ENO1 inside the new hypoxia module.
#   Section 16 requires a gene to be pooled over two thirds of the units.
#   Figures: axis ranges no longer delete error bars (coord_cartesian),
#   Figure 3 axis label, Figure 4 "effects differ", Figure 5 without the
#   circular CollagenMaturation module, Figure 7 wording on P4HA3 (an active
#   collagen P4H, Kukkola et al., J Biol Chem 2003, doi:10.1074/jbc.M306806200).
#
# ---- HOW TO RUN THE DEPTH CHECK ----------------------------------------------
#   1. Run once with ADJUST_DEPTH <- FALSE (block 12b)  (files have no suffix)
#   2. Set ADJUST_DEPTH <- TRUE and re-run, fresh session, from the top
#   3. Section 17 prints the comparison
#
# ---- KNOWN LIMITATIONS -------------------------------------------------------
#   Spatial autocorrelation between spots is not modelled, so every model CI
#   is too narrow. Sections 16 and 26 measure how large that is. Quote
#   empirical percentiles rather than model p values.
#   The fixed-effect collapse inside a patient treats serial sections as
#   independent, so each patient's variance is too small. Section 21 reports
#   a three-level model when metafor is installed.
#   The Buffa hypoxia metagene contains proliferation genes (e.g. KIF20A,
#   ANLN, KIF4A), so part of its score can track proliferation.
#   The mCAF module has three genes; POSTN and MMP11 are also in the
#   descriptive ECM modules (never in the same model).
#
# Run top to bottom in a fresh R session.
################################################################################

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject)
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
})
if (requireNamespace("viridis", quietly = TRUE)) suppressPackageStartupMessages(library(viridis))

## ============================================================================
## 0. CONFIG
## ============================================================================
## Expected layout under the root:
##   <root>/data/seurat/      the *-SeuratObj.rds files
##   <root>/data/manifest.csv optional Synapse manifest
##   <root>/data/upr_sets/    MSigDB UPR gene set TSVs (section 24)
##   <root>/data/estimate_audit/ESTIMATE_signatures_used.csv (section 20)
PROJECT_ROOT <- Sys.getenv("P4HA1_BRCA_ROOT", unset = path.expand("~/p4ha1_brca"))

SAMPLE_TABLE <- Sys.getenv("P4HA1_BRCA_SAMPLES", unset = file.path("docs", "brca_visium_samples.csv"))

## main run = primary tumours only. TRUE adds the 5 metastatic sections
## (one patient) and writes to a different output folder.
INCLUDE_METASTASIS <- FALSE

## pooling unit. "patient" is the primary analysis, "section" is the
## sensitivity analysis that pools sections directly.
POOL_LEVEL <- "patient"
stopifnot(POOL_LEVEL %in% c("patient", "section"))
K_UNIT <- if (POOL_LEVEL == "patient") "patients" else "sections"

## a downloaded .rds that is not in the sample table halts the run unless TRUE
ALLOW_UNLISTED_FILES <- FALSE

INPUT_MODE   <- "folder"   # "folder" = every *-SeuratObj.rds in DOWNLOAD_DIR
# "manifest" = drive from MANIFEST_CSV + entityIds

DOWNLOAD_DIR <- file.path(PROJECT_ROOT, "data", "seurat")
OUT_ROOT     <- file.path(PROJECT_ROOT, paste0(
  if (INCLUDE_METASTASIS) "cohort_with_metastasis" else "cohort_primary",
  if (POOL_LEVEL == "section") "_sectionpool" else "",
  Sys.getenv("P4HA1_BRCA_OUT_SUFFIX", unset = "")))   # e.g. "_v5" to keep earlier runs
MANIFEST_CSV <- file.path(PROJECT_ROOT, "data", "manifest.csv")  # skipped if absent

# only used in manifest mode
INCLUDE_SAMPLE_IDS <- character(0)
INCLUDE_ENTITY_IDS <- character(0)
USE_SYNAPSER       <- FALSE

SAVE_PER_SAMPLE_FIGS <- TRUE    # FALSE on re-runs: hundreds of PNGs, slow
SAVE_ANALYZED_RDS    <- FALSE
SPOT_SIZE            <- 250

## One seed for the whole run. AddModuleScore samples its control genes, so
## set.seed(SEED) is called immediately before EVERY AddModuleScore.
SEED <- 1234
set.seed(SEED)

## NOTE: ADJUST_DEPTH, RUN_ID, stamp() and outf() are defined in block 12b,
## immediately before section 13. Do not define them here as well.


## ============================================================================
## 1. Gene panel + modules
## ============================================================================
## v6: breast-specific cell-type and hypoxia modules; sources in the header.
BUFFA_HYPOXIA_M34030 <- c(   # MSigDB BUFFA_HYPOXIA_METAGENE, verbatim, 50 genes
  "NDRG1","LDHA","PGK1","TPI1","PGAM1","PFKP","ENO1","SLC2A1","PNP","BNIP3",
  "PSRC1","SLC16A1","DDIT4","ADM","HK2","SEC61G","AK4","CA9","P4HA1","ACOT7",
  "GPI","TUBB6","TUBA1C","CDKN3","CTSV","VEGFA","TUBA1A","LRRC42","PSMA7","GAPDH",
  "CHCHD2","YKT6","MIF","MAP7D1","MRPL15","MRPL13","MCTS1","UTP11","KIF4A","HILPDA",
  "MRGBP","KIF20A","MRPS17","ESRP1","SHCBP1","SLC25A32","CORO1C","ANLN","MAD2L2","ANKRD37")
stopifnot(length(BUFFA_HYPOXIA_M34030) == 50, !anyDuplicated(BUFFA_HYPOXIA_M34030))

MODULES <- list(
  CollagenSynthesis  = c("COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3"),
  CollagenMaturation = c("P4HA1","P4HA2","P4HB","PLOD1","PLOD2","PLOD3"),
  ECMCrosslink       = c("LOX","LOXL1","LOXL2","LOXL3","LOXL4"),
  ECMGlycoprotein    = c("FN1","POSTN","SPARC","FBN1","THBS2"),
  ECMRemodel         = c("MMP2","MMP11","MMP14","TIMP1","TIMP2"),
  mCAF               = c("MMP11","CDH11","POSTN"),                            # Cords 2023, Table 1
  iCAF               = c("PLA2G2A","CFD","C3","CD34","CXCL12","CD248"),       # Cords 2023, Table 1
  Hypoxia            = setdiff(BUFFA_HYPOXIA_M34030, "P4HA1"),                # P4HA1 is the outcome
  Epithelial         = c("EPCAM","KRT8","KRT18","KRT19","KRT5","KRT14"))      # Reed 2024

## genes whose presence is logged per section: P4HA1 plus every module gene
gene_panel <- unique(c("P4HA1", unlist(MODULES, use.names = FALSE)))

TARGET_GENES   <- c("COL1A1","COL1A2","COL3A1","COL5A1","COL6A1","FN1","POSTN","SPARC","LOX","LOXL2")
MODULE_TARGETS <- c("CollagenSynthesis","CollagenMaturation","ECMCrosslink",
                    "ECMGlycoprotein","Hypoxia","mCAF","iCAF","Epithelial")
LOGIT_MODULES  <- c("Hypoxia","Epithelial","mCAF","CollagenSynthesis","iCAF")

## CollagenMaturation contains P4HA1 (and P4HA2), so it is circular as a
## predictor. Scored and plotted descriptively only; never in a model of
## P4HA1 or P4HA2. Section 15.1 re-checks this and halts on failure.

## ============================================================================
## 2. Helpers
## ============================================================================
dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)
LOG_FILE <- file.path(OUT_ROOT, "run_log.txt")
log_msg <- function(...) {
  m <- paste0(format(Sys.time(), "%H:%M:%S"), "  ", sprintf(...))
  cat(m, "\n"); cat(m, "\n", file = LOG_FILE, append = TRUE)
}
fig <- function(expr, file, w, h) {
  tryCatch({ p <- force(expr); ggsave(file, p, width = w, height = h, dpi = 300) },
           error = function(e) log_msg("  fig failed (%s): %s", basename(file), conditionMessage(e)))
}
get_image_name <- function(obj) { im <- Images(obj); if (length(im)) im[1] else NA_character_ }

set_spot_size <- function(obj, value = SPOT_SIZE) {
  img <- get_image_name(obj)
  if (!is.na(img)) try(obj@images[[img]]@scale.factors$spot <- value, silent = TRUE)
  obj
}
get_coords_aligned <- function(obj, img) {   # works on Seurat v4 and v5
  co <- GetTissueCoordinates(obj, image = img)
  if ("cell" %in% colnames(co)) rownames(co) <- co$cell
  if (all(c("imagerow","imagecol") %in% colnames(co)))
    data.frame(x = co$imagecol, y = co$imagerow, row.names = rownames(co))
  else if (all(c("x","y") %in% colnames(co)))
    data.frame(x = co$x, y = co$y, row.names = rownames(co))
  else { num <- co[, sapply(co, is.numeric), drop = FALSE]
  data.frame(x = num[[2]], y = num[[1]], row.names = rownames(co)) }
}

## AddModuleScore with a bin fallback. In a few sparse BRCA sections hundreds
## of genes share one identical mean and 20 bins cannot be formed
## ("Insufficient data values to produce 20 bins"). 20 bins are tried first
## (same seed, so every section that worked before gives the identical score);
## only on that error 15 and then 10 bins are tried. Every fallback is logged
## to <OUT_ROOT>/nbin_fallback_log.csv. Any other error is raised unchanged.
NBIN_STEPS <- c(20, 15, 10)
add_module_score_fallback <- function(object, features, name, assay, nbin = 20, ctrl = 50) {
  sid_here <- get0("sid", envir = parent.frame(),
                   ifnotfound = get0("sample_id", envir = parent.frame(), ifnotfound = NA_character_))
  last_err <- NULL
  for (nb in NBIN_STEPS[NBIN_STEPS <= nbin]) {
    set.seed(SEED)
    res <- tryCatch(AddModuleScore(object, features = features, name = name,
                                   assay = assay, nbin = nb, ctrl = ctrl),
                    error = function(e) e)
    if (!inherits(res, "error")) {
      if (nb != nbin) {
        log_msg("  %s / %s: %d bins impossible (tied mean expression), scored with %d bins",
                sid_here, name, nbin, nb)
        fb <- file.path(OUT_ROOT, "nbin_fallback_log.csv")
        suppressWarnings(write.table(data.frame(time = format(Sys.time()), sample_id = sid_here, module = name,
                                                nbin_requested = nbin, nbin_used = nb),
                                     fb, sep = ",", row.names = FALSE, col.names = !file.exists(fb), append = TRUE))
      }
      return(res)
    }
    last_err <- res
    if (!grepl("bins", conditionMessage(res))) break
  }
  log_msg("  %s / %s: AddModuleScore failed: %s", sid_here, name, conditionMessage(last_err))
  stop(last_err)
}

add_modules_safe <- function(obj, modules) {
  present <- character(0)
  for (nm in names(modules)) {
    g <- intersect(modules[[nm]], rownames(obj))
    if (length(g) < 2) next
    ok <- tryCatch({
      set.seed(SEED)
      obj <- add_module_score_fallback(obj, features = list(g), name = nm,
                                       assay = DefaultAssay(obj), nbin = 20, ctrl = 50); TRUE },
      error = function(e) FALSE)
    if (!ok) next
    old <- paste0(nm, "1")
    if (old %in% colnames(obj@meta.data)) {
      names(obj@meta.data)[names(obj@meta.data) == old] <- nm
      present <- c(present, nm)
    }
  }
  list(obj = obj, modules = present)
}

## The sample table. Columns: sample_id, patient_id, htan_biospecimen_id,
## htan_participant_id, primary_metastasis, cancer_type.
if (!file.exists(SAMPLE_TABLE))
  stop("Sample table not found: ", SAMPLE_TABLE,
       ". Source the script from the repo root or set P4HA1_BRCA_SAMPLES.")
SAMPLES <- read.csv(SAMPLE_TABLE, stringsAsFactors = FALSE)
stopifnot(all(c("sample_id", "patient_id", "primary_metastasis") %in% names(SAMPLES)),
          !anyDuplicated(SAMPLES$sample_id),
          all(SAMPLES$primary_metastasis %in% c("Primary", "Metastasis")))

## section -> patient. Extended in section 4 with the file stems.
SAMPLE_PATIENT <- setNames(SAMPLES$patient_id, SAMPLES$sample_id)
patient_of <- function(sid) {
  p <- unname(SAMPLE_PATIENT[as.character(sid)])
  if (anyNA(p)) stop("patient_of(): no patient for section(s) ",
                     paste(unique(sid[is.na(p)]), collapse = ", "))
  p
}

## v5 (C4): the unit a descriptive summary counts. With POOL_LEVEL "patient"
## medians, box plots and paired tests use ONE value per patient (the mean of
## that patient's sections). With "section" the unit is the section.
UNIT_LABEL  <- if (identical(POOL_LEVEL, "patient")) "patient" else "section"
unit_of     <- function(sid) if (identical(POOL_LEVEL, "patient")) patient_of(sid) else as.character(sid)
n_units     <- function(sid) length(unique(unit_of(sid)))
unit_values <- function(x, sid) {
  ok <- is.finite(x)
  if (!any(ok)) return(numeric(0))
  as.numeric(tapply(x[ok], unit_of(sid[ok]), mean))
}
unit_median <- function(x, sid) { v <- unit_values(x, sid); if (length(v)) median(v) else NA_real_ }

## fixed effect, inverse variance collapse of the sections of each patient
collapse_to_patients <- function(yi, vi, sid) {
  pid <- patient_of(sid)
  ok  <- is.finite(yi) & is.finite(vi) & vi > 0
  yi <- yi[ok]; vi <- vi[ok]; pid <- pid[ok]
  if (!length(yi)) return(list(yi = numeric(0), vi = numeric(0), n_sec = 0L))
  w  <- 1 / vi
  sw <- tapply(w, pid, sum)
  list(yi = as.numeric(tapply(w * yi, pid, sum) / sw),
       vi = as.numeric(1 / sw),
       n_sec = length(yi))
}

## DerSimonian-Laird random-effects pooling. sid (section IDs, same order as
## yi) is required when POOL_LEVEL is "patient". k counts patients (or
## sections); k_sections always counts the sections that went in.
meta_pool <- function(yi, vi, sid = NULL) {
  yi <- as.numeric(yi); vi <- as.numeric(vi)
  n_sec <- sum(is.finite(yi) & is.finite(vi) & vi > 0)
  if (identical(POOL_LEVEL, "patient")) {
    if (is.null(sid) || length(sid) != length(yi))
      stop("meta_pool(): POOL_LEVEL is 'patient' but section IDs were not passed ",
           "(or have the wrong length). Every pooling call must pass sample_id.")
    cp <- collapse_to_patients(yi, vi, as.character(sid)); yi <- cp$yi; vi <- cp$vi
  }
  keep <- is.finite(yi) & is.finite(vi) & vi > 0; yi <- yi[keep]; vi <- vi[keep]
  k <- length(yi)
  if (k == 0) return(data.frame(est=NA, lower=NA, upper=NA, k=0, tau2=NA, I2=NA, k_sections=n_sec))
  if (k == 1) return(data.frame(est=yi, lower=yi-1.96*sqrt(vi), upper=yi+1.96*sqrt(vi), k=1, tau2=0, I2=0,
                                k_sections=n_sec))
  w <- 1/vi; yb <- sum(w*yi)/sum(w); Q <- sum(w*(yi-yb)^2); df <- k-1
  C <- sum(w) - sum(w^2)/sum(w); tau2 <- max(0, (Q-df)/C)
  ws <- 1/(vi+tau2); est <- sum(ws*yi)/sum(ws); se <- sqrt(1/sum(ws))
  data.frame(est=est, lower=est-1.96*se, upper=est+1.96*se, k=k, tau2=tau2,
             I2=if (Q > 0) max(0,(Q-df)/Q)*100 else 0, k_sections=n_sec)
}

## ============================================================================
## 3. PER-SECTION analysis
## ============================================================================
analyze_sample <- function(rds_path, sample_id, preservation, sr_version,
                           entity_id, out_dir, save_figs) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  obj <- readRDS(rds_path)
  obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
  DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else Assays(obj)[1]
  obj <- set_spot_size(obj, SPOT_SIZE)
  log_msg("  class=%s assays=%s images=%s spots=%d feats=%d",
          paste(class(obj), collapse="/"), paste(Assays(obj), collapse=","),
          paste(Images(obj), collapse=","), ncol(obj), nrow(obj))
  
  if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
    alt <- grep("snn_res", colnames(obj@meta.data), value = TRUE)
    obj$seurat_clusters <- if (length(alt)) obj@meta.data[[alt[1]]] else factor(0)
  }
  if (!"P4HA1" %in% rownames(obj)) stop("P4HA1 not in this object")
  
  present <- gene_panel %in% rownames(obj)
  log_msg("  genes present: %d / %d  (missing: %s)", sum(present), length(gene_panel),
          paste(gene_panel[!present], collapse=", "))
  
  ms <- add_modules_safe(obj, MODULES); obj <- ms$obj; mods_present <- ms$modules
  
  genes_here <- intersect(c("P4HA1", TARGET_GENES), rownames(obj))
  mods_here  <- intersect(MODULE_TARGETS, mods_present)
  meta_here  <- intersect(c(mods_here, "seurat_clusters"), colnames(obj@meta.data))
  expr <- FetchData(obj, vars = c(genes_here, meta_here), layer = "data")
  targets_here <- intersect(TARGET_GENES, colnames(expr))
  
  ## --- spatial visualisation ---
  if (save_figs) {
    fig(SpatialDimPlot(obj, group.by="seurat_clusters", label=TRUE, label.size=4),
        file.path(out_dir,"01_clusters.png"), 7, 6)
    f <- intersect(c("P4HA1","COL1A1","COL1A2","COL3A1"), rownames(obj))
    if (length(f)) fig(SpatialFeaturePlot(obj, features=f, ncol=2, image.alpha=.5),
                       file.path(out_dir,"02_p4ha1_collagens.png"), 9, 8)
    f <- intersect(c("FN1","POSTN","SPARC","LOX"), rownames(obj))
    if (length(f)) fig(SpatialFeaturePlot(obj, features=f, ncol=2, image.alpha=.5),
                       file.path(out_dir,"03_ecm_glycoproteins.png"), 9, 8)
    f <- intersect(c("CollagenMaturation","CollagenSynthesis","ECMCrosslink","Hypoxia"), mods_present)
    if (length(f)) fig(SpatialFeaturePlot(obj, features=f, ncol=2, image.alpha=.5),
                       file.path(out_dir,"04_modules.png"), 9, 8)
    f <- intersect(c("P4HA1","mCAF","Hypoxia","Epithelial"), c(rownames(obj), mods_present))
    if (length(f)) fig(SpatialFeaturePlot(obj, features=f, ncol=2, image.alpha=.5),
                       file.path(out_dir,"05_p4ha1_compartments.png"), 9, 8)
  }
  
  ## --- per-spot correlations ---
  global_cor <- do.call(rbind, lapply(targets_here, function(g) {
    ct <- suppressWarnings(cor.test(expr$P4HA1, expr[[g]], method="spearman", exact=FALSE))
    data.frame(gene=g, rho=unname(ct$estimate), p=ct$p.value) }))
  global_cor <- global_cor[order(-global_cor$rho), ]
  write.csv(global_cor, file.path(out_dir,"06_global_correlation_P4HA1.csv"), row.names=FALSE)
  
  cor_by_cluster <- expr %>% group_by(seurat_clusters) %>%
    summarise(across(all_of(targets_here),
                     ~ suppressWarnings(cor(P4HA1, .x, method="spearman"))),
              n_spots = n(), .groups="drop")
  write.csv(cor_by_cluster, file.path(out_dir,"07_correlation_by_cluster.csv"), row.names=FALSE)
  
  if (save_figs) {
    cl_long <- cor_by_cluster %>% pivot_longer(all_of(targets_here), names_to="gene", values_to="rho")
    fig(ggplot(cl_long, aes(gene, factor(seurat_clusters), fill=rho)) + geom_tile() +
          scale_fill_gradient2(low="#3B4CC0", mid="white", high="#B40426", midpoint=0, limits=c(-1,1)) +
          labs(x=NULL, y="Cluster", title="P4HA1 vs ECM genes (per cluster)") +
          theme_minimal(base_size=11) + theme(axis.text.x=element_text(angle=45, hjust=1)),
        file.path(out_dir,"08_correlation_heatmap.png"), 8, 5)
  }
  
  ## --- P4HA1 vs module scores ---
  module_cor <- do.call(rbind, lapply(mods_here, function(m) {
    ct <- suppressWarnings(cor.test(expr$P4HA1, expr[[m]], method="spearman", exact=FALSE))
    data.frame(module=m, rho=unname(ct$estimate), p=ct$p.value) }))
  if (!is.null(module_cor)) module_cor <- module_cor[order(-module_cor$rho), ]
  write.csv(module_cor, file.path(out_dir,"09_P4HA1_vs_modules.csv"), row.names=FALSE)
  
  if (save_figs) {
    mk <- function(y) ggplot(expr, aes(P4HA1, .data[[y]], colour=factor(seurat_clusters))) +
      geom_point(size=.6, alpha=.6) + geom_smooth(method="loess", se=FALSE, colour="black", linewidth=.6) +
      labs(colour="Cluster", title=paste("P4HA1 vs", y)) + theme_minimal(base_size=11)
    sc <- intersect(c("CollagenSynthesis","Hypoxia","mCAF","Epithelial"), colnames(expr))
    if (length(sc) >= 1)
      fig(wrap_plots(lapply(sc, mk), ncol=2), file.path(out_dir,"10_scatter_P4HA1_modules.png"), 11, 8)
  }
  
  ## --- Moran's I ---
  moran <- NULL
  if (requireNamespace("spdep", quietly=TRUE)) {
    img <- get_image_name(obj)
    moran <- tryCatch({
      coords <- get_coords_aligned(obj, img)
      lw <- spdep::nb2listw(spdep::knn2nb(spdep::knearneigh(as.matrix(coords), k=6)), style="W")
      gm <- intersect(c("P4HA1","COL1A1","COL1A2","COL3A1","FN1","POSTN","SPARC",
                        "LOX","LOXL2","CA9","VEGFA","EPCAM","ACTA2"), rownames(obj))
      ## v5.1: one gene failing (e.g. constant across spots) no longer drops the
      ## whole table; it is logged and skipped.
      tab <- do.call(rbind, lapply(gm, function(g) {
        v <- FetchData(obj, vars=g, layer="data")[,1]; names(v) <- colnames(obj)
        v <- v[rownames(coords)]
        if (anyNA(v) || stats::var(v) <= 0) {
          log_msg("  Moran: %s skipped (constant or missing values)", g); return(NULL) }
        mi <- tryCatch(spdep::moran.test(v, lw, zero.policy=TRUE, randomisation=TRUE),
                       error = function(e) { log_msg("  Moran: %s skipped: %s", g, conditionMessage(e)); NULL })
        if (is.null(mi)) return(NULL)
        data.frame(gene=g, morans_I=unname(mi$estimate[1]), expected=unname(mi$estimate[2]),
                   variance=unname(mi$estimate[3]), p_value=mi$p.value) }))
      if (is.null(tab)) stop("no gene could be tested")
      tab[order(-tab$morans_I), ]
    }, error=function(e){ log_msg("  Moran skipped: %s", conditionMessage(e)); NULL })
    if (!is.null(moran)) write.csv(moran, file.path(out_dir,"15_morans_I.csv"), row.names=FALSE)
  }
  
  ## --- P4HA1-high vs -low DE ---
  p4 <- expr$P4HA1; pos <- p4 > 0
  st <- rep("low", length(p4))
  if (sum(pos) >= 20) { thr <- median(p4[pos]); st[pos & p4 < thr] <- "mid"; st[pos & p4 >= thr] <- "high" }
  else st[pos] <- "high"
  obj$P4HA1_status <- factor(st, levels=c("low","mid","high"))
  if (save_figs)
    fig(SpatialDimPlot(obj, group.by="P4HA1_status",
                       cols=c(low="grey85", mid="#FDB863", high="#B2182B")),
        file.path(out_dir,"11_p4ha1_status_map.png"), 7, 6)
  if (sum(st=="high") >= 10 && sum(st=="low") >= 10) {
    Idents(obj) <- "P4HA1_status"
    de <- tryCatch({
      if (DefaultAssay(obj)=="SCT") obj <- PrepSCTFindMarkers(obj)
      d <- FindMarkers(obj, ident.1="high", ident.2="low", assay=DefaultAssay(obj),
                       logfc.threshold=0.25, min.pct=0.1)
      d$gene <- rownames(d); d[order(-d$avg_log2FC), ]
    }, error=function(e){ log_msg("  DE skipped: %s", conditionMessage(e)); NULL })
    if (!is.null(de)) write.csv(de, file.path(out_dir,"12_DE_P4HA1high_vs_low.csv"), row.names=FALSE)
  }
  
  ## --- within-compartment correlations ---
  if (all(c("Epithelial","mCAF") %in% colnames(expr))) {
    expr$Compartment <- with(expr, ifelse(
      Epithelial > quantile(Epithelial,.66,na.rm=TRUE) & mCAF < quantile(mCAF,.5,na.rm=TRUE),
      "Epithelial-dominant",
      ifelse(mCAF > quantile(mCAF,.66,na.rm=TRUE) & Epithelial < quantile(Epithelial,.5,na.rm=TRUE),
             "CAF-dominant","Mixed")))
  } else expr$Compartment <- "Mixed"
  comp_cor <- expr %>% group_by(Compartment) %>%
    summarise(across(all_of(targets_here),
                     ~ suppressWarnings(cor(P4HA1, .x, method="spearman"))),
              n_spots_compartment = n(), .groups="drop")
  write.csv(comp_cor, file.path(out_dir,"16_correlation_by_compartment.csv"), row.names=FALSE)
  comp_long <- comp_cor %>% pivot_longer(all_of(targets_here), names_to="gene", values_to="rho")
  if (save_figs)
    fig(ggplot(comp_long, aes(gene, Compartment, fill=rho)) + geom_tile() +
          scale_fill_gradient2(low="#3B4CC0", mid="white", high="#B40426", midpoint=0, limits=c(-.5,.5)) +
          labs(x=NULL, y=NULL, title="P4HA1 vs ECM genes within compartments") +
          theme_minimal(base_size=11) + theme(axis.text.x=element_text(angle=45, hjust=1)),
        file.path(out_dir,"17_compartment_heatmap.png"), 8, 3.5)
  
  ## --- partial correlations ---
  partial_tbl <- NULL
  if (requireNamespace("ppcor", quietly=TRUE)) {
    pct <- function(x,y,z){ d<-data.frame(x,y,z); d<-d[complete.cases(d),]
    tryCatch(ppcor::pcor.test(d$x,d$y,d$z, method="spearman")$estimate, error=function(e) NA_real_) }
    partial_tbl <- do.call(rbind, lapply(targets_here, function(g) data.frame(
      gene=g,
      rho_unadjusted       = suppressWarnings(cor(expr$P4HA1, expr[[g]], method="spearman")),
      rho_given_Hypoxia    = if ("Hypoxia"    %in% colnames(expr)) pct(expr$P4HA1, expr[[g]], expr$Hypoxia)    else NA,
      rho_given_Epithelial = if ("Epithelial" %in% colnames(expr)) pct(expr$P4HA1, expr[[g]], expr$Epithelial) else NA,
      rho_given_mCAF       = if ("mCAF"      %in% colnames(expr)) pct(expr$P4HA1, expr[[g]], expr$mCAF)      else NA)))
    write.csv(partial_tbl, file.path(out_dir,"18_partial_correlations.csv"), row.names=FALSE)
  }
  
  ## --- compartment map ---
  obj$Compartment <- expr$Compartment[match(colnames(obj), rownames(expr))]
  if (save_figs)
    fig(SpatialDimPlot(obj, group.by="Compartment",
                       cols=c("Epithelial-dominant"="#2166AC","CAF-dominant"="#B2182B","Mixed"="grey80")),
        file.path(out_dir,"19_compartment_map.png"), 7, 6)
  
  ## --- per-cluster summary ---
  cl_summary <- expr %>% group_by(seurat_clusters) %>%
    summarise(n_spots=n(), P4HA1_mean=mean(P4HA1), P4HA1_pct_pos=mean(P4HA1>0)*100,
              across(all_of(mods_here), ~ mean(.x, na.rm=TRUE)), .groups="drop") %>%
    arrange(desc(P4HA1_mean))
  write.csv(cl_summary, file.path(out_dir,"20_cluster_summary.csv"), row.names=FALSE)
  
  ## --- dot plot ---
  if (save_figs) {
    Idents(obj) <- "seurat_clusters"
    mk <- intersect(c("P4HA1","P4HA2","COL1A1","COL1A2","COL3A1","FN1","POSTN","SPARC","LOX","LOXL2",
                      "ACTA2","TAGLN","VIM","CA9","VEGFA","SLC2A1","BNIP3","EPCAM","KRT8","KRT18"),
                    rownames(obj))
    if (length(mk))
      fig(DotPlot(obj, features=mk, assay=DefaultAssay(obj), group.by="seurat_clusters") +
            RotatedAxis() + labs(title="P4HA1 / ECM / Hypoxia / Epithelial markers by cluster"),
          file.path(out_dir,"14_dotplot.png"), 12, 6)
  }
  
  ## --- per-section logistic regression + forest ---
  y <- as.integer(expr$P4HA1 > 0); logit <- NULL
  if (length(unique(y)) == 2) {
    logit <- do.call(rbind, lapply(intersect(LOGIT_MODULES, colnames(expr)), function(m){
      dd <- data.frame(y = y, x = expr[[m]])
      f <- tryCatch(glm(y ~ x, data = dd, family = binomial()), error=function(e) NULL)
      if (is.null(f)) return(NULL)
      s <- summary(f)$coefficients; if (nrow(s) < 2) return(NULL)
      data.frame(module=m, log_OR=s[2,"Estimate"], SE=s[2,"Std. Error"], p=s[2,"Pr(>|z|)"]) }))
    if (!is.null(logit) && nrow(logit) > 0) {
      write.csv(logit, file.path(out_dir,"25_logistic_regression_results.csv"), row.names=FALSE)
      if (save_figs) {
        lf <- logit
        lf$OR <- exp(lf$log_OR); lf$lo <- exp(lf$log_OR-1.96*lf$SE); lf$hi <- exp(lf$log_OR+1.96*lf$SE)
        lf$direction <- ifelse(lf$OR>1, "Promotes P4HA1", "Suppresses P4HA1")
        lf$module <- factor(lf$module, levels=lf$module[order(lf$OR)])
        lf$label  <- sprintf("OR = %.2f  (p = %.0e)", lf$OR, lf$p)
        fig(ggplot(lf, aes(OR, module, colour=direction)) +
              geom_vline(xintercept=1, linetype="dashed", colour="grey50") +
              geom_errorbarh(aes(xmin=lo, xmax=hi), height=.25, linewidth=.9) +
              geom_point(size=5) +
              geom_text(aes(label=label), hjust=-.15, vjust=-1, size=3.6, colour="black", show.legend=FALSE) +
              scale_x_log10(breaks=c(.2,.5,1,2,5)) + coord_cartesian(xlim=c(.1,10)) +
              scale_colour_manual(values=c("Promotes P4HA1"="#B2182B","Suppresses P4HA1"="#2166AC")) +
              labs(x="Odds ratio for P4HA1 detection (log scale, 95% CI)", y=NULL, colour=NULL,
                   title="What predicts P4HA1 expression?",
                   subtitle=sprintf("%s, n = %s spots", sample_id, format(nrow(expr), big.mark=",")),
                   caption="Logistic regression: P(P4HA1 > 0) ~ module score") +
              theme_minimal(base_size=13) +
              theme(legend.position="top", panel.grid.minor=element_blank(),
                    plot.title=element_text(face="bold")),
            file.path(out_dir,"26_forest_plot.png"), 8.5, 4.5)
      }
    }
  }
  
  if (SAVE_ANALYZED_RDS) saveRDS(obj, file.path(out_dir, paste0(sample_id,"_analyzed.rds")))
  
  tag <- function(df){ if (is.null(df)||nrow(df)==0) return(NULL)
    df$sample_id<-sample_id; df$entity_id<-entity_id; df$preservation<-preservation
    df$sr_version<-sr_version; df$n_spots_sample<-nrow(expr); df }
  list(global_cor=tag(global_cor), module_cor=tag(module_cor),
       comp_cor=tag(comp_long), partial_cor=tag(partial_tbl),
       logit=tag(logit), cl_summary=tag(cl_summary), moran=tag(moran))
}

## ============================================================================
## 4. Build the work list, then RUN THE LOOP
## ============================================================================
mani <- NULL
if (nzchar(MANIFEST_CSV) && file.exists(MANIFEST_CSV)) {
  mani <- read.csv(MANIFEST_CSV, check.names=FALSE, stringsAsFactors=FALSE)
  mani$sample_id    <- sub("-SeuratObj\\.rds$", "", basename(mani$Filename))
  mani$preservation <- ifelse(grepl("/FFPE/", mani$Path), "FFPE",
                              ifelse(grepl("/OCT/", mani$Path), "OCT", "other"))
  mani$sr_version   <- if ("Workflow Version" %in% colnames(mani)) mani[["Workflow Version"]] else NA
}

if (INPUT_MODE == "folder") {
  files <- list.files(DOWNLOAD_DIR, pattern="-SeuratObj\\.rds$", full.names=TRUE)
  if (!length(files)) stop("No *-SeuratObj.rds files found in ", DOWNLOAD_DIR)
  work <- data.frame(path=files, sample_id=sub("-SeuratObj\\.rds$","",basename(files)),
                     stringsAsFactors=FALSE)
  if (!is.null(mani)) { m <- mani[match(work$sample_id, mani$sample_id), ]
  work$entity_id<-m$entityId; work$preservation<-m$preservation; work$sr_version<-m$sr_version
  } else { work$entity_id<-NA; work$preservation<-NA; work$sr_version<-NA }
  log_msg("FOLDER mode: %d Seurat objects in %s", nrow(work), DOWNLOAD_DIR)
} else {
  stopifnot(!is.null(mani), all(c("entityId","Filename","Path") %in% colnames(mani)))
  if (length(INCLUDE_SAMPLE_IDS) > 0) mani <- mani[mani$sample_id %in% INCLUDE_SAMPLE_IDS, ]
  if (length(INCLUDE_ENTITY_IDS) > 0) mani <- mani[mani$entityId  %in% INCLUDE_ENTITY_IDS, ]
  resolve_path <- function(eid, fn){
    if (USE_SYNAPSER){ if(!requireNamespace("synapser",quietly=TRUE)) stop("install synapser")
      return(synapser::synGet(eid, downloadLocation=DOWNLOAD_DIR)$path) }
    cand <- file.path(DOWNLOAD_DIR, basename(fn)); if (file.exists(cand)) return(cand)
    h <- list.files(DOWNLOAD_DIR, pattern=paste0("^",basename(fn),"$"), recursive=TRUE, full.names=TRUE)
    if (length(h)) h[1] else NA_character_ }
  work <- data.frame(
    path=vapply(seq_len(nrow(mani)), function(i) tryCatch(resolve_path(mani$entityId[i], mani$Filename[i]),
                                                          error=function(e) NA_character_), character(1)),
    sample_id=mani$sample_id, entity_id=mani$entityId,
    preservation=mani$preservation, sr_version=mani$sr_version, stringsAsFactors=FALSE)
  log_msg("MANIFEST mode: %d rows", nrow(work))
}

## ---------------------------------------------------------------------------
## Match every file to the sample table, attach patient and status, then keep
## primary tumours only unless INCLUDE_METASTASIS is TRUE.
##   1. exact match of the file stem to sample_id
##   2. otherwise a table ID that is a prefix of the file stem, or vice versa,
##      AND is followed by a delimiter (- _ .) in the longer one.
##      v5 (C6): without the delimiter rule "X-U1" matched "X-U12".
##   Anything still unmatched halts the run unless ALLOW_UNLISTED_FILES.
## ---------------------------------------------------------------------------
prefix_at_boundary <- function(long, short)
  startsWith(long, short) & nchar(long) > nchar(short) &
  grepl("^[-_.]", substring(long, nchar(short) + 1L))

match_section <- function(stem) {
  if (stem %in% SAMPLES$sample_id) return(stem)
  ids <- SAMPLES$sample_id
  hit <- ids[prefix_at_boundary(stem, ids) | prefix_at_boundary(ids, stem)]
  if (length(hit) == 1) return(hit)
  if (length(hit) > 1)
    stop("File ", stem, " matches several table sections: ", paste(hit, collapse = ", "))
  NA_character_
}
work$table_id <- vapply(work$sample_id, match_section, character(1))
unlisted <- work$sample_id[is.na(work$table_id)]
if (length(unlisted)) {
  log_msg("Files not in %s: %s", SAMPLE_TABLE, paste(unlisted, collapse = ", "))
  if (!ALLOW_UNLISTED_FILES)
    stop(length(unlisted), " file(s) are not in the BRCA sample table (see log). ",
         "Remove them from DOWNLOAD_DIR or set ALLOW_UNLISTED_FILES <- TRUE.")
  work <- work[!is.na(work$table_id), , drop = FALSE]
}
inexact <- work$sample_id != work$table_id
if (any(inexact))
  log_msg("Matched by prefix (check these): %s",
          paste(sprintf("%s -> %s", work$sample_id[inexact], work$table_id[inexact]), collapse = "; "))
if (anyDuplicated(work$table_id))
  stop("Two files map to the same table section: ",
       paste(unique(work$table_id[duplicated(work$table_id)]), collapse = ", "))
tbl <- SAMPLES[match(work$table_id, SAMPLES$sample_id), ]
work$patient_id         <- tbl$patient_id
work$primary_metastasis <- tbl$primary_metastasis
if ("cancer_type" %in% names(tbl) && !all(tbl$cancer_type == "BRCA"))
  stop("Non-BRCA sections in the work list: ",
       paste(work$sample_id[tbl$cancer_type != "BRCA"], collapse = ", "))

## the analysis keys on the file stem, so map the file stem to the patient too
SAMPLE_PATIENT <- c(SAMPLE_PATIENT, setNames(work$patient_id, work$sample_id))
SAMPLE_PATIENT <- SAMPLE_PATIENT[!duplicated(names(SAMPLE_PATIENT))]

n_met <- sum(work$primary_metastasis == "Metastasis")
if (!INCLUDE_METASTASIS && n_met) {
  log_msg("Primary-only run: %d metastatic section(s) left out (%s)", n_met,
          paste(unique(work$patient_id[work$primary_metastasis == "Metastasis"]), collapse = ", "))
  work <- work[work$primary_metastasis == "Primary", , drop = FALSE]
}

exp_tbl <- if (INCLUDE_METASTASIS) SAMPLES else SAMPLES[SAMPLES$primary_metastasis == "Primary", ]
missing_secs <- setdiff(exp_tbl$sample_id, work$table_id)
log_msg("Cohort: %d of %d expected sections, %d of %d expected patients.",
        nrow(work), nrow(exp_tbl), length(unique(work$patient_id)), length(unique(exp_tbl$patient_id)))
if (length(missing_secs)) log_msg("  Not found (not downloaded or no Level 4 object): %s",
                                  paste(missing_secs, collapse = ", "))
if (!nrow(work)) stop("No BRCA sections left to analyse.")

acc <- list(global_cor=list(), module_cor=list(), comp_cor=list(),
            partial_cor=list(), logit=list(), cl_summary=list(), moran=list())
status <- data.frame()
for (i in seq_len(nrow(work))) {
  sid<-work$sample_id[i]; eid<-work$entity_id[i]; path<-work$path[i]
  log_msg("[%d/%d] %s", i, nrow(work), sid)
  if (is.na(path) || !file.exists(path)) {
    log_msg("  FILE NOT FOUND -> skipped")
    status<-rbind(status, data.frame(sample_id=sid, entity_id=eid, status="file_not_found")); next }
  res <- tryCatch(analyze_sample(path, sid, work$preservation[i], work$sr_version[i], eid,
                                 out_dir=file.path(OUT_ROOT,"per_sample",sid), save_figs=SAVE_PER_SAMPLE_FIGS),
                  error=function(e){ log_msg("  ERROR: %s", conditionMessage(e)); NULL })
  if (is.null(res)) { status<-rbind(status, data.frame(sample_id=sid, entity_id=eid, status="error")); next }
  for (k in names(acc)) acc[[k]][[length(acc[[k]])+1]] <- res[[k]]
  status <- rbind(status, data.frame(sample_id=sid, entity_id=eid, status="ok"))
}

bind <- function(lst) if (length(lst)) dplyr::bind_rows(lst) else data.frame()
all_global<-bind(acc$global_cor); all_module<-bind(acc$module_cor); all_comp<-bind(acc$comp_cor)
all_partial<-bind(acc$partial_cor); all_logit<-bind(acc$logit)
all_cluster<-bind(acc$cl_summary); all_moran<-bind(acc$moran)

## ============================================================================
## 5. Save combined per-section tables + run status
## ============================================================================
agg <- file.path(OUT_ROOT,"aggregate"); dir.create(agg, showWarnings=FALSE, recursive=TRUE)
write.csv(all_global,  file.path(agg,"per_sample_global_correlations.csv"), row.names=FALSE)
write.csv(all_module,  file.path(agg,"per_sample_module_correlations.csv"), row.names=FALSE)
write.csv(all_comp,    file.path(agg,"per_sample_compartment_correlations.csv"), row.names=FALSE)
write.csv(all_partial, file.path(agg,"per_sample_partial_correlations.csv"), row.names=FALSE)
write.csv(all_logit,   file.path(agg,"per_sample_logistic.csv"), row.names=FALSE)
write.csv(all_cluster, file.path(agg,"per_sample_cluster_summary.csv"), row.names=FALSE)
write.csv(all_moran,   file.path(agg,"per_sample_morans_I.csv"), row.names=FALSE)
write.csv(status,      file.path(agg,"run_status.csv"), row.names=FALSE)
## the section -> patient map used for pooling, read by sections 21 and 22
write.csv(work[, c("sample_id", "table_id", "patient_id", "primary_metastasis")],
          file.path(agg, "cohort_sections.csv"), row.names = FALSE)
log_msg("Sections ok: %d / %d", sum(status$status=="ok"), nrow(status))

## A section that failed the per-section analysis (for example an object with
## no P4HA1) has no outcome for any later model. Drop it here, once and visibly.
failed_secs <- status$sample_id[status$status != "ok"]
if (length(failed_secs)) {
  log_msg("Removed from the work list for sections 6 onward (failed in section 3): %s",
          paste(failed_secs, collapse = ", "))
  work <- work[!work$sample_id %in% failed_secs, , drop = FALSE]
}
log_msg("Work list for the rest of the run: %d sections, %d patients.",
        nrow(work), length(unique(work$patient_id)))

## ============================================================================
## 6. CROSS-SECTION pooled forest
## ============================================================================
if (nrow(all_logit) > 0) {
  pooled <- all_logit %>% group_by(module) %>%
    group_modify(~ meta_pool(.x$log_OR, .x$SE^2, .x$sample_id)) %>% ungroup() %>%
    mutate(OR=exp(est), lower=exp(lower), upper=exp(upper),
           direction=ifelse(OR>1,"Promotes P4HA1","Suppresses P4HA1"),
           label=sprintf("OR %.2f  (%d %s, I2=%.0f%%)", OR, as.integer(k), K_UNIT, I2)) %>%
    arrange(OR) %>% mutate(module=factor(module, levels=module))
  write.csv(pooled, file.path(agg,"POOLED_logistic_forest.csv"), row.names=FALSE)
  fig(ggplot(pooled, aes(OR, module, colour=direction)) +
        geom_vline(xintercept=1, linetype="dashed", colour="grey50") +
        geom_errorbarh(aes(xmin=lower, xmax=upper), height=.25, linewidth=.9) +
        geom_point(size=5) +
        geom_text(aes(label=label), hjust=-.1, vjust=-1.1, size=3.4, colour="black", show.legend=FALSE) +
        scale_x_log10(breaks=c(.2,.5,1,2,5)) + coord_cartesian(xlim=c(.1,10)) +
        scale_colour_manual(values=c("Promotes P4HA1"="#B2182B","Suppresses P4HA1"="#2166AC")) +
        labs(x="Pooled odds ratio for P4HA1 detection (log scale, 95% CI)", y=NULL, colour=NULL,
             title="What predicts P4HA1 across tumours?",
             subtitle=sprintf("Random-effects meta-analysis of %d patients (%d sections)",
                              length(unique(patient_of(all_logit$sample_id))), length(unique(all_logit$sample_id))),
             caption="Per-section logistic regression, DerSimonian-Laird pooled") +
        theme_minimal(base_size=13) +
        theme(legend.position="top", panel.grid.minor=element_blank(), plot.title=element_text(face="bold")),
      file.path(agg,"FOREST_pooled.png"), 9, 5)
}

## ============================================================================
## 7. Consistency boxplot: P4HA1 <-> module correlations, one point per unit
##    The only figure showing SPREAD across tumours rather than a pooled
##    average. v5 (C4): one value per patient (mean over its sections).
## ============================================================================
if (nrow(all_module) > 0) {
  mod_unit <- all_module %>% filter(is.finite(rho)) %>%
    mutate(unit = unit_of(sample_id)) %>%
    group_by(module, unit) %>%
    summarise(rho = mean(rho), n_sections = dplyr::n(), .groups = "drop")
  mod_summary <- mod_unit %>% group_by(module) %>%
    summarise(median_rho = median(rho), frac_positive = mean(rho > 0),
              n_units = dplyr::n(), n_sections = sum(n_sections), .groups = "drop") %>%
    arrange(desc(median_rho))
  write.csv(mod_summary, file.path(agg,"module_correlation_summary.csv"), row.names=FALSE)
  ord <- mod_summary$module[order(mod_summary$median_rho)]
  ## v5.2: CollagenMaturation contains P4HA1 (circular); in the CSV, not drawn
  fig(ggplot(mod_unit[mod_unit$module != "CollagenMaturation", ],
             aes(factor(module, levels=setdiff(ord, "CollagenMaturation")), rho)) +
        geom_hline(yintercept=0, linetype="dashed", colour="grey60") +
        geom_boxplot(outlier.size=.5, fill="grey90") +
        geom_jitter(width=.15, size=.9, alpha=.5) + coord_flip() +
        labs(x=NULL, y=sprintf("Per-%s Spearman rho (P4HA1 vs module)", UNIT_LABEL),
             title="Consistency of P4HA1 associations across tumours",
             subtitle=sprintf("One point per %s (mean of its sections), n = %d",
                              UNIT_LABEL, dplyr::n_distinct(mod_unit$unit))) +
        theme_minimal(base_size=12),
      file.path(agg,"BOX_module_correlations.png"), 8, 5)
}

log_msg("Core cohort analysis done. Aggregate outputs in: %s", normalizePath(agg))

## ============================================================================
## 8. COHORT within-compartment analysis
##    v5 (C4): summaries and the paired CAF vs epithelial test use one value
##    per patient, so serial sections are not counted as independent pairs.
## ============================================================================
COMPARTMENT_MIN_SPOTS <- 20
COMPARTMENT_LEVELS    <- c("Epithelial-dominant","Mixed","CAF-dominant")
CONTRAST              <- c("CAF-dominant","Epithelial-dominant")

comp_csv <- file.path(agg, "per_sample_compartment_correlations.csv")

if (!file.exists(comp_csv)) {
  log_msg("Cohort compartment step skipped: %s not found", comp_csv)
} else {
  cc <- tryCatch(read.csv(comp_csv, stringsAsFactors = FALSE, check.names = FALSE),
                 error = function(e) data.frame())
  needed <- c("sample_id","Compartment","gene","rho")
  if (!all(needed %in% colnames(cc))) {
    log_msg("Cohort compartment step skipped: missing columns (%s)",
            paste(setdiff(needed, colnames(cc)), collapse = ", "))
  } else {
    cc <- cc[!is.na(cc$rho), , drop = FALSE]
    if ("n_spots_compartment" %in% colnames(cc)) {
      before <- nrow(cc)
      cc <- cc[is.finite(cc$n_spots_compartment) &
                 cc$n_spots_compartment >= COMPARTMENT_MIN_SPOTS, , drop = FALSE]
      log_msg("Compartment spot filter (>= %d spots): kept %d / %d rows",
              COMPARTMENT_MIN_SPOTS, nrow(cc), before)
    } else {
      log_msg("WARNING: n_spots_compartment not in CSV; rho from tiny compartments is NOT screened.")
    }
    cc$Compartment <- factor(cc$Compartment, levels = COMPARTMENT_LEVELS)
    cc <- cc[!is.na(cc$Compartment), , drop = FALSE]
    
    if (nrow(cc) == 0) {
      log_msg("Cohort compartment step skipped: no rows left after cleaning/filtering")
    } else {
      ## one value per unit, gene and compartment
      cc_u <- cc %>% mutate(unit = unit_of(sample_id)) %>%
        group_by(gene, Compartment, unit) %>%
        summarise(rho = mean(rho), n_sections = dplyr::n(), .groups = "drop")
      cc_u$gene <- factor(cc_u$gene, levels = sort(unique(as.character(cc_u$gene))))
      
      comp_summary <- cc_u %>%
        group_by(gene, Compartment) %>%
        summarise(n_units       = dplyr::n(),
                  n_sections    = sum(n_sections),
                  median_rho    = median(rho, na.rm = TRUE),
                  q25           = unname(quantile(rho, .25, na.rm = TRUE)),
                  q75           = unname(quantile(rho, .75, na.rm = TRUE)),
                  frac_positive = mean(rho > 0, na.rm = TRUE),
                  .groups = "drop") %>%
        arrange(gene, Compartment)
      write.csv(comp_summary, file.path(agg, "COHORT_compartment_summary.csv"), row.names = FALSE)
      
      fig(ggplot(cc_u, aes(Compartment, rho, fill = Compartment)) +
            geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
            geom_boxplot(outlier.size = .4, alpha = .85) +
            geom_jitter(width = .15, size = .6, alpha = .4) +
            facet_wrap(~ gene, ncol = 5) +
            scale_fill_manual(values = c("Epithelial-dominant" = "#2166AC",
                                         "Mixed"               = "grey75",
                                         "CAF-dominant"        = "#B2182B"), drop = FALSE) +
            labs(x = NULL, y = sprintf("Per-%s Spearman rho (P4HA1 vs gene)", UNIT_LABEL),
                 title = "P4HA1 to ECM correlation by tissue compartment",
                 subtitle = sprintf("One point per %s (%d); compartments with < %d spots removed",
                                    UNIT_LABEL, dplyr::n_distinct(cc_u$unit), COMPARTMENT_MIN_SPOTS)) +
            theme_minimal(base_size = 11) +
            theme(legend.position = "none",
                  axis.text.x = element_text(angle = 30, hjust = 1),
                  strip.text  = element_text(face = "bold")),
          file.path(agg, "BOX_compartment_correlations.png"), 12, 7)
      
      test_rows <- lapply(levels(cc_u$gene), function(g) {
        sub <- cc_u[cc_u$gene == g & as.character(cc_u$Compartment) %in% CONTRAST,
                    c("unit","Compartment","rho")]
        empty <- data.frame(gene = g, n_pairs = 0L, median_CAF = NA_real_,
                            median_Epi = NA_real_, median_diff = NA_real_,
                            W = NA_real_, p = NA_real_)
        if (nrow(sub) == 0) return(empty)
        sub$Compartment <- factor(as.character(sub$Compartment), levels = CONTRAST)
        m <- tapply(sub$rho, list(sub$unit, sub$Compartment), function(x) mean(x, na.rm = TRUE))
        if (!all(CONTRAST %in% colnames(m))) return(empty)
        caf <- m[, "CAF-dominant"]; epi <- m[, "Epithelial-dominant"]
        keep <- is.finite(caf) & is.finite(epi)
        caf <- caf[keep]; epi <- epi[keep]; n_pairs <- length(caf)
        if (n_pairs < 3)
          return(data.frame(gene = g, n_pairs = n_pairs,
                            median_CAF  = if (n_pairs) median(caf) else NA_real_,
                            median_Epi  = if (n_pairs) median(epi) else NA_real_,
                            median_diff = if (n_pairs) median(caf - epi) else NA_real_,
                            W = NA_real_, p = NA_real_))
        wt <- suppressWarnings(stats::wilcox.test(caf, epi, paired = TRUE))
        data.frame(gene = g, n_pairs = n_pairs,
                   median_CAF = median(caf), median_Epi = median(epi),
                   median_diff = median(caf - epi),
                   W = unname(wt$statistic), p = wt$p.value)
      })
      comp_test <- do.call(rbind, test_rows)
      comp_test$pair_unit <- UNIT_LABEL
      comp_test$p_adj_BH <- p.adjust(comp_test$p, method = "BH")
      comp_test <- comp_test[order(comp_test$p_adj_BH, na.last = TRUE), ]
      write.csv(comp_test, file.path(agg, "COHORT_compartment_test.csv"), row.names = FALSE)
      
      log_msg("Cohort compartment analysis done: %d summary rows, %d genes tested (pairs = %ss)",
              nrow(comp_summary), nrow(comp_test), UNIT_LABEL)
    }
  }
}

## ============================================================================
## 9. ENO1 SENSITIVITY CHECK   [single pass over the cohort]
##    ENO1 is glycolytic and ubiquitous; the hypoxia module is fitted with and
##    without it. The object is loaded once and both sets are scored from it.
## ============================================================================
## v6: the module in use, with and without ENO1
HYPOXIA_7 <- MODULES$Hypoxia                 # name kept so nothing below changes
HYPOXIA_6 <- setdiff(HYPOXIA_7, "ENO1")

hypoxia_logit_both <- function(rds_path, sid) {
  obj <- readRDS(rds_path)
  obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
  DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else Assays(obj)[1]
  if (!"P4HA1" %in% rownames(obj)) return(NULL)
  
  fit_variant <- function(hyp_genes, tag) {
    g <- intersect(hyp_genes, rownames(obj))
    if (length(g) < 2) return(NULL)
    nm  <- paste0("HypTmp_", tag)
    set.seed(SEED)
    o <- tryCatch(add_module_score_fallback(obj, features = list(g), name = nm,
                                            assay = DefaultAssay(obj), nbin = 20, ctrl = 50),
                  error = function(e) NULL)
    col <- paste0(nm, "1")
    if (is.null(o) || !col %in% colnames(o@meta.data)) return(NULL)
    d <- FetchData(o, vars = c("P4HA1", col), layer = "data")
    y <- as.integer(d$P4HA1 > 0)
    if (length(unique(y)) != 2) return(NULL)
    dd <- data.frame(y = y, x = d[[col]])
    f <- tryCatch(glm(y ~ x, data = dd, family = binomial()), error = function(e) NULL)
    if (is.null(f)) return(NULL)
    s <- summary(f)$coefficients
    if (nrow(s) < 2) return(NULL)
    data.frame(sample_id = sid, variant = tag, n_genes_used = length(g),
               log_OR = s[2, "Estimate"], SE = s[2, "Std. Error"])
  }
  rbind(fit_variant(HYPOXIA_7, "with_ENO1"),
        fit_variant(HYPOXIA_6, "without_ENO1"))
}

sens <- lapply(seq_len(nrow(work)), function(i) {
  path <- work$path[i]; sid <- work$sample_id[i]
  if (is.na(path) || !file.exists(path)) return(NULL)
  log_msg("[eno1 %d/%d] %s", i, nrow(work), sid)
  tryCatch(hypoxia_logit_both(path, sid),
           error = function(e) { log_msg("  eno1 skipped: %s", conditionMessage(e)); NULL })
})
sens <- dplyr::bind_rows(sens)

if (nrow(sens) == 0) {                                   # v5 (C8)
  log_msg("ENO1 check: no section produced a fit.")
  pooled_sens <- data.frame()
} else {
  pooled_sens <- sens %>%
    dplyr::group_by(variant) %>%
    dplyr::group_modify(~ meta_pool(.x$log_OR, .x$SE^2, .x$sample_id)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(OR = exp(est), lower = exp(lower), upper = exp(upper))
}
write.csv(sens,        file.path(agg, "ENO1_sensitivity_per_sample.csv"), row.names = FALSE)
write.csv(pooled_sens, file.path(agg, "ENO1_sensitivity_pooled.csv"),     row.names = FALSE)

cat("\n=== Hypoxia pooled OR: ENO1 sensitivity ===\n")
for (v in c("with_ENO1", "without_ENO1")) {
  r <- if (nrow(pooled_sens)) pooled_sens[pooled_sens$variant == v, ] else pooled_sens
  if (nrow(r) == 1)
    cat(sprintf("  %-13s OR = %.2f  (95%% CI %.2f-%.2f)  k=%d %s  I2=%.0f%%\n",
                v, r$OR, r$lower, r$upper, as.integer(r$k), K_UNIT, r$I2))
}
cat("\nIf the two ORs are close, ENO1 is not distorting the hypoxia result.\n")

## ============================================================================
## 10. STRATIFIED MULTIVARIABLE MODEL
##     At the SAME level of hypoxia, does collagen (+ iCAF) still track P4HA1?
##       (A) ADJUSTED     P(P4HA1>0) ~ Hypoxia + CollagenSynthesis + iCAF
##       (B) STRATIFIED   within within-section hypoxia tertiles (relative)
##       (C) INTERACTION  Hypoxia * CollagenSynthesis
##     Predictors here are NOT standardised, so these numbers are not
##     comparable to sections 13, 15 and 16, which are.
## ============================================================================
STRAT_MODULES   <- c("Hypoxia", "CollagenSynthesis", "iCAF")
HYPOXIA_N_BINS  <- 3
STRAT_MIN_SPOTS <- 30
BIN_LABELS      <- c("Hyp_low", "Hyp_mid", "Hyp_high")

strat_expr_one <- function(rds_path) {
  obj <- readRDS(rds_path)
  obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
  DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else Assays(obj)[1]
  if (!"P4HA1" %in% rownames(obj)) return(NULL)
  
  score_cols <- character(0)
  for (nm in STRAT_MODULES) {
    g <- intersect(MODULES[[nm]], rownames(obj))
    if (length(g) < 2) return(NULL)
    set.seed(SEED)
    obj <- tryCatch(
      add_module_score_fallback(obj, features = list(g), name = nm,
                                assay = DefaultAssay(obj), nbin = 20, ctrl = 50),
      error = function(e) NULL)
    if (is.null(obj) || !paste0(nm, "1") %in% colnames(obj@meta.data)) return(NULL)
    names(obj@meta.data)[names(obj@meta.data) == paste0(nm, "1")] <- nm
    score_cols <- c(score_cols, nm)
  }
  d <- FetchData(obj, vars = c("P4HA1", score_cols), layer = "data")
  d$y <- as.integer(d$P4HA1 > 0)
  d
}

fit_strat_one <- function(d, sid) {
  out <- list(adjusted = NULL, stratified = NULL, interaction = NULL, collin = NULL)
  if (is.null(d) || nrow(d) < STRAT_MIN_SPOTS || length(unique(d$y)) != 2) return(out)
  
  cc <- suppressWarnings(cor(d[, STRAT_MODULES], method = "spearman"))
  out$collin <- data.frame(sample_id = sid,
                           Hyp_Col  = cc["Hypoxia", "CollagenSynthesis"],
                           Hyp_iCAF = cc["Hypoxia", "iCAF"],
                           Col_iCAF = cc["CollagenSynthesis", "iCAF"])
  
  ## (A) adjusted, hypoxia continuous
  fa <- tryCatch(glm(y ~ Hypoxia + CollagenSynthesis + iCAF, data = d, family = binomial()),
                 error = function(e) NULL)
  if (!is.null(fa)) {
    s  <- summary(fa)$coefficients
    tm <- intersect(c("Hypoxia", "CollagenSynthesis", "iCAF"), rownames(s))
    if (length(tm))
      out$adjusted <- data.frame(sample_id = sid, term = tm,
                                 log_OR = s[tm, "Estimate"], SE = s[tm, "Std. Error"],
                                 p = s[tm, "Pr(>|z|)"], row.names = NULL)
  }
  
  ## (B) stratified within within-section hypoxia tertiles
  br <- tryCatch(quantile(d$Hypoxia, probs = seq(0, 1, length.out = HYPOXIA_N_BINS + 1),
                          na.rm = TRUE), error = function(e) NULL)
  if (!is.null(br) && length(unique(br)) == HYPOXIA_N_BINS + 1) {
    d$hyp_bin <- cut(d$Hypoxia, breaks = br, labels = BIN_LABELS, include.lowest = TRUE)
    rows <- lapply(BIN_LABELS, function(b) {
      db <- d[which(d$hyp_bin == b), , drop = FALSE]
      if (nrow(db) < STRAT_MIN_SPOTS || length(unique(db$y)) != 2) return(NULL)
      fb <- tryCatch(glm(y ~ CollagenSynthesis + iCAF, data = db, family = binomial()),
                     error = function(e) NULL)
      if (is.null(fb)) return(NULL)
      s  <- summary(fb)$coefficients
      tm <- intersect(c("CollagenSynthesis", "iCAF"), rownames(s))
      if (!length(tm)) return(NULL)
      data.frame(sample_id = sid, hyp_bin = b, term = tm,
                 log_OR = s[tm, "Estimate"], SE = s[tm, "Std. Error"],
                 p = s[tm, "Pr(>|z|)"], n_spots = nrow(db), row.names = NULL)
    })
    out$stratified <- do.call(rbind, rows)
  }
  
  ## (C) interaction. Exports ALL coefficients so beta_C from THIS model is
  ##     available for the crossover z_H* = -beta_C/beta_HC.
  fi <- tryCatch(glm(y ~ Hypoxia * CollagenSynthesis + iCAF, data = d, family = binomial()),
                 error = function(e) NULL)
  if (!is.null(fi)) {
    s  <- summary(fi)$coefficients
    tm <- setdiff(rownames(s), "(Intercept)")
    if (length(tm))
      out$interaction <- data.frame(sample_id = sid, term = tm,
                                    log_OR = s[tm, "Estimate"], SE = s[tm, "Std. Error"],
                                    p = s[tm, "Pr(>|z|)"], row.names = NULL)
  }
  out
}

strat_acc <- list(adjusted = list(), stratified = list(), interaction = list(), collin = list())
for (i in seq_len(nrow(work))) {
  path <- work$path[i]; sid <- work$sample_id[i]
  if (is.na(path) || !file.exists(path)) next
  log_msg("[strat %d/%d] %s", i, nrow(work), sid)
  d <- tryCatch(strat_expr_one(path), error = function(e) { log_msg("  strat skipped: %s", conditionMessage(e)); NULL })
  r <- fit_strat_one(d, sid)
  for (k in names(strat_acc)) if (!is.null(r[[k]])) strat_acc[[k]][[length(strat_acc[[k]]) + 1]] <- r[[k]]
}

adj_all <- bind(strat_acc$adjusted)
str_all <- bind(strat_acc$stratified)
int_all <- bind(strat_acc$interaction)
col_all <- bind(strat_acc$collin)

write.csv(adj_all, file.path(agg, "STRAT_adjusted_per_sample.csv"),     row.names = FALSE)
write.csv(str_all, file.path(agg, "STRAT_stratified_per_sample.csv"),   row.names = FALSE)
write.csv(int_all, file.path(agg, "STRAT_interaction_per_sample.csv"),  row.names = FALSE)
write.csv(col_all, file.path(agg, "STRAT_predictor_collinearity.csv"),  row.names = FALSE)

if (nrow(col_all) > 0) {
  cat(sprintf("\n=== Predictor collinearity (median Spearman across %ss) ===\n", UNIT_LABEL))
  cat(sprintf("  Hypoxia vs Collagen : %.2f\n", unit_median(col_all$Hyp_Col,  col_all$sample_id)))
  cat(sprintf("  Hypoxia vs iCAF     : %.2f\n", unit_median(col_all$Hyp_iCAF, col_all$sample_id)))
  cat(sprintf("  Collagen vs iCAF    : %.2f\n", unit_median(col_all$Col_iCAF, col_all$sample_id)))
  cat("  (|rho| approaching ~0.7+ means the adjusted ORs are hard to separate.)\n")
}

if (nrow(adj_all) > 0) {
  adj_pooled <- adj_all %>% group_by(term) %>%
    group_modify(~ meta_pool(.x$log_OR, .x$SE^2, .x$sample_id)) %>% ungroup() %>%
    mutate(OR = exp(est), lower = exp(lower), upper = exp(upper)) %>%
    arrange(desc(OR))
  write.csv(adj_pooled, file.path(agg, "STRAT_adjusted_pooled.csv"), row.names = FALSE)
}

if (nrow(str_all) > 0) {
  str_pooled <- str_all %>% group_by(term, hyp_bin) %>%
    group_modify(~ meta_pool(.x$log_OR, .x$SE^2, .x$sample_id)) %>% ungroup() %>%
    mutate(OR = exp(est), lower = exp(lower), upper = exp(upper),
           hyp_bin = factor(hyp_bin, levels = BIN_LABELS))
  write.csv(str_pooled, file.path(agg, "STRAT_stratified_pooled.csv"), row.names = FALSE)
}

if (nrow(int_all) > 0) {
  int_pooled <- int_all %>% group_by(term) %>%
    group_modify(~ meta_pool(.x$log_OR, .x$SE^2, .x$sample_id)) %>% ungroup() %>%
    mutate(OR = exp(est), lower = exp(lower), upper = exp(upper))
  write.csv(int_pooled, file.path(agg, "STRAT_interaction_pooled.csv"), row.names = FALSE)
  
  it_row <- int_pooled[grepl(":", int_pooled$term), ]
  cat("\n=== Hypoxia x Collagen interaction (pooled) ===\n")
  if (nrow(it_row) == 1) {
    cat(sprintf("  interaction coefficient = %+.4f  (OR scale %.2f, 95%% CI %.2f-%.2f)  k=%d %s  I2=%.0f%%\n",
                it_row$est, it_row$OR, it_row$lower, it_row$upper, as.integer(it_row$k), K_UNIT, it_row$I2))
    if (it_row$lower > 1 || it_row$upper < 1)
      cat(sprintf("  CI EXCLUDES 1 -> the collagen effect DOES change with hypoxia (%s).\n",
                  if (it_row$est > 0) "stronger as hypoxia rises" else "weaker as hypoxia rises"))
    else
      cat("  CI includes 1 -> no evidence the collagen effect changes with hypoxia.\n")
  }
  
  bC  <- int_pooled$est[int_pooled$term == "CollagenSynthesis"]
  bHC <- it_row$est
  if (length(bC) == 1 && length(bHC) == 1 && is.finite(bC) && is.finite(bHC) && bHC != 0) {
    cat(sprintf("\n  beta_C (interaction model) = %+.4f ; beta_HC = %+.4f\n", bC, bHC))
    cat("  Collagen slope = beta_C + beta_HC * z_H\n")
    cat(sprintf("  Crossover z_H* = %+.3f  (ratio of two pooled estimates; no CI)\n", -bC / bHC))
    cat("  Check against HYPOXIA_SCORE_RANGE.csv: if z_H* lies outside the observed\n")
    cat("  range, the crossover is an extrapolation and should not be reported.\n")
  }
}

log_msg("Stratified multivariable step done.")

## ============================================================================
## 11. PROTEOSTASIS MODULES  (exploratory)
##
##     HYPOTHESIS (not established): P4H may also retain the molecule in the
##     ER for quality control, so P4HA1 may track proteostasis rather than
##     collagen.
##
##     MODULES (the audit re-runs at run time and can halt the section):
##       ERstress       27 core ER-stress/UPR genes hand-selected from
##                      HALLMARK_UNFOLDED_PROTEIN_RESPONSE (MSigDB M5922).
##                      VEGFA removed: it is also in the Hypoxia module.
##       ERstressStrict the above minus HYOU1, DDIT4, ERO1A (hypoxia-
##                      responsive). Every model is fitted with BOTH; the
##                      strict one is primary. Section 14 reports whether the
##                      objects carry ERO1A or the legacy ERO1L.
##       ERAD           GOBP_ERAD_PATHWAY (MSigDB M15534), with the genes it
##                      shares with ERstress stripped.
##
##     OPEN: HYOU1, DDIT4, ERO1A are removed as hypoxia-responsive but ATF4
##     and CHAC1 (PERK-responsive) are kept. Section 23 reports whether that
##     choice moves the result.
##
##     NOT INCLUDED: the ER chaperone complex (GO:0034663) is defined to
##     include P4HB, the beta subunit of P4HA1's own enzyme.
##
##     Gene threshold: max(2, 30% of the set), so a 6-gene and a 104-gene set
##     are held to the same proportion.
##
##     Figures live in section 12 and use the STRICT variant only.
## ============================================================================
PROT_DISJOINT   <- TRUE   # strip the genes shared with ERstress out of ERAD
PROT_MIN_SPOTS  <- 30
PROT_MIN_FRAC   <- 0.30
PROT_MIN_ABS    <- 2
PROT_N_BINS     <- 3
PROT_BIN_LABELS <- c("Hyp_low", "Hyp_mid", "Hyp_high")

min_genes_for <- function(gene_set) max(PROT_MIN_ABS, ceiling(PROT_MIN_FRAC * length(gene_set)))

ERSTRESS_GENES <- c("ATF3","ATF4","ATF6","ERN1","EIF2AK3","XBP1","HSPA5","HSP90B1",
                    "CALR","HERPUD1","EDEM1","DNAJB9","DNAJC3","DNAJA4","PDIA5",
                    "PDIA6","SERP1","WFS1","HYOU1","KDELR3","SSR1","ERO1A","PREB",
                    "BAG3","HSPA9","CHAC1","DDIT4")

HYPOXIA_RESPONSIVE_IN_ERSTRESS <- c("HYOU1",   # hypoxia up-regulated 1
                                    "DDIT4",   # HIF-1 target (REDD1)
                                    "ERO1A")   # hypoxia-inducible
ERSTRESS_STRICT_GENES <- setdiff(ERSTRESS_GENES, HYPOXIA_RESPONSIVE_IN_ERSTRESS)

ERAD_GENES_FULL <- c("AFG2B","AMFR","ANKZF1","AQP11","ATF6","ATXN3","ATXN3L","AUP1",
                     "BAG6","BCAP31","BRSK2","CALR","CALR3","CANX","CAV1","CCDC47","CLGN","DERL1",
                     "DERL2","DERL3","DNAJB12","DNAJB2","DNAJB9","DNAJC10","ECPAS","EDEM1","EDEM2",
                     "EDEM3","ERLEC1","ERLIN1","ERLIN2","FAF1","FAF2","FAM8A1","FBXO17","FBXO2",
                     "FBXO27","FBXO44","FBXO6","FOXRED2","GET4","HERPUD1","HM13","HSP90B1","HSPA5",
                     "JKAMP","MAN1A1","MAN1A2","MAN1B1","MAN1C1","MARCHF6","NCCRP1","NFE2L2","NPLOC4",
                     "OS9","PRKN","PSMC6","RCN3","RHBDD1","RHBDD2","RNF103","RNF121","RNF139","RNF145",
                     "RNF175","RNF185","RNF5","RNFT1","RNFT2","SDF2L1","SEC61B","SEL1L","SEL1L2",
                     "SELENOS","SGTA","SIK2","STT3B","STUB1","SVIP","SYVN1","TMEM129","TMEM259",
                     "TMEM67","TMUB1","TMUB2","TMX1","TOR1A","TRIM13","TRIM25","UBAC2","UBE2G2",
                     "UBE2J1","UBE2J2","UBE4A","UBE4B","UBQLN1","UBQLN2","UBXN1","UBXN10","UBXN2A",
                     "UBXN4","UBXN6","UBXN8","UFD1","UMOD","USP13","USP14","USP19","USP25","VCP",
                     "WFS1","XBP1","YOD1")

## strip shared genes using the FULL ERstress list so ERAD is identical
## regardless of which ERstress variant is being fitted
ERAD_GENES <- if (PROT_DISJOINT) setdiff(ERAD_GENES_FULL, ERSTRESS_GENES) else ERAD_GENES_FULL

MODULES$ERstress       <- ERSTRESS_GENES
MODULES$ERstressStrict <- ERSTRESS_STRICT_GENES
MODULES$ERAD           <- ERAD_GENES

ERSTRESS_VARIANTS <- c("ERstress", "ERstressStrict")
BASE_MODULES      <- c("Hypoxia", "CollagenSynthesis", "iCAF")
SCORE_MODULES     <- c(BASE_MODULES, ERSTRESS_VARIANTS, "ERAD")

## ---------------------------------------------------------------------------
## 11.0  CIRCULARITY / OVERLAP AUDIT  -- halts the section on failure
## ---------------------------------------------------------------------------
DANGER_GENES <- c("P4HA1","P4HA2","P4HA3","P4HB","PLOD1","PLOD2","PLOD3")
NEW_MODULES  <- c(ERSTRESS_VARIANTS, "ERAD")

cat("\n", strrep("=", 72), "\n", sep = "")
cat("SECTION 11 GENE-SET AUDIT\n")
cat(strrep("=", 72), "\n", sep = "")

audit_fail <- FALSE
for (nm in NEW_MODULES) {
  g   <- MODULES[[nm]]
  bad <- intersect(DANGER_GENES, g)
  cat(sprintf("\n%s  (n = %d genes, min detected required = %d)\n",
              nm, length(g), min_genes_for(g)))
  cat(sprintf("  circularity (P4HA1 + enzyme partners): %s\n",
              if (length(bad)) paste0("*** FAIL: ", paste(bad, collapse=", "), " ***") else "clean"))
  if (length(bad)) audit_fail <- TRUE
  for (other in setdiff(names(MODULES), c(nm, ERSTRESS_VARIANTS))) {
    ov <- intersect(g, MODULES[[other]])
    if (length(ov))
      cat(sprintf("  overlap with %-20s %s\n", other, paste(ov, collapse=", ")))
  }
}
cat(sprintf("\nERstress vs ERAD shared genes: %d  (PROT_DISJOINT = %s)\n",
            length(intersect(MODULES$ERstress, MODULES$ERAD)), PROT_DISJOINT))
cat(sprintf("Sensitivity variant removes %d hypoxia-responsive genes: %s\n",
            length(HYPOXIA_RESPONSIVE_IN_ERSTRESS),
            paste(HYPOXIA_RESPONSIVE_IN_ERSTRESS, collapse = ", ")))
cat("\nMinimum detected genes required per module:\n")
for (nm in SCORE_MODULES)
  cat(sprintf("  %-18s %3d genes -> need %d\n", nm, length(MODULES[[nm]]), min_genes_for(MODULES[[nm]])))
if (audit_fail)
  stop("Section 11 halted: a proteostasis module contains P4HA1 or an enzyme partner.")
cat("\nAudit passed.\n\n")

## ---------------------------------------------------------------------------
## 11.1  score all modules per spot
## ---------------------------------------------------------------------------
prot_expr_one <- function(rds_path, sid) {
  obj <- readRDS(rds_path)
  obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
  DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else Assays(obj)[1]
  if (!"P4HA1" %in% rownames(obj)) return(NULL)
  
  cov <- list()
  for (nm in SCORE_MODULES) {
    full <- MODULES[[nm]]
    g    <- intersect(full, rownames(obj))
    need <- min_genes_for(full)
    cov[[nm]] <- c(n_found = length(g), n_total = length(full), n_required = need)
    if (length(g) < need) {
      log_msg("  %s: only %d/%d genes of %s detected (need %d) -> section skipped",
              sid, length(g), length(full), nm, need)
      return(NULL)
    }
    set.seed(SEED)
    obj <- tryCatch(
      add_module_score_fallback(obj, features = list(g), name = nm,
                                assay = DefaultAssay(obj), nbin = 20, ctrl = 50),
      error = function(e) NULL)
    if (is.null(obj) || !paste0(nm, "1") %in% colnames(obj@meta.data)) return(NULL)
    names(obj@meta.data)[names(obj@meta.data) == paste0(nm, "1")] <- nm
  }
  d <- FetchData(obj, vars = c("P4HA1", SCORE_MODULES), layer = "data")
  d$y <- as.integer(d$P4HA1 > 0)
  attr(d, "coverage") <- do.call(rbind, lapply(names(cov), function(k)
    data.frame(sample_id = sid, module = k, n_found = unname(cov[[k]]["n_found"]),
               n_total = unname(cov[[k]]["n_total"]), n_required = unname(cov[[k]]["n_required"]))))
  d
}

## variance inflation factors, no car dependency
vif_calc <- function(d, vars) {
  do.call(rbind, lapply(vars, function(v) {
    others <- setdiff(vars, v)
    f  <- as.formula(paste(v, "~", paste(others, collapse = " + ")))
    r2 <- tryCatch(summary(lm(f, data = d))$r.squared, error = function(e) NA_real_)
    data.frame(term = v, R2 = r2, VIF = if (is.finite(r2) && r2 < 1) 1/(1-r2) else NA_real_)
  }))
}

coefs_of <- function(fit, keep = NULL) {
  if (is.null(fit)) return(NULL)
  s  <- summary(fit)$coefficients
  tm <- setdiff(rownames(s), "(Intercept)")
  if (!is.null(keep)) tm <- intersect(keep, tm)
  if (!length(tm)) return(NULL)
  data.frame(term = tm, log_OR = s[tm,"Estimate"], SE = s[tm,"Std. Error"],
             p = s[tm,"Pr(>|z|)"], row.names = NULL)
}

## ---------------------------------------------------------------------------
## 11.2  fit every model for one section, for BOTH ERstress variants
## ---------------------------------------------------------------------------
fit_prot_one <- function(d, sid) {
  out <- list(marginal=NULL, adjusted=NULL, stratified=NULL,
              interaction=NULL, vif=NULL, collin=NULL, coverage=NULL)
  if (is.null(d) || nrow(d) < PROT_MIN_SPOTS || length(unique(d$y)) != 2) return(out)
  out$coverage <- attr(d, "coverage")
  
  var_cor <- suppressWarnings(cor(d$ERstress, d$ERstressStrict, method = "spearman"))
  
  marg <- list(); adj <- list(); strat <- list(); int <- list()
  vifs <- list(); cols <- list()
  
  for (ES in ERSTRESS_VARIANTS) {
    preds <- c(BASE_MODULES, ES, "ERAD")
    
    cc <- suppressWarnings(cor(d[, preds], method = "spearman"))
    cols[[ES]] <- data.frame(sample_id = sid, variant = ES,
                             Hyp_ERstress  = cc["Hypoxia", ES],
                             Hyp_ERAD      = cc["Hypoxia","ERAD"],
                             Col_ERstress  = cc["CollagenSynthesis", ES],
                             Col_ERAD      = cc["CollagenSynthesis","ERAD"],
                             ERstress_ERAD = cc[ES,"ERAD"],
                             ERstress_variants_cor = var_cor)
    v <- vif_calc(d, preds); v$sample_id <- sid; v$variant <- ES; vifs[[ES]] <- v
    
    ## (A) MARGINAL, module alone (ERAD fitted under both variants; identical)
    mods_alone <- c(ES, "ERAD")
    marg[[ES]] <- do.call(rbind, lapply(mods_alone, function(m) {
      f <- tryCatch(glm(as.formula(paste("y ~", m)), data = d, family = binomial()),
                    error = function(e) NULL)
      r <- coefs_of(f, keep = m); if (is.null(r)) return(NULL)
      cbind(sample_id = sid, variant = ES, r)
    }))
    
    ## (B) ADJUSTED
    fa <- tryCatch(glm(as.formula(paste("y ~", paste(preds, collapse = " + "))),
                       data = d, family = binomial()), error = function(e) NULL)
    r <- coefs_of(fa, keep = preds)
    if (!is.null(r)) adj[[ES]] <- cbind(sample_id = sid, variant = ES, r)
    
    ## (C) STRATIFIED by within-section hypoxia tertiles
    br <- tryCatch(quantile(d$Hypoxia, probs = seq(0,1,length.out = PROT_N_BINS+1), na.rm = TRUE),
                   error = function(e) NULL)
    if (!is.null(br) && length(unique(br)) == PROT_N_BINS + 1) {
      d$hyp_bin <- cut(d$Hypoxia, breaks = br, labels = PROT_BIN_LABELS, include.lowest = TRUE)
      instrat <- c("CollagenSynthesis","iCAF", ES, "ERAD")
      strat[[ES]] <- do.call(rbind, lapply(PROT_BIN_LABELS, function(b) {
        db <- d[which(d$hyp_bin == b), , drop = FALSE]
        if (nrow(db) < PROT_MIN_SPOTS || length(unique(db$y)) != 2) return(NULL)
        fb <- tryCatch(glm(as.formula(paste("y ~", paste(instrat, collapse = " + "))),
                           data = db, family = binomial()), error = function(e) NULL)
        r <- coefs_of(fb, keep = instrat); if (is.null(r)) return(NULL)
        cbind(sample_id = sid, variant = ES, hyp_bin = b, r, n_spots = nrow(db))
      }))
    }
    
    ## (D) INTERACTION with hypoxia
    fi <- tryCatch(glm(as.formula(paste0("y ~ Hypoxia * ", ES,
                                         " + CollagenSynthesis + iCAF + ERAD")),
                       data = d, family = binomial()), error = function(e) NULL)
    r <- coefs_of(fi)
    if (!is.null(r)) {
      r$term[grepl(":", r$term)] <- "Hypoxia:ERstress"
      r$term[r$term == ES] <- "ERstress"
      int[[ES]] <- cbind(sample_id = sid, variant = ES, r)
    }
  }
  
  out$marginal   <- bind(marg);  out$adjusted    <- bind(adj)
  out$stratified <- bind(strat); out$interaction <- bind(int)
  out$vif        <- bind(vifs);  out$collin      <- bind(cols)
  for (nm in c("marginal","adjusted","stratified")) {
    if (!is.null(out[[nm]]) && nrow(out[[nm]]))
      out[[nm]]$term[out[[nm]]$term %in% ERSTRESS_VARIANTS] <- "ERstress"
  }
  if (!is.null(out$vif) && nrow(out$vif))
    out$vif$term[out$vif$term %in% ERSTRESS_VARIANTS] <- "ERstress"
  out
}

## ---------------------------------------------------------------------------
## 11.3  run across the cohort
## ---------------------------------------------------------------------------
prot_acc <- list(marginal=list(), adjusted=list(), stratified=list(),
                 interaction=list(), vif=list(), collin=list(), coverage=list())
for (i in seq_len(nrow(work))) {
  path <- work$path[i]; sid <- work$sample_id[i]
  if (is.na(path) || !file.exists(path)) next
  log_msg("[prot %d/%d] %s", i, nrow(work), sid)
  d <- tryCatch(prot_expr_one(path, sid),
                error = function(e) { log_msg("  prot skipped: %s", conditionMessage(e)); NULL })
  r <- fit_prot_one(d, sid)
  for (k in names(prot_acc))
    if (!is.null(r[[k]]) && nrow(r[[k]])) prot_acc[[k]][[length(prot_acc[[k]])+1]] <- r[[k]]
}

marg_all <- bind(prot_acc$marginal);    padj_all <- bind(prot_acc$adjusted)
pstr_all <- bind(prot_acc$stratified);  pint_all <- bind(prot_acc$interaction)
vif_all  <- bind(prot_acc$vif);         pcol_all <- bind(prot_acc$collin)
cov_all  <- bind(prot_acc$coverage)

write.csv(marg_all, file.path(agg,"PROT_marginal_per_sample.csv"),    row.names=FALSE)
write.csv(padj_all, file.path(agg,"PROT_adjusted_per_sample.csv"),    row.names=FALSE)
write.csv(pstr_all, file.path(agg,"PROT_stratified_per_sample.csv"),  row.names=FALSE)
write.csv(pint_all, file.path(agg,"PROT_interaction_per_sample.csv"), row.names=FALSE)
write.csv(vif_all,  file.path(agg,"PROT_VIF_per_sample.csv"),         row.names=FALSE)
write.csv(pcol_all, file.path(agg,"PROT_collinearity_per_sample.csv"),row.names=FALSE)
write.csv(cov_all,  file.path(agg,"PROT_gene_coverage.csv"),          row.names=FALSE)

n_in <- if (nrow(padj_all)) dplyr::n_distinct(padj_all$sample_id) else 0
cat(sprintf("\n=== Sections contributing: %d of %d (%d of %d patients) ===\n", n_in, nrow(work),
            if (n_in) length(unique(patient_of(unique(padj_all$sample_id)))) else 0L,
            length(unique(work$patient_id))))
if (n_in < nrow(work)) {
  missing <- setdiff(work$sample_id, unique(padj_all$sample_id))
  cat("  excluded: ", paste(missing, collapse = ", "), "\n")
  cat("  (reason is in the log above, per section)\n")
}

if (nrow(cov_all) > 0) {
  cat("\n=== Gene coverage per module (across contributing sections) ===\n")
  cs <- cov_all %>% group_by(module) %>%
    summarise(median_found = median(n_found), n_total = max(n_total),
              min_found = min(n_found), required = max(n_required), .groups = "drop")
  ## v5 (C5): median of an even number of sections can be x.5, so %d crashed
  for (i in seq_len(nrow(cs)))
    cat(sprintf("  %-18s %5.1f / %3.0f detected (min %3.0f, required %3.0f)\n", cs$module[i],
                cs$median_found[i], cs$n_total[i], cs$min_found[i], cs$required[i]))
}

if (nrow(pcol_all) > 0) {
  cat(sprintf("\n=== Predictor collinearity (median Spearman across %ss) ===\n", UNIT_LABEL))
  for (ES in ERSTRESS_VARIANTS) {
    p <- pcol_all[pcol_all$variant == ES, ]
    if (!nrow(p)) next
    cat(sprintf("  [%s]\n", ES))
    cat(sprintf("     Hypoxia vs ERstress : %+.2f\n", unit_median(p$Hyp_ERstress,  p$sample_id)))
    cat(sprintf("     Hypoxia vs ERAD     : %+.2f\n", unit_median(p$Hyp_ERAD,      p$sample_id)))
    cat(sprintf("     Collagen vs ERstress: %+.2f\n", unit_median(p$Col_ERstress,  p$sample_id)))
    cat(sprintf("     ERstress vs ERAD    : %+.2f\n", unit_median(p$ERstress_ERAD, p$sample_id)))
  }
  p1 <- pcol_all[pcol_all$variant == ERSTRESS_VARIANTS[1], ]
  cat(sprintf("  Correlation between the two ERstress variants: %+.3f (median)\n",
              unit_median(p1$ERstress_variants_cor, p1$sample_id)))
}
if (nrow(vif_all) > 0) {
  cat("\n=== Variance inflation factors (median across sections) ===\n")
  vs <- vif_all %>% group_by(variant, term) %>%
    summarise(median_VIF = median(VIF, na.rm=TRUE),
              max_VIF = suppressWarnings(max(VIF, na.rm=TRUE)), .groups="drop") %>%
    arrange(variant, desc(median_VIF))
  for (i in seq_len(nrow(vs)))
    cat(sprintf("  %-16s %-18s VIF %5.2f (max %6.2f)%s\n", vs$variant[i], vs$term[i],
                vs$median_VIF[i], vs$max_VIF[i],
                if (is.finite(vs$median_VIF[i]) && vs$median_VIF[i] > 5) "  <-- HIGH" else ""))
  cat("  VIF > 5: hard to separate.  VIF > 10: do not interpret alone.\n")
}

## ---------------------------------------------------------------------------
## 11.4  pool and report  (figures live in section 12)
## ---------------------------------------------------------------------------
pool_and_write <- function(df, group_vars, csv_name) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  out <- df %>% group_by(across(all_of(group_vars))) %>%
    group_modify(~ meta_pool(.x$log_OR, .x$SE^2, .x$sample_id)) %>% ungroup() %>%
    mutate(OR = exp(est), lower = exp(lower), upper = exp(upper))
  write.csv(out, file.path(agg, csv_name), row.names = FALSE)
  out
}

marg_pooled <- pool_and_write(marg_all, c("variant","term"),           "PROT_marginal_pooled.csv")
padj_pooled <- pool_and_write(padj_all, c("variant","term"),           "PROT_adjusted_pooled.csv")
pstr_pooled <- pool_and_write(pstr_all, c("variant","term","hyp_bin"), "PROT_stratified_pooled.csv")
pint_pooled <- pool_and_write(pint_all, c("variant","term"),           "PROT_interaction_pooled.csv")

excl1 <- function(lo, hi) isTRUE((lo > 1) | (hi < 1))
star  <- function(lo, hi) if (excl1(lo, hi)) "  *" else "   "

if (!is.null(marg_pooled)) {
  cat("\n=== MARGINAL (module alone, no adjustment) ===\n")
  for (i in seq_len(nrow(marg_pooled))) with(marg_pooled[i,],
                                             cat(sprintf("  %-16s %-10s OR %6.2f (95%% CI %5.2f-%7.2f) k=%d I2=%2.0f%% tau=%.2f%s\n",
                                                         variant, term, OR, lower, upper, as.integer(k), I2, sqrt(tau2), star(lower, upper))))
  cat(sprintf("  k counts %s. ERAD is identical under both variants by construction.\n", K_UNIT))
}

if (!is.null(padj_pooled)) {
  cat("\n=== ADJUSTED (all five predictors) ===\n")
  for (i in seq_len(nrow(padj_pooled))) with(padj_pooled[i,],
                                             cat(sprintf("  %-16s %-18s OR %6.2f (95%% CI %5.2f-%7.2f) k=%d I2=%2.0f%% tau=%.2f%s\n",
                                                         variant, term, OR, lower, upper, as.integer(k), I2, sqrt(tau2), star(lower, upper))))
  cat("  * = 95% CI excludes 1\n")
  
  a <- padj_pooled[padj_pooled$term == "ERstress" & padj_pooled$variant == "ERstress", ]
  b <- padj_pooled[padj_pooled$term == "ERstress" & padj_pooled$variant == "ERstressStrict", ]
  if (nrow(a) == 1 && nrow(b) == 1) {
    cat("\n  --- ERstress sensitivity (removing HYOU1, DDIT4, ERO1A) ---\n")
    cat(sprintf("     full 27 genes  : OR %.2f (%.2f-%.2f)\n", a$OR, a$lower, a$upper))
    cat(sprintf("     strict 24 genes: OR %.2f (%.2f-%.2f)\n", b$OR, b$lower, b$upper))
    cat(sprintf("     log-scale shift: %+.4f   ratio of ORs: %.3f   effect reduced %.0f%%\n",
                b$est - a$est, b$OR / a$OR, 100*(1 - b$OR/a$OR)))
    cat(sprintf("     VERDICT: %s\n",
                if (excl1(b$lower, b$upper) && excl1(a$lower, a$upper))
                  "both variants exclude 1 - not driven by the hypoxia-responsive genes"
                else if (excl1(a$lower, a$upper) && !excl1(b$lower, b$upper))
                  "*** DOES NOT SURVIVE - significant only with the hypoxia-responsive genes ***"
                else if (!excl1(a$lower, a$upper) && excl1(b$lower, b$upper))
                  "significant only in the strict variant"
                else "neither variant is significant"))
    cat("     Report the STRICT estimate, and the size of the drop as well as the verdict.\n")
  }
}

if (!is.null(pstr_pooled)) {
  cat("\n=== STRATIFIED by hypoxia tertile ===\n")
  ps <- pstr_pooled[order(pstr_pooled$variant, pstr_pooled$term,
                          match(pstr_pooled$hyp_bin, PROT_BIN_LABELS)), ]
  for (i in seq_len(nrow(ps))) with(ps[i,],
                                    cat(sprintf("  %-16s %-18s %-9s OR %6.2f (95%% CI %5.2f-%7.2f) k=%d%s\n",
                                                variant, term, hyp_bin, OR, lower, upper, as.integer(k), star(lower, upper))))
}

if (!is.null(pint_pooled)) {
  cat("\n=== Hypoxia x ERstress interaction ===\n")
  for (ES in ERSTRESS_VARIANTS) {
    sub <- pint_pooled[pint_pooled$variant == ES, ]
    it  <- sub[sub$term == "Hypoxia:ERstress", ]
    bE  <- sub$est[sub$term == "ERstress"]
    if (nrow(it) != 1) next
    cat(sprintf("\n  [%s]\n", ES))
    cat(sprintf("     interaction coefficient = %+.4f  (OR scale 95%% CI %.2f-%.2f)\n",
                it$est, it$lower, it$upper))
    if (excl1(it$lower, it$upper)) {
      cat(sprintf("     CI EXCLUDES 1 -> the ER-stress effect DOES change with hypoxia (%s).\n",
                  if (it$est > 0) "stronger as hypoxia rises" else "weaker as hypoxia rises"))
    } else {
      cat("     CI includes 1 -> no evidence the ER-stress effect changes with hypoxia.\n")
    }
    cat("     Do NOT quote the raw per-unit interaction OR: module scores span far\n")
    cat("     less than 1 unit. Use the slopes below.\n")
    if (length(bE) == 1 && is.finite(bE) && is.finite(it$est) && it$est != 0) {
      cat(sprintf("     ERstress log-odds slope = %+.4f %+.4f * z_Hypoxia\n", bE, it$est))
      cat(sprintf("     crossover at z_Hypoxia = %+.3f (ratio of two pooled estimates; no CI)\n",
                  -bE / it$est))
      cat("     implied ERstress OR at selected hypoxia scores:\n")
      for (zh in c(-0.3, -0.1, 0, 0.1, 0.3))
        cat(sprintf("        z_Hypoxia %+.1f -> OR %8.2f\n", zh, exp(bE + it$est * zh)))
    }
  }
}

cat("\n", strrep("=", 72), "\n", sep = "")
cat("SECTION 11 COMPLETE - READ IN THIS ORDER:\n")
cat("  1. sections contributing - is anyone excluded, and why?\n")
cat("  2. gene coverage         - is each score built on enough genes?\n")
cat("  3. collinearity + VIF    - are the predictors separable?\n")
cat("  4. marginal              - any signal before conditioning?\n")
cat("  5. adjusted + SENSITIVITY VERDICT\n")
cat("  6. stratified + interaction\n")
cat("EXPLORATORY pending sign-off on the gene sets.\n")
cat("Associational, mRNA-level. Cannot demonstrate a clamp.\n")
cat(strrep("=", 72), "\n", sep = "")

log_msg("Proteostasis step done.")

################################################################################
## 12. FIGURE EXPORT  --  journal style AND presentation style
##     Reads the pooled CSVs written above, so nothing is re-analysed. Every
##     panel is defined once and rendered twice:
##       journal : 8 pt, 89 / 183 mm columns, vector PDF + 600 dpi PNG
##       slide   : 16 pt, 13.33 x 7.5 in (16:9), one panel per slide
##     Filled point = CI excludes 1, open point = CI includes 1.
################################################################################

## ============================================================================
## 12.0  OBSERVED HYPOXIA SCORE RANGE   (cached; recomputed if stale)
##       Figure 4 is drawn only over the observed score range. Every section
##       of the work list gets a row (NA quantiles if it could not be scored),
##       so the cache can be checked against the work list.
## ============================================================================
HYPRANGE_CSV <- file.path(agg, "HYPOXIA_SCORE_RANGE.csv")

hr_fresh <- FALSE
if (file.exists(HYPRANGE_CSV)) {
  hr_old <- read.csv(HYPRANGE_CSV, stringsAsFactors = FALSE)
  ## v5 (C8): a cache built on another set of sections is not reused
  hr_fresh <- setequal(hr_old$sample_id, work$sample_id)
  if (!hr_fresh) log_msg("HYPOXIA_SCORE_RANGE.csv was made from other sections - recomputing.")
}

if (!hr_fresh) {
  log_msg("Computing HYPOXIA_SCORE_RANGE.csv (one extra pass over the cohort).")
  na_row <- function(sid, n = NA_integer_)
    data.frame(sample_id = sid, n_spots = n, q01 = NA_real_, q05 = NA_real_, q25 = NA_real_,
               q50 = NA_real_, q75 = NA_real_, q95 = NA_real_, q99 = NA_real_, row.names = NULL)
  hr <- do.call(rbind, lapply(seq_len(nrow(work)), function(i) {
    p <- work$path[i]; sid <- work$sample_id[i]
    if (is.na(p) || !file.exists(p)) return(na_row(sid))
    obj <- readRDS(p)
    obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
    DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else Assays(obj)[1]
    g <- intersect(MODULES$Hypoxia, rownames(obj))
    if (length(g) < 2) return(na_row(sid, ncol(obj)))
    set.seed(SEED)
    obj2 <- tryCatch(add_module_score_fallback(obj, features = list(g), name = "H",
                                               assay = DefaultAssay(obj), nbin = 20, ctrl = 50),
                     error = function(e) { log_msg("  %s: hypoxia range skipped: %s", sid, conditionMessage(e)); NULL })
    if (is.null(obj2)) return(na_row(sid, ncol(obj)))
    v <- obj2@meta.data$H1
    qs <- function(q) unname(quantile(v, q, na.rm = TRUE))
    data.frame(sample_id = sid, n_spots = length(v),
               q01 = qs(.01), q05 = qs(.05), q25 = qs(.25), q50 = qs(.50),
               q75 = qs(.75), q95 = qs(.95), q99 = qs(.99), row.names = NULL)
  }))
  write.csv(hr, HYPRANGE_CSV, row.names = FALSE)
  log_msg("Hypoxia range written: median 5th pct %.3f, median 95th pct %.3f",
          median(hr$q05, na.rm = TRUE), median(hr$q95, na.rm = TRUE))
} else {
  log_msg("Hypoxia score range already cached for these sections - skipping recomputation.")
}

## ============================================================================
## 12.1  STYLE DEFINITIONS
## ============================================================================
PRIMARY_VARIANT <- "ERstressStrict"
BIN_NICE        <- c("Low", "Mid", "High")

MM <- function(mm) mm / 25.4

CLR <- list(pos  = "#0F4C81",
            neg  = "#B2182B",
            grey = "#8C8C8C",
            full = "#9E9E9E")

## capabilities("cairo") can be TRUE while the cairo library fails to load
## when a device is opened. Open a throwaway device once and fall back to
## grDevices::pdf if it does not work.
PDF_DEV <- local({
  tf <- tempfile(fileext = ".pdf")
  ok <- isTRUE(capabilities("cairo")) &&
    !inherits(try({ grDevices::cairo_pdf(tf); grDevices::dev.off() }, silent = TRUE), "try-error")
  unlink(tf)
  if (ok) grDevices::cairo_pdf else grDevices::pdf
})
log_msg("PDF device in use: %s", if (identical(PDF_DEV, grDevices::cairo_pdf)) "cairo_pdf" else "pdf (cairo unavailable)")

FIG_STYLES <- list(
  journal = list(
    dir = "figures_journal", base = 8, pt = 1.9, stroke = 0.45, lwd = 0.45, lab = 2.1,
    grid = 0.25, border = 0.40, key = 3, mar = 3, tags = TRUE, titles = FALSE,
    split = FALSE, devices = c("pdf", "png"), dpi = 600, w1 = MM(89), w2 = MM(183)),
  slide = list(
    dir = "figures_slides", base = 16, pt = 4.4, stroke = 0.9, lwd = 1.0, lab = 4.8,
    grid = 0.5, border = 0.7, key = 6, mar = 10, tags = FALSE, titles = TRUE,
    split = TRUE, devices = c("png", "pdf"), dpi = 200, w1 = 13.333, w2 = 13.333, h_full = 7.5)
)

theme_fig <- function(S) {
  b <- S$base
  theme_bw(base_size = b, base_family = "sans") +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(linewidth = S$grid, colour = "grey90"),
      panel.border       = element_rect(linewidth = S$border, colour = "black", fill = NA),
      axis.text          = element_text(colour = "black", size = b),
      axis.title         = element_text(colour = "black", size = b),
      axis.ticks         = element_line(linewidth = S$border, colour = "black"),
      strip.background   = element_blank(),
      strip.text         = element_text(face = "bold", size = b, hjust = 0),
      plot.title         = element_text(face = "bold", size = b + 2, hjust = 0),
      plot.subtitle      = element_text(size = b - 1, hjust = 0, colour = "grey25"),
      plot.caption       = element_text(size = b - 2, hjust = 0, colour = "grey35"),
      plot.tag           = element_text(face = "bold", size = b + 3),
      plot.tag.position  = c(0.005, 0.99),
      legend.key.size    = unit(S$key, "mm"),
      legend.text        = element_text(size = b - 1),
      legend.title       = element_text(size = b - 1),
      legend.margin      = margin(0, 0, 0, 0),
      plot.background    = element_rect(fill = "white", colour = NA),
      plot.margin        = margin(S$mar, S$mar, S$mar, S$mar)
    )
}

## one save routine for every figure in sections 12 to 19
save_fig <- function(p, S, name, w, h) {
  outdir <- file.path(agg, S$dir); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  for (dv in S$devices) {
    f <- file.path(outdir, paste0(name, ".", dv))
    tryCatch({
      if (dv == "pdf") ggsave(f, p, width = w, height = h, device = PDF_DEV, bg = "white")
      else             ggsave(f, p, width = w, height = h, dpi = S$dpi, bg = "white")
    }, error = function(e) log_msg("  fig failed (%s): %s", basename(f), conditionMessage(e)))
  }
  cat(sprintf("  %-40s %5.1f x %4.1f in  -> %s\n", name, w, h, S$dir))
}

## ============================================================================
## 12.2  SHARED HELPERS
## ============================================================================
rd <- function(f) {
  p <- file.path(agg, f)
  if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else NULL
}

nice <- c(Hypoxia = "Hypoxia", CollagenSynthesis = "Collagen synthesis",
          iCAF = "iCAF", mCAF = "mCAF", Epithelial = "Epithelial",
          CollagenMaturation = "Collagen maturation", ECMCrosslink = "ECM crosslinking",
          ECMGlycoprotein = "ECM glycoprotein", ECMRemodel = "ECM remodelling",
          ERstress = "ER stress", ERAD = "ERAD")
relabel <- function(x) unname(ifelse(x %in% names(nice), nice[x], x))

XLAB_OR <- "Odds ratio for P4HA1 detection (95% CI)"

sig_flag <- function(df) ifelse(df$lower > 1 | df$upper < 1, "excl", "incl")

mk_forest <- function(df, S, xlim, breaks, ycol = "term", colour = CLR$pos,
                      show_stats = TRUE, xlab = NULL, ylab = NULL) {
  df$ypos <- df[[ycol]]
  df$sig  <- sig_flag(df)
  p <- ggplot(df, aes(OR, ypos)) +
    geom_vline(xintercept = 1, linetype = "22", linewidth = S$lwd * 0.7, colour = "grey40") +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0,
                   linewidth = S$lwd, colour = "black") +
    geom_point(aes(fill = sig), shape = 21, size = S$pt, stroke = S$stroke, colour = "black") +
    scale_fill_manual(values = c(excl = colour, incl = "white"), guide = "none") +
    ## v5.2: coord_cartesian, not scale limits. Scale limits DELETE any error
    ## bar that reaches past them (ER stress in Figure 3, ERAD High in Fig 2B).
    scale_x_log10(breaks = breaks, labels = function(x) sub("^0\\.", ".", x)) +
    coord_cartesian(xlim = xlim) +
    labs(x = xlab, y = ylab) +
    theme_fig(S)
  if (show_stats) {
    lab_df <- transform(df, .xpos = xlim[2])
    p <- p + geom_text(data = lab_df,
                       aes(x = .xpos, label = sprintf("%.2f (%.2f-%.2f)", OR, lower, upper)),
                       hjust = 1, vjust = -1.15, size = S$lab, colour = "grey20")
  }
  p
}

mk_strat <- function(df, S, xlim, breaks, colour = CLR$pos, xlab = NULL,
                     ylab = "Hypoxia tertile", nrow_facet = 1) {
  df$sig <- sig_flag(df)
  ggplot(df, aes(OR, hyp_bin)) +
    geom_vline(xintercept = 1, linetype = "22", linewidth = S$lwd * 0.7, colour = "grey40") +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, linewidth = S$lwd) +
    geom_point(aes(fill = sig), shape = 21, size = S$pt, stroke = S$stroke, colour = "black") +
    scale_fill_manual(values = c(excl = colour, incl = "white"), guide = "none") +
    scale_x_log10(breaks = breaks, labels = function(x) sub("^0\\.", ".", x)) +
    coord_cartesian(xlim = xlim) +                                   # v5.2, see mk_forest
    facet_wrap(~ term, nrow = nrow_facet) +
    geom_text(aes(label = sprintf("%.2f", OR)), vjust = -1.35,
              size = S$lab, colour = "grey20") +
    labs(x = xlab, y = ylab) +
    theme_fig(S)
}

## ============================================================================
## 12.3  READ THE POOLED RESULTS  (no re-analysis below this line)
## ============================================================================
pooled_fig <- rd("POOLED_logistic_forest.csv")         # section 6, marginal
sadj_fig   <- rd("STRAT_adjusted_pooled.csv")          # section 10 A
sstr_fig   <- rd("STRAT_stratified_pooled.csv")        # section 10 B
sint_fig   <- rd("STRAT_interaction_pooled.csv")       # section 10 C
padj_fig   <- rd("PROT_adjusted_pooled.csv")           # section 11 B
pstr_fig   <- rd("PROT_stratified_pooled.csv")         # section 11 C
pint_fig   <- rd("PROT_interaction_pooled.csv")        # section 11 D
permod_fig <- rd("per_sample_module_correlations.csv") # section 5
hrange     <- rd("HYPOXIA_SCORE_RANGE.csv")

k_of <- function(df) if (!is.null(df) && "k" %in% names(df)) max(df$k, na.rm = TRUE) else NA_integer_
sub_k <- function(df) {
  k <- k_of(df)
  if (is.na(k)) NULL else sprintf("Random-effects meta-analysis, k = %d %s", as.integer(k), K_UNIT)
}

## ============================================================================
## 12.4  RENDER EVERY FIGURE IN ONE STYLE
## ============================================================================
render_style <- function(S) {
  
  ttl <- function(p, t, s = NULL) if (isTRUE(S$titles)) p + labs(title = t, subtitle = s) else p
  tg  <- function(p, t) if (isTRUE(S$tags)) p + labs(tag = t) else p
  
  ## -- FIGURE 1: what predicts P4HA1 ----------------------------------------
  if (!is.null(pooled_fig) && !is.null(sadj_fig) && !is.null(sstr_fig)) {
    XL <- c(0.28, 4.2); BR <- c(0.5, 1, 2, 4)
    
    a <- pooled_fig %>% mutate(term = relabel(module)) %>% arrange(OR) %>%
      mutate(term = factor(term, levels = term))
    pA <- mk_forest(a, S, XL, BR, xlab = if (S$split) XLAB_OR else NULL)
    pA <- tg(ttl(pA, "Each predictor on its own", sub_k(pooled_fig)), "A")
    
    b <- sadj_fig %>% mutate(term = relabel(term)) %>% arrange(OR) %>%
      mutate(term = factor(term, levels = term))
    pB <- mk_forest(b, S, XL, BR, xlab = if (S$split) XLAB_OR else NULL)
    pB <- tg(ttl(pB, "All predictors in one model", sub_k(sadj_fig)), "B")
    
    cdf <- sstr_fig %>% mutate(term = relabel(term),
                               hyp_bin = factor(hyp_bin, BIN_LABELS, BIN_NICE))
    pC <- mk_strat(cdf, S, XL, BR, xlab = XLAB_OR)
    pC <- tg(ttl(pC, "Within hypoxia tertiles (relative, within section)",
                 "Holding hypoxia level constant"), "C")
    
    if (isTRUE(S$split)) {
      save_fig(pA, S, "Fig1A_marginal",   S$w2, S$h_full)
      save_fig(pB, S, "Fig1B_adjusted",   S$w2, S$h_full)
      save_fig(pC, S, "Fig1C_stratified", S$w2, S$h_full)
    } else {
      save_fig((pA / pB / pC) + plot_layout(heights = c(1.05, 0.85, 1.15)),
               S, "Figure1_main", S$w2, MM(135))
    }
  } else {
    cat("  Figure 1 skipped: a pooled CSV from section 6 or 10 is missing.\n")
  }
  
  ## -- FIGURE 2: proteostasis, strict ER stress module ----------------------
  if (!is.null(padj_fig) && !is.null(pstr_fig) && "variant" %in% names(padj_fig)) {
    XL <- c(0.25, 8); BR <- c(0.5, 1, 2, 4, 8)
    
    a <- padj_fig %>% filter(variant == PRIMARY_VARIANT) %>%
      mutate(term = relabel(term)) %>% arrange(OR) %>%
      mutate(term = factor(term, levels = term))
    pA <- mk_forest(a, S, XL, BR, xlab = if (S$split) XLAB_OR else NULL)
    pA <- tg(ttl(pA, "Adjusted model with proteostasis modules",
                 "Strict 24-gene ER stress module. Exploratory."), "A")
    
    lv <- c("ER stress", "ERAD", "Collagen synthesis", "iCAF")
    b <- pstr_fig %>% filter(variant == PRIMARY_VARIANT) %>%
      mutate(term = relabel(term),
             hyp_bin = factor(hyp_bin, BIN_LABELS, BIN_NICE)) %>%
      mutate(term = factor(term, levels = intersect(lv, unique(term))))
    pB <- mk_strat(b, S, XL, BR, xlab = XLAB_OR)
    pB <- tg(ttl(pB, "Within hypoxia tertiles",
                 "Strict 24-gene ER stress module. Exploratory."), "B")
    
    if (isTRUE(S$split)) {
      save_fig(pA, S, "Fig2A_proteostasis_adjusted",   S$w2, S$h_full)
      save_fig(pB, S, "Fig2B_proteostasis_stratified", S$w2, S$h_full)
    } else {
      save_fig((pA / pB) + plot_layout(heights = c(1, 1.05)),
               S, "Figure2_proteostasis", S$w2, MM(100))
    }
  } else {
    cat("  Figure 2 skipped: a pooled CSV from section 11 is missing.\n")
  }
  
  ## -- FIGURE 3: ER stress gene-set sensitivity -----------------------------
  if (!is.null(padj_fig) && dplyr::n_distinct(padj_fig$variant) == 2) {
    ord <- padj_fig %>% filter(variant == PRIMARY_VARIANT) %>% arrange(OR) %>%
      pull(term) %>% relabel()
    s <- padj_fig %>%
      mutate(term = factor(relabel(term), levels = ord),
             variant = factor(variant, c("ERstress", "ERstressStrict"),
                              c("Full (27 genes)", "Strict (24 genes)")))
    s$sig <- sig_flag(s)
    
    p3 <- ggplot(s, aes(OR, term, colour = variant)) +
      geom_vline(xintercept = 1, linetype = "22", linewidth = S$lwd * 0.7, colour = "grey40") +
      geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, linewidth = S$lwd,
                     position = position_dodge(width = .55)) +
      geom_point(aes(shape = sig), fill = "white", size = S$pt, stroke = S$stroke,
                 position = position_dodge(width = .55)) +
      scale_shape_manual(values = c(excl = 16, incl = 21), guide = "none") +
      scale_colour_manual(values = c("Full (27 genes)"   = CLR$full,
                                     "Strict (24 genes)" = CLR$neg), name = NULL) +
      scale_x_log10(breaks = c(0.5, 1, 2, 4, 8),
                    labels = function(x) sub("^0\\.", ".", x)) +
      coord_cartesian(xlim = c(0.35, 8)) +                           # v5.2: CIs were deleted at 5
      labs(x = "Adjusted odds ratio for P4HA1 detection (95% CI)", y = NULL) +
      theme_fig(S) + theme(legend.position = "top")
    p3 <- ttl(p3, "ER stress module with and without hypoxia-responsive genes",
              "Removing HYOU1, DDIT4, ERO1A. The strict estimate is the one to report.")
    
    save_fig(p3, S, "Figure3_sensitivity",
             if (S$split) S$w2 else S$w1 * 1.6,
             if (S$split) S$h_full else MM(62))
  } else {
    cat("  Figure 3 skipped: PROT_adjusted_pooled.csv does not contain both variants.\n")
  }
  
  ## -- FIGURE 4 (old): effect modification by hypoxia -----------------------
  ##    SUPERSEDED by section 13 (one joint model, correct bands, a p value).
  ##    Kept only to show what changed. Do not present it on its own.
  if (is.null(hrange) || !any(is.finite(hrange$q05))) {
    cat("  Figure 4 (old) skipped: HYPOXIA_SCORE_RANGE.csv not available.\n")
  } else if (!is.null(sint_fig) && !is.null(pint_fig)) {
    lo <- median(hrange$q05, na.rm = TRUE); hi <- median(hrange$q95, na.rm = TRUE)
    zz <- seq(lo, hi, length.out = 200)
    
    band <- function(df, main_term, variant = NULL, label) {
      d <- df
      if (!is.null(variant) && "variant" %in% names(d)) d <- d[d$variant == variant, ]
      im <- grepl(":", d$term)
      b0 <- d$est[d$term == main_term]; bi <- d$est[im]
      if (length(b0) != 1 || length(bi) != 1) return(NULL)
      se0 <- (log(d$upper[d$term == main_term]) - log(d$lower[d$term == main_term])) / (2 * 1.96)
      sei <- (log(d$upper[im]) - log(d$lower[im])) / (2 * 1.96)
      est <- b0 + bi * zz
      se  <- sqrt(se0^2 + (zz^2) * sei^2)
      data.frame(z = zz, OR = exp(est), lo = exp(est - 1.96 * se),
                 hi = exp(est + 1.96 * se), panel = label, cross = -b0 / bi)
    }
    
    sl <- bind_rows(
      band(sint_fig, "CollagenSynthesis", label = "Collagen synthesis"),
      band(pint_fig, "ERstress", variant = PRIMARY_VARIANT, label = "ER stress"))
    
    if (!is.null(sl) && nrow(sl)) {
      cr <- sl %>% group_by(panel) %>%
        summarise(cross = dplyr::first(cross), .groups = "drop") %>%
        filter(cross >= lo, cross <= hi)
      
      p4 <- ggplot(sl, aes(z, OR)) +
        geom_hline(yintercept = 1, linetype = "22", linewidth = S$lwd * 0.7, colour = "grey40") +
        geom_ribbon(aes(ymin = lo, ymax = hi, fill = panel), alpha = .16, colour = NA) +
        geom_line(aes(colour = panel), linewidth = S$lwd * 1.4) +
        scale_colour_manual(values = c("Collagen synthesis" = CLR$neg, "ER stress" = CLR$pos), name = NULL) +
        scale_fill_manual(values   = c("Collagen synthesis" = CLR$neg, "ER stress" = CLR$pos), name = NULL) +
        scale_y_log10() + coord_cartesian(xlim = c(lo, hi)) +
        labs(x = "Hypoxia module score (spot level, 5th to 95th percentile)",
             y = XLAB_OR,
             caption = "Superseded by section 13. Bands ignore the covariance between the two pooled coefficients.") +
        theme_fig(S) + theme(legend.position = "top")
      if (nrow(cr))
        p4 <- p4 + geom_vline(data = cr, aes(xintercept = cross),
                              linetype = "dotted", linewidth = S$lwd * 0.7, colour = "grey30")
      p4 <- ttl(p4, "How the collagen and ER stress effects change with hypoxia (superseded)",
                "Slopes from two separate models, over the observed score range only")
      
      save_fig(p4, S, "Figure4_effect_modification_OLD",
               if (S$split) S$w2 else S$w1 * 1.5,
               if (S$split) S$h_full else MM(62))
    }
  }
  
  ## -- FIGURE 5: per-patient consistency ------------------------------------
  ##    v5 (C4): one point per patient (mean of its sections).
  if (!is.null(permod_fig) && all(c("module", "rho", "sample_id") %in% names(permod_fig))) {
    ## v5.2: CollagenMaturation contains P4HA1 itself (section 1), so its rho
    ## with P4HA1 is partly P4HA1 against itself. Kept in the CSV, not drawn.
    pm <- permod_fig[is.finite(permod_fig$rho) & permod_fig$module != "CollagenMaturation", , drop = FALSE] %>%
      mutate(unit = unit_of(sample_id)) %>%
      group_by(module, unit) %>% summarise(rho = mean(rho), .groups = "drop")
    pm$module_nice <- relabel(pm$module)
    ordm <- pm %>% group_by(module_nice) %>%
      summarise(med = median(rho, na.rm = TRUE), .groups = "drop") %>%
      arrange(med) %>% pull(module_nice)
    pm$module_nice <- factor(pm$module_nice, levels = ordm)
    n_un <- dplyr::n_distinct(pm$unit)
    
    p5 <- ggplot(pm, aes(module_nice, rho)) +
      geom_hline(yintercept = 0, linetype = "22", linewidth = S$lwd * 0.7, colour = "grey40") +
      geom_boxplot(outlier.shape = NA, fill = "grey92", linewidth = S$lwd * 0.8, width = .6) +
      geom_jitter(width = .13, size = S$pt * 0.55, alpha = .5, colour = CLR$pos) +
      coord_flip() +
      labs(x = NULL, y = sprintf("Per-%s Spearman rho (P4HA1 vs module score)", UNIT_LABEL),
           caption = sprintf("One point per %s (mean of its sections), n = %d.", UNIT_LABEL, n_un)) +
      theme_fig(S) +
      theme(panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(linewidth = S$grid, colour = "grey90"))
    p5 <- ttl(p5, "How consistent is each association across tumours?",
              sprintf("One point per %s, n = %d", UNIT_LABEL, n_un))
    
    save_fig(p5, S, "Figure5_consistency",
             if (S$split) S$w2 else S$w1 * 1.6,
             if (S$split) S$h_full else MM(70))
  } else {
    cat("  Figure 5 skipped: per_sample_module_correlations.csv missing or malformed.\n")
  }
  
  invisible(TRUE)
}

## ============================================================================
## 12.5  RUN BOTH STYLES
## ============================================================================
for (sn in names(FIG_STYLES)) {
  S <- FIG_STYLES[[sn]]
  cat(sprintf("\n--- %s figures -> %s ---\n", toupper(sn),
              normalizePath(file.path(agg, S$dir), mustWork = FALSE)))
  render_style(S)
}

cat("\nfigures_journal : PDFs are the versions to submit.\n")
cat("figures_slides  : 16:9 PNGs, one panel per slide, titles baked in.\n")
cat("Section 11 output remains EXPLORATORY pending sign-off on the gene sets.\n")

log_msg("SECTIONS 0 TO 12 COMPLETE. Everything under: %s", normalizePath(OUT_ROOT))

################################################################################
## 14. DETECTION AUDIT FOR P4HA1, P4HA2 AND P4HA3
##     P4HA1, P4HA2 and P4HA3 all encode collagen P4H alpha subunits; P4HA3
##     also forms an active enzyme that hydroxylates collagen (Kukkola et al.,
##     J Biol Chem 2003, doi:10.1074/jbc.M306806200). A null result for
##     P4HA3 only means "no effect" if P4HA3 is detected often enough to be
##     modelled. This section counts and does nothing else: no module scores,
##     no glm, no pooling. It writes unsuffixed files deliberately (it fits
##     nothing, so ADJUST_DEPTH does not apply).
################################################################################

stopifnot(exists("work"), exists("agg"), exists("log_msg"), exists("bind"))

## ============================================================================
## 14.0  SETTINGS
## ============================================================================
P4HA_GENES <- c("P4HA1", "P4HA2", "P4HA3")

## legacy symbol aliases to look for if the primary symbol is absent
ALIASES <- list(
  P4HA1  = character(0),
  P4HA2  = character(0),
  P4HA3  = character(0),
  ERO1A  = "ERO1L")

## ~10 events per predictor; the joint model has 6 predictors + intercept.
## Sections below this are flagged, not deleted.
MIN_EVENTS_PER_TUMOUR <- 70

JOINT_CACHE_DETECTION <- TRUE   # cache per-gene detection rates (sections 16, 19, 25)
DETECT_CACHE_MIN_RATE <- 0.005  # do not cache genes under this rate

## ============================================================================
## 14.1  CIRCULARITY RE-AUDIT FOR THE NEW OUTCOMES
## ============================================================================
cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 14.1  CIRCULARITY RE-AUDIT FOR THE NEW OUTCOMES\n")
cat(strrep("=", 74), "\n", sep = "")

PREDICTOR_MODULES <- c("Hypoxia", "CollagenSynthesis", "iCAF",
                       intersect(c("ERstressStrict", "ERstress", "ERAD"), names(MODULES)))
for (g in P4HA_GENES) {
  hits <- names(MODULES)[vapply(MODULES, function(m) g %in% m, logical(1))]
  used <- intersect(hits, PREDICTOR_MODULES)
  cat(sprintf("\n  %s appears in modules: %s\n", g,
              if (length(hits)) paste(hits, collapse = ", ") else "none"))
  cat(sprintf("     of which used as PREDICTORS: %s\n",
              if (length(used)) paste0("*** CIRCULAR: ", paste(used, collapse = ", "), " ***")
              else "none, safe"))
}
cat("\n  CollagenMaturation contains P4HA1 and P4HA2 and is excluded from every model.\n")

## ============================================================================
## 14.2  COUNT, ONE SECTION AT A TIME
## ============================================================================
detect_one <- function(rds_path, sid) {
  obj <- readRDS(rds_path)
  obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
  DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else Assays(obj)[1]
  rn <- rownames(obj)
  n  <- ncol(obj)
  
  resolve <- function(g) {
    if (g %in% rn) return(g)
    alt <- intersect(ALIASES[[g]], rn)
    if (length(alt)) return(alt[1])
    NA_character_
  }
  
  cnt_col <- grep("^nCount_", colnames(obj@meta.data), value = TRUE)
  depth   <- if (length(cnt_col)) obj@meta.data[[cnt_col[1]]] else rep(NA_real_, n)
  feat_col <- grep("^nFeature_", colnames(obj@meta.data), value = TRUE)
  nfeat    <- if (length(feat_col)) obj@meta.data[[feat_col[1]]] else rep(NA_real_, n)
  
  qs <- function(v, q) if (all(is.na(v))) NA_real_ else unname(quantile(v, q, na.rm = TRUE))
  depth_row <- data.frame(
    sample_id = sid, n_spots = n,
    depth_col = if (length(cnt_col)) cnt_col[1] else NA_character_,
    depth_q25 = qs(depth, .25), depth_median = qs(depth, .50), depth_q75 = qs(depth, .75),
    nFeature_median = qs(nfeat, .50), row.names = NULL)
  
  pos_list <- list()
  rows <- do.call(rbind, lapply(P4HA_GENES, function(g) {
    sym <- resolve(g)
    if (is.na(sym)) {
      return(data.frame(sample_id = sid, gene = g, symbol_found = NA_character_,
                        present = FALSE, n_spots = n, n_detected = NA_integer_,
                        detect_rate = NA_real_, mean_expr_pos = NA_real_,
                        depth_rho = NA_real_, row.names = NULL))
    }
    v   <- as.numeric(FetchData(obj, vars = sym, layer = "data")[, 1])
    pos <- v > 0
    pos_list[[g]] <<- pos
    rho <- if (all(is.na(depth)) || length(unique(pos)) < 2) NA_real_ else
      suppressWarnings(cor(as.numeric(pos), depth, method = "spearman"))
    data.frame(sample_id = sid, gene = g, symbol_found = sym, present = TRUE,
               n_spots = n, n_detected = sum(pos), detect_rate = mean(pos),
               mean_expr_pos = if (any(pos)) mean(v[pos]) else NA_real_,
               depth_rho = rho, row.names = NULL)
  }))
  
  co <- data.frame(sample_id = sid, jaccard_A1_A2 = NA_real_, both = NA_integer_)
  if (!is.null(pos_list$P4HA1) && !is.null(pos_list$P4HA2)) {
    a <- pos_list$P4HA1; b <- pos_list$P4HA2
    u <- sum(a | b)
    co$both <- sum(a & b)
    co$jaccard_A1_A2 <- if (u > 0) sum(a & b) / u else NA_real_
  }
  
  ero <- data.frame(sample_id = sid,
                    ERO1A_present = "ERO1A" %in% rn,
                    ERO1L_present = "ERO1L" %in% rn)
  
  rates <- NULL
  if (JOINT_CACHE_DETECTION) {
    m <- tryCatch(GetAssayData(obj, assay = DefaultAssay(obj), layer = "counts"),
                  error = function(e) NULL)
    if (is.null(m) || !nrow(m))
      m <- tryCatch(GetAssayData(obj, assay = DefaultAssay(obj), layer = "data"),
                    error = function(e) NULL)
    if (!is.null(m) && nrow(m)) {
      r <- Matrix::rowSums(m > 0) / ncol(m)
      r <- r[r >= DETECT_CACHE_MIN_RATE]
      if (length(r))
        rates <- data.frame(sample_id = sid, gene = names(r),
                            detect_rate = as.numeric(r), row.names = NULL)
    }
  }
  
  list(genes = rows, depth = depth_row, codetect = co, ero = ero, rates = rates)
}

acc14 <- list(genes = list(), depth = list(), codetect = list(), ero = list(), rates = list())
for (i in seq_len(nrow(work))) {
  path <- work$path[i]; sid <- work$sample_id[i]
  if (is.na(path) || !file.exists(path)) { log_msg("[detect %d/%d] %s FILE MISSING", i, nrow(work), sid); next }
  log_msg("[detect %d/%d] %s", i, nrow(work), sid)
  r <- tryCatch(detect_one(path, sid),
                error = function(e) { log_msg("  detect failed: %s", conditionMessage(e)); NULL })
  if (is.null(r)) next
  for (k in names(acc14))
    if (!is.null(r[[k]]) && nrow(r[[k]])) acc14[[k]][[length(acc14[[k]]) + 1]] <- r[[k]]
}

det_all   <- bind(acc14$genes)
depth_all <- bind(acc14$depth)
co_all    <- bind(acc14$codetect)
ero_all   <- bind(acc14$ero)

if (nrow(det_all) == 0) stop("Section 14: no section could be read.")

write.csv(det_all,   file.path(agg, "DETECT_p4ha_per_sample.csv"), row.names = FALSE)
write.csv(depth_all, file.path(agg, "DETECT_depth_per_sample.csv"), row.names = FALSE)

if (JOINT_CACHE_DETECTION && length(acc14$rates)) {
  rates_all <- bind(acc14$rates)
  write.csv(rates_all, file.path(agg, "DETECT_cohort_gene_rates.csv"), row.names = FALSE)
  log_msg("Cached detection rates for %d gene-by-section combinations.", nrow(rates_all))
}

## ============================================================================
## 14.3  REPORT
## ============================================================================
cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 14  DETECTION AUDIT\n")
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("Sections read: %d (%d patients)\n", dplyr::n_distinct(det_all$sample_id),
            length(unique(patient_of(unique(det_all$sample_id))))))

cat("\n--- is each gene in the object at all? ---\n")
for (g in P4HA_GENES) {
  s <- det_all[det_all$gene == g, ]
  cat(sprintf("  %-6s present in %2d / %2d sections", g, sum(s$present), nrow(s)))
  syms <- unique(na.omit(s$symbol_found))
  cat(sprintf("   symbol(s): %s\n", if (length(syms)) paste(syms, collapse = ", ") else "NOT FOUND"))
}

cat("\n--- detection rate per gene (fraction of spots with count > 0) ---\n")
det_ok <- det_all[det_all$present, , drop = FALSE]
det_ok$usable <- det_ok$n_detected >= MIN_EVENTS_PER_TUMOUR &
  (det_ok$n_spots - det_ok$n_detected) >= MIN_EVENTS_PER_TUMOUR
det_summary <- det_ok %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    k_sections       = dplyr::n(),
    k_patients       = length(unique(patient_of(sample_id))),
    median_rate      = median(detect_rate, na.rm = TRUE),
    min_rate         = min(detect_rate, na.rm = TRUE),
    max_rate         = max(detect_rate, na.rm = TRUE),
    median_events    = median(n_detected, na.rm = TRUE),
    min_events       = min(n_detected, na.rm = TRUE),
    n_sections_ok    = sum(usable, na.rm = TRUE),
    n_patients_ok    = length(unique(patient_of(sample_id[usable %in% TRUE]))),
    median_depth_rho = median(depth_rho, na.rm = TRUE),
    .groups = "drop")
write.csv(det_summary, file.path(agg, "DETECT_p4ha_summary.csv"), row.names = FALSE)

for (i in seq_len(nrow(det_summary))) with(det_summary[i, ],
                                           cat(sprintf("  %-6s median %5.1f%% (range %4.1f%% to %5.1f%%)  median %5.0f events (min %4.0f)  usable in %2d/%2d sections, %2d/%2d patients\n",
                                                       gene, 100*median_rate, 100*min_rate, 100*max_rate,
                                                       median_events, min_events, n_sections_ok, k_sections, n_patients_ok, k_patients)))
cat(sprintf("  \"usable\" = at least %d detected AND %d undetected spots (rule of thumb).\n",
            MIN_EVENTS_PER_TUMOUR, MIN_EVENTS_PER_TUMOUR))

cat("\n--- per section, side by side ---\n")
wide <- det_all %>%
  dplyr::select(sample_id, gene, detect_rate, n_detected) %>%
  tidyr::pivot_wider(names_from = gene, values_from = c(detect_rate, n_detected))
if ("detect_rate_P4HA1" %in% names(wide)) wide <- dplyr::arrange(wide, dplyr::desc(detect_rate_P4HA1))
gcol <- function(r, nm) if (nm %in% names(r)) r[[nm]] else NA_real_
cat(sprintf("  %-30s %-10s %8s %8s %8s %10s\n", "section", "patient", "A1 %", "A2 %", "A3 %", "A3 events"))
for (i in seq_len(nrow(wide))) {
  r <- wide[i, ]
  f  <- function(x) if (is.null(x) || length(x) == 0 || is.na(x)) "       ." else sprintf("%8.1f", 100*x)
  fe <- function(x) if (is.null(x) || length(x) == 0 || is.na(x)) "         ." else sprintf("%10.0f", x)
  a3 <- gcol(r, "n_detected_P4HA3")
  low3 <- !is.na(a3) && a3 < MIN_EVENTS_PER_TUMOUR
  cat(sprintf("  %-30s %-10s %s %s %s %s%s\n", r$sample_id, patient_of(r$sample_id),
              f(gcol(r, "detect_rate_P4HA1")), f(gcol(r, "detect_rate_P4HA2")),
              f(gcol(r, "detect_rate_P4HA3")), fe(a3),
              if (low3) "   <-- too few for P4HA3" else ""))
}

cat("\n--- does detection just follow sequencing depth? ---\n")
for (g in P4HA_GENES) {
  s <- det_all[det_all$gene == g & det_all$present & is.finite(det_all$depth_rho), ]
  if (!nrow(s)) next
  cat(sprintf("  %-6s median Spearman rho(detected, nCount) = %+.2f  (range %+.2f to %+.2f)\n",
              g, median(s$depth_rho), min(s$depth_rho), max(s$depth_rho)))
}
cat("  The higher this is, the more depth adjustment matters. Set ADJUST_DEPTH\n")
cat("  <- TRUE in block 12b to answer it directly.\n")

if (nrow(co_all) && any(is.finite(co_all$jaccard_A1_A2))) {
  cat("\n--- do P4HA1 and P4HA2 mark the same spots? ---\n")
  cat(sprintf("  median Jaccard overlap = %.2f  (range %.2f to %.2f)\n",
              median(co_all$jaccard_A1_A2, na.rm = TRUE),
              min(co_all$jaccard_A1_A2, na.rm = TRUE),
              max(co_all$jaccard_A1_A2, na.rm = TRUE)))
  cat("  High overlap: P4HA2 is close to a repeat of P4HA1, so agreement between\n")
  cat("  them is weak evidence. Low overlap: closer to independent tests.\n")
}

cat("\n--- ERO1A vs ERO1L ---\n")
if (nrow(ero_all)) {
  cat(sprintf("  ERO1A present in %d / %d sections\n", sum(ero_all$ERO1A_present), nrow(ero_all)))
  cat(sprintf("  ERO1L present in %d / %d sections\n", sum(ero_all$ERO1L_present), nrow(ero_all)))
  if (all(!ero_all$ERO1A_present) && any(ero_all$ERO1L_present)) {
    cat("  *** The pipeline asks for ERO1A and the objects carry ERO1L. The strict\n")
    cat("      ER stress module then removes TWO hypoxia-responsive genes, not three.\n")
    cat("      Fix the gene list before quoting the sensitivity analysis. ***\n")
  } else if (all(ero_all$ERO1A_present)) {
    cat("  ERO1A is the symbol in use; the strict module drops three genes as described.\n")
  } else {
    cat("  ERO1A is missing from some sections; the strict module differs between them.\n")
  }
}

## ============================================================================
## 14.4  VERDICT
## ============================================================================
cat("\n", strrep("=", 74), "\n", sep = "")
cat("VERDICT: CAN THE P4HA2 / P4HA3 TEST BE RUN?\n")
cat(strrep("=", 74), "\n", sep = "")
for (g in P4HA_GENES) {
  s <- det_summary[det_summary$gene == g, ]
  if (!nrow(s)) { cat(sprintf("  %-6s NOT PRESENT in any section. Test impossible.\n", g)); next }
  frac <- s$n_patients_ok / s$k_patients
  verdict <- if (frac >= 0.8) "fine, model as usual" else
    if (frac >= 0.5) "marginal, report k and per-section events alongside every estimate" else
      "UNDERPOWERED, a null result here would mean nothing"
  cat(sprintf("  %-6s usable in %2d/%2d patients (%.0f%%)  ->  %s\n",
              g, s$n_patients_ok, s$k_patients, 100*frac, verdict))
}
cat("\n  If P4HA3 is underpowered, report its CI and whether it excludes the P4HA1\n")
cat("  effect size (an equivalence argument), or state plainly that P4HA3 is too\n")
cat("  sparse in Visium to serve as a control. Never report \"no effect in P4HA3\"\n")
cat("  without the power statement.\n")
cat(strrep("=", 74), "\n", sep = "")

log_msg("Section 14 complete.")

################################################################################
## 12b. RUN STATE, PROVENANCE AND SHARED HELPERS FOR SECTIONS 13 TO 26
##
##   The filename suffix and the model formula must both come from
##   ADJUST_DEPTH at the same moment. Everything that guards against them
##   disagreeing is here, ahead of every section that could be affected.
##   DO NOT redefine any of these by hand in the console. Change ADJUST_DEPTH
##   here and re-run from the top of the file in a fresh session.
################################################################################

## ADJUST_DEPTH is the ONLY place this is set.
ADJUST_DEPTH <- FALSE

OUT_SUFFIX <- if (isTRUE(ADJUST_DEPTH)) "_depthadj" else ""
outf <- function(base, ext = "csv") sprintf("%s%s.%s", base, OUT_SUFFIX, ext)

RUN_ID <- paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), "-",
                 if (isTRUE(ADJUST_DEPTH)) "depthadj" else "plain")

stamp <- function(df) {
  if (is.null(df) || !nrow(df)) return(df)
  df$adjust_depth <- isTRUE(ADJUST_DEPTH)
  df$run_id       <- RUN_ID
  df
}

write_out <- function(df, base) {
  write.csv(stamp(df), file.path(agg, outf(base)), row.names = FALSE)
  invisible(NULL)
}

assert_run_state <- function(where) {
  problems <- character(0)
  
  if (!exists("ADJUST_DEPTH")) problems <- c(problems, "ADJUST_DEPTH is not defined")
  if (!exists("OUT_SUFFIX"))   problems <- c(problems, "OUT_SUFFIX is not defined")
  if (!exists("outf"))         problems <- c(problems, "outf() is not defined")
  if (!exists("RUN_ID"))       problems <- c(problems, "RUN_ID is not defined")
  if (length(problems)) stop(sprintf("%s: %s. Re-run block 12b.",
                                     where, paste(problems, collapse = "; ")))
  
  expect_suffix <- if (isTRUE(ADJUST_DEPTH)) "_depthadj" else ""
  if (!identical(OUT_SUFFIX, expect_suffix))
    problems <- c(problems, sprintf(
      "OUT_SUFFIX is '%s' but ADJUST_DEPTH = %s implies '%s'",
      OUT_SUFFIX, ADJUST_DEPTH, expect_suffix))
  
  test_name <- outf("__probe__")
  if (!identical(test_name, paste0("__probe__", expect_suffix, ".csv")))
    problems <- c(problems, sprintf(
      "outf() returns '%s', which does not match ADJUST_DEPTH = %s. It was probably redefined by hand.",
      test_name, ADJUST_DEPTH))
  
  flag_names <- c("JOINT_ADJUST_DEPTH", "ISO_ADJUST_DEPTH", "RND_ADJUST_DEPTH", "METAB_ADJUST_DEPTH",
                  "SC_ADJUST_DEPTH", "EL_ADJUST_DEPTH", "UPR_ADJUST_DEPTH", "SP_ADJUST_DEPTH")
  local_flags <- unlist(lapply(flag_names, function(f) if (exists(f)) setNames(get(f), f) else NULL))
  bad <- names(local_flags)[local_flags != ADJUST_DEPTH]
  if (length(bad))
    problems <- c(problems, sprintf(
      "%s disagree(s) with ADJUST_DEPTH = %s. A section was re-entered below its settings block.",
      paste(bad, collapse = ", "), ADJUST_DEPTH))
  
  if (length(problems)) {
    cat("\n", strrep("!", 74), "\n", sep = "")
    cat(sprintf("%s: RUN STATE IS INCONSISTENT. NOTHING WAS FITTED.\n", where))
    for (p in problems) cat("  - ", p, "\n", sep = "")
    cat("\n  Fix: start a fresh R session and run from section 0.\n")
    cat(strrep("!", 74), "\n", sep = "")
    stop("Run state inconsistent; see the message above.")
  }
  
  cat(sprintf("\n[run state OK] %s | ADJUST_DEPTH = %s | suffix = '%s' | run_id = %s\n",
              where, ADJUST_DEPTH, OUT_SUFFIX, RUN_ID))
  invisible(TRUE)
}

write_run_manifest <- function() {
  gene_hash <- function(g) paste0(length(g), ":", substr(paste(sort(g), collapse = ","), 1, 60))
  mf <- data.frame(
    key = c("run_id", "timestamp", "adjust_depth", "seed", "out_suffix", "pool_level",
            "include_metastasis", "R_version", "Seurat_version",
            "n_sections_in_worklist", "n_patients_in_worklist",
            "module_Hypoxia", "module_CollagenSynthesis", "module_iCAF",
            "module_ERstressStrict", "module_ERAD", "joint_formula"),
    value = c(RUN_ID, format(Sys.time()), as.character(ADJUST_DEPTH),
              as.character(SEED), OUT_SUFFIX, POOL_LEVEL, as.character(INCLUDE_METASTASIS),
              paste(R.version$major, R.version$minor, sep = "."),
              as.character(utils::packageVersion("Seurat")),
              as.character(nrow(work)), as.character(length(unique(work$patient_id))),
              gene_hash(MODULES$Hypoxia), gene_hash(MODULES$CollagenSynthesis),
              gene_hash(MODULES$iCAF), gene_hash(MODULES$ERstressStrict),
              gene_hash(MODULES$ERAD),
              paste0("y ~ zHyp*zCol + zHyp*zER + ziCAF + zERAD",
                     if (ADJUST_DEPTH) " + zDepth" else "")),
    stringsAsFactors = FALSE)
  write.csv(mf, file.path(agg, outf("RUN_manifest")), row.names = FALSE)
  cat("\n--- run manifest ---\n")
  for (i in seq_len(nrow(mf))) cat(sprintf("  %-26s %s\n", mf$key[i], mf$value[i]))
  invisible(mf)
}

## meta_pool returns est / lower / upper; recover SE and a Wald p from the
## interval. The statistic is called zstat because in 13.5 a grouping
## variable is called z.
pool_with_p <- function(yi, vi, sid = NULL) {
  m <- meta_pool(yi, vi, sid)
  if (!is.finite(m$est) || !is.finite(m$lower)) {
    m$SE <- NA_real_; m$zstat <- NA_real_; m$p <- NA_real_; return(m)
  }
  se <- (m$upper - m$lower) / (2 * 1.96)
  m$SE    <- se
  m$zstat <- if (is.finite(se) && se > 0) m$est / se else NA_real_
  m$p     <- if (is.finite(m$zstat)) 2 * pnorm(-abs(m$zstat)) else NA_real_
  m
}

fmt_p <- function(p) {
  if (length(p) != 1 || !is.finite(p)) return("NA")
  if (p < 2e-16) return("< 2e-16")
  if (p < 0.001)  return(sprintf("%.1e", p))
  sprintf("%.3f", p)
}

## k label used on every figure and report from here on
k_lab <- function(r) sprintf("k = %d %s (%d sections)", as.integer(r$k), K_UNIT, as.integer(r$k_sections))

## ---------------------------------------------------------------------------
## v5 (C3): SAFETY CHECK AGAINST SECTION 15, SECTION BY SECTION
##   Sections 20, 23, 24 and 26 refit the section 15 model (or a part of it)
##   and must reproduce it. They are compared per section on the sections
##   both fitted, and the pooled reference is recomputed on those sections.
##   A section fitted here but NOT in section 15 halts: the inclusion rules
##   then differ. A section fitted in 15 but not here (no ESTIMATE table, for
##   example) is listed and the comparison continues without it.
## ---------------------------------------------------------------------------
s15_fitted <- function(outcome) {
  if (!exists("iso_stat")) stop("Section 15 has not run in this session (iso_stat missing).")
  unique(iso_stat$sample_id[iso_stat$outcome == outcome & iso_stat$fitted])
}

check_against_s15 <- function(per_section, outcome, where, tol = 1e-6) {
  if (!exists("iso_con") || !exists("iso_con_pooled"))
    stop(where, " needs section 15 (iso_con, iso_con_pooled) in the same session for its safety check.")
  ref  <- iso_con[iso_con$outcome == outcome & iso_con$contrast == "interaction_diff",
                  c("sample_id", "est", "var")]
  here <- per_section[, c("sample_id", "est", "var")]
  extra <- setdiff(here$sample_id, ref$sample_id)
  lost  <- setdiff(ref$sample_id, here$sample_id)
  m <- merge(here, ref, by = "sample_id", suffixes = c("", "_s15"))
  dmax <- if (nrow(m)) max(abs(m$est - m$est_s15)) else NA_real_
  
  cat("\n", strrep("=", 74), "\n", sep = "")
  cat(sprintf("%s  SAFETY CHECK against section 15 (%s)\n", where, outcome))
  cat(strrep("=", 74), "\n", sep = "")
  cat(sprintf("  sections fitted here %d, in section 15 %d, in both %d\n",
              nrow(here), nrow(ref), nrow(m)))
  if (length(lost))
    cat(sprintf("  fitted in section 15 but not here (reason in the log): %s\n", paste(lost, collapse = ", ")))
  if (length(extra))
    stop(where, " halted: section(s) fitted here but not in section 15: ", paste(extra, collapse = ", "),
         ". The inclusion rules differ, so nothing below is comparable.")
  if (!nrow(m) || !is.finite(dmax) || dmax > tol)
    stop(sprintf(paste0("%s halted: per-section estimates differ from section 15 (largest difference %s). ",
                        "Scoring or model differs; nothing below is interpretable."), where, format(dmax)))
  ref_pool  <- pool_with_p(m$est_s15, m$var_s15, m$sample_id)
  here_pool <- pool_with_p(m$est, m$var, m$sample_id)
  cat(sprintf("  largest per-section difference %.1e\n", dmax))
  cat(sprintf("  pooled here %+.6f   section 15 on the same sections %+.6f   (%s)\n",
              here_pool$est, ref_pool$est, k_lab(ref_pool)))
  if (length(lost)) {
    full <- iso_con_pooled$est[iso_con_pooled$outcome == outcome & iso_con_pooled$contrast == "interaction_diff"]
    cat(sprintf("  section 15 on all its sections: %+.6f (differs only by the section(s) listed above)\n", full))
  }
  cat("  PASSED\n")
  invisible(ref_pool)
}


################################################################################
## 13. FIGURE 4, REBUILT FROM ONE JOINT MODEL, WITH THE P VALUE ON IT
##
##   One model per section containing both interactions:
##       y ~ zHyp * zCol + zHyp * zER + ziCAF + zERAD  [+ zDepth]
##   so every quantity the figure needs is a contrast on a single fit and
##   vcov() supplies all covariances. All predictors are z-scored within the
##   section, so every effect is per 1 SD.
##   Converging lines are not the same as lines that swap order; the report
##   and the figure say which one the data support.
################################################################################

stopifnot(exists("MODULES"), exists("meta_pool"), exists("work"), exists("agg"))
if (!all(c("ERstressStrict", "ERAD") %in% names(MODULES)))
  stop("Section 13 needs section 11 to have run first (ERstressStrict / ERAD not defined).")

## ============================================================================
## 13.0  SETTINGS
## ============================================================================
JOINT_ES           <- "ERstressStrict"
JOINT_MIN_SPOTS    <- 30
JOINT_ADJUST_DEPTH <- ADJUST_DEPTH       # derived, never set by hand
JOINT_Z_GRID       <- seq(-2, 2, length.out = 41)
JOINT_ALPHA        <- 0.05

JOINT_SCORE_MODULES <- c("Hypoxia", "CollagenSynthesis", "iCAF", JOINT_ES, "ERAD")

JOINT_RENAME <- c(Hypoxia = "zHyp", CollagenSynthesis = "zCol",
                  iCAF = "ziCAF", ERAD = "zERAD")
JOINT_RENAME[JOINT_ES] <- "zER"

T_COL <- "zCol"; T_ER <- "zER"; T_HYP <- "zHyp"
I_COL <- paste0(T_HYP, ":", T_COL)
I_ER  <- paste0(T_HYP, ":", T_ER)

LAB_COL <- "Collagen synthesis"
LAB_ER  <- "ER stress"
LAB_DIF <- "difference"

assert_run_state("Section 13")
write_run_manifest()

## ============================================================================
## 13.1  SCORE, THEN z-STANDARDISE WITHIN THE SECTION
##       Also returns a fingerprint of each module score so two runs can be
##       compared and any drift traced to scoring or to modelling.
## ============================================================================
joint_expr_one <- function(rds_path, sid) {
  obj <- readRDS(rds_path)
  obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
  DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else Assays(obj)[1]
  if (!"P4HA1" %in% rownames(obj)) return(NULL)
  
  fp <- list()
  for (nm in JOINT_SCORE_MODULES) {
    full <- MODULES[[nm]]
    g    <- intersect(full, rownames(obj))
    need <- min_genes_for(full)
    if (length(g) < need) {
      log_msg("  %s: %s has only %d/%d genes (need %d) -> section skipped",
              sid, nm, length(g), length(full), need)
      return(NULL)
    }
    set.seed(SEED)
    obj <- tryCatch(add_module_score_fallback(obj, features = list(g), name = nm,
                                              assay = DefaultAssay(obj), nbin = 20, ctrl = 50),
                    error = function(e) NULL)
    if (is.null(obj) || !paste0(nm, "1") %in% colnames(obj@meta.data)) return(NULL)
    names(obj@meta.data)[names(obj@meta.data) == paste0(nm, "1")] <- nm
    v <- obj@meta.data[[nm]]
    fp[[nm]] <- data.frame(sample_id = sid, module = nm,
                           n_genes_in_set = length(full), n_genes_used = length(g),
                           score_mean = mean(v, na.rm = TRUE),
                           score_sd   = sd(v,   na.rm = TRUE),
                           score_min  = min(v,  na.rm = TRUE),
                           score_max  = max(v,  na.rm = TRUE),
                           row.names = NULL)
  }
  
  vars <- c("P4HA1", JOINT_SCORE_MODULES)
  cnt_col <- grep("^nCount_", colnames(obj@meta.data), value = TRUE)
  if (JOINT_ADJUST_DEPTH) {
    if (!length(cnt_col)) {
      log_msg("  %s: no nCount_* column, depth adjustment impossible -> skipped", sid)
      return(NULL)
    }
    vars <- c(vars, cnt_col[1])
  }
  d <- FetchData(obj, vars = vars, layer = "data")
  d$y <- as.integer(d$P4HA1 > 0)
  
  zz <- function(x) { s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) return(NULL)
  as.numeric((x - mean(x, na.rm = TRUE)) / s) }
  for (nm in JOINT_SCORE_MODULES) {
    v <- zz(d[[nm]])
    if (is.null(v)) { log_msg("  %s: %s has zero variance -> section skipped", sid, nm); return(NULL) }
    d[[ JOINT_RENAME[[nm]] ]] <- v
  }
  if (JOINT_ADJUST_DEPTH) {
    v <- zz(log1p(d[[cnt_col[1]]]))   # log1p: counts are strongly right skewed
    if (is.null(v)) return(NULL)
    d$zDepth <- v
  }
  
  attr(d, "fingerprint") <- do.call(rbind, fp)
  attr(d, "n_spots")     <- nrow(d)
  d
}

## ============================================================================
## 13.2  FIT ONE SECTION
## ============================================================================
fit_joint_one <- function(d, sid) {
  out <- list(coef = NULL, contrast = NULL, curve = NULL, modelcheck = NULL)
  if (is.null(d) || nrow(d) < JOINT_MIN_SPOTS || length(unique(d$y)) != 2) return(out)
  
  rhs <- paste0(T_HYP, " * ", T_COL, " + ", T_HYP, " * ", T_ER, " + ziCAF + zERAD",
                if (JOINT_ADJUST_DEPTH) " + zDepth" else "")
  f <- tryCatch(glm(as.formula(paste("y ~", rhs)), data = d, family = binomial()),
                error = function(e) NULL)
  if (is.null(f) || !f$converged) return(out)
  
  b <- coef(f); V <- vcov(f)
  
  ## read the model, not the flag
  has_depth <- "zDepth" %in% names(b)
  if (!identical(has_depth, isTRUE(JOINT_ADJUST_DEPTH)))
    stop(sprintf(paste0("Section 13, section %s: JOINT_ADJUST_DEPTH = %s but the fitted model %s ",
                        "a zDepth term. Start a fresh session from section 0."),
                 sid, JOINT_ADJUST_DEPTH, if (has_depth) "DOES contain" else "does NOT contain"))
  out$modelcheck <- data.frame(sample_id = sid, has_zDepth = has_depth,
                               n_coef = length(b), formula = rhs, row.names = NULL)
  
  need <- c(T_COL, T_ER, I_COL, I_ER)
  if (!all(need %in% names(b)) || any(!is.finite(b[need]))) return(out)
  
  s  <- summary(f)$coefficients
  tm <- setdiff(rownames(s), "(Intercept)")
  out$coef <- data.frame(sample_id = sid, term = tm,
                         log_OR = s[tm, "Estimate"], SE = s[tm, "Std. Error"],
                         p = s[tm, "Pr(>|z|)"], row.names = NULL)
  
  ## estimate L b, variance L V L' : every covariance is included
  Lvec <- function(w) { L <- setNames(numeric(length(b)), names(b))
  for (nm in names(w)) L[nm] <- w[[nm]]; L }
  est_of <- function(L) as.numeric(L %*% b)
  var_of <- function(L) as.numeric(t(L) %*% V %*% L)
  
  L_int  <- Lvec(setNames(list(1, -1), c(I_COL, I_ER)))
  L_main <- Lvec(setNames(list(1, -1), c(T_COL, T_ER)))
  out$contrast <- data.frame(
    sample_id = sid,
    contrast  = c("interaction_diff", "main_diff"),
    est       = c(est_of(L_int),  est_of(L_main)),
    var       = c(var_of(L_int),  var_of(L_main)),
    n_spots   = nrow(d),
    n_events  = sum(d$y),
    row.names = NULL)
  
  out$curve <- do.call(rbind, lapply(JOINT_Z_GRID, function(z) {
    Lc <- Lvec(setNames(list(1, z),          c(T_COL, I_COL)))
    Le <- Lvec(setNames(list(1, z),          c(T_ER,  I_ER)))
    Ld <- Lvec(setNames(list(1, -1, z, -z),  c(T_COL, T_ER, I_COL, I_ER)))
    data.frame(sample_id = sid, z = z,
               series = c(LAB_COL, LAB_ER, LAB_DIF),
               est    = c(est_of(Lc), est_of(Le), est_of(Ld)),
               var    = c(var_of(Lc), var_of(Le), var_of(Ld)),
               row.names = NULL)
  }))
  
  out
}

## ============================================================================
## 13.3  RUN ACROSS THE COHORT
## ============================================================================
log_msg("Section 13 start. ERstress = %s. ADJUST_DEPTH = %s. run_id = %s.",
        JOINT_ES, JOINT_ADJUST_DEPTH, RUN_ID)

joint_acc <- list(coef = list(), contrast = list(), curve = list(),
                  fingerprint = list(), modelcheck = list())
for (i in seq_len(nrow(work))) {
  path <- work$path[i]; sid <- work$sample_id[i]
  if (is.na(path) || !file.exists(path)) next
  log_msg("[joint %d/%d] %s", i, nrow(work), sid)
  d <- tryCatch(joint_expr_one(path, sid),
                error = function(e) { log_msg("  joint skipped: %s", conditionMessage(e)); NULL })
  if (!is.null(d) && !is.null(attr(d, "fingerprint")))
    joint_acc$fingerprint[[length(joint_acc$fingerprint) + 1]] <- attr(d, "fingerprint")
  r <- fit_joint_one(d, sid)
  for (k in c("coef", "contrast", "curve", "modelcheck"))
    if (!is.null(r[[k]]) && nrow(r[[k]])) joint_acc[[k]][[length(joint_acc[[k]]) + 1]] <- r[[k]]
}

jcoef_all <- bind(joint_acc$coef)
jcon_all  <- bind(joint_acc$contrast)
jcur_all  <- bind(joint_acc$curve)
jfp_all   <- bind(joint_acc$fingerprint)
jmc_all   <- bind(joint_acc$modelcheck)

if (nrow(jcon_all) == 0) stop("Section 13: no section produced a usable joint fit.")

cat("\n--- model verification (reads the fitted models, not the flag) ---\n")
cat(sprintf("  sections fitted           : %d\n", nrow(jmc_all)))
cat(sprintf("  models containing zDepth  : %d\n", sum(jmc_all$has_zDepth)))
cat(sprintf("  formula used              : %s\n", unique(jmc_all$formula)[1]))
if (length(unique(jmc_all$formula)) > 1)
  stop("Section 13: more than one formula was fitted across sections. Aborting.")
cat(sprintf("  -> this run IS %s\n",
            if (all(jmc_all$has_zDepth)) "DEPTH ADJUSTED" else "UNADJUSTED"))

write_out(jcoef_all, "JOINT_coefficients_per_sample")
write_out(jcon_all,  "JOINT_contrasts_per_sample")
write_out(jfp_all,   "JOINT_module_fingerprint")
write_out(jmc_all,   "JOINT_model_check")

k_used <- dplyr::n_distinct(jcon_all$sample_id)
log_msg("Section 13: %d sections (%d patients) contributed.", k_used,
        length(unique(patient_of(unique(jcon_all$sample_id)))))
if (k_used < nrow(work)) {
  miss <- setdiff(work$sample_id, unique(jcon_all$sample_id))
  cat(sprintf("\n  sections excluded from section 13 (%d): %s\n", length(miss),
              paste(miss, collapse = ", ")))
}

## ============================================================================
## 13.4  POOL THE COEFFICIENTS AND THE SCALAR CONTRASTS
## ============================================================================
jcoef_pooled <- jcoef_all %>%
  group_by(term) %>% group_modify(~ pool_with_p(.x$log_OR, .x$SE^2, .x$sample_id)) %>% ungroup() %>%
  mutate(OR = exp(est), lower_OR = exp(lower), upper_OR = exp(upper))
write_out(jcoef_pooled, "JOINT_coefficients_pooled")

jcon_pooled <- jcon_all %>%
  group_by(contrast) %>% group_modify(~ pool_with_p(.x$est, .x$var, .x$sample_id)) %>% ungroup() %>%
  mutate(ratio_of_ORs = exp(est), lower_ratio = exp(lower), upper_ratio = exp(upper))
write_out(jcon_pooled, "JOINT_slope_difference_pooled")

## ============================================================================
## 13.5  POOL THE CURVES, one meta-analysis per grid point per series
## ============================================================================
jcur_pooled <- jcur_all %>%
  group_by(z, series) %>% group_modify(~ pool_with_p(.x$est, .x$var, .x$sample_id)) %>% ungroup() %>%
  arrange(series, z) %>%
  mutate(OR = exp(est), lower_OR = exp(lower), upper_OR = exp(upper))
write_out(jcur_pooled, "JOINT_slope_curves_pooled")

jdif <- jcur_pooled[jcur_pooled$series == LAB_DIF, ]
jdif$sig <- is.finite(jdif$p) & jdif$p < JOINT_ALPHA
jlines <- jcur_pooled[jcur_pooled$series != LAB_DIF, ]
jlines$series <- factor(jlines$series, levels = c(LAB_COL, LAB_ER))

cross_z <- NA_real_; cross_sig <- NA
{
  s  <- sign(jdif$est)
  ix <- which(diff(s) != 0)
  if (length(ix)) {
    i1 <- ix[1]
    x1 <- jdif$z[i1];   x2 <- jdif$z[i1 + 1]
    y1 <- jdif$est[i1]; y2 <- jdif$est[i1 + 1]
    cross_z   <- if (y2 != y1) x1 - y1 * (x2 - x1) / (y2 - y1) else x1
    cross_sig <- jdif$sig[which.min(abs(jdif$z - cross_z))]
  }
}

sig_lo <- sig_hi <- NA_real_
if (any(jdif$sig)) {
  r      <- rle(jdif$sig)
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  w      <- which(r$values)
  bigg   <- w[which.max(r$lengths[w])]
  sig_lo <- jdif$z[starts[bigg]]
  sig_hi <- jdif$z[ends[bigg]]
}

## ============================================================================
## 13.6  CONSOLE REPORT
## ============================================================================
ri  <- jcon_pooled[jcon_pooled$contrast == "interaction_diff", ]
rm_ <- jcon_pooled[jcon_pooled$contrast == "main_diff", ]

cat("\n", strrep("=", 74), "\n", sep = "")
cat(sprintf("SECTION 13  [%s]  run_id %s\n",
            if (JOINT_ADJUST_DEPTH) "DEPTH ADJUSTED" else "UNADJUSTED", RUN_ID))
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("Model : y ~ %s\n", unique(jmc_all$formula)[1]))
cat(sprintf("ER stress module : %s   |   sections fitted : %d   |   pooled over %s\n", JOINT_ES, k_used, K_UNIT))
cat("All predictors z scored WITHIN section, so every effect is per 1 SD.\n")

cat("\n--- pooled coefficients from the joint model ---\n")
for (i in seq_len(nrow(jcoef_pooled))) with(jcoef_pooled[i, ],
                                            cat(sprintf("  %-16s beta %+7.4f  OR %5.2f (%.2f-%.2f)  k=%2d  I2=%3.0f%%  p=%s\n",
                                                        term, est, OR, lower_OR, upper_OR, as.integer(k), I2, fmt_p(p))))

cat("\n--- do the collagen and ER stress slopes respond differently to hypoxia? ---\n")
if (nrow(ri) == 1) {
  cat("  (collagen x hypoxia) minus (ER stress x hypoxia):\n")
  cat(sprintf("     estimate %+.4f   95%% CI %+.4f to %+.4f   %s   I2=%.0f%%   p = %s\n",
              ri$est, ri$lower, ri$upper, k_lab(ri), ri$I2, fmt_p(ri$p)))
  if (is.finite(ri$p) && ri$p < JOINT_ALPHA) {
    cat(sprintf("     -> the two lines DO respond to hypoxia at different rates (%s).\n",
                if (ri$est < 0) "collagen falls relative to ER stress as hypoxia rises"
                else            "collagen rises relative to ER stress as hypoxia rises"))
  } else {
    cat("     -> no evidence the two lines respond to hypoxia differently.\n")
  }
}
if (nrow(rm_) == 1)
  cat(sprintf("\n  Gap at average hypoxia (z = 0): %+.4f (%+.4f to %+.4f) p = %s\n",
              rm_$est, rm_$lower, rm_$upper, fmt_p(rm_$p)))

cat("\n--- where are the two lines actually distinguishable? ---\n")
for (zt in c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)) {
  r <- jdif[which.min(abs(jdif$z - zt)), ]
  cat(sprintf("   %+5.1f   %+7.4f  (%+7.4f to %+7.4f)   %s%s\n",
              r$z, r$est, r$lower, r$upper, fmt_p(r$p), if (r$sig) "  *" else ""))
}
if (is.finite(sig_lo)) {
  cat(sprintf("\n  Lines separable from z = %+.2f to z = %+.2f SD (%d of %d grid points).\n",
              sig_lo, sig_hi, sum(jdif$sig), nrow(jdif)))
} else {
  cat("\n  The two lines are not separable anywhere on the grid.\n")
}

cat("\n--- HOW TO SAY THIS OUT LOUD ---\n")
if (is.finite(cross_z)) {
  cat(sprintf("  The two lines meet at z_Hyp = %+.2f SD, and at that point the\n", cross_z))
  if (isTRUE(cross_sig)) {
    cat("  difference IS significant, so they do swap order.\n")
  } else {
    cat("  difference is NOT significant.\n")
    cat("  Supported     : the two effects converge as hypoxia changes.\n")
    cat("  NOT supported : that they swap order.\n")
  }
} else {
  cat("  The two lines do not meet anywhere within +/- 2 SD of hypoxia.\n")
}

cat("\nCAVEATS\n")
cat("  1. Spatial autocorrelation is not modelled. Every CI here is too narrow.\n")
cat("     Sections 16 and 26 measure how large that is. Quote empirical percentiles.\n")
cat("  2. Per 1 SD within section. Not comparable to sections 10 and 11.\n")
cat(sprintf("  3. Depth adjustment is %s (verified against the fitted models).\n",
            if (JOINT_ADJUST_DEPTH) "ON" else "OFF"))
cat("  4. Associational, mRNA level.\n")
cat(strrep("=", 74), "\n", sep = "")

## ============================================================================
## 13.7  FIGURE 4, REBUILT
##       v5 (C2): title and annotations are computed from this run.
## ============================================================================
if (!exists("FIG_STYLES") || !exists("theme_fig")) {
  cat("\n  Figure skipped: section 12 has not run.\n")
} else if (nrow(ri) == 1) {
  stat_lab <- sprintf("Difference between the two slopes  %+.3f  (95%% CI %+.3f to %+.3f)\np = %s      %s      I2 = %.0f%%%s",
                      ri$est, ri$lower, ri$upper, fmt_p(ri$p), k_lab(ri), ri$I2,
                      if (JOINT_ADJUST_DEPTH) "\nsequencing depth adjusted" else "")
  fig4_title <- if (is.finite(ri$p) && ri$p < JOINT_ALPHA)
    sprintf("The collagen effect %s relative to ER stress as hypoxia rises",
            if (ri$est < 0) "weakens" else "strengthens")
  else "No evidence that the collagen and ER stress effects respond differently to hypoxia"
  
  y_top <- max(jlines$upper_OR, na.rm = TRUE)
  y_bot <- min(jlines$lower_OR, na.rm = TRUE)
  y_hi  <- y_top * (y_top / y_bot)^0.30
  
  for (sn in names(FIG_STYLES)) {
    S <- FIG_STYLES[[sn]]
    p4 <- ggplot(jlines, aes(z, OR))
    
    if (is.finite(sig_lo)) {
      p4 <- p4 +
        annotate("rect", xmin = sig_lo, xmax = sig_hi,
                 ymin = -Inf, ymax = Inf, fill = "grey40", alpha = 0.14) +
        annotate("text", x = (sig_lo + sig_hi) / 2, y = y_top, vjust = -0.4,
                 size = S$lab * 0.95, colour = "grey40", lineheight = 1.1,
                 label = sprintf("effects differ (p < %.2f)\n%+.1f to %+.1f SD", JOINT_ALPHA, sig_lo, sig_hi))
    }
    
    p4 <- p4 +
      geom_hline(yintercept = 1, linetype = "22", linewidth = S$lwd * 0.7, colour = "grey40") +
      geom_ribbon(aes(ymin = lower_OR, ymax = upper_OR, fill = series),
                  alpha = .16, colour = NA) +
      geom_line(aes(colour = series), linewidth = S$lwd * 1.4) +
      scale_colour_manual(values = setNames(c(CLR$neg, CLR$pos), c(LAB_COL, LAB_ER)), name = NULL) +
      scale_fill_manual(  values = setNames(c(CLR$neg, CLR$pos), c(LAB_COL, LAB_ER)), name = NULL)
    
    if (is.finite(cross_z)) {
      p4 <- p4 +
        geom_vline(xintercept = cross_z, linetype = "dotted",
                   linewidth = S$lwd * 0.7, colour = "grey35") +
        annotate("text", x = cross_z, y = y_bot, hjust = 1.06, vjust = -0.25,
                 size = S$lab * 0.95, colour = "grey35", lineheight = 1.1,
                 label = if (isTRUE(cross_sig))
                   sprintf("lines meet at %+.1f SD", cross_z)
                 else
                   sprintf("lines meet at %+.1f SD\n(difference not significant here)", cross_z))
    }
    
    p4 <- p4 +
      annotate("text", x = min(jlines$z), y = y_hi, label = stat_lab,
               hjust = 0, vjust = 1, size = S$lab * 1.1, colour = "grey15",
               lineheight = 1.15) +
      scale_y_log10() +
      coord_cartesian(xlim = range(jlines$z), ylim = c(y_bot, y_hi)) +
      labs(x = "Hypoxia score (SD within section)", y = XLAB_OR,
           caption = paste0("Both slopes from one joint model per section; all covariances included. ",
                            if (is.finite(sig_lo))
                              sprintf("Shaded: hypoxia range where the two effects differ (p < %.2f).", JOINT_ALPHA)
                            else "")) +
      theme_fig(S) + theme(legend.position = "top")
    
    if (isTRUE(S$titles))
      p4 <- p4 + labs(
        title = fig4_title,
        subtitle = paste0("Effect of each module on P4HA1 detection, as a function of hypoxia",
                          if (JOINT_ADJUST_DEPTH) " (sequencing depth adjusted)" else ""))
    
    save_fig(p4, S, sub("\\.csv$", "", outf("Figure4_joint_effect_modification")),
             if (isTRUE(S$split)) S$w2 else S$w1 * 1.5,
             if (isTRUE(S$split)) S$h_full else MM(68))
  }
}

log_msg("Section 13 complete.")

################################################################################
## 15. THE SPECIFICITY TEST: P4HA1 vs P4HA2 vs P4HA3
##
##   v5.2: the old premise ("P4HA3 does not hydroxylate collagen") was wrong.
##   All three genes encode collagen P4H alpha subunits, and P4HA3 forms an
##   active collagen hydroxylase (Kukkola et al., J Biol Chem 2003,
##   doi:10.1074/jbc.M306806200). This section therefore asks only whether
##   the pattern is specific to P4HA1 among the three alpha subunits; it is
##   NOT a test of enzyme function versus gene family.
##   The model does not change between outcomes: identical predictors, spots,
##   standardisation and pooling; only y changes.
################################################################################

stopifnot(exists("pool_with_p"), exists("fmt_p"))

## ============================================================================
## 15.0  SETTINGS
## ============================================================================
ISO_OUTCOMES     <- c("P4HA1", "P4HA2", "P4HA3")
ISO_ES           <- "ERstressStrict"
ISO_MIN_SPOTS    <- 30
ISO_MIN_EVENTS   <- 70
ISO_ADJUST_DEPTH <- ADJUST_DEPTH    # derived, never set by hand
ISO_Z_GRID       <- seq(-2, 2, length.out = 41)
ISO_ALPHA        <- 0.05

ISO_SCORE_MODULES <- c("Hypoxia", "CollagenSynthesis", "iCAF", ISO_ES, "ERAD")
ISO_RENAME <- c(Hypoxia = "zHyp", CollagenSynthesis = "zCol",
                iCAF = "ziCAF", ERAD = "zERAD")
ISO_RENAME[ISO_ES] <- "zER"

IT_COL <- "zCol"; IT_ER <- "zER"; IT_HYP <- "zHyp"
II_COL <- paste0(IT_HYP, ":", IT_COL)
II_ER  <- paste0(IT_HYP, ":", IT_ER)

assert_run_state("Section 15")

## ============================================================================
## 15.1  CIRCULARITY GATE  -- hard stop
## ============================================================================
cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 15.1  CIRCULARITY GATE\n")
cat(strrep("=", 74), "\n", sep = "")
iso_fail <- FALSE
for (out_gene in ISO_OUTCOMES) {
  bad <- ISO_SCORE_MODULES[vapply(ISO_SCORE_MODULES,
                                  function(m) out_gene %in% MODULES[[m]], logical(1))]
  cat(sprintf("  outcome %-6s : %s\n", out_gene,
              if (length(bad)) paste0("*** CIRCULAR, in ", paste(bad, collapse = ", "), " ***")
              else "clean, no predictor module contains it"))
  if (length(bad)) iso_fail <- TRUE
}
if (iso_fail) stop("Section 15 halted: an outcome gene appears in a predictor module.")

## ============================================================================
## 15.2  SCORE AND STANDARDISE ONCE PER SECTION (reused for all outcomes)
## ============================================================================
iso_expr_one <- function(rds_path, sid) {
  obj <- readRDS(rds_path)
  obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
  DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else Assays(obj)[1]
  
  missing_out <- setdiff(ISO_OUTCOMES, rownames(obj))
  if (length(missing_out)) {
    log_msg("  %s: outcome gene(s) absent (%s) -> section skipped",
            sid, paste(missing_out, collapse = ", "))
    return(NULL)
  }
  
  fp <- list()
  for (nm in ISO_SCORE_MODULES) {
    full <- MODULES[[nm]]
    g    <- intersect(full, rownames(obj))
    need <- min_genes_for(full)
    if (length(g) < need) {
      log_msg("  %s: %s has only %d/%d genes (need %d) -> section skipped",
              sid, nm, length(g), length(full), need)
      return(NULL)
    }
    set.seed(SEED)
    obj <- tryCatch(add_module_score_fallback(obj, features = list(g), name = nm,
                                              assay = DefaultAssay(obj), nbin = 20, ctrl = 50),
                    error = function(e) NULL)
    if (is.null(obj) || !paste0(nm, "1") %in% colnames(obj@meta.data)) return(NULL)
    names(obj@meta.data)[names(obj@meta.data) == paste0(nm, "1")] <- nm
    v <- obj@meta.data[[nm]]
    fp[[nm]] <- data.frame(sample_id = sid, module = nm,
                           n_genes_in_set = length(full), n_genes_used = length(g),
                           score_mean = mean(v, na.rm = TRUE),
                           score_sd   = sd(v,   na.rm = TRUE),
                           row.names = NULL)
  }
  
  vars <- c(ISO_OUTCOMES, ISO_SCORE_MODULES)
  cnt_col <- grep("^nCount_", colnames(obj@meta.data), value = TRUE)
  if (ISO_ADJUST_DEPTH) {
    if (!length(cnt_col)) return(NULL)
    vars <- c(vars, cnt_col[1])
  }
  d <- FetchData(obj, vars = vars, layer = "data")
  
  zz <- function(x) { s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) return(NULL)
  as.numeric((x - mean(x, na.rm = TRUE)) / s) }
  for (nm in ISO_SCORE_MODULES) {
    v <- zz(d[[nm]])
    if (is.null(v)) { log_msg("  %s: %s zero variance -> section skipped", sid, nm); return(NULL) }
    d[[ ISO_RENAME[[nm]] ]] <- v
  }
  if (ISO_ADJUST_DEPTH) {
    v <- zz(log1p(d[[cnt_col[1]]])); if (is.null(v)) return(NULL); d$zDepth <- v
  }
  attr(d, "fingerprint") <- do.call(rbind, fp)
  d
}

## ============================================================================
## 15.3  FIT THE SAME MODEL FOR ONE OUTCOME
## ============================================================================
fit_iso_one <- function(d, sid, out_gene) {
  out <- list(coef = NULL, contrast = NULL, curve = NULL, status = NULL)
  
  y  <- as.integer(d[[out_gene]] > 0)
  ev <- sum(y); nev <- length(y) - ev
  status <- data.frame(sample_id = sid, outcome = out_gene, n_spots = length(y),
                       n_events = ev, n_nonevents = nev, fitted = FALSE,
                       has_zDepth = NA, reason = NA_character_, row.names = NULL)
  
  if (length(y) < ISO_MIN_SPOTS || ev == 0 || nev == 0) {
    status$reason <- "no variation in outcome"; out$status <- status; return(out)
  }
  if (ev < ISO_MIN_EVENTS || nev < ISO_MIN_EVENTS) {
    status$reason <- sprintf("too few events (%d) or non-events (%d), threshold %d",
                             ev, nev, ISO_MIN_EVENTS)
    out$status <- status; return(out)
  }
  
  dd <- d; dd$y <- y
  rhs <- paste0(IT_HYP, " * ", IT_COL, " + ", IT_HYP, " * ", IT_ER, " + ziCAF + zERAD",
                if (ISO_ADJUST_DEPTH) " + zDepth" else "")
  f <- tryCatch(glm(as.formula(paste("y ~", rhs)), data = dd, family = binomial()),
                error = function(e) NULL)
  if (is.null(f) || !f$converged) {
    status$reason <- "glm did not converge"; out$status <- status; return(out)
  }
  
  b <- coef(f); V <- vcov(f)
  
  has_depth <- "zDepth" %in% names(b)
  if (!identical(has_depth, isTRUE(ISO_ADJUST_DEPTH)))
    stop(sprintf("Section 15, %s / %s: ISO_ADJUST_DEPTH = %s but the fitted model %s zDepth.",
                 sid, out_gene, ISO_ADJUST_DEPTH,
                 if (has_depth) "DOES contain" else "does NOT contain"))
  status$has_zDepth <- has_depth
  
  need <- c(IT_COL, IT_ER, II_COL, II_ER)
  if (!all(need %in% names(b)) || any(!is.finite(b[need]))) {
    status$reason <- "coefficients not estimable"; out$status <- status; return(out)
  }
  status$fitted <- TRUE; out$status <- status
  
  s  <- summary(f)$coefficients
  tm <- setdiff(rownames(s), "(Intercept)")
  out$coef <- data.frame(sample_id = sid, outcome = out_gene, term = tm,
                         log_OR = s[tm, "Estimate"], SE = s[tm, "Std. Error"],
                         p = s[tm, "Pr(>|z|)"], row.names = NULL)
  
  Lvec <- function(w) { L <- setNames(numeric(length(b)), names(b))
  for (nm in names(w)) L[nm] <- w[[nm]]; L }
  est_of <- function(L) as.numeric(L %*% b)
  var_of <- function(L) as.numeric(t(L) %*% V %*% L)
  
  L_int  <- Lvec(setNames(list(1, -1), c(II_COL, II_ER)))
  L_main <- Lvec(setNames(list(1, -1), c(IT_COL, IT_ER)))
  out$contrast <- data.frame(
    sample_id = sid, outcome = out_gene,
    contrast = c("interaction_diff", "main_diff"),
    est = c(est_of(L_int), est_of(L_main)),
    var = c(var_of(L_int), var_of(L_main)),
    n_spots = nrow(dd), n_events = ev, row.names = NULL)
  
  out$curve <- do.call(rbind, lapply(ISO_Z_GRID, function(z) {
    Lc <- Lvec(setNames(list(1, z),         c(IT_COL, II_COL)))
    Le <- Lvec(setNames(list(1, z),         c(IT_ER,  II_ER)))
    Ld <- Lvec(setNames(list(1, -1, z, -z), c(IT_COL, IT_ER, II_COL, II_ER)))
    data.frame(sample_id = sid, outcome = out_gene, z = z,
               series = c("Collagen synthesis", "ER stress", "difference"),
               est = c(est_of(Lc), est_of(Le), est_of(Ld)),
               var = c(var_of(Lc), var_of(Le), var_of(Ld)), row.names = NULL)
  }))
  out
}

## ============================================================================
## 15.4  RUN
## ============================================================================
log_msg("Section 15 start. ADJUST_DEPTH = %s. run_id = %s.", ISO_ADJUST_DEPTH, RUN_ID)

iso_acc <- list(coef = list(), contrast = list(), curve = list(),
                status = list(), fingerprint = list())
for (i in seq_len(nrow(work))) {
  path <- work$path[i]; sid <- work$sample_id[i]
  if (is.na(path) || !file.exists(path)) next
  log_msg("[iso %d/%d] %s", i, nrow(work), sid)
  d <- tryCatch(iso_expr_one(path, sid),
                error = function(e) { log_msg("  iso skipped: %s", conditionMessage(e)); NULL })
  if (is.null(d)) next
  if (!is.null(attr(d, "fingerprint")))
    iso_acc$fingerprint[[length(iso_acc$fingerprint) + 1]] <- attr(d, "fingerprint")
  for (og in ISO_OUTCOMES) {
    r <- fit_iso_one(d, sid, og)
    for (k in c("coef", "contrast", "curve", "status"))
      if (!is.null(r[[k]]) && nrow(r[[k]])) iso_acc[[k]][[length(iso_acc[[k]]) + 1]] <- r[[k]]
  }
}

iso_coef <- bind(iso_acc$coef)
iso_con  <- bind(iso_acc$contrast)
iso_cur  <- bind(iso_acc$curve)
iso_stat <- bind(iso_acc$status)
iso_fp   <- bind(iso_acc$fingerprint)

if (nrow(iso_con) == 0) stop("Section 15: no outcome produced a usable fit.")

cat("\n--- model verification ---\n")
cat(sprintf("  fits with zDepth: %d of %d fitted models\n",
            sum(iso_stat$has_zDepth, na.rm = TRUE), sum(iso_stat$fitted)))
cat(sprintf("  -> this run IS %s\n",
            if (all(iso_stat$has_zDepth[iso_stat$fitted])) "DEPTH ADJUSTED" else "UNADJUSTED"))

write_out(iso_coef, "ISOFORM_coefficients_per_sample")
write_out(iso_con,  "ISOFORM_contrasts_per_sample")
write_out(iso_stat, "ISOFORM_fit_status")
write_out(iso_fp,   "ISOFORM_module_fingerprint")

## ============================================================================
## 15.5  POOL
## ============================================================================
iso_coef_pooled <- iso_coef %>%
  group_by(outcome, term) %>% group_modify(~ pool_with_p(.x$log_OR, .x$SE^2, .x$sample_id)) %>% ungroup() %>%
  mutate(OR = exp(est), lower_OR = exp(lower), upper_OR = exp(upper))
write_out(iso_coef_pooled, "ISOFORM_coefficients_pooled")

iso_con_pooled <- iso_con %>%
  group_by(outcome, contrast) %>% group_modify(~ pool_with_p(.x$est, .x$var, .x$sample_id)) %>% ungroup()
write_out(iso_con_pooled, "ISOFORM_slope_difference_pooled")

iso_cur_pooled <- iso_cur %>%
  group_by(outcome, series, z) %>% group_modify(~ pool_with_p(.x$est, .x$var, .x$sample_id)) %>% ungroup() %>%
  mutate(OR = exp(est), lower_OR = exp(lower), upper_OR = exp(upper)) %>%
  arrange(outcome, series, z)
write_out(iso_cur_pooled, "ISOFORM_slope_curves_pooled")

## ============================================================================
## 15.6  REPORT
## ============================================================================
cat("\n", strrep("=", 74), "\n", sep = "")
cat(sprintf("SECTION 15  [%s]  run_id %s\n",
            if (ISO_ADJUST_DEPTH) "DEPTH ADJUSTED" else "UNADJUSTED", RUN_ID))
cat(strrep("=", 74), "\n", sep = "")

cat("\n--- which sections contributed, per outcome ---\n")
for (og in ISO_OUTCOMES) {
  s <- iso_stat[iso_stat$outcome == og, ]
  fs <- s$sample_id[s$fitted]
  cat(sprintf("  %-6s fitted in %2d / %2d sections (%d patients)   median events %5.0f\n",
              og, sum(s$fitted), nrow(s), length(unique(patient_of(fs))),
              if (length(fs)) median(s$n_events[s$fitted]) else NA_real_))
  drop <- s[!s$fitted, ]
  if (nrow(drop))
    for (j in seq_len(nrow(drop)))
      cat(sprintf("         excluded: %-30s %s\n", drop$sample_id[j], drop$reason[j]))
}
cat("  Name these sections whenever you quote k.\n")

cat("\n--- the two interaction terms, side by side ---\n")
cat(sprintf("  %-18s %-8s %9s %19s %6s %6s %11s\n",
            "term", "outcome", "beta", "95% CI", "k", "I2", "p"))
for (tm in c(II_COL, II_ER)) {
  for (og in ISO_OUTCOMES) {
    r <- iso_coef_pooled[iso_coef_pooled$term == tm & iso_coef_pooled$outcome == og, ]
    if (!nrow(r)) next
    cat(sprintf("  %-18s %-8s %+9.4f  %+8.4f to %+8.4f %5d %5.0f%% %10s\n",
                tm, og, r$est, r$lower, r$upper, as.integer(r$k), r$I2, fmt_p(r$p)))
  }
  cat("\n")
}

cat("--- THE HEADLINE: difference between the two slopes, per outcome ---\n")
for (og in ISO_OUTCOMES) {
  r <- iso_con_pooled[iso_con_pooled$outcome == og &
                        iso_con_pooled$contrast == "interaction_diff", ]
  if (!nrow(r)) next
  cat(sprintf("  %-6s %+.4f  (95%% CI %+.4f to %+.4f)  %s  I2=%3.0f%%  p = %-10s %s\n",
              og, r$est, r$lower, r$upper, k_lab(r), r$I2, fmt_p(r$p),
              if (is.finite(r$p) && r$p < ISO_ALPHA) "SIGNIFICANT" else "not significant"))
}

## v5 (C2): the P4HA3 statement is decided by THIS run's numbers, including
## whether the P4HA3 interval is narrow enough to rule out the P4HA1 effect.
{
  a3 <- iso_con_pooled[iso_con_pooled$outcome == "P4HA3" &
                         iso_con_pooled$contrast == "interaction_diff", ]
  a1 <- iso_con_pooled[iso_con_pooled$outcome == "P4HA1" &
                         iso_con_pooled$contrast == "interaction_diff", ]
  cat("\n--- WHAT YOU MAY AND MAY NOT CLAIM ABOUT P4HA3 ---\n")
  if (nrow(a3) != 1 || nrow(a1) != 1) {
    cat("  P4HA3 or P4HA1 could not be pooled; no claim is possible.\n")
  } else {
    sig1 <- is.finite(a1$p) && a1$p < ISO_ALPHA
    sig3 <- is.finite(a3$p) && a3$p < ISO_ALPHA
    same_sign <- sign(a3$est) == sign(a1$est)
    a3_covers_a1 <- is.finite(a3$lower) && a3$lower <= a1$est && a1$est <= a3$upper
    cat(sprintf("  P4HA1 %+.4f (%s, p = %s)   P4HA3 %+.4f (95%% CI %+.4f to %+.4f, %s, p = %s)\n",
                a1$est, k_lab(a1), fmt_p(a1$p), a3$est, a3$lower, a3$upper, k_lab(a3), fmt_p(a3$p)))
    if (!sig1) {
      cat("  P4HA1 itself shows no significant slope difference, so there is no\n")
      cat("  pattern for P4HA3 to confirm or deny. No specificity claim is possible.\n")
    } else if (!sig3 && !a3_covers_a1) {
      cat("  P4HA3 shows no interaction AND its interval excludes the P4HA1 effect\n")
      cat("  size. That is a real negative, not just a lack of power: the\n")
      cat("  specificity reading is supported.\n")
    } else if (!sig3) {
      cat("  P4HA3 is not significant, BUT its interval still contains the P4HA1\n")
      cat("  effect size. This is a lack of power, not evidence of no effect.\n")
      cat("  Do NOT describe P4HA3 as 'showing no effect'.\n")
    } else if (same_sign) {
      cat("  *** P4HA3 shows the SAME pattern and it IS significant. The pattern is\n")
      cat("  then shared by the collagen P4H alpha subunits, not specific to P4HA1. ***\n")
    } else {
      cat("  *** P4HA3 is SIGNIFICANT but in the OPPOSITE direction to P4HA1. This\n")
      cat("  is neither the clean negative the specificity argument assumes nor\n")
      cat("  co-regulation. It needs its own explanation before anything is claimed. ***\n")
    }
  }
}

cat("\nCAVEATS\n")
cat("  Spatial autocorrelation not modelled; all CIs too narrow.\n")
cat("  Per 1 SD within section; not comparable to sections 10 and 11.\n")
cat("  Coefficient magnitudes are not comparable across genes with different\n")
cat("  detection rates. Presence or absence of an interaction is.\n")
cat("  No multiple-comparison correction applied.\n")
cat(strrep("=", 74), "\n", sep = "")

## ============================================================================
## 15.7  FIGURE 7
## ============================================================================
if (exists("FIG_STYLES") && exists("theme_fig")) {
  lines_df <- iso_cur_pooled[iso_cur_pooled$series != "difference", ]
  lines_df$series  <- factor(lines_df$series, c("Collagen synthesis", "ER stress"))
  lines_df$outcome <- factor(lines_df$outcome, ISO_OUTCOMES)
  
  lab_df <- do.call(rbind, lapply(ISO_OUTCOMES, function(og) {
    r <- iso_con_pooled[iso_con_pooled$outcome == og &
                          iso_con_pooled$contrast == "interaction_diff", ]
    if (!nrow(r)) return(NULL)
    data.frame(outcome = factor(og, ISO_OUTCOMES),
               lab = sprintf("diff %+.3f, p = %s\n%s", r$est, fmt_p(r$p), k_lab(r)),
               row.names = NULL)
  }))
  
  y_top <- max(lines_df$upper_OR, na.rm = TRUE)
  y_bot <- min(lines_df$lower_OR, na.rm = TRUE)
  y_hi  <- y_top * (y_top / y_bot)^0.28
  
  for (sn in names(FIG_STYLES)) {
    S <- FIG_STYLES[[sn]]
    p7 <- ggplot(lines_df, aes(z, OR)) +
      geom_hline(yintercept = 1, linetype = "22", linewidth = S$lwd * 0.7, colour = "grey40") +
      geom_ribbon(aes(ymin = lower_OR, ymax = upper_OR, fill = series), alpha = .16, colour = NA) +
      geom_line(aes(colour = series), linewidth = S$lwd * 1.3) +
      scale_colour_manual(values = setNames(c(CLR$neg, CLR$pos),
                                            c("Collagen synthesis", "ER stress")), name = NULL) +
      scale_fill_manual(values   = setNames(c(CLR$neg, CLR$pos),
                                            c("Collagen synthesis", "ER stress")), name = NULL) +
      facet_wrap(~ outcome, nrow = 1) + scale_y_log10() +
      coord_cartesian(ylim = c(y_bot, y_hi)) +
      labs(x = "Hypoxia score (SD within section)",
           y = "Odds ratio for detection (95% CI)",
           caption = paste0("Identical model and predictors; only the outcome gene differs. Each gene is fitted only in sections with >= ",
                            ISO_MIN_EVENTS, " detected spots, so k differs between panels.",
                            if (ISO_ADJUST_DEPTH) " Sequencing depth adjusted." else "")) +
      theme_fig(S) + theme(legend.position = "top")
    
    if (!is.null(lab_df) && nrow(lab_df))
      p7 <- p7 + geom_text(data = lab_df, inherit.aes = FALSE,
                           aes(x = min(lines_df$z), y = y_hi, label = lab),
                           hjust = 0, vjust = 1, size = S$lab, colour = "grey20",
                           lineheight = 1.1)
    
    if (isTRUE(S$titles))
      p7 <- p7 + labs(title = "Is the pattern specific to P4HA1?",
                      subtitle = "The three collagen prolyl 4-hydroxylase alpha subunits, same model")
    
    save_fig(p7, S, sub("\\.csv$", "", outf("Figure7_isoform_comparison")),
             S$w2, if (isTRUE(S$split)) S$h_full else MM(72))
  }
}

log_msg("Section 15 complete.")


################################################################################
## 16. THE RANDOM-GENE CONTROL
##
##   Model p values assume independent spots. A reference built from real
##   genes in the same tissue carries the same spatial inflation, so the
##   comparison is immune to that assumption. Genes are matched to P4HA1 on
##   detection rate.
################################################################################

RATES_CSV <- file.path(agg, "DETECT_cohort_gene_rates.csv")
if (!file.exists(RATES_CSV))
  stop("Section 16 needs DETECT_cohort_gene_rates.csv. Run section 14 first.")

## ============================================================================
## 16.0  SETTINGS
## ============================================================================
RND_N_GENES      <- 1000
RND_SEED         <- 20260831
RND_CALIPER      <- 0.25
RND_MIN_EVENTS   <- 70
RND_MIN_SPOTS    <- 30
RND_MIN_K_FRAC   <- 2 / 3     # share of pooling units a gene must be pooled over
RND_ES           <- "ERstressStrict"
RND_ADJUST_DEPTH <- ADJUST_DEPTH    # derived, never set by hand

RND_SCORE_MODULES <- c("Hypoxia", "CollagenSynthesis", "iCAF", RND_ES, "ERAD")
RND_RENAME <- c(Hypoxia = "zHyp", CollagenSynthesis = "zCol",
                iCAF = "ziCAF", ERAD = "zERAD")
RND_RENAME[RND_ES] <- "zER"

RT_COL <- "zCol"; RT_ER <- "zER"; RT_HYP <- "zHyp"
RI_COL <- paste0(RT_HYP, ":", RT_COL)
RI_ER  <- paste0(RT_HYP, ":", RT_ER)

assert_run_state("Section 16")

## ============================================================================
## 16.1  CHOOSE THE SECTIONS AND THE GENES
##       Built on the sections where P4HA1 itself can be modelled (at least
##       RND_MIN_EVENTS positive and negative spots). A section with P4HA1 in a
##       handful of spots falls under the detection cache threshold, and the
##       "present in every section" filter would then drop P4HA1 itself.
## ============================================================================
rates <- read.csv(RATES_CSV, stringsAsFactors = FALSE)

det16 <- read.csv(file.path(agg, "DETECT_p4ha_per_sample.csv"), stringsAsFactors = FALSE)
det16 <- det16[det16$gene == "P4HA1" & as.logical(det16$present), , drop = FALSE]
RND_SECTIONS <- det16$sample_id[det16$n_detected >= RND_MIN_EVENTS &
                                  (det16$n_spots - det16$n_detected) >= RND_MIN_EVENTS]
RND_SECTIONS <- intersect(RND_SECTIONS, work$sample_id)
rnd_dropped <- setdiff(work$sample_id, RND_SECTIONS)
if (length(rnd_dropped))
  log_msg("Section 16: %d section(s) left out, P4HA1 below %d events: %s",
          length(rnd_dropped), RND_MIN_EVENTS, paste(rnd_dropped, collapse = ", "))
if (length(RND_SECTIONS) < 2) stop("Section 16: fewer than 2 sections where P4HA1 can be modelled.")
rates <- rates[rates$sample_id %in% RND_SECTIONS, , drop = FALSE]

## a gene enters the reference only if pooled over two thirds of the units
RND_MIN_K <- max(2L, as.integer(ceiling(RND_MIN_K_FRAC * n_units(RND_SECTIONS))))

gene_med <- rates %>%
  group_by(gene) %>%
  summarise(median_rate = median(detect_rate), n_tumours = dplyr::n(), .groups = "drop") %>%
  filter(n_tumours == dplyr::n_distinct(rates$sample_id))

target <- gene_med$median_rate[gene_med$gene == "P4HA1"]
if (!length(target)) stop("P4HA1 not in the cached detection rates.")

## every outcome gene and every gene inside a predictor module is excluded
EXCLUDE <- unique(c("P4HA1", "P4HA2", "P4HA3", "P4HB", "PLOD1", "PLOD2", "PLOD3",
                    unlist(MODULES[RND_SCORE_MODULES], use.names = FALSE)))

pool_genes <- gene_med %>%
  filter(!gene %in% EXCLUDE,
         median_rate >= target * (1 - RND_CALIPER),
         median_rate <= target * (1 + RND_CALIPER))

cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 16.1  GENE SELECTION\n")
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  sections used               : %d (%d %s)\n",
            length(RND_SECTIONS), n_units(RND_SECTIONS), K_UNIT))
cat(sprintf("  minimum k per pooled gene   : %d %s\n", RND_MIN_K, K_UNIT))
cat(sprintf("  P4HA1 median detection rate : %.3f\n", target))
cat(sprintf("  caliper                     : %.0f%% relative (%.3f to %.3f)\n",
            100 * RND_CALIPER, target * (1 - RND_CALIPER), target * (1 + RND_CALIPER)))
cat(sprintf("  genes present in all sections: %d\n", nrow(gene_med)))
cat(sprintf("  eligible inside the caliper  : %d\n", nrow(pool_genes)))

if (nrow(pool_genes) < RND_N_GENES) {
  cat(sprintf("  *** only %d eligible. Using all of them. Describe the reference as\n", nrow(pool_genes)))
  cat("      \"all comparable genes\", not \"1000 random genes\". Do NOT widen the\n")
  cat("      caliper to reach 1000: that breaks the matching. ***\n")
  sel <- pool_genes
} else {
  set.seed(RND_SEED)
  sel <- pool_genes[sample(nrow(pool_genes), RND_N_GENES), ]
}
if (!nrow(sel)) stop("Section 16: no gene inside the caliper.")
sel <- sel[order(sel$median_rate), ]
write_out(sel, "RANDOM_gene_selection")

cat(sprintf("\n  selected: %d genes, detection median %.3f (range %.3f to %.3f)\n",
            nrow(sel), median(sel$median_rate), min(sel$median_rate), max(sel$median_rate)))
cat(sprintf("  P4HA1 sits at the %.0fth percentile of the selected set.\n",
            100 * mean(sel$median_rate <= target)))

RND_GENES <- sel$gene

## ============================================================================
## 16.2  FIT, ONE SECTION AT A TIME (model matrix built once per section)
## ============================================================================
rnd_contrast <- function(b, V) {
  need <- c(RT_COL, RT_ER, RI_COL, RI_ER)
  if (!all(need %in% names(b)) || any(!is.finite(b[need]))) return(NULL)
  if (!all(need %in% rownames(V))) return(NULL)
  L <- setNames(numeric(length(b)), names(b))
  L[RI_COL] <- 1; L[RI_ER] <- -1
  Lv <- L[rownames(V)]
  est <- sum(Lv * b[rownames(V)])
  var <- as.numeric(t(Lv) %*% V %*% Lv)
  if (!is.finite(est) || !is.finite(var) || var <= 0) return(NULL)
  c(est = est, var = var)
}

vcov_from_glmfit <- function(fitg, X) {
  p  <- fitg$rank
  Qr <- fitg$qr
  if (is.null(Qr) || p < 1) return(NULL)
  cv <- tryCatch(chol2inv(Qr$qr[1:p, 1:p, drop = FALSE]), error = function(e) NULL)
  if (is.null(cv)) return(NULL)
  nm <- colnames(X)[Qr$pivot[1:p]]
  dimnames(cv) <- list(nm, nm)
  cv                                  # binomial: dispersion is 1
}

rnd_one_tumour <- function(rds_path, sid) {
  obj <- readRDS(rds_path)
  obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
  DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else Assays(obj)[1]
  if (!"P4HA1" %in% rownames(obj)) return(NULL)
  
  for (nm in RND_SCORE_MODULES) {
    full <- MODULES[[nm]]
    g    <- intersect(full, rownames(obj))
    if (length(g) < min_genes_for(full)) {
      log_msg("  %s: %s coverage too low -> section skipped", sid, nm); return(NULL)
    }
    set.seed(SEED)
    obj <- tryCatch(add_module_score_fallback(obj, features = list(g), name = nm,
                                              assay = DefaultAssay(obj), nbin = 20, ctrl = 50),
                    error = function(e) NULL)
    if (is.null(obj) || !paste0(nm, "1") %in% colnames(obj@meta.data)) return(NULL)
    names(obj@meta.data)[names(obj@meta.data) == paste0(nm, "1")] <- nm
  }
  
  md <- obj@meta.data
  zz <- function(x) { s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) return(NULL)
  as.numeric((x - mean(x, na.rm = TRUE)) / s) }
  Z <- list()
  for (nm in RND_SCORE_MODULES) {
    v <- zz(md[[nm]]); if (is.null(v)) return(NULL)
    Z[[ RND_RENAME[[nm]] ]] <- v
  }
  dz <- as.data.frame(Z)
  if (RND_ADJUST_DEPTH) {
    cc <- grep("^nCount_", colnames(md), value = TRUE)
    if (!length(cc)) return(NULL)
    v <- zz(log1p(md[[cc[1]]])); if (is.null(v)) return(NULL)
    dz$zDepth <- v
  }
  if (nrow(dz) < RND_MIN_SPOTS) return(NULL)
  
  rhs <- paste0("~ ", RT_HYP, " * ", RT_COL, " + ", RT_HYP, " * ", RT_ER,
                " + ziCAF + zERAD", if (RND_ADJUST_DEPTH) " + zDepth" else "")
  X <- model.matrix(as.formula(rhs), data = dz)
  
  has_depth <- "zDepth" %in% colnames(X)
  if (!identical(has_depth, isTRUE(RND_ADJUST_DEPTH)))
    stop(sprintf("Section 16, %s: RND_ADJUST_DEPTH = %s but the design matrix %s zDepth.",
                 sid, RND_ADJUST_DEPTH, if (has_depth) "DOES contain" else "does NOT contain"))
  
  genes_here <- intersect(RND_GENES, rownames(obj))
  if (!length(genes_here)) return(NULL)
  E <- GetAssayData(obj, assay = DefaultAssay(obj), layer = "data")[genes_here, , drop = FALSE]
  
  fam <- binomial()
  out <- vector("list", length(genes_here))
  for (j in seq_along(genes_here)) {
    y <- as.integer(E[j, ] > 0)
    ev <- sum(y); nev <- length(y) - ev
    if (ev < RND_MIN_EVENTS || nev < RND_MIN_EVENTS) next
    fitg <- tryCatch(suppressWarnings(glm.fit(X, y, family = fam)), error = function(e) NULL)
    if (is.null(fitg) || !isTRUE(fitg$converged)) next
    V <- vcov_from_glmfit(fitg, X); if (is.null(V)) next
    b <- fitg$coefficients; names(b) <- colnames(X)
    cc <- rnd_contrast(b, V); if (is.null(cc)) next
    out[[j]] <- data.frame(sample_id = sid, gene = genes_here[j],
                           est = unname(cc["est"]), var = unname(cc["var"]),
                           n_events = ev, row.names = NULL)
  }
  do.call(rbind, out)
}

## ============================================================================
## 16.3  RUN
## ============================================================================
log_msg("Section 16 start. %d genes. ADJUST_DEPTH = %s. run_id = %s.",
        length(RND_GENES), RND_ADJUST_DEPTH, RUN_ID)

rnd_acc <- list()
for (i in seq_len(nrow(work))) {
  path <- work$path[i]; sid <- work$sample_id[i]
  if (is.na(path) || !file.exists(path)) next
  if (!sid %in% RND_SECTIONS) next
  log_msg("[rnd %d/%d] %s", i, nrow(work), sid)
  r <- tryCatch(rnd_one_tumour(path, sid),
                error = function(e) { log_msg("  rnd skipped: %s", conditionMessage(e)); NULL })
  if (!is.null(r) && nrow(r)) rnd_acc[[length(rnd_acc) + 1]] <- r
}
rnd_all <- bind(rnd_acc)
if (nrow(rnd_all) == 0) stop("Section 16: no random gene produced a usable fit.")

## ============================================================================
## 16.4  POOL EACH GENE ACROSS SECTIONS
## ============================================================================
rnd_pooled <- rnd_all %>%
  group_by(gene) %>% group_modify(~ pool_with_p(.x$est, .x$var, .x$sample_id)) %>% ungroup()

n_before <- nrow(rnd_pooled)
rnd_pooled <- rnd_pooled[is.finite(rnd_pooled$est) & rnd_pooled$k >= RND_MIN_K, ]
if (!nrow(rnd_pooled)) stop(sprintf("Section 16: no random gene was pooled over >= %d %s.", RND_MIN_K, K_UNIT))
write_out(rnd_pooled, "RANDOM_pooled_per_gene")

cat("\n", strrep("=", 74), "\n", sep = "")
cat(sprintf("SECTION 16  [%s]  run_id %s\n",
            if (RND_ADJUST_DEPTH) "DEPTH ADJUSTED" else "UNADJUSTED", RUN_ID))
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  genes attempted / pooled / kept (k >= %d) : %d / %d / %d\n",
            RND_MIN_K, length(RND_GENES), n_before, nrow(rnd_pooled)))

cat("\n  distribution of the pooled slope difference:\n")
qq <- quantile(rnd_pooled$est, c(0, .025, .25, .5, .75, .975, 1))
for (nmq in names(qq)) cat(sprintf("     %-6s %+.4f\n", nmq, qq[[nmq]]))
cat(sprintf("     significant at 0.05 in %d of %d (%.1f%%)\n",
            sum(rnd_pooled$p < 0.05, na.rm = TRUE), nrow(rnd_pooled),
            100 * mean(rnd_pooled$p < 0.05, na.rm = TRUE)))
cat("     About 5% would be expected. Well above that means the model p values\n")
cat("     are inflated, most likely by spatial autocorrelation.\n")

## ============================================================================
## 16.5  WHERE DO THE THREE ISOFORMS SIT?
##       v5 (C7): each isoform is re-pooled on the sections the reference was
##       built on, so observed and reference share the same sections.
## ============================================================================
obs <- NULL
if (exists("iso_con")) {
  ic <- iso_con[iso_con$contrast == "interaction_diff" & iso_con$sample_id %in% RND_SECTIONS, ]
  if (nrow(ic))
    obs <- as.data.frame(ic %>% group_by(outcome) %>%
                           group_modify(~ pool_with_p(.x$est, .x$var, .x$sample_id)) %>% ungroup())
} else {
  iso_file <- file.path(agg, outf("ISOFORM_slope_difference_pooled"))
  if (file.exists(iso_file)) {
    io <- read.csv(iso_file, stringsAsFactors = FALSE)
    obs <- io[io$contrast == "interaction_diff", ]
    cat("\n  NOTE: section 15 is not in this session; isoforms are pooled on all of\n")
    cat("  their own sections, not only the reference sections.\n")
  }
}

ref_summary <- NULL
if (!is.null(obs) && nrow(obs)) {
  cat("\n--- where the three isoforms sit (same sections as the reference) ---\n")
  ref_summary <- do.call(rbind, lapply(seq_len(nrow(obs)), function(i) {
    o <- obs[i, ]
    p_dir <- mean(rnd_pooled$est <= o$est)
    p_mag <- mean(abs(rnd_pooled$est) >= abs(o$est))
    data.frame(outcome = o$outcome, est = o$est, model_p = o$p, k = o$k,
               percentile = 100 * p_dir, emp_p_directional = p_dir,
               emp_p_magnitude = p_mag, n_reference = nrow(rnd_pooled), row.names = NULL)
  }))
  write_out(ref_summary, "RANDOM_reference_summary")
  
  for (i in seq_len(nrow(ref_summary))) with(ref_summary[i, ],
                                             cat(sprintf("  %-6s effect %+.4f  k=%2d  percentile %5.1f  empirical p: directional %.4f, magnitude %.4f%s\n",
                                                         outcome, est, as.integer(k), percentile, emp_p_directional, emp_p_magnitude,
                                                         if (k < RND_MIN_K) "   (k below the reference minimum)" else "")))
  cat("\n  Directional = fraction of random genes at least as negative. Quote it\n")
  cat("  and show the other; do not pick whichever is smaller. P4HA2 and P4HA3\n")
  cat("  are compared with a reference matched on P4HA1's detection rate, so\n")
  cat("  their percentiles are only indicative.\n")
  
  a1 <- ref_summary[ref_summary$outcome == "P4HA1", ]
  if (nrow(a1) == 1) {
    cat("\n--- HOW TO SAY IT ---\n")
    if (a1$emp_p_directional <= 0.01) {
      cat(sprintf("  P4HA1 is in the extreme %.1f%% of %d detection-matched genes. The\n",
                  100 * a1$emp_p_directional, a1$n_reference))
      cat("  behaviour is not generic, and since the reference carries the same\n")
      cat("  spatial structure, autocorrelation cannot explain the position.\n")
    } else if (a1$emp_p_directional <= 0.05) {
      cat(sprintf("  P4HA1 is in the top %.1f%%. Real but not dramatic. Report the\n",
                  100 * a1$emp_p_directional))
      cat("  percentile alongside the model p and do not lean on the model p alone.\n")
    } else {
      cat(sprintf("  P4HA1 sits at the %.0fth percentile, empirical p = %.3f. It is NOT\n",
                  a1$percentile, a1$emp_p_directional))
      cat("  unusual against this reference. The behaviour is largely generic and\n")
      cat("  the model p is measuring the wrong null. Report this plainly.\n")
    }
  }
}

cat("\nWHAT THIS SETTLES AND WHAT IT DOES NOT\n")
cat("  Does: whether the behaviour is specific to P4HA1 or common in this tissue.\n")
cat("  Does not: fix the CIs, which stay too narrow for every gene here.\n")
cat("  Matching is on detection rate only; the reference may contain fibroblast-\n")
cat("  or hypoxia-responsive genes that behave the same way for real reasons.\n")
cat(strrep("=", 74), "\n", sep = "")

## ============================================================================
## 16.6  FIGURE 8
## ============================================================================
if (exists("FIG_STYLES") && exists("theme_fig")) {
  for (sn in names(FIG_STYLES)) {
    S <- FIG_STYLES[[sn]]
    p8 <- ggplot(rnd_pooled, aes(est)) +
      geom_histogram(bins = 50, fill = "grey72", colour = "white", linewidth = 0.2) +
      geom_vline(xintercept = 0, linetype = "22", linewidth = S$lwd * 0.7, colour = "grey40") +
      labs(x = "Pooled slope difference (collagen minus ER stress, log odds per SD)",
           y = "Number of random genes",
           caption = sprintf("%d genes matched to P4HA1 on detection rate, same model, sections and pooling.%s",
                             nrow(rnd_pooled),
                             if (RND_ADJUST_DEPTH) " Sequencing depth adjusted." else "")) +
      theme_fig(S)
    
    if (!is.null(obs) && nrow(obs)) {
      cols <- c(P4HA1 = CLR$neg, P4HA2 = CLR$pos, P4HA3 = "grey25")
      od <- obs[obs$outcome %in% names(cols), ]
      ymax <- max(ggplot_build(p8)$data[[1]]$count, na.rm = TRUE)
      p8 <- p8 +
        geom_vline(data = od, aes(xintercept = est, colour = outcome), linewidth = S$lwd * 1.5) +
        geom_text(data = od, aes(x = est, y = ymax, label = outcome, colour = outcome),
                  angle = 90, hjust = 1.1, vjust = -0.4, size = S$lab, show.legend = FALSE) +
        scale_colour_manual(values = cols, guide = "none")
    }
    
    if (isTRUE(S$titles))
      p8 <- p8 + labs(title = "Is the behaviour specific to P4HA1?",
                      subtitle = "P4HA1 against a reference of detection-matched random genes")
    
    save_fig(p8, S, sub("\\.csv$", "", outf("Figure8_random_gene_control")),
             if (isTRUE(S$split)) S$w2 else S$w1 * 1.5,
             if (isTRUE(S$split)) S$h_full else MM(66))
  }
}

log_msg("Section 16 complete.")

################################################################################
## 17. DEPTH ADJUSTED vs UNADJUSTED, AND RUN-TO-RUN REPRODUCIBILITY
##     Run sections 13, 15, 16 once with ADJUST_DEPTH FALSE and once with
##     TRUE, each in a FRESH session from section 0. Then run this. It
##     verifies the two files are different runs before comparing them, and
##     diffs the module fingerprints so a change can be traced.
################################################################################

cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 17  COMPARISON AND REPRODUCIBILITY CHECK\n")
cat(strrep("=", 74), "\n", sep = "")

cmp_read <- function(base) {
  f0 <- file.path(agg, paste0(base, ".csv"))
  f1 <- file.path(agg, paste0(base, "_depthadj.csv"))
  if (!file.exists(f0) || !file.exists(f1)) return(NULL)
  a <- read.csv(f0, stringsAsFactors = FALSE)
  b <- read.csv(f1, stringsAsFactors = FALSE)
  
  if ("adjust_depth" %in% names(a) && "adjust_depth" %in% names(b)) {
    if (isTRUE(all(a$adjust_depth == b$adjust_depth))) {
      cat(sprintf("\n  *** %s: BOTH files report adjust_depth = %s. These are not two\n",
                  base, unique(a$adjust_depth)[1]))
      cat("  different runs. Re-run from a fresh session with ADJUST_DEPTH <- TRUE.\n")
      cat("  Comparison skipped. ***\n")
      return(NULL)
    }
    if ("run_id" %in% names(a) && identical(unique(a$run_id), unique(b$run_id))) {
      cat(sprintf("\n  *** %s: both files carry the same run_id. Skipped. ***\n", base))
      return(NULL)
    }
  } else {
    cat(sprintf("\n  NOTE: %s has no provenance stamp; check the run manifests by hand.\n", base))
  }
  list(unadj = a, adj = b)
}

cmp_line <- function(label, e0, l0, u0, p0, e1, l1, u1, p1) {
  shift <- e1 - e0
  pct   <- if (is.finite(e0) && abs(e0) > 1e-12) 100 * shift / abs(e0) else NA_real_
  sig0  <- is.finite(p0) && p0 < 0.05
  sig1  <- is.finite(p1) && p1 < 0.05
  cat(sprintf("  %-22s\n", label))
  cat(sprintf("     unadjusted %+8.4f (%+.4f to %+.4f)  p = %s\n", e0, l0, u0, fmt_p(p0)))
  cat(sprintf("     adjusted   %+8.4f (%+.4f to %+.4f)  p = %s\n", e1, l1, u1, fmt_p(p1)))
  cat(sprintf("     shift      %+8.4f  (%s)%s\n", shift,
              if (is.finite(pct)) sprintf("%+.1f%%", pct) else "n/a",
              if (sig0 != sig1) "   <-- SIGNIFICANCE CHANGED" else ""))
}

## ---- 17.1  section 13 headline --------------------------------------------
j <- cmp_read("JOINT_slope_difference_pooled")
if (is.null(j)) {
  cat("\n  17.1 unavailable.\n")
} else {
  cat("\n--- 17.1  section 13, the slope difference for P4HA1 ---\n")
  for (ct in intersect(j$unadj$contrast, j$adj$contrast)) {
    a <- j$unadj[j$unadj$contrast == ct, ]; b <- j$adj[j$adj$contrast == ct, ]
    if (nrow(a) && nrow(b)) cmp_line(ct, a$est, a$lower, a$upper, a$p, b$est, b$lower, b$upper, b$p)
  }
}

## ---- 17.2  section 15 per outcome -----------------------------------------
s <- cmp_read("ISOFORM_slope_difference_pooled")
if (is.null(s)) {
  cat("\n  17.2 unavailable.\n")
} else {
  cat("\n--- 17.2  section 15, per outcome ---\n")
  a0 <- s$unadj[s$unadj$contrast == "interaction_diff", ]
  a1 <- s$adj[s$adj$contrast == "interaction_diff", ]
  for (og in intersect(a0$outcome, a1$outcome)) {
    x <- a0[a0$outcome == og, ]; y <- a1[a1$outcome == og, ]
    cmp_line(og, x$est, x$lower, x$upper, x$p, y$est, y$lower, y$upper, y$p)
    if (x$k != y$k)
      cat(sprintf("     NOTE: k changed from %d to %d. Say so if it happens.\n", as.integer(x$k), as.integer(y$k)))
  }
}

## ---- 17.3  section 16 percentiles -----------------------------------------
r <- cmp_read("RANDOM_reference_summary")
if (is.null(r)) {
  cat("\n  17.3 unavailable.\n")
} else {
  cat("\n--- 17.3  section 16, empirical percentile ---\n")
  cat(sprintf("  %-8s %14s %14s %10s\n", "outcome", "unadj pctile", "adj pctile", "shift"))
  for (og in intersect(r$unadj$outcome, r$adj$outcome)) {
    x <- r$unadj[r$unadj$outcome == og, ]; y <- r$adj[r$adj$outcome == og, ]
    cat(sprintf("  %-8s %13.1f %14.1f %+10.1f\n",
                og, x$percentile, y$percentile, y$percentile - x$percentile))
  }
  cat("\n  The percentile is the number to quote publicly.\n")
}

## ---- 17.4  module fingerprint diff ----------------------------------------
{
  f0 <- file.path(agg, "JOINT_module_fingerprint.csv")
  f1 <- file.path(agg, "JOINT_module_fingerprint_depthadj.csv")
  if (file.exists(f0) && file.exists(f1)) {
    p0 <- read.csv(f0, stringsAsFactors = FALSE)
    p1 <- read.csv(f1, stringsAsFactors = FALSE)
    key <- c("sample_id", "module")
    m <- merge(p0[, c(key, "n_genes_used", "score_mean", "score_sd")],
               p1[, c(key, "n_genes_used", "score_mean", "score_sd")],
               by = key, suffixes = c("_a", "_b"))
    m$dmean  <- m$score_mean_b - m$score_mean_a
    m$dsd    <- m$score_sd_b   - m$score_sd_a
    m$dgenes <- m$n_genes_used_b - m$n_genes_used_a
    changed <- m[abs(m$dmean) > 1e-9 | abs(m$dsd) > 1e-9 | m$dgenes != 0, ]
    cat("\n--- 17.4  did the module scores themselves change between runs? ---\n")
    if (!nrow(changed)) {
      cat("  No. All module scores are identical across the two runs, so any\n")
      cat("  difference in results comes from the model, not from scoring.\n")
    } else {
      cat(sprintf("  YES: %d of %d module-by-section scores differ. Something upstream\n", nrow(changed), nrow(m)))
      cat("  changed (gene lists, Seurat version, or a contaminated session).\n")
      changed <- changed[order(-abs(changed$dmean)), ]
      for (i in seq_len(min(10, nrow(changed)))) with(changed[i, ],
                                                      cat(sprintf("     %-30s %-18s d_mean %+.6f  d_sd %+.6f  d_genes %+d\n",
                                                                  sample_id, module, dmean, dsd, as.integer(dgenes))))
    }
  } else {
    cat("\n  17.4 unavailable: needs both fingerprint files from section 13.\n")
  }
}

## ---- 17.5  inflation both ways --------------------------------------------
{
  f0 <- file.path(agg, "RANDOM_pooled_per_gene.csv")
  f1 <- file.path(agg, "RANDOM_pooled_per_gene_depthadj.csv")
  if (file.exists(f0) && file.exists(f1)) {
    g0 <- read.csv(f0, stringsAsFactors = FALSE)
    g1 <- read.csv(f1, stringsAsFactors = FALSE)
    cat("\n--- 17.5  how many random genes are themselves significant? ---\n")
    cat(sprintf("  unadjusted : %d of %d (%.1f%%)\n",
                sum(g0$p < 0.05, na.rm = TRUE), nrow(g0), 100*mean(g0$p < 0.05, na.rm = TRUE)))
    cat(sprintf("  adjusted   : %d of %d (%.1f%%)\n",
                sum(g1$p < 0.05, na.rm = TRUE), nrow(g1), 100*mean(g1$p < 0.05, na.rm = TRUE)))
    cat("  About 5% expected. Depth adjustment is not expected to fix this; the\n")
    cat("  cause is spatial autocorrelation.\n")
  }
}
cat(strrep("=", 74), "\n", sep = "")

log_msg("Section 17 complete. run_id %s, ADJUST_DEPTH %s.", RUN_ID, ADJUST_DEPTH)

################################################################################
## 18. THE METABOLIC AXIS: 2-OXOGLUTARATE SUPPLY INSTEAD OF HYPOXIA  (Figure 9)
##
##   P4H consumes oxygen AND 2-oxoglutarate in the same reaction. Sections 13
##   and 15 came in through the oxygen door; this section comes in through the
##   cofactor door, with the same model and zMet in place of zHyp.
##
##   v5 (C1) -- WHAT CHANGED AND WHY
##   * The colorectal-cancer survival signature (ASRGL1, MECP2, NOS2, NOS3)
##     is GONE. It was selected in another cancer to predict survival, not to
##     measure 2-oxoglutarate supply, so it had no meaning as a BRCA measure.
##   * Every module here is a named standard gene set (or a subset of one,
##     clearly labelled), and its detection in THIS cohort is printed at run
##     time. No coverage figure from another cohort is quoted.
##   * Each metabolic module is scored on its own. Before, a section that
##     could not score one of them was removed from ALL of them.
##   * The hypoxia-axis result (section 13) is re-pooled on exactly the same
##     BRCA sections as the metabolic result before they are compared, and is
##     drawn as the first panel of Figure 9.
##
##   GENE SETS
##   Glutamine  GOBP_GLUTAMINE_FAMILY_AMINO_ACID_METABOLIC_PROCESS (M11431,
##              GO:0009064, 71 genes), verbatim. Primary.
##   OGcore     12 enzymes that make or consume 2-oxoglutarate, a HAND-PICKED
##              subset of the set above. Secondary.
##   Weaknesses, to state before anyone asks: the GO term also covers the
##   urea cycle, proline metabolism and NOS genes; and transcript abundance
##   is a weak proxy for metabolite level. Section 19 tests whether the
##   module's genes even move together in this tissue.
##
##   THE GATE: if the metabolic module largely mirrors hypoxia, agreement
##   between the two axes is one finding reported twice. 18.5 measures this
##   before any result is shown.
################################################################################

stopifnot(exists("MODULES"), exists("meta_pool"), exists("work"), exists("agg"))
if (!exists("pool_with_p") || !exists("fmt_p") || !exists("assert_run_state"))
  stop("Section 18 needs the 12b helper block (pool_with_p, fmt_p, assert_run_state).")

## ============================================================================
## 18.0  GENE SETS AND SETTINGS
## ============================================================================
## MSigDB M11431, GO:0009064, release 2025.1.Hs, verbatim. Do not edit by hand.
GLUTAMINE_GENES <- c(
  "GLYATL1B","MTHFS","FTCD","AZIN2","CLN3","NOXRED1","UROC1","CPS1",
  "AMDHD1","DAO","NAGS","FABP1","FAH","FOLH1B","SIRT4","FOLH1",
  "FPGS","DDAH1","SLC7A11","GAD1","GAD2","PHGDH","BLOC1S6","GGT1",
  "GLS2","GCLC","GCLM","GLS","GLUD1","GLUD2","GLUL","GOT1",
  "GOT2","RIMKLA","PYCR2","HAL","ARG1","MIR21","MECP2","ASL",
  "ASS1","MTTP","NOS1","NOS2","NOS3","ATP2B4","OAT","ODC1",
  "OTC","AADAT","LGSN","AZIN1","PFAS","PRODH","NIT2","RIMKLB",
  "PYCR1","ALDH18A1","PRODH2","PYCR3","TAT","UCP2","CAD","ALDH5A1",
  "AGMAT","DGLUCY","ATCAY","ALDH4A1","ARHGAP11B","GLYATL1","NR1H4")

## Enzymes that produce or consume 2-oxoglutarate, or sit one step away:
##   GLS, GLS2 (glutamine -> glutamate); GLUD1, GLUD2 (glutamate -> 2-OG);
##   GOT1, GOT2 (transamination); GLUL (reverse direction); OAT; ALDH5A1;
##   ALDH4A1; CAD (competing glutamine sink); PHGDH (transaminates to 2-OG).
## HAND-SELECTED, hence secondary only.
OG_CORE_GENES <- c("GLS","GLS2","GLUD1","GLUD2","GOT1","GOT2","GLUL","OAT",
                   "ALDH5A1","ALDH4A1","CAD","PHGDH")

MODULES$Glutamine <- GLUTAMINE_GENES
MODULES$OGcore    <- OG_CORE_GENES
MODULES$PaperSig  <- NULL          # v5 (C1): removed, see header

METAB_VARIANTS     <- c("Glutamine", "OGcore")
METAB_PRIMARY      <- "Glutamine"
METAB_ES           <- "ERstressStrict"
METAB_BASE_MODULES <- c("Hypoxia", "CollagenSynthesis", "iCAF", METAB_ES, "ERAD")
METAB_MIN_SPOTS    <- 30
METAB_Z_GRID       <- seq(-2, 2, length.out = 41)
METAB_ALPHA        <- 0.05
METAB_ADJUST_DEPTH <- ADJUST_DEPTH     # derived, never set by hand

METAB_COLLIN_WARN <- 0.50
METAB_COLLIN_STOP <- 0.80

MT_COL <- "zCol"; MT_ER <- "zER"; MT_MET <- "zMet"
MI_COL <- paste0(MT_MET, ":", MT_COL)
MI_ER  <- paste0(MT_MET, ":", MT_ER)

LAB_COL <- "Collagen synthesis"
LAB_ER  <- "ER stress"
LAB_DIF <- "difference"

NICE_AXIS <- c(Hypoxia   = "Hypoxia axis (section 13)",
               Glutamine = "Glutamine family (GO, 71 genes)",
               OGcore    = "2-OG core (12 genes, hand-picked)")

assert_run_state("Section 18")

## ============================================================================
## 18.1  GENE-SET AUDIT, halts on failure
## ============================================================================
cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 18.1  GENE-SET AUDIT\n")
cat(strrep("=", 74), "\n", sep = "")

DANGER <- c("P4HA1","P4HA2","P4HA3","P4HB","PLOD1","PLOD2","PLOD3","PPIB")

metab_fail <- FALSE
for (nm in METAB_VARIANTS) {
  g   <- MODULES[[nm]]
  bad <- intersect(DANGER, g)
  cat(sprintf("\n  %s  (n = %d genes, need %d detected)\n", nm, length(g), min_genes_for(g)))
  cat(sprintf("    circularity : %s\n",
              if (length(bad)) paste0("*** FAIL: ", paste(bad, collapse = ", "), " ***") else "clean"))
  if (length(bad)) metab_fail <- TRUE
  ov_any <- FALSE
  for (other in METAB_BASE_MODULES) {
    ov <- intersect(g, MODULES[[other]])
    if (length(ov)) {
      cat(sprintf("    *** overlaps %s: %s ***\n", other, paste(ov, collapse = ", ")))
      ov_any <- TRUE; metab_fail <- TRUE
    }
  }
  if (!ov_any) cat("    overlap with the predictor modules : none\n")
}
if (metab_fail)
  stop("Section 18 halted: a metabolic module is circular or overlaps another predictor.")
cat("\n  Audit passed.\n")

## ============================================================================
## 18.2  SCORE, STANDARDISE, AND MEASURE THE COLLINEARITY GATE
##       The base modules must all score (they are in every model). Each
##       metabolic module is scored on its own and skipped on its own.
## ============================================================================
metab_expr_one <- function(rds_path, sid) {
  obj <- readRDS(rds_path)
  obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
  DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else Assays(obj)[1]
  if (!"P4HA1" %in% rownames(obj)) return(NULL)
  
  score_one <- function(obj, nm) {
    set.seed(SEED)
    g <- intersect(MODULES[[nm]], rownames(obj))
    o <- tryCatch(add_module_score_fallback(obj, features = list(g), name = nm,
                                            assay = DefaultAssay(obj), nbin = 20, ctrl = 50),
                  error = function(e) NULL)
    if (is.null(o) || !paste0(nm, "1") %in% colnames(o@meta.data)) return(NULL)
    names(o@meta.data)[names(o@meta.data) == paste0(nm, "1")] <- nm
    o
  }
  
  cov <- list()
  for (nm in c(METAB_BASE_MODULES, METAB_VARIANTS)) {
    full <- MODULES[[nm]]; g <- intersect(full, rownames(obj))
    cov[[nm]] <- data.frame(sample_id = sid, module = nm, n_total = length(full),
                            n_found = length(g), n_required = min_genes_for(full), row.names = NULL)
  }
  
  for (nm in METAB_BASE_MODULES) {
    if (cov[[nm]]$n_found < cov[[nm]]$n_required) {
      log_msg("  %s: %s has %d/%d genes (need %d) -> section skipped",
              sid, nm, cov[[nm]]$n_found, cov[[nm]]$n_total, cov[[nm]]$n_required)
      return(NULL)
    }
    o <- score_one(obj, nm); if (is.null(o)) return(NULL); obj <- o
  }
  have <- character(0)
  for (nm in METAB_VARIANTS) {
    if (cov[[nm]]$n_found < cov[[nm]]$n_required) {
      log_msg("  %s: %s has %d/%d genes (need %d) -> this module skipped here",
              sid, nm, cov[[nm]]$n_found, cov[[nm]]$n_total, cov[[nm]]$n_required)
      next
    }
    o <- score_one(obj, nm); if (is.null(o)) next
    obj <- o; have <- c(have, nm)
  }
  if (!length(have)) return(NULL)
  
  scored <- c(METAB_BASE_MODULES, have)
  vars <- c("P4HA1", scored)
  cnt_col <- grep("^nCount_", colnames(obj@meta.data), value = TRUE)
  if (METAB_ADJUST_DEPTH) {
    if (!length(cnt_col)) return(NULL)
    vars <- c(vars, cnt_col[1])
  }
  d <- FetchData(obj, vars = vars, layer = "data")
  d$y <- as.integer(d$P4HA1 > 0)
  
  ## collinearity on the RAW scores, one row per metabolic module
  cc <- suppressWarnings(cor(d[, scored], method = "spearman"))
  collin <- data.frame(
    sample_id      = sid,
    variant        = have,
    rho_vs_Hypoxia = unname(cc[have, "Hypoxia"]),
    rho_vs_Col     = unname(cc[have, "CollagenSynthesis"]),
    rho_vs_ER      = unname(cc[have, METAB_ES]),
    rho_vs_primary = if (METAB_PRIMARY %in% have) unname(cc[have, METAB_PRIMARY]) else NA_real_,
    row.names = NULL)
  
  zz <- function(x) { s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) return(NULL)
  as.numeric((x - mean(x, na.rm = TRUE)) / s) }
  for (nm in METAB_BASE_MODULES) {
    v <- zz(d[[nm]]); if (is.null(v)) return(NULL)
    d[[paste0("z", nm)]] <- v
  }
  for (nm in have) {
    v <- zz(d[[nm]])
    if (is.null(v)) { have <- setdiff(have, nm); next }
    d[[paste0("z", nm)]] <- v
  }
  if (!length(have)) return(NULL)
  if (METAB_ADJUST_DEPTH) {
    v <- zz(log1p(d[[cnt_col[1]]])); if (is.null(v)) return(NULL); d$zDepth <- v
  }
  
  attr(d, "have")     <- have
  attr(d, "collin")   <- collin[collin$variant %in% have, , drop = FALSE]
  attr(d, "coverage") <- do.call(rbind, cov)
  d
}

## ============================================================================
## 18.3  FIT ONE SECTION, ONE METABOLIC MODULE
##       Same shape as sections 13 and 15, with zMet replacing zHyp.
## ============================================================================
fit_metab_one <- function(d, sid, variant) {
  out <- list(coef = NULL, contrast = NULL, curve = NULL)
  if (is.null(d) || !variant %in% attr(d, "have")) return(out)
  if (nrow(d) < METAB_MIN_SPOTS || length(unique(d$y)) != 2) return(out)
  
  dd <- d
  dd$zMet  <- d[[paste0("z", variant)]]
  dd$zCol  <- d$zCollagenSynthesis
  dd$zER   <- d[[paste0("z", METAB_ES)]]
  
  rhs <- paste0(MT_MET, " * ", MT_COL, " + ", MT_MET, " * ", MT_ER, " + ziCAF + zERAD",
                if (METAB_ADJUST_DEPTH) " + zDepth" else "")
  f <- tryCatch(glm(as.formula(paste("y ~", rhs)), data = dd, family = binomial()),
                error = function(e) NULL)
  if (is.null(f) || !f$converged) return(out)
  
  b <- coef(f); V <- vcov(f)
  has_depth <- "zDepth" %in% names(b)
  if (!identical(has_depth, isTRUE(METAB_ADJUST_DEPTH)))
    stop(sprintf("Section 18, %s / %s: flag and fitted model disagree about zDepth.", sid, variant))
  
  need <- c(MT_COL, MT_ER, MI_COL, MI_ER)
  if (!all(need %in% names(b)) || any(!is.finite(b[need]))) return(out)
  
  s  <- summary(f)$coefficients
  tm <- setdiff(rownames(s), "(Intercept)")
  out$coef <- data.frame(sample_id = sid, variant = variant, term = tm,
                         log_OR = s[tm, "Estimate"], SE = s[tm, "Std. Error"],
                         p = s[tm, "Pr(>|z|)"], row.names = NULL)
  
  Lvec <- function(w) { L <- setNames(numeric(length(b)), names(b))
  for (nm in names(w)) L[nm] <- w[[nm]]; L }
  est_of <- function(L) as.numeric(L %*% b)
  var_of <- function(L) as.numeric(t(L) %*% V %*% L)
  
  L_int  <- Lvec(setNames(list(1, -1), c(MI_COL, MI_ER)))
  L_main <- Lvec(setNames(list(1, -1), c(MT_COL, MT_ER)))
  out$contrast <- data.frame(
    sample_id = sid, variant = variant,
    contrast = c("interaction_diff", "main_diff"),
    est = c(est_of(L_int), est_of(L_main)),
    var = c(var_of(L_int), var_of(L_main)),
    n_spots = nrow(dd), n_events = sum(dd$y), row.names = NULL)
  
  out$curve <- do.call(rbind, lapply(METAB_Z_GRID, function(z) {
    Lc <- Lvec(setNames(list(1, z),         c(MT_COL, MI_COL)))
    Le <- Lvec(setNames(list(1, z),         c(MT_ER,  MI_ER)))
    Ld <- Lvec(setNames(list(1, -1, z, -z), c(MT_COL, MT_ER, MI_COL, MI_ER)))
    data.frame(sample_id = sid, variant = variant, z = z,
               series = c(LAB_COL, LAB_ER, LAB_DIF),
               est = c(est_of(Lc), est_of(Le), est_of(Ld)),
               var = c(var_of(Lc), var_of(Le), var_of(Ld)), row.names = NULL)
  }))
  out
}

## ============================================================================
## 18.4  RUN THE COHORT
## ============================================================================
log_msg("Section 18 start. Modules: %s. ADJUST_DEPTH = %s.",
        paste(METAB_VARIANTS, collapse = ", "), METAB_ADJUST_DEPTH)

metab_acc <- list(coef = list(), contrast = list(), curve = list(),
                  collin = list(), coverage = list())
for (i in seq_len(nrow(work))) {
  path <- work$path[i]; sid <- work$sample_id[i]
  if (is.na(path) || !file.exists(path)) next
  log_msg("[metab %d/%d] %s", i, nrow(work), sid)
  d <- tryCatch(metab_expr_one(path, sid),
                error = function(e) { log_msg("  metab skipped: %s", conditionMessage(e)); NULL })
  if (is.null(d)) next
  metab_acc$collin[[length(metab_acc$collin) + 1]]     <- attr(d, "collin")
  metab_acc$coverage[[length(metab_acc$coverage) + 1]] <- attr(d, "coverage")
  for (vr in METAB_VARIANTS) {
    r <- fit_metab_one(d, sid, vr)
    for (k in c("coef", "contrast", "curve"))
      if (!is.null(r[[k]]) && nrow(r[[k]])) metab_acc[[k]][[length(metab_acc[[k]]) + 1]] <- r[[k]]
  }
}

mcoef <- bind(metab_acc$coef)
mcon  <- bind(metab_acc$contrast)
mcur  <- bind(metab_acc$curve)
mcol  <- bind(metab_acc$collin)
mcov  <- bind(metab_acc$coverage)

if (nrow(mcon) == 0) stop("Section 18: no section produced a usable fit.")

write_out(mcoef, "METAB_coefficients_per_sample")
write_out(mcon,  "METAB_contrasts_per_sample")
write_out(mcol,  "METAB_collinearity")
write_out(mcov,  "METAB_gene_coverage")

## ============================================================================
## 18.5  THE COLLINEARITY GATE, reported before any result
## ============================================================================
cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 18.5  IS THIS REALLY A SECOND, INDEPENDENT AXIS?\n")
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  Spearman correlation with the hypoxia module (median over %ss):\n\n", UNIT_LABEL))

gate_verdict <- list()
for (vr in METAB_VARIANTS) {
  cv <- mcol[mcol$variant == vr, ]
  if (!nrow(cv)) { gate_verdict[[vr]] <- "not scored"; next }
  med <- unit_median(cv$rho_vs_Hypoxia, cv$sample_id)
  status <- if (abs(med) >= METAB_COLLIN_STOP) "TOO HIGH"
  else if (abs(med) >= METAB_COLLIN_WARN) "HIGH"
  else "acceptable"
  gate_verdict[[vr]] <- status
  cat(sprintf("    %-10s median %+.2f   range %+.2f to %+.2f   -> %s\n",
              vr, med, min(cv$rho_vs_Hypoxia, na.rm = TRUE), max(cv$rho_vs_Hypoxia, na.rm = TRUE), status))
}
cat(sprintf("\n  thresholds: |rho| >= %.2f is HIGH, >= %.2f is TOO HIGH\n",
            METAB_COLLIN_WARN, METAB_COLLIN_STOP))
cat("    acceptable : the axes are largely separate; agreement is independent support.\n")
cat("    HIGH       : substantial shared variance; agreement is weaker than it looks.\n")
cat("    TOO HIGH   : the metabolic module is largely hypoxia in another form.\n")

if (identical(gate_verdict[[METAB_PRIMARY]], "TOO HIGH")) {
  cat("\n  *** The PRIMARY module fails the gate. Everything below is still\n")
  cat("  computed, but it cannot be described as an independent axis. ***\n")
}

pc <- mcol[mcol$variant == METAB_PRIMARY, ]
if (nrow(pc)) {
  cat("\n  other correlations of the primary module (median):\n")
  oc <- mcol[mcol$variant == "OGcore", ]
  if (nrow(oc)) cat(sprintf("    vs the 12-gene core subset : %+.2f\n", unit_median(oc$rho_vs_primary, oc$sample_id)))
  cat(sprintf("    vs collagen synthesis      : %+.2f\n", unit_median(pc$rho_vs_Col, pc$sample_id)))
  cat(sprintf("    vs ER stress               : %+.2f\n", unit_median(pc$rho_vs_ER, pc$sample_id)))
}

cat("\n  gene coverage per module in THIS cohort (across sections):\n")
cs <- mcov %>% group_by(module) %>%
  summarise(median_found = median(n_found), total = max(n_total),
            min_found = min(n_found), required = max(n_required),
            n_below = sum(n_found < n_required), n = dplyr::n(), .groups = "drop")
## v5 (C5): %d on a median crashed with an even number of sections
for (i in seq_len(nrow(cs))) with(cs[i, ],
                                  cat(sprintf("    %-20s %5.1f / %3.0f detected (min %3.0f, need %3.0f)  below threshold in %d of %d sections\n",
                                              module, median_found, total, min_found, required, as.integer(n_below), as.integer(n))))

## ============================================================================
## 18.6  POOL AND REPORT
## ============================================================================
mcoef_pooled <- mcoef %>%
  group_by(variant, term) %>% group_modify(~ pool_with_p(.x$log_OR, .x$SE^2, .x$sample_id)) %>% ungroup() %>%
  mutate(OR = exp(est), lower_OR = exp(lower), upper_OR = exp(upper))
write_out(mcoef_pooled, "METAB_coefficients_pooled")

mcon_pooled <- mcon %>%
  group_by(variant, contrast) %>% group_modify(~ pool_with_p(.x$est, .x$var, .x$sample_id)) %>% ungroup()
write_out(mcon_pooled, "METAB_slope_difference_pooled")

mcur_pooled <- mcur %>%
  group_by(variant, series, z) %>% group_modify(~ pool_with_p(.x$est, .x$var, .x$sample_id)) %>% ungroup() %>%
  mutate(OR = exp(est), lower_OR = exp(lower), upper_OR = exp(upper)) %>%
  arrange(variant, series, z)
write_out(mcur_pooled, "METAB_slope_curves_pooled")

## ---- the hypoxia axis on EXACTLY the same sections -------------------------
## v5 (C1): the two axes are compared on one set of BRCA sections. Sections
## fitted on only one axis are named, not silently mixed in.
prim_secs <- unique(mcon$sample_id[mcon$variant == METAB_PRIMARY & mcon$contrast == "interaction_diff"])
hyp_ref <- NULL; hyp_curve <- NULL; hyp_note <- NULL
if (exists("jcon_all") && exists("jcur_all")) {
  only_met <- setdiff(prim_secs, jcon_all$sample_id)
  only_hyp <- setdiff(unique(jcon_all$sample_id), prim_secs)
  common   <- intersect(prim_secs, jcon_all$sample_id)
  if (length(common)) {
    jc <- jcon_all[jcon_all$contrast == "interaction_diff" & jcon_all$sample_id %in% common, ]
    hyp_ref <- pool_with_p(jc$est, jc$var, jc$sample_id)
    hyp_curve <- jcur_all[jcur_all$sample_id %in% common, ] %>%
      group_by(series, z) %>% group_modify(~ pool_with_p(.x$est, .x$var, .x$sample_id)) %>% ungroup() %>%
      mutate(OR = exp(est), lower_OR = exp(lower), upper_OR = exp(upper), variant = "Hypoxia")
    ## the metabolic result on the same common sections, for the side by side
    mc <- mcon[mcon$variant == METAB_PRIMARY & mcon$contrast == "interaction_diff" & mcon$sample_id %in% common, ]
    met_common <- pool_with_p(mc$est, mc$var, mc$sample_id)
  }
  hyp_note <- sprintf("common sections %d; only on the metabolic axis %d; only on the hypoxia axis %d",
                      length(common), length(only_met), length(only_hyp))
} else {
  hyp_note <- "section 13 is not in this session; no like-for-like comparison is possible"
}

cat("\n", strrep("=", 74), "\n", sep = "")
cat(sprintf("SECTION 18  THE METABOLIC AXIS  [%s]\n",
            if (METAB_ADJUST_DEPTH) "DEPTH ADJUSTED" else "UNADJUSTED"))
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("Model : y ~ %s\n",
            paste0(MT_MET, "*", MT_COL, " + ", MT_MET, "*", MT_ER, " + ziCAF + zERAD",
                   if (METAB_ADJUST_DEPTH) " + zDepth" else "")))
cat("zMet replaces zHyp as the axis. Everything else is identical to section 13.\n")

cat("\n--- the two interaction terms ---\n")
cat(sprintf("  %-12s %-14s %9s %21s %5s %6s %11s\n",
            "module", "term", "beta", "95% CI", "k", "I2", "p"))
for (vr in METAB_VARIANTS) {
  for (tm in c(MI_COL, MI_ER)) {
    r <- mcoef_pooled[mcoef_pooled$variant == vr & mcoef_pooled$term == tm, ]
    if (!nrow(r)) next
    cat(sprintf("  %-12s %-14s %+9.4f  %+8.4f to %+8.4f %5d %5.0f%% %10s\n",
                vr, tm, r$est, r$lower, r$upper, as.integer(r$k), r$I2, fmt_p(r$p)))
  }
}

cat("\n--- THE HEADLINE: difference between the two slopes, along the metabolic axis ---\n")
for (vr in METAB_VARIANTS) {
  r <- mcon_pooled[mcon_pooled$variant == vr & mcon_pooled$contrast == "interaction_diff", ]
  if (!nrow(r)) next
  cat(sprintf("  %-12s %+.4f  (95%% CI %+.4f to %+.4f)  %s  I2=%3.0f%%  p = %-10s %s\n",
              vr, r$est, r$lower, r$upper, k_lab(r), r$I2, fmt_p(r$p),
              if (is.finite(r$p) && r$p < METAB_ALPHA) "SIGNIFICANT" else "not significant"))
}

cat("\n--- SIDE BY SIDE WITH THE HYPOXIA AXIS, SAME BRCA SECTIONS ---\n")
cat(sprintf("  %s\n", hyp_note))
if (!is.null(hyp_ref) && is.finite(hyp_ref$est)) {
  cat(sprintf("  hypoxia axis   : %+.4f (%+.4f to %+.4f)  p = %s  %s\n",
              hyp_ref$est, hyp_ref$lower, hyp_ref$upper, fmt_p(hyp_ref$p), k_lab(hyp_ref)))
  cat(sprintf("  metabolic axis : %+.4f (%+.4f to %+.4f)  p = %s  %s\n",
              met_common$est, met_common$lower, met_common$upper, fmt_p(met_common$p), k_lab(met_common)))
  same_sign <- sign(met_common$est) == sign(hyp_ref$est)
  sig_met   <- is.finite(met_common$p) && met_common$p < METAB_ALPHA
  sig_hyp   <- is.finite(hyp_ref$p) && hyp_ref$p < METAB_ALPHA
  cat("\n  VERDICT\n")
  if (!sig_hyp) {
    cat("    The hypoxia axis itself shows no significant slope difference on these\n")
    cat("    sections, so there is nothing for the metabolic axis to corroborate.\n")
  } else if (sig_met && same_sign) {
    cat("    The same decoupling appears along BOTH axes, in the same direction.\n")
    cat(sprintf("    How much that counts depends on the gate above: %s.\n",
                gate_verdict[[METAB_PRIMARY]]))
  } else if (sig_met && !same_sign) {
    cat("    *** The metabolic axis is significant but in the OPPOSITE direction.\n")
    cat("    That is neither corroboration nor a clean negative. Report it as it is. ***\n")
  } else {
    cat("    The decoupling does NOT appear along the metabolic axis. Either it is\n")
    cat("    specific to oxygen limitation, or the transcript-level module is too\n")
    cat("    poor a proxy for 2-oxoglutarate to detect it (section 19 tests the\n")
    cat("    module's coherence). Do NOT present this as evidence against the mechanism.\n")
  }
}

cat("\n--- broad set vs focused subset ---\n")
{
  a <- mcon_pooled[mcon_pooled$variant == "Glutamine" & mcon_pooled$contrast == "interaction_diff", ]
  b <- mcon_pooled[mcon_pooled$variant == "OGcore"    & mcon_pooled$contrast == "interaction_diff", ]
  if (nrow(a) == 1 && nrow(b) == 1) {
    cat(sprintf("  71-gene GOBP set : %+.4f  p = %s  %s\n", a$est, fmt_p(a$p), k_lab(a)))
    cat(sprintf("  12-gene 2-OG core: %+.4f  p = %s  %s\n", b$est, fmt_p(b$p), k_lab(b)))
    cat("  The subset was chosen by hand, so it is suggestive only. Report the\n")
    cat("  GOBP set as primary.\n")
  }
}

cat("\nCAVEATS TO STATE BEFORE ANYONE ASKS\n")
cat("  1. The GO term covers arginine, glutamate, glutamine AND proline, so it\n")
cat("     includes urea-cycle, proline and NOS genes that do not feed 2-OG.\n")
cat("  2. Transcript abundance is a weak proxy for metabolite concentration.\n")
cat("  3. Spatial autocorrelation not modelled, so all CIs are too narrow.\n")
cat("  4. Per 1 SD within section, not comparable to sections 10 and 11.\n")
cat(strrep("=", 74), "\n", sep = "")

## ============================================================================
## 18.7  FIGURE 9
##       Panel 1: the hypoxia axis (section 13), re-pooled on the SAME BRCA
##       sections. Panels 2 and 3: the two metabolic modules. One cancer, one
##       set of sections, one model.
## ============================================================================
if (!exists("FIG_STYLES") || !exists("theme_fig")) {
  cat("\n  Figure 9 skipped: section 12 has not run.\n")
} else {
  ld <- mcur_pooled[mcur_pooled$series != LAB_DIF, c("variant", "series", "z", "OR", "lower_OR", "upper_OR")]
  if (!is.null(hyp_curve))
    ld <- rbind(as.data.frame(hyp_curve[hyp_curve$series != LAB_DIF,
                                        c("variant", "series", "z", "OR", "lower_OR", "upper_OR")]),
                as.data.frame(ld))
  ld <- ld[ld$variant %in% names(NICE_AXIS), ]
  lv <- unname(NICE_AXIS[intersect(names(NICE_AXIS), unique(ld$variant))])
  ld$series  <- factor(ld$series, c(LAB_COL, LAB_ER))
  ld$variant <- factor(unname(NICE_AXIS[ld$variant]), levels = lv)
  
  lab_rows <- list()
  if (!is.null(hyp_ref) && is.finite(hyp_ref$est))
    lab_rows[["Hypoxia"]] <- data.frame(variant = NICE_AXIS[["Hypoxia"]],
                                        lab = sprintf("diff %+.3f, p = %s\n%s", hyp_ref$est, fmt_p(hyp_ref$p), k_lab(hyp_ref)))
  for (vr in METAB_VARIANTS) {
    r <- mcon_pooled[mcon_pooled$variant == vr & mcon_pooled$contrast == "interaction_diff", ]
    if (nrow(r)) lab_rows[[vr]] <- data.frame(variant = NICE_AXIS[[vr]],
                                              lab = sprintf("diff %+.3f, p = %s\n%s", r$est, fmt_p(r$p), k_lab(r)))
  }
  lab_df <- do.call(rbind, lab_rows)
  if (!is.null(lab_df)) lab_df$variant <- factor(lab_df$variant, levels = lv)
  
  y_top <- max(ld$upper_OR, na.rm = TRUE); y_bot <- min(ld$lower_OR, na.rm = TRUE)
  y_hi  <- y_top * (y_top / y_bot)^0.28
  med_rho <- if (nrow(pc)) unit_median(pc$rho_vs_Hypoxia, pc$sample_id) else NA_real_
  
  for (sn in names(FIG_STYLES)) {
    S <- FIG_STYLES[[sn]]
    p9 <- ggplot(ld, aes(z, OR)) +
      geom_hline(yintercept = 1, linetype = "22", linewidth = S$lwd * 0.7, colour = "grey40") +
      geom_ribbon(aes(ymin = lower_OR, ymax = upper_OR, fill = series), alpha = .16, colour = NA) +
      geom_line(aes(colour = series), linewidth = S$lwd * 1.3) +
      scale_colour_manual(values = setNames(c(CLR$neg, CLR$pos), c(LAB_COL, LAB_ER)), name = NULL) +
      scale_fill_manual(values   = setNames(c(CLR$neg, CLR$pos), c(LAB_COL, LAB_ER)), name = NULL) +
      facet_wrap(~ variant, nrow = 1) + scale_y_log10() +
      coord_cartesian(ylim = c(y_bot, y_hi)) +
      labs(x = "Axis module score (SD within section)",
           y = "Odds ratio for P4HA1 detection (95% CI)",
           caption = sprintf(paste0("BRCA only. Same model in every panel; only the axis module differs.\n",
                                    "Hypoxia panel re-pooled on the sections of the glutamine panel.\n",
                                    "Median correlation of the glutamine module with hypoxia: %+.2f."),
                             med_rho)) +
      theme_fig(S) + theme(legend.position = "top")
    
    if (!is.null(lab_df) && nrow(lab_df))
      p9 <- p9 + geom_text(data = lab_df, inherit.aes = FALSE,
                           aes(x = min(ld$z), y = y_hi, label = lab),
                           hjust = 0, vjust = 1, size = S$lab, colour = "grey20", lineheight = 1.1)
    
    if (isTRUE(S$titles))
      p9 <- p9 + labs(title = "Does the same decoupling appear along the cofactor axis?",
                      subtitle = "2-oxoglutarate is a substrate of P4H alongside oxygen")
    
    save_fig(p9, S, sub("\\.csv$", "", outf("Figure9_metabolic_axis")),
             S$w2, if (isTRUE(S$split)) S$h_full else MM(72))
  }
}

log_msg("Section 18 complete.")

################################################################################
## 19. IS EACH MODULE INTERNALLY COHERENT?   (Figure 10)
##
##   A module score is a mean over member genes. It only means something if
##   the members move together. This asks a question about the module alone,
##   so it needs no external yardstick (AddModuleScore already subtracts the
##   generic activity such yardsticks measure).
##
##   Two measures: mean pairwise Spearman between member genes, and the share
##   of variance taken by PC1. Visium is sparse, so absolute values are close
##   to meaningless; each module is compared with 30 random gene sets of the
##   same size matched on detection rate, and reported as a PERCENTILE.
##
##   Reference modules are tested alongside, so "coherent" can be judged
##   against modules that are expected to work in THIS cohort. Whether they do
##   is read from this run, not assumed.
##
##   v5 (C4): the summary is a median over patients (each patient = mean of
##   its sections), not over sections.
################################################################################

stopifnot(exists("MODULES"), exists("work"), exists("agg"))
if (!exists("fmt_p") || !exists("assert_run_state") || !exists("write_out"))
  stop("Section 19 needs the 12b helper block.")
if (!"Glutamine" %in% names(MODULES))
  stop("Section 19 needs section 18 to have defined MODULES$Glutamine.")

RATES_CSV <- file.path(agg, "DETECT_cohort_gene_rates.csv")
if (!file.exists(RATES_CSV))
  stop("Section 19 needs DETECT_cohort_gene_rates.csv from section 14.")

## ============================================================================
## 19.0b  CANDIDATE METABOLIC GENE SETS (all standard sets, verbatim)
##   GlnNarrow    GOBP_GLUTAMINE_METABOLIC_PROCESS, M11960, 21 genes. All 21
##                are inside the 71-gene set, so it is a slice of the same
##                instrument. It also carries liver genes (CPS1, FABP1, MTTP,
##                NR1H4) that are not expected in breast tissue; the coverage
##                printed below shows how many are actually detected here.
##   TCA_KEGG     KEGG_CITRATE_CYCLE_TCA_CYCLE, M3985, 32 genes.
##   TCA_Reactome REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE, M784, 34 genes.
##                The two TCA sets share only 18 genes, so if both cohere the
##                conclusion is robust to the definition.
##   OXPHOS       HALLMARK_OXIDATIVE_PHOSPHORYLATION, M5936, 200 genes, minus
##                genes already in another predictor module.
##   All four are clean of P4HA1, P4HA2, P4HA3, P4HB, PLOD1-3 and PPIB.
## ============================================================================
CAND_GLNNARROW <- c(
  "GLYATL1B", "CPS1", "FABP1", "SIRT4", "PHGDH", "BLOC1S6", "GLS2", "GLS",
  "GLUD1", "GLUL", "MECP2", "ASL", "MTTP", "LGSN", "PFAS", "NIT2",
  "UCP2", "CAD", "ARHGAP11B", "GLYATL1", "NR1H4")

CAND_TCA_KEGG <- c(
  "CS", "DLAT", "DLD", "DLST", "FH", "SUCLG2P2", "IDH1", "IDH2",
  "IDH3A", "IDH3B", "IDH3G", "MDH1", "MDH2", "ACLY", "ACO1", "OGDH",
  "ACO2", "PC", "PCK1", "PCK2", "PDHA1", "PDHA2", "PDHB", "OGDHL",
  "SDHA", "SDHB", "SDHC", "SDHD", "SDHCP5", "SUCLG2", "SUCLG1", "SUCLA2")

CAND_TCA_REACTOME <- c(
  "CS", "IDH3G", "SDHA", "ACAT1", "DLD", "FH", "ACO2", "IDH3B",
  "OGDH", "NNT", "SDHB", "DLST", "TRAP1", "KGD4", "ISCA1", "ISCU",
  "SUCLA2", "SIRT3", "SDHC", "MDH2", "SDHAF4", "SUCLG1", "FXN", "ISCA2",
  "IDH3A", "SDHAF2", "SUCLG2", "IDH2", "SDHAF3", "SDHD", "SDHAF1", "LYRM4",
  "CSKMT", "NFS1")

CAND_OXPHOS <- c(
  "ABCB7", "ACAA1", "ACAA2", "ACADM", "ACADSB", "ACADVL", "ACAT1", "ACO2",
  "AFG3L2", "AIFM1", "ALAS1", "ALDH6A1", "ATP1B1", "ATP5F1A", "ATP5F1B", "ATP5F1C",
  "ATP5F1D", "ATP5F1E", "ATP5PB", "ATP5MC1", "ATP5MC2", "ATP5MC3", "ATP5PD", "ATP5ME",
  "ATP5PF", "ATP5MF", "ATP5MG", "ATP5PO", "ATP6AP1", "ATP6V0B", "ATP6V0C", "ATP6V0E1",
  "ATP6V1C1", "ATP6V1D", "ATP6V1E1", "ATP6V1F", "ATP6V1G1", "ATP6V1H", "BAX", "BCKDHA",
  "BDH2", "MPC1", "CASP7", "COX10", "COX11", "COX15", "COX17", "COX4I1",
  "COX5A", "COX5B", "COX6A1", "COX6B1", "COX6C", "COX7A2", "COX7A2L", "COX7B",
  "COX7C", "COX8A", "CPT1A", "CS", "CYB5A", "CYB5R3", "CYC1", "CYCS",
  "DECR1", "DLAT", "DLD", "DLST", "ECH1", "ECHS1", "ECI1", "ETFA",
  "ETFB", "ETFDH", "FDX1", "FH", "FXN", "GLUD1", "GOT2", "GPI",
  "GPX4", "GRPEL1", "HADHA", "HADHB", "HCCS", "HSD17B10", "HSPA9", "HTRA2",
  "IDH1", "IDH2", "IDH3A", "IDH3B", "IDH3G", "IMMT", "ISCA1", "ISCU",
  "LDHA", "LDHB", "LRPPRC", "MAOB", "MDH1", "MDH2", "MFN2", "MGST3",
  "MRPL11", "MRPL15", "MRPL34", "MRPL35", "MRPS11", "MRPS12", "MRPS15", "MRPS22",
  "MRPS30", "MTRF1", "MTRR", "MTX2", "NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4",
  "NDUFA5", "NDUFA6", "NDUFA7", "NDUFA8", "NDUFA9", "NDUFAB1", "NDUFB1", "NDUFB2",
  "NDUFB3", "NDUFB4", "NDUFB5", "NDUFB6", "NDUFB7", "NDUFB8", "NDUFC1", "NDUFC2",
  "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS6", "NDUFS7", "NDUFS8", "NDUFV1",
  "NDUFV2", "NNT", "NQO2", "OAT", "OGDH", "OPA1", "OXA1L", "PDHA1",
  "PDHB", "PDHX", "PDK4", "PDP1", "PHB2", "PHYH", "PMPCA", "POLR2F",
  "POR", "PRDX3", "RETSAT", "RHOT1", "RHOT2", "SDHA", "SDHB", "SDHC",
  "SDHD", "SLC25A11", "SLC25A12", "SLC25A20", "SLC25A3", "SLC25A4", "SLC25A5", "SLC25A6",
  "SUCLA2", "SUCLG1", "SUPV3L1", "SURF1", "TCIRG1", "TIMM10", "TIMM13", "TIMM17A",
  "TIMM50", "TIMM8B", "TIMM9", "TOMM22", "TOMM70", "UQCR10", "UQCR11", "UQCRB",
  "UQCRC1", "UQCRC2", "UQCRFS1", "UQCRH", "UQCRQ", "VDAC1", "VDAC2", "VDAC3")

## genes already in another predictor module are removed (LDHA, HSPA9)
CAND_OXPHOS <- setdiff(CAND_OXPHOS,
                       unique(unlist(MODULES[c("Hypoxia", "CollagenSynthesis", "iCAF",
                                               "ERstressStrict", "ERAD")], use.names = FALSE)))

MODULES$GlnNarrow    <- CAND_GLNNARROW
MODULES$TCA_KEGG     <- CAND_TCA_KEGG
MODULES$TCA_Reactome <- CAND_TCA_REACTOME
MODULES$OXPHOS       <- CAND_OXPHOS

cat("\n  candidate sets added:\n")
for (nm in c("GlnNarrow", "TCA_KEGG", "TCA_Reactome", "OXPHOS")) {
  bad <- intersect(MODULES[[nm]], c("P4HA1","P4HA2","P4HA3","P4HB","PLOD1","PLOD2","PLOD3","PPIB"))
  if (length(bad)) stop(sprintf("Section 19: %s contains %s", nm, paste(bad, collapse=", ")))
  cat(sprintf("    %-14s %3d genes\n", nm, length(MODULES[[nm]])))
}

## ============================================================================
## 19.0  SETTINGS
## ============================================================================
COH_REFERENCES <- c("Hypoxia", "CollagenSynthesis")
COH_MODULES <- c("Hypoxia", "CollagenSynthesis", "ERstressStrict", "ERAD",
                 "Glutamine", "OGcore",
                 "GlnNarrow", "TCA_KEGG", "TCA_Reactome", "OXPHOS")

COH_N_NULL     <- 30
COH_MIN_GENES  <- 4
COH_MIN_SPOTS  <- 30
COH_SEED       <- 90125
COH_CALIPER    <- 0.35

assert_run_state("Section 19")

## ============================================================================
## 19.1  NULL POOL, MATCHED ON DETECTION RATE
## ============================================================================
rates <- read.csv(RATES_CSV, stringsAsFactors = FALSE)
## gene pool from the sections where P4HA1 can be modelled, as in 16.1
if (exists("RND_SECTIONS")) rates <- rates[rates$sample_id %in% RND_SECTIONS, , drop = FALSE]
gene_med <- rates %>%
  group_by(gene) %>%
  summarise(median_rate = median(detect_rate), n_tum = dplyr::n(), .groups = "drop") %>%
  filter(n_tum == dplyr::n_distinct(rates$sample_id))

ALL_MODULE_GENES <- unique(unlist(MODULES[COH_MODULES], use.names = FALSE))
null_pool <- gene_med[!gene_med$gene %in% ALL_MODULE_GENES, ]

cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 19.1  SETUP\n")
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  genes present in all sections : %d\n", nrow(gene_med)))
cat(sprintf("  excluded (module members)     : %d\n", length(ALL_MODULE_GENES)))
cat(sprintf("  null pool                     : %d genes\n", nrow(null_pool)))
cat("\n  modules under test (detection in THIS cohort):\n")
for (nm in COH_MODULES) {
  det <- intersect(MODULES[[nm]], gene_med$gene)
  mr  <- if (length(det)) median(gene_med$median_rate[gene_med$gene %in% det]) else NA_real_
  cat(sprintf("    %-18s %3d genes, %3d detected cohort-wide, median rate %.3f\n",
              nm, length(MODULES[[nm]]), length(det), mr))
}

draw_null_sets <- function(module_genes, n_sets) {
  det <- intersect(module_genes, gene_med$gene)
  if (length(det) < COH_MIN_GENES) return(NULL)
  target <- median(gene_med$median_rate[gene_med$gene %in% det])
  elig <- null_pool$gene[null_pool$median_rate >= target * (1 - COH_CALIPER) &
                           null_pool$median_rate <= target * (1 + COH_CALIPER)]
  if (length(elig) < length(det) * 2) return(NULL)
  lapply(seq_len(n_sets), function(i) sample(elig, length(det)))
}

set.seed(COH_SEED)
NULL_SETS <- lapply(COH_MODULES, function(nm) draw_null_sets(MODULES[[nm]], COH_N_NULL))
names(NULL_SETS) <- COH_MODULES

## v5.1 FALLBACK. The objects store a filtered feature list (~7,500 genes), so a
## module gene can be absent from some sections although it is expressed where
## present (LDHA: 33 of 48 sections, detection 0.97; CA9: 9 of 48). The strict
## rule above (present in EVERY section) then leaves Hypoxia with 3 genes and no
## null. For modules that failed only that way, the target rate and the set size
## come from member genes present in >= half of the sections. The null genes are
## still drawn from the all-sections pool. A separate seed keeps the null sets of
## every other module identical to the strict run.
COH_FALLBACK <- character(0)
gene_any <- rates %>%
  group_by(gene) %>%
  summarise(median_rate = median(detect_rate), n_tum = dplyr::n(), .groups = "drop") %>%
  filter(n_tum >= 0.5 * dplyr::n_distinct(rates$sample_id))
set.seed(COH_SEED + 1L)
for (nm in COH_MODULES[vapply(NULL_SETS, is.null, logical(1))]) {
  det <- intersect(MODULES[[nm]], gene_any$gene)
  if (length(det) < COH_MIN_GENES) next
  target <- median(gene_any$median_rate[gene_any$gene %in% det])
  elig <- null_pool$gene[null_pool$median_rate >= target * (1 - COH_CALIPER) &
                           null_pool$median_rate <= target * (1 + COH_CALIPER)]
  if (length(elig) < length(det) * 2) next
  NULL_SETS[[nm]] <- lapply(seq_len(COH_N_NULL), function(i) sample(elig, length(det)))
  COH_FALLBACK <- c(COH_FALLBACK, nm)
}

cat("\n  null sets drawn:\n")
for (nm in COH_MODULES)
  cat(sprintf("    %-18s %s%s\n", nm,
              if (is.null(NULL_SETS[[nm]])) "*** could not match, module NOT TESTED ***"
              else sprintf("%d sets of %d genes", length(NULL_SETS[[nm]]), length(NULL_SETS[[nm]][[1]])),
              if (nm %in% COH_FALLBACK) "   (fallback: genes in >= half of the sections)" else ""))

## ============================================================================
## 19.2  THE TWO COHERENCE MEASURES
## ============================================================================
coh_meanrho <- function(mat) {
  if (is.null(mat) || nrow(mat) < COH_MIN_GENES) return(NA_real_)
  keep <- apply(mat, 1, function(v) sd(v) > 0)
  if (sum(keep) < COH_MIN_GENES) return(NA_real_)
  m <- mat[keep, , drop = FALSE]
  cc <- suppressWarnings(cor(t(as.matrix(m)), method = "spearman"))
  if (!is.matrix(cc)) return(NA_real_)
  mean(cc[upper.tri(cc)], na.rm = TRUE)
}

coh_pc1 <- function(mat) {
  if (is.null(mat) || nrow(mat) < COH_MIN_GENES) return(NA_real_)
  keep <- apply(mat, 1, function(v) sd(v) > 0)
  if (sum(keep) < COH_MIN_GENES) return(NA_real_)
  m <- scale(t(as.matrix(mat[keep, , drop = FALSE])))
  m <- m[, apply(m, 2, function(v) all(is.finite(v))), drop = FALSE]
  if (ncol(m) < COH_MIN_GENES) return(NA_real_)
  ev <- tryCatch(svd(m, nu = 0, nv = 0)$d^2, error = function(e) NULL)
  if (is.null(ev) || !length(ev)) return(NA_real_)
  ev[1] / sum(ev)
}

## ============================================================================
## 19.3  RUN, ONE SECTION AT A TIME
## ============================================================================
coh_one <- function(rds_path, sid) {
  obj <- readRDS(rds_path)
  obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
  DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else Assays(obj)[1]
  if (ncol(obj) < COH_MIN_SPOTS) return(NULL)
  
  E <- tryCatch(GetAssayData(obj, assay = DefaultAssay(obj), layer = "data"),
                error = function(e) NULL)
  if (is.null(E) || !nrow(E)) return(NULL)
  rn <- rownames(E)
  
  grab <- function(genes) {
    g <- intersect(genes, rn)
    if (length(g) < COH_MIN_GENES) return(NULL)
    E[g, , drop = FALSE]
  }
  
  rows <- list()
  for (nm in COH_MODULES) {
    real <- grab(MODULES[[nm]])
    if (is.null(real)) {
      log_msg("  %s: %s has too few genes here, skipped", sid, nm); next
    }
    r_rho <- coh_meanrho(real); r_pc1 <- coh_pc1(real)
    
    ns <- NULL_SETS[[nm]]
    null_rho <- null_pc1 <- rep(NA_real_, if (is.null(ns)) 0 else length(ns))
    if (!is.null(ns)) {
      for (j in seq_along(ns)) {
        nm_mat <- grab(ns[[j]])
        if (is.null(nm_mat)) next
        null_rho[j] <- coh_meanrho(nm_mat)
        null_pc1[j] <- coh_pc1(nm_mat)
      }
    }
    
    pct_rho <- if (sum(is.finite(null_rho)) >= 10)
      mean(null_rho[is.finite(null_rho)] < r_rho) else NA_real_
    pct_pc1 <- if (sum(is.finite(null_pc1)) >= 10)
      mean(null_pc1[is.finite(null_pc1)] < r_pc1) else NA_real_
    
    rows[[length(rows) + 1]] <- data.frame(
      sample_id = sid, module = nm,
      n_genes_used = nrow(real),
      mean_rho = r_rho, pc1_var = r_pc1,
      null_mean_rho = if (length(null_rho)) median(null_rho, na.rm = TRUE) else NA_real_,
      null_pc1 = if (length(null_pc1)) median(null_pc1, na.rm = TRUE) else NA_real_,
      pct_rho = pct_rho, pct_pc1 = pct_pc1,
      n_null_ok = sum(is.finite(null_rho)),
      row.names = NULL)
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

log_msg("Section 19 start. %d modules, %d null sets each.", length(COH_MODULES), COH_N_NULL)

coh_acc <- list()
for (i in seq_len(nrow(work))) {
  path <- work$path[i]; sid <- work$sample_id[i]
  if (is.na(path) || !file.exists(path)) next
  log_msg("[coh %d/%d] %s", i, nrow(work), sid)
  r <- tryCatch(coh_one(path, sid),
                error = function(e) { log_msg("  coh skipped: %s", conditionMessage(e)); NULL })
  if (!is.null(r) && nrow(r)) coh_acc[[length(coh_acc) + 1]] <- r
}
coh_all <- bind(coh_acc)
if (nrow(coh_all) == 0) stop("Section 19: no section produced coherence estimates.")
write_out(coh_all, "COH_per_sample")

## ============================================================================
## 19.4  SUMMARISE  (one value per unit, then the median over units)
## ============================================================================
COH_COLS <- c("n_genes_used", "mean_rho", "null_mean_rho", "pct_rho", "pc1_var", "null_pc1", "pct_pc1")
coh_unit <- coh_all %>%
  mutate(unit = unit_of(sample_id)) %>%
  group_by(module, unit) %>%
  summarise(across(all_of(COH_COLS), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

coh_sum <- coh_unit %>%
  group_by(module) %>%
  summarise(k = dplyr::n(),
            median_genes = median(n_genes_used),
            median_rho = median(mean_rho, na.rm = TRUE),
            median_null_rho = median(null_mean_rho, na.rm = TRUE),
            median_pct_rho = median(pct_rho, na.rm = TRUE),
            median_pc1 = median(pc1_var, na.rm = TRUE),
            median_null_pc1 = median(null_pc1, na.rm = TRUE),
            median_pct_pc1 = median(pct_pc1, na.rm = TRUE),
            .groups = "drop")
coh_sum <- coh_sum[match(COH_MODULES, coh_sum$module), ]
coh_sum <- coh_sum[!is.na(coh_sum$module), ]
write_out(coh_sum, "COH_summary")

cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 19  DO THE GENES INSIDE EACH MODULE BEHAVE AS A UNIT?\n")
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  Medians over %ss. The percentile against the module's own null is the\n", UNIT_LABEL))
cat("  number that matters.\n\n")

cat(sprintf("  %-18s %4s %6s %8s %9s %9s %8s %9s\n",
            "module", "k", "genes", "rho", "null rho", "pctile", "PC1", "pctile"))
for (i in seq_len(nrow(coh_sum))) with(coh_sum[i, ],
                                       cat(sprintf("  %-18s %4d %6.0f %8.3f %9.3f %8.0f%% %8.3f %8.0f%%\n",
                                                   module, as.integer(k), median_genes, median_rho, median_null_rho, 100*median_pct_rho,
                                                   median_pc1, 100*median_pct_pc1)))

## ============================================================================
## 19.5  VERDICT
## ============================================================================
COH_PCT_MIN <- 0.90

passes <- function(mod) {
  r <- coh_sum[coh_sum$module == mod, ]
  if (!nrow(r)) return(FALSE)
  isTRUE(r$median_pct_rho >= COH_PCT_MIN) || isTRUE(r$median_pct_pc1 >= COH_PCT_MIN)
}
## v5.1: a module without a matched null was never tested; say so instead of
## reporting it as "not coherent".
coh_tested <- function(mod) {
  r <- coh_sum[coh_sum$module == mod, ]
  nrow(r) > 0 && (is.finite(r$median_pct_rho) || is.finite(r$median_pct_pc1))
}
coh_label <- function(mod, caps = TRUE) {
  if (!coh_tested(mod)) return("NOT TESTED (no matched null)")
  if (passes(mod)) (if (caps) "COHERENT" else "coherent")
  else (if (caps) "not coherent" else "NOT coherent")
}

cat("\n", strrep("-", 74), "\n", sep = "")
cat(sprintf("  criterion: median percentile >= %.0f%% on either measure\n\n", 100*COH_PCT_MIN))
for (nm in coh_sum$module)
  cat(sprintf("    %-18s %s%s\n", nm, coh_label(nm),
              if (exists("COH_FALLBACK") && nm %in% COH_FALLBACK) "   [fallback null]" else ""))

hyp_ok  <- passes("Hypoxia")
col_ok  <- passes("CollagenSynthesis")
glut_ok <- passes("Glutamine")
og_ok   <- passes("OGcore")
ctrl_ok <- hyp_ok || col_ok

cat("\n", strrep("=", 74), "\n", sep = "")
cat("VERDICT\n")
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  Reference modules in this cohort: hypoxia %s, collagen synthesis %s.\n",
            coh_label("Hypoxia", caps = FALSE), coh_label("CollagenSynthesis", caps = FALSE)))
if (coh_tested("Hypoxia") && !hyp_ok) {
  cat("  NOTE: the hypoxia module is the axis of sections 10-16. If it is not\n")
  cat("  coherent, the hypoxia score is a composite of its genes, not one\n")
  cat("  coordinated programme; read the main result with that in mind.\n")
}

if (!ctrl_ok) {
  cat("\n  NEITHER REFERENCE MODULE IS COHERENT, SO THIS SECTION PROVES NOTHING.\n")
  cat("  The measure does not separate known modules from random genes in this\n")
  cat("  data. Do not conclude anything about the metabolic modules from it.\n")
} else if (!glut_ok && !og_ok) {
  cat("\n  THE GLUTAMINE MODULE IS NOT A COHERENT UNIT IN THIS COHORT.\n")
  cat("  Its member genes do not move together, so the score is an average of\n")
  cat("  unrelated quantities. Section 18's result for it says something about\n")
  cat("  the instrument, not the biology: the cofactor hypothesis has not been\n")
  cat("  tested by it, in either direction.\n")
} else {
  cat(sprintf("\n  THE GLUTAMINE MODULE IS COHERENT (71-gene set %s, 12-gene core %s).\n",
              if (glut_ok) "coherent" else "not coherent", if (og_ok) "coherent" else "not coherent"))
  cat("  Its score represents something, so section 18's result cannot be\n")
  cat("  dismissed as an incoherent instrument. Two readings remain: the effect\n")
  cat("  is specific to oxygen, or the module tracks transcriptional\n")
  cat("  co-regulation rather than metabolic flux.\n")
}

CANDIDATES <- c("GlnNarrow", "TCA_KEGG", "TCA_Reactome", "OXPHOS")
cand_ok <- CANDIDATES[vapply(CANDIDATES, passes, logical(1))]

cat("\n", strrep("-", 74), "\n", sep = "")
cat("  CANDIDATE SETS: is there a usable metabolic instrument at all?\n\n")
for (nm in intersect(CANDIDATES, coh_sum$module))
  cat(sprintf("    %-14s %s\n", nm, if (passes(nm)) "COHERENT" else "not coherent"))

if (!ctrl_ok) {
  cat("\n    References failed, so ignore the above.\n")
} else if (!length(cand_ok)) {
  cat("\n    None of the candidates is coherent either, while at least one\n")
  cat("    reference module is. Metabolic genes do not move together in this\n")
  cat("    data, consistent with flux being controlled mostly after\n")
  cat("    transcription. The cofactor axis is then not testable with these\n")
  cat("    Visium transcript data. Say that plainly.\n")
} else {
  cat(sprintf("\n    %d candidate(s) coherent: %s\n", length(cand_ok),
              paste(cand_ok, collapse = ", ")))
  cat("    Before rerunning section 18 on one, check its correlation with the\n")
  cat("    hypoxia module and whether it plausibly reports on 2-OG supply\n")
  cat("    (OXPHOS is downstream of the TCA cycle). Pick ONE on biological grounds.\n")
}

cat("\n  CAVEATS\n")
cat("    Coherence is necessary for a module score to mean one thing. It is not\n")
cat("    sufficient: a coherent module can still measure the wrong thing.\n")
cat(strrep("=", 74), "\n", sep = "")

## ============================================================================
## 19.6  FIGURE 10  (one point per unit)
## ============================================================================
if (!exists("FIG_STYLES") || !exists("theme_fig")) {
  cat("\n  Figure 10 skipped: section 12 has not run.\n")
} else {
  pd <- coh_unit[is.finite(coh_unit$pct_rho), ]
  lv <- intersect(COH_MODULES, unique(pd$module))
  nice_m <- c(Hypoxia = "Hypoxia", CollagenSynthesis = "Collagen synth.",
              ERstressStrict = "ER stress", ERAD = "ERAD",
              Glutamine = "Glutamine (71)", OGcore = "2-OG core (12)",
              GlnNarrow = "Glutamine (21)", TCA_KEGG = "TCA KEGG",
              TCA_Reactome = "TCA Reactome", OXPHOS = "OXPHOS")
  pd$is_metab <- ifelse(pd$module %in% c("Glutamine", "OGcore", CANDIDATES), "metabolic", "reference")
  pd$module <- factor(nice_m[pd$module], levels = unname(nice_m[lv]))
  
  for (sn in names(FIG_STYLES)) {
    S <- FIG_STYLES[[sn]]
    p10 <- ggplot(pd, aes(module, 100 * pct_rho)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 100 * COH_PCT_MIN, ymax = Inf,
               fill = CLR$pos, alpha = 0.07) +
      geom_hline(yintercept = 100 * COH_PCT_MIN, linetype = "22",
                 linewidth = S$lwd * 0.7, colour = "grey40") +
      geom_hline(yintercept = 50, linetype = "solid",
                 linewidth = S$lwd * 0.5, colour = "grey80") +
      geom_boxplot(aes(fill = is_metab), outlier.shape = NA,
                   linewidth = S$lwd * 0.8, width = .55) +
      geom_jitter(width = .13, size = S$pt * 0.5, alpha = .55, colour = "grey25") +
      scale_fill_manual(values = c(reference = "grey92", metabolic = "#F4C9C4"), guide = "none") +
      labs(x = NULL, y = "Percentile within matched random sets",
           caption = sprintf("One point per %s (mean of its sections). Shaded: above %.0f%%, where a coherent module should sit. Grey line: chance.",
                             UNIT_LABEL, 100 * COH_PCT_MIN)) +
      theme_fig(S) +
      theme(panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(linewidth = S$grid, colour = "grey90"),
            axis.text.x = element_text(angle = 20, hjust = 1))
    
    if (isTRUE(S$titles))
      p10 <- p10 + labs(title = "Do the genes inside each module behave as a unit?",
                        subtitle = "Compared against random gene sets of the same size and detection rate")
    
    save_fig(p10, S, sub("\\.csv$", "", outf("Figure10_module_coherence")),
             S$w2, if (isTRUE(S$split)) S$h_full else MM(74))
  }
}

log_msg("Section 19 complete.")

################################################################################
## 20. DOES CELL COMPOSITION EXPLAIN THE RESULT?
##
##   A spot reads several cells. If oxygenated spots are mostly fibroblasts
##   and hypoxic spots mostly tumour cells, the decoupling can appear without
##   any cell changing its behaviour.
##
##   Composition measure: stroma and immune scores REBUILT from the ESTIMATE
##   signatures with every gene shared with a predictor module, an outcome or
##   the hydroxylase family removed (the stromal signature contains COL1A2,
##   COL3A1, COL6A3 and CXCL12). The per-spot ESTIMATE table stored in each
##   object (obj@misc$Tumor_purity_estimate$table) is used only to validate
##   the rebuilt scores; a section whose table is missing or cannot be
##   aligned to its spots is skipped and named in the safety check.
##
##   THREE MODELS, SAME SPOTS, SAME POOLING
##     M0 baseline    : the section 15 model unchanged; must reproduce it.
##     M1 adjusted    : M0 + zStroma + zImmune.
##     M2 low stroma  : M0 in the lowest third of spots by stroma score.
##   Read the size of the slope difference and how much survives, not p.
################################################################################

stopifnot(exists("MODULES"), exists("work"), exists("agg"), exists("SEED"),
          exists("ADJUST_DEPTH"), exists("pool_with_p"), exists("fmt_p"),
          exists("write_out"), exists("assert_run_state"), exists("log_msg"),
          exists("bind"), exists("min_genes_for"), exists("check_against_s15"))

## ============================================================================
## 20.0  SETTINGS
## ============================================================================
SC_SIG_CSV       <- file.path(PROJECT_ROOT, "data", "estimate_audit", "ESTIMATE_signatures_used.csv")
SC_OUTCOMES      <- c("P4HA1", "P4HA2", "P4HA3")
SC_ES            <- "ERstressStrict"
SC_BASE_MODULES  <- c("Hypoxia", "CollagenSynthesis", "iCAF", SC_ES, "ERAD")
SC_NEW_MODULES   <- c("StromaESTc", "ImmuneESTc")
SC_MIN_SPOTS     <- 30
SC_MIN_EVENTS    <- 70          # same as section 15
SC_ADJUST_DEPTH  <- ADJUST_DEPTH
SC_STRAT_Q       <- 1/3
SC_FAMILY        <- c("P4HA1", "P4HA2", "P4HA3", "P4HB", "PLOD1", "PLOD2", "PLOD3", "PPIB")

SC_RENAME <- c(Hypoxia = "zHyp", CollagenSynthesis = "zCol", iCAF = "ziCAF",
               ERAD = "zERAD", StromaESTc = "zStroma", ImmuneESTc = "zImmune")
SC_RENAME[SC_ES] <- "zER"

SC_II_COL <- "zHyp:zCol"
SC_II_ER  <- "zHyp:zER"
SC_BASE_RHS <- "zHyp * zCol + zHyp * zER + ziCAF + zERAD"
SC_MODELS <- c("M0_baseline", "M1_adjusted", "M2_low_stroma")

assert_run_state("Section 20")

## ============================================================================
## 20.1  BUILD THE CLEAN STROMA AND IMMUNE GENE SETS (computed, never typed)
## ============================================================================
if (!file.exists(SC_SIG_CSV))
  stop("Section 20 needs ", SC_SIG_CSV, ". Run estimate_overlap_audit.R first.")
sig <- read.csv(SC_SIG_CSV, stringsAsFactors = FALSE)
STROMA_RAW <- unique(sig$gene[sig$signature == "stromal"])
IMMUNE_RAW <- unique(sig$gene[sig$signature == "immune"])
stopifnot(length(STROMA_RAW) == 141, length(IMMUNE_RAW) == 141)

SC_BLOCK <- unique(c(unlist(MODULES[SC_BASE_MODULES], use.names = FALSE),
                     SC_FAMILY, SC_OUTCOMES))
both_sigs <- intersect(STROMA_RAW, IMMUNE_RAW)
STROMA_CLEAN <- setdiff(STROMA_RAW, c(SC_BLOCK, both_sigs))
IMMUNE_CLEAN <- setdiff(IMMUNE_RAW, c(SC_BLOCK, both_sigs))

SC_MOD <- MODULES
SC_MOD$StromaESTc <- STROMA_CLEAN
SC_MOD$ImmuneESTc <- IMMUNE_CLEAN

cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 20.1  CLEAN COMPOSITION GENE SETS\n")
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  stromal: %d -> %d genes, removed: %s\n", length(STROMA_RAW), length(STROMA_CLEAN),
            paste(sort(setdiff(STROMA_RAW, STROMA_CLEAN)), collapse = ", ")))
cat(sprintf("  immune : %d -> %d genes, removed: %s\n", length(IMMUNE_RAW), length(IMMUNE_CLEAN),
            { r <- sort(setdiff(IMMUNE_RAW, IMMUNE_CLEAN)); if (length(r)) paste(r, collapse = ", ") else "none" }))
if (length(both_sigs))
  cat(sprintf("  genes in BOTH signatures, removed from both: %s\n", paste(both_sigs, collapse = ", ")))

write_out(data.frame(set  = c(rep("StromaESTc", length(STROMA_CLEAN)),
                              rep("ImmuneESTc", length(IMMUNE_CLEAN))),
                     gene = c(STROMA_CLEAN, IMMUNE_CLEAN)),
          "SC_signature_genes_used")

for (og in SC_OUTCOMES) {
  bad <- names(SC_MOD[c(SC_BASE_MODULES, SC_NEW_MODULES)])[vapply(
    SC_MOD[c(SC_BASE_MODULES, SC_NEW_MODULES)], function(g) og %in% g, logical(1))]
  if (length(bad)) stop("Section 20 halted: ", og, " is inside ", paste(bad, collapse = ", "))
}
shared_new <- lapply(SC_BASE_MODULES, function(m) intersect(SC_MOD[[m]], c(STROMA_CLEAN, IMMUNE_CLEAN)))
if (any(lengths(shared_new) > 0)) stop("Section 20 halted: a clean set still shares genes with a module.")
cat("  circularity gate: clean\n")

## ============================================================================
## 20.2  SCORE, STANDARDISE AND VALIDATE, ONE SECTION AT A TIME
##       Scoring follows section 15 exactly, so the five original scores are
##       identical to section 15's.
## ============================================================================
sc_expr_one <- function(rds_path, sid) {
  obj <- readRDS(rds_path)
  obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
  DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else Assays(obj)[1]
  
  if (length(setdiff(SC_OUTCOMES, rownames(obj)))) {
    log_msg("  %s: outcome gene absent -> skipped", sid); return(NULL)
  }
  
  ## v5 (C3): the ESTIMATE table must be alignable to the spots. Before, a
  ## missing meta column made max(abs(numeric(0))) = -Inf and the check passed.
  est <- tryCatch(obj@misc$Tumor_purity_estimate$table, error = function(e) NULL)
  if (is.null(est) || !nrow(est)) { log_msg("  %s: no ESTIMATE table -> skipped", sid); return(NULL) }
  if (!all(c("stromal", "immune") %in% colnames(est))) {
    log_msg("  %s: ESTIMATE table has no stromal/immune columns -> skipped", sid); return(NULL)
  }
  if (!is.null(rownames(est)) && all(colnames(obj) %in% rownames(est))) {
    est <- est[colnames(obj), , drop = FALSE]
  } else {
    meta_pur <- obj@meta.data$Tumor_purity_estimate
    aligned <- !is.null(meta_pur) && "purity" %in% colnames(est) &&
      nrow(est) == ncol(obj) && length(meta_pur) == nrow(est) &&
      isTRUE(max(abs(est$purity - meta_pur), na.rm = FALSE) <= 1e-6)
    if (!aligned) {
      log_msg("  %s: ESTIMATE table cannot be aligned to the spots -> skipped", sid); return(NULL)
    }
  }
  
  for (nm in c(SC_BASE_MODULES, SC_NEW_MODULES)) {
    full <- SC_MOD[[nm]]
    g    <- intersect(full, rownames(obj))
    if (length(g) < min_genes_for(full)) {
      log_msg("  %s: %s has %d/%d genes -> skipped", sid, nm, length(g), length(full))
      return(NULL)
    }
    set.seed(SEED)
    obj <- tryCatch(add_module_score_fallback(obj, features = list(g), name = nm,
                                              assay = DefaultAssay(obj), nbin = 20, ctrl = 50),
                    error = function(e) NULL)
    if (is.null(obj) || !paste0(nm, "1") %in% colnames(obj@meta.data)) return(NULL)
    names(obj@meta.data)[names(obj@meta.data) == paste0(nm, "1")] <- nm
  }
  
  vars <- c(SC_OUTCOMES, SC_BASE_MODULES, SC_NEW_MODULES)
  cnt_col <- grep("^nCount_", colnames(obj@meta.data), value = TRUE)
  if (SC_ADJUST_DEPTH) { if (!length(cnt_col)) return(NULL); vars <- c(vars, cnt_col[1]) }
  d <- FetchData(obj, vars = vars, layer = "data")
  
  zz <- function(x) { s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) return(NULL)
  as.numeric((x - mean(x, na.rm = TRUE)) / s) }
  for (nm in c(SC_BASE_MODULES, SC_NEW_MODULES)) {
    v <- zz(d[[nm]]); if (is.null(v)) { log_msg("  %s: %s zero variance -> skipped", sid, nm); return(NULL) }
    d[[ SC_RENAME[[nm]] ]] <- v
  }
  if (SC_ADJUST_DEPTH) { v <- zz(log1p(d[[cnt_col[1]]])); if (is.null(v)) return(NULL); d$zDepth <- v }
  
  sp <- function(a, b) suppressWarnings(cor(a, b, method = "spearman"))
  attr(d, "validation") <- data.frame(
    sample_id                = sid,
    n_spots                  = nrow(d),
    rho_stroma_vs_ESTIMATE   = sp(d$StromaESTc, est$stromal),
    rho_immune_vs_ESTIMATE   = sp(d$ImmuneESTc, est$immune),
    rho_stroma_vs_collagen   = sp(d$zStroma, d$zCol),
    rho_stroma_vs_hypoxia    = sp(d$zStroma, d$zHyp),
    rho_stroma_vs_ERstress   = sp(d$zStroma, d$zER),
    rho_immune_vs_hypoxia    = sp(d$zImmune, d$zHyp),
    row.names = NULL)
  d
}

## ============================================================================
## 20.3  FIT ONE MODEL FOR ONE OUTCOME
##       M2 keeps the whole-section z-scores, so its contrast is in the same
##       units as M0.
## ============================================================================
fit_sc_one <- function(d, sid, out_gene, model) {
  out <- list(coef = NULL, contrast = NULL, status = NULL)
  
  if (model == "M2_low_stroma") d <- d[d$zStroma <= quantile(d$zStroma, SC_STRAT_Q), , drop = FALSE]
  
  y  <- as.integer(d[[out_gene]] > 0)
  ev <- sum(y); nev <- length(y) - ev
  status <- data.frame(sample_id = sid, outcome = out_gene, model = model,
                       n_spots = length(y), n_events = ev, n_nonevents = nev,
                       fitted = FALSE, has_zDepth = NA, reason = NA_character_, row.names = NULL)
  
  if (length(y) < SC_MIN_SPOTS || ev == 0 || nev == 0) {
    status$reason <- "no variation in outcome"; out$status <- status; return(out) }
  if (ev < SC_MIN_EVENTS || nev < SC_MIN_EVENTS) {
    status$reason <- sprintf("too few events (%d) or non-events (%d), threshold %d", ev, nev, SC_MIN_EVENTS)
    out$status <- status; return(out) }
  
  rhs <- SC_BASE_RHS
  if (model == "M1_adjusted") rhs <- paste(rhs, "+ zStroma + zImmune")
  if (SC_ADJUST_DEPTH)        rhs <- paste(rhs, "+ zDepth")
  
  dd <- d; dd$y <- y
  f <- tryCatch(glm(as.formula(paste("y ~", rhs)), data = dd, family = binomial()),
                error = function(e) NULL)
  if (is.null(f) || !f$converged) { status$reason <- "glm did not converge"; out$status <- status; return(out) }
  
  b <- coef(f); V <- vcov(f)
  has_depth <- "zDepth" %in% names(b)
  if (!identical(has_depth, isTRUE(SC_ADJUST_DEPTH)))
    stop(sprintf("Section 20, %s / %s: ADJUST_DEPTH = %s but the fitted model %s zDepth.",
                 sid, model, SC_ADJUST_DEPTH, if (has_depth) "DOES contain" else "does NOT contain"))
  status$has_zDepth <- has_depth
  
  need <- c("zCol", "zER", SC_II_COL, SC_II_ER)
  if (model == "M1_adjusted") need <- c(need, "zStroma", "zImmune")
  if (!all(need %in% names(b)) || any(!is.finite(b[need]))) {
    status$reason <- "coefficients not estimable"; out$status <- status; return(out) }
  status$fitted <- TRUE; out$status <- status
  
  s  <- summary(f)$coefficients
  tm <- setdiff(rownames(s), "(Intercept)")
  out$coef <- data.frame(sample_id = sid, outcome = out_gene, model = model, term = tm,
                         log_OR = s[tm, "Estimate"], SE = s[tm, "Std. Error"],
                         p = s[tm, "Pr(>|z|)"], row.names = NULL)
  
  Lvec <- function(w) { L <- setNames(numeric(length(b)), names(b)); for (nm in names(w)) L[nm] <- w[[nm]]; L }
  L_int  <- Lvec(setNames(list(1, -1), c(SC_II_COL, SC_II_ER)))
  L_main <- Lvec(setNames(list(1, -1), c("zCol", "zER")))
  out$contrast <- data.frame(
    sample_id = sid, outcome = out_gene, model = model,
    contrast = c("interaction_diff", "main_diff"),
    est = c(as.numeric(L_int %*% b), as.numeric(L_main %*% b)),
    var = c(as.numeric(t(L_int) %*% V %*% L_int), as.numeric(t(L_main) %*% V %*% L_main)),
    n_spots = nrow(dd), n_events = ev, row.names = NULL)
  out
}

## ============================================================================
## 20.4  RUN
## ============================================================================
log_msg("Section 20 start. ADJUST_DEPTH = %s. run_id = %s.", SC_ADJUST_DEPTH, RUN_ID)

sc_acc <- list(coef = list(), contrast = list(), status = list(), validation = list())
for (i in seq_len(nrow(work))) {
  path <- work$path[i]; sid <- work$sample_id[i]
  if (is.na(path) || !file.exists(path)) next
  log_msg("[sc %d/%d] %s", i, nrow(work), sid)
  d <- tryCatch(sc_expr_one(path, sid),
                error = function(e) { log_msg("  sc skipped: %s", conditionMessage(e)); NULL })
  if (is.null(d)) next
  sc_acc$validation[[length(sc_acc$validation) + 1]] <- attr(d, "validation")
  for (og in SC_OUTCOMES) for (md in SC_MODELS) {
    r <- fit_sc_one(d, sid, og, md)
    for (k in c("coef", "contrast", "status"))
      if (!is.null(r[[k]]) && nrow(r[[k]])) sc_acc[[k]][[length(sc_acc[[k]]) + 1]] <- r[[k]]
  }
}

sc_coef <- bind(sc_acc$coef)
sc_con  <- bind(sc_acc$contrast)
sc_stat <- bind(sc_acc$status)
sc_val  <- bind(sc_acc$validation)
if (nrow(sc_con) == 0) stop("Section 20: no model produced a usable fit.")

write_out(sc_val,  "SC_validation_per_sample")
write_out(sc_coef, "SC_coefficients_per_sample")
write_out(sc_con,  "SC_contrasts_per_sample")
write_out(sc_stat, "SC_fit_status")

## ============================================================================
## 20.5  POOL
## ============================================================================
sc_coef_pooled <- sc_coef %>%
  group_by(outcome, model, term) %>% group_modify(~ pool_with_p(.x$log_OR, .x$SE^2, .x$sample_id)) %>% ungroup()
write_out(sc_coef_pooled, "SC_coefficients_pooled")

sc_con_pooled <- sc_con %>%
  group_by(outcome, model, contrast) %>% group_modify(~ pool_with_p(.x$est, .x$var, .x$sample_id)) %>% ungroup()
write_out(sc_con_pooled, "SC_slope_difference_pooled")

## ============================================================================
## 20.6  SAFETY CHECK: M0 MUST REPRODUCE SECTION 15, SECTION BY SECTION
## ============================================================================
m0_est <- setNames(rep(NA_real_, length(SC_OUTCOMES)), SC_OUTCOMES)
for (og in SC_OUTCOMES) {
  m0_ps <- sc_con[sc_con$outcome == og & sc_con$model == "M0_baseline" &
                    sc_con$contrast == "interaction_diff", ]
  if (!nrow(m0_ps)) { cat(sprintf("\n  %s: no M0 fit, outcome left out of section 20.\n", og)); next }
  check_against_s15(m0_ps, og, "Section 20")
  m0_est[og] <- sc_con_pooled$est[sc_con_pooled$outcome == og & sc_con_pooled$model == "M0_baseline" &
                                    sc_con_pooled$contrast == "interaction_diff"]
}

## ============================================================================
## 20.7  REPORT
## ============================================================================
cat("\n", strrep("=", 74), "\n", sep = "")
cat(sprintf("SECTION 20  [%s]  run_id %s\n", if (SC_ADJUST_DEPTH) "DEPTH ADJUSTED" else "UNADJUSTED", RUN_ID))
cat(strrep("=", 74), "\n", sep = "")

um <- function(x) unit_median(x, sc_val$sample_id)
cat(sprintf("\n--- is the rebuilt composition measure the same thing as ESTIMATE? (median over %ss) ---\n", UNIT_LABEL))
cat(sprintf("  stroma vs ESTIMATE stromal : median rho %.2f (min %.2f)\n",
            um(sc_val$rho_stroma_vs_ESTIMATE), min(sc_val$rho_stroma_vs_ESTIMATE, na.rm = TRUE)))
cat(sprintf("  immune vs ESTIMATE immune  : median rho %.2f (min %.2f)\n",
            um(sc_val$rho_immune_vs_ESTIMATE), min(sc_val$rho_immune_vs_ESTIMATE, na.rm = TRUE)))

cat("\n--- how strongly stroma travels with the model's predictors ---\n")
cat(sprintf("  stroma vs collagen synthesis : median rho %+.2f\n", um(sc_val$rho_stroma_vs_collagen)))
cat(sprintf("  stroma vs hypoxia            : median rho %+.2f\n", um(sc_val$rho_stroma_vs_hypoxia)))
cat(sprintf("  stroma vs ER stress          : median rho %+.2f\n", um(sc_val$rho_stroma_vs_ERstress)))
cat(sprintf("  immune vs hypoxia            : median rho %+.2f\n", um(sc_val$rho_immune_vs_hypoxia)))

for (og in SC_OUTCOMES) {
  if (!is.finite(m0_est[og])) next
  cat("\n", strrep("-", 74), "\n", sep = "")
  cat(sprintf("  OUTCOME %s\n", og))
  cat(strrep("-", 74), "\n", sep = "")
  
  cat("\n  sections per model\n")
  for (md in SC_MODELS) {
    s <- sc_stat[sc_stat$outcome == og & sc_stat$model == md, ]
    cat(sprintf("    %-14s fitted in %2d / %2d   median spots %5.0f   median events %5.0f\n",
                md, sum(s$fitted), nrow(s),
                if (any(s$fitted)) median(s$n_spots[s$fitted]) else NA_real_,
                if (any(s$fitted)) median(s$n_events[s$fitted]) else NA_real_))
    drop <- s[!s$fitted, ]
    for (j in seq_len(nrow(drop)))
      cat(sprintf("                   excluded: %-30s %s\n", drop$sample_id[j], drop$reason[j]))
  }
  
  cat("\n  THE HEADLINE: difference between the two slopes\n")
  for (md in SC_MODELS) {
    r <- sc_con_pooled[sc_con_pooled$outcome == og & sc_con_pooled$model == md &
                         sc_con_pooled$contrast == "interaction_diff", ]
    if (!nrow(r)) next
    cat(sprintf("    %-14s %+.4f  (95%% CI %+.4f to %+.4f)  %s  I2=%3.0f%%  p = %-10s  retained %4.0f%% of M0\n",
                md, r$est, r$lower, r$upper, k_lab(r), r$I2, fmt_p(r$p), 100 * r$est / m0_est[og]))
  }
  
  cat("\n  the individual terms\n")
  cat(sprintf("    %-10s %-14s %9s %19s %6s %11s\n", "term", "model", "beta", "95% CI", "I2", "p"))
  for (tm in c("zCol", "zER", SC_II_COL, SC_II_ER, "zStroma", "zImmune")) {
    for (md in SC_MODELS) {
      r <- sc_coef_pooled[sc_coef_pooled$outcome == og & sc_coef_pooled$model == md &
                            sc_coef_pooled$term == tm, ]
      if (!nrow(r)) next
      cat(sprintf("    %-10s %-14s %+9.4f  %+8.4f to %+8.4f %5.0f%% %10s\n",
                  tm, md, r$est, r$lower, r$upper, r$I2, fmt_p(r$p)))
    }
  }
}

cat("\nCAVEATS\n")
cat("  Model p values are inflated by spatial autocorrelation; read effect sizes.\n")
cat(sprintf("  k counts %s (POOL_LEVEL = %s); section 21 gives the other unit.\n", K_UNIT, POOL_LEVEL))
cat("  ESTIMATE signatures were built for bulk tumours, not single spots.\n")
cat(sprintf("  Stroma and collagen synthesis travel together here (median rho %+.2f), so\n",
            um(sc_val$rho_stroma_vs_collagen)))
cat("  M1 asks a narrower question than 'does composition explain it'. M2 is cleaner.\n")
cat("  M2 may lose sections (fewer events); compare k before comparing estimates.\n")
cat("  Compare models WITHIN a gene, not across genes with different detection.\n")
cat(strrep("=", 74), "\n", sep = "")

log_msg("Section 20 complete.")

################################################################################
## 21. SECTION-LEVEL SENSITIVITY
##
##   The pipeline pools over patients (POOL_LEVEL). This section reports the
##   other direction: the result with sections as the unit.
##   1. Reproduces the pipeline's pooled results from the per-section CSVs
##      with its own code; stops if it cannot.
##   2. Pools the same per-section estimates with sections as the unit.
##   3. With metafor: a three-level model (sections nested in patients).
##   4. Repeats the patient pooling treating HT339B1 and HT339B2 as one
##      person (different HTAN participant IDs, shared HT339 prefix; not
##      confirmed whether they are the same person).
##   No model is refitted. Standalone: needs only the aggregate CSVs.
################################################################################

if (!exists("OUT_ROOT")) {
  PROJECT_ROOT       <- Sys.getenv("P4HA1_BRCA_ROOT", unset = path.expand("~/p4ha1_brca"))
  INCLUDE_METASTASIS <- FALSE      # keep in step with section 0
  POOL_LEVEL         <- "patient"  # keep in step with section 0
  OUT_ROOT <- file.path(PROJECT_ROOT, paste0(
    if (INCLUDE_METASTASIS) "cohort_with_metastasis" else "cohort_primary",
    if (POOL_LEVEL == "section") "_sectionpool" else ""))
}
AGG_DIR <- file.path(OUT_ROOT, "aggregate")
OUT_DIR <- AGG_DIR
TOL <- 1e-6

if (!dir.exists(AGG_DIR)) stop("AGG_DIR not found: ", AGG_DIR, ". Set P4HA1_BRCA_ROOT.")

## ============================================================================
## 21.0  READ (unadjusted rows only, one run_id per file)
## ============================================================================
read_one <- function(name, required = TRUE) {
  f <- file.path(AGG_DIR, paste0(name, ".csv"))
  if (!file.exists(f)) { if (required) stop("Missing file: ", f) else return(NULL) }
  d <- read.csv(f, stringsAsFactors = FALSE)
  if ("adjust_depth" %in% names(d)) d <- d[!as.logical(d$adjust_depth), , drop = FALSE]
  if ("run_id" %in% names(d)) {
    ids <- unique(d$run_id)
    if (length(ids) != 1) stop(name, " holds several run_ids: ", paste(ids, collapse = ", "))
    cat(sprintf("  %-34s run_id %s, %d rows\n", name, ids, nrow(d)))
  } else {
    cat(sprintf("  %-34s (no run_id column), %d rows\n", name, nrow(d)))
  }
  d
}

cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 21  SECTION-LEVEL SENSITIVITY\n")
cat(strrep("=", 74), "\n", sep = "")
cat("\nInputs:\n")
s21_iso    <- read_one("ISOFORM_contrasts_per_sample")
s21_iso_pl <- read_one("ISOFORM_slope_difference_pooled")
s21_sc     <- read_one("SC_contrasts_per_sample", required = FALSE)
s21_sc_pl  <- read_one("SC_slope_difference_pooled", required = FALSE)
cmap       <- read_one("cohort_sections")

s21_iso <- s21_iso[s21_iso$contrast == "interaction_diff", ]
if (!is.null(s21_sc)) s21_sc <- s21_sc[s21_sc$contrast == "interaction_diff", ]

## patient ID from the cohort table (BRCA IDs do not follow one pattern)
pid_of <- function(sid) {
  p <- cmap$patient_id[match(sid, cmap$sample_id)]
  if (anyNA(p)) stop("No patient for: ", paste(unique(sid[is.na(p)]), collapse = ", "))
  p
}
map <- unique(data.frame(sample_id = s21_iso$sample_id, patient = pid_of(s21_iso$sample_id)))
map <- map[order(map$patient, map$sample_id), ]
cat(sprintf("\nSections %d, patients %d\n", nrow(map), length(unique(map$patient))))
print(table(map$patient))

## ============================================================================
## 21.2  POOLING
## ============================================================================
dl_pool <- function(y, v) {
  ok <- is.finite(y) & is.finite(v) & v > 0; y <- y[ok]; v <- v[ok]; k <- length(y)
  if (k == 0) return(c(est = NA, lower = NA, upper = NA, p = NA, k = 0, I2 = NA))
  w  <- 1 / v; fe <- sum(w * y) / sum(w)
  Q  <- sum(w * (y - fe)^2); df <- k - 1
  C  <- sum(w) - sum(w^2) / sum(w)
  t2 <- if (k > 1 && C > 0) max(0, (Q - df) / C) else 0
  ws <- 1 / (v + t2); est <- sum(ws * y) / sum(ws); se <- sqrt(1 / sum(ws))
  I2 <- if (k > 1 && Q > 0) max(0, (Q - df) / Q) * 100 else 0
  c(est = est, lower = est - 1.96 * se, upper = est + 1.96 * se,
    p = 2 * pnorm(-abs(est / se)), k = k, I2 = I2)
}

fe_collapse_pid <- function(d, pid) {
  d$patient <- pid
  do.call(rbind, lapply(split(d, d$patient), function(x) {
    ok <- is.finite(x$est) & is.finite(x$var) & x$var > 0; x <- x[ok, ]
    if (!nrow(x)) return(NULL)
    w <- 1 / x$var
    data.frame(patient = x$patient[1], n_sections = nrow(x),
               est = sum(w * x$est) / sum(w), var = 1 / sum(w))
  }))
}
pat_pool_pid <- function(d, pid) { pc <- fe_collapse_pid(d, pid); dl_pool(pc$est, pc$var) }
run_pool <- function(d) if (identical(POOL_LEVEL, "section")) dl_pool(d$est, d$var) else pat_pool_pid(d, pid_of(d$sample_id))

## ============================================================================
## 21.3  CHECK: REPRODUCE THE PIPELINE'S POOLED RESULTS
## ============================================================================
cat(sprintf("\nCheck: does this pooling reproduce the pipeline's (POOL_LEVEL = %s)?\n", POOL_LEVEL))
chk <- list()
for (og in unique(s21_iso$outcome)) {
  mine <- run_pool(s21_iso[s21_iso$outcome == og, ])
  ref  <- s21_iso_pl$est[s21_iso_pl$outcome == og & s21_iso_pl$contrast == "interaction_diff"]
  chk[[length(chk) + 1]] <- data.frame(label = paste("S15", og), mine = mine[["est"]],
                                       ref = if (length(ref) == 1) ref else NA_real_)
}
if (!is.null(s21_sc)) for (og in unique(s21_sc$outcome)) for (md in unique(s21_sc$model)) {
  s <- s21_sc[s21_sc$outcome == og & s21_sc$model == md, ]
  mine <- run_pool(s)
  ref  <- s21_sc_pl$est[s21_sc_pl$outcome == og & s21_sc_pl$model == md & s21_sc_pl$contrast == "interaction_diff"]
  chk[[length(chk) + 1]] <- data.frame(label = paste("S20", og, md), mine = mine[["est"]],
                                       ref = if (length(ref) == 1) ref else NA_real_)
}
chk <- do.call(rbind, chk)
chk$ok <- is.finite(chk$mine) & is.finite(chk$ref) & abs(chk$mine - chk$ref) < TOL
for (i in seq_len(nrow(chk)))
  cat(sprintf("  %-26s mine %+.6f  pipeline %+.6f  %s\n",
              chk$label[i], chk$mine[i], chk$ref[i], if (isTRUE(chk$ok[i])) "ok" else "*** MISMATCH ***"))
if (!all(chk$ok)) stop("Section 21 halted: this pooling does not reproduce the pipeline's. Nothing below is interpretable.")
cat("  PASSED\n")

## ============================================================================
## 21.4  PATIENT LEVEL vs SECTION LEVEL
## ============================================================================
analyses <- list()
for (og in unique(s21_iso$outcome))
  analyses[[paste("S15", og)]] <- s21_iso[s21_iso$outcome == og, ]
if (!is.null(s21_sc)) for (og in unique(s21_sc$outcome)) for (md in c("M1_adjusted", "M2_low_stroma"))
  analyses[[paste("S20", og, md)]] <- s21_sc[s21_sc$outcome == og & s21_sc$model == md, ]
analyses <- analyses[vapply(analyses, nrow, integer(1)) > 0]

pid_merged <- function(sid) { p <- pid_of(sid); p[p == "HT339B2"] <- "HT339B1"; p }

has_metafor <- requireNamespace("metafor", quietly = TRUE)
res21 <- do.call(rbind, lapply(names(analyses), function(lb) {
  d   <- analyses[[lb]]
  s   <- dl_pool(d$est, d$var)
  pp  <- pat_pool_pid(d, pid_of(d$sample_id))
  pm  <- pat_pool_pid(d, pid_merged(d$sample_id))
  ml <- c(est = NA, lower = NA, upper = NA)
  if (has_metafor && nrow(d) >= 3) {
    fit <- tryCatch(metafor::rma.mv(yi = d$est, V = d$var,
                                    random = ~ 1 | patient / sample_id,
                                    data = data.frame(d, patient = pid_of(d$sample_id)),
                                    method = "REML"),
                    error = function(e) NULL)
    if (!is.null(fit)) ml <- c(est = as.numeric(fit$b), lower = fit$ci.lb, upper = fit$ci.ub)
  }
  data.frame(analysis = lb,
             k_pat = pp[["k"]], pat_est = pp[["est"]], pat_lo = pp[["lower"]], pat_hi = pp[["upper"]], pat_p = pp[["p"]], pat_I2 = pp[["I2"]],
             k_sec = s[["k"]],  sec_est = s[["est"]],  sec_lo = s[["lower"]],  sec_hi = s[["upper"]],  sec_p = s[["p"]],
             k_patm = pm[["k"]], patm_est = pm[["est"]], patm_lo = pm[["lower"]], patm_hi = pm[["upper"]],
             ml_est = ml[["est"]], ml_lo = ml[["lower"]], ml_hi = ml[["upper"]],
             row.names = NULL)
}))
res21$ci_excludes_0_pat <- res21$pat_hi < 0 | res21$pat_lo > 0
res21$ci_excludes_0_sec <- res21$sec_hi < 0 | res21$sec_lo > 0
write.csv(res21, file.path(OUT_DIR, "SECTION_level_sensitivity.csv"), row.names = FALSE)

cat("\n", strrep("-", 74), "\n", sep = "")
cat("SLOPE DIFFERENCE: PATIENTS (primary) vs SECTIONS\n")
cat(strrep("-", 74), "\n", sep = "")
fmt <- function(e, l, h) if (is.na(e)) "            n/a            " else sprintf("%+.4f (%+.4f to %+.4f)", e, l, h)
for (i in seq_len(nrow(res21))) with(res21[i, ], {
  cat(sprintf("\n  %s\n", analysis))
  cat(sprintf("    patients     k=%2d  %s  p = %.2g  I2 = %.0f%%\n", as.integer(k_pat), fmt(pat_est, pat_lo, pat_hi), pat_p, pat_I2))
  cat(sprintf("    sections     k=%2d  %s  p = %.2g%s\n", as.integer(k_sec), fmt(sec_est, sec_lo, sec_hi), sec_p,
              if (isTRUE(ci_excludes_0_sec != ci_excludes_0_pat)) "   <-- CONCLUSION CHANGES" else ""))
  cat(sprintf("    HT339 merged k=%2d  %s\n", as.integer(k_patm), fmt(patm_est, patm_lo, patm_hi)))
  if (has_metafor) cat(sprintf("    3-level            %s\n", fmt(ml_est, ml_lo, ml_hi)))
})

cat("\nNOTES\n")
if (!has_metafor) cat("  metafor not installed, so the three-level check was skipped.\n")
if (is.null(s21_sc)) cat("  Section 20 outputs not found, so only section 15 contrasts were checked.\n")
cat("  Within a patient, sections are combined as fixed effect, assuming they are\n")
cat("  independent; serial sections are not, so each patient's variance is too\n")
cat("  small. The three-level model does not make that assumption.\n")
cat("  Output: ", file.path(OUT_DIR, "SECTION_level_sensitivity.csv"), "\n", sep = "")
cat(strrep("=", 74), "\n", sep = "")

################################################################################
## 22. INFLUENCE OF SINGLE PATIENTS, AND THE METASTASIS
##   22.A  Leave one PATIENT out (with all of their sections).
##   22.B  With vs without the metastatic patient (HT268B1). Only possible on
##         the output of a run with INCLUDE_METASTASIS <- TRUE.
##   22.C  Sign concordance across patients and sections, before pooling.
##   No model is refitted. Standalone: needs only the aggregate CSVs.
################################################################################

if (!identical(POOL_LEVEL, "patient")) {
  cat("\nSection 22 skipped: it is written for POOL_LEVEL = 'patient'. Section 21 covers the section level.\n")
} else {
  
  cat("\n", strrep("=", 74), "\n", sep = "")
  cat("SECTION 22  SINGLE-PATIENT INFLUENCE AND THE METASTASIS\n")
  cat(strrep("=", 74), "\n", sep = "")
  cat("\nInputs:\n")
  s22_iso    <- read_one("ISOFORM_contrasts_per_sample")
  s22_iso_pl <- read_one("ISOFORM_slope_difference_pooled")
  s22_sc     <- read_one("SC_contrasts_per_sample", required = FALSE)
  s22_sc_pl  <- read_one("SC_slope_difference_pooled", required = FALSE)
  cmap       <- read_one("cohort_sections")
  
  s22_iso <- s22_iso[s22_iso$contrast == "interaction_diff", ]
  if (!is.null(s22_sc)) s22_sc <- s22_sc[s22_sc$contrast == "interaction_diff", ]
  
  met_patients <- unique(cmap$patient_id[cmap$primary_metastasis == "Metastasis"])
  
  fe_collapse22 <- function(d) fe_collapse_pid(d, pid_of(d$sample_id))
  pat_pool22    <- function(pc) dl_pool(pc$est, pc$var)
  excl0 <- function(lo, hi) (is.finite(lo) && is.finite(hi)) && (hi < 0 || lo > 0)
  
  analyses22 <- list()
  for (og in unique(s22_iso$outcome)) analyses22[[paste("S15", og)]] <- s22_iso[s22_iso$outcome == og, ]
  if (!is.null(s22_sc))
    for (og in unique(s22_sc$outcome)) for (md in c("M1_adjusted", "M2_low_stroma"))
      analyses22[[paste("S20", og, md)]] <- s22_sc[s22_sc$outcome == og & s22_sc$model == md, ]
  analyses22 <- analyses22[vapply(analyses22, nrow, integer(1)) > 0]
  
  ## ---- check: reproduce the pipeline -------------------------------------
  cat("\nCheck: reproduce the pipeline's pooled estimates\n")
  ok_all <- TRUE
  for (lb in names(analyses22)) {
    d <- analyses22[[lb]]
    mine <- pat_pool22(fe_collapse22(d))[["est"]]
    ref <- if (grepl("^S15", lb)) {
      og <- sub("^S15 ", "", lb)
      s22_iso_pl$est[s22_iso_pl$outcome == og & s22_iso_pl$contrast == "interaction_diff"]
    } else {
      p <- strsplit(lb, " ")[[1]]
      s22_sc_pl$est[s22_sc_pl$outcome == p[2] & s22_sc_pl$model == p[3] & s22_sc_pl$contrast == "interaction_diff"]
    }
    good <- length(ref) == 1 && isTRUE(abs(mine - ref) < TOL)
    ok_all <- ok_all && good
    cat(sprintf("  %-26s mine %+.6f  pipeline %+.6f  %s\n", lb, mine,
                if (length(ref) == 1) ref else NA_real_, if (good) "ok" else "*** MISMATCH ***"))
  }
  if (!ok_all) stop("Section 22 halted: pooling does not reproduce the pipeline. Nothing below is interpretable.")
  cat("  PASSED\n")
  
  ## ---- leave one out, metastasis, concordance ----------------------------
  res22 <- list(); loo_rows <- list()
  for (lb in names(analyses22)) {
    d    <- analyses22[[lb]]
    pc   <- fe_collapse22(d)
    pall <- pat_pool22(pc)
    
    loo <- do.call(rbind, lapply(pc$patient, function(s) {
      r <- pat_pool22(pc[pc$patient != s, ])
      data.frame(analysis = lb, dropped = s, n_sections_dropped = pc$n_sections[pc$patient == s],
                 est = r[["est"]], lower = r[["lower"]], upper = r[["upper"]], p = r[["p"]], k = r[["k"]],
                 shift = r[["est"]] - pall[["est"]],
                 still_excludes_0 = excl0(r[["lower"]], r[["upper"]]), row.names = NULL)
    }))
    loo_rows[[lb]] <- loo
    
    has_met <- any(pc$patient %in% met_patients)
    prim <- if (has_met) pat_pool22(pc[!pc$patient %in% met_patients, ])
    else c(est = NA, lower = NA, upper = NA, p = NA, k = NA, I2 = NA)
    
    res22[[lb]] <- data.frame(
      analysis = lb,
      k = pall[["k"]], est = pall[["est"]], lower = pall[["lower"]], upper = pall[["upper"]],
      loo_min = min(loo$est, na.rm = TRUE), loo_max = max(loo$est, na.rm = TRUE),
      loo_max_abs_shift = max(abs(loo$shift), na.rm = TRUE),
      most_influential = loo$dropped[which.max(abs(loo$shift))],
      loo_any_conclusion_change = any(loo$still_excludes_0 != excl0(pall[["lower"]], pall[["upper"]])),
      has_metastasis = has_met,
      primary_k = prim[["k"]], primary_est = prim[["est"]],
      primary_lower = prim[["lower"]], primary_upper = prim[["upper"]],
      same_sign_patients = sum(sign(pc$est) == sign(pall[["est"]]), na.rm = TRUE), n_patients = nrow(pc),
      same_sign_sections = sum(sign(d$est)  == sign(pall[["est"]]), na.rm = TRUE), n_sections = nrow(d),
      row.names = NULL)
  }
  res_df <- do.call(rbind, res22)
  loo_df <- do.call(rbind, loo_rows)
  write.csv(res_df, file.path(OUT_DIR, "LOO_summary.csv"), row.names = FALSE)
  write.csv(loo_df, file.path(OUT_DIR, "LOO_per_dropped_patient.csv"), row.names = FALSE)
  
  for (i in seq_len(nrow(res_df))) with(res_df[i, ], {
    cat("\n", strrep("-", 74), "\n", sep = "")
    cat(sprintf("  %s\n", analysis))
    cat(strrep("-", 74), "\n", sep = "")
    cat(sprintf("    all patients     k=%2d  %+.4f (%+.4f to %+.4f)\n", as.integer(k), est, lower, upper))
    cat(sprintf("    leave one out          range %+.4f to %+.4f   largest shift %+.4f when %s is dropped\n",
                loo_min, loo_max, loo_max_abs_shift, most_influential))
    cat(sprintf("    conclusion            %s\n",
                if (loo_any_conclusion_change) "*** CHANGES when at least one patient is dropped ***"
                else "unchanged by dropping any single patient"))
    if (has_metastasis)
      cat(sprintf("    without metastasis k=%2d  %+.4f (%+.4f to %+.4f)%s\n",
                  as.integer(primary_k), primary_est, primary_lower, primary_upper,
                  if (excl0(primary_lower, primary_upper) != excl0(lower, upper))
                    "   <-- CONCLUSION CHANGES" else ""))
    else
      cat("    metastasis            not in this run (primary-only); rerun with INCLUDE_METASTASIS <- TRUE\n")
    cat(sprintf("    same direction as the pooled estimate in %d of %d patients (%d of %d sections)\n",
                as.integer(same_sign_patients), as.integer(n_patients),
                as.integer(same_sign_sections), as.integer(n_sections)))
  })
  
  hl <- "S15 P4HA1"
  if (hl %in% names(loo_rows)) {
    l <- loo_rows[[hl]]
    l <- l[order(-abs(l$shift)), ]
    cat("\n", strrep("-", 74), "\n", sep = "")
    cat("  MOST INFLUENTIAL PATIENTS FOR THE HEADLINE (S15 P4HA1)\n")
    cat(strrep("-", 74), "\n", sep = "")
    for (j in seq_len(min(3, nrow(l)))) with(l[j, ],
                                             cat(sprintf("    drop %-12s (%d sections) -> %+.4f  (shift %+.4f)\n",
                                                         dropped, as.integer(n_sections_dropped), est, shift)))
  }
  
  cat("\nNOTES\n")
  cat("  No model was refitted; every number comes from the stored per-section fits.\n")
  cat("  Model p values remain inflated by spatial autocorrelation.\n")
  cat(strrep("=", 74), "\n", sep = "")
}

################################################################################
## 23. IS THE RESULT CARRIED BY ONE ER-STRESS GENE?
##
##   The ER stress arm rests on one hand-built 24-gene module (section 19
##   reports whether it is coherent in this cohort). If one or two members
##   produce the whole finding, "P4HA1 tracks ER stress" is really "P4HA1
##   tracks that gene". The section 15 model is refitted once with the full
##   module and once with each gene removed in turn; nothing else changes.
##   The full-module fit must reproduce section 15, section by section.
##   Also answers whether keeping ATF4 / CHAC1 (PERK-responsive) while
##   dropping HYOU1 / DDIT4 / ERO1A matters. Outcome: P4HA1 only.
##   Run after the main pipeline in the same FRESH session.
################################################################################

stopifnot(exists("MODULES"), exists("work"), exists("agg"), exists("SEED"),
          exists("ADJUST_DEPTH"), exists("pool_with_p"), exists("fmt_p"),
          exists("write_out"), exists("assert_run_state"), exists("log_msg"),
          exists("bind"), exists("min_genes_for"), exists("check_against_s15"))

## ============================================================================
## 23.0  SETTINGS
## ============================================================================
EL_OUTCOME      <- "P4HA1"
EL_ES           <- "ERstressStrict"
EL_OTHER        <- c("Hypoxia", "CollagenSynthesis", "iCAF", "ERAD")
EL_MIN_SPOTS    <- 30
EL_MIN_EVENTS   <- 70
EL_ADJUST_DEPTH <- ADJUST_DEPTH

EL_GENES    <- MODULES[[EL_ES]]
EL_VARIANTS <- c("FULL", EL_GENES)      # FULL first, then one per dropped gene
EL_FLAGGED  <- intersect(c("ATF4", "CHAC1"), EL_GENES)

EL_RENAME <- c(Hypoxia = "zHyp", CollagenSynthesis = "zCol", iCAF = "ziCAF", ERAD = "zERAD")
EL_II_COL <- "zHyp:zCol"
EL_II_ER  <- "zHyp:zER"
EL_RHS <- paste0("zHyp * zCol + zHyp * zER + ziCAF + zERAD",
                 if (EL_ADJUST_DEPTH) " + zDepth" else "")

assert_run_state("Section 23")

## v5 (C3): only the sections section 15 fitted for this outcome
EL_SECTIONS <- s15_fitted(EL_OUTCOME)

cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 23  LEAVE ONE ER-STRESS GENE OUT\n")
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  module %s has %d genes, so %d model fits per section, on %d sections\n",
            EL_ES, length(EL_GENES), length(EL_VARIANTS), length(EL_SECTIONS)))
cat(sprintf("  genes: %s\n", paste(EL_GENES, collapse = ", ")))

for (m in c(EL_OTHER, EL_ES))
  if (EL_OUTCOME %in% MODULES[[m]])
    stop("Section 23 halted: ", EL_OUTCOME, " is inside ", m)

## ============================================================================
## 23.1  ONE SECTION: SCORE THE FIXED MODULES ONCE, THEN EVERY ER VARIANT
## ============================================================================
el_one <- function(rds_path, sid) {
  obj <- readRDS(rds_path)
  obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
  DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else Assays(obj)[1]
  if (!EL_OUTCOME %in% rownames(obj)) { log_msg("  %s: outcome absent -> skipped", sid); return(NULL) }
  
  for (nm in EL_OTHER) {
    full <- MODULES[[nm]]; g <- intersect(full, rownames(obj))
    if (length(g) < min_genes_for(full)) {
      log_msg("  %s: %s has %d/%d genes -> skipped", sid, nm, length(g), length(full)); return(NULL)
    }
    set.seed(SEED)
    obj <- tryCatch(add_module_score_fallback(obj, features = list(g), name = nm,
                                              assay = DefaultAssay(obj), nbin = 20, ctrl = 50),
                    error = function(e) NULL)
    if (is.null(obj) || !paste0(nm, "1") %in% colnames(obj@meta.data)) return(NULL)
    names(obj@meta.data)[names(obj@meta.data) == paste0(nm, "1")] <- nm
  }
  
  for (v in EL_VARIANTS) {
    set_v <- if (identical(v, "FULL")) EL_GENES else setdiff(EL_GENES, v)
    g <- intersect(set_v, rownames(obj))
    ## the FULL variant uses the full module's threshold, exactly as section 15
    if (length(g) < min_genes_for(if (identical(v, "FULL")) EL_GENES else set_v)) {
      log_msg("  %s: ER variant %s has %d/%d genes -> variant skipped", sid, v, length(g), length(set_v))
      next
    }
    nm <- paste0("ERv_", v)
    set.seed(SEED)
    o <- tryCatch(add_module_score_fallback(obj, features = list(g), name = nm,
                                            assay = DefaultAssay(obj), nbin = 20, ctrl = 50),
                  error = function(e) NULL)
    if (is.null(o) || !paste0(nm, "1") %in% colnames(o@meta.data)) next
    obj <- o
    names(obj@meta.data)[names(obj@meta.data) == paste0(nm, "1")] <- nm
  }
  
  er_cols <- grep("^ERv_", colnames(obj@meta.data), value = TRUE)
  if (!length(er_cols)) return(NULL)
  
  vars <- c(EL_OUTCOME, EL_OTHER, er_cols)
  cnt_col <- grep("^nCount_", colnames(obj@meta.data), value = TRUE)
  if (EL_ADJUST_DEPTH) { if (!length(cnt_col)) return(NULL); vars <- c(vars, cnt_col[1]) }
  d <- FetchData(obj, vars = vars, layer = "data")
  
  zz <- function(x) { s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) return(NULL)
  as.numeric((x - mean(x, na.rm = TRUE)) / s) }
  for (nm in EL_OTHER) {
    v <- zz(d[[nm]]); if (is.null(v)) return(NULL); d[[ EL_RENAME[[nm]] ]] <- v
  }
  if (EL_ADJUST_DEPTH) { v <- zz(log1p(d[[cnt_col[1]]])); if (is.null(v)) return(NULL); d$zDepth <- v }
  d
}

## ============================================================================
## 23.2  FIT ONE VARIANT
## ============================================================================
el_fit <- function(d, sid, variant) {
  col <- paste0("ERv_", variant)
  if (!col %in% names(d)) return(NULL)
  zer <- (function(x) { s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) NULL
  else as.numeric((x - mean(x, na.rm = TRUE)) / s) })(d[[col]])
  if (is.null(zer)) return(NULL)
  
  dd <- d; dd$zER <- zer
  y <- as.integer(dd[[EL_OUTCOME]] > 0); dd$y <- y
  ev <- sum(y); nev <- length(y) - ev
  if (length(y) < EL_MIN_SPOTS || ev < EL_MIN_EVENTS || nev < EL_MIN_EVENTS) return(NULL)
  
  f <- tryCatch(glm(as.formula(paste("y ~", EL_RHS)), data = dd, family = binomial()),
                error = function(e) NULL)
  if (is.null(f) || !f$converged) return(NULL)
  
  b <- coef(f); V <- vcov(f)
  has_depth <- "zDepth" %in% names(b)
  if (!identical(has_depth, isTRUE(EL_ADJUST_DEPTH)))
    stop(sprintf("Section 23, %s / %s: ADJUST_DEPTH = %s but the fitted model %s zDepth.",
                 sid, variant, EL_ADJUST_DEPTH, if (has_depth) "DOES contain" else "does NOT contain"))
  
  need <- c("zCol", "zER", EL_II_COL, EL_II_ER)
  if (!all(need %in% names(b)) || any(!is.finite(b[need]))) return(NULL)
  
  L <- setNames(numeric(length(b)), names(b)); L[EL_II_COL] <- 1; L[EL_II_ER] <- -1
  data.frame(sample_id = sid, variant = variant,
             est = as.numeric(L %*% b), var = as.numeric(t(L) %*% V %*% L),
             beta_ER = unname(b["zER"]), beta_hypER = unname(b[EL_II_ER]),
             n_spots = nrow(dd), n_events = ev, row.names = NULL)
}

## ============================================================================
## 23.3  RUN
## ============================================================================
log_msg("Section 23 start. %d variants. ADJUST_DEPTH = %s. run_id = %s.",
        length(EL_VARIANTS), EL_ADJUST_DEPTH, RUN_ID)

el_acc <- list()
for (i in seq_len(nrow(work))) {
  path <- work$path[i]; sid <- work$sample_id[i]
  if (is.na(path) || !file.exists(path)) next
  if (!sid %in% EL_SECTIONS) next
  log_msg("[elo %d/%d] %s", i, nrow(work), sid)
  d <- tryCatch(el_one(path, sid),
                error = function(e) { log_msg("  elo skipped: %s", conditionMessage(e)); NULL })
  if (is.null(d)) next
  for (v in EL_VARIANTS) {
    r <- el_fit(d, sid, v)
    if (!is.null(r)) el_acc[[length(el_acc) + 1]] <- r
  }
  rm(d); gc(verbose = FALSE)
}
el_all <- bind(el_acc)
if (nrow(el_all) == 0) stop("Section 23: nothing was fitted.")
write_out(el_all, "ERLOO_per_sample")

## ============================================================================
## 23.4  SAFETY CHECK: THE FULL MODULE MUST REPRODUCE SECTION 15
## ============================================================================
check_against_s15(el_all[el_all$variant == "FULL", ], EL_OUTCOME, "Section 23")

## every variant is pooled on the sections where the FULL module was fitted,
## so a removal is never confounded with a change of sections
el_secs <- unique(el_all$sample_id[el_all$variant == "FULL"])
el_use  <- el_all[el_all$sample_id %in% el_secs, ]
el_pooled <- el_use %>%
  group_by(variant) %>% group_modify(~ pool_with_p(.x$est, .x$var, .x$sample_id)) %>% ungroup()
write_out(el_pooled, "ERLOO_pooled")

## ============================================================================
## 23.5  REPORT
## ============================================================================
excl0 <- function(lo, hi) (is.finite(lo) && is.finite(hi)) && (hi < 0 || lo > 0)
full <- el_pooled[el_pooled$variant == "FULL", ]
base_sig <- excl0(full$lower, full$upper)

el_rep <- el_pooled[el_pooled$variant != "FULL", ]
el_rep$shift          <- el_rep$est - full$est
el_rep$pct_retained   <- 100 * el_rep$est / full$est
el_rep$still_excl0    <- vapply(seq_len(nrow(el_rep)),
                                function(i) excl0(el_rep$lower[i], el_rep$upper[i]), logical(1))
el_rep <- el_rep[order(-abs(el_rep$shift)), ]
write_out(el_rep, "ERLOO_ranked")

cat("\n", strrep("=", 74), "\n", sep = "")
cat(sprintf("SECTION 23  [%s]  run_id %s\n",
            if (EL_ADJUST_DEPTH) "DEPTH ADJUSTED" else "UNADJUSTED", RUN_ID))
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("\n  full module (%d genes) : %+.4f  (95%% CI %+.4f to %+.4f)  %s  p = %s\n",
            length(EL_GENES), full$est, full$lower, full$upper, k_lab(full), fmt_p(full$p)))

cat("\n  each gene removed in turn, ordered by how much it moves the result\n\n")
cat(sprintf("  %-10s %9s %21s %5s %9s %8s\n", "dropped", "estimate", "95% CI", "k", "shift", "kept %"))
for (i in seq_len(nrow(el_rep))) with(el_rep[i, ],
                                      cat(sprintf("  %-10s %+9.4f  %+8.4f to %+8.4f %5d %+9.4f %7.0f%%%s%s\n",
                                                  variant, est, lower, upper, as.integer(k), shift, pct_retained,
                                                  if (still_excl0 != base_sig) "   <-- CONCLUSION CHANGES" else "",
                                                  if (variant %in% EL_FLAGGED) "   [flagged gene]" else "")))

cat(sprintf("\n  range across all %d single-gene removals: %+.4f to %+.4f\n",
            nrow(el_rep), min(el_rep$est), max(el_rep$est)))

changed <- el_rep$variant[el_rep$still_excl0 != base_sig]
cat("\n", strrep("-", 74), "\n", sep = "")
cat("  VERDICT\n")
cat(strrep("-", 74), "\n", sep = "")
if (!base_sig) {
  cat("  The full-module slope difference does not exclude zero in this cohort,\n")
  cat("  so there is no result for a single gene to carry. Removals are shown\n")
  cat("  for completeness only.\n")
} else if (!length(changed)) {
  cat("  NO SINGLE GENE CARRIES THE RESULT. Removing any one gene leaves the\n")
  cat("  slope difference in the same direction and still excluding zero.\n")
} else {
  cat(sprintf("  *** THE RESULT DEPENDS ON %d GENE(S): %s ***\n",
              length(changed), paste(changed, collapse = ", ")))
  cat("  Removing it changes whether the interval excludes zero. The claim has\n")
  cat("  to be rewritten around that gene rather than around 'ER stress'.\n")
}

if (length(EL_FLAGGED)) {
  cat("\n  THE ATF4 / CHAC1 QUESTION\n")
  for (g in EL_FLAGGED) {
    r <- el_rep[el_rep$variant == g, ]
    if (nrow(r) == 1)
      cat(sprintf("    without %-6s %+.4f  (%.0f%% of the full estimate)\n", g, r$est, r$pct_retained))
  }
  cat("    If these are close to the full estimate, the inconsistent exclusion\n")
  cat("    rule does not affect the result.\n")
}

cat("\nCAVEATS\n")
cat("  Only P4HA1 was refitted. Removing one of 24 genes is a small perturbation;\n")
cat("  section 24 tests whether different ER stress DEFINITIONS agree.\n")
cat("  Model p values remain inflated by spatial autocorrelation.\n")
cat(strrep("=", 74), "\n", sep = "")

log_msg("Section 23 complete.")

################################################################################
## 24. DOES THE ER-STRESS RESULT HOLD WITH INDEPENDENT UPR GENE SETS?
##
##   Section 23 tested dependence on MEMBERS; this tests dependence on the
##   DEFINITION. The section 15 model is refitted with zER built from each
##   published set in turn; ERstressStrict is refitted as the reference and
##   must reproduce section 15.
##   Diagnostics printed before any estimate: overlap with ERstressStrict
##   (Jaccard) and correlation with the hypoxia module. Genes shared with
##   Hypoxia, CollagenSynthesis or iCAF are stripped automatically.
##   INPUT: a folder of MSigDB TSV files (data/upr_sets), e.g.
##     HALLMARK_UNFOLDED_PROTEIN_RESPONSE (the parent of the hand-built
##       module, so NOT independent), GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_
##       STRESS, and the IRE1 / PERK / ATF6 arms of the UPR.
################################################################################

stopifnot(exists("MODULES"), exists("work"), exists("agg"), exists("SEED"),
          exists("ADJUST_DEPTH"), exists("pool_with_p"), exists("fmt_p"),
          exists("write_out"), exists("assert_run_state"), exists("log_msg"),
          exists("bind"), exists("min_genes_for"), exists("check_against_s15"))

## ============================================================================
## 24.0  SETTINGS
## ============================================================================
UPR_DIR          <- file.path(PROJECT_ROOT, "data", "upr_sets")
UPR_OUTCOME      <- "P4HA1"
UPR_REF          <- "ERstressStrict"
UPR_OTHER        <- c("Hypoxia", "CollagenSynthesis", "iCAF", "ERAD")
UPR_MIN_GENES    <- 8
UPR_MIN_SPOTS    <- 30
UPR_MIN_EVENTS   <- 70
UPR_ADJUST_DEPTH <- ADJUST_DEPTH
UPR_COLLIN_WARN  <- 0.50
UPR_JACCARD_WARN <- 0.50
UPR_FAMILY_MAX   <- 0.05

UPR_FAMILY <- c("P4HA1", "P4HA2", "P4HA3", "P4HB", "PLOD1", "PLOD2", "PLOD3", "PPIB")
UPR_RENAME <- c(Hypoxia = "zHyp", CollagenSynthesis = "zCol", iCAF = "ziCAF", ERAD = "zERAD")
UPR_II_COL <- "zHyp:zCol"
UPR_II_ER  <- "zHyp:zER"
UPR_RHS <- paste0("zHyp * zCol + zHyp * zER + ziCAF + zERAD",
                  if (UPR_ADJUST_DEPTH) " + zDepth" else "")

assert_run_state("Section 24")
UPR_SECTIONS <- s15_fitted(UPR_OUTCOME)

## ============================================================================
## 24.1  READ THE CANDIDATE SETS FROM THE TSVs
## ============================================================================
if (!dir.exists(UPR_DIR)) stop("Section 24 needs the folder ", UPR_DIR, " with MSigDB TSVs in it.")
tsvs <- list.files(UPR_DIR, pattern = "\\.tsv$", full.names = TRUE, ignore.case = TRUE)
if (!length(tsvs)) stop("No .tsv files in ", UPR_DIR)

read_msigdb <- function(path) {
  x <- readLines(path, warn = FALSE)
  key <- grep("^(MAPPED_SYMBOLS|GENE_SYMBOLS)\t", x, value = TRUE)
  if (!length(key)) stop("No GENE_SYMBOLS or MAPPED_SYMBOLS row in ", basename(path))
  g <- strsplit(sub("^[A-Z_]+\t", "", key[1]), "[,\t ]+")[[1]]
  unique(g[nzchar(g)])
}
short_name <- function(path) {
  nm <- sub("\\.tsv$", "", basename(path), ignore.case = TRUE)
  nm <- sub("\\.v[0-9].*$", "", nm)
  nm <- sub("^HALLMARK_", "HALLMARK.", nm)
  nm <- sub("^GOBP_", "GO.", nm)
  nm <- sub("^REACTOME_", "REACT.", nm)
  substr(nm, 1, 34)
}

cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 24.1  CANDIDATE SETS\n")
cat(strrep("=", 74), "\n", sep = "")

raw_sets <- list()
for (f in tsvs) {
  g <- tryCatch(read_msigdb(f), error = function(e) { cat("  skipped ", basename(f), ": ",
                                                          conditionMessage(e), "\n", sep = ""); NULL })
  if (!is.null(g)) raw_sets[[short_name(f)]] <- g
}
if (!length(raw_sets)) stop("Section 24: no usable gene set was read.")

## ============================================================================
## 24.2  AUDIT AND CLEAN
##       Family genes (e.g. P4HB) are stripped and named; a set made of more
##       than UPR_FAMILY_MAX family genes is dropped. Genes shared with
##       Hypoxia, collagen or iCAF are stripped.
## ============================================================================
STRIP_AGAINST <- c("Hypoxia", "CollagenSynthesis", "iCAF")
strip_pool <- unique(unlist(MODULES[STRIP_AGAINST], use.names = FALSE))
ref_genes  <- MODULES[[UPR_REF]]

jacc <- function(a, b) { u <- length(union(a, b)); if (u) length(intersect(a, b)) / u else NA_real_ }

UPR_SETS <- list(); audit_rows <- list()
for (nm in names(raw_sets)) {
  g0  <- raw_sets[[nm]]
  fam <- intersect(g0, c(UPR_FAMILY, UPR_OUTCOME))
  stripped <- intersect(g0, strip_pool)
  g <- setdiff(g0, c(strip_pool, fam))
  fam_frac <- if (length(g0)) length(fam) / length(g0) else 0
  keep <- length(g) >= UPR_MIN_GENES && fam_frac <= UPR_FAMILY_MAX
  
  cat(sprintf("\n  %-34s %3d genes\n", nm, length(g0)))
  if (length(fam))
    cat(sprintf("     *** CIRCULAR GENES REMOVED: %s  (%.1f%% of the set) ***\n",
                paste(sort(fam), collapse = ", "), 100 * fam_frac))
  if (length(stripped))
    cat(sprintf("     stripped (shared with %s): %s\n",
                paste(STRIP_AGAINST, collapse = "/"), paste(sort(stripped), collapse = ", ")))
  cat(sprintf("     kept %3d genes | overlap with %s: %d shared, Jaccard %.2f%s\n",
              length(g), UPR_REF, length(intersect(g, ref_genes)), jacc(g, ref_genes),
              if (!keep) sprintf("   *** DROPPED: %s ***",
                                 if (fam_frac > UPR_FAMILY_MAX) "built around family genes"
                                 else "too small") else ""))
  audit_rows[[nm]] <- data.frame(set = nm, n_raw = length(g0), n_kept = length(g),
                                 family_removed = paste(sort(fam), collapse = ";"),
                                 family_fraction = fam_frac,
                                 stripped = paste(sort(stripped), collapse = ";"),
                                 n_shared_with_ref = length(intersect(g, ref_genes)),
                                 jaccard_with_ref = jacc(g, ref_genes),
                                 used = keep, row.names = NULL)
  if (keep) UPR_SETS[[nm]] <- g
}
audit <- bind(audit_rows)
write_out(audit, "UPR_set_audit")
if (!length(UPR_SETS)) stop("Section 24: every candidate set was dropped.")
write_out(data.frame(set = rep(names(UPR_SETS), lengths(UPR_SETS)),
                     gene = unlist(UPR_SETS, use.names = FALSE)), "UPR_genes_used")

UPR_VARIANTS <- c(UPR_REF, names(UPR_SETS))
for (nm in names(UPR_SETS)) MODULES[[nm]] <- UPR_SETS[[nm]]   # left in place for section 25

## ============================================================================
## 24.3  ONE SECTION: SCORE EVERYTHING ONCE
## ============================================================================
upr_one <- function(rds_path, sid) {
  obj <- readRDS(rds_path)
  obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
  DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else Assays(obj)[1]
  if (!UPR_OUTCOME %in% rownames(obj)) { log_msg("  %s: outcome absent -> skipped", sid); return(NULL) }
  
  score <- function(obj, nm, gset) {
    g <- intersect(gset, rownames(obj))
    if (length(g) < min_genes_for(gset)) return(NULL)
    set.seed(SEED)
    o <- tryCatch(add_module_score_fallback(obj, features = list(g), name = nm,
                                            assay = DefaultAssay(obj), nbin = 20, ctrl = 50),
                  error = function(e) NULL)
    if (is.null(o) || !paste0(nm, "1") %in% colnames(o@meta.data)) return(NULL)
    names(o@meta.data)[names(o@meta.data) == paste0(nm, "1")] <- nm
    o
  }
  
  for (nm in UPR_OTHER) {
    o <- score(obj, nm, MODULES[[nm]])
    if (is.null(o)) { log_msg("  %s: %s could not be scored -> skipped", sid, nm); return(NULL) }
    obj <- o
  }
  have <- character(0)
  for (v in UPR_VARIANTS) {
    o <- score(obj, paste0("UPRv_", v), if (identical(v, UPR_REF)) MODULES[[UPR_REF]] else UPR_SETS[[v]])
    if (is.null(o)) { log_msg("  %s: variant %s too few genes here -> variant skipped", sid, v); next }
    obj <- o; have <- c(have, v)
  }
  if (!length(have)) return(NULL)
  
  vars <- c(UPR_OUTCOME, UPR_OTHER, paste0("UPRv_", have))
  cnt_col <- grep("^nCount_", colnames(obj@meta.data), value = TRUE)
  if (UPR_ADJUST_DEPTH) { if (!length(cnt_col)) return(NULL); vars <- c(vars, cnt_col[1]) }
  d <- FetchData(obj, vars = vars, layer = "data")
  
  zz <- function(x) { s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) return(NULL)
  as.numeric((x - mean(x, na.rm = TRUE)) / s) }
  for (nm in UPR_OTHER) { v <- zz(d[[nm]]); if (is.null(v)) return(NULL); d[[ UPR_RENAME[[nm]] ]] <- v }
  if (UPR_ADJUST_DEPTH) { v <- zz(log1p(d[[cnt_col[1]]])); if (is.null(v)) return(NULL); d$zDepth <- v }
  
  sp <- function(a, b) suppressWarnings(cor(a, b, method = "spearman"))
  ref_col <- paste0("UPRv_", UPR_REF)
  attr(d, "diag") <- do.call(rbind, lapply(have, function(v) data.frame(
    sample_id = sid, variant = v,
    rho_vs_hypoxia = sp(d[[paste0("UPRv_", v)]], d$Hypoxia),
    rho_vs_ref     = if (ref_col %in% names(d)) sp(d[[paste0("UPRv_", v)]], d[[ref_col]]) else NA_real_,
    row.names = NULL)))
  attr(d, "have") <- have
  d
}

## ============================================================================
## 24.4  FIT ONE VARIANT
## ============================================================================
upr_fit <- function(d, sid, variant) {
  col <- paste0("UPRv_", variant)
  if (!col %in% names(d)) return(NULL)
  zer <- (function(x) { s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) NULL
  else as.numeric((x - mean(x, na.rm = TRUE)) / s) })(d[[col]])
  if (is.null(zer)) return(NULL)
  
  dd <- d; dd$zER <- zer
  y <- as.integer(dd[[UPR_OUTCOME]] > 0); dd$y <- y
  ev <- sum(y); nev <- length(y) - ev
  if (length(y) < UPR_MIN_SPOTS || ev < UPR_MIN_EVENTS || nev < UPR_MIN_EVENTS) return(NULL)
  
  f <- tryCatch(glm(as.formula(paste("y ~", UPR_RHS)), data = dd, family = binomial()),
                error = function(e) NULL)
  if (is.null(f) || !f$converged) return(NULL)
  
  b <- coef(f); V <- vcov(f)
  if (!identical("zDepth" %in% names(b), isTRUE(UPR_ADJUST_DEPTH)))
    stop(sprintf("Section 24, %s / %s: flag and fitted model disagree about zDepth.", sid, variant))
  need <- c("zCol", "zER", UPR_II_COL, UPR_II_ER)
  if (!all(need %in% names(b)) || any(!is.finite(b[need]))) return(NULL)
  
  L <- setNames(numeric(length(b)), names(b)); L[UPR_II_COL] <- 1; L[UPR_II_ER] <- -1
  data.frame(sample_id = sid, variant = variant,
             est = as.numeric(L %*% b), var = as.numeric(t(L) %*% V %*% L),
             beta_ER = unname(b["zER"]), beta_hypER = unname(b[UPR_II_ER]),
             n_spots = nrow(dd), n_events = ev, row.names = NULL)
}

## ============================================================================
## 24.5  RUN
## ============================================================================
log_msg("Section 24 start. %d variants. ADJUST_DEPTH = %s. run_id = %s.",
        length(UPR_VARIANTS), UPR_ADJUST_DEPTH, RUN_ID)

upr_acc <- list(fit = list(), diag = list())
for (i in seq_len(nrow(work))) {
  path <- work$path[i]; sid <- work$sample_id[i]
  if (is.na(path) || !file.exists(path)) next
  if (!sid %in% UPR_SECTIONS) next
  log_msg("[upr %d/%d] %s", i, nrow(work), sid)
  d <- tryCatch(upr_one(path, sid),
                error = function(e) { log_msg("  upr skipped: %s", conditionMessage(e)); NULL })
  if (is.null(d)) next
  upr_acc$diag[[length(upr_acc$diag) + 1]] <- attr(d, "diag")
  for (v in attr(d, "have")) {
    r <- upr_fit(d, sid, v)
    if (!is.null(r)) upr_acc$fit[[length(upr_acc$fit) + 1]] <- r
  }
  rm(d); gc(verbose = FALSE)
}
upr_all  <- bind(upr_acc$fit)
upr_diag <- bind(upr_acc$diag)
if (nrow(upr_all) == 0) stop("Section 24: nothing was fitted.")
write_out(upr_all,  "UPR_per_sample")
write_out(upr_diag, "UPR_diagnostics_per_sample")

## ============================================================================
## 24.6  SAFETY CHECK (reference vs section 15, section by section)
## ============================================================================
check_against_s15(upr_all[upr_all$variant == UPR_REF, ], UPR_OUTCOME, "Section 24")

## every set pooled on the sections where the reference was fitted
upr_secs <- unique(upr_all$sample_id[upr_all$variant == UPR_REF])
upr_pooled <- upr_all[upr_all$sample_id %in% upr_secs, ] %>%
  group_by(variant) %>% group_modify(~ pool_with_p(.x$est, .x$var, .x$sample_id)) %>% ungroup()
write_out(upr_pooled, "UPR_pooled")
refr <- upr_pooled[upr_pooled$variant == UPR_REF, ]

## ============================================================================
## 24.7  REPORT
## ============================================================================
excl0 <- function(lo, hi) (is.finite(lo) && is.finite(hi)) && (hi < 0 || lo > 0)
ref_sig <- excl0(refr$lower, refr$upper)

cat("\n", strrep("=", 74), "\n", sep = "")
cat(sprintf("SECTION 24  [%s]  run_id %s\n",
            if (UPR_ADJUST_DEPTH) "DEPTH ADJUSTED" else "UNADJUSTED", RUN_ID))
cat(strrep("=", 74), "\n", sep = "")

cat(sprintf("\n--- diagnostics, read these before the estimates (medians over %ss) ---\n", UNIT_LABEL))
cat(sprintf("  %-34s %10s %12s %12s\n", "set", "Jaccard", "rho hypoxia", "rho with ref"))
for (v in UPR_VARIANTS) {
  dg <- upr_diag[upr_diag$variant == v, ]
  jc <- if (identical(v, UPR_REF)) 1 else audit$jaccard_with_ref[audit$set == v]
  if (!nrow(dg)) next
  mh <- unit_median(dg$rho_vs_hypoxia, dg$sample_id); mr <- unit_median(dg$rho_vs_ref, dg$sample_id)
  ## v5: standard flag order ("%+12.2f")
  cat(sprintf("  %-34s %10.2f %+12.2f %+12.2f%s%s\n", v, jc, mh, mr,
              if (is.finite(mh) && abs(mh) >= UPR_COLLIN_WARN) "  <-- tracks hypoxia" else "",
              if (!identical(v, UPR_REF) && is.finite(jc) && jc >= UPR_JACCARD_WARN)
                "  <-- not independent" else ""))
}
cat("\n  Jaccard = overlap with the reference module AFTER stripping. A set that\n")
cat("  tracks hypoxia can manufacture the interaction being tested.\n")

cat("\n--- THE HEADLINE: slope difference under each definition of ER stress ---\n")
cat(sprintf("  %-34s %9s %21s %5s %11s\n", "set", "estimate", "95% CI", "k", "p"))
for (v in UPR_VARIANTS) {
  r <- upr_pooled[upr_pooled$variant == v, ]
  if (!nrow(r)) next
  cat(sprintf("  %-34s %+9.4f  %+8.4f to %+8.4f %5d %10s%s%s\n",
              v, r$est, r$lower, r$upper, as.integer(r$k), fmt_p(r$p),
              if (identical(v, UPR_REF)) "   [reference]" else "",
              if (!identical(v, UPR_REF) && excl0(r$lower, r$upper) == ref_sig &&
                  sign(r$est) == sign(refr$est)) "   agrees" else ""))
}

UPR_VERDICT <- NA_character_
cand <- upr_pooled[upr_pooled$variant != UPR_REF, ]
if (nrow(cand)) {
  indep <- audit$set[audit$used & is.finite(audit$jaccard_with_ref) &
                       audit$jaccard_with_ref < UPR_JACCARD_WARN]
  ci <- cand[cand$variant %in% indep, ]
  agree <- ci$variant[sign(ci$est) == sign(refr$est) &
                        vapply(seq_len(nrow(ci)), function(i) excl0(ci$lower[i], ci$upper[i]), logical(1))]
  same_dir <- ci$variant[sign(ci$est) == sign(refr$est)]
  
  cat("\n", strrep("-", 74), "\n", sep = "")
  cat("  VERDICT\n")
  cat(strrep("-", 74), "\n", sep = "")
  cat(sprintf("  independent sets tested (Jaccard < %.2f): %d\n", UPR_JACCARD_WARN, nrow(ci)))
  cat(sprintf("  same direction as the reference          : %d\n", length(same_dir)))
  cat(sprintf("  same direction AND excluding zero        : %d\n", length(agree)))
  if (!ref_sig) {
    UPR_VERDICT <- "reference not significant"
    cat("\n  The reference module itself does not exclude zero in this cohort, so\n")
    cat("  there is no result for other definitions to reproduce.\n")
  } else if (length(agree) >= 2) {
    UPR_VERDICT <- "replicated"
    cat("\n  THE DEFINITION IS NOT THE RESULT. Several independently published\n")
    cat("  definitions of ER stress give the same answer.\n")
  } else if (length(same_dir) >= 2) {
    UPR_VERDICT <- "same direction only"
    cat("\n  PARTIAL. The independent sets point the same way but not all exclude\n")
    cat("  zero: consistent in direction, weaker in precision.\n")
  } else {
    UPR_VERDICT <- "not replicated"
    cat("\n  *** THE RESULT IS SPECIFIC TO THE HAND-BUILT MODULE. Independent\n")
    cat("  definitions of ER stress do not reproduce it. ***\n")
  }
  cat("\n  If the IRE1 / PERK / ATF6 arms disagree, that is a biological result:\n")
  cat("  name which arm carries it.\n")
}

cat("\nCAVEATS\n")
cat("  Only P4HA1 was refitted. Sets were stripped of hydroxylase-family genes\n")
cat("  and of genes shared with Hypoxia, collagen and iCAF (see UPR_set_audit),\n")
cat("  so they are not the published sets as-is. Section 25 tests their coherence.\n")
cat(strrep("=", 74), "\n", sep = "")

log_msg("Section 24 complete.")

################################################################################
## 25. ARE THE PUBLISHED ER-STRESS SETS ANY MORE COHERENT THAN OURS?
##
##   Same method as section 19, applied to the reference modules, our
##   ERstressStrict and every UPR set section 24 read. If no ER-stress
##   definition is coherent, incoherence is a property of ER-stress gene sets
##   in Visium data, not of the hand-built module. If published sets cohere
##   and ours does not, ours is the weak one.
##   The null pool differs from section 19 (different modules are excluded),
##   so percentiles are not directly comparable; hypoxia and collagen are
##   included as positive controls and must still land high.
##   Needs section 24 in the same session (UPR_SETS).
################################################################################

stopifnot(exists("MODULES"), exists("work"), exists("agg"),
          exists("write_out"), exists("assert_run_state"), exists("log_msg"), exists("bind"))
if (!exists("UPR_SETS"))
  stop("Section 25 needs section 24 to have run in this session (UPR_SETS not found).")
if (!exists("coh_meanrho") || !exists("coh_pc1") || !exists("COH_COLS"))
  stop("Section 25 needs section 19 to have run in this session (coherence functions not found).")

## ============================================================================
## 25.0  SETTINGS
## ============================================================================
UC_REFS      <- c("Hypoxia", "CollagenSynthesis", "ERstressStrict", "ERAD")
UC_MODULES   <- c(UC_REFS, names(UPR_SETS))
UC_N_NULL    <- 30
UC_MIN_GENES <- 4
UC_MIN_SPOTS <- 30
UC_SEED      <- 90125
UC_CALIPER   <- 0.35
UC_PCT_MIN   <- 0.90

RATES_CSV <- file.path(agg, "DETECT_cohort_gene_rates.csv")
if (!file.exists(RATES_CSV)) stop("Section 25 needs DETECT_cohort_gene_rates.csv from section 14.")

assert_run_state("Section 25")

miss <- setdiff(UC_MODULES, names(MODULES))
if (length(miss)) stop("Not in MODULES: ", paste(miss, collapse = ", "))

## ============================================================================
## 25.1  NULL POOL, MATCHED ON DETECTION RATE
## ============================================================================
rates <- read.csv(RATES_CSV, stringsAsFactors = FALSE)
if (exists("RND_SECTIONS")) rates <- rates[rates$sample_id %in% RND_SECTIONS, , drop = FALSE]
gene_med <- rates %>%
  group_by(gene) %>%
  summarise(median_rate = median(detect_rate), n_tum = dplyr::n(), .groups = "drop") %>%
  filter(n_tum == dplyr::n_distinct(rates$sample_id))

UC_ALL_GENES <- unique(unlist(MODULES[UC_MODULES], use.names = FALSE))
uc_pool <- gene_med[!gene_med$gene %in% UC_ALL_GENES, ]

cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 25.1  SETUP\n")
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  genes present in all sections : %d\n", nrow(gene_med)))
cat(sprintf("  null pool after exclusions    : %d\n", nrow(uc_pool)))
for (nm in UC_MODULES) {
  det <- intersect(MODULES[[nm]], gene_med$gene)
  cat(sprintf("    %-34s %3d genes, %3d detected cohort-wide\n",
              nm, length(MODULES[[nm]]), length(det)))
}

uc_draw <- function(genes, n_sets) {
  det <- intersect(genes, gene_med$gene)
  if (length(det) < UC_MIN_GENES) return(NULL)
  target <- median(gene_med$median_rate[gene_med$gene %in% det])
  elig <- uc_pool$gene[uc_pool$median_rate >= target * (1 - UC_CALIPER) &
                         uc_pool$median_rate <= target * (1 + UC_CALIPER)]
  if (length(elig) < length(det) * 2) return(NULL)
  lapply(seq_len(n_sets), function(i) sample(elig, length(det)))
}
set.seed(UC_SEED)
UC_NULLS <- lapply(UC_MODULES, function(nm) uc_draw(MODULES[[nm]], UC_N_NULL))
names(UC_NULLS) <- UC_MODULES
for (nm in UC_MODULES)
  if (is.null(UC_NULLS[[nm]]))
    cat(sprintf("  *** %s: could not match null sets, it will be skipped ***\n", nm))

## ============================================================================
## 25.2  RUN  (the two measures are the section 19 functions)
## ============================================================================
uc_one <- function(rds_path, sid) {
  obj <- readRDS(rds_path)
  obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
  DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else Assays(obj)[1]
  if (ncol(obj) < UC_MIN_SPOTS) return(NULL)
  E <- tryCatch(GetAssayData(obj, assay = DefaultAssay(obj), layer = "data"),
                error = function(e) NULL)
  if (is.null(E) || !nrow(E)) return(NULL)
  rn <- rownames(E)
  grab <- function(g) { g <- intersect(g, rn); if (length(g) < UC_MIN_GENES) NULL else E[g, , drop = FALSE] }
  
  rows <- list()
  for (nm in UC_MODULES) {
    ns <- UC_NULLS[[nm]]; if (is.null(ns)) next
    real <- grab(MODULES[[nm]]); if (is.null(real)) next
    r_rho <- coh_meanrho(real); r_pc1 <- coh_pc1(real)
    n_rho <- n_pc1 <- rep(NA_real_, length(ns))
    for (j in seq_along(ns)) {
      mm <- grab(ns[[j]]); if (is.null(mm)) next
      n_rho[j] <- coh_meanrho(mm); n_pc1[j] <- coh_pc1(mm)
    }
    pr <- if (sum(is.finite(n_rho)) >= 10) mean(n_rho[is.finite(n_rho)] < r_rho) else NA_real_
    pp <- if (sum(is.finite(n_pc1)) >= 10) mean(n_pc1[is.finite(n_pc1)] < r_pc1) else NA_real_
    rows[[length(rows) + 1]] <- data.frame(
      sample_id = sid, module = nm, n_genes_used = nrow(real),
      mean_rho = r_rho, pc1_var = r_pc1,
      null_mean_rho = median(n_rho, na.rm = TRUE), null_pc1 = median(n_pc1, na.rm = TRUE),
      pct_rho = pr, pct_pc1 = pp, row.names = NULL)
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

log_msg("Section 25 start. %d modules, %d null sets each.", length(UC_MODULES), UC_N_NULL)
uc_acc <- list()
for (i in seq_len(nrow(work))) {
  path <- work$path[i]; sid <- work$sample_id[i]
  if (is.na(path) || !file.exists(path)) next
  log_msg("[uc %d/%d] %s", i, nrow(work), sid)
  r <- tryCatch(uc_one(path, sid),
                error = function(e) { log_msg("  uc skipped: %s", conditionMessage(e)); NULL })
  if (!is.null(r) && nrow(r)) uc_acc[[length(uc_acc) + 1]] <- r
}
uc_all <- bind(uc_acc)
if (nrow(uc_all) == 0) stop("Section 25: nothing computed.")
write_out(uc_all, "UPRCOH_per_sample")

## one value per unit, then the median over units (v5, C4)
uc_unit <- uc_all %>%
  mutate(unit = unit_of(sample_id)) %>%
  group_by(module, unit) %>%
  summarise(across(all_of(COH_COLS), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
uc_sum <- uc_unit %>%
  group_by(module) %>%
  summarise(k = dplyr::n(), median_genes = median(n_genes_used),
            median_rho = median(mean_rho, na.rm = TRUE),
            median_null_rho = median(null_mean_rho, na.rm = TRUE),
            median_pct_rho = median(pct_rho, na.rm = TRUE),
            median_pc1 = median(pc1_var, na.rm = TRUE),
            median_null_pc1 = median(null_pc1, na.rm = TRUE),
            median_pct_pc1 = median(pct_pc1, na.rm = TRUE), .groups = "drop")
uc_sum <- uc_sum[match(UC_MODULES, uc_sum$module), ]
uc_sum <- uc_sum[!is.na(uc_sum$module), ]
write_out(uc_sum, "UPRCOH_summary")

## ============================================================================
## 25.3  REPORT
## ============================================================================
cat("\n", strrep("=", 74), "\n", sep = "")
cat(sprintf("SECTION 25  COHERENCE OF THE PUBLISHED ER-STRESS SETS  run_id %s\n", RUN_ID))
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("\n  %-34s %4s %6s %8s %9s %9s %8s %9s\n",
            "module", "k", "genes", "rho", "null rho", "pctile", "PC1", "pctile"))
for (i in seq_len(nrow(uc_sum))) with(uc_sum[i, ],
                                      cat(sprintf("  %-34s %4d %6.0f %8.3f %9.3f %8.0f%% %8.3f %8.0f%%\n",
                                                  module, as.integer(k), median_genes, median_rho, median_null_rho, 100 * median_pct_rho,
                                                  median_pc1, 100 * median_pct_pc1)))

uc_passes <- function(mod) {
  r <- uc_sum[uc_sum$module == mod, ]
  if (!nrow(r)) return(FALSE)
  isTRUE(r$median_pct_rho >= UC_PCT_MIN) || isTRUE(r$median_pct_pc1 >= UC_PCT_MIN)
}
cat(sprintf("\n  criterion: median percentile >= %.0f%% on either measure\n\n", 100 * UC_PCT_MIN))
for (nm in uc_sum$module)
  cat(sprintf("    %-34s %s\n", nm, if (uc_passes(nm)) "COHERENT" else "not coherent"))

ctrl_ok_25 <- uc_passes("Hypoxia") || uc_passes("CollagenSynthesis")
cat("\n", strrep("-", 74), "\n", sep = "")
cat("  VERDICT\n")
cat(strrep("-", 74), "\n", sep = "")
if (!ctrl_ok_25) {
  cat("  The positive controls failed in this run, so nothing here is readable.\n")
  if (exists("coh_sum")) {
    s19 <- vapply(c("Hypoxia", "CollagenSynthesis"), passes, logical(1))
    cat(sprintf("  In section 19: hypoxia %s, collagen %s.%s\n",
                if (s19[1]) "coherent" else "not coherent", if (s19[2]) "coherent" else "not coherent",
                if (any(s19)) " So the null pool or matching here is at fault, not the candidate sets." else ""))
  }
} else {
  upr_names <- names(UPR_SETS)
  coh  <- upr_names[vapply(upr_names, uc_passes, logical(1))]
  ours <- uc_passes("ERstressStrict")
  cat(sprintf("  our module (ERstressStrict) : %s\n", if (ours) "COHERENT" else "not coherent"))
  cat(sprintf("  published sets coherent     : %d of %d\n", length(coh), length(upr_names)))
  s24 <- if (exists("UPR_VERDICT")) UPR_VERDICT else NA_character_
  if (!ours && !length(coh)) {
    cat("\n  NO ER-STRESS DEFINITION IS COHERENT IN THIS DATA, while at least one\n")
    cat("  positive control is. Incoherence is then a property of ER-stress gene\n")
    cat("  sets in Visium data, not a defect of the hand-built module: the score is\n")
    cat("  a composite readout, not one coordinated programme.\n")
    cat(sprintf("  Section 24 verdict on whether the definitions agree: %s.\n",
                if (is.na(s24)) "not available" else s24))
  } else if (!ours && length(coh)) {
    cat(sprintf("\n  *** SOME PUBLISHED SETS COHERE AND OURS DOES NOT: %s ***\n",
                paste(coh, collapse = ", ")))
    cat(sprintf("  Section 24 verdict: %s. If those sets also reproduce the result, the\n",
                if (is.na(s24)) "not available" else s24))
    cat("  primary definition should probably move to one of them, with the hand-\n")
    cat("  built module kept as a sensitivity analysis.\n")
  } else {
    cat("\n  Our module coheres in this run. Compare with section 19 before using\n")
    cat("  either number: the null pool is not the same.\n")
  }
}

cat("\nCAVEATS\n")
cat("  Coherence is necessary for a module score to mean one thing, not\n")
cat("  sufficient, and its absence does not by itself invalidate a result that\n")
cat("  replicates across patients, gene sets and definitions.\n")
cat(strrep("=", 74), "\n", sep = "")

log_msg("Section 25 complete.")

################################################################################
## 26. SPATIAL PERMUTATION TEST
##
##   Neighbouring spots are not independent, so model CIs are too optimistic.
##   For every section the map of P4HA1 detection is moved as one rigid sheet
##   (a random shift on the Visium grid with wrap-around, and a random 180
##   degree turn). Patches stay patches, but no longer sit on their true
##   positions, so any real link with the predictors is broken. The section
##   15 model is refitted and pooled exactly as in section 15. Repeating this
##   SP_N_PERM times gives the null distribution of the pooled slope
##   difference under "same spatial structure, no real link".
##   The unpermuted fit on all spots must reproduce section 15, section by
##   section. The smallest possible empirical p is 1 / (SP_N_PERM + 1).
################################################################################

stopifnot(exists("MODULES"), exists("work"), exists("SEED"), exists("ADJUST_DEPTH"),
          exists("meta_pool"), exists("write_out"), exists("assert_run_state"),
          exists("log_msg"), exists("min_genes_for"), exists("check_against_s15"))

## ============================================================================
## 26.0  SETTINGS
## ============================================================================
SP_OUTCOME      <- "P4HA1"
SP_N_PERM       <- 1000
SP_SEED         <- 20261005
SP_PROBE        <- 20
SP_MIN_SHIFT    <- 10          # smallest allowed move, in spot spacings
SP_MODS         <- c("Hypoxia", "CollagenSynthesis", "ERstressStrict", "iCAF", "ERAD")
SP_RENAME       <- c(Hypoxia = "zHyp", CollagenSynthesis = "zCol", ERstressStrict = "zER",
                     iCAF = "ziCAF", ERAD = "zERAD")
SP_MIN_SPOTS    <- 30
SP_MIN_EVENTS   <- 70
SP_ADJUST_DEPTH <- ADJUST_DEPTH
SP_II_COL       <- "zHyp:zCol"
SP_II_ER        <- "zHyp:zER"
SP_RHS <- paste0("zHyp * zCol + zHyp * zER + ziCAF + zERAD",
                 if (SP_ADJUST_DEPTH) " + zDepth" else "")

assert_run_state("Section 26")
for (m in SP_MODS) if (SP_OUTCOME %in% MODULES[[m]])
  stop("Section 26 halted: ", SP_OUTCOME, " is inside ", m)
SP_SECTIONS <- s15_fitted(SP_OUTCOME)

cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 26  SPATIAL PERMUTATION TEST\n")
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  outcome %s, %d permutations, seed %d, %d sections\n",
            SP_OUTCOME, SP_N_PERM, SP_SEED, length(SP_SECTIONS)))

## ============================================================================
## 26.1  ONE SECTION: SCORES, OUTCOME AND GRID POSITION OF EVERY SPOT
## ============================================================================
sp_load <- function(rds_path, sid) {
  obj <- readRDS(rds_path)
  obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
  DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else Assays(obj)[1]
  if (!SP_OUTCOME %in% rownames(obj)) { log_msg("  %s: outcome absent -> skipped", sid); return(NULL) }
  
  for (nm in SP_MODS) {
    full <- MODULES[[nm]]; g <- intersect(full, rownames(obj))
    if (length(g) < min_genes_for(full)) {
      log_msg("  %s: %s has %d/%d genes -> skipped", sid, nm, length(g), length(full)); return(NULL)
    }
    set.seed(SEED)
    obj <- tryCatch(add_module_score_fallback(obj, features = list(g), name = nm,
                                              assay = DefaultAssay(obj), nbin = 20, ctrl = 50),
                    error = function(e) NULL)
    if (is.null(obj) || !paste0(nm, "1") %in% colnames(obj@meta.data)) return(NULL)
    names(obj@meta.data)[names(obj@meta.data) == paste0(nm, "1")] <- nm
  }
  
  vars <- c(SP_OUTCOME, SP_MODS)
  cnt_col <- grep("^nCount_", colnames(obj@meta.data), value = TRUE)
  if (SP_ADJUST_DEPTH) { if (!length(cnt_col)) return(NULL); vars <- c(vars, cnt_col[1]) }
  d <- FetchData(obj, vars = vars, layer = "data")
  
  zz <- function(x) { s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) return(NULL)
  as.numeric((x - mean(x, na.rm = TRUE)) / s) }
  for (nm in SP_MODS) { v <- zz(d[[nm]]); if (is.null(v)) return(NULL); d[[ SP_RENAME[[nm]] ]] <- v }
  if (SP_ADJUST_DEPTH) { v <- zz(log1p(d[[cnt_col[1]]])); if (is.null(v)) return(NULL); d$zDepth <- v }
  d$y <- as.integer(d[[SP_OUTCOME]] > 0)
  
  ## Visium array row and col, matched by barcode. Spots without a position
  ## are kept here (NA) so the safety check uses exactly section 15's spots.
  im <- Images(obj)
  if (!length(im)) { log_msg("  %s: no image -> skipped", sid); return(NULL) }
  co <- tryCatch(obj@images[[im[1]]]@coordinates, error = function(e) NULL)
  if (is.null(co) || !all(c("row", "col") %in% colnames(co))) {
    log_msg("  %s: no row/col grid coordinates -> skipped", sid); return(NULL)
  }
  idx <- match(rownames(d), rownames(co))
  d$grow <- as.integer(co$row[idx]); d$gcol <- as.integer(co$col[idx])
  if (any((d$grow + d$gcol) %% 2L != 0L, na.rm = TRUE))
    stop("Section 26, ", sid, ": row + col is not even for every spot. The grid is not the ",
         "expected Visium layout, so the shift rule below would be wrong.")
  d
}

## ============================================================================
## 26.2  FITS
## ============================================================================
sp_contrast <- function(b, V) {
  need <- c(SP_II_COL, SP_II_ER)
  if (!all(need %in% names(b)) || any(!is.finite(b[need]))) return(NULL)
  L <- setNames(numeric(length(b)), names(b)); L[SP_II_COL] <- 1; L[SP_II_ER] <- -1
  c(est = as.numeric(L %*% b), var = as.numeric(t(L) %*% V %*% L),
    b_col = unname(b[SP_II_COL]), b_er = unname(b[SP_II_ER]))
}
sp_usable <- function(y) {
  ev <- sum(y); length(y) >= SP_MIN_SPOTS && ev >= SP_MIN_EVENTS && (length(y) - ev) >= SP_MIN_EVENTS
}
sp_fit_obs <- function(d) {
  if (!sp_usable(d$y)) return(NULL)
  f <- tryCatch(glm(as.formula(paste("y ~", SP_RHS)), data = d, family = binomial()),
                error = function(e) NULL)
  if (is.null(f) || !f$converged) return(NULL)
  sp_contrast(coef(f), vcov(f))
}
sp_fit_fast <- function(X, y) {
  if (!sp_usable(y)) return(NULL)
  fg <- tryCatch(suppressWarnings(glm.fit(X, y, family = binomial())), error = function(e) NULL)
  if (is.null(fg) || !isTRUE(fg$converged)) return(NULL)
  p <- fg$rank; if (p < ncol(X)) return(NULL)
  cv <- tryCatch(chol2inv(fg$qr$qr[1:p, 1:p, drop = FALSE]), error = function(e) NULL)
  if (is.null(cv)) return(NULL)
  nm <- colnames(X)[fg$qr$pivot[1:p]]; dimnames(cv) <- list(nm, nm)
  b <- fg$coefficients; names(b) <- colnames(X)
  sp_contrast(b, cv[names(b), names(b)])
}

## ============================================================================
## 26.3  THE GRID OF ONE SECTION AND THE MOVE
## ============================================================================
sp_grid <- function(d) {
  R <- max(d$grow) + 1L; if (R %% 2L) R <- R + 1L
  C <- max(d$gcol) + 1L; if (C %% 2L) C <- C + 1L
  sx <- d$gcol * 0.5; sy <- d$grow * sqrt(3) / 2
  near <- matrix(NA_integer_, R, C); direct <- matrix(FALSE, R, C)
  direct[cbind(d$grow + 1L, d$gcol + 1L)] <- TRUE
  for (r in 0:(R - 1L)) {
    cc <- seq.int(r %% 2L, C - 1L, by = 2L)
    dx <- outer(cc * 0.5, sx, "-"); dy <- r * sqrt(3) / 2 - sy
    d2 <- dx^2 + matrix(dy^2, nrow = length(cc), ncol = length(sx), byrow = TRUE)
    near[r + 1L, cc + 1L] <- max.col(-d2, ties.method = "first")
  }
  list(R = R, C = C, near = near, direct = direct)
}
sp_move <- function(d, G) {
  repeat {
    dr <- sample.int(G$R, 1L) - 1L
    dc <- sample.int(G$C, 1L) - 1L
    if ((dr + dc) %% 2L) dc <- (dc + 1L) %% G$C
    flip <- runif(1) < 0.5
    wr <- min(dr, G$R - dr) * sqrt(3) / 2; wc <- min(dc, G$C - dc) * 0.5
    if (flip || sqrt(wr^2 + wc^2) >= SP_MIN_SHIFT) break
  }
  r <- d$grow; k <- d$gcol
  if (flip) { r <- G$R - 1L - r; k <- G$C - 1L - k }
  r <- (r + dr) %% G$R; k <- (k + dc) %% G$C
  src <- G$near[cbind(r + 1L, k + 1L)]
  list(src = src, direct = mean(G$direct[cbind(r + 1L, k + 1L)]))
}

## ============================================================================
## 26.4  LOAD EVERY SECTION ONCE AND FIT THE UNPERMUTED MODEL
## ============================================================================
log_msg("Section 26 start. ADJUST_DEPTH = %s. run_id = %s.", SP_ADJUST_DEPTH, RUN_ID)
sp_dat <- list(); sp_obs <- list()
for (i in seq_len(nrow(work))) {
  path <- work$path[i]; sid <- work$sample_id[i]
  if (is.na(path) || !file.exists(path)) next
  if (!sid %in% SP_SECTIONS) next
  log_msg("[sp load %d/%d] %s", i, nrow(work), sid)
  d <- tryCatch(sp_load(path, sid),
                error = function(e) { if (grepl("^Section 26", conditionMessage(e))) stop(e)
                  log_msg("  skipped: %s", conditionMessage(e)); NULL })
  if (is.null(d)) next
  o_all <- sp_fit_obs(d)                      # all spots, as section 15
  if (is.null(o_all)) { log_msg("  %s: model not fitted -> skipped", sid); next }
  n_all <- nrow(d); has_pos <- !is.na(d$grow) & !is.na(d$gcol)
  n_nopos <- sum(!has_pos)
  if (n_nopos) { log_msg("  %s: %d of %d spots have no grid position and are set aside",
                         sid, n_nopos, n_all); d <- d[has_pos, , drop = FALSE] }
  o <- if (n_nopos) sp_fit_obs(d) else o_all  # positioned spots: what is permuted
  if (is.null(o)) { log_msg("  %s: model not fitted on positioned spots -> skipped", sid); next }
  X <- model.matrix(as.formula(paste("~", SP_RHS)), data = d)
  sp_dat[[sid]] <- list(X = X, y = d$y, d = d[, c("grow", "gcol")], G = sp_grid(d))
  sp_obs[[sid]] <- data.frame(sample_id = sid, est = o[["est"]], var = o[["var"]],
                              b_col = o[["b_col"]], b_er = o[["b_er"]],
                              est_all = o_all[["est"]], var_all = o_all[["var"]],
                              n_spots_all = n_all, n_no_position = n_nopos,
                              n_spots = nrow(d), n_events = sum(d$y), row.names = NULL)
  gc(verbose = FALSE)
}
sp_obs <- do.call(rbind, sp_obs)
if (is.null(sp_obs) || !nrow(sp_obs)) stop("Section 26: nothing was fitted.")
write_out(sp_obs, "SPPERM_observed_per_tumour")

## ---- safety check, section by section --------------------------------------
check_against_s15(data.frame(sample_id = sp_obs$sample_id, est = sp_obs$est_all, var = sp_obs$var_all),
                  SP_OUTCOME, "Section 26")
obs_pool <- meta_pool(sp_obs$est, sp_obs$var, sp_obs$sample_id)
cat(sprintf("  spots without a grid position: %d in %d section(s)\n",
            sum(sp_obs$n_no_position), sum(sp_obs$n_no_position > 0)))
cat(sprintf("  estimate on positioned spots = %+.6f  <- the value the permutations are compared with\n",
            obs_pool$est))
chk <- vapply(names(sp_dat), function(s) {
  f <- sp_fit_fast(sp_dat[[s]]$X, sp_dat[[s]]$y)
  if (is.null(f)) NA_real_ else abs(f[["est"]] - sp_obs$est[sp_obs$sample_id == s]) }, numeric(1))
if (anyNA(chk) || max(chk) > 1e-6)
  stop("Section 26 halted: the fast fitter does not match glm on the unpermuted data.")
cat("  fast fitter PASSED (agrees with glm, largest difference ", format(max(chk), digits = 2), ")\n", sep = "")

## ============================================================================
## 26.5  PERMUTATIONS
## ============================================================================
set.seed(SP_SEED)
sids <- names(sp_dat); K <- length(sids)
null_est <- matrix(NA_real_, SP_N_PERM, K, dimnames = list(NULL, sids))
null_var <- null_est
null_pool <- data.frame(perm = seq_len(SP_N_PERM), est = NA_real_, k = NA_integer_,
                        pooled_b_col = NA_real_, pooled_b_er = NA_real_, direct = NA_real_)
t0 <- Sys.time()
for (b in seq_len(SP_N_PERM)) {
  bc <- rep(NA_real_, K); be <- bc; dv <- bc
  for (j in seq_len(K)) {
    S <- sp_dat[[j]]
    mv <- sp_move(S$d, S$G)
    f <- sp_fit_fast(S$X, S$y[mv$src])
    dv[j] <- mv$direct
    if (is.null(f)) next
    null_est[b, j] <- f[["est"]]; null_var[b, j] <- f[["var"]]
    bc[j] <- f[["b_col"]]; be[j] <- f[["b_er"]]
  }
  p <- meta_pool(null_est[b, ], null_var[b, ], sids)
  null_pool$est[b] <- p$est; null_pool$k[b] <- p$k
  null_pool$pooled_b_col[b] <- mean(bc, na.rm = TRUE)   # simple means, description only
  null_pool$pooled_b_er[b]  <- mean(be, na.rm = TRUE)
  null_pool$direct[b] <- mean(dv)
  if (b == SP_PROBE) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    log_msg("  %d permutations took %.0f s. Estimated total: %.0f minutes.",
            SP_PROBE, el, el / SP_PROBE * SP_N_PERM / 60)
  }
  if (b %% 100 == 0) log_msg("  permutation %d / %d", b, SP_N_PERM)
}
write_out(null_pool, "SPPERM_null_pooled")

## ============================================================================
## 26.6  REPORT
## ============================================================================
nl <- null_pool$est[is.finite(null_pool$est)]; B <- length(nl)
obs <- obs_pool$est
p_one <- (1 + sum(nl <= obs)) / (B + 1)                 # as negative or more
p_two <- (1 + sum(abs(nl) >= abs(obs))) / (B + 1)
q <- quantile(nl, c(0.025, 0.5, 0.975))

per_t <- do.call(rbind, lapply(sids, function(s) {
  v <- null_est[, s]; v <- v[is.finite(v)]; o <- sp_obs$est[sp_obs$sample_id == s]
  data.frame(sample_id = s, patient = patient_of(s), observed = o,
             null_mean = mean(v), null_sd = sd(v),
             model_se = sqrt(sp_obs$var[sp_obs$sample_id == s]),
             p_one_sided = (1 + sum(v <= o)) / (length(v) + 1), n_perm_ok = length(v))
}))
per_t$se_ratio <- per_t$null_sd / per_t$model_se
write_out(per_t, "SPPERM_per_tumour")

sp_summary <- data.frame(outcome = SP_OUTCOME, observed = obs, k = obs_pool$k,
                         k_sections = obs_pool$k_sections,
                         model_lower = obs_pool$lower, model_upper = obs_pool$upper,
                         n_perm = B, null_mean = mean(nl), null_sd = sd(nl),
                         null_q025 = q[[1]], null_median = q[[2]], null_q975 = q[[3]],
                         p_one_sided = p_one, p_two_sided = p_two,
                         median_k_in_null = median(null_pool$k, na.rm = TRUE),
                         mean_direct_share = mean(null_pool$direct, na.rm = TRUE),
                         median_se_ratio = unit_median(per_t$se_ratio, per_t$sample_id),
                         seed = SP_SEED)
write_out(sp_summary, "SPPERM_summary")

cat("\n", strrep("=", 74), "\n", sep = "")
cat(sprintf("SECTION 26  RESULT  run_id %s\n", RUN_ID))
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  observed pooled slope difference : %+.4f  (model CI %+.4f to %+.4f, %s)\n",
            obs, obs_pool$lower, obs_pool$upper, k_lab(obs_pool)))
cat(sprintf("  permutations usable              : %d of %d\n", B, SP_N_PERM))
cat(sprintf("  null: mean %+.4f, sd %.4f, middle 95%% %+.4f to %+.4f\n", mean(nl), sd(nl), q[[1]], q[[3]]))
cat(sprintf("  permuted estimates as negative or more : %d\n", sum(nl <= obs)))
cat(sprintf("  empirical p, one sided           : %.4f\n", p_one))
cat(sprintf("  empirical p, two sided           : %.4f\n", p_two))
cat(sprintf("  smallest p this run can give     : %.4f\n", 1 / (B + 1)))
cat(sprintf("  %s used per permutation (median)   : %.0f of %d\n",
            K_UNIT, median(null_pool$k, na.rm = TRUE), n_units(sids)))
cat(sprintf("  spots that kept a true tissue value (mean share) : %.0f%%\n",
            100 * mean(null_pool$direct, na.rm = TRUE)))
cat(sprintf("  null sd / model SE, median over %ss : %.2f  (1 = model SE honest, above 1 = too small)\n",
            UNIT_LABEL, unit_median(per_t$se_ratio, per_t$sample_id)))
cat(sprintf("  sections with own one sided p < 0.05 : %d of %d\n", sum(per_t$p_one_sided < 0.05), K))

cat("\n  VERDICT\n")
if (p_two < 0.05) {
  cat("  The observed estimate lies outside what moved maps produce.\n")
  cat("  Spatial structure of the outcome alone does not explain it.\n")
} else {
  cat("  *** The observed estimate is NOT distinguishable from moved maps. ***\n")
  cat("  Spatial structure of the outcome alone can produce a value this large.\n")
}
cat("\n  CAVEATS\n")
cat("  Only the outcome map is moved; the predictors keep their own structure.\n")
cat("  Wrap-around and nearest-spot filling make the kept structure approximate.\n")
cat("  The model CI above is the uncorrected one and is shown for comparison only.\n")
cat(strrep("=", 74), "\n", sep = "")
log_msg("Section 26 complete.")
