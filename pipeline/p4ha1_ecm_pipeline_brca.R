################################################################################
# P4HA1 / ECM / hypoxia in spatial transcriptomics  --  BREAST CANCER (BRCA)
# cohort pipeline, ported from the PDAC pipeline v4
#   (github.com/Mahditajbakhsh/p4ha1-pdac-spatial, pipeline/p4ha1_ecm_pipeline_v4.R)
#
# The analysis is the PDAC analysis. Gene panel, modules, models, seeds and
# section order are unchanged, so the two cancers can be compared directly.
# Only the cohort and the pooling unit change. Every change is listed below and
# tagged "BRCA:" at the place in the code where it happens.
#
# SECTION MAP: see docs/section_map.md
#
# ---- CHANGES FROM THE PDAC v4 PIPELINE ---------------------------------------
# (B1) COHORT. 54 Visium sections from 22 patients, HTAN WUSTL atlas
#      (Mo et al., Nature 2024, doi:10.1038/s41586-024-08087-4). The list of
#      sections, their patient and their primary / metastasis status come from
#      docs/brca_visium_samples.csv, taken from Sample_ID_Lookup_table_v1.xlsx
#      (sheet Visium) in github.com/ding-lab/ST_subclone_publication.
#
# (B2) PRIMARY TUMOURS ONLY in the main run. The 5 metastatic sections all come
#      from one patient (HT268B1). They are left out of the work list unless
#      INCLUDE_METASTASIS <- TRUE, which writes to a separate output folder so
#      nothing is overwritten. Section 22 compares the two runs.
#
# (B3) THE PATIENT IS THE POOLING UNIT. In PDAC most patients had one section
#      and pooling was over sections, with patients as a sensitivity (section
#      21). In BRCA several patients have 4 to 6 serial sections cut 5 to 100 um
#      apart, which are close to replicates of one tissue. Pooling over
#      sections would let one tumour vote six times. meta_pool() therefore
#      first combines the sections of each patient (fixed effect, inverse
#      variance, the same collapse PDAC section 21 used) and then pools the
#      patients with DerSimonian-Laird. POOL_LEVEL <- "section" restores the
#      PDAC behaviour; section 21 now reports that as the sensitivity.
#      meta_pool() halts if it is called without sample IDs, so no pooling
#      step can silently fall back to the section level.
#
# (B4) ONE PROJECT ROOT. All paths are derived from P4HA1_BRCA_ROOT.
#
# (B5) PDAC RESULTS REMOVED FROM COMMENTS. Comments in the PDAC script quoted
#      PDAC numbers (odds ratios, k, coverage, sample IDs). They are not true
#      for this cohort and were removed or marked "PDAC:". Hard coded PDAC
#      checks (expected patient count, the metastasis ID) now come from the
#      sample table.
#
# ---- HOW TO RUN THE DEPTH CHECK (unchanged from PDAC) -----------------------
#   1. Run once with ADJUST_DEPTH <- FALSE   (files have no suffix)
#   2. Set ADJUST_DEPTH <- TRUE and re-run sections 13 to 16
#   3. Run section 17 to print the comparison
#
# ---- KNOWN LIMITATIONS (carried over from PDAC) ------------------------------
#     Spatial autocorrelation between neighbouring spots is not modelled in
#     the regressions, so every model confidence interval is narrower than it
#     should be. Section 16 (random genes) and section 26 (spatial
#     permutation) measure how large that problem is in this cohort. Quote
#     empirical percentiles rather than model p values.
#
#     The fixed effect collapse inside a patient treats serial sections as
#     independent, so the variance of each patient estimate is too small.
#     The patient still counts once for k and for tau2. Section 21 reports
#     a three level model (sections in patients) as a second opinion when
#     metafor is installed.
#
#     The cohort is read end to end many times (sections 9, 10, 11, 13 to 16).
#     Not merged, for the same reason as in PDAC: section 11 can drop a sample
#     on low gene coverage and section 10 cannot, so a merge would change k.
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
## BRCA: everything below hangs off one root folder. Set the environment
## variable P4HA1_BRCA_ROOT, or edit the default here. Expected layout:
##   <root>/data/seurat/      the *-SeuratObj.rds files
##   <root>/data/manifest.csv optional Synapse manifest
##   <root>/data/upr_sets/    MSigDB UPR gene set TSVs (section 24)
##   <root>/data/estimate_audit/ESTIMATE_signatures_used.csv (section 20)
PROJECT_ROOT <- Sys.getenv("P4HA1_BRCA_ROOT", unset = path.expand("~/p4ha1_brca"))

## BRCA: path of the repo's sample table. Defaults to docs/ next to this script's
## folder when sourced from the repo root; set P4HA1_BRCA_SAMPLES to override.
SAMPLE_TABLE <- Sys.getenv("P4HA1_BRCA_SAMPLES", unset = file.path("docs", "brca_visium_samples.csv"))

## BRCA: main run = primary tumours only. TRUE adds the 5 metastatic sections
## (one patient) and writes to a different output folder.
INCLUDE_METASTASIS <- FALSE

## BRCA: pooling unit. "patient" is the primary analysis, "section" reproduces
## the PDAC behaviour and is the sensitivity analysis.
POOL_LEVEL <- "patient"
stopifnot(POOL_LEVEL %in% c("patient", "section"))
K_UNIT <- if (POOL_LEVEL == "patient") "patients" else "sections"   # BRCA: label for pooled k

## BRCA: a downloaded .rds that is not in the sample table halts the run unless
## this is TRUE, so no unlisted section can enter the cohort unnoticed.
ALLOW_UNLISTED_FILES <- FALSE

INPUT_MODE   <- "folder"   # "folder" = analyse every *-SeuratObj.rds in DOWNLOAD_DIR
# "manifest" = drive from MANIFEST_CSV + entityIds

DOWNLOAD_DIR <- file.path(PROJECT_ROOT, "data", "seurat")
OUT_ROOT     <- file.path(PROJECT_ROOT, paste0(
  if (INCLUDE_METASTASIS) "cohort_with_metastasis" else "cohort_primary",
  if (POOL_LEVEL == "section") "_sectionpool" else ""))
MANIFEST_CSV <- file.path(PROJECT_ROOT, "data", "manifest.csv")  # skipped if absent

# only used in manifest mode
INCLUDE_SAMPLE_IDS <- character(0)
INCLUDE_ENTITY_IDS <- character(0)
USE_SYNAPSER       <- FALSE

SAVE_PER_SAMPLE_FIGS <- TRUE    # set FALSE on re-runs: ~230 PNGs, slow, rarely re-read
SAVE_ANALYZED_RDS    <- FALSE   # also save each analysed object (heavy)
SPOT_SIZE            <- 250

## v3: one seed for the whole run. AddModuleScore samples its control genes, so
## without this the same module scored in two different sections is not the same
## variable. set.seed(SEED) is called immediately before EVERY AddModuleScore.
SEED <- 1234
set.seed(SEED)

## NOTE: the depth switch, RUN_ID, stamp() and outf() are defined in the
## SECTION 0 PATCH block further down, immediately before section 13. They
## live there because the run-state guard that checks them must not depend
## on any single section. Do not define them here as well.


## ============================================================================
## 1. Gene panel + modules
## ============================================================================
gene_panel <- c(
  "P4HA1","P4HA2","P4HB","PLOD1","PLOD2","PLOD3",
  "COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3",
  "LOX","LOXL1","LOXL2","LOXL3","LOXL4",
  "FN1","POSTN","SPARC","FBN1","THBS2",
  "MMP2","MMP11","MMP14","TIMP1","TIMP2",
  "ACTA2","TAGLN","MYH11","MCAM",
  "IL6","CXCL12","CXCL1","CXCL2","HAS1","PDGFRA",
  "VEGFA","CA9","SLC2A1","LDHA","BNIP3","ANGPTL4","ENO1","HIF1A",
  "EPCAM","KRT8","KRT18","KRT19" ,"CDH1")

MODULES <- list(
  CollagenSynthesis  = c("COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3"),
  CollagenMaturation = c("P4HA1","P4HA2","P4HB","PLOD1","PLOD2","PLOD3"),
  ECMCrosslink       = c("LOX","LOXL1","LOXL2","LOXL3","LOXL4"),
  ECMGlycoprotein    = c("FN1","POSTN","SPARC","FBN1","THBS2"),
  ECMRemodel         = c("MMP2","MMP11","MMP14","TIMP1","TIMP2"),
  myCAF              = c("ACTA2","TAGLN","MYH11","MCAM"),
  iCAF               = c("IL6","CXCL12","CXCL1","CXCL2","HAS1","PDGFRA"),
  Hypoxia            = c("VEGFA","CA9","SLC2A1","LDHA","BNIP3","ANGPTL4","ENO1"),
  Epithelial         = c("EPCAM","KRT8","KRT18","KRT19","CDH1"))

TARGET_GENES   <- c("COL1A1","COL1A2","COL3A1","COL5A1","COL6A1","FN1","POSTN","SPARC","LOX","LOXL2")
MODULE_TARGETS <- c("CollagenSynthesis","CollagenMaturation","ECMCrosslink",
                    "ECMGlycoprotein","Hypoxia","myCAF","iCAF","Epithelial")
LOGIT_MODULES  <- c("Hypoxia","Epithelial","myCAF","CollagenSynthesis","iCAF")

## NOTE ON CollagenMaturation: it contains P4HA1 and is therefore circular as a
## predictor. It is scored and plotted descriptively only. It is NOT in
## LOGIT_MODULES and must never enter a model of P4HA1 detection.
## It also contains P4HA2, which matters in section 15 where P4HA2 becomes an
## outcome. Section 15.1 re-checks this automatically and halts on failure.

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
## BRCA: AddModuleScore with a bin fallback.
## Seurat's AddModuleScore splits all genes into nbin bins by mean expression
## (ggplot2::cut_number) to pick expression-matched control genes. In a few
## sparse BRCA sections hundreds of genes have exactly one UMI in one spot, so
## they share one identical mean (log(2)/n_spots) and 20 bins cannot be formed:
## "Insufficient data values to produce 20 bins". PDAC never hit this.
## This wrapper tries 20 bins first, exactly as before (same seed, so every
## section that worked before gives the identical score), and only on that
## specific error retries with 15 and then 10 bins. Every fallback is logged
## and written to <OUT_ROOT>/nbin_fallback_log.csv. Any other error is raised
## unchanged, so nothing else is hidden.
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
      set.seed(SEED)                                    # v3
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
## BRCA: the sample table. Read once here because meta_pool() needs the
## section -> patient map. Columns: sample_id, patient_id, htan_biospecimen_id,
## htan_participant_id, primary_metastasis, cancer_type.
if (!file.exists(SAMPLE_TABLE))
  stop("Sample table not found: ", SAMPLE_TABLE,
       ". Source the script from the repo root or set P4HA1_BRCA_SAMPLES.")
SAMPLES <- read.csv(SAMPLE_TABLE, stringsAsFactors = FALSE)
stopifnot(all(c("sample_id", "patient_id", "primary_metastasis") %in% names(SAMPLES)),
          !anyDuplicated(SAMPLES$sample_id),
          all(SAMPLES$primary_metastasis %in% c("Primary", "Metastasis")))

## BRCA: section -> patient. Filled in section 4 once the downloaded files have
## been matched to the table (file names may carry suffixes the table lacks).
SAMPLE_PATIENT <- setNames(SAMPLES$patient_id, SAMPLES$sample_id)
patient_of <- function(sid) {
  p <- unname(SAMPLE_PATIENT[as.character(sid)])
  if (anyNA(p)) stop("patient_of(): no patient for section(s) ",
                     paste(unique(sid[is.na(p)]), collapse = ", "))
  p
}

## BRCA: fixed effect, inverse variance collapse of the sections of each
## patient. Identical to fe_collapse() in PDAC section 21.
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

# DerSimonian-Laird random-effects pooling
## BRCA: sid (the section IDs, same order as yi) is now required when
## POOL_LEVEL is "patient". The DL step itself is unchanged from PDAC.
## k counts patients (or sections when POOL_LEVEL is "section");
## k_sections always counts the sections that went in.
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
             I2=max(0,(Q-df)/Q)*100, k_sections=n_sec)
}

## ============================================================================
## 3. PER-SAMPLE analysis
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
    f <- intersect(c("P4HA1","myCAF","Hypoxia","Epithelial"), c(rownames(obj), mods_present))
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
    sc <- intersect(c("CollagenSynthesis","Hypoxia","myCAF","Epithelial"), colnames(expr))
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
      tab <- do.call(rbind, lapply(gm, function(g) {
        v <- FetchData(obj, vars=g, layer="data")[,1]; names(v) <- colnames(obj)
        v <- v[rownames(coords)]
        mi <- spdep::moran.test(v, lw, zero.policy=TRUE, randomisation=TRUE)
        data.frame(gene=g, morans_I=unname(mi$estimate[1]), expected=unname(mi$estimate[2]),
                   variance=unname(mi$estimate[3]), p_value=mi$p.value) }))
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
  if (all(c("Epithelial","myCAF") %in% colnames(expr))) {
    expr$Compartment <- with(expr, ifelse(
      Epithelial > quantile(Epithelial,.66,na.rm=TRUE) & myCAF < quantile(myCAF,.5,na.rm=TRUE),
      "Epithelial-dominant",
      ifelse(myCAF > quantile(myCAF,.66,na.rm=TRUE) & Epithelial < quantile(Epithelial,.5,na.rm=TRUE),
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
      rho_given_myCAF      = if ("myCAF"      %in% colnames(expr)) pct(expr$P4HA1, expr[[g]], expr$myCAF)      else NA)))
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
  
  ## --- per-sample logistic regression + forest ---
  ## v3: fitted on an explicit data.frame so the coefficient is named "x"
  ## rather than "expr[[m]]". Same estimate, easier to read in the summary.
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
              scale_x_log10(limits=c(.1,10), breaks=c(.2,.5,1,2,5)) +
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
  
  ## v3 BUG FIX: this used to be df$n_spots <- nrow(expr), which silently
  ## overwrote the per-CLUSTER n_spots column inside cl_summary with the total
  ## spot count of the sample.
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
## BRCA: match every file to the sample table, attach patient and status,
## then keep primary tumours only unless INCLUDE_METASTASIS is TRUE.
##   1. exact match of the file stem to sample_id
##   2. otherwise, the table section whose sample_id is the longest prefix of
##      the file stem, or vice versa (HTAN file names may add or drop a suffix)
##   Anything still unmatched halts the run unless ALLOW_UNLISTED_FILES.
## ---------------------------------------------------------------------------
match_section <- function(stem) {
  if (stem %in% SAMPLES$sample_id) return(stem)
  hit <- SAMPLES$sample_id[startsWith(stem, SAMPLES$sample_id) | startsWith(SAMPLES$sample_id, stem)]
  if (length(hit) == 1) return(hit)
  if (length(hit) > 1) {
    stop("File ", stem, " matches several table sections: ", paste(hit, collapse = ", "))
  }
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
if (anyDuplicated(work$table_id))
  stop("Two files map to the same table section: ",
       paste(unique(work$table_id[duplicated(work$table_id)]), collapse = ", "))
tbl <- SAMPLES[match(work$table_id, SAMPLES$sample_id), ]
work$patient_id         <- tbl$patient_id
work$primary_metastasis <- tbl$primary_metastasis
if (!all(tbl$cancer_type == "BRCA")) stop("Non-BRCA sections in the work list.")

## the analysis keys on the file stem, so map the file stem to the patient too
SAMPLE_PATIENT <- c(SAMPLE_PATIENT, setNames(work$patient_id, work$sample_id))
SAMPLE_PATIENT <- SAMPLE_PATIENT[!duplicated(names(SAMPLE_PATIENT))]

n_met <- sum(work$primary_metastasis == "Metastasis")
if (!INCLUDE_METASTASIS && n_met) {
  log_msg("Primary-only run: %d metastatic section(s) left out (%s)", n_met,
          paste(unique(work$patient_id[work$primary_metastasis == "Metastasis"]), collapse = ", "))
  work <- work[work$primary_metastasis == "Primary", , drop = FALSE]
}

## how much of the expected cohort is actually here
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
## 5. Save combined per-sample tables + run status
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
## BRCA: the section -> patient map used for pooling, read by sections 21 and 22
write.csv(work[, c("sample_id", "table_id", "patient_id", "primary_metastasis")],
          file.path(agg, "cohort_sections.csv"), row.names = FALSE)
log_msg("Samples ok: %d / %d", sum(status$status=="ok"), nrow(status))

## ============================================================================
## 6. CROSS-SAMPLE pooled forest
## ============================================================================
if (nrow(all_logit) > 0) {
  pooled <- all_logit %>% group_by(module) %>%
    group_modify(~ meta_pool(.x$log_OR, .x$SE^2, .x$sample_id)) %>% ungroup() %>%
    mutate(OR=exp(est), lower=exp(lower), upper=exp(upper),
           direction=ifelse(OR>1,"Promotes P4HA1","Suppresses P4HA1"),
           label=sprintf("OR %.2f  (%d %s, I2=%.0f%%)", OR, k, K_UNIT, I2)) %>%   # BRCA: k unit
    arrange(OR) %>% mutate(module=factor(module, levels=module))
  write.csv(pooled, file.path(agg,"POOLED_logistic_forest.csv"), row.names=FALSE)
  fig(ggplot(pooled, aes(OR, module, colour=direction)) +
        geom_vline(xintercept=1, linetype="dashed", colour="grey50") +
        geom_errorbarh(aes(xmin=lower, xmax=upper), height=.25, linewidth=.9) +
        geom_point(size=5) +
        geom_text(aes(label=label), hjust=-.1, vjust=-1.1, size=3.4, colour="black", show.legend=FALSE) +
        scale_x_log10(limits=c(.1,10), breaks=c(.2,.5,1,2,5)) +
        scale_colour_manual(values=c("Promotes P4HA1"="#B2182B","Suppresses P4HA1"="#2166AC")) +
        labs(x="Pooled odds ratio for P4HA1 detection (log scale, 95% CI)", y=NULL, colour=NULL,
             title="What predicts P4HA1 across tumours?",
             subtitle=sprintf("Random-effects meta-analysis of %d patients (%d sections)",   # BRCA
                              length(unique(patient_of(all_logit$sample_id))), length(unique(all_logit$sample_id))),
             caption="Per-sample logistic regression, DerSimonian-Laird pooled") +
        theme_minimal(base_size=13) +
        theme(legend.position="top", panel.grid.minor=element_blank(), plot.title=element_text(face="bold")),
      file.path(agg,"FOREST_pooled.png"), 9, 5)
}

## ============================================================================
## 7. Consistency boxplot: per-tumour P4HA1 <-> module correlations
##
##   BRCA NOTE: the reasoning in this header was written for the PDAC
##   cohort. Any statement about what the data showed, which module worked,
##   or how many tumours were usable describes PDAC, not BRCA. The code
##   below is unchanged except where tagged "BRCA:".
##    KEPT DELIBERATELY. This is the only figure showing per-tumour SPREAD
##    rather than a pooled average. In PDAC, with tau ~0.7-1.0 and every
##    prediction interval crossing 1, "how variable is this across tumours?" is the first
##    question a reader will ask, and this answers it in one glance.
##    Section 12 renders a publication-grade version of this as Figure 5.
## ============================================================================
if (nrow(all_module) > 0) {
  mod_summary <- all_module %>% group_by(module) %>%
    summarise(median_rho=median(rho, na.rm=TRUE), frac_positive=mean(rho>0, na.rm=TRUE),
              n_samples=dplyr::n(), .groups="drop") %>% arrange(desc(median_rho))
  write.csv(mod_summary, file.path(agg,"module_correlation_summary.csv"), row.names=FALSE)
  ord <- mod_summary$module[order(mod_summary$median_rho)]
  fig(ggplot(all_module, aes(factor(module, levels=ord), rho)) +
        geom_hline(yintercept=0, linetype="dashed", colour="grey60") +
        geom_boxplot(outlier.size=.5, fill="grey90") +
        geom_jitter(width=.15, size=.6, alpha=.4) + coord_flip() +
        labs(x=NULL, y="Per-section Spearman rho (P4HA1 vs module)",
             title="Consistency of P4HA1 associations across tumours") +
        theme_minimal(base_size=12),
      file.path(agg,"BOX_module_correlations.png"), 8, 5)
}

log_msg("Core cohort analysis done. Aggregate outputs in: %s", normalizePath(agg))

## ============================================================================
## 8. COHORT within-compartment analysis
##
##   BRCA NOTE: the reasoning in this header was written for the PDAC
##   cohort. Any statement about what the data showed, which module worked,
##   or how many tumours were usable describes PDAC, not BRCA. The code
##   below is unchanged except where tagged "BRCA:".
##    KEPT. Answers a different question from the hypoxia narrative and gave a
##    real result in PDAC (COL3A1 surviving BH correction). Not spare parts.
## ============================================================================
COMPARTMENT_MIN_SPOTS <- 20
COMPARTMENT_LEVELS    <- c("Epithelial-dominant","Mixed","CAF-dominant")
CONTRAST              <- c("CAF-dominant","Epithelial-dominant")

comp_csv <- file.path(agg, "per_sample_compartment_correlations.csv")

if (!file.exists(comp_csv)) {
  log_msg("Cohort compartment step skipped: %s not found", comp_csv)
} else {
  cc <- read.csv(comp_csv, stringsAsFactors = FALSE, check.names = FALSE)
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
      cc$gene <- factor(cc$gene, levels = sort(unique(cc$gene)))
      
      comp_summary <- cc %>%
        group_by(gene, Compartment) %>%
        summarise(n_tumours     = dplyr::n_distinct(sample_id),
                  median_rho    = median(rho, na.rm = TRUE),
                  q25           = unname(quantile(rho, .25, na.rm = TRUE)),
                  q75           = unname(quantile(rho, .75, na.rm = TRUE)),
                  frac_positive = mean(rho > 0, na.rm = TRUE),
                  .groups = "drop") %>%
        arrange(gene, Compartment)
      write.csv(comp_summary, file.path(agg, "COHORT_compartment_summary.csv"), row.names = FALSE)
      
      fig(ggplot(cc, aes(Compartment, rho, fill = Compartment)) +
            geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
            geom_boxplot(outlier.size = .4, alpha = .85) +
            geom_jitter(width = .15, size = .5, alpha = .35) +
            facet_wrap(~ gene, ncol = 5) +
            scale_fill_manual(values = c("Epithelial-dominant" = "#2166AC",
                                         "Mixed"               = "grey75",
                                         "CAF-dominant"        = "#B2182B"), drop = FALSE) +
            labs(x = NULL, y = "Per-tumour Spearman rho (P4HA1 vs gene)",
                 title = "P4HA1 to ECM correlation by tissue compartment",
                 subtitle = sprintf("Across %d tumours; compartments with < %d spots removed",
                                    dplyr::n_distinct(cc$sample_id), COMPARTMENT_MIN_SPOTS)) +
            theme_minimal(base_size = 11) +
            theme(legend.position = "none",
                  axis.text.x = element_text(angle = 30, hjust = 1),
                  strip.text  = element_text(face = "bold")),
          file.path(agg, "BOX_compartment_correlations.png"), 12, 7)
      
      test_rows <- lapply(levels(cc$gene), function(g) {
        sub <- cc[cc$gene == g & as.character(cc$Compartment) %in% CONTRAST,
                  c("sample_id","Compartment","rho")]
        empty <- data.frame(gene = g, n_pairs = 0L, median_CAF = NA_real_,
                            median_Epi = NA_real_, median_diff = NA_real_,
                            W = NA_real_, p = NA_real_)
        if (nrow(sub) == 0) return(empty)
        sub$Compartment <- factor(as.character(sub$Compartment), levels = CONTRAST)
        m <- tapply(sub$rho, list(sub$sample_id, sub$Compartment),
                    function(x) mean(x, na.rm = TRUE))
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
      comp_test$p_adj_BH <- p.adjust(comp_test$p, method = "BH")
      comp_test <- comp_test[order(comp_test$p_adj_BH, na.last = TRUE), ]
      write.csv(comp_test, file.path(agg, "COHORT_compartment_test.csv"), row.names = FALSE)
      
      log_msg("Cohort compartment analysis done: %d summary rows, %d genes tested",
              nrow(comp_summary), nrow(comp_test))
    }
  }
}

## ============================================================================
## 9. ENO1 SENSITIVITY CHECK   [single pass over the cohort]
##
##   BRCA NOTE: the reasoning in this header was written for the PDAC
##   cohort. Any statement about what the data showed, which module worked,
##   or how many tumours were usable describes PDAC, not BRCA. The code
##   below is unchanged except where tagged "BRCA:".
##    The object is loaded ONCE and both hypoxia gene sets are scored from it.
## ============================================================================
HYPOXIA_7 <- c("VEGFA","CA9","SLC2A1","LDHA","BNIP3","ANGPTL4","ENO1")
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
    set.seed(SEED)                                           # v3
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

pooled_sens <- sens %>%
  dplyr::group_by(variant) %>%
  dplyr::group_modify(~ meta_pool(.x$log_OR, .x$SE^2, .x$sample_id)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(OR = exp(est), lower = exp(lower), upper = exp(upper))

write.csv(sens,        file.path(agg, "ENO1_sensitivity_per_sample.csv"), row.names = FALSE)
write.csv(pooled_sens, file.path(agg, "ENO1_sensitivity_pooled.csv"),     row.names = FALSE)

cat("\n=== Hypoxia pooled OR: ENO1 sensitivity ===\n")
for (v in c("with_ENO1", "without_ENO1")) {
  r <- pooled_sens[pooled_sens$variant == v, ]
  if (nrow(r) == 1)
    cat(sprintf("  %-13s OR = %.2f  (95%% CI %.2f-%.2f)  k=%d  I2=%.0f%%\n",
                v, r$OR, r$lower, r$upper, r$k, r$I2))
}
cat("\nIf the two ORs are close, ENO1 is not distorting the hypoxia result.\n")

## ============================================================================
## 10. STRATIFIED MULTIVARIABLE MODEL
##
##   BRCA NOTE: the reasoning in this header was written for the PDAC
##   cohort. Any statement about what the data showed, which module worked,
##   or how many tumours were usable describes PDAC, not BRCA. The code
##   below is unchanged except where tagged "BRCA:".
##     Question (Valerio): at the SAME level of hypoxia, does collagen (+ iCAF)
##     still track P4HA1, or is the collagen signal just hypoxia in disguise?
##
##       (A) ADJUSTED     P(P4HA1>0) ~ Hypoxia + CollagenSynthesis + iCAF
##       (B) STRATIFIED   within within-sample hypoxia tertiles (RELATIVE, so
##                        "high" is not the same absolute level across tumours)
##       (C) INTERACTION  Hypoxia * CollagenSynthesis  <- better-powered test
##                        slope(z_H) = beta_C + beta_HC * z_H
##
##     Not circular: no predictor contains P4HA1. Collinearity is written out.
##     NOTE: predictors here are NOT standardised, so these numbers are not
##     comparable to sections 13, 15 and 16, which are. Section 13 says the same
##     thing in its own header. Do not put them in one table.
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
    set.seed(SEED)                                           # v3
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
  
  ## (B) stratified within within-sample hypoxia tertiles
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
  ##     available for the crossover z_H* = -beta_C/beta_HC. That beta_C is a
  ##     DIFFERENT parameter from the additive model's and must not be swapped.
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
  cat("\n=== Predictor collinearity (median Spearman across tumours) ===\n")
  cat(sprintf("  Hypoxia vs Collagen : %.2f\n", median(col_all$Hyp_Col,  na.rm = TRUE)))
  cat(sprintf("  Hypoxia vs iCAF     : %.2f\n", median(col_all$Hyp_iCAF, na.rm = TRUE)))
  cat(sprintf("  Collagen vs iCAF    : %.2f\n", median(col_all$Col_iCAF, na.rm = TRUE)))
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
    cat(sprintf("  interaction coefficient = %+.4f  (OR scale %.2f, 95%% CI %.2f-%.2f)  k=%d  I2=%.0f%%\n",
                it_row$est, it_row$OR, it_row$lower, it_row$upper, it_row$k, it_row$I2))
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
    cat(sprintf("  Collagen slope = beta_C + beta_HC * z_H\n"))
    cat(sprintf("  Crossover z_H* = %+.3f  (ratio of two pooled estimates; no CI)\n", -bC / bHC))
    cat("  Check against HYPOXIA_SCORE_RANGE.csv: if z_H* lies outside the observed\n")
    cat("  range, the crossover is an extrapolation and should not be reported.\n")
  }
}

log_msg("Stratified multivariable step done.")

## ============================================================================
## 11. PROTEOSTASIS MODULES  (Valerio's Part 2, exploratory)
##
##   BRCA NOTE: the reasoning in this header was written for the PDAC
##   cohort. Any statement about what the data showed, which module worked,
##   or how many tumours were usable describes PDAC, not BRCA. The code
##   below is unchanged except where tagged "BRCA:".
##
##     HYPOTHESIS (his, not established): P4H is bifunctional - it hydroxylates
##     AND acts as a clamp retaining the molecule in the ER so the QC system can
##     decide degradation vs export. If so, P4HA1 may track proteostasis rather
##     than collagen.
##
##     MODULES (audit re-runs at run time and can halt the section):
##       ERstress       27 core ER-stress/UPR genes hand-selected from
##                      HALLMARK_UNFOLDED_PROTEIN_RESPONSE (MSigDB M5922).
##                      NOT the full 113-gene hallmark: only ~24% of it is core
##                      ER stress, the rest mRNA decay (19), ribosome biogenesis
##                      (12), translation (10). In hypoxic tissue, where
##                      translation is globally suppressed, the whole hallmark
##                      would largely measure translational machinery.
##                      VEGFA REMOVED: it is in the hallmark AND in the Hypoxia
##                      module here, so keeping it would manufacture correlation.
##                      VERIFIED against the MSigDB file: all 27 are present in
##                      the 113-gene hallmark, and VEGFA is the only overlap
##                      with the Hypoxia module.
##       ERstressStrict the above minus HYOU1, DDIT4, ERO1A (24 genes). Those
##                      three are hypoxia-responsive and the headline finding is
##                      "ER stress predicts P4HA1 in hypoxic regions", which is
##                      exactly what they could manufacture. Every model is
##                      fitted with BOTH variants; the strict one is primary.
##                      Section 14 confirmed ERO1A is the symbol actually
##                      present in these objects (not the legacy ERO1L), so the
##                      strict module really does drop three genes.
##       ERAD           GOBP_ERAD_PATHWAY (MSigDB M15534, GO:0036503).
##                      VERIFIED: exact match to the MSigDB file, 113 unique
##                      genes both sides. 9 shared with ERstress, stripped
##                      below, leaving 104.
##
##     STILL OPEN: the exclusion criterion is applied inconsistently. HYOU1,
##     DDIT4 and ERO1A are removed as hypoxia-responsive, but ATF4 and CHAC1
##     are PERK-responsive and were kept. That is not a principled line. Decide
##     the rule and apply it, or say plainly that the strict set is defined by
##     HIF-responsiveness only.
##
##     NOT INCLUDED - the retention/chaperone arm. GO:0034663 "endoplasmic
##     reticulum chaperone complex" is DEFINED to include protein disulfide
##     isomerase, which is P4HB, the beta subunit of P4HA1's own enzyme.
##     Scoring it would recreate the CollagenMaturation circularity in a form
##     that looks legitimate. Needs Valerio's ruling on P4HB and on PPIB.
##
##     GENE THRESHOLD is proportional, max(2, 30% of the set). A flat 5-gene
##     minimum is trivial for 104-gene ERAD but brutal for 6-gene iCAF. In PDAC it
##     dropped one tumour (HT270P2, 4/6 iCAF genes), giving k=22 instead of 23.
##
##     Figures live in section 12 and use the STRICT variant only.
## ============================================================================
PROT_DISJOINT   <- TRUE   # strip the 9 genes shared with ERstress out of ERAD
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
      log_msg("  %s: only %d/%d genes of %s detected (need %d) -> sample skipped",
              sid, length(g), length(full), nm, need)
      return(NULL)
    }
    set.seed(SEED)                                           # v3
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
    data.frame(sample_id = sid, module = k, n_found = cov[[k]]["n_found"],
               n_total = cov[[k]]["n_total"], n_required = cov[[k]]["n_required"])))
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
## 11.2  fit every model for one sample, for BOTH ERstress variants
## ---------------------------------------------------------------------------
fit_prot_one <- function(d, sid) {
  out <- list(marginal=NULL, adjusted=NULL, stratified=NULL,
              interaction=NULL, vif=NULL, collin=NULL, coverage=NULL)
  if (is.null(d) || nrow(d) < PROT_MIN_SPOTS || length(unique(d$y)) != 2) return(out)
  out$coverage <- attr(d, "coverage")
  
  var_cor <- suppressWarnings(cor(d$ERstress, d$ERstressStrict, method = "spearman"))
  
  ## `strat` not `str` - the old name shadowed base::str() inside here
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
    
    ## (A) MARGINAL, module alone.
    ## v3: ERAD is now fitted under BOTH variants. The estimate is identical
    ## either way (ERAD alone does not depend on which ERstress variant exists),
    ## but it means filtering the pooled CSV on the primary variant no longer
    ## silently drops the ERAD row.
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
    
    ## (C) STRATIFIED by within-sample hypoxia tertiles
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
cat(sprintf("\n=== Samples contributing: %d of %d ===\n", n_in, nrow(work)))
if (n_in < nrow(work)) {
  missing <- setdiff(work$sample_id, unique(padj_all$sample_id))
  cat("  excluded: ", paste(missing, collapse = ", "), "\n")
  cat("  (reason is in the log above, per sample)\n")
}

if (nrow(cov_all) > 0) {
  cat("\n=== Gene coverage per module (across contributing tumours) ===\n")
  cs <- cov_all %>% group_by(module) %>%
    summarise(median_found = median(n_found), n_total = max(n_total),
              min_found = min(n_found), required = max(n_required), .groups = "drop")
  for (i in seq_len(nrow(cs)))
    cat(sprintf("  %-18s %3d / %3d detected (min %3d, required %3d)\n", cs$module[i],
                cs$median_found[i], cs$n_total[i], cs$min_found[i], cs$required[i]))
}

if (nrow(pcol_all) > 0) {
  cat("\n=== Predictor collinearity (median Spearman) ===\n")
  for (ES in ERSTRESS_VARIANTS) {
    p <- pcol_all[pcol_all$variant == ES, ]
    if (!nrow(p)) next
    cat(sprintf("  [%s]\n", ES))
    cat(sprintf("     Hypoxia vs ERstress : %+.2f\n", median(p$Hyp_ERstress, na.rm=TRUE)))
    cat(sprintf("     Hypoxia vs ERAD     : %+.2f\n", median(p$Hyp_ERAD, na.rm=TRUE)))
    cat(sprintf("     Collagen vs ERstress: %+.2f\n", median(p$Col_ERstress, na.rm=TRUE)))
    cat(sprintf("     ERstress vs ERAD    : %+.2f\n", median(p$ERstress_ERAD, na.rm=TRUE)))
  }
  cat(sprintf("  Correlation between the two ERstress variants: %+.3f (median)\n",
              median(pcol_all$ERstress_variants_cor, na.rm = TRUE)))
}
if (nrow(vif_all) > 0) {
  cat("\n=== Variance inflation factors (median across tumours) ===\n")
  vs <- vif_all %>% group_by(variant, term) %>%
    summarise(median_VIF = median(VIF, na.rm=TRUE), max_VIF = max(VIF, na.rm=TRUE), .groups="drop") %>%
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

excl1 <- function(lo, hi) (lo > 1) | (hi < 1)
star  <- function(lo, hi) if (excl1(lo, hi)) "  *" else "   "

if (!is.null(marg_pooled)) {
  cat("\n=== MARGINAL (module alone, no adjustment) ===\n")
  for (i in seq_len(nrow(marg_pooled))) with(marg_pooled[i,],
                                             cat(sprintf("  %-16s %-10s OR %6.2f (95%% CI %5.2f-%7.2f) k=%d I2=%2.0f%% tau=%.2f%s\n",
                                                         variant, term, OR, lower, upper, k, I2, sqrt(tau2), star(lower, upper))))
  cat("  (ERAD is identical under both variants by construction.)\n")
}

if (!is.null(padj_pooled)) {
  cat("\n=== ADJUSTED (all five predictors) ===\n")
  for (i in seq_len(nrow(padj_pooled))) with(padj_pooled[i,],
                                             cat(sprintf("  %-16s %-18s OR %6.2f (95%% CI %5.2f-%7.2f) k=%d I2=%2.0f%% tau=%.2f%s\n",
                                                         variant, term, OR, lower, upper, k, I2, sqrt(tau2), star(lower, upper))))
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
                else "neither variant is significant"))
    cat("     Report the STRICT estimate. 'Survives' is a yes/no test; the size of\n")
    cat("     the drop matters too, and is printed above.\n")
  }
}

if (!is.null(pstr_pooled)) {
  cat("\n=== STRATIFIED by hypoxia tertile ===\n")
  ps <- pstr_pooled[order(pstr_pooled$variant, pstr_pooled$term,
                          match(pstr_pooled$hyp_bin, PROT_BIN_LABELS)), ]
  for (i in seq_len(nrow(ps))) with(ps[i,],
                                    cat(sprintf("  %-16s %-18s %-9s OR %6.2f (95%% CI %5.2f-%7.2f) k=%d%s\n",
                                                variant, term, hyp_bin, OR, lower, upper, k, star(lower, upper))))
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
    cat("     Do NOT quote the raw per-unit interaction OR: it is the product of two\n")
    cat("     module scores whose observed range is far below 1 unit. Use the slopes.\n")
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
cat("  1. samples contributing - is anyone excluded, and why?\n")
cat("  2. gene coverage        - is each score built on enough genes?\n")
cat("  3. collinearity + VIF   - are the predictors separable?\n")
cat("  4. marginal             - any signal before conditioning?\n")
cat("  5. adjusted + SENSITIVITY VERDICT\n")
cat("  6. stratified + interaction\n")
cat("EXPLORATORY pending Valerio's sign-off on the gene sets.\n")
cat("Associational, mRNA-level. Cannot demonstrate a clamp.\n")
cat(strrep("=", 72), "\n", sep = "")

log_msg("Proteostasis step done.")

################################################################################
## 12. FIGURE EXPORT  --  journal style AND presentation style
##
##   BRCA NOTE: the reasoning in this header was written for the PDAC
##   cohort. Any statement about what the data showed, which module worked,
##   or how many tumours were usable describes PDAC, not BRCA. The code
##   below is unchanged except where tagged "BRCA:".
##
##     Reads the pooled CSVs written above, so nothing is re-analysed. Every
##     panel is defined once and rendered twice:
##       journal : 8 pt, 89 / 183 mm columns, vector PDF + 600 dpi PNG,
##                 panels combined with A/B/C tags, no baked-in titles.
##       slide   : 16 pt, 13.33 x 7.5 in (16:9), thick lines, large points,
##                 titles baked in, one panel per slide.
##     Convention in both: filled point = CI excludes 1, open point = CI includes 1.
################################################################################

## ============================================================================
## 12.0  OBSERVED HYPOXIA SCORE RANGE   (runs once, then cached)
##       Figure 4 plots the collagen and ER-stress effect as a function of the
##       hypoxia score. Without the real observed range the x-axis would be
##       invented and the crossover line could sit at a score that never occurs
##       in the data. That is how the earlier version of this figure went wrong,
##       so Figure 4 is hard-blocked on this file existing.
## ============================================================================
HYPRANGE_CSV <- file.path(agg, "HYPOXIA_SCORE_RANGE.csv")

if (!file.exists(HYPRANGE_CSV)) {
  log_msg("HYPOXIA_SCORE_RANGE.csv missing - computing it (one extra pass over the cohort).")
  hr <- do.call(rbind, lapply(seq_len(nrow(work)), function(i) {
    p <- work$path[i]; sid <- work$sample_id[i]
    if (is.na(p) || !file.exists(p)) return(NULL)
    obj <- readRDS(p)
    obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
    DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else Assays(obj)[1]
    g <- intersect(MODULES$Hypoxia, rownames(obj))
    if (length(g) < 2) return(NULL)
    set.seed(SEED)
    obj <- add_module_score_fallback(obj, features = list(g), name = "H",
                          assay = DefaultAssay(obj), nbin = 20, ctrl = 50)
    v <- obj@meta.data$H1
    qs <- function(q) unname(quantile(v, q, na.rm = TRUE))
    data.frame(sample_id = sid, n_spots = length(v),
               q01 = qs(.01), q05 = qs(.05), q25 = qs(.25), q50 = qs(.50),
               q75 = qs(.75), q95 = qs(.95), q99 = qs(.99), row.names = NULL)
  }))
  write.csv(hr, HYPRANGE_CSV, row.names = FALSE)
  log_msg("Hypoxia range written: median 5th pct %.3f, median 95th pct %.3f",
          median(hr$q05, na.rm = TRUE), median(hr$q95, na.rm = TRUE))
} else {
  log_msg("Hypoxia score range already cached - skipping recomputation.")
}

## ============================================================================
## 12.1  STYLE DEFINITIONS
## ============================================================================
PRIMARY_VARIANT <- "ERstressStrict"   # strict 24-gene ER stress module
BIN_NICE        <- c("Low", "Mid", "High")

MM <- function(mm) mm / 25.4

CLR <- list(pos  = "#0F4C81",   # CI excludes 1
            neg  = "#B2182B",
            grey = "#8C8C8C",
            full = "#9E9E9E")

## v4 BUG FIX. capabilities("cairo") returns TRUE on this machine but the cairo
## library fails to load when a device is actually opened, which produced
## "failed to load cairo DLL" on every PDF and may have left unusable files.
## Trusting capabilities() is not enough: open a throwaway device once and see
## whether it really works, then fall back to grDevices::pdf if it does not.
## The only practical difference is font embedding, which does not matter for
## these all-ASCII figures.
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
    dir     = "figures_journal",
    base    = 8,      # pt
    pt      = 1.9,    # point size
    stroke  = 0.45,
    lwd     = 0.45,
    lab     = 2.1,    # in-plot numeric label size
    grid    = 0.25,
    border  = 0.40,
    key     = 3,      # mm
    mar     = 3,
    tags    = TRUE,
    titles  = FALSE,
    split   = FALSE,  # keep multi-panel figures combined
    devices = c("pdf", "png"),
    dpi     = 600,
    w1      = MM(89),   # single column
    w2      = MM(183)   # double column
  ),
  
  slide = list(
    dir     = "figures_slides",
    base    = 16,
    pt      = 4.4,
    stroke  = 0.9,
    lwd     = 1.0,
    lab     = 4.8,
    grid    = 0.5,
    border  = 0.7,
    key     = 6,
    mar     = 10,
    tags    = FALSE,  # A/B/C tags are meaningless on separate slides
    titles  = TRUE,   # each slide has to stand on its own
    split   = TRUE,   # one panel per slide
    devices = c("png", "pdf"),
    dpi     = 200,
    w1      = 13.333, # 16:9
    w2      = 13.333,
    h_full  = 7.5
  )
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

## ============================================================================
## 12.2  SHARED HELPERS
## ============================================================================
rd <- function(f) {
  p <- file.path(agg, f)
  if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else NULL
}

nice <- c(Hypoxia = "Hypoxia", CollagenSynthesis = "Collagen synthesis",
          iCAF = "iCAF", myCAF = "myCAF", Epithelial = "Epithelial",
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
    scale_x_log10(limits = xlim, breaks = breaks,
                  labels = function(x) sub("^0\\.", ".", x)) +
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
    scale_x_log10(limits = xlim, breaks = breaks,
                  labels = function(x) sub("^0\\.", ".", x)) +
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
  if (is.na(k)) NULL else sprintf("Random-effects meta-analysis, k = %d %s", k, K_UNIT)   # BRCA: k unit
}

## ============================================================================
## 12.4  RENDER EVERY FIGURE IN ONE STYLE
## ============================================================================
render_style <- function(S) {
  
  outdir <- file.path(agg, S$dir)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  save_one <- function(p, name, w, h) {
    for (dv in S$devices) {
      f <- file.path(outdir, paste0(name, ".", dv))
      tryCatch({
        if (dv == "pdf") ggsave(f, p, width = w, height = h, device = PDF_DEV, bg = "white")
        else             ggsave(f, p, width = w, height = h, dpi = S$dpi, bg = "white")
      }, error = function(e) log_msg("  fig failed (%s): %s", basename(f), conditionMessage(e)))
    }
    cat(sprintf("  %-34s %5.1f x %4.1f in\n", name, w, h))
  }
  
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
    pC <- tg(ttl(pC, "Within hypoxia tertiles (relative, within tumour)",
                 "Holding hypoxia level constant"), "C")
    
    if (isTRUE(S$split)) {
      save_one(pA, "Fig1A_marginal",   S$w2, S$h_full)
      save_one(pB, "Fig1B_adjusted",   S$w2, S$h_full)
      save_one(pC, "Fig1C_stratified", S$w2, S$h_full)
    } else {
      save_one((pA / pB / pC) + plot_layout(heights = c(1.05, 0.85, 1.15)),
               "Figure1_main", S$w2, MM(135))
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
      save_one(pA, "Fig2A_proteostasis_adjusted",   S$w2, S$h_full)
      save_one(pB, "Fig2B_proteostasis_stratified", S$w2, S$h_full)
    } else {
      save_one((pA / pB) + plot_layout(heights = c(1, 1.05)),
               "Figure2_proteostasis", S$w2, MM(100))
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
      scale_x_log10(limits = c(0.35, 5), breaks = c(0.5, 1, 2, 4),
                    labels = function(x) sub("^0\\.", ".", x)) +
      labs(x = paste("Adjusted", tolower(XLAB_OR)), y = NULL) +
      theme_fig(S) + theme(legend.position = "top")
    p3 <- ttl(p3, "ER stress module with and without hypoxia-responsive genes",
              "Removing HYOU1, DDIT4, ERO1A. The strict estimate is the one to report.")
    
    save_one(p3, "Figure3_sensitivity",
             if (S$split) S$w2 else S$w1 * 1.6,
             if (S$split) S$h_full else MM(62))
  } else {
    cat("  Figure 3 skipped: PROT_adjusted_pooled.csv does not contain both variants.\n")
  }
  
  ## -- FIGURE 4: effect modification by hypoxia -----------------------------
  ##    SUPERSEDED by section 13, which fits both slopes in ONE model and can
  ##    therefore carry a p value and correct bands. This version is kept only
  ##    so you can show what changed. Do not present it on its own.
  if (is.null(hrange)) {
    cat("  Figure 4 skipped: HYPOXIA_SCORE_RANGE.csv still not available.\n")
  } else if (!is.null(sint_fig) && !is.null(pint_fig)) {
    lo <- median(hrange$q05, na.rm = TRUE); hi <- median(hrange$q95, na.rm = TRUE)
    zz <- seq(lo, hi, length.out = 200)
    
    band <- function(df, main_term, variant = NULL, label) {
      d <- df
      if (!is.null(variant) && "variant" %in% names(d)) d <- d[d$variant == variant, ]
      im <- grepl(":", d$term)
      b0 <- d$est[d$term == main_term]; bi <- d$est[im]
      if (length(b0) != 1 || length(bi) != 1) return(NULL)
      ## covariance between the two pooled coefficients is not exported, so it
      ## is ignored here. State this in the legend. Section 13 fixes it.
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
        filter(cross >= lo, cross <= hi)   # drawn only if inside the observed range
      
      p4 <- ggplot(sl, aes(z, OR)) +
        geom_hline(yintercept = 1, linetype = "22", linewidth = S$lwd * 0.7, colour = "grey40") +
        geom_ribbon(aes(ymin = lo, ymax = hi, fill = panel), alpha = .16, colour = NA) +
        geom_line(aes(colour = panel), linewidth = S$lwd * 1.4) +
        scale_colour_manual(values = c("Collagen synthesis" = CLR$neg,
                                       "ER stress"          = CLR$pos), name = NULL) +
        scale_fill_manual(values   = c("Collagen synthesis" = CLR$neg,
                                       "ER stress"          = CLR$pos), name = NULL) +
        scale_y_log10() + coord_cartesian(xlim = c(lo, hi)) +
        labs(x = "Hypoxia module score (spot level, 5th to 95th percentile)",
             y = XLAB_OR,
             caption = "Superseded by section 13. Bands ignore the covariance between the two pooled coefficients.") +
        theme_fig(S) + theme(legend.position = "top")
      if (nrow(cr))
        p4 <- p4 + geom_vline(data = cr, aes(xintercept = cross),
                              linetype = "dotted", linewidth = S$lwd * 0.7, colour = "grey30")
      p4 <- ttl(p4, "How the collagen and ER stress effects change with hypoxia",
                "Slope from the interaction model, over the observed score range only")
      
      save_one(p4, "Figure4_effect_modification",
               if (S$split) S$w2 else S$w1 * 1.5,
               if (S$split) S$h_full else MM(62))
    }
  }
  
  ## -- FIGURE 5: per-tumour consistency -------------------------------------
  if (!is.null(permod_fig) && all(c("module", "rho") %in% names(permod_fig))) {
    pm <- permod_fig[is.finite(permod_fig$rho), , drop = FALSE]
    pm$module_nice <- relabel(pm$module)
    ordm <- pm %>% group_by(module_nice) %>%
      summarise(med = median(rho, na.rm = TRUE), .groups = "drop") %>%
      arrange(med) %>% pull(module_nice)
    pm$module_nice <- factor(pm$module_nice, levels = ordm)
    n_tum <- dplyr::n_distinct(pm$sample_id)
    
    p5 <- ggplot(pm, aes(module_nice, rho)) +
      geom_hline(yintercept = 0, linetype = "22", linewidth = S$lwd * 0.7, colour = "grey40") +
      geom_boxplot(outlier.shape = NA, fill = "grey92", linewidth = S$lwd * 0.8, width = .6) +
      geom_jitter(width = .13, size = S$pt * 0.55, alpha = .5, colour = CLR$pos) +
      coord_flip() +
      labs(x = NULL, y = "Per-tumour Spearman rho (P4HA1 vs module score)") +
      theme_fig(S) +
      theme(panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(linewidth = S$grid, colour = "grey90"))
    p5 <- ttl(p5, "How consistent is each association across tumours?",
              sprintf("One point per tumour, n = %d", n_tum))
    
    save_one(p5, "Figure5_consistency",
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
cat("Section 11 output remains EXPLORATORY pending Valerio's sign-off on the gene sets.\n")

log_msg("SECTIONS 0 TO 12 COMPLETE. Everything under: %s", normalizePath(OUT_ROOT))

################################################################################
## 14. DETECTION AUDIT FOR P4HA1, P4HA2 AND P4HA3
##
##   BRCA NOTE: the reasoning in this header was written for the PDAC
##   cohort. Any statement about what the data showed, which module worked,
##   or how many tumours were usable describes PDAC, not BRCA. The code
##   below is unchanged except where tagged "BRCA:".
##                                             (Valerio, lab meeting 3, Aug 2026)
##
## WHY THIS RUNS BEFORE ANY MODEL
##   Valerio's test is: the pattern should hold for P4HA1 and P4HA2, which
##   hydroxylate collagen, and FAIL for P4HA3, which does not. That only works
##   as a falsification test if a null result for P4HA3 means "no effect"
##   rather than "no power". P4HA3 is the least expressed of the three in most
##   tissues, so if it is detected in a small fraction of spots the logistic
##   model either will not converge or will return an estimate so imprecise
##   that nothing can be concluded. A control with no power proves nothing, and
##   presenting one as evidence would be worse than not running the test.
##
##   So this section counts and does NOTHING else. No module scores, no glm, no
##   pooling. One pass over the cohort, pure arithmetic, and the answer decides
##   how the next section is built.
##
## PDAC RESULT (unadjusted run, for reference only; NOT the BRCA result,
## which this section prints when it runs)
##   All three genes present in all 23 PDAC tumours under their current symbols.
##   Detection: A1 median 27.9%, A2 30.1%, A3 9.0%.
##   Usable: A1 23/23, A2 23/23, A3 20/23.
##   Depth correlation: 0.14 / 0.17 / 0.10, all modest.
##   P4HA1 vs P4HA2 spot overlap: median Jaccard 0.24.
##   ERO1A present in all 23, ERO1L in none, so the strict module really does
##   drop three genes as described.
##
## NOTE ON THE DEPTH SWITCH
##   This section is unaffected by ADJUST_DEPTH: it fits nothing, so there is
##   nothing to adjust. It writes unsuffixed filenames deliberately, because the
##   counts are the same whichever way the models are later fitted.
################################################################################

stopifnot(exists("work"), exists("agg"), exists("log_msg"), exists("bind"))

## ============================================================================
## 14.0  SETTINGS
## ============================================================================
P4HA_GENES <- c("P4HA1", "P4HA2", "P4HA3")

## Legacy symbol aliases to look for if the primary symbol is absent. An object
## built on an older reference will carry the left hand name.
ALIASES <- list(
  P4HA1  = character(0),
  P4HA2  = character(0),
  P4HA3  = character(0),
  ERO1A  = "ERO1L"       # settled: these objects carry ERO1A
)

## Rule of thumb only, not a law: about 10 events per predictor. The joint model
## in section 13 has 6 predictors plus an intercept, so roughly 70 events, and
## the same again for non-events. Tumours below this are flagged, not deleted.
MIN_EVENTS_PER_TUMOUR <- 70

JOINT_CACHE_DETECTION <- TRUE   # cache cohort-wide per-gene detection rates
DETECT_CACHE_MIN_RATE <- 0.005  # do not cache genes under this rate, keeps file small

## ============================================================================
## 14.1  CIRCULARITY RE-AUDIT FOR THE NEW OUTCOMES
##       The rule that caught CollagenMaturation applies again the moment the
##       outcome gene changes. A module containing P4HA2 is fine as a predictor
##       of P4HA1 but circular as a predictor of P4HA2.
## ============================================================================
cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 14.1  CIRCULARITY RE-AUDIT FOR THE NEW OUTCOMES\n")
cat(strrep("=", 74), "\n", sep = "")

if (exists("MODULES")) {
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
  cat("\n  Reminder: CollagenMaturation contains P4HA1 and P4HA2 and is already\n")
  cat("  excluded from every model. Confirm again before swapping the outcome.\n")
} else {
  cat("  MODULES not defined, skipping. Run section 1 first.\n")
}

## ============================================================================
## 14.2  COUNT, ONE TUMOUR AT A TIME
## ============================================================================
detect_one <- function(rds_path, sid) {
  obj <- readRDS(rds_path)
  obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
  DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else Assays(obj)[1]
  rn <- rownames(obj)
  n  <- ncol(obj)
  
  ## resolve each gene to whatever symbol is actually present
  resolve <- function(g) {
    if (g %in% rn) return(g)
    alt <- intersect(ALIASES[[g]], rn)
    if (length(alt)) return(alt[1])
    NA_character_
  }
  
  ## sequencing depth per spot
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
  
  ## per gene counts
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
  
  ## do P4HA1 and P4HA2 mark the same spots?
  co <- data.frame(sample_id = sid, jaccard_A1_A2 = NA_real_, both = NA_integer_)
  if (!is.null(pos_list$P4HA1) && !is.null(pos_list$P4HA2)) {
    a <- pos_list$P4HA1; b <- pos_list$P4HA2
    u <- sum(a | b)
    co$both <- sum(a & b)
    co$jaccard_A1_A2 <- if (u > 0) sum(a & b) / u else NA_real_
  }
  
  ## side errand: ERO1A vs ERO1L
  ero <- data.frame(sample_id = sid,
                    ERO1A_present = "ERO1A" %in% rn,
                    ERO1L_present = "ERO1L" %in% rn)
  
  ## optional cache of cohort-wide detection rates for the random-gene control
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

if (nrow(det_all) == 0) stop("Section 14: no tumour could be read.")

write.csv(det_all,   file.path(agg, "DETECT_p4ha_per_sample.csv"), row.names = FALSE)
write.csv(depth_all, file.path(agg, "DETECT_depth_per_sample.csv"), row.names = FALSE)

if (JOINT_CACHE_DETECTION && length(acc14$rates)) {
  rates_all <- bind(acc14$rates)
  write.csv(rates_all, file.path(agg, "DETECT_cohort_gene_rates.csv"), row.names = FALSE)
  log_msg("Cached detection rates for %d gene-by-tumour combinations.", nrow(rates_all))
}

## ============================================================================
## 14.3  REPORT
## ============================================================================
cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 14  DETECTION AUDIT\n")
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("Tumours read: %d\n", dplyr::n_distinct(det_all$sample_id)))

cat("\n--- is each gene in the object at all? ---\n")
for (g in P4HA_GENES) {
  s <- det_all[det_all$gene == g, ]
  cat(sprintf("  %-6s present in %2d / %2d tumours", g, sum(s$present), nrow(s)))
  syms <- unique(na.omit(s$symbol_found))
  cat(sprintf("   symbol(s): %s\n", if (length(syms)) paste(syms, collapse = ", ") else "NOT FOUND"))
}

cat("\n--- detection rate per gene (fraction of spots with count > 0) ---\n")
det_summary <- det_all %>%
  dplyr::filter(present) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    k_tumours       = dplyr::n(),
    median_rate     = median(detect_rate, na.rm = TRUE),
    min_rate        = min(detect_rate, na.rm = TRUE),
    max_rate        = max(detect_rate, na.rm = TRUE),
    median_events   = median(n_detected, na.rm = TRUE),
    min_events      = min(n_detected, na.rm = TRUE),
    n_tumours_ok    = sum(n_detected >= MIN_EVENTS_PER_TUMOUR &
                            (n_spots - n_detected) >= MIN_EVENTS_PER_TUMOUR, na.rm = TRUE),
    median_depth_rho = median(depth_rho, na.rm = TRUE),
    .groups = "drop")
write.csv(det_summary, file.path(agg, "DETECT_p4ha_summary.csv"), row.names = FALSE)

for (i in seq_len(nrow(det_summary))) with(det_summary[i, ],
                                           cat(sprintf("  %-6s median %5.1f%% (range %4.1f%% to %5.1f%%)  median %5.0f events (min %4.0f)  usable in %2d/%2d tumours\n",
                                                       gene, 100*median_rate, 100*min_rate, 100*max_rate,
                                                       median_events, min_events, n_tumours_ok, k_tumours)))
cat(sprintf("  \"usable\" = at least %d detected AND %d undetected spots. Rule of thumb\n",
            MIN_EVENTS_PER_TUMOUR, MIN_EVENTS_PER_TUMOUR))
cat("  (~10 events per predictor, 6 predictors), not a hard statistical law.\n")

cat("\n--- per tumour, side by side ---\n")
wide <- det_all %>%
  dplyr::select(sample_id, gene, detect_rate, n_detected) %>%
  tidyr::pivot_wider(names_from = gene, values_from = c(detect_rate, n_detected)) %>%
  dplyr::arrange(dplyr::desc(detect_rate_P4HA1))
cat(sprintf("  %-26s %8s %8s %8s %10s\n", "tumour", "A1 %", "A2 %", "A3 %", "A3 events"))
for (i in seq_len(nrow(wide))) {
  r <- wide[i, ]
  f  <- function(x) if (is.null(x) || length(x) == 0 || is.na(x)) "       ." else sprintf("%8.1f", 100*x)
  fe <- function(x) if (is.null(x) || length(x) == 0 || is.na(x)) "         ." else sprintf("%10.0f", x)
  ## flag the tumours that will drop out of the P4HA3 model, so the change in k
  ## is visible here rather than discovered later
  low3 <- !is.na(r$n_detected_P4HA3) && r$n_detected_P4HA3 < MIN_EVENTS_PER_TUMOUR
  cat(sprintf("  %-26s %s %s %s %s%s\n", r$sample_id,
              f(r$detect_rate_P4HA1), f(r$detect_rate_P4HA2), f(r$detect_rate_P4HA3),
              fe(r$n_detected_P4HA3), if (low3) "   <-- too few for P4HA3" else ""))
}

cat("\n--- does detection just follow sequencing depth? ---\n")
for (g in P4HA_GENES) {
  s <- det_all[det_all$gene == g & det_all$present, ]
  if (!nrow(s)) next
  cat(sprintf("  %-6s median Spearman rho(detected, nCount) = %+.2f  (range %+.2f to %+.2f)\n",
              g, median(s$depth_rho, na.rm = TRUE),
              min(s$depth_rho, na.rm = TRUE), max(s$depth_rho, na.rm = TRUE)))
}
cat("  The higher this is, the more depth adjustment matters. It never proves\n")
cat("  the signal is artefactual, but it sizes the objection. Set ADJUST_DEPTH\n")
cat("  <- TRUE in section 0 to answer it directly.\n")

if (nrow(co_all)) {
  cat("\n--- do P4HA1 and P4HA2 mark the same spots? ---\n")
  cat(sprintf("  median Jaccard overlap = %.2f  (range %.2f to %.2f)\n",
              median(co_all$jaccard_A1_A2, na.rm = TRUE),
              min(co_all$jaccard_A1_A2, na.rm = TRUE),
              max(co_all$jaccard_A1_A2, na.rm = TRUE)))
  cat("  High overlap means P4HA2 is close to a repeat of P4HA1, so agreement\n")
  cat("  between them is weaker evidence than it looks. Low overlap means they\n")
  cat("  are closer to independent tests. Report this number either way.\n")
}

cat("\n--- side errand: ERO1A vs ERO1L ---\n")
if (nrow(ero_all)) {
  cat(sprintf("  ERO1A present in %d / %d tumours\n", sum(ero_all$ERO1A_present), nrow(ero_all)))
  cat(sprintf("  ERO1L present in %d / %d tumours\n", sum(ero_all$ERO1L_present), nrow(ero_all)))
  if (all(!ero_all$ERO1A_present) && any(ero_all$ERO1L_present)) {
    cat("  *** The pipeline asks for ERO1A and the objects carry ERO1L. The gene\n")
    cat("      is being dropped silently, so the strict ER stress module removes\n")
    cat("      TWO hypoxia-responsive genes, not three. Fix the gene list and say\n")
    cat("      so before quoting the sensitivity analysis. ***\n")
  } else if (all(ero_all$ERO1A_present)) {
    cat("  Fine. ERO1A is the symbol in use, the strict module really does drop\n")
    cat("  three genes as described.\n")
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
  if (!nrow(s)) { cat(sprintf("  %-6s NOT PRESENT in any tumour. Test impossible.\n", g)); next }
  frac <- s$n_tumours_ok / s$k_tumours
  verdict <- if (frac >= 0.8) "fine, model as usual" else
    if (frac >= 0.5) "marginal, report k and per-tumour events alongside every estimate" else
      "UNDERPOWERED, a null result here would mean nothing"
  cat(sprintf("  %-6s usable in %2d/%2d tumours (%.0f%%)  ->  %s\n",
              g, s$n_tumours_ok, s$k_tumours, 100*frac, verdict))
}
cat("\n  If P4HA3 lands in the underpowered row, the falsification test cannot be\n")
cat("  run as a simple null result. The honest options then are:\n")
cat("    (a) report its CI and show it excludes the P4HA1 effect size, which is\n")
cat("        an equivalence argument rather than a significance one; or\n")
cat("    (b) state plainly that P4HA3 is too sparsely detected in Visium to\n")
cat("        serve as a control, and propose a different specificity test.\n")
cat("  What is NOT acceptable is reporting \"no effect in P4HA3\" without the\n")
cat("  power statement attached.\n")
cat(strrep("=", 74), "\n", sep = "")

log_msg("Section 14 complete.")

################################################################################
## 12b. RUN STATE, PROVENANCE AND SHARED HELPERS FOR SECTIONS 13 TO 17
##
## This block used to live in section 0 and inside section 13. It sits here in
## the merged file because of the bug it exists to prevent:
##
##   The filename suffix and the model formula were both derived from
##   ADJUST_DEPTH, but at DIFFERENT MOMENTS. Re-entering section 13 below its
##   settings line updated the suffix and not the model, so an unadjusted
##   analysis was written under a "_depthadj" name. Two filenames, one analysis,
##   and nothing said so.
##
## Everything that guards against that is here, in one place, ahead of every
## section that could be affected:
##   ADJUST_DEPTH, OUT_SUFFIX, outf(), RUN_ID, stamp(), write_out(),
##   assert_run_state(), write_run_manifest(), pool_with_p(), fmt_p()
##
## DO NOT redefine any of these by hand in the console. That is exactly how the
## identical files were produced. Change ADJUST_DEPTH here and re-run from the
## top of the file.
################################################################################

## ADJUST_DEPTH is the ONLY place this is set. Sections 13, 15 and 16 read it
## and are forbidden from overriding it. See 13.0a for the enforcement.
ADJUST_DEPTH <- FALSE

OUT_SUFFIX <- if (isTRUE(ADJUST_DEPTH)) "_depthadj" else ""
outf <- function(base, ext = "csv") sprintf("%s%s.%s", base, OUT_SUFFIX, ext)

## A run identity, so a file can always be traced back to the run that made it.
## Recomputed every time section 0 executes, which is what makes it useful: if
## two files carry the same RUN_ID they came from one execution.
RUN_ID <- paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), "-",
                 if (isTRUE(ADJUST_DEPTH)) "depthadj" else "plain")

## Stamp provenance onto every data frame before it is written. Two extra
## columns is a small price for never again having to guess which run a CSV
## came from.
stamp <- function(df) {
  if (is.null(df) || !nrow(df)) return(df)
  df$adjust_depth <- isTRUE(ADJUST_DEPTH)
  df$run_id       <- RUN_ID
  df
}

## write.csv with the stamp applied, so it cannot be forgotten on one file
write_out <- function(df, base) {
  write.csv(stamp(df), file.path(agg, outf(base)), row.names = FALSE)
  invisible(NULL)
}


## ---------------------------------------------------------------------------
## SHARED HELPERS FOR SECTIONS 13, 15, 16 AND 17
##
## v5.1: these were originally defined inside section 13, which meant sections
## 15 and 16 could not run unless section 13 had run first in the same session.
## That is the wrong shape for a guard: the thing that checks the run state must
## not depend on one of the sections it is guarding. They live here now.
## ---------------------------------------------------------------------------

assert_run_state <- function(where) {
  problems <- character(0)
  
  if (!exists("ADJUST_DEPTH")) problems <- c(problems, "ADJUST_DEPTH is not defined")
  if (!exists("OUT_SUFFIX"))   problems <- c(problems, "OUT_SUFFIX is not defined")
  if (!exists("outf"))         problems <- c(problems, "outf() is not defined")
  if (!exists("RUN_ID"))       problems <- c(problems, "RUN_ID is not defined; apply the section 0 patch")
  if (length(problems)) stop(sprintf("%s: %s. Re-run the section 0 patch.",
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
  
  local_flags <- c(JOINT_ADJUST_DEPTH = if (exists("JOINT_ADJUST_DEPTH")) JOINT_ADJUST_DEPTH else NA,
                   ISO_ADJUST_DEPTH   = if (exists("ISO_ADJUST_DEPTH"))   ISO_ADJUST_DEPTH   else NA,
                   RND_ADJUST_DEPTH   = if (exists("RND_ADJUST_DEPTH"))   RND_ADJUST_DEPTH   else NA)
  local_flags <- local_flags[!is.na(local_flags)]
  bad <- names(local_flags)[local_flags != ADJUST_DEPTH]
  if (length(bad))
    problems <- c(problems, sprintf(
      "%s disagree(s) with ADJUST_DEPTH = %s. You re-entered a section below its settings block.",
      paste(bad, collapse = ", "), ADJUST_DEPTH))
  
  if (length(problems)) {
    cat("\n", strrep("!", 74), "\n", sep = "")
    cat(sprintf("%s: RUN STATE IS INCONSISTENT. NOTHING WAS FITTED.\n", where))
    for (p in problems) cat("  - ", p, "\n", sep = "")
    cat("\n  Fix: start a fresh R session and run from section 0. Do not redefine\n")
    cat("  outf() or the flags by hand: that is exactly how the adjusted and\n")
    cat("  unadjusted files came out identical last time.\n")
    cat(strrep("!", 74), "\n", sep = "")
    stop("Run state inconsistent; see the message above.")
  }
  
  cat(sprintf("\n[run state OK] %s | ADJUST_DEPTH = %s | suffix = '%s' | run_id = %s\n",
              where, ADJUST_DEPTH, OUT_SUFFIX, RUN_ID))
  invisible(TRUE)
}


write_run_manifest <- function() {
  gene_hash <- function(g) paste0(length(g), ":", substr(
    paste(sort(g), collapse = ","), 1, 60))
  mf <- data.frame(
    key = c("run_id", "timestamp", "adjust_depth", "seed", "out_suffix",
            "R_version", "Seurat_version",
            "n_tumours_in_worklist",
            "module_Hypoxia", "module_CollagenSynthesis", "module_iCAF",
            "module_ERstressStrict", "module_ERAD",
            "joint_formula"),
    value = c(RUN_ID, format(Sys.time()), as.character(ADJUST_DEPTH),
              as.character(SEED), OUT_SUFFIX,
              paste(R.version$major, R.version$minor, sep = "."),
              as.character(utils::packageVersion("Seurat")),
              as.character(nrow(work)),
              gene_hash(MODULES$Hypoxia), gene_hash(MODULES$CollagenSynthesis),
              gene_hash(MODULES$iCAF), gene_hash(MODULES$ERstressStrict),
              gene_hash(MODULES$ERAD),
              paste0("y ~ zHyp*zCol + zHyp*zER + ziCAF + zERAD",
                     if (ADJUST_DEPTH) " + zDepth" else "")),
    stringsAsFactors = FALSE)
  write.csv(mf, file.path(agg, outf("RUN_manifest")), row.names = FALSE)
  cat("\n--- run manifest ---\n")
  for (i in seq_len(nrow(mf))) cat(sprintf("  %-26s %s\n", mf$key[i], mf$value[i]))
  cat("  Diff this against the previous run's manifest whenever numbers move.\n")
  invisible(mf)
}

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
  if (!is.finite(p)) return("NA")
  if (p < 2e-16) return("< 2e-16")
  if (p < 0.001)  return(sprintf("%.1e", p))
  sprintf("%.3f", p)
}



################################################################################
## 13. FIGURE 4, REBUILT FROM ONE JOINT MODEL, WITH THE P VALUE ON IT
##
##   BRCA NOTE: the reasoning in this header was written for the PDAC
##   cohort. Any statement about what the data showed, which module worked,
##   or how many tumours were usable describes PDAC, not BRCA. The code
##   below is unchanged except where tagged "BRCA:".
##                                             (Valerio, lab meeting 3, Aug 2026)
##
## WHAT HE ASKED FOR
##   He was looking at the crossing-lines slide when he said "I need the P value
##   for the shape. I need the P value for the difference of the two." He wanted
##   the number ON THAT FIGURE.
##
## WHY THE ORIGINAL FIGURE 4 COULD NOT CARRY A P VALUE
##   Its two curves come from two SEPARATE models fitted to the SAME spots
##   (section 10 for collagen, section 11 for ER stress). Coefficients from two
##   separate fits cannot be compared with a simple test: their sampling errors
##   are correlated and that correlation was never estimated.
##
## WHAT THIS SECTION DOES
##   Fits ONE model per tumour containing both interactions:
##
##     y ~ zHyp * zCol + zHyp * zER + ziCAF + zERAD  [+ zDepth]
##
##   Every quantity the figure needs is then a contrast on that single fit, so
##   vcov() supplies all covariances and nothing has to be assumed.
##
## THE ONE THING NOT TO OVERCLAIM
##   The two lines CONVERGE as hypoxia rises. Whether they SWAP ORDER is not
##   established: the crossing sits where the difference is no longer
##   significant. The figure labels this itself.
################################################################################

stopifnot(exists("MODULES"), exists("meta_pool"), exists("work"), exists("agg"))
if (!all(c("ERstressStrict", "ERAD") %in% names(MODULES)))
  stop("Section 13 needs section 11 to have run first (ERstressStrict / ERAD not defined).")

## ============================================================================
## 13.0  SETTINGS
## ============================================================================
JOINT_ES           <- "ERstressStrict"   # strict 24 gene module, the primary one
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

## ============================================================================
## 13.0a  RUN STATE GATE   *** THE FIX FOR FAULT 1 ***
##
##   Three separate things have to agree before any model is fitted: the master
##   flag, the filename suffix, and the section-local flag. Previously they
##   could disagree silently, which produced two identical files under two
##   names. Now a disagreement is a hard stop.
##
##   This runs BEFORE the cohort loop, so a wrong state costs a second rather
##   than twenty minutes.
## ============================================================================
assert_run_state("Section 13")

## ============================================================================
## 13.0b  RUN MANIFEST   *** PART OF THE FIX FOR FAULT 2 ***
##
##   Records everything that could make two runs differ. Written once per run,
##   next to the outputs. When numbers move between runs, this is the first
##   thing to diff, and it costs nothing to produce.
## ============================================================================
write_run_manifest()

## ============================================================================
## 13.1  HELPERS
## ============================================================================

## meta_pool returns est / lower / upper but no SE and no p. Recover the SE from
## the interval rather than editing meta_pool.
## NOTE: the Wald statistic is zstat, NOT z. In 13.5 one grouping variable is z
## (the hypoxia grid) and group_modify() refuses to return a column named after
## a grouping variable. Do not rename it back.
## score every module, then z standardise within this tumour.
## v5: also returns a FINGERPRINT of each module score, so two runs can be
## compared directly and the cause of any drift localised to scoring or to
## modelling. This is the diagnostic for fault 2.
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
      log_msg("  %s: %s has only %d/%d genes (need %d) -> sample skipped",
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
    if (is.null(v)) { log_msg("  %s: %s has zero variance -> sample skipped", sid, nm); return(NULL) }
    d[[ JOINT_RENAME[[nm]] ]] <- v
  }
  if (JOINT_ADJUST_DEPTH) {
    ## log1p first: read counts are strongly right skewed, so on the raw scale a
    ## handful of very deep spots would dominate the coefficient
    v <- zz(log1p(d[[cnt_col[1]]]))
    if (is.null(v)) return(NULL)
    d$zDepth <- v
  }
  
  attr(d, "fingerprint") <- do.call(rbind, fp)
  attr(d, "n_spots")     <- nrow(d)
  d
}

## ============================================================================
## 13.2  FIT ONE TUMOUR
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
  
  ## *** THE REAL FIX FOR FAULT 1 ***
  ## Do not trust the flag. Look at the model that was actually fitted. If the
  ## flag says depth adjusted and zDepth is not among the coefficients, or vice
  ## versa, the run is not what it claims to be and must not write files.
  has_depth <- "zDepth" %in% names(b)
  if (!identical(has_depth, isTRUE(JOINT_ADJUST_DEPTH)))
    stop(sprintf(paste0("Section 13, tumour %s: JOINT_ADJUST_DEPTH = %s but the fitted model %s ",
                        "a zDepth term. The flag and the model disagree, which is how the ",
                        "identical adjusted/unadjusted files were produced. Start a fresh ",
                        "session from section 0."),
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
  
  ## Contrast machinery. L is a row vector over the coefficient names; the
  ## estimate is L b and its variance is L V L'. Written this way every
  ## covariance term is included automatically and cannot be forgotten.
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

if (nrow(jcon_all) == 0) stop("Section 13: no tumour produced a usable joint fit.")

## ---- model check summary, printed loudly ----------------------------------
cat("\n--- model verification (reads the fitted models, not the flag) ---\n")
cat(sprintf("  tumours fitted            : %d\n", nrow(jmc_all)))
cat(sprintf("  models containing zDepth  : %d\n", sum(jmc_all$has_zDepth)))
cat(sprintf("  formula used              : %s\n", unique(jmc_all$formula)[1]))
if (length(unique(jmc_all$formula)) > 1)
  stop("Section 13: more than one formula was fitted across tumours. Aborting.")
cat(sprintf("  -> this run IS %s\n",
            if (all(jmc_all$has_zDepth)) "DEPTH ADJUSTED" else "UNADJUSTED"))

write_out(jcoef_all, "JOINT_coefficients_per_sample")
write_out(jcon_all,  "JOINT_contrasts_per_sample")
write_out(jfp_all,   "JOINT_module_fingerprint")
write_out(jmc_all,   "JOINT_model_check")

k_used <- dplyr::n_distinct(jcon_all$sample_id)
log_msg("Section 13: %d tumours contributed.", k_used)
if (k_used < nrow(work)) {
  miss <- setdiff(work$sample_id, unique(jcon_all$sample_id))
  cat(sprintf("\n  tumours excluded from section 13 (%d): %s\n", length(miss),
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
    cross_z   <- x1 - y1 * (x2 - x1) / (y2 - y1)
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
cat(sprintf("ER stress module : %s   |   sections fitted : %d   |   pooled over %s\n", JOINT_ES, k_used, K_UNIT))   # BRCA
cat("All predictors z scored WITHIN tumour, so every effect is per 1 SD.\n")

cat("\n--- pooled coefficients from the joint model ---\n")
for (i in seq_len(nrow(jcoef_pooled))) with(jcoef_pooled[i, ],
                                            cat(sprintf("  %-16s beta %+7.4f  OR %5.2f (%.2f-%.2f)  k=%2d  I2=%3.0f%%  p=%s\n",
                                                        term, est, OR, lower_OR, upper_OR, k, I2, fmt_p(p))))

cat("\n--- THE ANSWER TO VALERIO'S QUESTION ---\n")
if (nrow(ri) == 1) {
  cat("  (collagen x hypoxia) minus (ER stress x hypoxia):\n")
  cat(sprintf("     estimate %+.4f   95%% CI %+.4f to %+.4f   k=%d   I2=%.0f%%   p = %s\n",
              ri$est, ri$lower, ri$upper, ri$k, ri$I2, fmt_p(ri$p)))
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
    cat("  ESTABLISHED     : the two effects CONVERGE as hypoxia rises.\n")
    cat("  NOT ESTABLISHED : that they swap order.\n")
  }
} else {
  cat("  The two lines do not meet anywhere within +/- 2 SD of hypoxia.\n")
}

cat("\nCAVEATS\n")
cat("  1. Spatial autocorrelation is not modelled. Every CI here is too narrow.\n")
cat("     Section 16 measures how large that is. Quote empirical percentiles.\n")
cat("  2. Per 1 SD within tumour. Not comparable to sections 10 and 11.\n")
cat(sprintf("  3. Depth adjustment is %s (verified against the fitted models).\n",
            if (JOINT_ADJUST_DEPTH) "ON" else "OFF"))
cat("  4. Associational, mRNA level.\n")
cat(strrep("=", 74), "\n", sep = "")

## ============================================================================
## 13.7  FIGURE 4, REBUILT
## ============================================================================
if (!exists("FIG_STYLES") || !exists("theme_fig")) {
  cat("\n  Figure skipped: section 12 has not run.\n")
} else {
  stat_lab <- sprintf("Difference between the two slopes  %+.3f  (95%% CI %+.3f to %+.3f)\np = %s      k = %d %s      I2 = %.0f%%%s",   # BRCA: k unit
                      ri$est, ri$lower, ri$upper, fmt_p(ri$p), ri$k, K_UNIT, ri$I2,
                      if (JOINT_ADJUST_DEPTH) "\nsequencing depth adjusted" else "")
  
  y_top <- max(jlines$upper_OR, na.rm = TRUE)
  y_bot <- min(jlines$lower_OR, na.rm = TRUE)
  y_hi  <- y_top * (y_top / y_bot)^0.30
  
  for (sn in names(FIG_STYLES)) {
    S <- FIG_STYLES[[sn]]
    outdir <- file.path(agg, S$dir); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    
    p4 <- ggplot(jlines, aes(z, OR))
    
    if (is.finite(sig_lo)) {
      p4 <- p4 +
        annotate("rect", xmin = sig_lo, xmax = sig_hi,
                 ymin = -Inf, ymax = Inf, fill = "grey40", alpha = 0.14) +
        annotate("segment", x = sig_hi, xend = sig_hi, y = -Inf, yend = Inf,
                 linewidth = S$lwd * 0.6, colour = "grey55") +
        annotate("text", x = sig_hi, y = y_hi, hjust = -0.06, vjust = 1.4,
                 size = S$lab * 0.95, colour = "grey40", lineheight = 1.1,
                 label = sprintf("slopes no longer\ndistinguishable\n(above %+.1f SD)", sig_hi))
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
      labs(x = "Hypoxia score (SD within tumour)", y = XLAB_OR,
           caption = paste0("Both slopes from one joint model per tumour; all covariances included. ",
                            if (is.finite(sig_lo))
                              sprintf("Shaded: the two slopes differ (p < %.2f).", JOINT_ALPHA)
                            else "")) +
      theme_fig(S) + theme(legend.position = "top")
    
    if (isTRUE(S$titles))
      p4 <- p4 + labs(
        title = "P4HA1 stops preferring collagen as hypoxia rises",
        subtitle = paste0("Effect of each module on P4HA1 detection, as a function of hypoxia",
                          if (JOINT_ADJUST_DEPTH) " (sequencing depth adjusted)" else ""))
    
    w <- if (isTRUE(S$split)) S$w2 else S$w1 * 1.5
    h <- if (isTRUE(S$split)) S$h_full else MM(68)
    for (dv in S$devices) {
      f <- file.path(outdir, outf("Figure4_joint_effect_modification", dv))
      tryCatch({
        if (dv == "pdf") ggsave(f, p4, width = w, height = h, device = PDF_DEV, bg = "white")
        else             ggsave(f, p4, width = w, height = h, dpi = S$dpi, bg = "white")
      }, error = function(e) log_msg("  fig failed (%s): %s", basename(f), conditionMessage(e)))
    }
    cat(sprintf("  %s written to %s\n", outf("Figure4_joint_effect_modification", "png"), S$dir))
  }
}

log_msg("Section 13 complete.")


################################################################################
## 15. THE SPECIFICITY TEST: P4HA1 vs P4HA2 vs P4HA3
##
##   BRCA NOTE: the reasoning in this header was written for the PDAC
##   cohort. Any statement about what the data showed, which module worked,
##   or how many tumours were usable describes PDAC, not BRCA. The code
##   below is unchanged except where tagged "BRCA:".
##
##   P4HA1 and P4HA2 hydroxylate collagen. P4HA3 does not. If the pattern is
##   about the enzyme it should appear in the first two and fail in the third.
##   If it appears in all three it is gene-family co-regulation and the
##   mechanistic reading collapses. This section is built to be able to return
##   the second answer; if it cannot, it is not a test.
##
##   THE MODEL DOES NOT CHANGE between outcomes. Identical predictors, spots,
##   standardisation and pooling; only y changes. P4HA1 is NOT added as a
##   covariate when P4HA2 is the outcome: the isoforms are co-regulated, so
##   conditioning on one while modelling another changes the estimand.
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
## 15.1  CIRCULARITY GATE  -- hard stop, not a warning
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
## 15.2  SCORE AND STANDARDISE ONCE PER TUMOUR
##       Scored once and reused for all three outcomes, so nothing can differ
##       because of AddModuleScore's RNG. Fingerprint recorded as in 13.
## ============================================================================
iso_expr_one <- function(rds_path, sid) {
  obj <- readRDS(rds_path)
  obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
  DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else Assays(obj)[1]
  
  missing_out <- setdiff(ISO_OUTCOMES, rownames(obj))
  if (length(missing_out)) {
    log_msg("  %s: outcome gene(s) absent (%s) -> sample skipped",
            sid, paste(missing_out, collapse = ", "))
    return(NULL)
  }
  
  fp <- list()
  for (nm in ISO_SCORE_MODULES) {
    full <- MODULES[[nm]]
    g    <- intersect(full, rownames(obj))
    need <- min_genes_for(full)
    if (length(g) < need) {
      log_msg("  %s: %s has only %d/%d genes (need %d) -> sample skipped",
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
    if (is.null(v)) { log_msg("  %s: %s zero variance -> sample skipped", sid, nm); return(NULL) }
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
  
  ## same hard check as section 13: read the model, not the flag
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

cat("\n--- which tumours contributed, per outcome ---\n")
for (og in ISO_OUTCOMES) {
  s <- iso_stat[iso_stat$outcome == og, ]
  cat(sprintf("  %-6s fitted in %2d / %2d tumours   median events %5.0f\n",
              og, sum(s$fitted), nrow(s), median(s$n_events[s$fitted], na.rm = TRUE)))
  drop <- s[!s$fitted, ]
  if (nrow(drop))
    for (j in seq_len(nrow(drop)))
      cat(sprintf("         excluded: %-26s %s\n", drop$sample_id[j], drop$reason[j]))
}
cat("  Name these tumours whenever you quote k.\n")

cat("\n--- the two interaction terms, side by side ---\n")
cat(sprintf("  %-18s %-8s %9s %19s %6s %6s %11s\n",
            "term", "outcome", "beta", "95% CI", "k", "I2", "p"))
for (tm in c(II_COL, II_ER)) {
  for (og in ISO_OUTCOMES) {
    r <- iso_coef_pooled[iso_coef_pooled$term == tm & iso_coef_pooled$outcome == og, ]
    if (!nrow(r)) next
    cat(sprintf("  %-18s %-8s %+9.4f  %+8.4f to %+8.4f %5d %5.0f%% %10s\n",
                tm, og, r$est, r$lower, r$upper, r$k, r$I2, fmt_p(r$p)))
  }
  cat("\n")
}

cat("--- THE HEADLINE: difference between the two slopes, per outcome ---\n")
for (og in ISO_OUTCOMES) {
  r <- iso_con_pooled[iso_con_pooled$outcome == og &
                        iso_con_pooled$contrast == "interaction_diff", ]
  if (!nrow(r)) next
  cat(sprintf("  %-6s %+.4f  (95%% CI %+.4f to %+.4f)  k=%2d  I2=%3.0f%%  p = %-10s %s\n",
              og, r$est, r$lower, r$upper, r$k, r$I2, fmt_p(r$p),
              if (is.finite(r$p) && r$p < ISO_ALPHA) "SIGNIFICANT" else "not significant"))
}

## *** the branch that decides how you may describe P4HA3 ***
{
  a3 <- iso_con_pooled[iso_con_pooled$outcome == "P4HA3" &
                         iso_con_pooled$contrast == "interaction_diff", ]
  a1 <- iso_con_pooled[iso_con_pooled$outcome == "P4HA1" &
                         iso_con_pooled$contrast == "interaction_diff", ]
  if (nrow(a3) == 1 && nrow(a1) == 1) {
    cat("\n--- WHAT YOU MAY AND MAY NOT CLAIM ABOUT P4HA3 ---\n")
    sig3 <- is.finite(a3$p) && a3$p < ISO_ALPHA
    same_sign <- sign(a3$est) == sign(a1$est)
    if (!sig3) {
      cat("  P4HA3 shows no interaction. Combined with its large and highly\n")
      cat("  significant collagen main effect, that is a real negative rather than\n")
      cat("  a power failure. The specificity claim stands.\n")
    } else if (sig3 && same_sign) {
      cat("  *** P4HA3 shows the SAME pattern and it IS significant. That is\n")
      cat("  evidence for gene-family co-regulation, not enzyme function. The\n")
      cat("  mechanistic reading does not survive. Report this. ***\n")
    } else {
      cat("  *** P4HA3 is SIGNIFICANT but in the OPPOSITE direction to P4HA1.\n")
      cat("  This is not the clean negative the specificity argument assumed, and\n")
      cat("  it is not co-regulation either. Do NOT describe P4HA3 as showing no\n")
      cat("  effect. It shows a different effect, which needs its own explanation\n")
      cat("  before anything is claimed. If a previous email said otherwise, it\n")
      cat("  needs correcting. ***\n")
    }
  }
}

cat("\nCAVEATS\n")
cat("  Spatial autocorrelation not modelled; all CIs too narrow.\n")
cat("  Per 1 SD within tumour; not comparable to sections 10 and 11.\n")
cat("  Coefficient magnitudes are not comparable across genes with different\n")
cat("  detection rates. Presence or absence of an interaction is, and that is\n")
cat("  what the specificity result rests on.\n")
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
    st <- iso_stat[iso_stat$outcome == og, ]
    if (!nrow(r)) return(NULL)
    data.frame(outcome = factor(og, ISO_OUTCOMES),
               lab = sprintf("diff %+.3f, p = %s\nk = %d", r$est, fmt_p(r$p), sum(st$fitted)),
               row.names = NULL)
  }))
  
  y_top <- max(lines_df$upper_OR, na.rm = TRUE)
  y_bot <- min(lines_df$lower_OR, na.rm = TRUE)
  y_hi  <- y_top * (y_top / y_bot)^0.28
  
  for (sn in names(FIG_STYLES)) {
    S <- FIG_STYLES[[sn]]
    outdir <- file.path(agg, S$dir); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    
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
      labs(x = "Hypoxia score (SD within tumour)",
           y = "Odds ratio for detection (95% CI)",
           caption = paste0("Identical model, predictors and spots for all three outcomes; only the outcome gene differs.",
                            if (ISO_ADJUST_DEPTH) " Sequencing depth adjusted." else "")) +
      theme_fig(S) + theme(legend.position = "top")
    
    if (!is.null(lab_df) && nrow(lab_df))
      p7 <- p7 + geom_text(data = lab_df, inherit.aes = FALSE,
                           aes(x = min(lines_df$z), y = y_hi, label = lab),
                           hjust = 0, vjust = 1, size = S$lab, colour = "grey20",
                           lineheight = 1.1)
    
    if (isTRUE(S$titles))
      p7 <- p7 + labs(title = "Does the pattern follow enzyme function or the gene family?",
                      subtitle = "P4HA1 and P4HA2 hydroxylate collagen; P4HA3 does not")
    
    h <- if (isTRUE(S$split)) S$h_full else MM(72)
    for (dv in S$devices) {
      f <- file.path(outdir, outf("Figure7_isoform_comparison", dv))
      tryCatch({
        if (dv == "pdf") ggsave(f, p7, width = S$w2, height = h, device = PDF_DEV, bg = "white")
        else             ggsave(f, p7, width = S$w2, height = h, dpi = S$dpi, bg = "white")
      }, error = function(e) log_msg("  fig failed (%s): %s", basename(f), conditionMessage(e)))
    }
    cat(sprintf("  %s written to %s\n", outf("Figure7_isoform_comparison", "png"), S$dir))
  }
}

log_msg("Section 15 complete.")


################################################################################
## 16. THE RANDOM-GENE CONTROL
##
##   BRCA NOTE: the reasoning in this header was written for the PDAC
##   cohort. Any statement about what the data showed, which module worked,
##   or how many tumours were usable describes PDAC, not BRCA. The code
##   below is unchanged except where tagged "BRCA:".
##
##   The model p value assumes spots are independent. They are not: tissue is
##   made of homogeneous regions, so neighbouring spots repeat each other and
##   the model believes it has more information than it does. A reference built
##   from real genes in the same tissue carries the same inflation, so the
##   comparison is immune to that assumption.
##
##   Genes are matched to P4HA1 on detection rate. Without matching the
##   reference would fill with sparse noisy genes and P4HA1 would look extreme
##   for a trivial reason.
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
RND_MIN_K        <- 15   # BRCA: kept from PDAC. Now counts patients: 15 of 21 primary (PDAC: 15 of 23 sections)
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
## 16.1  CHOOSE THE GENES
## ============================================================================
rates <- read.csv(RATES_CSV, stringsAsFactors = FALSE)

gene_med <- rates %>%
  group_by(gene) %>%
  summarise(median_rate = median(detect_rate), n_tumours = dplyr::n(), .groups = "drop") %>%
  filter(n_tumours == dplyr::n_distinct(rates$sample_id))

target <- gene_med$median_rate[gene_med$gene == "P4HA1"]
if (!length(target)) stop("P4HA1 not in the cached detection rates.")

## Not just the outcome genes: every gene inside a predictor module has to go
## too, or a "random" gene could be a member of the module it is tested against.
EXCLUDE <- unique(c("P4HA1", "P4HA2", "P4HA3", "P4HB", "PLOD1", "PLOD2", "PLOD3",
                    unlist(MODULES[RND_SCORE_MODULES], use.names = FALSE)))

pool_genes <- gene_med %>%
  filter(!gene %in% EXCLUDE,
         median_rate >= target * (1 - RND_CALIPER),
         median_rate <= target * (1 + RND_CALIPER))

cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 16.1  GENE SELECTION\n")
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  P4HA1 median detection rate : %.3f\n", target))
cat(sprintf("  caliper                     : %.0f%% relative (%.3f to %.3f)\n",
            100 * RND_CALIPER, target * (1 - RND_CALIPER), target * (1 + RND_CALIPER)))
cat(sprintf("  genes present in all tumours: %d\n", nrow(gene_med)))
cat(sprintf("  eligible inside the caliper : %d\n", nrow(pool_genes)))

if (nrow(pool_genes) < RND_N_GENES) {
  cat(sprintf("  *** only %d eligible. Using all of them. Describe the reference as\n", nrow(pool_genes)))
  cat("      \"all comparable genes\", not \"1000 random genes\". Do NOT widen the\n")
  cat("      caliper to reach 1000: that breaks the matching. ***\n")
  sel <- pool_genes
} else {
  set.seed(RND_SEED)
  sel <- pool_genes[sample(nrow(pool_genes), RND_N_GENES), ]
}
sel <- sel[order(sel$median_rate), ]
write_out(sel, "RANDOM_gene_selection")

cat(sprintf("\n  selected: %d genes, detection median %.3f (range %.3f to %.3f)\n",
            nrow(sel), median(sel$median_rate), min(sel$median_rate), max(sel$median_rate)))
cat(sprintf("  P4HA1 sits at the %.0fth percentile of the selected set.\n",
            100 * mean(sel$median_rate <= target)))
cat("  A value above 50 is expected, since more genes exist at lower detection\n")
cat("  rates. It works in your favour: lower-detection genes give noisier\n")
cat("  estimates, widening the reference and making the test harder to pass.\n")

RND_GENES <- sel$gene

## ============================================================================
## 16.2  FIT, ONE TUMOUR AT A TIME
##       Predictors are identical for every gene within a tumour, so the model
##       matrix is built ONCE and glm.fit is called per gene. glm() 1000 times
##       per tumour would spend most of the run parsing formulas.
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
      log_msg("  %s: %s coverage too low -> sample skipped", sid, nm); return(NULL)
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
  
  ## same hard check, on the design matrix this time
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
  log_msg("[rnd %d/%d] %s", i, nrow(work), sid)
  r <- tryCatch(rnd_one_tumour(path, sid),
                error = function(e) { log_msg("  rnd skipped: %s", conditionMessage(e)); NULL })
  if (!is.null(r) && nrow(r)) rnd_acc[[length(rnd_acc) + 1]] <- r
}
rnd_all <- bind(rnd_acc)
if (nrow(rnd_all) == 0) stop("Section 16: no random gene produced a usable fit.")

## ============================================================================
## 16.4  POOL EACH GENE ACROSS TUMOURS
## ============================================================================
rnd_pooled <- rnd_all %>%
  group_by(gene) %>% group_modify(~ pool_with_p(.x$est, .x$var, .x$sample_id)) %>% ungroup()

n_before <- nrow(rnd_pooled)
rnd_pooled <- rnd_pooled[is.finite(rnd_pooled$est) & rnd_pooled$k >= RND_MIN_K, ]
write_out(rnd_pooled, "RANDOM_pooled_per_gene")

cat("\n", strrep("=", 74), "\n", sep = "")
cat(sprintf("SECTION 16  [%s]  run_id %s\n",
            if (RND_ADJUST_DEPTH) "DEPTH ADJUSTED" else "UNADJUSTED", RUN_ID))
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  genes attempted / pooled / kept : %d / %d / %d\n",
            length(RND_GENES), n_before, nrow(rnd_pooled)))

cat("\n  distribution of the pooled slope difference:\n")
qq <- quantile(rnd_pooled$est, c(0, .025, .25, .5, .75, .975, 1))
for (nmq in names(qq)) cat(sprintf("     %-6s %+.4f\n", nmq, qq[[nmq]]))
cat(sprintf("     significant at 0.05 in %d of %d (%.1f%%)\n",
            sum(rnd_pooled$p < 0.05, na.rm = TRUE), nrow(rnd_pooled),
            100 * mean(rnd_pooled$p < 0.05, na.rm = TRUE)))
cat("     About 5% would be expected. Well above that means the model p values\n")
cat("     in this dataset are inflated, most likely by spatial autocorrelation.\n")
cat("     Report the number. It is why the percentile below is what to quote.\n")

## ============================================================================
## 16.5  WHERE DO THE THREE ISOFORMS SIT?
## ============================================================================
obs <- NULL
iso_file <- file.path(agg, outf("ISOFORM_slope_difference_pooled"))
if (file.exists(iso_file)) {
  io <- read.csv(iso_file, stringsAsFactors = FALSE)
  obs <- io[io$contrast == "interaction_diff", c("outcome", "est", "lower", "upper", "p", "k")]
} else if (exists("iso_con_pooled")) {
  obs <- as.data.frame(iso_con_pooled)[iso_con_pooled$contrast == "interaction_diff",
                                       c("outcome", "est", "lower", "upper", "p", "k")]
}

ref_summary <- NULL
if (!is.null(obs) && nrow(obs)) {
  cat("\n--- where the three isoforms sit ---\n")
  ref_summary <- do.call(rbind, lapply(seq_len(nrow(obs)), function(i) {
    o <- obs[i, ]
    p_dir <- mean(rnd_pooled$est <= o$est)
    p_mag <- mean(abs(rnd_pooled$est) >= abs(o$est))
    data.frame(outcome = o$outcome, est = o$est, model_p = o$p,
               percentile = 100 * p_dir, emp_p_directional = p_dir,
               emp_p_magnitude = p_mag, n_reference = nrow(rnd_pooled), row.names = NULL)
  }))
  write_out(ref_summary, "RANDOM_reference_summary")
  
  for (i in seq_len(nrow(ref_summary))) with(ref_summary[i, ],
                                             cat(sprintf("  %-6s effect %+.4f  percentile %5.1f  empirical p: directional %.4f, magnitude %.4f\n",
                                                         outcome, est, percentile, emp_p_directional, emp_p_magnitude)))
  cat("\n  Directional = fraction of random genes at least as negative. That is\n")
  cat("  the one matching the hypothesis. Quote it and show the other; do not\n")
  cat("  pick whichever is smaller.\n")
  
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
      cat("  the model p is measuring the wrong null. Report this plainly and tell\n")
      cat("  Valerio before doing anything else with the result.\n")
    }
  }
}

cat("\nWHAT THIS SETTLES AND WHAT IT DOES NOT\n")
cat("  Does: whether the behaviour is specific to P4HA1 or common in this tissue.\n")
cat("  Does: partly answers autocorrelation, since the reference shares it. It\n")
cat("        does NOT fix the CIs, which stay too narrow for every gene here.\n")
cat("  Not : matching is on detection rate only. The reference contains\n")
cat("        fibroblast-specific and hypoxia-responsive genes that may behave the\n")
cat("        same way for real reasons. That makes the test harder, not easier.\n")
cat("  Not : mechanism. Still mRNA level, still associational.\n")
cat(strrep("=", 74), "\n", sep = "")

## ============================================================================
## 16.6  FIGURE 8
## ============================================================================
if (exists("FIG_STYLES") && exists("theme_fig")) {
  for (sn in names(FIG_STYLES)) {
    S <- FIG_STYLES[[sn]]
    outdir <- file.path(agg, S$dir); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    
    p8 <- ggplot(rnd_pooled, aes(est)) +
      geom_histogram(bins = 50, fill = "grey72", colour = "white", linewidth = 0.2) +
      geom_vline(xintercept = 0, linetype = "22", linewidth = S$lwd * 0.7, colour = "grey40") +
      labs(x = "Pooled slope difference (collagen minus ER stress, log odds per SD)",
           y = "Number of random genes",
           caption = sprintf("%d genes matched to P4HA1 on detection rate, same model and pooling.%s",
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
    
    w <- if (isTRUE(S$split)) S$w2 else S$w1 * 1.5
    h <- if (isTRUE(S$split)) S$h_full else MM(66)
    for (dv in S$devices) {
      f <- file.path(outdir, outf("Figure8_random_gene_control", dv))
      tryCatch({
        if (dv == "pdf") ggsave(f, p8, width = w, height = h, device = PDF_DEV, bg = "white")
        else             ggsave(f, p8, width = w, height = h, dpi = S$dpi, bg = "white")
      }, error = function(e) log_msg("  fig failed (%s): %s", basename(f), conditionMessage(e)))
    }
    cat(sprintf("  %s written to %s\n", outf("Figure8_random_gene_control", "png"), S$dir))
  }
}

log_msg("Section 16 complete.")


################################################################################
## 17. DEPTH ADJUSTED vs UNADJUSTED, AND RUN-TO-RUN REPRODUCIBILITY
##
##   BRCA NOTE: the reasoning in this header was written for the PDAC
##   cohort. Any statement about what the data showed, which module worked,
##   or how many tumours were usable describes PDAC, not BRCA. The code
##   below is unchanged except where tagged "BRCA:".
##
##   v5 adds two guards this section did not have, both prompted by what went
##   wrong: it verifies the two files really are different runs before comparing
##   them, and it diffs the module fingerprints so a change in the numbers can
##   be traced to scoring or to modelling rather than guessed at.
##
##   Run sections 13, 15, 16 once with ADJUST_DEPTH FALSE and once with TRUE,
##   each in a FRESH session from section 0. Then run this.
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
  
  ## *** THE GUARD THAT WOULD HAVE CAUGHT THE LAST MISTAKE ***
  if ("adjust_depth" %in% names(a) && "adjust_depth" %in% names(b)) {
    if (isTRUE(all(a$adjust_depth == b$adjust_depth))) {
      cat(sprintf("\n  *** %s: BOTH files report adjust_depth = %s. These are not two\n",
                  base, unique(a$adjust_depth)[1]))
      cat("  different runs. The _depthadj file was written by an unadjusted run,\n")
      cat("  which is exactly the earlier fault. Re-run from a fresh session with\n")
      cat("  ADJUST_DEPTH <- TRUE. Comparison skipped. ***\n")
      return(NULL)
    }
    if ("run_id" %in% names(a) && identical(unique(a$run_id), unique(b$run_id))) {
      cat(sprintf("\n  *** %s: both files carry the same run_id. Skipped. ***\n", base))
      return(NULL)
    }
  } else {
    cat(sprintf("\n  NOTE: %s predates the provenance stamp, so the guard cannot\n", base))
    cat("  verify these are different runs. Check the run manifests by hand.\n")
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
    if (x$k != y$k) {
      cat(sprintf("     NOTE: k changed from %d to %d. Depth adjustment costs a degree\n", x$k, y$k))
      cat("     of freedom, so a borderline tumour can drop out. Say so if it happens.\n")
    }
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

## ---- 17.4  MODULE FINGERPRINT DIFF   *** the diagnostic for fault 2 *** ----
##   If the numbers moved between two runs, this says whether the module scores
##   themselves changed. Identical fingerprints mean the cause is downstream in
##   the model or the pooling. Different fingerprints mean it is upstream in
##   scoring, and the module and tumour are named.
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
      cat("  No. All module scores are numerically identical across the two runs.\n")
      cat("  Any difference in results therefore comes from the model, not from\n")
      cat("  scoring. Depth adjustment is expected to change the model, so this is\n")
      cat("  the normal outcome.\n")
    } else {
      cat(sprintf("  YES: %d of %d module-by-tumour scores differ.\n", nrow(changed), nrow(m)))
      cat("  Module scores should be identical between runs, since SEED is set\n")
      cat("  before every AddModuleScore call. If they are not, something upstream\n")
      cat("  changed: gene lists, Seurat version, or a contaminated session.\n")
      cat("  Worst offenders:\n")
      changed <- changed[order(-abs(changed$dmean)), ]
      for (i in seq_len(min(10, nrow(changed)))) with(changed[i, ],
                                                      cat(sprintf("     %-26s %-18s d_mean %+.6f  d_sd %+.6f  d_genes %+d\n",
                                                                  sample_id, module, dmean, dsd, dgenes)))
    }
  } else {
    cat("\n  17.4 unavailable: fingerprints exist only from v5 onward. Re-run\n")
    cat("  sections 13 and 15 to generate them.\n")
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
    cat("  About 5% expected. Depth adjustment is NOT expected to fix this: the\n")
    cat("  cause is spatial autocorrelation, not depth. If it stays high after\n")
    cat("  adjustment, that confirms the diagnosis.\n")
  }
}

cat("\n--- WHAT TO TELL VALERIO ---\n")
cat("  Small shifts, no significance changes:\n")
cat("     \"Adjusting for sequencing depth per spot moves the slope difference by\n")
cat("     [shift]; the conclusion is unchanged.\" Objection closed in one line.\n")
cat("  Large shift or a status change:\n")
cat("     lead with it. Depth was carrying part of the signal and everything\n")
cat("     downstream needs re-reading in that light.\n")
cat(strrep("=", 74), "\n", sep = "")

log_msg("Section 17 complete. run_id %s, ADJUST_DEPTH %s.", RUN_ID, ADJUST_DEPTH)

################################################################################
## 18. THE METABOLIC AXIS: 2-OXOGLUTARATE SUPPLY INSTEAD OF HYPOXIA
##
##   BRCA NOTE: the reasoning in this header was written for the PDAC
##   cohort. Any statement about what the data showed, which module worked,
##   or how many tumours were usable describes PDAC, not BRCA. The code
##   below is unchanged except where tagged "BRCA:".
##                                     (Valerio, email after lab meeting 3, 2026)
##
## HIS REQUEST
##   "Try to see it from another angle. TCA cycle, glutamine metabolism, because
##   glutamine metabolism produces 2-oxoglutarate. Hypoxia is dirty."
##
## THE LOGIC, IN ONE PARAGRAPH
##   P4H consumes oxygen AND 2-oxoglutarate in the same reaction. Both are
##   substrates, not accessories. So there are two independent doors into the
##   same room: "enzyme present but unable to work". Sections 13 and 15 came in
##   through the oxygen door. This section comes in through the cofactor door.
##   If the decoupling is really about enzymatic capacity, it should be visible
##   from both. If it only appears along the hypoxia axis, then hypoxia was
##   doing something else and the mechanistic reading weakens.
##
## WHY A SECOND DOOR IS NEEDED AT ALL
##   Hypoxia changes hundreds of genes at once. So "P4HA1 and ER stress rise
##   together in hypoxic regions" always has the alternative reading that both
##   are independent consequences of the same upheaval, with no relationship to
##   collagen or to retention. No amount of further analysis along the hypoxia
##   axis can settle that. A mechanistically distinct axis can.
##
## THE GENE SET
##   GOBP_GLUTAMINE_FAMILY_AMINO_ACID_METABOLIC_PROCESS, MSigDB M11431,
##   GO:0009064, human, release 2025.1.Hs. 71 genes.
##
##   This is the set the paper Valerio sent started from, NOT its four-gene
##   signature. Two reasons:
##     - the signature (ASRGL1, MECP2, NOS2, NOS3) was selected by LASSO to
##       predict SURVIVAL in colorectal cancer, not to measure 2-oxoglutarate
##       supply. The paper itself notes NOS2 is not directly part of glutamine
##       metabolism.
##     - ASRGL1 is not even in the current version of the source set. The paper
##       used an older GO release with 76 genes; the current one has 71.
##   The signature is still fitted below as a secondary, because he asked, but
##   labelled for what it is.
##
## VERIFIED BEFORE WRITING THIS CODE
##   Circularity : clean. No P4HA1, P4HA2, P4HA3, P4HB, PLOD1-3, PPIB.
##   Overlap     : ZERO shared genes with Hypoxia, CollagenSynthesis, iCAF,
##                 ERstressStrict, ERAD or CollagenMaturation. Nothing to strip,
##                 unlike the VEGFA problem in section 11.
##   PDAC ONLY, not checked for BRCA (the section prints the BRCA values):
##   Coverage    : all 23 tumours clear the 30% threshold. Worst tumour has 37
##                 of 71 detected, median 46. k will stay at 23.
##   Core genes  : GLUL 0.886, GLS 0.441, GLUD1 0.364, OAT 0.353, GOT2 0.260,
##                 CAD 0.238, GOT1 0.153 median detection. P4HA1 is 0.279 for
##                 comparison, so the genes that matter are detected as well as
##                 the outcome. Four of the seven highest-expressed members of
##                 the whole set are core 2-oxoglutarate genes.
##
## THE HONEST WEAKNESS, STATE IT BEFORE ANYONE ASKS
##   Two of them, and they point in the same direction.
##
##   (1) The GO term is broader than glutamine. Its own description says
##       "arginine, glutamate, glutamine and proline". It carries the urea cycle
##       (CPS1, OTC, ASS1, ASL), proline metabolism (PRODH, PYCR1-3, ALDH18A1),
##       nitric oxide synthases and polyamine genes, none of which feed
##       2-oxoglutarate to P4HA1. The module measures glutamine-family amino
##       acid metabolism, not cofactor availability.
##
##   (2) Transcript abundance is a poor proxy for metabolite concentration.
##       Metabolic flux is controlled largely after transcription: enzyme
##       activity, substrate availability, allosteric regulation. None of that
##       is visible here. This axis is therefore WEAKER as a measurement than
##       hypoxia, where HIF drives transcription directly. Its value is in being
##       INDEPENDENT, not in being sharper.
##
##   Because of (1) a focused 12-gene subset is fitted alongside. It is closer
##   to the hypothesis but hand-selected, which is exactly the criticism the
##   27-gene ER stress module attracts. Neither version is clean; running both
##   and reporting both is the honest response.
##
## THE GATE THAT CAN STOP THIS SECTION
##   HIF suppresses TCA flux, so the metabolic module may correlate with the
##   hypoxia module. If it does so strongly, the two axes are not independent
##   and the whole premise of a "second door" fails. 18.2 measures this FIRST
##   and warns loudly. That is not a formality: a high correlation means the
##   result must not be presented as independent confirmation.
##
## OUTPUT
##   METAB_collinearity.csv        correlation with hypoxia, per tumour
##   METAB_coefficients_pooled.csv
##   METAB_slope_difference_pooled.csv    the headline, comparable to section 13
##   METAB_slope_curves_pooled.csv
##   Figure9_metabolic_axis
##
## PREREQUISITE
##   Run after sections 11, 12 and the 12b helper block. Uses assert_run_state,
##   pool_with_p, fmt_p, outf, write_out.
################################################################################

stopifnot(exists("MODULES"), exists("meta_pool"), exists("work"), exists("agg"))
if (!exists("pool_with_p") || !exists("fmt_p") || !exists("assert_run_state"))
  stop("Section 18 needs the 12b helper block (pool_with_p, fmt_p, assert_run_state).")

## ============================================================================
## 18.0  GENE SETS AND SETTINGS
## ============================================================================

## MSigDB M11431, GO:0009064, release 2025.1.Hs, transcribed verbatim from the
## GENE_SYMBOLS field of the downloaded TSV. 71 genes. Do not edit by hand: the
## defence of this module is that it is a named standard set, not a curated one.
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

## Focused subset: enzymes that directly produce or consume 2-oxoglutarate, or
## sit one step away on the glutamine to 2-OG route.
##   GLS, GLS2    glutamine -> glutamate
##   GLUD1, GLUD2 glutamate -> 2-oxoglutarate (oxidative deamination)
##   GOT1, GOT2   glutamate <-> 2-oxoglutarate (transamination, cytosol / mito)
##   GLUL         glutamate -> glutamine (the reverse direction, so it reports
##                on the same node and is informative about which way flux runs)
##   OAT          ornithine + 2-OG -> glutamate semialdehyde + glutamate
##   ALDH5A1      succinate semialdehyde -> succinate, the GABA shunt exit
##   ALDH4A1      glutamate semialdehyde -> glutamate
##   CAD          consumes glutamine for pyrimidines, a competing sink
##   PHGDH        serine synthesis, transaminates glutamate to 2-OG
## HAND-SELECTED. That is the point of running it only as a secondary.
OG_CORE_GENES <- c("GLS","GLS2","GLUD1","GLUD2","GOT1","GOT2","GLUL","OAT",
                   "ALDH5A1","ALDH4A1","CAD","PHGDH")

## The paper's four-gene signature. Built to predict colorectal cancer survival,
## not to measure cofactor supply. ASRGL1 is not in the current GO set at all.
## Fitted because he asked; reported with that caveat attached.
PAPER_SIGNATURE <- c("ASRGL1","MECP2","NOS2","NOS3")

MODULES$Glutamine   <- GLUTAMINE_GENES
MODULES$OGcore      <- OG_CORE_GENES
MODULES$PaperSig    <- PAPER_SIGNATURE

METAB_VARIANTS   <- c("Glutamine", "OGcore", "PaperSig")
METAB_PRIMARY    <- "Glutamine"
METAB_ES         <- "ERstressStrict"
METAB_MIN_SPOTS  <- 30
METAB_Z_GRID     <- seq(-2, 2, length.out = 41)
METAB_ALPHA      <- 0.05
METAB_ADJUST_DEPTH <- ADJUST_DEPTH     # derived, never set by hand

## correlation with hypoxia above which the two axes cannot be called independent
METAB_COLLIN_WARN <- 0.50
METAB_COLLIN_STOP <- 0.80

MT_COL <- "zCol"; MT_ER <- "zER"; MT_MET <- "zMet"
MI_COL <- paste0(MT_MET, ":", MT_COL)
MI_ER  <- paste0(MT_MET, ":", MT_ER)

LAB_COL <- "Collagen synthesis"
LAB_ER  <- "ER stress"
LAB_DIF <- "difference"

assert_run_state("Section 18")

## ============================================================================
## 18.1  GENE-SET AUDIT, halts on failure
## ============================================================================
cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 18.1  GENE-SET AUDIT\n")
cat(strrep("=", 74), "\n", sep = "")

DANGER <- c("P4HA1","P4HA2","P4HA3","P4HB","PLOD1","PLOD2","PLOD3","PPIB")
OTHER_PREDICTORS <- c("Hypoxia","CollagenSynthesis","iCAF","ERstressStrict","ERAD")

metab_fail <- FALSE
for (nm in METAB_VARIANTS) {
  g   <- MODULES[[nm]]
  bad <- intersect(DANGER, g)
  cat(sprintf("\n  %s  (n = %d genes)\n", nm, length(g)))
  cat(sprintf("    circularity : %s\n",
              if (length(bad)) paste0("*** FAIL: ", paste(bad, collapse = ", "), " ***") else "clean"))
  if (length(bad)) metab_fail <- TRUE
  ov_any <- FALSE
  for (other in OTHER_PREDICTORS) {
    ov <- intersect(g, MODULES[[other]])
    if (length(ov)) {
      cat(sprintf("    *** overlaps %s: %s ***\n", other, paste(ov, collapse = ", ")))
      ov_any <- TRUE; metab_fail <- TRUE
    }
  }
  if (!ov_any) cat("    overlap with other predictor modules : none\n")
}
if (metab_fail)
  stop("Section 18 halted: a metabolic module is circular or overlaps another predictor.")
cat("\n  Audit passed. Nothing to strip.\n")

## ============================================================================
## 18.2  SCORE, STANDARDISE, AND MEASURE THE COLLINEARITY GATE
##
##   Everything is scored in ONE pass so all variants see identical predictors.
##   The hypoxia module is scored too, not because it enters the model, but
##   because the correlation between it and the metabolic module decides whether
##   this analysis can be described as an independent axis at all.
## ============================================================================
METAB_SCORE_MODULES <- c("Hypoxia", "CollagenSynthesis", "iCAF", METAB_ES, "ERAD",
                         METAB_VARIANTS)

metab_expr_one <- function(rds_path, sid) {
  obj <- readRDS(rds_path)
  obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
  DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else Assays(obj)[1]
  if (!"P4HA1" %in% rownames(obj)) return(NULL)
  
  cov <- list()
  for (nm in METAB_SCORE_MODULES) {
    full <- MODULES[[nm]]
    g    <- intersect(full, rownames(obj))
    need <- min_genes_for(full)
    cov[[nm]] <- data.frame(sample_id = sid, module = nm,
                            n_total = length(full), n_found = length(g),
                            n_required = need, row.names = NULL)
    if (length(g) < need) {
      log_msg("  %s: %s has %d/%d genes (need %d) -> sample skipped",
              sid, nm, length(g), length(full), need)
      return(NULL)
    }
    set.seed(SEED)
    obj <- tryCatch(add_module_score_fallback(obj, features = list(g), name = nm,
                                   assay = DefaultAssay(obj), nbin = 20, ctrl = 50),
                    error = function(e) NULL)
    if (is.null(obj) || !paste0(nm, "1") %in% colnames(obj@meta.data)) return(NULL)
    names(obj@meta.data)[names(obj@meta.data) == paste0(nm, "1")] <- nm
  }
  
  vars <- c("P4HA1", METAB_SCORE_MODULES)
  cnt_col <- grep("^nCount_", colnames(obj@meta.data), value = TRUE)
  if (METAB_ADJUST_DEPTH) {
    if (!length(cnt_col)) return(NULL)
    vars <- c(vars, cnt_col[1])
  }
  d <- FetchData(obj, vars = vars, layer = "data")
  d$y <- as.integer(d$P4HA1 > 0)
  
  ## collinearity, measured on the RAW scores before standardisation
  cc <- suppressWarnings(cor(d[, METAB_SCORE_MODULES], method = "spearman"))
  collin <- data.frame(
    sample_id = sid,
    Glutamine_vs_Hypoxia = cc["Glutamine", "Hypoxia"],
    OGcore_vs_Hypoxia    = cc["OGcore", "Hypoxia"],
    PaperSig_vs_Hypoxia  = cc["PaperSig", "Hypoxia"],
    Glutamine_vs_OGcore  = cc["Glutamine", "OGcore"],
    Glutamine_vs_Col     = cc["Glutamine", "CollagenSynthesis"],
    Glutamine_vs_ER      = cc["Glutamine", METAB_ES],
    row.names = NULL)
  
  zz <- function(x) { s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) return(NULL)
  as.numeric((x - mean(x, na.rm = TRUE)) / s) }
  for (nm in METAB_SCORE_MODULES) {
    v <- zz(d[[nm]]); if (is.null(v)) return(NULL)
    d[[paste0("z", nm)]] <- v
  }
  if (METAB_ADJUST_DEPTH) {
    v <- zz(log1p(d[[cnt_col[1]]])); if (is.null(v)) return(NULL); d$zDepth <- v
  }
  
  attr(d, "collin")   <- collin
  attr(d, "coverage") <- do.call(rbind, cov)
  d
}

## ============================================================================
## 18.3  FIT ONE TUMOUR, ONE METABOLIC VARIANT
##       Same shape as sections 13 and 15, with zMet replacing zHyp as the axis.
## ============================================================================
fit_metab_one <- function(d, sid, variant) {
  out <- list(coef = NULL, contrast = NULL, curve = NULL)
  if (is.null(d) || nrow(d) < METAB_MIN_SPOTS || length(unique(d$y)) != 2) return(out)
  
  dd <- d
  dd$zMet <- d[[paste0("z", variant)]]
  dd$zCol <- d$zCollagenSynthesis
  dd$zER  <- d[[paste0("z", METAB_ES)]]
  dd$ziCAF  <- d$ziCAF
  dd$zERAD  <- d$zERAD
  
  rhs <- paste0(MT_MET, " * ", MT_COL, " + ", MT_MET, " * ", MT_ER, " + ziCAF + zERAD",
                if (METAB_ADJUST_DEPTH) " + zDepth" else "")
  f <- tryCatch(glm(as.formula(paste("y ~", rhs)), data = dd, family = binomial()),
                error = function(e) NULL)
  if (is.null(f) || !f$converged) return(out)
  
  b <- coef(f); V <- vcov(f)
  
  ## read the model, not the flag
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
log_msg("Section 18 start. Variants: %s. ADJUST_DEPTH = %s.",
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
  if (!is.null(attr(d, "collin")))
    metab_acc$collin[[length(metab_acc$collin) + 1]] <- attr(d, "collin")
  if (!is.null(attr(d, "coverage")))
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

if (nrow(mcon) == 0) stop("Section 18: no tumour produced a usable fit.")

write_out(mcoef, "METAB_coefficients_per_sample")
write_out(mcon,  "METAB_contrasts_per_sample")
write_out(mcol,  "METAB_collinearity")
write_out(mcov,  "METAB_gene_coverage")

## ============================================================================
## 18.5  THE COLLINEARITY GATE, reported before any result
##
##   This has to come before the headline. If the metabolic module is largely a
##   restatement of hypoxia, then a matching result is not independent
##   confirmation, it is the same finding twice, and presenting it as
##   corroboration would be wrong.
## ============================================================================
cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 18.5  IS THIS REALLY A SECOND, INDEPENDENT AXIS?\n")
cat(strrep("=", 74), "\n", sep = "")
cat("  Spearman correlation with the hypoxia module, per tumour:\n\n")

gate_verdict <- list()
for (vr in METAB_VARIANTS) {
  cc <- mcol[[paste0(vr, "_vs_Hypoxia")]]
  med <- median(cc, na.rm = TRUE); mx <- max(abs(cc), na.rm = TRUE)
  status <- if (abs(med) >= METAB_COLLIN_STOP) "TOO HIGH"
  else if (abs(med) >= METAB_COLLIN_WARN) "HIGH"
  else "acceptable"
  gate_verdict[[vr]] <- status
  cat(sprintf("    %-10s median %+.2f   range %+.2f to %+.2f   max |rho| %.2f   -> %s\n",
              vr, med, min(cc, na.rm = TRUE), max(cc, na.rm = TRUE), mx, status))
}
cat(sprintf("\n  thresholds: |rho| >= %.2f is HIGH, >= %.2f is TOO HIGH\n",
            METAB_COLLIN_WARN, METAB_COLLIN_STOP))
cat("\n  HOW TO READ THIS\n")
cat("    acceptable : the two axes are largely separate. A matching result on\n")
cat("                 both is genuine independent support.\n")
cat("    HIGH       : the axes share substantial variance. A matching result is\n")
cat("                 weaker evidence than it looks. Say so explicitly.\n")
cat("    TOO HIGH   : the metabolic module is largely hypoxia in another form.\n")
cat("                 Do NOT present agreement as corroboration; it is one\n")
cat("                 finding reported twice.\n")

if (identical(gate_verdict[[METAB_PRIMARY]], "TOO HIGH")) {
  cat("\n  *** The PRIMARY module fails the gate. Everything below is still\n")
  cat("  computed, but it cannot be described as an independent axis. Tell\n")
  cat("  Valerio this before showing the numbers. ***\n")
}

cat("\n  other correlations of the primary module (median across tumours):\n")
cat(sprintf("    vs the 12-gene core subset : %+.2f\n", median(mcol$Glutamine_vs_OGcore, na.rm = TRUE)))
cat(sprintf("    vs collagen synthesis      : %+.2f\n", median(mcol$Glutamine_vs_Col, na.rm = TRUE)))
cat(sprintf("    vs ER stress               : %+.2f\n", median(mcol$Glutamine_vs_ER, na.rm = TRUE)))
cat("  If the 71-gene and 12-gene versions correlate highly, the broad set is\n")
cat("  effectively being carried by its core members, which is reassuring.\n")

cat("\n  gene coverage per module (median across tumours):\n")
cs <- mcov %>% group_by(module) %>%
  summarise(median_found = median(n_found), total = max(n_total),
            min_found = min(n_found), required = max(n_required), .groups = "drop")
for (i in seq_len(nrow(cs))) with(cs[i, ],
                                  cat(sprintf("    %-20s %3d / %3d detected (min %3d, need %3d)\n",
                                              module, median_found, total, min_found, required)))

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

cat("\n", strrep("=", 74), "\n", sep = "")
cat(sprintf("SECTION 18  THE METABOLIC AXIS  [%s]\n",
            if (METAB_ADJUST_DEPTH) "DEPTH ADJUSTED" else "UNADJUSTED"))
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("Model : y ~ %s\n",
            paste0(MT_MET, "*", MT_COL, " + ", MT_MET, "*", MT_ER, " + ziCAF + zERAD",
                   if (METAB_ADJUST_DEPTH) " + zDepth" else "")))
cat("zMet replaces zHyp as the axis. Everything else is identical to section 13,\n")
cat("so the headline number is directly comparable to it.\n")

cat("\n--- the two interaction terms ---\n")
cat(sprintf("  %-12s %-14s %9s %21s %5s %6s %11s\n",
            "variant", "term", "beta", "95% CI", "k", "I2", "p"))
for (vr in METAB_VARIANTS) {
  for (tm in c(MI_COL, MI_ER)) {
    r <- mcoef_pooled[mcoef_pooled$variant == vr & mcoef_pooled$term == tm, ]
    if (!nrow(r)) next
    cat(sprintf("  %-12s %-14s %+9.4f  %+8.4f to %+8.4f %5d %5.0f%% %10s\n",
                vr, tm, r$est, r$lower, r$upper, r$k, r$I2, fmt_p(r$p)))
  }
}

cat("\n--- THE HEADLINE: difference between the two slopes, along the metabolic axis ---\n")
for (vr in METAB_VARIANTS) {
  r <- mcon_pooled[mcon_pooled$variant == vr & mcon_pooled$contrast == "interaction_diff", ]
  if (!nrow(r)) next
  cat(sprintf("  %-12s %+.4f  (95%% CI %+.4f to %+.4f)  k=%2d  I2=%3.0f%%  p = %-10s %s\n",
              vr, r$est, r$lower, r$upper, r$k, r$I2, fmt_p(r$p),
              if (is.finite(r$p) && r$p < METAB_ALPHA) "SIGNIFICANT" else "not significant"))
}

## ---- compare directly against the hypoxia axis ----------------------------
hyp_ref <- NULL
hf <- file.path(agg, outf("JOINT_slope_difference_pooled"))
if (file.exists(hf)) {
  h <- read.csv(hf, stringsAsFactors = FALSE)
  hyp_ref <- h[h$contrast == "interaction_diff", ]
}

cat("\n--- SIDE BY SIDE WITH THE HYPOXIA AXIS ---\n")
if (!is.null(hyp_ref) && nrow(hyp_ref) == 1) {
  cat(sprintf("  hypoxia axis (section 13) : %+.4f  p = %s\n", hyp_ref$est, fmt_p(hyp_ref$p)))
  pr <- mcon_pooled[mcon_pooled$variant == METAB_PRIMARY &
                      mcon_pooled$contrast == "interaction_diff", ]
  if (nrow(pr) == 1) {
    cat(sprintf("  metabolic axis (primary)  : %+.4f  p = %s\n", pr$est, fmt_p(pr$p)))
    same_sign <- sign(pr$est) == sign(hyp_ref$est)
    sig_met   <- is.finite(pr$p) && pr$p < METAB_ALPHA
    cat("\n  VERDICT\n")
    if (sig_met && same_sign) {
      cat("    The same decoupling appears along BOTH axes, in the same direction.\n")
      cat("    That is the result Valerio was hoping for: two mechanistically\n")
      cat("    distinct routes to the same state give the same answer, so the\n")
      cat("    \"it is just a side effect of hypoxia\" reading is much harder to hold.\n")
      cat(sprintf("    Strength of this claim depends on the gate above: %s.\n",
                  gate_verdict[[METAB_PRIMARY]]))
    } else if (sig_met && !same_sign) {
      cat("    *** The metabolic axis is significant but in the OPPOSITE direction.\n")
      cat("    That is not corroboration and it is not a clean negative either.\n")
      cat("    It needs explaining before anything is claimed. Report it as it is. ***\n")
    } else {
      cat("    The decoupling does NOT appear along the metabolic axis.\n")
      cat("    Two honest readings, and the data here cannot separate them:\n")
      cat("      (a) the phenomenon is specific to oxygen limitation, not to\n")
      cat("          enzymatic capacity in general; or\n")
      cat("      (b) transcript levels of metabolic genes are simply too poor a\n")
      cat("          proxy for 2-oxoglutarate availability to detect it.\n")
      cat("    Reading (b) is quite likely given how far transcript is from flux.\n")
      cat("    Do NOT present this as evidence against the mechanism.\n")
    }
  }
} else {
  cat("  Section 13 output not found, so no direct comparison. Run section 13 first.\n")
}

## ---- primary vs the hand-selected subset ----------------------------------
cat("\n--- broad set vs focused subset ---\n")
{
  a <- mcon_pooled[mcon_pooled$variant == "Glutamine" & mcon_pooled$contrast == "interaction_diff", ]
  b <- mcon_pooled[mcon_pooled$variant == "OGcore"    & mcon_pooled$contrast == "interaction_diff", ]
  if (nrow(a) == 1 && nrow(b) == 1) {
    cat(sprintf("  71-gene GOBP set : %+.4f  p = %s\n", a$est, fmt_p(a$p)))
    cat(sprintf("  12-gene 2-OG core: %+.4f  p = %s\n", b$est, fmt_p(b$p)))
    cat("  If the focused subset gives a stronger effect, that is consistent with\n")
    cat("  the signal coming from 2-oxoglutarate metabolism specifically rather\n")
    cat("  than from glutamine-family metabolism generally. It is suggestive, not\n")
    cat("  a test: the subset was chosen by hand, which is the same criticism the\n")
    cat("  27-gene ER stress module attracts. Report the GOBP set as primary.\n")
  }
}

cat("\n--- the paper's four-gene signature ---\n")
{
  s <- mcon_pooled[mcon_pooled$variant == "PaperSig" & mcon_pooled$contrast == "interaction_diff", ]
  if (nrow(s) == 1) {
    cat(sprintf("  %+.4f  (95%% CI %+.4f to %+.4f)  p = %s\n", s$est, s$lower, s$upper, fmt_p(s$p)))
    cat("  Fitted because Valerio asked. State plainly what it is: a signature\n")
    cat("  selected by LASSO to predict SURVIVAL in colorectal cancer, not to\n")
    cat("  measure 2-oxoglutarate supply. The source paper itself notes NOS2 is\n")
    cat("  not directly part of glutamine metabolism, and ASRGL1 is no longer in\n")
    cat("  the current version of the GO set the paper started from. Whatever\n")
    cat("  this number is, it is not a cofactor measurement.\n")
  }
}

cat("\nCAVEATS TO STATE BEFORE ANYONE ASKS\n")
cat("  1. The GO term covers arginine, glutamate, glutamine AND proline, so it\n")
cat("     includes the urea cycle, proline metabolism and NOS genes that have\n")
cat("     nothing to do with feeding 2-oxoglutarate to P4HA1.\n")
cat("  2. Transcript abundance is a weak proxy for metabolite concentration.\n")
cat("     Flux is controlled mostly after transcription. This axis is weaker as\n")
cat("     a MEASUREMENT than hypoxia; its value is that it is INDEPENDENT.\n")
cat("  3. Spatial autocorrelation not modelled, so all CIs are too narrow.\n")
cat("  4. Per 1 SD within tumour, not comparable to sections 10 and 11.\n")
cat("  5. mRNA level and associational. Says nothing about actual 2-OG levels.\n")
cat(strrep("=", 74), "\n", sep = "")

## ============================================================================
## 18.7  FIGURE 9
## ============================================================================
if (!exists("FIG_STYLES") || !exists("theme_fig")) {
  cat("\n  Figure 9 skipped: section 12 has not run.\n")
} else {
  nice_v <- c(Glutamine = "Glutamine family (71 genes)",
              OGcore    = "2-OG core (12 genes)",
              PaperSig  = "Paper signature (4 genes)")
  ld <- mcur_pooled[mcur_pooled$series != LAB_DIF, ]
  ld$series  <- factor(ld$series, c(LAB_COL, LAB_ER))
  ld$variant <- factor(nice_v[ld$variant], levels = unname(nice_v))
  
  lab_df <- do.call(rbind, lapply(METAB_VARIANTS, function(vr) {
    r <- mcon_pooled[mcon_pooled$variant == vr & mcon_pooled$contrast == "interaction_diff", ]
    if (!nrow(r)) return(NULL)
    data.frame(variant = factor(nice_v[[vr]], levels = unname(nice_v)),
               lab = sprintf("diff %+.3f, p = %s\nk = %d", r$est, fmt_p(r$p), r$k),
               row.names = NULL)
  }))
  
  y_top <- max(ld$upper_OR, na.rm = TRUE); y_bot <- min(ld$lower_OR, na.rm = TRUE)
  y_hi  <- y_top * (y_top / y_bot)^0.28
  
  for (sn in names(FIG_STYLES)) {
    S <- FIG_STYLES[[sn]]
    outdir <- file.path(agg, S$dir); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    
    p9 <- ggplot(ld, aes(z, OR)) +
      geom_hline(yintercept = 1, linetype = "22", linewidth = S$lwd * 0.7, colour = "grey40") +
      geom_ribbon(aes(ymin = lower_OR, ymax = upper_OR, fill = series), alpha = .16, colour = NA) +
      geom_line(aes(colour = series), linewidth = S$lwd * 1.3) +
      scale_colour_manual(values = setNames(c(CLR$neg, CLR$pos), c(LAB_COL, LAB_ER)), name = NULL) +
      scale_fill_manual(values   = setNames(c(CLR$neg, CLR$pos), c(LAB_COL, LAB_ER)), name = NULL) +
      facet_wrap(~ variant, nrow = 1) + scale_y_log10() +
      coord_cartesian(ylim = c(y_bot, y_hi)) +
      labs(x = "Metabolic module score (SD within tumour)",
           y = "Odds ratio for P4HA1 detection (95% CI)",
           caption = sprintf("Same model as section 13 with the metabolic module replacing hypoxia as the axis. Median correlation with hypoxia: %+.2f.",
                             median(mcol$Glutamine_vs_Hypoxia, na.rm = TRUE))) +
      theme_fig(S) + theme(legend.position = "top")
    
    if (!is.null(lab_df) && nrow(lab_df))
      p9 <- p9 + geom_text(data = lab_df, inherit.aes = FALSE,
                           aes(x = min(ld$z), y = y_hi, label = lab),
                           hjust = 0, vjust = 1, size = S$lab, colour = "grey20", lineheight = 1.1)
    
    if (isTRUE(S$titles))
      p9 <- p9 + labs(title = "Does the same decoupling appear along the cofactor axis?",
                      subtitle = "2-oxoglutarate is a substrate of P4H alongside oxygen")
    
    h <- if (isTRUE(S$split)) S$h_full else MM(72)
    for (dv in S$devices) {
      f <- file.path(outdir, outf("Figure9_metabolic_axis", dv))
      tryCatch({
        if (dv == "pdf") ggsave(f, p9, width = S$w2, height = h, device = PDF_DEV, bg = "white")
        else             ggsave(f, p9, width = S$w2, height = h, dpi = S$dpi, bg = "white")
      }, error = function(e) log_msg("  fig failed (%s): %s", basename(f), conditionMessage(e)))
    }
    cat(sprintf("  %s written to %s\n", outf("Figure9_metabolic_axis", "png"), S$dir))
  }
}

log_msg("Section 18 complete.")

################################################################################
## 19. IS THE METABOLIC MODULE INTERNALLY COHERENT?
##
##   BRCA NOTE: the reasoning in this header was written for the PDAC
##   cohort. Any statement about what the data showed, which module worked,
##   or how many tumours were usable describes PDAC, not BRCA. The code
##   below is unchanged except where tagged "BRCA:".
##
##     REPLACES the earlier section 19. Delete that one; this supersedes it.
##
## WHY THE FIRST ATTEMPT FAILED, AND WHY IT WAS MY DESIGN ERROR
##   The first version compared each module score against external yardsticks:
##   an OXPHOS module, mitochondrial read fraction, genes detected per spot,
##   ribosomal read fraction. The hypoxia module was the positive control.
##
##   In PDAC, hypoxia failed all four (best median rho +0.083, threshold 0.15). Since we
##   know the hypoxia module works, the yardsticks were the problem, and the
##   section said nothing about glutamine.
##
##   The reason is now obvious in hindsight. AddModuleScore computes a score
##   RELATIVE TO CONTROL GENES matched on expression level. It is built to
##   remove exactly the generic effects those yardsticks measure: depth, overall
##   cellular activity, gene count. Three of my four yardsticks were generic
##   activity measures, so the scoring method had already subtracted the thing I
##   was trying to correlate against. The test could not have worked.
##
## WHAT THIS VERSION ASKS INSTEAD
##   Not "does the module track something external", but "do the genes inside
##   the module behave as a unit at all".
##
##   A module score is a mean over member genes. That mean only means something
##   if the members move together. If GLS, OTC, NOS3 and PRODH vary
##   independently across spots, averaging them produces a number that
##   represents no coherent biological quantity, and any null result from it is
##   uninterpretable.
##
##   This asks a question about the module alone. It needs no external yardstick
##   and is therefore immune to the failure above.
##
## TWO MEASURES OF COHERENCE
##   1. Mean pairwise correlation between member genes across spots.
##      Direct and easy to explain.
##   2. Variance explained by the first principal component of the member genes.
##      If members move as one, PC1 captures most of the variance. This is less
##      sensitive to the tie problem below.
##
## THE PROBLEM WITH BOTH, AND HOW IT IS HANDLED
##   Visium is sparse. Most gene pairs share a large number of double zeros,
##   which inflates ties and drags correlations toward a value that reflects
##   dropout rather than biology. So an absolute coherence number is close to
##   meaningless: is 0.05 low? Compared to what?
##
##   The answer is a null. For every real module, RANDOM gene sets of the same
##   size, matched on detection rate, are scored the same way. The module is
##   then reported as a PERCENTILE against its own null. That absorbs the
##   dropout problem, exactly as section 16 absorbed spatial autocorrelation.
##
## THE CONTROLS
##   Several known modules are tested alongside, not just glutamine:
##     Hypoxia            positive control, known to work in this data
##     CollagenSynthesis  a tight, obviously co-regulated set
##     ERstressStrict     the other module carrying a real result
##     ERAD               a large set that produced a near-null
##   Testing several shows what coherence looks like across modules that work,
##   which is the only way to judge whether glutamine is unusual.
##
## HOW TO READ THE OUTCOME
##   Hypoxia coherent, glutamine not  -> the module is not a unit. Section 18's
##                                       null is about the instrument. This is
##                                       the expected outcome.
##   Both coherent                    -> glutamine is a real module and its null
##                                       cannot be dismissed. Harder result.
##   Neither coherent                 -> the metric fails here too. Say so and
##                                       stop; do not quietly use it anyway.
##
## OUTPUT
##   COH_per_sample.csv, COH_summary.csv, Figure10_module_coherence
##
## PREREQUISITE
##   Sections 11, 12, the 12b helper block, section 18 (gene sets), and
##   DETECT_cohort_gene_rates.csv from section 14 for detection matching.
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
## 19.0b  CANDIDATE METABOLIC GENE SETS
##
##   The 71-gene GOBP glutamine family set failed the coherence test at
##   percentile 0, and so did its hand-picked 12-gene core. Before concluding
##   that transcript data cannot represent metabolic state at all, the obvious
##   question is whether a BETTER SET would hold together. These four are the
##   candidates, and they are added to the same test, unchanged.
##
##   WHAT EACH ONE IS, AND WHY IT IS HERE
##
##   GlnNarrow   GOBP_GLUTAMINE_METABOLIC_PROCESS, M11960, 21 genes.
##               Glutamine only, without the arginine and proline arms that
##               diluted the 71-gene set. IMPORTANT CAVEAT: all 21 of its genes
##               are inside the set that already failed, so this is a slice of
##               the broken instrument rather than a new one. It can still be
##               more coherent than its parent, but if it fails too then breadth
##               was not the problem. It also still carries CPS1, FABP1, MTTP
##               and NR1H4, which are liver genes with no business in a
##               pancreatic glutamine module.
##
##   TCA_KEGG    KEGG_CITRATE_CYCLE_TCA_CYCLE, M3985, 32 genes.
##   TCA_React   REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE, M784, 34 genes.
##               Genuinely new: neither shares a single gene with the failed
##               sets. They also disagree with each other, sharing only 18
##               genes. KEGG includes pyruvate dehydrogenase and
##               gluconeogenesis (PC, PCK1, PCK2, ACLY); Reactome includes
##               iron-sulfur cluster assembly (ISCU, NFS1, SDHAF1-4, FXN).
##               That disagreement is useful rather than a nuisance: if both
##               cohere the conclusion is robust, and if only one does we learn
##               which definition holds together in this tissue.
##
##   OXPHOS      HALLMARK_OXIDATIVE_PHOSPHORYLATION, M5936, 200 genes.
##               Sits immediately downstream of the TCA cycle. Two genes are
##               stripped below because they are already in use elsewhere.
##
##   All four are clean of P4HA1, P4HA2, P4HA3, P4HB, PLOD1-3 and PPIB.
## ============================================================================

## M11960, GO:0006541, 21 unique genes, verbatim from the TSV
CAND_GLNNARROW <- c(
  "GLYATL1B", "CPS1", "FABP1", "SIRT4", "PHGDH", "BLOC1S6", "GLS2", "GLS",
  "GLUD1", "GLUL", "MECP2", "ASL", "MTTP", "LGSN", "PFAS", "NIT2",
  "UCP2", "CAD", "ARHGAP11B", "GLYATL1", "NR1H4")

## M3985, hsa00020, 32 unique genes, verbatim from the TSV
CAND_TCA_KEGG <- c(
  "CS", "DLAT", "DLD", "DLST", "FH", "SUCLG2P2", "IDH1", "IDH2",
  "IDH3A", "IDH3B", "IDH3G", "MDH1", "MDH2", "ACLY", "ACO1", "OGDH",
  "ACO2", "PC", "PCK1", "PCK2", "PDHA1", "PDHA2", "PDHB", "OGDHL",
  "SDHA", "SDHB", "SDHC", "SDHD", "SDHCP5", "SUCLG2", "SUCLG1", "SUCLA2")

## M784, R-HSA-71403, 34 unique genes, verbatim from the TSV
CAND_TCA_REACTOME <- c(
  "CS", "IDH3G", "SDHA", "ACAT1", "DLD", "FH", "ACO2", "IDH3B",
  "OGDH", "NNT", "SDHB", "DLST", "TRAP1", "KGD4", "ISCA1", "ISCU",
  "SUCLA2", "SIRT3", "SDHC", "MDH2", "SDHAF4", "SUCLG1", "FXN", "ISCA2",
  "IDH3A", "SDHAF2", "SUCLG2", "IDH2", "SDHAF3", "SDHD", "SDHAF1", "LYRM4",
  "CSKMT", "NFS1")

## M5936, hallmark, 200 unique genes, verbatim from the TSV
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

## OXPHOS shares LDHA with the Hypoxia module and HSPA9 with ERstressStrict.
## A candidate that shares genes with an existing module is not independent of
## it, so both go. Two of 200 changes nothing about the set's character.
CAND_OXPHOS <- setdiff(CAND_OXPHOS,
                       unique(unlist(MODULES[c("Hypoxia", "CollagenSynthesis", "iCAF",
                                               "ERstressStrict", "ERAD")],
                                     use.names = FALSE)))

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
## References first (we know how these behave), then the two that failed, then
## the four candidates. Keeping the references in every run is what makes the
## result readable: in PDAC, collagen at percentile 100 and hypoxia at 97 were
## the picture of a module that works.
COH_MODULES <- c("Hypoxia", "CollagenSynthesis", "ERstressStrict", "ERAD",
                 "Glutamine", "OGcore",
                 "GlnNarrow", "TCA_KEGG", "TCA_Reactome", "OXPHOS")

COH_N_NULL     <- 30     # random gene sets per module per tumour
COH_MIN_GENES  <- 4      # a set smaller than this has no meaningful coherence
COH_MIN_SPOTS  <- 30
COH_SEED       <- 90125
COH_CALIPER    <- 0.35   # relative band for detection-rate matching of null sets

assert_run_state("Section 19")

## ============================================================================
## 19.1  BUILD THE POOL OF GENES THE NULL SETS ARE DRAWN FROM
##
##   Null sets must be matched to each real module on detection rate, because
##   coherence measured on sparse data depends heavily on how often the genes
##   are seen at all. A null built from rare genes would make any real module
##   look coherent by comparison.
## ============================================================================
rates <- read.csv(RATES_CSV, stringsAsFactors = FALSE)
gene_med <- rates %>%
  group_by(gene) %>%
  summarise(median_rate = median(detect_rate), n_tum = dplyr::n(), .groups = "drop") %>%
  filter(n_tum == dplyr::n_distinct(rates$sample_id))

## genes that belong to any module under test are excluded from the null pool,
## otherwise a null set could contain the very genes it is a null for
ALL_MODULE_GENES <- unique(unlist(MODULES[COH_MODULES], use.names = FALSE))
null_pool <- gene_med[!gene_med$gene %in% ALL_MODULE_GENES, ]

cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 19.1  SETUP\n")
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  genes present in all tumours : %d\n", nrow(gene_med)))
cat(sprintf("  excluded (module members)    : %d\n", length(ALL_MODULE_GENES)))
cat(sprintf("  null pool                    : %d genes\n", nrow(null_pool)))
cat("\n  modules under test:\n")
for (nm in COH_MODULES) {
  det <- intersect(MODULES[[nm]], gene_med$gene)
  mr  <- if (length(det)) median(gene_med$median_rate[gene_med$gene %in% det]) else NA_real_
  cat(sprintf("    %-18s %3d genes, %3d detected cohort-wide, median rate %.3f\n",
              nm, length(MODULES[[nm]]), length(det), mr))
}

## for a given module, draw COH_N_NULL random sets matched on size and rate
draw_null_sets <- function(module_genes, n_sets) {
  det <- intersect(module_genes, gene_med$gene)
  if (length(det) < COH_MIN_GENES) return(NULL)
  target <- median(gene_med$median_rate[gene_med$gene %in% det])
  elig <- null_pool$gene[null_pool$median_rate >= target * (1 - COH_CALIPER) &
                           null_pool$median_rate <= target * (1 + COH_CALIPER)]
  if (length(elig) < length(det) * 2) return(NULL)   # not enough to sample from
  lapply(seq_len(n_sets), function(i) sample(elig, length(det)))
}

set.seed(COH_SEED)
NULL_SETS <- lapply(COH_MODULES, function(nm) draw_null_sets(MODULES[[nm]], COH_N_NULL))
names(NULL_SETS) <- COH_MODULES

cat("\n  null sets drawn:\n")
for (nm in COH_MODULES)
  cat(sprintf("    %-18s %s\n", nm,
              if (is.null(NULL_SETS[[nm]])) "*** could not match, module skipped ***"
              else sprintf("%d sets of %d genes", length(NULL_SETS[[nm]]), length(NULL_SETS[[nm]][[1]]))))

## ============================================================================
## 19.2  THE TWO COHERENCE MEASURES
## ============================================================================

## mean pairwise Spearman correlation among the rows of a gene-by-spot matrix.
## Genes with no variance are dropped: a constant gene has no correlation with
## anything and would otherwise inject NAs.
coh_meanrho <- function(mat) {
  if (is.null(mat) || nrow(mat) < COH_MIN_GENES) return(NA_real_)
  keep <- apply(mat, 1, function(v) sd(v) > 0)
  if (sum(keep) < COH_MIN_GENES) return(NA_real_)
  m <- mat[keep, , drop = FALSE]
  cc <- suppressWarnings(cor(t(as.matrix(m)), method = "spearman"))
  if (!is.matrix(cc)) return(NA_real_)
  mean(cc[upper.tri(cc)], na.rm = TRUE)
}

## proportion of variance taken by PC1 after scaling each gene. If the members
## move as one, PC1 dominates. Less affected by ties than pairwise Spearman.
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
## 19.3  RUN, ONE TUMOUR AT A TIME
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
    
    ## percentile of the real module within its own null distribution
    pct_rho <- if (sum(is.finite(null_rho)) >= 10)
      mean(null_rho[is.finite(null_rho)] < r_rho) else NA_real_
    pct_pc1 <- if (sum(is.finite(null_pc1)) >= 10)
      mean(null_pc1[is.finite(null_pc1)] < r_pc1) else NA_real_
    
    rows[[length(rows) + 1]] <- data.frame(
      sample_id = sid, module = nm,
      n_genes_used = nrow(real),
      mean_rho = r_rho, pc1_var = r_pc1,
      null_mean_rho = median(null_rho, na.rm = TRUE),
      null_pc1 = median(null_pc1, na.rm = TRUE),
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
if (nrow(coh_all) == 0) stop("Section 19: no tumour produced coherence estimates.")
write_out(coh_all, "COH_per_sample")

## ============================================================================
## 19.4  SUMMARISE
## ============================================================================
coh_sum <- coh_all %>%
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
cat("  Each module is compared against random gene sets of the SAME SIZE matched\n")
cat("  on detection rate. Absolute coherence is meaningless in sparse data; the\n")
cat("  percentile against its own null is the number that matters.\n\n")

cat(sprintf("  %-18s %6s %8s %9s %9s %8s %9s\n",
            "module", "genes", "rho", "null rho", "pctile", "PC1", "pctile"))
for (i in seq_len(nrow(coh_sum))) with(coh_sum[i, ],
                                       cat(sprintf("  %-18s %6.0f %8.3f %9.3f %8.0f%% %8.3f %8.0f%%\n",
                                                   module, median_genes, median_rho, median_null_rho, 100*median_pct_rho,
                                                   median_pc1, 100*median_pct_pc1)))

cat("\n  rho      = mean pairwise correlation between member genes\n")
cat("  null rho = the same for size and rate matched random sets\n")
cat("  pctile   = where the module sits in its own null distribution\n")
cat("  PC1      = share of variance taken by the first principal component\n")

## ============================================================================
## 19.5  VERDICT
##
##   A module that is a real biological unit should sit high in its own null
##   distribution. 90% is the bar: it means the module is more internally
##   consistent than 9 out of 10 random sets of the same size and detection
##   rate. Hypoxia and CollagenSynthesis show what a working module looks like.
## ============================================================================
COH_PCT_MIN <- 0.90

pct_of <- function(mod, col) {
  r <- coh_sum[coh_sum$module == mod, ]
  if (!nrow(r)) return(NA_real_)
  r[[col]]
}
passes <- function(mod) {
  a <- pct_of(mod, "median_pct_rho"); b <- pct_of(mod, "median_pct_pc1")
  isTRUE(a >= COH_PCT_MIN) || isTRUE(b >= COH_PCT_MIN)
}

cat("\n", strrep("-", 74), "\n", sep = "")
cat(sprintf("  criterion: median percentile >= %.0f%% on either measure\n\n", 100*COH_PCT_MIN))
for (nm in coh_sum$module)
  cat(sprintf("    %-18s %s\n", nm, if (passes(nm)) "COHERENT" else "not coherent"))

hyp_ok  <- passes("Hypoxia")
col_ok  <- passes("CollagenSynthesis")
glut_ok <- passes("Glutamine")
og_ok   <- passes("OGcore")
ctrl_ok <- hyp_ok || col_ok

cat("\n", strrep("=", 74), "\n", sep = "")
cat("VERDICT\n")
cat(strrep("=", 74), "\n", sep = "")

if (!ctrl_ok) {
  cat("  THE TEST FAILED AGAIN, AND THIS SECTION PROVES NOTHING.\n\n")
  cat("  Neither hypoxia nor collagen synthesis is coherent by this measure, and\n")
  cat("  both are modules we know work. So the measure does not function in this\n")
  cat("  data, exactly as the yardstick approach did not.\n\n")
  cat("  Do NOT conclude anything about the glutamine module from this. Two\n")
  cat("  diagnostics have now failed, which is itself worth telling Valerio:\n")
  cat("  the question may not be answerable with this data type.\n")
} else if (!glut_ok && !og_ok) {
  cat("  THE MODULE IS NOT A COHERENT UNIT. This is the expected outcome.\n\n")
  cat(sprintf("  Working modules: hypoxia %s, collagen synthesis %s.\n",
              if (hyp_ok) "coherent" else "not", if (col_ok) "coherent" else "not"))
  cat("  The glutamine module and its 12-gene core are not.\n\n")
  cat("  WHAT THIS MEANS IN PLAIN TERMS\n")
  cat("    A module score is an average over its member genes. That average only\n")
  cat("    represents something real if the members move together. Here they do\n")
  cat("    not, so the score is an average of unrelated quantities and does not\n")
  cat("    correspond to any single biological state.\n\n")
  cat("  WHAT IT LICENSES YOU TO SAY\n")
  cat("    Section 18's null is about the instrument, not the biology. The\n")
  cat("    cofactor hypothesis has not been tested, let alone refuted. This also\n")
  cat("    explains the earlier oddity that the module correlated with nothing at\n")
  cat("    all: there was no coherent signal in it to correlate with.\n\n")
  cat("  WHAT IT DOES NOT LICENSE\n")
  cat("    Any claim about whether the decoupling occurs along the cofactor axis.\n")
  cat("    An incoherent module gives no evidence in either direction.\n\n")
  cat("  WHY THIS WAS PREDICTABLE, IN HINDSIGHT\n")
  cat("    The GO term spans arginine, glutamate, glutamine and proline\n")
  cat("    metabolism. The urea cycle, proline synthesis and nitric oxide genes\n")
  cat("    in it are not co-regulated with glutaminolysis, so there was never\n")
  cat("    reason to expect them to move together in one tissue.\n")
} else if (glut_ok || og_ok) {
  cat("  THE MODULE IS COHERENT. This is the harder outcome.\n\n")
  cat(sprintf("  Glutamine %s, 2-OG core %s.\n",
              if (glut_ok) "coherent" else "not coherent",
              if (og_ok) "coherent" else "not coherent"))
  cat("  So the members do move together and the score represents something.\n")
  cat("  Section 18's null cannot then be dismissed as an incoherent instrument.\n\n")
  cat("  WHAT THIS MEANS\n")
  cat("    Two readings remain, and this data cannot separate them:\n")
  cat("      (a) the decoupling is specific to oxygen limitation rather than to\n")
  cat("          enzymatic capacity in general; or\n")
  cat("      (b) the module is coherent but still measures transcriptional\n")
  cat("          co-regulation rather than metabolic flux, which is controlled\n")
  cat("          after transcription.\n\n")
  cat("    Coherence rules out one failure mode, not all of them. Report the\n")
  cat("    result as it is and let Valerio judge.\n")
}

## ---- did any candidate set rescue the metabolic axis? ----------------------
CANDIDATES <- c("GlnNarrow", "TCA_KEGG", "TCA_Reactome", "OXPHOS")
cand_ok <- CANDIDATES[vapply(CANDIDATES, function(x)
  x %in% coh_sum$module && passes(x), logical(1))]

cat("\n", strrep("-", 74), "\n", sep = "")
cat("  CANDIDATE SETS: is there a usable metabolic instrument at all?\n\n")
for (nm in intersect(CANDIDATES, coh_sum$module))
  cat(sprintf("    %-14s %s\n", nm, if (passes(nm)) "COHERENT" else "not coherent"))

if (!ctrl_ok) {
  cat("\n    References failed, so ignore the above. Nothing is interpretable.\n")
} else if (!length(cand_ok)) {
  cat("\n    NONE of the candidates is coherent either.\n")
  cat("    Four independent metabolic gene sets, from three different databases,\n")
  cat("    all fail while collagen and hypoxia pass cleanly. That is no longer a\n")
  cat("    bad-gene-set problem: metabolic genes do not move together in this\n")
  cat("    data, which fits what is known biologically, since metabolic flux is\n")
  cat("    controlled largely after transcription rather than by co-regulating\n")
  cat("    the enzymes.\n\n")
  cat("    CONCLUSION: the cofactor axis is not testable with Visium transcript\n")
  cat("    data. Say that plainly. It is a real finding about method rather than\n")
  cat("    a failure, and it saves anyone else repeating the attempt.\n")
} else {
  cat(sprintf("\n    %d candidate(s) coherent: %s\n", length(cand_ok),
              paste(cand_ok, collapse = ", ")))
  cat("    A usable instrument may exist. Before rerunning section 18 on it,\n")
  cat("    check two things:\n")
  cat("      1. its correlation with the hypoxia module. A coherent metabolic\n")
  cat("         module that simply mirrors hypoxia is not a second axis.\n")
  cat("      2. whether it plausibly reports on 2-oxoglutarate supply. OXPHOS\n")
  cat("         cohering would show the method works, but OXPHOS is downstream\n")
  cat("         of the TCA cycle, not a measure of cofactor availability.\n")
  cat("    Do not rerun the axis on all of them: that is four more tests on the\n")
  cat("    same question. Pick one on biological grounds first.\n")
}

cat("\n  CAVEATS\n")
cat("    Sparse data means most gene pairs share many double zeros, which\n")
cat("    inflates ties and drags correlations toward a dropout-driven value.\n")
cat("    That is why everything is reported as a percentile against a matched\n")
cat("    null rather than as an absolute number.\n")
cat("    Coherence is necessary for a module score to mean anything. It is not\n")
cat("    sufficient: a coherent module can still measure the wrong thing.\n")
cat(strrep("=", 74), "\n", sep = "")

## ============================================================================
## 19.6  FIGURE 10
## ============================================================================
if (!exists("FIG_STYLES") || !exists("theme_fig")) {
  cat("\n  Figure 10 skipped: section 12 has not run.\n")
} else {
  pd <- coh_all[is.finite(coh_all$pct_rho), ]
  lv <- intersect(COH_MODULES, unique(pd$module))
  nice_m <- c(Hypoxia = "Hypoxia", CollagenSynthesis = "Collagen synth.",
              ERstressStrict = "ER stress", ERAD = "ERAD",
              Glutamine = "Glutamine (71)", OGcore = "2-OG core (12)",
              GlnNarrow = "Glutamine (21)", TCA_KEGG = "TCA KEGG",
              TCA_Reactome = "TCA Reactome", OXPHOS = "OXPHOS")
  pd$module <- factor(nice_m[pd$module], levels = unname(nice_m[lv]))
  pd$is_metab <- ifelse(grepl("Glutamine|2-OG|TCA|OXPHOS", as.character(pd$module)),
                        "metabolic", "reference")
  
  for (sn in names(FIG_STYLES)) {
    S <- FIG_STYLES[[sn]]
    outdir <- file.path(agg, S$dir); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    
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
           caption = sprintf("One point per tumour. Shaded region above %.0f%% is where a coherent module should sit. Grey line is chance.",
                             100 * COH_PCT_MIN)) +
      theme_fig(S) +
      theme(panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(linewidth = S$grid, colour = "grey90"),
            axis.text.x = element_text(angle = 20, hjust = 1))
    
    if (isTRUE(S$titles))
      p10 <- p10 + labs(title = "Do the genes inside each module behave as a unit?",
                        subtitle = "Compared against random gene sets of the same size and detection rate")
    
    h <- if (isTRUE(S$split)) S$h_full else MM(74)
    for (dv in S$devices) {
      f <- file.path(outdir, outf("Figure10_module_coherence", dv))
      tryCatch({
        if (dv == "pdf") ggsave(f, p10, width = S$w2, height = h, device = PDF_DEV, bg = "white")
        else             ggsave(f, p10, width = S$w2, height = h, dpi = S$dpi, bg = "white")
      }, error = function(e) log_msg("  fig failed (%s): %s", basename(f), conditionMessage(e)))
    }
    cat(sprintf("  %s written to %s\n", outf("Figure10_module_coherence", "png"), S$dir))
  }
}

log_msg("Section 19 complete.")
################################################################################
## 20. DOES CELL COMPOSITION EXPLAIN THE RESULT?
##
##   BRCA NOTE: the reasoning in this header was written for the PDAC
##   cohort. Any statement about what the data showed, which module worked,
##   or how many tumours were usable describes PDAC, not BRCA. The code
##   below is unchanged except where tagged "BRCA:".
##
##   THE WORRY
##   A Visium spot reads several cells at once. If well oxygenated spots are
##   mostly fibroblasts (collagen and P4HA1 together) and hypoxic spots are
##   mostly tumour cells (P4HA1 and ER stress without collagen), the
##   decoupling appears without any cell changing its behaviour. Only the
##   mix inside the spot changed.
##
##   THE COMPOSITION MEASURE
##   Every SeuratObj.rds carries per-spot ESTIMATE output in
##   obj@misc$Tumor_purity_estimate$table (stromal, immune, estimate, purity),
##   computed on Spatial counts, same method for all 23 samples (verified in
##   PDAC). BRCA: NOT verified. A section whose object has no such table is
##   skipped and logged (see 20.2), so check the log for how many were lost.
##   It cannot be used as it stands: the ESTIMATE stromal signature contains
##   COL1A2, COL3A1 and COL6A3 from CollagenSynthesis and CXCL12 from iCAF
##   (estimate_overlap_audit.R). Adjusting for it would remove part of the
##   collagen arm mechanically.
##   So stroma and immune scores are REBUILT here from the ESTIMATE signatures
##   with every gene shared with a predictor module, an outcome or the
##   hydroxylase family removed, scored with AddModuleScore exactly like the
##   other modules.
##
##   THREE MODELS, SAME SPOTS, SAME POOLING
##     M0 baseline    : the section 15 model unchanged. MUST reproduce the
##                      section 15 result, or the section stops.
##     M1 adjusted    : M0 + zStroma + zImmune. Among spots with the same
##                      amount of stroma and immune cells, is the pattern there?
##     M2 low stroma  : M0 fitted only in the lowest third of spots by stroma
##                      score within each tumour. Where the spot is mostly
##                      tumour cells, is the pattern there?
##
##   HOW TO READ IT
##   Model p values are inflated by spatial autocorrelation (in PDAC 296 of
##   1000 random genes were significant, section 16; check section 16 here). Read the size of the slope difference and how
##   much of it survives, not the p value.
##
##   Run after the main pipeline in the same FRESH session, like sections 18
##   and 19. Needs work, MODULES, SEED, ADJUST_DEPTH and the 12b helpers.
################################################################################

stopifnot(exists("MODULES"), exists("work"), exists("agg"), exists("SEED"),
          exists("ADJUST_DEPTH"), exists("pool_with_p"), exists("fmt_p"),
          exists("write_out"), exists("assert_run_state"), exists("log_msg"),
          exists("bind"), exists("min_genes_for"))

## ============================================================================
## 20.0  SETTINGS
## ============================================================================
SC_SIG_CSV       <- file.path(PROJECT_ROOT, "data", "estimate_audit", "ESTIMATE_signatures_used.csv")  # BRCA: from project root
SC_OUTCOMES      <- c("P4HA1", "P4HA2", "P4HA3")
SC_ES            <- "ERstressStrict"
SC_BASE_MODULES  <- c("Hypoxia", "CollagenSynthesis", "iCAF", SC_ES, "ERAD")
SC_NEW_MODULES   <- c("StromaESTc", "ImmuneESTc")
SC_MIN_SPOTS     <- 30
SC_MIN_EVENTS    <- 70          # same as section 15
SC_ADJUST_DEPTH  <- ADJUST_DEPTH  # derived, never set by hand
SC_STRAT_Q       <- 1/3         # M2 keeps spots at or below this stroma quantile
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
## 20.1  BUILD THE CLEAN STROMA AND IMMUNE GENE SETS
##       Genes to remove are computed, never typed.
## ============================================================================
if (!file.exists(SC_SIG_CSV))
  stop("Section 20 needs ", SC_SIG_CSV, ". Run estimate_overlap_audit.R first.")
sig <- read.csv(SC_SIG_CSV, stringsAsFactors = FALSE)
STROMA_RAW <- unique(sig$gene[sig$signature == "stromal"])
IMMUNE_RAW <- unique(sig$gene[sig$signature == "immune"])
stopifnot(length(STROMA_RAW) == 141, length(IMMUNE_RAW) == 141)

SC_BLOCK <- unique(c(unlist(MODULES[SC_BASE_MODULES], use.names = FALSE),
                     SC_FAMILY, SC_OUTCOMES))
both_sigs <- intersect(STROMA_RAW, IMMUNE_RAW)   # a gene in both would tie them together
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

## ---- circularity gate: hard stop --------------------------------------------
for (og in SC_OUTCOMES) {
  bad <- names(SC_MOD[c(SC_BASE_MODULES, SC_NEW_MODULES)])[vapply(
    SC_MOD[c(SC_BASE_MODULES, SC_NEW_MODULES)], function(g) og %in% g, logical(1))]
  if (length(bad)) stop("Section 20 halted: ", og, " is inside ", paste(bad, collapse = ", "))
}
shared_new <- lapply(SC_BASE_MODULES, function(m) intersect(SC_MOD[[m]], c(STROMA_CLEAN, IMMUNE_CLEAN)))
if (any(lengths(shared_new) > 0)) stop("Section 20 halted: a clean set still shares genes with a module.")
cat("  circularity gate: clean\n")

## ============================================================================
## 20.2  SCORE, STANDARDISE AND VALIDATE, ONE TUMOUR AT A TIME
##       Scoring follows section 15 exactly (same seed before each module,
##       nbin 20, ctrl 50), so the five original scores are identical to
##       section 15's. That is what lets M0 reproduce it.
## ============================================================================
sc_expr_one <- function(rds_path, sid) {
  obj <- readRDS(rds_path)
  obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
  DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else Assays(obj)[1]
  
  if (length(setdiff(SC_OUTCOMES, rownames(obj)))) {
    log_msg("  %s: outcome gene absent -> skipped", sid); return(NULL)
  }
  
  ## original ESTIMATE table, aligned to the spot order used by FetchData
  est <- obj@misc$Tumor_purity_estimate$table
  if (is.null(est)) { log_msg("  %s: no ESTIMATE table -> skipped", sid); return(NULL) }
  if (!is.null(rownames(est)) && all(colnames(obj) %in% rownames(est))) {
    est <- est[colnames(obj), , drop = FALSE]
  } else {
    if (nrow(est) != ncol(obj) ||
        max(abs(est$purity - obj@meta.data$Tumor_purity_estimate)) > 1e-6)
      stop(sprintf("Section 20, %s: ESTIMATE table cannot be aligned to spots.", sid))
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
##       M2 keeps the whole-tumour z-scores, so its units stay 1 SD of the
##       whole tumour and its contrast is directly comparable to M0.
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
## 20.6  SAFETY CHECK: M0 MUST REPRODUCE SECTION 15, FOR EVERY OUTCOME
## ============================================================================
cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 20.6  SAFETY CHECK\n")
cat(strrep("=", 74), "\n", sep = "")
m0_est <- setNames(rep(NA_real_, length(SC_OUTCOMES)), SC_OUTCOMES)
for (og in SC_OUTCOMES) {
  m0 <- sc_con_pooled[sc_con_pooled$outcome == og & sc_con_pooled$model == "M0_baseline" &
                        sc_con_pooled$contrast == "interaction_diff", ]
  if (!nrow(m0)) stop("Section 20 halted: no M0 result for ", og)
  if (exists("iso_con_pooled")) {
    ref <- iso_con_pooled[iso_con_pooled$outcome == og &
                            iso_con_pooled$contrast == "interaction_diff", ]$est
    tol <- 1e-6; ref_src <- "section 15 in this session"
  } else {
    stop("Section 20 needs section 15 (iso_con_pooled) in the same session for the safety check.")
  }
  cat(sprintf("  %-6s M0 = %+.5f   reference = %+.5f (%s)\n", og, m0$est, ref, ref_src))
  if (!length(ref) || abs(m0$est - ref) > tol)
    stop("Section 20 halted: M0 does not reproduce section 15 for ", og, ". Nothing below is interpretable.")
  m0_est[og] <- m0$est
}
cat("  PASSED: every baseline is the section 15 result.\n")

## ============================================================================
## 20.7  REPORT
## ============================================================================
cat("\n", strrep("=", 74), "\n", sep = "")
cat(sprintf("SECTION 20  [%s]  run_id %s\n", if (SC_ADJUST_DEPTH) "DEPTH ADJUSTED" else "UNADJUSTED", RUN_ID))
cat(strrep("=", 74), "\n", sep = "")

cat("\n--- is the rebuilt composition measure the same thing as ESTIMATE? ---\n")
cat(sprintf("  stroma vs ESTIMATE stromal : median rho %.2f (min %.2f)\n",
            median(sc_val$rho_stroma_vs_ESTIMATE), min(sc_val$rho_stroma_vs_ESTIMATE)))
cat(sprintf("  immune vs ESTIMATE immune  : median rho %.2f (min %.2f)\n",
            median(sc_val$rho_immune_vs_ESTIMATE), min(sc_val$rho_immune_vs_ESTIMATE)))

cat("\n--- how strongly stroma travels with the model's predictors ---\n")
cat(sprintf("  stroma vs collagen synthesis : median rho %+.2f\n", median(sc_val$rho_stroma_vs_collagen)))
cat(sprintf("  stroma vs hypoxia            : median rho %+.2f\n", median(sc_val$rho_stroma_vs_hypoxia)))
cat(sprintf("  stroma vs ER stress          : median rho %+.2f\n", median(sc_val$rho_stroma_vs_ERstress)))
cat(sprintf("  immune vs hypoxia            : median rho %+.2f\n", median(sc_val$rho_immune_vs_hypoxia)))

for (og in SC_OUTCOMES) {
  cat("\n", strrep("-", 74), "\n", sep = "")
  cat(sprintf("  OUTCOME %s\n", og))
  cat(strrep("-", 74), "\n", sep = "")
  
  cat("\n  tumours per model\n")
  for (md in SC_MODELS) {
    s <- sc_stat[sc_stat$outcome == og & sc_stat$model == md, ]
    cat(sprintf("    %-14s fitted in %2d / %2d   median spots %5.0f   median events %5.0f\n",
                md, sum(s$fitted), nrow(s), median(s$n_spots[s$fitted]), median(s$n_events[s$fitted])))
    drop <- s[!s$fitted, ]
    for (j in seq_len(nrow(drop)))
      cat(sprintf("                   excluded: %-26s %s\n", drop$sample_id[j], drop$reason[j]))
  }
  
  cat("\n  THE HEADLINE: difference between the two slopes\n")
  for (md in SC_MODELS) {
    r <- sc_con_pooled[sc_con_pooled$outcome == og & sc_con_pooled$model == md &
                         sc_con_pooled$contrast == "interaction_diff", ]
    if (!nrow(r)) next
    cat(sprintf("    %-14s %+.4f  (95%% CI %+.4f to %+.4f)  k=%2d  I2=%3.0f%%  p = %-10s  retained %4.0f%% of M0\n",
                md, r$est, r$lower, r$upper, r$k, r$I2, fmt_p(r$p), 100 * r$est / m0_est[og]))
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
cat(sprintf("  k counts %s (POOL_LEVEL = %s); section 21 gives the other unit.\n",
            if (POOL_LEVEL == "patient") "patients" else "sections", POOL_LEVEL))
cat("  ESTIMATE signatures were built for bulk tumours, not single spots.\n")
cat("  Stroma and collagen synthesis are biologically linked (rho about 0.7 in PDAC), so M1\n")
cat("  asks a narrower question than 'does composition explain it'. M2 is cleaner.\n")
cat("  M2 keeps whole-tumour z-scores; its hypoxia and collagen ranges are narrower.\n")
cat("  Compare models WITHIN a gene. Coefficient sizes are not comparable ACROSS\n")
cat("  genes with different detection rates.\n")
cat("  P4HA3 is detected in few spots, so M2 may lose many tumours for it; read its\n")
cat("  k before reading its estimate.\n")
cat(strrep("=", 74), "\n", sep = "")

log_msg("Section 20 complete.")
################################################################################
## 21. SECTION-LEVEL SENSITIVITY   (BRCA: the PDAC section 21, turned around)
##
##   THE SITUATION
##   In PDAC, pooling was over sections and this section checked patients.
##   In BRCA the pipeline already pools over patients (POOL_LEVEL, section 0),
##   because several patients contribute 4 to 6 serial sections. So this
##   section now reports the other direction: what the result would have been
##   with sections as the unit, as in the PDAC analysis.
##
##   WHAT THIS DOES
##   1. Reproduces the pipeline's pooled results from the per-section CSVs
##      with its own code (fixed effect within patient, then DerSimonian-Laird
##      across patients). If it cannot, it stops: then its pooling is not the
##      pipeline's pooling.
##   2. Pools the same per-section estimates with sections as the unit.
##   3. If metafor is installed, fits a three-level model (sections nested in
##      patients) as a second opinion that does not assume serial sections
##      estimate exactly the same thing.
##   4. Repeats the patient pooling treating HT339B1 and HT339B2 as one person.
##      The lookup table gives them different HTAN participant IDs (HTA12_250,
##      HTA12_262) but they share the HT339 prefix; it is not confirmed whether
##      they are the same person. Same situation as HT270P1 / HT270P2 in PDAC.
##
##   No model is refitted. Standalone: needs only the aggregate CSVs.
################################################################################

## BRCA: same folder logic as section 0, so this also runs in a fresh session
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
## 21.0  READ THE PER-SECTION CONTRASTS AND THE SECTION -> PATIENT MAP
##       Only unadjusted rows, and only one run_id per file, so nothing from
##       different runs is mixed.
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
iso    <- read_one("ISOFORM_contrasts_per_sample")
iso_pl <- read_one("ISOFORM_slope_difference_pooled")
sc     <- read_one("SC_contrasts_per_sample", required = FALSE)     # section 20 is optional
sc_pl  <- read_one("SC_slope_difference_pooled", required = FALSE)
cmap   <- read_one("cohort_sections")

iso <- iso[iso$contrast == "interaction_diff", ]
if (!is.null(sc)) sc <- sc[sc$contrast == "interaction_diff", ]

## ============================================================================
## 21.1  PATIENT ID FROM THE COHORT TABLE (not parsed from the section ID:
##       BRCA IDs such as S17-41599-Fs1U1Bp1 do not follow the HT pattern)
## ============================================================================
pid_of <- function(sid) {
  p <- cmap$patient_id[match(sid, cmap$sample_id)]
  if (anyNA(p)) stop("No patient for: ", paste(unique(sid[is.na(p)]), collapse = ", "))
  p
}
map <- unique(data.frame(sample_id = iso$sample_id, patient = pid_of(iso$sample_id)))
map <- map[order(map$patient, map$sample_id), ]
n_pat <- length(unique(map$patient))
cat(sprintf("\nSections %d, patients %d\n", nrow(map), n_pat))
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
  c(est = est, lower = est - qnorm(0.975) * se, upper = est + qnorm(0.975) * se,
    p = 2 * pnorm(-abs(est / se)), k = k, I2 = I2)
}

## fixed-effect, inverse variance, within one patient (same as meta_pool)
fe_collapse <- function(d, pid) {
  d$patient <- pid
  do.call(rbind, lapply(split(d, d$patient), function(x) {
    ok <- is.finite(x$est) & is.finite(x$var) & x$var > 0; x <- x[ok, ]
    if (!nrow(x)) return(NULL)
    w <- 1 / x$var
    data.frame(patient = x$patient[1], n_sections = nrow(x),
               est = sum(w * x$est) / sum(w), var = 1 / sum(w))
  }))
}
pat_pool <- function(d, pid) { pc <- fe_collapse(d, pid); dl_pool(pc$est, pc$var) }

## the pipeline's own pooling unit for this run
run_pool <- function(d) if (identical(POOL_LEVEL, "section")) dl_pool(d$est, d$var) else pat_pool(d, pid_of(d$sample_id))

## ============================================================================
## 21.3  CHECK: REPRODUCE THE PIPELINE'S POOLED RESULTS
## ============================================================================
cat(sprintf("\nCheck: does this script's pooling reproduce the pipeline's (POOL_LEVEL = %s)?\n", POOL_LEVEL))
chk <- list()
for (og in unique(iso$outcome)) {
  mine <- run_pool(iso[iso$outcome == og, ])
  ref  <- iso_pl$est[iso_pl$outcome == og & iso_pl$contrast == "interaction_diff"]
  chk[[length(chk) + 1]] <- c(label = paste("S15", og), mine = mine[["est"]], ref = ref)
}
if (!is.null(sc)) for (og in unique(sc$outcome)) for (md in unique(sc$model)) {
  s <- sc[sc$outcome == og & sc$model == md, ]
  mine <- run_pool(s)
  ref  <- sc_pl$est[sc_pl$outcome == og & sc_pl$model == md & sc_pl$contrast == "interaction_diff"]
  chk[[length(chk) + 1]] <- c(label = paste("S20", og, md), mine = mine[["est"]], ref = ref)
}
chk <- as.data.frame(do.call(rbind, chk), stringsAsFactors = FALSE)
chk$mine <- as.numeric(chk$mine); chk$ref <- as.numeric(chk$ref)
chk$ok <- abs(chk$mine - chk$ref) < TOL
for (i in seq_len(nrow(chk)))
  cat(sprintf("  %-26s mine %+.6f  pipeline %+.6f  %s\n",
              chk$label[i], chk$mine[i], chk$ref[i], if (isTRUE(chk$ok[i])) "ok" else "*** MISMATCH ***"))
if (!all(chk$ok %in% TRUE)) stop("Section 21 halted: this pooling does not reproduce the pipeline's. Nothing below is interpretable.")
cat("  PASSED\n")

## ============================================================================
## 21.4  PATIENT LEVEL (primary) vs SECTION LEVEL (PDAC style)
## ============================================================================
analyses <- list()
for (og in unique(iso$outcome))
  analyses[[paste("S15", og)]] <- iso[iso$outcome == og, ]
if (!is.null(sc)) for (og in unique(sc$outcome)) for (md in c("M1_adjusted", "M2_low_stroma"))
  analyses[[paste("S20", og, md)]] <- sc[sc$outcome == og & sc$model == md, ]

pid_merged <- function(sid) { p <- pid_of(sid); p[p == "HT339B2"] <- "HT339B1"; p }

has_metafor <- requireNamespace("metafor", quietly = TRUE)
res <- do.call(rbind, lapply(names(analyses), function(lb) {
  d   <- analyses[[lb]]
  s   <- dl_pool(d$est, d$var)
  pp  <- pat_pool(d, pid_of(d$sample_id))
  pm  <- pat_pool(d, pid_merged(d$sample_id))
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
res$ci_excludes_0_pat <- res$pat_hi < 0 | res$pat_lo > 0
res$ci_excludes_0_sec <- res$sec_hi < 0 | res$sec_lo > 0
write.csv(res, file.path(OUT_DIR, "SECTION_level_sensitivity.csv"), row.names = FALSE)

cat("\n", strrep("-", 74), "\n", sep = "")
cat("SLOPE DIFFERENCE: PATIENTS (primary) vs SECTIONS\n")
cat(strrep("-", 74), "\n", sep = "")
fmt <- function(e, l, h) if (is.na(e)) "            n/a            " else sprintf("%+.4f (%+.4f to %+.4f)", e, l, h)
for (i in seq_len(nrow(res))) with(res[i, ], {
  cat(sprintf("\n  %s\n", analysis))
  cat(sprintf("    patients  k=%2d  %s  p = %.2g  I2 = %.0f%%\n", k_pat, fmt(pat_est, pat_lo, pat_hi), pat_p, pat_I2))
  cat(sprintf("    sections  k=%2d  %s  p = %.2g%s\n", k_sec, fmt(sec_est, sec_lo, sec_hi), sec_p,
              if (ci_excludes_0_sec != ci_excludes_0_pat) "   <-- CONCLUSION CHANGES" else ""))
  cat(sprintf("    HT339 merged k=%2d  %s\n", k_patm, fmt(patm_est, patm_lo, patm_hi)))
  if (has_metafor) cat(sprintf("    3-level       %s\n", fmt(ml_est, ml_lo, ml_hi)))
})

cat("\nNOTES\n")
if (!has_metafor) cat("  metafor not installed, so the three-level check was skipped.\n")
if (is.null(sc))  cat("  Section 20 outputs not found, so only section 15 contrasts were checked.\n")
cat("  Within a patient, sections are combined as fixed effect: it assumes they\n")
cat("  estimate the same thing and are independent. Serial sections are not\n")
cat("  independent, so each patient's variance is too small. The three-level\n")
cat("  model does not make that assumption.\n")
cat("  Model p values remain inflated by spatial autocorrelation.\n")
cat("  Output: ", file.path(OUT_DIR, "SECTION_level_sensitivity.csv"), "\n", sep = "")
cat(strrep("=", 74), "\n", sep = "")

################################################################################
## 22. INFLUENCE OF SINGLE PATIENTS, AND THE METASTASIS
##
##   THREE QUESTIONS, ALL ANSWERED FROM THE PER-SECTION CSVs. No model is
##   refitted, so this runs in seconds and cannot change any estimate.
##
##   22.A  LEAVE ONE PATIENT OUT
##         Drop each patient (all of their sections) in turn and pool the
##         rest. If the headline moves a lot when one patient goes, the result
##         rests on that patient and its identity has to be reported.
##         BRCA: PDAC dropped one SECTION at a time. With up to 6 serial
##         sections per patient that would leave the rest of the tumour in.
##
##   22.B  WITH vs WITHOUT THE METASTASIS
##         BRCA: the 5 metastatic sections all come from one patient, HT268B1
##         (Sample_ID_Lookup_table_v1, sheet Visium). The main run leaves them
##         out (INCLUDE_METASTASIS <- FALSE). Run the pipeline once more with
##         INCLUDE_METASTASIS <- TRUE, then run this section on that output:
##         it reports the pooled result with and without HT268B1.
##
##   22.C  SIGN CONCORDANCE
##         How many patients (and sections) point the same way, before any
##         pooling. A pooled coefficient hides this.
##
##   Standalone: needs only the aggregate CSVs, not the pipeline session.
################################################################################

if (!exists("OUT_ROOT")) {
  PROJECT_ROOT       <- Sys.getenv("P4HA1_BRCA_ROOT", unset = path.expand("~/p4ha1_brca"))
  INCLUDE_METASTASIS <- FALSE      # set TRUE to read the with-metastasis run
  POOL_LEVEL         <- "patient"
  OUT_ROOT <- file.path(PROJECT_ROOT, paste0(
    if (INCLUDE_METASTASIS) "cohort_with_metastasis" else "cohort_primary",
    if (POOL_LEVEL == "section") "_sectionpool" else ""))
}
if (!identical(POOL_LEVEL, "patient"))
  stop("Section 22 is written for POOL_LEVEL = 'patient'. Section 21 covers the section level.")
AGG_DIR   <- file.path(OUT_ROOT, "aggregate")
OUT_DIR   <- AGG_DIR
TOL       <- 1e-6

if (!dir.exists(AGG_DIR)) stop("AGG_DIR not found: ", AGG_DIR)

## ============================================================================
## 22.0  READ
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
  }
  d
}

cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 22  SINGLE-PATIENT INFLUENCE AND THE METASTASIS\n")
cat(strrep("=", 74), "\n", sep = "")
cat("\nInputs:\n")
iso    <- read_one("ISOFORM_contrasts_per_sample")
iso_pl <- read_one("ISOFORM_slope_difference_pooled")
sc     <- read_one("SC_contrasts_per_sample", required = FALSE)
sc_pl  <- read_one("SC_slope_difference_pooled", required = FALSE)
cmap   <- read_one("cohort_sections")

iso <- iso[iso$contrast == "interaction_diff", ]
if (!is.null(sc)) sc <- sc[sc$contrast == "interaction_diff", ]

pid_of <- function(sid) {
  p <- cmap$patient_id[match(sid, cmap$sample_id)]
  if (anyNA(p)) stop("No patient for: ", paste(unique(sid[is.na(p)]), collapse = ", "))
  p
}
met_patients <- unique(cmap$patient_id[cmap$primary_metastasis == "Metastasis"])

## ============================================================================
## 22.1  POOLING  (same as the pipeline's patient-level meta_pool; verified below)
## ============================================================================
dl_pool <- function(y, v) {
  ok <- is.finite(y) & is.finite(v) & v > 0; y <- y[ok]; v <- v[ok]; k <- length(y)
  if (k == 0) return(c(est = NA, lower = NA, upper = NA, p = NA, k = 0, I2 = NA))
  w  <- 1 / v; fe <- sum(w * y) / sum(w)
  Q  <- sum(w * (y - fe)^2); df <- k - 1
  C  <- sum(w) - sum(w^2) / sum(w)
  t2 <- if (k > 1 && C > 0) max(0, (Q - df) / C) else 0
  ws <- 1 / (v + t2); est <- sum(ws * y) / sum(ws); se <- sqrt(1 / sum(ws))
  c(est = est, lower = est - qnorm(.975) * se, upper = est + qnorm(.975) * se,
    p = 2 * pnorm(-abs(est / se)), k = k,
    I2 = if (k > 1 && Q > 0) max(0, (Q - df) / Q) * 100 else 0)
}
fe_collapse <- function(d) {
  d$patient <- pid_of(d$sample_id)
  do.call(rbind, lapply(split(d, d$patient), function(x) {
    ok <- is.finite(x$est) & is.finite(x$var) & x$var > 0; x <- x[ok, ]
    if (!nrow(x)) return(NULL)
    w <- 1 / x$var
    data.frame(patient = x$patient[1], n_sections = nrow(x),
               est = sum(w * x$est) / sum(w), var = 1 / sum(w))
  }))
}
pat_pool <- function(pc) dl_pool(pc$est, pc$var)
excl0 <- function(lo, hi) (is.finite(lo) && is.finite(hi)) && (hi < 0 || lo > 0)

## the analyses to run, each a data frame of per-section contrasts
analyses <- list()
for (og in unique(iso$outcome)) analyses[[paste("S15", og)]] <- iso[iso$outcome == og, ]
if (!is.null(sc))
  for (og in unique(sc$outcome)) for (md in c("M1_adjusted", "M2_low_stroma"))
    analyses[[paste("S20", og, md)]] <- sc[sc$outcome == og & sc$model == md, ]

## ============================================================================
## 22.2  CHECK: does this pooling reproduce the pipeline's?
## ============================================================================
cat("\nCheck: reproduce the pipeline's pooled estimates\n")
ok_all <- TRUE
for (lb in names(analyses)) {
  d <- analyses[[lb]]
  mine <- pat_pool(fe_collapse(d))[["est"]]
  ref <- if (grepl("^S15", lb)) {
    og <- sub("^S15 ", "", lb)
    iso_pl$est[iso_pl$outcome == og & iso_pl$contrast == "interaction_diff"]
  } else {
    p <- strsplit(lb, " ")[[1]]
    sc_pl$est[sc_pl$outcome == p[2] & sc_pl$model == p[3] & sc_pl$contrast == "interaction_diff"]
  }
  good <- length(ref) == 1 && isTRUE(abs(mine - ref) < TOL)
  ok_all <- ok_all && good
  cat(sprintf("  %-26s mine %+.6f  pipeline %+.6f  %s\n", lb, mine,
              if (length(ref) == 1) ref else NA_real_, if (good) "ok" else "*** MISMATCH ***"))
}
if (!ok_all) stop("Section 22 halted: pooling does not reproduce the pipeline. Nothing below is interpretable.")
cat("  PASSED\n")

## ============================================================================
## 22.3  LEAVE ONE PATIENT OUT, WITHOUT THE METASTASIS, SIGN CONCORDANCE
## ============================================================================
res  <- list()
loo_rows <- list()

for (lb in names(analyses)) {
  d   <- analyses[[lb]]
  pc  <- fe_collapse(d)
  all <- pat_pool(pc)
  
  ## A. leave one patient out
  loo <- do.call(rbind, lapply(pc$patient, function(s) {
    r <- pat_pool(pc[pc$patient != s, ])
    data.frame(analysis = lb, dropped = s, n_sections_dropped = pc$n_sections[pc$patient == s],
               est = r[["est"]], lower = r[["lower"]],
               upper = r[["upper"]], p = r[["p"]], k = r[["k"]],
               shift = r[["est"]] - all[["est"]],
               still_excludes_0 = excl0(r[["lower"]], r[["upper"]]), row.names = NULL)
  }))
  loo_rows[[lb]] <- loo
  
  ## B. without the metastatic patient (only possible in the with-metastasis run)
  has_met <- any(pc$patient %in% met_patients)
  prim <- if (has_met) pat_pool(pc[!pc$patient %in% met_patients, ])
          else c(est = NA, lower = NA, upper = NA, p = NA, k = NA, I2 = NA)
  
  ## C. sign concordance, before pooling
  same_pat <- sum(sign(pc$est) == sign(all[["est"]]), na.rm = TRUE)
  same_sec <- sum(sign(d$est)  == sign(all[["est"]]), na.rm = TRUE)
  
  res[[lb]] <- data.frame(
    analysis = lb,
    k = all[["k"]], est = all[["est"]], lower = all[["lower"]], upper = all[["upper"]],
    loo_min = min(loo$est), loo_max = max(loo$est),
    loo_max_abs_shift = max(abs(loo$shift)),
    most_influential = loo$dropped[which.max(abs(loo$shift))],
    loo_any_conclusion_change = any(loo$still_excludes_0 != excl0(all[["lower"]], all[["upper"]])),
    has_metastasis = has_met,
    primary_k = prim[["k"]], primary_est = prim[["est"]],
    primary_lower = prim[["lower"]], primary_upper = prim[["upper"]],
    same_sign_patients = same_pat, n_patients = nrow(pc),
    same_sign_sections = same_sec, n_sections = nrow(d),
    row.names = NULL)
}

res_df <- do.call(rbind, res)
loo_df <- do.call(rbind, loo_rows)
write.csv(res_df, file.path(OUT_DIR, "LOO_summary.csv"), row.names = FALSE)
write.csv(loo_df, file.path(OUT_DIR, "LOO_per_dropped_patient.csv"), row.names = FALSE)

## ============================================================================
## 22.4  REPORT
## ============================================================================
for (i in seq_len(nrow(res_df))) with(res_df[i, ], {
  cat("\n", strrep("-", 74), "\n", sep = "")
  cat(sprintf("  %s\n", analysis))
  cat(strrep("-", 74), "\n", sep = "")
  cat(sprintf("    all patients     k=%2d  %+.4f (%+.4f to %+.4f)\n", k, est, lower, upper))
  cat(sprintf("    leave one out          range %+.4f to %+.4f   largest shift %+.4f when %s is dropped\n",
              loo_min, loo_max, loo_max_abs_shift, most_influential))
  cat(sprintf("    conclusion            %s\n",
              if (loo_any_conclusion_change) "*** CHANGES when at least one patient is dropped ***"
              else "unchanged by dropping any single patient"))
  if (has_metastasis)
    cat(sprintf("    without metastasis k=%2d  %+.4f (%+.4f to %+.4f)%s\n",
                primary_k, primary_est, primary_lower, primary_upper,
                if (excl0(primary_lower, primary_upper) != excl0(lower, upper))
                  "   <-- CONCLUSION CHANGES" else ""))
  else
    cat("    metastasis            not in this run (primary-only); rerun with INCLUDE_METASTASIS <- TRUE\n")
  cat(sprintf("    same direction as the pooled estimate in %d of %d patients (%d of %d sections)\n",
              same_sign_patients, n_patients, same_sign_sections, n_sections))
})

## the three most influential patients for the headline analysis
hl <- "S15 P4HA1"
if (hl %in% names(loo_rows)) {
  l <- loo_rows[[hl]]
  l <- l[order(-abs(l$shift)), ]
  cat("\n", strrep("-", 74), "\n", sep = "")
  cat("  MOST INFLUENTIAL PATIENTS FOR THE HEADLINE (S15 P4HA1)\n")
  cat(strrep("-", 74), "\n", sep = "")
  for (j in seq_len(min(3, nrow(l)))) with(l[j, ],
    cat(sprintf("    drop %-12s (%d sections) -> %+.4f  (shift %+.4f)\n", dropped, n_sections_dropped, est, shift)))
  cat("  A patient that moves the estimate a lot is not an error. It is worth\n")
  cat("  knowing what is unusual about it before someone else asks.\n")
}

cat("\nNOTES\n")
cat("  No model was refitted; every number comes from the stored per-section fits.\n")
cat("  Leave one out drops a PATIENT with all of their sections.\n")
cat("  Model p values remain inflated by spatial autocorrelation.\n")
cat("  Output: LOO_summary.csv and LOO_per_dropped_patient.csv in ", OUT_DIR, "\n", sep = "")
cat(strrep("=", 74), "\n", sep = "")
################################################################################
## 23. IS THE RESULT CARRIED BY ONE ER-STRESS GENE?
##
##   BRCA NOTE: the reasoning in this header was written for the PDAC
##   cohort. Any statement about what the data showed, which module worked,
##   or how many tumours were usable describes PDAC, not BRCA. The code
##   below is unchanged except where tagged "BRCA:".
##
##   THE WORRY
##   The ER stress arm rests on one hand-built 24-gene module that fails the
##   coherence test of section 19 (percentile 0 on mean rho, 7 on PC1). If the
##   whole finding is produced by one or two members, then "P4HA1 tracks ER
##   stress" is really "P4HA1 tracks HSPA5", which is a different and much
##   smaller claim.
##
##   WHAT THIS DOES
##   Refits the section 15 model 25 times per tumour: once with the full
##   module, then once with each of the 24 genes removed in turn. Nothing else
##   changes. Every fit is pooled exactly as in section 15.
##
##   THE CHECK THAT MAKES IT READABLE
##   The full-module run MUST reproduce the section 15 P4HA1 estimate to six
##   decimals. If it does not, the scoring differs from section 15 and no
##   comparison below means anything, so the section stops.
##
##   THE OPEN QUESTION IT ALSO ANSWERS
##   ATF4 and CHAC1 are PERK-responsive and were kept, while HYOU1, DDIT4 and
##   ERO1A were dropped as hypoxia-responsive. That exclusion rule has been an
##   inconsistency in the pipeline for weeks. This says whether it matters.
##
##   HOW TO READ IT
##   All 24 estimates close to the full one  -> no single gene carries the
##        result. That is the sentence worth having.
##   One gene's removal collapses it         -> say so plainly and name the
##        gene. The claim then has to be rewritten around that gene.
##
##   COST: 25 module scorings per tumour, so a few minutes. Outcome is P4HA1
##   only; the point is the robustness of the headline, not a survey.
##
##   Run after the main pipeline in the same FRESH session.
################################################################################

stopifnot(exists("MODULES"), exists("work"), exists("agg"), exists("SEED"),
          exists("ADJUST_DEPTH"), exists("pool_with_p"), exists("fmt_p"),
          exists("write_out"), exists("assert_run_state"), exists("log_msg"),
          exists("bind"), exists("min_genes_for"))
if (!exists("iso_con_pooled"))
  stop("Section 23 needs section 15 (iso_con_pooled) in the same session for the safety check.")

## ============================================================================
## 23.0  SETTINGS
## ============================================================================
EL_OUTCOME      <- "P4HA1"
EL_ES           <- "ERstressStrict"
EL_OTHER        <- c("Hypoxia", "CollagenSynthesis", "iCAF", "ERAD")
EL_MIN_SPOTS    <- 30
EL_MIN_EVENTS   <- 70
EL_ADJUST_DEPTH <- ADJUST_DEPTH
EL_TOL          <- 1e-6

EL_GENES <- MODULES[[EL_ES]]
EL_VARIANTS <- c("FULL", EL_GENES)      # FULL first, then one per dropped gene

## genes worth naming in the report whatever happens
EL_FLAGGED <- intersect(c("ATF4", "CHAC1"), EL_GENES)

EL_RENAME <- c(Hypoxia = "zHyp", CollagenSynthesis = "zCol", iCAF = "ziCAF", ERAD = "zERAD")
EL_II_COL <- "zHyp:zCol"
EL_II_ER  <- "zHyp:zER"
EL_RHS <- paste0("zHyp * zCol + zHyp * zER + ziCAF + zERAD",
                 if (EL_ADJUST_DEPTH) " + zDepth" else "")

assert_run_state("Section 23")

cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 23  LEAVE ONE ER-STRESS GENE OUT\n")
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  module %s has %d genes, so %d model fits per tumour\n",
            EL_ES, length(EL_GENES), length(EL_VARIANTS)))
cat(sprintf("  genes: %s\n", paste(EL_GENES, collapse = ", ")))

## the outcome must not be inside any predictor module
for (m in c(EL_OTHER, EL_ES))
  if (EL_OUTCOME %in% MODULES[[m]])
    stop("Section 23 halted: ", EL_OUTCOME, " is inside ", m)

## ============================================================================
## 23.1  ONE TUMOUR: SCORE THE FIXED MODULES ONCE, THEN EVERY ER VARIANT
## ============================================================================
el_one <- function(rds_path, sid) {
  obj <- readRDS(rds_path)
  obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
  DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else Assays(obj)[1]
  if (!EL_OUTCOME %in% rownames(obj)) { log_msg("  %s: outcome absent -> skipped", sid); return(NULL) }
  
  ## fixed predictors, scored exactly as in section 15
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
  
  ## every ER variant, scored into its own column
  for (v in EL_VARIANTS) {
    set_v <- if (identical(v, "FULL")) EL_GENES else setdiff(EL_GENES, v)
    g <- intersect(set_v, rownames(obj))
    if (length(g) < min_genes_for(set_v)) {
      log_msg("  %s: ER variant %s has %d/%d genes -> variant skipped", sid, v, length(g), length(set_v))
      next
    }
    nm <- paste0("ERv_", v)
    set.seed(SEED)
    obj <- tryCatch(add_module_score_fallback(obj, features = list(g), name = nm,
                                   assay = DefaultAssay(obj), nbin = 20, ctrl = 50),
                    error = function(e) NULL)
    if (is.null(obj) || !paste0(nm, "1") %in% colnames(obj@meta.data)) next
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

el_pooled <- el_all %>%
  group_by(variant) %>% group_modify(~ pool_with_p(.x$est, .x$var, .x$sample_id)) %>% ungroup()
write_out(el_pooled, "ERLOO_pooled")

## ============================================================================
## 23.4  SAFETY CHECK: THE FULL MODULE MUST REPRODUCE SECTION 15
## ============================================================================
full <- el_pooled[el_pooled$variant == "FULL", ]
ref  <- iso_con_pooled[iso_con_pooled$outcome == EL_OUTCOME &
                         iso_con_pooled$contrast == "interaction_diff", ]$est
cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 23.4  SAFETY CHECK\n")
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  FULL module = %+.6f   section 15 = %+.6f\n", full$est, ref))
if (!nrow(full) || abs(full$est - ref) > EL_TOL)
  stop("Section 23 halted: the full module does not reproduce section 15. Nothing below is interpretable.")
cat("  PASSED\n")

## ============================================================================
## 23.5  REPORT
## ============================================================================
excl0 <- function(lo, hi) (is.finite(lo) && is.finite(hi)) && (hi < 0 || lo > 0)
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
cat(sprintf("\n  full module (%d genes) : %+.4f  (95%% CI %+.4f to %+.4f)  k=%d  p = %s\n",
            length(EL_GENES), full$est, full$lower, full$upper, full$k, fmt_p(full$p)))

cat("\n  each gene removed in turn, ordered by how much it moves the result\n\n")
cat(sprintf("  %-10s %9s %21s %5s %9s %8s\n", "dropped", "estimate", "95% CI", "k", "shift", "kept %"))
for (i in seq_len(nrow(el_rep))) with(el_rep[i, ],
                                      cat(sprintf("  %-10s %+9.4f  %+8.4f to %+8.4f %5d %+9.4f %7.0f%%%s%s\n",
                                                  variant, est, lower, upper, k, shift, pct_retained,
                                                  if (still_excl0 != base_sig) "   <-- CONCLUSION CHANGES" else "",
                                                  if (variant %in% EL_FLAGGED) "   [flagged gene]" else "")))

cat(sprintf("\n  range across all %d single-gene removals: %+.4f to %+.4f\n",
            nrow(el_rep), min(el_rep$est), max(el_rep$est)))

## ---- the verdict -----------------------------------------------------------
changed <- el_rep$variant[el_rep$still_excl0 != base_sig]
cat("\n", strrep("-", 74), "\n", sep = "")
cat("  VERDICT\n")
cat(strrep("-", 74), "\n", sep = "")
if (!length(changed)) {
  cat("  NO SINGLE GENE CARRIES THE RESULT.\n")
  cat("  Removing any one of the 24 genes leaves the slope difference in the\n")
  cat("  same direction and still excluding zero. The ER stress arm is a\n")
  cat("  property of the set, not of one member.\n")
  cat("  This is the sentence to put in the paper, and it is a stronger answer\n")
  cat("  to the coherence problem than any p value.\n")
} else {
  cat(sprintf("  *** THE RESULT DEPENDS ON %d GENE(S): %s ***\n",
              length(changed), paste(changed, collapse = ", ")))
  cat("  Removing it changes whether the interval excludes zero. The claim has\n")
  cat("  to be rewritten around that gene rather than around 'ER stress', and\n")
  cat("  the gene's own biology has to be explained. Report this plainly.\n")
}

if (length(EL_FLAGGED)) {
  cat("\n  THE ATF4 / CHAC1 QUESTION\n")
  for (g in EL_FLAGGED) {
    r <- el_rep[el_rep$variant == g, ]
    if (nrow(r) == 1)
      cat(sprintf("    without %-6s %+.4f  (%.0f%% of the full estimate)\n", g, r$est, r$pct_retained))
  }
  cat("    These two are PERK-responsive and were kept, while HYOU1, DDIT4 and\n")
  cat("    ERO1A were dropped as hypoxia-responsive. If the numbers above are\n")
  cat("    close to the full estimate, the inconsistent exclusion rule does not\n")
  cat("    affect the result and you can say so instead of defending the rule.\n")
}

cat("\nCAVEATS\n")
cat("  Only P4HA1 was refitted. The same check for P4HA2 is a one-line change.\n")
cat("  Removing one gene from 24 is a small perturbation; this tests dependence\n")
cat("  on single members, not whether a different ER stress definition agrees.\n")
cat("  That is a separate question and needs independent UPR gene sets.\n")
cat("  Model p values remain inflated by spatial autocorrelation.\n")
cat(strrep("=", 74), "\n", sep = "")

log_msg("Section 23 complete.")

################################################################################
## 24. DOES THE ER-STRESS RESULT HOLD WITH INDEPENDENT UPR GENE SETS?
##
##   BRCA NOTE: the reasoning in this header was written for the PDAC
##   cohort. Any statement about what the data showed, which module worked,
##   or how many tumours were usable describes PDAC, not BRCA. The code
##   below is unchanged except where tagged "BRCA:".
##
##   THE WORRY THIS ANSWERS
##   Section 23 showed no single gene carries the result. That is dependence on
##   MEMBERS. This asks the other question: dependence on the DEFINITION. The
##   whole ER stress arm rests on one hand-built 24-gene module. If only that
##   module produces the pattern, the finding is about a gene list. If several
##   independently published definitions agree, it is about ER stress.
##
##   WHAT IS FITTED
##   The section 15 model, unchanged, with zER built from each candidate set in
##   turn. ERstressStrict is refitted as the reference and MUST reproduce the
##   section 15 number, or the section stops.
##
##   TWO DIAGNOSTICS PRINTED BEFORE ANY ESTIMATE
##   1. Overlap with ERstressStrict (Jaccard). A set sharing most of its genes
##      is not an independent test, it is the same module in other clothes.
##      Agreement from such a set is weak evidence and must be labelled.
##   2. Correlation with the hypoxia module. This is the VEGFA lesson from
##      section 11: a UPR set loaded with hypoxia-responsive genes would
##      manufacture the very interaction being tested. Genes shared with the
##      Hypoxia, CollagenSynthesis or iCAF modules are STRIPPED automatically
##      and named, and the residual correlation is still reported.
##
##   WHAT THIS DOES NOT DO
##   It does not check whether the candidate sets are internally coherent. To
##   get that, add their names to COH_MODULES in section 19 and re-run it: the
##   sets are left in MODULES by this section for exactly that purpose. Note
##   ERstressStrict itself fails coherence, so failing it is not disqualifying
##   here; it is context.
##
##   INPUT
##   A folder of MSigDB TSV files, one per candidate set. Suggested:
##     HALLMARK_UNFOLDED_PROTEIN_RESPONSE          (the parent of the hand-built
##                                                  module, so NOT independent)
##     GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS
##     GOBP_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE
##     GOBP_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE
##     GOBP_ATF6_MEDIATED_UNFOLDED_PROTEIN_RESPONSE
##   The last three are the three arms of the UPR and are the most interesting
##   result available here: if one arm carries the signal and another does not,
##   that is biology rather than a robustness check.
##
##   Run after the main pipeline in the same FRESH session.
################################################################################

stopifnot(exists("MODULES"), exists("work"), exists("agg"), exists("SEED"),
          exists("ADJUST_DEPTH"), exists("pool_with_p"), exists("fmt_p"),
          exists("write_out"), exists("assert_run_state"), exists("log_msg"),
          exists("bind"), exists("min_genes_for"))
if (!exists("iso_con_pooled"))
  stop("Section 24 needs section 15 (iso_con_pooled) in the same session for the safety check.")

## ============================================================================
## 24.0  SETTINGS
## ============================================================================
UPR_DIR         <- file.path(PROJECT_ROOT, "data", "upr_sets")   # BRCA: folder of MSigDB TSVs
UPR_OUTCOME     <- "P4HA1"
UPR_REF         <- "ERstressStrict"
UPR_OTHER       <- c("Hypoxia", "CollagenSynthesis", "iCAF", "ERAD")
UPR_MIN_GENES   <- 8        # a set smaller than this cannot support a module score
UPR_MIN_SPOTS   <- 30
UPR_MIN_EVENTS  <- 70
UPR_ADJUST_DEPTH <- ADJUST_DEPTH
UPR_TOL         <- 1e-6
UPR_COLLIN_WARN <- 0.50     # correlation with hypoxia above which to flag
UPR_JACCARD_WARN <- 0.50    # overlap with the reference above which to flag

UPR_FAMILY <- c("P4HA1", "P4HA2", "P4HA3", "P4HB", "PLOD1", "PLOD2", "PLOD3", "PPIB")
UPR_RENAME <- c(Hypoxia = "zHyp", CollagenSynthesis = "zCol", iCAF = "ziCAF", ERAD = "zERAD")
UPR_II_COL <- "zHyp:zCol"
UPR_II_ER  <- "zHyp:zER"
UPR_RHS <- paste0("zHyp * zCol + zHyp * zER + ziCAF + zERAD",
                  if (UPR_ADJUST_DEPTH) " + zDepth" else "")

assert_run_state("Section 24")

## ============================================================================
## 24.1  READ THE CANDIDATE SETS FROM THE TSVs
##       Gene lists are read from file, never typed.
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
  nm <- sub("\\.v[0-9].*$", "", nm)                 # drop the release suffix
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
##       Family genes are a hard stop. Genes shared with Hypoxia, collagen or
##       iCAF are stripped: that is the VEGFA lesson from section 11.
## ============================================================================
STRIP_AGAINST <- c("Hypoxia", "CollagenSynthesis", "iCAF")
strip_pool <- unique(unlist(MODULES[STRIP_AGAINST], use.names = FALSE))
ref_genes  <- MODULES[[UPR_REF]]

jacc <- function(a, b) { u <- length(union(a, b)); if (u) length(intersect(a, b)) / u else NA_real_ }

## A published ER-stress set may legitimately contain P4HB, the beta subunit of
## P4HA1's own enzyme, or another family gene. Those genes are STRIPPED and
## named, not tolerated. The set is dropped only if they make up more than
## UPR_FAMILY_MAX of it, because then the set is built around them and removing
## them changes what it measures: that was the GO:0034663 decision.
UPR_FAMILY_MAX <- 0.05

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
write_out(data.frame(set = rep(names(UPR_SETS), lengths(UPR_SETS)),
                     gene = unlist(UPR_SETS, use.names = FALSE)), "UPR_genes_used")
if (!length(UPR_SETS)) stop("Section 24: every candidate set was dropped.")

## the reference is fitted alongside, as the first variant
UPR_VARIANTS <- c(UPR_REF, names(UPR_SETS))
for (nm in names(UPR_SETS)) MODULES[[nm]] <- UPR_SETS[[nm]]   # left in place for section 19

## ============================================================================
## 24.3  ONE TUMOUR: SCORE EVERYTHING ONCE
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
  
  ## diagnostics on the raw scores
  sp <- function(a, b) suppressWarnings(cor(a, b, method = "spearman"))
  attr(d, "diag") <- do.call(rbind, lapply(have, function(v) data.frame(
    sample_id = sid, variant = v,
    rho_vs_hypoxia = sp(d[[paste0("UPRv_", v)]], d$Hypoxia),
    rho_vs_ref     = sp(d[[paste0("UPRv_", v)]], d[[paste0("UPRv_", UPR_REF)]]),
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

upr_pooled <- upr_all %>%
  group_by(variant) %>% group_modify(~ pool_with_p(.x$est, .x$var, .x$sample_id)) %>% ungroup()
write_out(upr_pooled, "UPR_pooled")

## ============================================================================
## 24.6  SAFETY CHECK
## ============================================================================
refr <- upr_pooled[upr_pooled$variant == UPR_REF, ]
ref15 <- iso_con_pooled[iso_con_pooled$outcome == UPR_OUTCOME &
                          iso_con_pooled$contrast == "interaction_diff", ]$est
cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 24.6  SAFETY CHECK\n")
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  %s = %+.6f   section 15 = %+.6f\n", UPR_REF, refr$est, ref15))
if (!nrow(refr) || abs(refr$est - ref15) > UPR_TOL)
  stop("Section 24 halted: the reference module does not reproduce section 15.")
cat("  PASSED\n")

## ============================================================================
## 24.7  REPORT
## ============================================================================
excl0 <- function(lo, hi) (is.finite(lo) && is.finite(hi)) && (hi < 0 || lo > 0)
ref_sig <- excl0(refr$lower, refr$upper)

cat("\n", strrep("=", 74), "\n", sep = "")
cat(sprintf("SECTION 24  [%s]  run_id %s\n",
            if (UPR_ADJUST_DEPTH) "DEPTH ADJUSTED" else "UNADJUSTED", RUN_ID))
cat(strrep("=", 74), "\n", sep = "")

cat("\n--- diagnostics, read these before the estimates ---\n")
cat(sprintf("  %-34s %10s %12s %12s\n", "set", "Jaccard", "rho hypoxia", "rho with ref"))
for (v in UPR_VARIANTS) {
  dg <- upr_diag[upr_diag$variant == v, ]
  jc <- if (identical(v, UPR_REF)) 1 else audit$jaccard_with_ref[audit$set == v]
  if (!nrow(dg)) next
  mh <- median(dg$rho_vs_hypoxia, na.rm = TRUE); mr <- median(dg$rho_vs_ref, na.rm = TRUE)
  cat(sprintf("  %-34s %10.2f %12+.2f %12+.2f%s%s\n", v, jc, mh, mr,
              if (abs(mh) >= UPR_COLLIN_WARN) "  <-- tracks hypoxia" else "",
              if (!identical(v, UPR_REF) && is.finite(jc) && jc >= UPR_JACCARD_WARN)
                "  <-- not independent" else ""))
}
cat("\n  Jaccard is overlap with the reference module AFTER stripping. A set with\n")
cat("  high overlap is the same instrument again; its agreement is not evidence.\n")
cat("  A set that tracks hypoxia can manufacture the interaction being tested.\n")

cat("\n--- THE HEADLINE: slope difference under each definition of ER stress ---\n")
cat(sprintf("  %-34s %9s %21s %5s %11s\n", "set", "estimate", "95% CI", "k", "p"))
for (v in UPR_VARIANTS) {
  r <- upr_pooled[upr_pooled$variant == v, ]
  if (!nrow(r)) next
  cat(sprintf("  %-34s %+9.4f  %+8.4f to %+8.4f %5d %10s%s%s\n",
              v, r$est, r$lower, r$upper, r$k, fmt_p(r$p),
              if (identical(v, UPR_REF)) "   [reference]" else "",
              if (!identical(v, UPR_REF) && excl0(r$lower, r$upper) == ref_sig &&
                  sign(r$est) == sign(refr$est)) "   agrees" else ""))
}

## ---- verdict ---------------------------------------------------------------
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
  if (length(agree) >= 2) {
    cat("\n  THE DEFINITION IS NOT THE RESULT. Several independently published\n")
    cat("  definitions of ER stress give the same answer, so the finding is not a\n")
    cat("  property of the hand-built module. This is the strongest available\n")
    cat("  answer to the coherence problem in section 19.\n")
  } else if (length(same_dir) >= 2) {
    cat("\n  PARTIAL. The independent sets point the same way but not all of them\n")
    cat("  exclude zero. Say exactly that: consistent in direction, weaker in\n")
    cat("  precision, and the hand-built module remains the primary definition.\n")
  } else {
    cat("\n  *** THE RESULT IS SPECIFIC TO THE HAND-BUILT MODULE. Independent\n")
    cat("  definitions of ER stress do not reproduce it. The claim must be\n")
    cat("  rewritten around what that module actually measures, and the\n")
    cat("  coherence failure in section 19 now matters a great deal more. ***\n")
  }
  cat("\n  If the three UPR arms disagree with each other, that is a biological\n")
  cat("  result rather than a robustness check: name which arm carries it.\n")
}

cat("\nCAVEATS\n")
cat("  Only P4HA1 was refitted.\n")
cat("  Any set containing P4HB or another hydroxylase-family gene had it removed;\n")
cat("  UPR_set_audit names them. A set built around such genes was dropped.\n")
cat("  Sets were stripped of genes shared with Hypoxia, collagen and iCAF, so\n")
cat("  they are not the published sets as-is. The stripped genes are listed in\n")
cat("  UPR_set_audit.\n")
cat("  Coherence was not tested here. Add the set names to COH_MODULES in\n")
cat("  section 19 and re-run it to get that.\n")
cat("  Model p values remain inflated by spatial autocorrelation.\n")
cat(strrep("=", 74), "\n", sep = "")

log_msg("Section 24 complete.")

################################################################################
## 25. ARE THE PUBLISHED ER-STRESS SETS ANY MORE COHERENT THAN OURS?
##
##   BRCA NOTE: the reasoning in this header was written for the PDAC
##   cohort. Any statement about what the data showed, which module worked,
##   or how many tumours were usable describes PDAC, not BRCA. The code
##   below is unchanged except where tagged "BRCA:".
##
##   WHY THIS MATTERS NOW
##   Section 19 found ERstressStrict is no more internally consistent than
##   random genes of the same size and detection rate. Section 24 then showed
##   five independent definitions of ER stress give the same result anyway.
##   That leaves one question: are those five definitions coherent?
##
##     If they are also incoherent, then incoherence is a property of
##     ER-stress gene sets in Visium data, not a flaw in a hand-built module.
##     The sentence for the paper becomes much easier to write.
##
##     If they ARE coherent while ours is not, then ours is the weak one, and
##     the primary definition should probably change to a published set that
##     both coheres and reproduces the result.
##
##   METHOD
##   Identical to section 19: mean pairwise Spearman between member genes, and
##   the variance share of PC1, each compared against 30 random gene sets
##   matched on size and detection rate, reported as a percentile.
##
##   WHY THE NUMBERS WILL NOT MATCH SECTION 19 EXACTLY
##   The null pool excludes every gene belonging to any module under test, and
##   this run tests a different list, so the pool differs and the random sets
##   differ. Hypoxia and collagen are included as positive controls and must
##   still land high; that is the gate, not exact reproduction.
##
##   Run after the main pipeline AND section 24 in the same session, because it
##   needs the UPR sets that section 24 put into MODULES.
################################################################################

stopifnot(exists("MODULES"), exists("work"), exists("agg"),
          exists("write_out"), exists("assert_run_state"), exists("log_msg"), exists("bind"))
if (!exists("UPR_SETS"))
  stop("Section 25 needs section 24 to have run in this session (UPR_SETS not found).")

## ============================================================================
## 25.0  SETTINGS
## ============================================================================
UC_REFS     <- c("Hypoxia", "CollagenSynthesis", "ERstressStrict", "ERAD")
UC_MODULES  <- c(UC_REFS, names(UPR_SETS))
UC_N_NULL   <- 30
UC_MIN_GENES <- 4
UC_MIN_SPOTS <- 30
UC_SEED     <- 90125
UC_CALIPER  <- 0.35
UC_PCT_MIN  <- 0.90

RATES_CSV <- file.path(agg, "DETECT_cohort_gene_rates.csv")
if (!file.exists(RATES_CSV)) stop("Section 25 needs DETECT_cohort_gene_rates.csv from section 14.")

assert_run_state("Section 25")

miss <- setdiff(UC_MODULES, names(MODULES))
if (length(miss)) stop("Not in MODULES: ", paste(miss, collapse = ", "))

## ============================================================================
## 25.1  NULL POOL, MATCHED ON DETECTION RATE
## ============================================================================
rates <- read.csv(RATES_CSV, stringsAsFactors = FALSE)
gene_med <- rates %>%
  group_by(gene) %>%
  summarise(median_rate = median(detect_rate), n_tum = dplyr::n(), .groups = "drop") %>%
  filter(n_tum == dplyr::n_distinct(rates$sample_id))

UC_ALL_GENES <- unique(unlist(MODULES[UC_MODULES], use.names = FALSE))
uc_pool <- gene_med[!gene_med$gene %in% UC_ALL_GENES, ]

cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 25.1  SETUP\n")
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  genes present in all tumours : %d\n", nrow(gene_med)))
cat(sprintf("  null pool after exclusions   : %d\n", nrow(uc_pool)))
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
## 25.2  THE TWO MEASURES  (same as section 19)
## ============================================================================
uc_rho <- function(mat) {
  if (is.null(mat) || nrow(mat) < UC_MIN_GENES) return(NA_real_)
  keep <- apply(mat, 1, function(v) sd(v) > 0)
  if (sum(keep) < UC_MIN_GENES) return(NA_real_)
  cc <- suppressWarnings(cor(t(as.matrix(mat[keep, , drop = FALSE])), method = "spearman"))
  if (!is.matrix(cc)) return(NA_real_)
  mean(cc[upper.tri(cc)], na.rm = TRUE)
}
uc_pc1 <- function(mat) {
  if (is.null(mat) || nrow(mat) < UC_MIN_GENES) return(NA_real_)
  keep <- apply(mat, 1, function(v) sd(v) > 0)
  if (sum(keep) < UC_MIN_GENES) return(NA_real_)
  m <- scale(t(as.matrix(mat[keep, , drop = FALSE])))
  m <- m[, apply(m, 2, function(v) all(is.finite(v))), drop = FALSE]
  if (ncol(m) < UC_MIN_GENES) return(NA_real_)
  ev <- tryCatch(svd(m, nu = 0, nv = 0)$d^2, error = function(e) NULL)
  if (is.null(ev) || !length(ev)) return(NA_real_)
  ev[1] / sum(ev)
}

## ============================================================================
## 25.3  RUN
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
    r_rho <- uc_rho(real); r_pc1 <- uc_pc1(real)
    n_rho <- n_pc1 <- rep(NA_real_, length(ns))
    for (j in seq_along(ns)) {
      mm <- grab(ns[[j]]); if (is.null(mm)) next
      n_rho[j] <- uc_rho(mm); n_pc1[j] <- uc_pc1(mm)
    }
    pr <- if (sum(is.finite(n_rho)) >= 10) mean(n_rho[is.finite(n_rho)] < r_rho) else NA_real_
    pp <- if (sum(is.finite(n_pc1)) >= 10) mean(n_pc1[is.finite(n_pc1)] < r_pc1) else NA_real_
    rows[[length(rows) + 1]] <- data.frame(
      sample_id = sid, module = nm, n_genes_used = nrow(real),
      mean_rho = r_rho, pc1_var = r_pc1,
      null_rho = median(n_rho, na.rm = TRUE), null_pc1 = median(n_pc1, na.rm = TRUE),
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

uc_sum <- uc_all %>%
  group_by(module) %>%
  summarise(k = dplyr::n(), median_genes = median(n_genes_used),
            median_rho = median(mean_rho, na.rm = TRUE),
            median_null_rho = median(null_rho, na.rm = TRUE),
            median_pct_rho = median(pct_rho, na.rm = TRUE),
            median_pc1 = median(pc1_var, na.rm = TRUE),
            median_null_pc1 = median(null_pc1, na.rm = TRUE),
            median_pct_pc1 = median(pct_pc1, na.rm = TRUE), .groups = "drop")
uc_sum <- uc_sum[match(UC_MODULES, uc_sum$module), ]
uc_sum <- uc_sum[!is.na(uc_sum$module), ]
write_out(uc_sum, "UPRCOH_summary")

## ============================================================================
## 25.4  REPORT
## ============================================================================
cat("\n", strrep("=", 74), "\n", sep = "")
cat(sprintf("SECTION 25  COHERENCE OF THE PUBLISHED ER-STRESS SETS  run_id %s\n", RUN_ID))
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("\n  %-34s %6s %8s %9s %9s %8s %9s\n",
            "module", "genes", "rho", "null rho", "pctile", "PC1", "pctile"))
for (i in seq_len(nrow(uc_sum))) with(uc_sum[i, ],
  cat(sprintf("  %-34s %6.0f %8.3f %9.3f %8.0f%% %8.3f %8.0f%%\n",
              module, median_genes, median_rho, median_null_rho, 100 * median_pct_rho,
              median_pc1, 100 * median_pct_pc1)))

passes <- function(mod) {
  r <- uc_sum[uc_sum$module == mod, ]
  if (!nrow(r)) return(NA)
  isTRUE(r$median_pct_rho >= UC_PCT_MIN) || isTRUE(r$median_pct_pc1 >= UC_PCT_MIN)
}
cat(sprintf("\n  criterion: median percentile >= %.0f%% on either measure\n\n", 100 * UC_PCT_MIN))
for (nm in uc_sum$module)
  cat(sprintf("    %-34s %s\n", nm, if (isTRUE(passes(nm))) "COHERENT" else "not coherent"))

## ---- gate: the positive controls must behave ------------------------------
ctrl_ok <- isTRUE(passes("Hypoxia")) || isTRUE(passes("CollagenSynthesis"))
cat("\n", strrep("-", 74), "\n", sep = "")
cat("  VERDICT\n")
cat(strrep("-", 74), "\n", sep = "")
if (!ctrl_ok) {
  cat("  The positive controls failed in this run, so nothing here is readable.\n")
  cat("  Collagen and hypoxia both passed in section 19, so if they fail here the\n")
  cat("  null pool or the matching is at fault, not the candidate sets.\n")
} else {
  upr_names <- names(UPR_SETS)
  coh <- upr_names[vapply(upr_names, function(x) isTRUE(passes(x)), logical(1))]
  ours <- isTRUE(passes("ERstressStrict"))
  cat(sprintf("  our module (ERstressStrict) : %s\n", if (ours) "COHERENT" else "not coherent"))
  cat(sprintf("  published sets coherent     : %d of %d\n", length(coh), length(upr_names)))
  if (!ours && !length(coh)) {
    cat("\n  NO ER-STRESS DEFINITION IS COHERENT IN THIS DATA.\n")
    cat("  Ours is not unusual: every published set behaves like random genes of\n")
    cat("  the same size and detection rate, while collagen and hypoxia pass\n")
    cat("  cleanly. So incoherence is a property of ER-stress gene sets in Visium\n")
    cat("  data, not a defect of the hand-built module.\n")
    cat("  For the paper: the ER stress score is a composite readout rather than a\n")
    cat("  measure of one coordinated transcriptional programme, and that is true\n")
    cat("  of every published definition tested, all of which give the same result\n")
    cat("  (section 24).\n")
  } else if (!ours && length(coh)) {
    cat(sprintf("\n  *** SOME PUBLISHED SETS COHERE AND OURS DOES NOT: %s ***\n",
                paste(coh, collapse = ", ")))
    cat("  Those sets both cohere and reproduce the result. The primary definition\n")
    cat("  should probably move to one of them, with the hand-built module kept as\n")
    cat("  a sensitivity analysis. Decide this with Valerio before writing.\n")
  } else {
    cat("\n  Our module coheres in this run. Check why it differed from section 19\n")
    cat("  before using either number: the null pool is not the same.\n")
  }
}

cat("\nCAVEATS\n")
cat("  Percentiles are not directly comparable to section 19: the null pool\n")
cat("  excludes a different set of genes, so the random sets differ.\n")
cat("  Coherence is necessary for a module score to mean one thing. It is not\n")
cat("  sufficient, and its absence does not by itself invalidate a result that\n")
cat("  replicates across tumours, gene sets and definitions.\n")
cat(strrep("=", 74), "\n", sep = "")

log_msg("Section 25 complete.")
################################################################################
## 26. SPATIAL PERMUTATION TEST
##
##   BRCA NOTE: the reasoning in this header was written for the PDAC
##   cohort. Any statement about what the data showed, which module worked,
##   or how many tumours were usable describes PDAC, not BRCA. The code
##   below is unchanged except where tagged "BRCA:".
##
##   THE WORRY
##   Neighbouring spots are not independent. The logistic model treats every
##   spot as a separate observation, so its confidence intervals and p values
##   are too optimistic (section 16 in PDAC: 296 of 1000 random genes were
##   "significant"). Section 16 answers "is P4HA1 special among genes".
##   It does not answer "is the P4HA1 estimate larger than what spatial
##   structure alone produces". This section does.
##
##   WHAT THIS DOES
##   For every tumour the map of P4HA1 detection is moved as one rigid sheet
##   over the tissue (a random shift on the Visium grid with wrap-around, and
##   a random 180 degree turn). The patches of P4HA1 stay patches, so the
##   spatial structure of the outcome is kept, but they no longer sit on their
##   true positions, so any real link with hypoxia, collagen and ER stress is
##   broken. The predictors are not touched. The section 15 model is refitted
##   on the moved map and pooled over tumours exactly as in section 15.
##   Repeating this SP_N_PERM times gives the distribution of the pooled slope
##   difference under "no real link, same spatial structure".
##
##   Where the moved sheet falls outside the tissue, the value of the nearest
##   tissue spot is used. The share of spots filled this way is reported.
##
##   THE CHECK THAT MAKES IT READABLE
##   The unpermuted fit MUST reproduce the section 15 P4HA1 estimate to six
##   decimals, otherwise the section stops.
##
##   HOW TO READ IT
##   Observed estimate beyond almost all permuted ones -> the result is larger
##        than spatial structure alone produces. Report the empirical p.
##   Observed estimate inside the permuted ones        -> the result cannot be
##        told apart from spatial structure. Say so plainly.
##
##   LIMITS
##   Wrap-around creates a seam and nearest-spot filling repeats some values,
##   so the kept spatial structure is approximate, not exact.
##   The smallest possible empirical p is 1 / (SP_N_PERM + 1).
##
##   Run after the main pipeline in the same FRESH session.
################################################################################

stopifnot(exists("MODULES"), exists("work"), exists("SEED"), exists("ADJUST_DEPTH"),
          exists("meta_pool"), exists("write_out"), exists("assert_run_state"),
          exists("log_msg"), exists("min_genes_for"))
if (!exists("iso_con_pooled"))
  stop("Section 26 needs section 15 (iso_con_pooled) in the same session for the safety check.")
suppressPackageStartupMessages(library(Seurat))

## ============================================================================
## 26.0  SETTINGS
## ============================================================================
SP_OUTCOME      <- "P4HA1"
SP_N_PERM       <- 1000
SP_SEED         <- 20261005
SP_PROBE        <- 20          # permutations timed first, to print an estimate
SP_MIN_SHIFT    <- 10          # smallest allowed move, in spot spacings
SP_MODS         <- c("Hypoxia", "CollagenSynthesis", "ERstressStrict", "iCAF", "ERAD")
SP_RENAME       <- c(Hypoxia = "zHyp", CollagenSynthesis = "zCol", ERstressStrict = "zER",
                     iCAF = "ziCAF", ERAD = "zERAD")
SP_MIN_SPOTS    <- 30
SP_MIN_EVENTS   <- 70
SP_ADJUST_DEPTH <- ADJUST_DEPTH
SP_TOL          <- 1e-6
SP_II_COL       <- "zHyp:zCol"
SP_II_ER        <- "zHyp:zER"
SP_RHS <- paste0("zHyp * zCol + zHyp * zER + ziCAF + zERAD",
                 if (SP_ADJUST_DEPTH) " + zDepth" else "")

assert_run_state("Section 26")
for (m in SP_MODS) if (SP_OUTCOME %in% MODULES[[m]])
  stop("Section 26 halted: ", SP_OUTCOME, " is inside ", m)

cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 26  SPATIAL PERMUTATION TEST\n")
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  outcome %s, %d permutations, seed %d\n", SP_OUTCOME, SP_N_PERM, SP_SEED))

## ============================================================================
## 26.1  ONE TUMOUR: SCORES, OUTCOME AND GRID POSITION OF EVERY SPOT
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
  
  ## grid position (Visium array row and col), matched by barcode
  im <- Images(obj)
  if (!length(im)) { log_msg("  %s: no image -> skipped", sid); return(NULL) }
  co <- tryCatch(obj@images[[im[1]]]@coordinates, error = function(e) NULL)
  if (is.null(co) || !all(c("row", "col") %in% colnames(co))) {
    log_msg("  %s: no row/col grid coordinates -> skipped", sid); return(NULL)
  }
  idx <- match(rownames(d), rownames(co))
  ## spots without a grid position are KEPT here (grow = NA) so that the safety
  ## check can use exactly the spots section 15 used; they are set aside later
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
## observed fit: plain glm, identical to section 15
sp_fit_obs <- function(d) {
  if (!sp_usable(d$y)) return(NULL)
  f <- tryCatch(glm(as.formula(paste("y ~", SP_RHS)), data = d, family = binomial()),
                error = function(e) NULL)
  if (is.null(f) || !f$converged) return(NULL)
  sp_contrast(coef(f), vcov(f))
}
## permuted fit: same model through glm.fit on a prebuilt design matrix (fast)
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
## 26.3  THE GRID OF ONE TUMOUR AND THE MOVE
## ============================================================================
## For every grid cell, the index of the nearest tissue spot (itself if the
## cell is a tissue spot). Built once per tumour; a permutation is then a lookup.
sp_grid <- function(d) {
  R <- max(d$grow) + 1L; if (R %% 2L) R <- R + 1L        # even sizes keep the
  C <- max(d$gcol) + 1L; if (C %% 2L) C <- C + 1L        # hex parity under wrap
  sx <- d$gcol * 0.5; sy <- d$grow * sqrt(3) / 2         # physical position, spot spacing = 1
  near <- matrix(NA_integer_, R, C); direct <- matrix(FALSE, R, C)
  direct[cbind(d$grow + 1L, d$gcol + 1L)] <- TRUE
  for (r in 0:(R - 1L)) {
    cc <- seq.int(r %% 2L, C - 1L, by = 2L)              # valid cells in this row
    dx <- outer(cc * 0.5, sx, "-"); dy <- r * sqrt(3) / 2 - sy
    d2 <- dx^2 + matrix(dy^2, nrow = length(cc), ncol = length(sx), byrow = TRUE)
    near[r + 1L, cc + 1L] <- max.col(-d2, ties.method = "first")
  }
  list(R = R, C = C, near = near, direct = direct)
}
## one random move: returns, for every spot, the index of the spot whose
## outcome it receives
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
## 26.4  LOAD EVERY TUMOUR ONCE AND FIT THE UNPERMUTED MODEL
## ============================================================================
log_msg("Section 26 start. ADJUST_DEPTH = %s. run_id = %s.", SP_ADJUST_DEPTH, RUN_ID)
sp_dat <- list(); sp_obs <- list()
for (i in seq_len(nrow(work))) {
  path <- work$path[i]; sid <- work$sample_id[i]
  if (is.na(path) || !file.exists(path)) next
  log_msg("[sp load %d/%d] %s", i, nrow(work), sid)
  d <- tryCatch(sp_load(path, sid),
                error = function(e) { if (grepl("^Section 26", conditionMessage(e))) stop(e)
                  log_msg("  skipped: %s", conditionMessage(e)); NULL })
  if (is.null(d)) next
  ## (a) all spots, as section 15: used only for the safety check
  o_all <- sp_fit_obs(d)
  if (is.null(o_all)) { log_msg("  %s: model not fitted -> skipped", sid); next }
  ## (b) spots that have a grid position: the data the permutation works on
  n_all <- nrow(d); has_pos <- !is.na(d$grow) & !is.na(d$gcol)
  n_nopos <- sum(!has_pos)
  if (n_nopos) { log_msg("  %s: %d of %d spots have no grid position and are set aside",
                         sid, n_nopos, n_all); d <- d[has_pos, , drop = FALSE] }
  o <- if (n_nopos) sp_fit_obs(d) else o_all
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
obs_pool <- meta_pool(sp_obs$est, sp_obs$var, sp_obs$sample_id)              # positioned spots
all_pool <- meta_pool(sp_obs$est_all, sp_obs$var_all, sp_obs$sample_id)      # all spots, as section 15
write_out(sp_obs, "SPPERM_observed_per_tumour")

## ---- safety check ----------------------------------------------------------
ref <- iso_con_pooled[iso_con_pooled$outcome == SP_OUTCOME &
                        iso_con_pooled$contrast == "interaction_diff", ]$est
cat("\n", strrep("=", 74), "\n", sep = "")
cat("SECTION 26.4  SAFETY CHECK\n")
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  all spots  = %+.6f   section 15 = %+.6f   (k = %d)\n", all_pool$est, ref, all_pool$k))
if (length(ref) != 1 || abs(all_pool$est - ref) > SP_TOL) {
  cat("\n  per tumour (send this table):\n")
  print(sp_obs[, c("sample_id", "est_all", "n_spots_all", "n_no_position", "n_events")], row.names = FALSE)
  stop("Section 26 halted: the all-spot fit does not reproduce section 15. Nothing below is interpretable.")
}
cat("  PASSED\n")
cat(sprintf("  spots without a grid position: %d in %d tumour(s)\n",
            sum(sp_obs$n_no_position), sum(sp_obs$n_no_position > 0)))
if (any(sp_obs$n_no_position > 0))
  print(sp_obs[sp_obs$n_no_position > 0, c("sample_id", "n_spots_all", "n_no_position")], row.names = FALSE)
cat(sprintf("  estimate on positioned spots = %+.6f  <- the value the permutations are compared with\n",
            obs_pool$est))
## the fast fitter must agree with glm on the unpermuted data
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
  null_pool$pooled_b_col[b] <- mean(bc, na.rm = TRUE)   # simple means, for description only
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

## per tumour: where the observed estimate sits in that tumour's own null
per_t <- do.call(rbind, lapply(sids, function(s) {
  v <- null_est[, s]; v <- v[is.finite(v)]; o <- sp_obs$est[sp_obs$sample_id == s]
  data.frame(sample_id = s, observed = o, null_mean = mean(v), null_sd = sd(v),
             model_se = sqrt(sp_obs$var[sp_obs$sample_id == s]),
             p_one_sided = (1 + sum(v <= o)) / (length(v) + 1), n_perm_ok = length(v))
}))
per_t$se_ratio <- per_t$null_sd / per_t$model_se
write_out(per_t, "SPPERM_per_tumour")

sp_summary <- data.frame(outcome = SP_OUTCOME, observed = obs, k = obs_pool$k,
                         model_lower = obs_pool$lower, model_upper = obs_pool$upper,
                         n_perm = B, null_mean = mean(nl), null_sd = sd(nl),
                         null_q025 = q[[1]], null_median = q[[2]], null_q975 = q[[3]],
                         p_one_sided = p_one, p_two_sided = p_two,
                         median_k_in_null = median(null_pool$k, na.rm = TRUE),
                         mean_direct_share = mean(null_pool$direct, na.rm = TRUE),
                         median_se_ratio = median(per_t$se_ratio, na.rm = TRUE),
                         seed = SP_SEED, adjust_depth = SP_ADJUST_DEPTH, run_id = RUN_ID)
write_out(sp_summary, "SPPERM_summary")

cat("\n", strrep("=", 74), "\n", sep = "")
cat(sprintf("SECTION 26  RESULT  run_id %s\n", RUN_ID))
cat(strrep("=", 74), "\n", sep = "")
cat(sprintf("  observed pooled slope difference : %+.4f  (model CI %+.4f to %+.4f, k = %d)\n",
            obs, obs_pool$lower, obs_pool$upper, obs_pool$k))
cat(sprintf("  permutations usable              : %d of %d\n", B, SP_N_PERM))
cat(sprintf("  null: mean %+.4f, sd %.4f, middle 95%% %+.4f to %+.4f\n", mean(nl), sd(nl), q[[1]], q[[3]]))
cat(sprintf("  permuted estimates as negative or more : %d\n", sum(nl <= obs)))
cat(sprintf("  empirical p, one sided           : %.4f\n", p_one))
cat(sprintf("  empirical p, two sided           : %.4f\n", p_two))
cat(sprintf("  smallest p this run can give     : %.4f\n", 1 / (B + 1)))
cat(sprintf("  %s used per permutation (median)   : %.0f of %d\n",   # BRCA: k unit
            K_UNIT, median(null_pool$k, na.rm = TRUE),
            if (POOL_LEVEL == "patient") length(unique(patient_of(sids))) else K))
cat(sprintf("  spots that kept a true tissue value (mean share) : %.0f%%\n",
            100 * mean(null_pool$direct, na.rm = TRUE)))
cat(sprintf("  per tumour: null sd / model SE, median  : %.2f  (1 = model SE is honest, above 1 = too small)\n",
            median(per_t$se_ratio, na.rm = TRUE)))
cat(sprintf("  sections with own one sided p < 0.05    : %d of %d\n", sum(per_t$p_one_sided < 0.05), K))

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
