# Section map

Top level sections of `pipeline/p4ha1_ecm_pipeline_brca.R`, in the order they
appear in the file. Line numbers refer to that file. Section order and numbering
are the same as in the PDAC pipeline v4, so section N here is section N there.

As in PDAC, sections 13 and 14 are physically out of numeric order, with a shared
helper block for sections 13, 15, 16 and 17 between them. Sourcing the file top to
bottom still runs them in the intended order.

Sections changed for BRCA beyond paths and labels: 0 (config), 2 (pooling helpers),
4 (work list), 5 (cohort table), 21 (now the section level sensitivity) and 22
(leave one patient out, metastasis comparison).

| Line | Section |
|---|---|
| 77 | 0. CONFIG |
| 136 | 1. Gene panel + modules |
| 172 | 2. Helpers |
| 282 | 3. PER-SAMPLE analysis |
| 524 | 4. Build the work list, then RUN THE LOOP |
| 639 | 5. Save combined per-sample tables + run status |
| 656 | 6. CROSS-SAMPLE pooled forest |
| 684 | 7. Consistency boxplot: per-tumour P4HA1 <-> module correlations |
| 715 | 8. COHORT within-compartment analysis |
| 824 | 9. ENO1 SENSITIVITY CHECK   [single pass over the cohort] |
| 894 | 10. STRATIFIED MULTIVARIABLE MODEL |
| 1078 | 11. PROTEOSTASIS MODULES  (Valerio's Part 2, exploratory) |
| 1539 | 12. FIGURE EXPORT  --  journal style AND presentation style |
| 1996 | 14. DETECTION AUDIT FOR P4HA1, P4HA2 AND P4HA3 |
| 2485 | 13. FIGURE 4, REBUILT FROM ONE JOINT MODEL, WITH THE P VALUE ON IT |
| 2977 | 15. THE SPECIFICITY TEST: P4HA1 vs P4HA2 vs P4HA3 |
| 3387 | 16. THE RANDOM-GENE CONTROL |
| 3745 | 17. DEPTH ADJUSTED vs UNADJUSTED, AND RUN-TO-RUN REPRODUCIBILITY |
| 3923 | 18. THE METABOLIC AXIS: 2-OXOGLUTARATE SUPPLY INSTEAD OF HYPOXIA |
| 4555 | 19. IS THE METABOLIC MODULE INTERNALLY COHERENT? |
| 5159 | 20. DOES CELL COMPOSITION EXPLAIN THE RESULT? |
| 5552 | 21. SECTION-LEVEL SENSITIVITY   (BRCA: the PDAC section 21, turned around) |
| 5760 | 22. INFLUENCE OF SINGLE PATIENTS, AND THE METASTASIS |
| 5992 | 23. IS THE RESULT CARRIED BY ONE ER-STRESS GENE? |
| 6287 | 24. DOES THE ER-STRESS RESULT HOLD WITH INDEPENDENT UPR GENE SETS? |
| 6684 | 25. ARE THE PUBLISHED ER-STRESS SETS ANY MORE COHERENT THAN OURS? |
| 6939 | 26. SPATIAL PERMUTATION TEST |

## Subsections

| Line | Subsection |
|---|---|
| 1180 | 11.0  CIRCULARITY / OVERLAP AUDIT  -- halts the section on failure |
| 1217 | 11.1  score all modules per spot |
| 1273 | 11.2  fit every model for one sample, for BOTH ERstress variants |
| 1360 | 11.3  run across the cohort |
| 1433 | 11.4  pool and report  (figures live in section 12) |
| 1556 | 12.0  OBSERVED HYPOXIA SCORE RANGE   (runs once, then cached) |
| 1592 | 12.1  STYLE DEFINITIONS |
| 1692 | 12.2  SHARED HELPERS |
| 1751 | 12.3  READ THE POOLED RESULTS  (no re-analysis below this line) |
| 1770 | 12.4  RENDER EVERY FIGURE IN ONE STYLE |
| 1980 | 12.5  RUN BOTH STYLES |
| 2037 | 14.0  SETTINGS |
| 2059 | 14.1  CIRCULARITY RE-AUDIT FOR THE NEW OUTCOMES |
| 2087 | 14.2  COUNT, ONE TUMOUR AT A TIME |
| 2201 | 14.3  REPORT |
| 2298 | 14.4  VERDICT |
| 2523 | 13.0  SETTINGS |
| 2546 | 13.0a  RUN STATE GATE   *** THE FIX FOR FAULT 1 *** |
| 2559 | 13.0b  RUN MANIFEST   *** PART OF THE FIX FOR FAULT 2 *** |
| 2568 | 13.1  HELPERS |
| 2645 | 13.2  FIT ONE TUMOUR |
| 2716 | 13.3  RUN ACROSS THE COHORT |
| 2768 | 13.4  POOL THE COEFFICIENTS AND THE SCALAR CONTRASTS |
| 2781 | 13.5  POOL THE CURVES, one meta-analysis per grid point per series |
| 2819 | 13.6  CONSOLE REPORT |
| 2891 | 13.7  FIGURE 4, REBUILT |
| 2999 | 15.0  SETTINGS |
| 3021 | 15.1  CIRCULARITY GATE  -- hard stop, not a warning |
| 3038 | 15.2  SCORE AND STANDARDISE ONCE PER TUMOUR |
| 3101 | 15.3  FIT THE SAME MODEL FOR ONE OUTCOME |
| 3179 | 15.4  RUN |
| 3221 | 15.5  POOL |
| 3239 | 15.6  REPORT |
| 3320 | 15.7  FIGURE 7 |
| 3410 | 16.0  SETTINGS |
| 3433 | 16.1  CHOOSE THE GENES |
| 3487 | 16.2  FIT, ONE TUMOUR AT A TIME |
| 3586 | 16.3  RUN |
| 3604 | 16.4  POOL EACH GENE ACROSS TUMOURS |
| 3631 | 16.5  WHERE DO THE THREE ISOFORMS SIT? |
| 3696 | 16.6  FIGURE 8 |
| 4026 | 18.0  GENE SETS AND SETTINGS |
| 4091 | 18.1  GENE-SET AUDIT, halts on failure |
| 4123 | 18.2  SCORE, STANDARDISE, AND MEASURE THE COLLINEARITY GATE |
| 4197 | 18.3  FIT ONE TUMOUR, ONE METABOLIC VARIANT |
| 4260 | 18.4  RUN THE COHORT |
| 4299 | 18.5  THE COLLINEARITY GATE, reported before any result |
| 4355 | 18.6  POOL AND REPORT |
| 4490 | 18.7  FIGURE 9 |
| 4649 | 19.0b  CANDIDATE METABOLIC GENE SETS |
| 4757 | 19.0  SETTINGS |
| 4776 | 19.1  BUILD THE POOL OF GENES THE NULL SETS ARE DRAWN FROM |
| 4830 | 19.2  THE TWO COHERENCE MEASURES |
| 4861 | 19.3  RUN, ONE TUMOUR AT A TIME |
| 4935 | 19.4  SUMMARISE |
| 4972 | 19.5  VERDICT |
| 5102 | 19.6  FIGURE 10 |
| 5212 | 20.0  SETTINGS |
| 5237 | 20.1  BUILD THE CLEAN STROMA AND IMMUNE GENE SETS |
| 5283 | 20.2  SCORE, STANDARDISE AND VALIDATE, ONE TUMOUR AT A TIME |
| 5351 | 20.3  FIT ONE MODEL FOR ONE OUTCOME |
| 5413 | 20.4  RUN |
| 5445 | 20.5  POOL |
| 5456 | 20.6  SAFETY CHECK: M0 MUST REPRODUCE SECTION 15, FOR EVERY OUTCOME |
| 5481 | 20.7  REPORT |
| 5594 | 21.0  READ THE PER-SECTION CONTRASTS AND THE SECTION -> PATIENT MAP |
| 5627 | 21.1  PATIENT ID FROM THE COHORT TABLE (not parsed from the section ID: |
| 5642 | 21.2  POOLING |
| 5674 | 21.3  CHECK: REPRODUCE THE PIPELINE'S POOLED RESULTS |
| 5699 | 21.4  PATIENT LEVEL (primary) vs SECTION LEVEL (PDAC style) |
| 5803 | 22.0  READ |
| 5839 | 22.1  POOLING  (same as the pipeline's patient-level meta_pool; verified below) |
| 5874 | 22.2  CHECK: does this pooling reproduce the pipeline's? |
| 5897 | 22.3  LEAVE ONE PATIENT OUT, WITHOUT THE METASTASIS, SIGN CONCORDANCE |
| 5948 | 22.4  REPORT |
| 6041 | 23.0  SETTINGS |
| 6078 | 23.1  ONE TUMOUR: SCORE THE FIXED MODULES ONCE, THEN EVERY ER VARIANT |
| 6135 | 23.2  FIT ONE VARIANT |
| 6170 | 23.3  RUN |
| 6198 | 23.4  SAFETY CHECK: THE FULL MODULE MUST REPRODUCE SECTION 15 |
| 6212 | 23.5  REPORT |
| 6346 | 24.0  SETTINGS |
| 6370 | 24.1  READ THE CANDIDATE SETS FROM THE TSVs |
| 6406 | 24.2  AUDIT AND CLEAN |
| 6464 | 24.3  ONE TUMOUR: SCORE EVERYTHING ONCE |
| 6519 | 24.4  FIT ONE VARIANT |
| 6551 | 24.5  RUN |
| 6582 | 24.6  SAFETY CHECK |
| 6596 | 24.7  REPORT |
| 6726 | 25.0  SETTINGS |
| 6746 | 25.1  NULL POOL, MATCHED ON DETECTION RATE |
| 6785 | 25.2  THE TWO MEASURES  (same as section 19) |
| 6808 | 25.3  RUN |
| 6871 | 25.4  REPORT |
| 6994 | 26.0  SETTINGS |
| 7023 | 26.1  ONE TUMOUR: SCORES, OUTCOME AND GRID POSITION OF EVERY SPOT |
| 7073 | 26.2  FITS |
| 7107 | 26.3  THE GRID OF ONE TUMOUR AND THE MOVE |
| 7144 | 26.4  LOAD EVERY TUMOUR ONCE AND FIT THE UNPERMUTED MODEL |
| 7209 | 26.5  PERMUTATIONS |
| 7244 | 26.6  REPORT |
