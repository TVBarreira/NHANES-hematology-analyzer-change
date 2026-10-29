###############################################################################
## M01_build_cbc_frame.R
##
## Builds the analytic frame for "Hematology Analyzer Change and Blood-Count
## Measures of Biological Aging and Inflammation in NHANES" (JAMA Research
## Letter). Downloads NHANES demographic, complete-blood-count and biochemistry
## files for the 11 cycles 1999-2000 to 2021-2023, applies the exclusions,
## defines the adolescent referent (12-19 years, nonpregnant), computes the
## composite indices (RAR, NLR, PLR, SII, SIRI), and saves data/m_cbc_frame.rds.
## Raw downloads are cached per cycle under data/.
##
## Notes that matter for reproduction:
##   * NHANES laboratory file names change across cycles (LAB25, L25_B, CBC_D...;
##     LAB18, L40_B, BIOPRO_D...). Each target is resolved against candidate
##     names; see grab() and the CYCLES table.
##   * Biochemistry variable names also change; each analyte is resolved against
##     candidate names, each with its unit multiplier.
##   * Alkaline phosphatase is kept both as released (alp_raw) and after the
##     NCHS bridging equations (alp), so the 2017-2018 external validation
##     compares the observed step against the NCHS split-sample prediction on the
##     released scale.
##   * The 2-year MEC weight (WTMEC2YR) is used for all cycles.
###############################################################################

rm(list = ls())

.nosleep_handle <- NULL
if (requireNamespace("NoSleepR", quietly = TRUE))
  .nosleep_handle <- try(NoSleepR::nosleep_on(keep_display = TRUE), silent = TRUE)

library(nhanesA); library(dplyr); library(tibble)
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

## Project folder: set the environment variable NHANES_BIOAGE_DIR, or run from
## the project folder (it must contain data/ and results/ subfolders).
base_dir <- Sys.getenv("NHANES_BIOAGE_DIR", unset = getwd())
if (!dir.exists(base_dir)) stop("Project folder not found: ", base_dir)
data_dir <- file.path(base_dir, "data"); results_dir <- file.path(base_dir, "results")
fig_dir <- file.path(results_dir, "figures")
for (.d in c(data_dir, results_dir, fig_dir))
  if (!dir.exists(.d)) dir.create(.d, recursive = TRUE)

## ---------------------------------------------------------------------
## Cycles, with EXPLICIT file names per cycle.
##
## The 2017-March 2020 pre-pandemic file (tag "P") is EXCLUDED: it re-releases
## 2017-2018 participants under new sequence numbers and would duplicate them.
##
## If nhanesA reports a name as unavailable, resolve it with
##   nhanesSearchTableNames("CBC"); nhanesSearchTableNames("BIOPRO")
## and correct the map here rather than letting the cycle drop.
## ---------------------------------------------------------------------
CYCLES <- tribble(
  ~key, ~year, ~label,        ~demo,     ~cbc,      ~bio,
  "A",  1999,  "1999-2000",   "DEMO",    "LAB25",   "LAB18",
  "B",  2001,  "2001-2002",   "DEMO_B",  "L25_B",   "L40_B",
  "C",  2003,  "2003-2004",   "DEMO_C",  "L25_C",   "L40_C",
  "D",  2005,  "2005-2006",   "DEMO_D",  "CBC_D",   "BIOPRO_D",
  "E",  2007,  "2007-2008",   "DEMO_E",  "CBC_E",   "BIOPRO_E",
  "F",  2009,  "2009-2010",   "DEMO_F",  "CBC_F",   "BIOPRO_F",
  "G",  2011,  "2011-2012",   "DEMO_G",  "CBC_G",   "BIOPRO_G",
  "H",  2013,  "2013-2014",   "DEMO_H",  "CBC_H",   "BIOPRO_H",
  "I",  2015,  "2015-2016",   "DEMO_I",  "CBC_I",   "BIOPRO_I",
  "J",  2017,  "2017-2018",   "DEMO_J",  "CBC_J",   "BIOPRO_J",
  "L",  2022,  "2021-2023",   "DEMO_L",  "CBC_L",   "BIOPRO_L")

## NEED = abort if absent. OPT = NA-fill and log.
DEMO_NEED <- c("SEQN", "RIDAGEYR", "RIAGENDR", "SDMVPSU", "SDMVSTRA", "RIDSTATR")
DEMO_OPT  <- c("RIDRETH1", "RIDEXPRG")
CBC_NEED  <- c("SEQN", "LBXRDW", "LBXWBCSI")
NEW_CBC   <- c("LBXHGB", "LBXMC", "LBXRBCSI")   # hemoglobin, MCHC, RBC count
CBC_OPT   <- c("LBXMCVSI", "LBXLYPCT", "LBXNEPCT", "LBXMOPCT",
               "LBXPLTSI", "LBXHCT", "LBXMCHSI",
               "LBXHGB", "LBXMC", "LBXRBCSI")   # 2026-10-18: Hgb, MCHC, RBC - NCHS named
                                                #   HCT, MCH and MCHC at the 2013 boundary
## BIOCHEMISTRY NAMES ARE NOT STABLE ACROSS CYCLES. L40_B (2001-2002) does not
## carry LBXSAPSI or LBXSCR. Each target is therefore resolved against
## candidate names, each carrying the multiplier that
## puts it in the conventional unit. The SI-named variants are listed because
## silently picking one up unconverted would corrupt PhenoAge with no error at
## all: creatinine in umol/L is ~88x its mg/dL value.
BIO_ALIAS <- list(
  albumin = list(c("LBXSAL", 1), c("LBDSAL", 1),
                 c("LBXSALSI", 0.1), c("LBDSALSI", 0.1)),            # g/L -> g/dL
  alp     = list(c("LBXSAPSI", 1), c("LBDSAPSI", 1),
                 c("LBXSAP", 1), c("LBDSAP", 1)),                    # U/L throughout
  creat   = list(c("LBXSCR", 1), c("LBDSCR", 1),
                 c("LBXSCRSI", 1/88.4017), c("LBDSCRSI", 1/88.4017)),# umol/L -> mg/dL
  glucose = list(c("LBXSGL", 1), c("LBDSGL", 1),
                 c("LBXSGLSI", 1/0.0555)),                           # mmol/L -> mg/dL
  totchol = list(c("LBXSCH", 1), c("LBDSCH", 1),
                 c("LBXSCHSI", 38.67)),                              # mmol/L -> mg/dL
  bun     = list(c("LBXSBU", 1), c("LBDSBU", 1),
                 c("LBXSBUSI", 2.8)),                                # mmol/L -> mg/dL
  uap     = list(c("LBXSUA", 1), c("LBDSUA", 1),
                 c("LBXSUASI", 1/59.48)))                            # umol/L -> mg/dL

## Returns SEQN plus the seven targets, already unit-converted, and records
## which source name was used. Any target it cannot resolve triggers a full
## dump of the file's column names, so the real name is visible in one run
## instead of after another round trip.
resolve_bio <- function(df, cyc) {
  if (is.null(df)) return(list(d = NULL, used = character(0)))
  out <- data.frame(SEQN = df$SEQN); used <- character(0)
  for (tgt in names(BIO_ALIAS)) {
    hit <- NULL
    for (cand in BIO_ALIAS[[tgt]])
      if (cand[1] %in% names(df)) { hit <- cand; break }
    if (is.null(hit)) {
      message("cycle ", cyc, ": cannot resolve ", tgt,
              ". Columns present in this file:\n  ",
              paste(names(df), collapse = ", "))
      out[[tgt]] <- NA_real_; used <- c(used, paste0(tgt, "=UNRESOLVED"))
    } else {
      mult <- as.numeric(hit[2])
      out[[tgt]] <- as.numeric(df[[hit[1]]]) * mult
      used <- c(used, paste0(tgt, "=", hit[1],
                             if (mult != 1) sprintf("*%.5g", mult) else ""))
    }
  }
  list(d = out, used = used)
}
## MEC weight: preferred name first. DEMO_L is expected to carry WTMEC2YR;
## WTMECPRP belongs to the excluded 2017-March 2020 file and is a fallback only.
WT_CANDIDATES <- c("WTMEC2YR", "WTMECPRP", "WTMEC4YR")

grab <- function(nm) {
  out <- try(nhanes(nm, translated = FALSE), silent = TRUE)
  if (inherits(out, "try-error") || is.null(out) || !nrow(out)) return(NULL)
  message("    ", nm, ": ", nrow(out), " rows, ", ncol(out), " cols")
  out
}
take <- function(df, need, opt, what, cyc) {
  if (is.null(df)) return(list(d = NULL, missing = character(0)))
  absent_need <- setdiff(need, names(df))
  if (length(absent_need))
    stop("cycle ", cyc, ", ", what, ": required column(s) absent: ",
         paste(absent_need, collapse = ", "),
         "\n  Inspect with nhanesTableVars() and update the NEED/OPT lists.")
  absent_opt <- setdiff(opt, names(df))
  for (v in absent_opt) df[[v]] <- NA_real_
  list(d = df %>% select(all_of(c(need, opt))), missing = absent_opt)
}
no_suffixes <- function(df, where) {
  bad <- grep("\\.(x|y)$", names(df), value = TRUE)
  if (length(bad)) stop("join created suffixed columns at ", where, ": ",
                        paste(bad, collapse = ", "),
                        "\n  A column was joined that already existed.")
  invisible(TRUE)
}

## ---------------------------------------------------------------------
## Download and assemble
## ---------------------------------------------------------------------
frames <- list(); blog <- list()
for (i in seq_len(nrow(CYCLES))) {
  cy <- CYCLES[i, ]
  ## Cache of the selected columns for this cycle. An older cache lacking the
  ## hemoglobin, MCHC and RBC columns is upgraded by fetching only the CBC table.
  cache <- file.path(data_dir, paste0("m_raw_v4_", cy$key, ".rds"))
  if (file.exists(cache)) {
    message("cycle ", cy$label, ": cached")
    frames[[cy$key]] <- readRDS(cache)
    blog[[cy$key]] <- data.frame(cycle = cy$label, source = "cache",
                                 n = nrow(frames[[cy$key]]),
                                 wt_var = NA_character_, missing = NA_character_)
    next
  }
  ## Incremental upgrade. An older cache holds everything except Hgb, MCHC and RBC,
  ## so fetch ONLY this cycle's CBC table and join those three columns on. The
  ## demographics, weights and biochemistry already in the cache are untouched,
  ## so every number built from them is identical. If anything fails, this
  ## cycle alone falls back to the full download below.
  cache3 <- file.path(data_dir, paste0("m_raw_v3_", cy$key, ".rds"))
  if (file.exists(cache3)) {
    up <- tryCatch({
      d3 <- readRDS(cache3); n0 <- nrow(d3)
      rc3 <- grab(cy$cbc)
      if (is.null(rc3) || !"SEQN" %in% names(rc3)) stop("CBC table unavailable")
      add <- rc3[, intersect(c("SEQN", NEW_CBC), names(rc3)), drop = FALSE]
      for (v in setdiff(NEW_CBC, names(add))) add[[v]] <- NA_real_
      d3 <- d3[, setdiff(names(d3), NEW_CBC), drop = FALSE]
      d3 <- left_join(d3, add, by = "SEQN")
      no_suffixes(d3, paste(cy$label, "v3 cache + new CBC columns"))
      if (nrow(d3) != n0) stop("row count changed (", n0, " -> ", nrow(d3), ")")
      saveRDS(d3, cache); d3
    }, error = function(e) {
      message("    incremental upgrade failed (", conditionMessage(e), ") - full download for this cycle")
      NULL })
    if (!is.null(up)) {
      message("cycle ", cy$label, ": v3 cache + CBC table only")
      frames[[cy$key]] <- up
      blog[[cy$key]] <- data.frame(cycle = cy$label, source = "v3 cache + CBC table",
                                   n = nrow(up), wt_var = NA_character_, missing = NA_character_)
      next
    }
  }
  message("cycle ", cy$label, " ...")
  rd <- grab(cy$demo); rc <- grab(cy$cbc); rb <- grab(cy$bio)
  if (is.null(rd)) stop("cycle ", cy$label, ": demographics file ", cy$demo,
                        " unavailable. Check the name in CYCLES.")
  if (is.null(rc)) stop("cycle ", cy$label, ": CBC file ", cy$cbc,
                        " unavailable. Check the name in CYCLES -- the early ",
                        "cycles are LAB25 / L25_B / L25_C, not CBC_<tag>.")

  wt_var <- WT_CANDIDATES[WT_CANDIDATES %in% names(rd)][1]
  if (is.na(wt_var)) stop("cycle ", cy$label, ": no MEC weight found among ",
                          paste(WT_CANDIDATES, collapse = "/"),
                          ". Columns present: ", paste(names(rd), collapse = ", "))
  if (wt_var != "WTMEC2YR")
    warning("cycle ", cy$label, ": using ", wt_var,
            " as the MEC weight. Confirm this is the 2-year MEC weight before ",
            "pooling -- see the pooling note in M02.")
  rd$WTMEC <- rd[[wt_var]]

  td <- take(rd, c(DEMO_NEED, "WTMEC"), DEMO_OPT, "DEMO", cy$label)
  tc <- take(rc, CBC_NEED, CBC_OPT, "CBC", cy$label)
  tb <- resolve_bio(rb, cy$label)

  d <- td$d %>% left_join(tc$d, by = "SEQN"); no_suffixes(d, paste(cy$label, "DEMO+CBC"))
  if (!is.null(tb$d)) {
    d <- d %>% left_join(tb$d, by = "SEQN"); no_suffixes(d, paste(cy$label, "+BIOPRO"))
  } else for (v in names(BIO_ALIAS)) d[[v]] <- NA_real_
  message("    biochem: ", paste(tb$used, collapse = "  "))

  ## a cycle contributing no usable haematology is a build failure, not a datum
  n_cbc <- sum(!is.na(d$LBXRDW) & !is.na(d$LBXWBCSI) & d$WTMEC > 0, na.rm = TRUE)
  if (n_cbc < 500)
    stop("cycle ", cy$label, ": only ", n_cbc,
         " participants with RDW, WBC and a positive MEC weight. ",
         "Expected several thousand -- the merge or the file name is wrong.")

  d$year <- cy$year; d$cycle <- cy$label; d$wt_var <- wt_var
  saveRDS(d, cache); frames[[cy$key]] <- d
  blog[[cy$key]] <- data.frame(
    cycle = cy$label, source = paste(cy$demo, cy$cbc, cy$bio, sep = " / "),
    n = nrow(d), wt_var = wt_var,
    missing = paste(c(td$missing, tc$missing), collapse = ";"),
    biochem = paste(tb$used, collapse = " "))
  message("    kept ", nrow(d), " rows; ", n_cbc, " with a usable CBC")
}
raw <- bind_rows(frames)
stopifnot(length(unique(raw$cycle)) == nrow(CYCLES))
stopifnot(!anyDuplicated(raw$SEQN))          # the "P" file is the usual culprit
message("assembled: ", nrow(raw), " rows across ", nrow(CYCLES), " cycles")
write.csv(bind_rows(blog), file.path(results_dir, "m_build_log.csv"), row.names = FALSE)

## ---------------------------------------------------------------------
## Rename, screen implausible values, apply the CDC chemistry bridging
##
## RIDEXPRG is asked only of women in the childbearing age range, so NA is the
## norm for men, children and older women. Coding NA and "cannot ascertain" as
## not-pregnant is therefore the intended behaviour, not a lapse -- but the
## per-cycle count is logged so a cycle that omits the variable entirely is
## visible rather than silently treated as a cycle with no pregnancies.
## ---------------------------------------------------------------------
d <- raw %>%
  transmute(
    seqn = SEQN, year, cycle,
    age = RIDAGEYR, sex = RIAGENDR, race = RIDRETH1,
    pregnant = ifelse(!is.na(RIDEXPRG) & RIDEXPRG == 1, 1L, 0L),
    preg_known = !is.na(RIDEXPRG),
    exam_status = RIDSTATR,
    psu = SDMVPSU, stratum = SDMVSTRA, wt2yr = WTMEC, wt_var,
    rdw   = ifelse(LBXRDW   <= 10   | LBXRDW   >= 40,  NA, LBXRDW),
    wbc   = ifelse(LBXWBCSI <= 0.05 | LBXWBCSI >= 400, NA, LBXWBCSI),
    mcv   = ifelse(LBXMCVSI <= 50   | LBXMCVSI >= 150, NA, LBXMCVSI),
    lymph = ifelse(LBXLYPCT <  0    | LBXLYPCT >  100, NA, LBXLYPCT),
    neut  = ifelse(LBXNEPCT <  0    | LBXNEPCT >  100, NA, LBXNEPCT),
    mono  = ifelse(LBXMOPCT <  0    | LBXMOPCT >  100, NA, LBXMOPCT),
    plt   = ifelse(LBXPLTSI <= 0    | LBXPLTSI >= 2000, NA, LBXPLTSI),
    hct = LBXHCT, mch = LBXMCHSI, hgb = LBXHGB, mchc = LBXMC, rbc = LBXRBCSI,
    albumin, alp, creat, glucose, totchol, bun, uap,
    ## ALP kept in its released state as well. The 2017-2018 external
    ## validation compares the observed step against the NCHS split-sample
    ## prediction; it must therefore be estimated on the released scale, before
    ## that equation is applied.
    alp_raw = alp)

stopifnot(all(c("albumin","alp","creat","glucose","alp_raw") %in% names(d)))
cat("\nPhenoAge input completeness by cycle (all 8 inputs present):\n")
print(tapply(with(d, !is.na(albumin) & !is.na(alp) & !is.na(creat) &
                     !is.na(glucose) & !is.na(lymph) & !is.na(mcv) &
                     !is.na(rdw) & !is.na(wbc)), d$cycle, sum))
cat("\nAny cycle at 0 above means a name is still unresolved -- see the dump.\n")
cat("\nRIDEXPRG availability by cycle (non-missing count):\n")
print(tapply(d$preg_known, d$cycle, sum))

## NCHS bridging equations as a two-stage chain. Stage 1 puts the Cobas 8000
## (2021-2023) onto the Cobas 6000 scale; stage 2 puts the Cobas 6000
## (2017-2018 and, after stage 1, 2021-2023) onto the DxC 660i scale used
## earlier. Where NCHS states adjustment is not recommended, none is applied --
## hence glucose is bridged only at 2021-2023, and albumin, creatinine and
## total cholesterol are not bridged at 2021-2023.
d$harmonised <- FALSE
i22 <- which(d$year == 2022)
if (length(i22)) {
  d$alp[i22]     <- 0.6098  + 0.9041 * d$alp[i22]
  d$glucose[i22] <- -0.7821 + 0.9669 * d$glucose[i22]
}
ico <- which(d$year %in% c(2017, 2022))
if (length(ico)) {
  d$albumin[ico] <- 1.044 * d$albumin[ico] + 0.01128
  d$alp[ico]     <- 10^(1.001 * log10(d$alp[ico]) - 0.04294)
  d$creat[ico]   <- 1.051 * d$creat[ico] - 0.06945
  d$totchol[ico] <- 1.046 * d$totchol[ico] - 2.203
  d$harmonised[ico] <- TRUE
}
message("chemistry harmonised for ", length(ico), " participants")

## ---------------------------------------------------------------------
## Composite indices
##
## NHANES reports the differential as PERCENTAGES; absolute counts are
## percentage/100 x WBC. Which indices retain WBC after simplification is the
## paper's organising point:
##   NLR  = neut% / lymph%                        WBC cancels
##   SII  = plt x neut_ct / lymph_ct
##        = plt x neut% / lymph%                  WBC cancels
##   PLR  = plt / lymph_ct                        WBC retained (denominator)
##   SIRI = neut_ct x mono_ct / lymph_ct
##        proportional to WBC x neut% x mono%/lymph%   WBC retained
##   RAR  = rdw / albumin                         RDW, undiluted
## ---------------------------------------------------------------------
d <- d %>% mutate(
  neut_ct  = neut/100  * wbc,
  lymph_ct = lymph/100 * wbc,
  mono_ct  = mono/100  * wbc,
  RAR  = rdw / albumin,
  NLR  = ifelse(lymph > 0, neut / lymph, NA_real_),
  SIRI = ifelse(lymph_ct > 0, neut_ct * mono_ct / lymph_ct, NA_real_),
  SII  = ifelse(lymph_ct > 0, plt * neut_ct / lymph_ct, NA_real_),
  PLR  = ifelse(lymph_ct > 0, plt / lymph_ct, NA_real_))

## ---------------------------------------------------------------------
## Analysis frame, then the ledger. The ledger's final row is nrow(ana) by
## construction so a STROBE flow diagram drawn from it reconciles.
## ---------------------------------------------------------------------
ana <- d %>% filter(!is.na(exam_status), exam_status == 2,
                    !is.na(rdw), !is.na(wbc),
                    pregnant == 0, !is.na(wt2yr), wt2yr > 0)

ledger <- tibble(
  step = c("released in these cycles",
           "MEC-examined",
           "  and with a valid RDW",
           "  and with a valid WBC",
           "  and non-pregnant",
           "  and with a positive MEC weight  [= ANALYSIS FRAME]",
           "of the frame: with the full differential",
           "of the frame: with albumin, so RAR computable",
           "of the frame: with platelets, so SII computable",
           "of the frame: adults 20+ (the clinical frame, for comparison)"),
  n = c(nrow(d),
        sum(!is.na(d$exam_status) & d$exam_status == 2),
        sum(!is.na(d$exam_status) & d$exam_status == 2 & !is.na(d$rdw)),
        sum(!is.na(d$exam_status) & d$exam_status == 2 & !is.na(d$rdw) & !is.na(d$wbc)),
        sum(!is.na(d$exam_status) & d$exam_status == 2 & !is.na(d$rdw) &
              !is.na(d$wbc) & d$pregnant == 0),
        nrow(ana),
        sum(!is.na(ana$lymph) & !is.na(ana$neut) & !is.na(ana$mono)),
        sum(!is.na(ana$RAR)),
        sum(!is.na(ana$SII)),
        sum(ana$age >= 20)))
stopifnot(ledger$n[6] == nrow(ana))
write.csv(ledger, file.path(results_dir, "m_sample_ledger.csv"), row.names = FALSE)
cat("\n===== SAMPLE LEDGER =====\n")
print(as.data.frame(ledger), row.names = FALSE)
cat("\nCompare the last row against row 6. The difference is what the clinical\n")
cat("paper's adult restriction was costing this analysis.\n")

## Pooled MEC weight: each 2-year weight divided by the number of cycles
## pooled. NOTE two things before quoting any weighted total. (a) 1999-2002 has
## 4-year weights (WTMEC4YR) that NCHS prefers for that span; using 2-year
## weights divided by 11 is the standard pooling workaround but is an
## approximation there. (b) 2021-2023 spans about three years, so equal
## weighting of cycles is likewise approximate. Neither affects a WITHIN-cycle
## mean or a difference between two adjacent cycles, which is what this paper
## estimates -- but it does affect any pooled prevalence.
ana$wt <- ana$wt2yr / length(unique(ana$cycle))

ana <- ana %>% mutate(ref_adol  = age >= 12 & age <= 19,
                      ref_young = age >= 20 & age <= 34)
cat("\nadolescent referent:", sum(ana$ref_adol),
    "| young adult referent:", sum(ana$ref_young), "\n")
for (r in c("ref_adol", "ref_young")) {
  cat("  ", r, "by cycle:\n"); print(tapply(ana[[r]], ana$cycle, sum))
}
if (min(tapply(ana$ref_adol, ana$cycle, sum)) < 300)
  warning("a cycle has fewer than 300 adolescents; the adolescent-referent ",
          "step will be imprecise in that cycle.")

saveRDS(ana, file.path(data_dir, "m_cbc_frame.rds"))
cat("\nanalysis frame:", nrow(ana), "participants, ages",
    min(ana$age, na.rm = TRUE), "to", max(ana$age, na.rm = TRUE), "\n")
cat("\nSaved: data/m_cbc_frame.rds, results/m_sample_ledger.csv,\n")
cat("       results/m_build_log.csv\n")
cat("NEXT: run M02_measurement_analyses.R\n")

if (!is.null(.nosleep_handle) && !inherits(.nosleep_handle, "try-error"))
  try(NoSleepR::nosleep_off(.nosleep_handle), silent = TRUE)
message("Done: CBC frame built.")
