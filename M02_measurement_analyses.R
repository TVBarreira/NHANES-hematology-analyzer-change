###############################################################################
## M02_measurement_analyses.R
##
## All estimation for "Hematology Analyzer Change and Blood-Count Measures of
## Biological Aging and Inflammation in NHANES" (JAMA Research Letter).
## Run M01_build_cbc_frame.R first.
##
## Sections
##   1  Analytic frame and referents
##   2  Change in each blood-count and chemistry analyte between adjacent cycles,
##      in the adolescent referent (2b: all 10 boundaries; 2c: change by
##      percentile, additive vs proportional; 2d: correction factors)
##   3  Composite indices (RAR, NLR, PLR, SII, SIRI) at the analyzer boundary
##   4  Phenotypic age (BioAge phenoage_calc trained in NHANES III, CRP omitted)
##      and its change at the analyzer boundary, overall and by age group
##   5  Counterfactual correction of RDW and WBC (and chemistry placebo tests);
##      5b: composite measures observed and corrected, adults, every boundary
##   6  Pooled-vs-within-cycle quartile reassignment
##   7  Figure panel series and the token-to-value file m_placeholders.csv
##
## Outputs (results/)
##   m_analyte_steps.csv, m_all_boundary_steps.csv, m_all_boundary_specificity.csv,
##   m_correction_form_quantiles.csv, m_correction_form_call.csv,
##   m_correction_factors.csv, m_index_steps.csv, m_pheno_steps.csv,
##   m_counterfactual.csv, m_derived_boundary_steps.csv, m_quartile_shift.csv,
##   m_figure_panels.csv, m_placeholders.csv
##
## Design notes
##   * Steps are estimated as the difference between adjacent cycles, never from
##     a trend model with a step term (which attributes any departure from the
##     assumed trend to the step).
##   * Phenotypic age is TRAINED in NHANES III with phenoage_calc(fit = NULL)
##     and projected with that fit; published coefficients with the CRP term
##     deleted are a different, uncalibrated measure (about 0.8 y higher).
##   * The phenotypic age step is estimated on phenotypic age itself with age as
##     a covariate, not on a pooled acceleration residual, whose level in any
##     cycle depends on the composition of the pool.
##   * Chemistry analytes, for which no platform changed at this boundary, are
##     corrected in the same way as placebo tests; the largest of them bounds the
##     method's noise floor.
###############################################################################

rm(list = ls())

.nosleep_handle <- NULL
if (requireNamespace("NoSleepR", quietly = TRUE))
  .nosleep_handle <- try(NoSleepR::nosleep_on(keep_display = TRUE), silent = TRUE)

library(dplyr); library(survey); library(tibble); library(tidyr)
library(BioAge)   # phenoage_calc + bundled NHANES3; must match 01_build_data
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
options(survey.lonely.psu = "adjust")

## Project folder: set the environment variable NHANES_BIOAGE_DIR, or run from
## the project folder (it must contain data/ and results/ subfolders).
base_dir <- Sys.getenv("NHANES_BIOAGE_DIR", unset = getwd())
if (!dir.exists(base_dir)) stop("Project folder not found: ", base_dir)
data_dir <- file.path(base_dir, "data"); results_dir <- file.path(base_dir, "results")
frame_path <- file.path(data_dir, "m_cbc_frame.rds")
if (!file.exists(frame_path)) stop("Run M01_build_cbc_frame first: ", frame_path)
ana <- readRDS(frame_path)
message("frame loaded: ", nrow(ana), " participants")

PH <- list()                                  # placeholder collector
put <- function(token, value, fmt = "%.2f") {
  PH[[token]] <<- if (is.character(value)) value else sprintf(fmt, value); invisible()
}
fmt_p <- function(p) if (is.na(p)) "NA" else if (p < .001) "<.001" else
  sub("^0", "", sprintf("%.2f", p))

## =====================================================================
## 0. Sample ledger tokens
## =====================================================================
led <- read.csv(file.path(results_dir, "m_sample_ledger.csv"))
gl <- function(pat) led$n[grep(pat, led$step, fixed = TRUE)][1]
put("N_RELEASED", format(gl("released in these cycles"), big.mark = ","))
put("N_CBC_NP",   format(gl("[= ANALYSIS FRAME]"),        big.mark = ","))
put("N_RAR",      format(gl("RAR computable"),            big.mark = ","))
put("N_SII",      format(gl("SII computable"),            big.mark = ","))

## =====================================================================
## 1. Designs
## =====================================================================
mk_design <- function(df) svydesign(id = ~psu, strata = ~stratum, weights = ~wt,
                                    nest = TRUE, data = df)

## Step between two adjacent cycles, within an optional logical subset column.
## Returns the difference, its CI, and P, from a single indicator regression --
## so the SE respects the survey design rather than being pooled by hand.
step_est <- function(df, yvar, y0, y1, subset_col = NULL, adjust_age = FALSE) {
  s <- df[df$year %in% c(y0, y1), , drop = FALSE]
  if (!is.null(subset_col)) s <- s[which(s[[subset_col]]), , drop = FALSE]
  s$.y <- s[[yvar]]; s <- s[!is.na(s$.y) & !is.na(s$wt) & s$wt > 0, , drop = FALSE]
  out <- list(n = nrow(s), est = NA_real_, lo = NA_real_, hi = NA_real_,
              p = NA_real_, m0 = NA_real_, m1 = NA_real_, sd = NA_real_)
  if (nrow(s) < 200 || length(unique(s$year)) < 2) return(out)
  s$.post <- as.integer(s$year == y1)
  d <- try(mk_design(s), silent = TRUE); if (inherits(d, "try-error")) return(out)
  f <- if (adjust_age) .y ~ .post + age else .y ~ .post
  m <- try(svyglm(f, design = d), silent = TRUE)
  if (inherits(m, "try-error")) return(out)
  co <- summary(m)$coefficients
  out$est <- co[".post", 1]; out$p <- co[".post", ncol(co)]
  ## CI built from the design df rather than confint(), whose parm matching and
  ## ddf argument name both vary across survey package versions
  tcrit <- qt(0.975, degf(d))
  out$lo <- out$est - tcrit * co[".post", 2]
  out$hi <- out$est + tcrit * co[".post", 2]
  mm <- try(svyby(~.y, ~year, d, svymean, na.rm = TRUE), silent = TRUE)
  if (!inherits(mm, "try-error")) { out$m0 <- mm[[2]][1]; out$m1 <- mm[[2]][2] }
  out$sd <- sd(df[[yvar]], na.rm = TRUE)
  out
}

## =====================================================================
## 2. Analyte steps at the haematology boundary, both referents
##    plus the external-validation step: ALP at the 2017-2018 chemistry
##    boundary, where NCHS's split-sample equation gives an independent value.
## =====================================================================
ANALYTES <- c(rdw = "RDW (%)", wbc = "WBC (10^3/uL)",
              hct = "Hematocrit (%)", hgb = "Hemoglobin (g/dL)",
              mch = "MCH (pg)", mchc = "MCHC (g/dL)", rbc = "RBC (10^6/uL)",
              mcv = "MCV (fL)",
              lymph = "Lymphocyte (%)", neut = "Neutrophil (%)",
              mono = "Monocyte (%)", plt = "Platelets (10^3/uL)",
              albumin = "Albumin (g/dL)", alp = "ALP (U/L)",
              creat = "Creatinine (mg/dL)", glucose = "Glucose (mg/dL)")
## analytes absent from the frame are skipped, so an older M01 frame still runs
ANALYTES <- ANALYTES[names(ANALYTES) %in% names(ana)]
## the five CBC analytes NCHS named as shifting at the 2013-2014 boundary
NCHS_NAMED <- c("rdw", "wbc", "hct", "mch", "mchc")
rows <- list()
for (v in names(ANALYTES)) for (fr in c("ref_adol", "ref_young", NA)) {
  sc <- if (is.na(fr)) NULL else fr
  r <- step_est(ana, v, 2011, 2013, sc)
  rows[[length(rows) + 1]] <- data.frame(
    analyte = v, label = ANALYTES[[v]], boundary = "2011-12 to 2013-14",
    frame = ifelse(is.na(fr), "all", fr), n = r$n,
    mean_pre = r$m0, mean_post = r$m1, step = r$est, lo = r$lo, hi = r$hi,
    sd_units = r$est / r$sd, p = r$p, row.names = NULL)
}
## EXTERNAL VALIDATION, chemistry boundary. Run on alp_raw, the released
## values, because the harmonised series has had the NCHS equation applied and
## therefore cannot be compared against that equation's own prediction. Both
## directions are reported: the raw step should approach the predicted +5.84,
## and the harmonised step should approach zero.
for (v in c("alp_raw", "alp")) for (fr in c("ref_adol", "ref_young", NA)) {
  sc <- if (is.na(fr)) NULL else fr
  r <- step_est(ana, v, 2015, 2017, sc)
  rows[[length(rows) + 1]] <- data.frame(
    analyte = v, label = paste0("ALP (U/L), ",
      ifelse(v == "alp_raw", "as released", "after NCHS bridging")),
    boundary = "2015-16 to 2017-18", frame = ifelse(is.na(fr), "all", fr),
    n = r$n, mean_pre = r$m0, mean_post = r$m1, step = r$est,
    lo = r$lo, hi = r$hi, sd_units = r$est / r$sd, p = r$p, row.names = NULL)
}
## FLANKING-BOUNDARY DIAGNOSTIC (kept for m_analyte_steps.csv continuity; the
## full ten-boundary series is section 2b). A calibration step should appear at the
## analyser boundary and NOT at the boundaries either side of it. This is what
## separates a step from ordinary secular drift, and it is what the raw
## adjacent-cycle estimate cannot tell you on its own -- referent glucose rose
## +2.96 at 2011-13 with no chemistry platform change there.
for (v in names(ANALYTES)) for (bd in list(c(2009, 2011), c(2013, 2015))) {
  r <- step_est(ana, v, bd[1], bd[2], "ref_adol")
  rows[[length(rows) + 1]] <- data.frame(
    analyte = v, label = ANALYTES[[v]],
    boundary = paste0(bd[1], " to ", bd[2], " (flanking)"),
    frame = "ref_adol", n = r$n, mean_pre = r$m0, mean_post = r$m1,
    step = r$est, lo = r$lo, hi = r$hi, sd_units = r$est / r$sd,
    p = r$p, row.names = NULL)
}
asteps <- bind_rows(rows)
write.csv(asteps, file.path(results_dir, "m_analyte_steps.csv"), row.names = FALSE)

gs <- function(v, fr, b = "2011-12 to 2013-14")
  asteps[asteps$analyte == v & asteps$frame == fr & asteps$boundary == b, ][1, ]
for (v in c("rdw", "wbc", "mcv", "lymph")) {
  a <- gs(v, "ref_adol"); y <- gs(v, "ref_young")
  put(paste0("STEP_", toupper(ifelse(v == "lymph", "lym", v)), "_ADOL"), a$step)
  put(paste0("SD_",   toupper(ifelse(v == "lymph", "lym", v)), "_ADOL"), a$sd_units)
  put(paste0("STEP_", toupper(ifelse(v == "lymph", "lym", v)), "_YOUNG"), y$step)
}
cat("\n===== ANALYTE STEPS, ADOLESCENT REFERENT =====\n")
print(asteps %>% filter(frame == "ref_adol", boundary == "2011-12 to 2013-14") %>%
        select(label, n, step, lo, hi, sd_units, p), row.names = FALSE)
cat("\nRDW: compare against the +0.71 NCHS reports for the full weighted sample.\n")
cat("ALP at 2015-16 to 2017-18 is the external check; NCHS predicts +5.84.\n")
print(asteps %>% filter(boundary == "2015-16 to 2017-18") %>%
        select(label, frame, n, step, lo, hi, p), row.names = FALSE)

## =====================================================================
## 2b. ALL-BOUNDARY PLACEBO SERIES
##     Same adolescent-referent contrast, every adjacent-cycle pair. For RDW
##     and WBC the analyser account predicts a step at 2011-2013 and nowhere
##     else. Chemistry analytes are EXPECTED to step where the chemistry
##     platform changed without a bridging study:
##       1999-2001  Hitachi 917 (1999-2000) to Beckman LX20     (first run 2026-10-17
##                  found albumin -0.61 SD, ALP -0.43 SD, creatinine +0.37 SD here)
##       2005-2007  LX20 to DxC 800
##     Those are confirmations, not noise. 2017-2022 spans the pandemic gap and
##     the Cobas transitions, which are harmonised by NCHS equations upstream.
##     A hit anywhere else is either real two-year change in adolescents
##     (implausible) or an undocumented method change worth reporting.
##
##     v2026-10-18: each row now carries the referent MEAN in the earlier cycle
##     and the PERCENT step, so a proportional correction can be applied where
##     an analyte's form calls for it (ALP: NCHS's own ALP bridging equations
##     are log-linear, and adolescent ALP is roughly twice adult ALP, so an
##     additive adolescent step would overcorrect adults).
## =====================================================================
yrs   <- sort(unique(ana$year))
pairs <- data.frame(from = head(yrs, -1), to = tail(yrs, -1))
rows  <- list()
for (v in names(ANALYTES)) for (i in seq_len(nrow(pairs))) {
  r <- step_est(ana, v, pairs$from[i], pairs$to[i], "ref_adol")
  ## referent mean in the EARLIER cycle, for the percent step
  .sub <- ana[ana$ref_adol & ana$year == pairs$from[i] & !is.na(ana[[v]]) &
              !is.na(ana$wt) & ana$wt > 0, , drop = FALSE]
  mf <- if (nrow(.sub) > 30)
          tryCatch(as.numeric(coef(svymean(as.formula(paste0("~", v)), mk_design(.sub)))),
                   error = function(e) NA_real_) else NA_real_
  rows[[length(rows) + 1]] <- data.frame(
    analyte = v, label = ANALYTES[[v]], from = pairs$from[i], to = pairs$to[i],
    boundary = paste0(pairs$from[i], "-", pairs$to[i]), n = r$n,
    step = r$est, lo = r$lo, hi = r$hi, sd_units = r$est / r$sd, p = r$p,
    mean_from = mf, pct = 100 * r$est / mf,
    pct_lo = 100 * r$lo / mf, pct_hi = 100 * r$hi / mf,
    row.names = NULL)
}
allb <- bind_rows(rows)
write.csv(allb, file.path(results_dir, "m_all_boundary_steps.csv"), row.names = FALSE)

spec <- allb %>% filter(!is.na(p)) %>% group_by(analyte, label) %>%
  summarise(boundaries = n(), hits_p05 = sum(p < .05), hits_p01 = sum(p < .01),
            largest_sd = max(abs(sd_units)),
            at = boundary[which.max(abs(sd_units))], .groups = "drop") %>%
  arrange(desc(largest_sd))
write.csv(spec, file.path(results_dir, "m_all_boundary_specificity.csv"), row.names = FALSE)

cat("\n===== ALL-BOUNDARY PLACEBO SERIES, ADOLESCENT REFERENT (SD units) =====\n")
cat("* P<.05  ** P<.01  *** P<.001. RDW and WBC should light up at 2011-2013 only.\n\n")
wide <- allb %>%
  mutate(cell = sprintf("%+.2f%s", sd_units,
         ifelse(is.na(p), "", ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "")))))) %>%
  select(label, boundary, cell) %>%
  tidyr::pivot_wider(names_from = boundary, values_from = cell)
print(as.data.frame(wide), row.names = FALSE)
cat("\n--- specificity ---\n"); print(as.data.frame(spec), row.names = FALSE)
.sp <- function(v, col) { r <- spec[spec$analyte == v, ]; if (nrow(r)) r[[col]][1] else NA }
put("N_BOUNDARIES", nrow(pairs), "%d")
put("HITS_RDW", .sp("rdw", "hits_p05"), "%d"); put("HITS_WBC", .sp("wbc", "hits_p05"), "%d")
put("RDW_HIT_AT", as.character(.sp("rdw", "at"))); put("WBC_HIT_AT", as.character(.sp("wbc", "at")))
cat(sprintf("\neMethods sentence: of %d adjacent-cycle boundaries, RDW stepped at %s (%d of %d significant at P<.05) and WBC at %s (%d of %d).\n",
            nrow(pairs), .sp("rdw","at"), .sp("rdw","hits_p05"), nrow(pairs),
            .sp("wbc","at"), .sp("wbc","hits_p05"), nrow(pairs)))

## =====================================================================
## 2c. CORRECTION FORM: additive or proportional?
##     A mean step cannot tell an additive offset from a proportional one,
##     and they diverge in the tails - where RAR cutpoints and top-quartile
##     PhenoAge live. NCHS's chemistry equations are mostly multiplicative, so
##     the question is fair. Weighted quantiles per cycle; cycles are
##     independent samples, so SE(diff) = sqrt(se1^2 + se2^2).
## =====================================================================
PROBS <- c(.10, .25, .50, .75, .90)
## analytes that get a correction factor: NCHS-named, present in the frame, and
## with a significant adolescent step at 2011-2013
CORR_SET <- intersect(NCHS_NAMED, names(ANALYTES))
CORR_SET <- CORR_SET[vapply(CORR_SET, function(v) { g <- gs(v, "ref_adol"); !is.na(g$p) && g$p < .05 }, logical(1))]
cat("\ncorrection factors will be offered for:", paste(CORR_SET, collapse = ", "), "\n")
cat("NCHS-named analytes WITHOUT a significant adolescent step (no factor offered):",
    paste(setdiff(intersect(NCHS_NAMED, names(ANALYTES)), CORR_SET), collapse = ", "), "\n")
svyq <- function(df, v, probs) {
  df <- df[!is.na(df[[v]]) & !is.na(df$wt) & df$wt > 0, , drop = FALSE]
  q <- svyquantile(as.formula(paste0("~", v)), mk_design(df), quantiles = probs,
                   ci = TRUE, se = TRUE)
  ## survey >= 4.1 returns a list of matrices (quantile, ci.2.5, ci.97.5, se);
  ## older versions differ, so fall back to the CI width when se is absent
  m  <- if (is.list(q) && !is.data.frame(q)) q[[1]] else as.matrix(q)
  cn <- colnames(m)
  qq <- if ("quantile" %in% cn) m[, "quantile"] else m[, 1]
  se <- if ("se" %in% cn) m[, "se"] else
        (m[, grep("97.5", cn)[1]] - m[, grep("2.5", cn)[1]]) / 3.92
  data.frame(prob = probs, q = as.numeric(qq), se = as.numeric(se))
}
cform <- list()
for (v in CORR_SET) {
  a <- svyq(ana[ana$ref_adol & ana$year == 2011, ], v, PROBS)
  b <- svyq(ana[ana$ref_adol & ana$year == 2013, ], v, PROBS)
  r <- data.frame(analyte = v, label = ANALYTES[[v]], prob = PROBS,
                  q2011 = a$q, q2013 = b$q, step = b$q - a$q, se = sqrt(a$se^2 + b$se^2))
  r$lo <- r$step - 1.96 * r$se; r$hi <- r$step + 1.96 * r$se
  r$pct <- 100 * r$step / r$q2011; r$pct_se <- 100 * r$se / r$q2011
  cform[[v]] <- r
}
cform <- bind_rows(cform)
write.csv(cform, file.path(results_dir, "m_correction_form_quantiles.csv"), row.names = FALSE)

callrows <- list()
for (v in CORR_SET) {
  r  <- cform[cform$analyte == v, ]
  fa <- summary(lm(step ~ q2011, data = r, weights = 1 / r$se^2))$coefficients
  fp <- summary(lm(pct  ~ q2011, data = r, weights = 1 / r$pct_se^2))$coefficients
  pa <- fa["q2011", 4]; pp <- fp["q2011", 4]
  form <- if (pa >= .05 && pp <  .05) "additive" else
          if (pp >= .05 && pa <  .05) "proportional" else
          if (pa >= .05 && pp >= .05) "indistinguishable - report additive" else
                                      "neither constant - report by percentile"
  callrows[[v]] <- data.frame(
    analyte = v, p_absolute_slope = pa, p_percent_slope = pp, form = form,
    additive_factor = weighted.mean(r$step, 1 / r$se^2),
    proportional_factor = 1 + weighted.mean(r$pct, 1 / r$pct_se^2) / 100, row.names = NULL)
}
cformcall <- bind_rows(callrows)
write.csv(cformcall, file.path(results_dir, "m_correction_form_call.csv"), row.names = FALSE)

cat("\n===== CORRECTION FORM: STEP BY PERCENTILE, ADOLESCENT REFERENT =====\n")
for (v in CORR_SET) {
  r <- cform[cform$analyte == v, ]
  cat(sprintf("\n%s\n  %-5s %8s %8s %18s %8s\n", ANALYTES[[v]], "pct", "2011-12", "2013-14", "step (95% CI)", "% step"))
  for (i in seq_len(nrow(r)))
    cat(sprintf("  P%-4d %8.2f %8.2f %7.2f (%.2f,%.2f) %7.1f%%\n", 100 * r$prob[i],
                r$q2011[i], r$q2013[i], r$step[i], r$lo[i], r$hi[i], r$pct[i]))
}
cat("\n--- form call (flat absolute step = additive; flat percent step = proportional) ---\n")
print(cformcall, row.names = FALSE)
cat("Five percentiles give the slope tests little power - read the table too.\n")
for (v in CORR_SET) {
  r <- cformcall[cformcall$analyte == v, ]
  put(paste0("FORM_", toupper(v)), r$form); put(paste0("PROP_", toupper(v)), r$proportional_factor, "%.4f")
}

## =====================================================================
## 2d. THE CORRECTION FACTORS - the letter's deliverable
##     Apply to raw analyte values for participants examined 2013-2014 onward,
##     BEFORE computing any composite. Additive: subtract. Proportional:
##     divide. NCHS's reported whole-sample shift sits beside ours.
## =====================================================================
NCHS_REPORTED <- c(rdw = 13.55 - 12.84, wbc = 7.41 - 7.03)   # CBC_H analytic note gives means for these two only
corr <- bind_rows(lapply(CORR_SET, function(v) {
  a <- gs(v, "ref_adol"); f <- cformcall[cformcall$analyte == v, ]
  data.frame(analyte = v, label = ANALYTES[[v]],
             applies_to = "cycles 2013-2014 onward (Coulter DxH 800)",
             additive_step = a$step, lo = a$lo, hi = a$hi, sd_units = a$sd_units,
             form_call = f$form, proportional_factor = f$proportional_factor,
             boundaries_tested = .sp(v, "boundaries"), boundaries_with_step = .sp(v, "hits_p05"),
             nchs_reported_shift = if (v %in% names(NCHS_REPORTED)) NCHS_REPORTED[[v]] else NA_real_, referent = "adolescents 12-19 y, non-pregnant",
             row.names = NULL)
}))
write.csv(corr, file.path(results_dir, "m_correction_factors.csv"), row.names = FALSE)
cat("\n===== CORRECTION FACTORS (apply to raw analyte values, 2013-2014 onward) =====\n")
print(corr %>% select(label, additive_step, lo, hi, form_call, proportional_factor,
                      boundaries_with_step, boundaries_tested, nchs_reported_shift),
      row.names = FALSE)
cat("\nHow to apply: for each participant examined in 2013-2014 or later, replace\n")
cat("  RDW with RDW - additive_step  (or RDW / proportional_factor if the form call says so)\n")
cat("  WBC with WBC - additive_step  (likewise)\n")
cat("then compute PhenoAge, RAR, SIRI or any other composite from the corrected values.\n")
for (v in CORR_SET) {
  r <- corr[corr$analyte == v, ]
  put(paste0("CORR_", toupper(v)), r$additive_step)
  put(paste0("CORR_", toupper(v), "_LO"), r$lo); put(paste0("CORR_", toupper(v), "_HI"), r$hi)
}

## =====================================================================
## 3. Composite index steps, and the algebra table
## =====================================================================
INDEX <- c(RAR = "RAR (RDW / albumin)", NLR = "NLR (neut % / lymph %)",
           SIRI = "SIRI (neut x mono / lymph counts)",
           SII = "SII (platelet x neut / lymph counts)",
           PLR = "PLR (platelet / lymph count)")
RETAINS <- c(RAR = "RDW, undiluted", NLR = "none - WBC cancels",
             SIRI = "WBC and differential", SII = "platelets only - WBC cancels",
             PLR = "platelets, WBC in denominator")
rows <- list()
for (v in names(INDEX)) for (fr in c("ref_adol", "ref_young", NA)) {
  sc <- if (is.na(fr)) NULL else fr
  r <- step_est(ana, v, 2011, 2013, sc)
  rows[[length(rows) + 1]] <- data.frame(
    index = v, label = INDEX[[v]], retains = RETAINS[[v]],
    frame = ifelse(is.na(fr), "all", fr), n = r$n,
    mean_pre = r$m0, mean_post = r$m1, step = r$est,
    pct = 100 * r$est / r$m0, sd_units = r$est / r$sd,
    lo = r$lo, hi = r$hi, p = r$p, row.names = NULL)
}
isteps <- bind_rows(rows)
write.csv(isteps, file.path(results_dir, "m_index_steps.csv"), row.names = FALSE)
for (v in names(INDEX)) {
  r <- isteps[isteps$index == v & isteps$frame == "ref_adol", ][1, ]
  put(paste0("STEP_", v, "_PCT"), r$pct, "%.1f")
  put(paste0("STEP_", v, "_SD"),  r$sd_units)
  put(paste0("P_", v), fmt_p(r$p))
}
cat("\n===== INDEX STEPS, ADOLESCENT REFERENT =====\n")
cat("`retains` is the organising point: exposure follows which affected\n")
cat("analytes survive the index's algebra.\n")
print(isteps %>% filter(frame == "ref_adol") %>%
        select(label, retains, n, step, pct, sd_units, p), row.names = FALSE)

## =====================================================================
## 4. PhenoAge -- computed EXACTLY as 01_build_data computes it
##
## Phenotypic age is computed by TRAINING the Gompertz model in NHANES III with
## phenoage_calc(fit = NULL) and projecting onto NHANES 1999-2023 with that fit;
## the CRP-free variant is RETRAINED on the 8-marker set rather than obtained by
## deleting the CRP term from the published 9-marker equation. The two are
## different measures on different scales (deleting a term leaves the surviving
## coefficients at their 9-marker values, so that measure is not calibrated).
##
## This section now makes the same call, so both manuscripts report ONE
## measure. Three details are matched deliberately:
##   1. training sample   NHANES3, ages 20-84, exactly as in 01_build_data
##   2. biomarker set     the _nocrp 8, in the same order
##   3. 5-SD trim         within cycle, mirroring BioAge's nhanes_all.R
##
## The age ceiling also now matches 01_build_data's B1a: drop the top-coded
## category only (age < 85 through 2005-2006, age < 80 from 2007-2008), which
## retains 84 and 79 respectively. B1b, the common 20-79 ceiling, is already
## written into 01_build_data as a sensitivity and is re-run in section 7b.
## =====================================================================
AGE_MIN <- 20
age_topcode <- function(year) ifelse(year <= 2005, 85, 80)   # wave <= 4 -> 85
PA_INPUTS <- c("albumin", "creat", "glucose", "lymph", "mcv", "rdw", "alp", "wbc")
PHENO_BIOMARKERS <- c("albumin_gL", "lymph", "mcv", "glucose_mmol",
                      "rdw", "creat_umol", "alp", "wbc")

## 5-SD trim, WITHIN cycle, as BioAge's own cleaning script does. Applied ONCE
## to the observed data and never re-applied after a counterfactual shift: a
## trim recomputed on shifted values would interact with the correction and
## contaminate the result. Note the window is cycle-specific, so for RDW -- 
## whose mean jumps 0.72 at the boundary -- the post-boundary window slides up
## with it. That is faithful to the clinical pipeline, not an endorsement.
trim5 <- function(df, cols) {
  df %>% group_by(cycle) %>%
    mutate(across(all_of(cols), ~ ifelse(
      . > mean(., na.rm = TRUE) + 5 * sd(., na.rm = TRUE) |
      . < mean(., na.rm = TRUE) - 5 * sd(., na.rm = TRUE), NA_real_, .))) %>%
    ungroup()
}

adu <- ana %>%
  filter(age >= AGE_MIN, age < age_topcode(year)) %>%
  trim5(PA_INPUTS)

## Put the frame in the units and column names phenoage_calc expects. These
## conversions mirror 01_build_data: albumin g/dL -> g/L, creatinine
## mg/dL -> umol/L, glucose mg/dL -> mmol/L.
as_bioage <- function(df, over = list()) {
  g <- function(v) if (!is.null(over[[v]])) over[[v]] else df[[v]]
  data.frame(
    sampleID     = paste0(df$year, "_", df$seqn),
    age          = df$age,
    gender       = df$sex,
    albumin_gL   = g("albumin") * 10,
    creat_umol   = g("creat")   * 88.4017,
    glucose_mmol = g("glucose") * 0.0555,
    lymph        = g("lymph"),
    mcv          = g("mcv"),
    rdw          = g("rdw"),
    alp          = g("alp"),
    wbc          = g("wbc"))
}

## Train once, in NHANES III, ages 20-84 -- the same fit 01_build_data uses.
message("training PhenoAge in NHANES III (same call as 01_build_data) ...")
train_pheno <- phenoage_calc(NHANES3 %>% filter(age >= AGE_MIN, age <= 84),
                             PHENO_BIOMARKERS, fit = NULL)
stopifnot(!is.null(train_pheno$fit))

## Project. Returns $data with a `phenoage` column; join back by sampleID
## rather than assuming row order is preserved.
project_pheno <- function(df, over = list()) {
  bd  <- as_bioage(df, over)
  out <- phenoage_calc(bd, PHENO_BIOMARKERS, fit = train_pheno$fit)$data
  if (!"phenoage" %in% names(out))
    stop("phenoage_calc returned no `phenoage` column. Columns: ",
         paste(names(out), collapse = ", "))
  as.numeric(out$phenoage[match(bd$sampleID, out$sampleID)])
}

adu$pa <- project_pheno(adu)
adu <- adu %>% filter(!is.na(pa))
put("N_PHENO", format(nrow(adu), big.mark = ","))
message("PhenoAge computable in ", nrow(adu), " adults below the top-code")

## Acceleration = residual of PhenoAge on chronological age, fitted ONCE in
## the pooled adult frame, as 01_build_data does (residuals() of a pooled lm).
## Used for the DISPLAY series only; every step below is estimated on PhenoAge
## itself with age adjustment, which is algebraically the same contrast.
accel_of <- function(pa, df) {
  ok <- !is.na(pa) & !is.na(df$age)
  dd <- df[ok, ]; dd$.pa <- pa[ok]
  m <- svyglm(.pa ~ age, design = mk_design(dd)); cf <- coef(m)
  out <- rep(NA_real_, length(pa))
  out[ok] <- dd$.pa - (cf[[1]] + cf[[2]] * dd$age)
  out
}
adu$accel <- accel_of(adu$pa, adu)

## Scale fingerprint. If these three do not match the clinical pipeline's
## phenoage_nocrp, the two runs are still not on one measure and the step in
## years is not comparable between the papers.
cat("\n===== PHENOAGE SCALE FINGERPRINT =====\n")
.d11 <- adu %>% filter(year == 2011)
cat(sprintf("SD of PhenoAge, 2011-2012        : %.3f\n", sd(.d11$pa)))
cat(sprintf("slope of PhenoAge on age, pooled : %.3f\n",
            coef(svyglm(pa ~ age, design = mk_design(adu)))[2]))
cat(sprintf("mean PhenoAge - mean age, 2011-12: %.3f\n",
            weighted.mean(.d11$pa - .d11$age, .d11$wt)))
cat("Compare against phenoage_nocrp in the clinical data. They should match.\n")

## ---- baseline step, age-group steps, slope test ----
## PRIMARY ESTIMATOR: on PhenoAge itself, age-adjusted, restricted to the two
## boundary cycles. NOT on the pooled acceleration residual. A residual's level
## in any cycle depends on the composition of the pool it was fitted in, so the
## a residual-based step changes when the pool changes (cycles included, the
## age window, complete-case vs imputation). A
## pure intercept shift cancels in a between-cycle difference, but a change in
## the fitted slope rotates the residuals, and cycles with different mean ages
## then move by different amounts. Estimating on PhenoAge is invariant to all
## of it. The residual is kept for the display series only.
pr <- list()
b  <- step_est(adu, "pa",    2011, 2013, adjust_age = TRUE)
br <- step_est(adu, "accel", 2011, 2013, adjust_age = TRUE)
pr[[1]] <- data.frame(what = "overall, on PhenoAge (PRIMARY)", n = b$n, est = b$est,
                      lo = b$lo, hi = b$hi, p = b$p)
pr[[2]] <- data.frame(what = "overall, via pooled residual (v2 route)", n = br$n,
                      est = br$est, lo = br$lo, hi = br$hi, p = br$p)
adu$agegrp <- cut(adu$age, c(19, 39, 59, Inf), labels = c("20-39", "40-59", "60+"))
for (g in levels(adu$agegrp)) {
  adu$.in <- adu$agegrp == g
  r <- step_est(adu, "pa", 2011, 2013, ".in", adjust_age = TRUE)
  pr[[length(pr) + 1]] <- data.frame(what = paste("age", g), n = r$n, est = r$est,
                                     lo = r$lo, hi = r$hi, p = r$p)
}
## slope of PhenoAge on chronological age, across cycles: an additive offset
## steps the intercept and leaves this interaction null
des_adu <- mk_design(adu)
slope_m <- svyglm(pa ~ age * factor(year), design = des_adu)
slope_p <- regTermTest(slope_m, "age:factor(year)")$p
pr[[length(pr) + 1]] <- data.frame(what = "age x cycle interaction (slope test)",
                                   n = nrow(adu), est = NA, lo = NA, hi = NA,
                                   p = as.numeric(slope_p))
psteps <- bind_rows(pr)
write.csv(psteps, file.path(results_dir, "m_pheno_steps.csv"), row.names = FALSE)
cat("\n===== PHENOTYPIC AGE ACCELERATION STEP =====\n")
print(psteps, row.names = FALSE)
cat("\nPrior run gave 2.61 (2.12-3.10) overall; 2.83 / 2.51 / 2.53 by age group;\n")
cat("slope P = .34. Equal steps across age groups plus a null slope test is the\n")
cat("signature of an additive offset: level breaks, within-cycle ranking does not.\n")
put("PA_STEP", b$est); put("PA_LO", b$lo); put("PA_HI", b$hi)
put("PA_SLOPE_P", fmt_p(as.numeric(slope_p)))

## =====================================================================
## 5. Counterfactual and leave-one-out
##    Shift each analyte back by its ADOLESCENT-referent step for everyone
##    examined from 2013-2014 on, leave all else untouched, recompute.
## =====================================================================
ADOL_STEP <- setNames(
  sapply(PA_INPUTS, function(v) { r <- gs(v, "ref_adol"); ifelse(is.na(r$step), 0, r$step) }),
  PA_INPUTS)
cat("\nadolescent-referent steps applied in the counterfactual:\n"); print(round(ADOL_STEP, 4))

shift_back <- function(df, vars) {
  ov <- list(); post <- df$year >= 2013
  for (v in vars) { x <- df[[v]]; x[post] <- x[post] - ADOL_STEP[[v]]; ov[[v]] <- x }
  ov
}
## The counterfactual is now measured on PhenoAge with the same primary
## estimator, so it no longer depends on the pooled residual fit at all.
run_cf <- function(vars, tag, kind) {
  dd <- adu; dd$pa_cf <- project_pheno(adu, shift_back(adu, vars))
  r <- step_est(dd, "pa_cf", 2011, 2013, adjust_age = TRUE)
  data.frame(correction = tag, kind = kind, vars = paste(vars, collapse = "+"),
             residual = r$est, lo = r$lo, hi = r$hi, p = r$p,
             pct_removed = 100 * (1 - r$est / b$est), row.names = NULL)
}
## HAEMATOLOGY vs CHEMISTRY. There is no chemistry platform change at the
## 2011-12 / 2013-14 boundary, so correcting a chemistry analyte there is
## unjustified in principle: its adjacent-cycle difference is secular change or
## noise, not calibration. Referent glucose moved +2.96 at this boundary and
## can "remove" a share of the step comparable to WBC -- which is not evidence
## about glucose, it is the method's own noise floor. Chemistry corrections are
## therefore reported as PLACEBO tests, and the largest of them bounds the
## precision of the whole counterfactual.
HAEM <- c("rdw", "wbc", "mcv", "lymph")
CHEM <- setdiff(PA_INPUTS, HAEM)
cf <- bind_rows(
  data.frame(correction = "none", kind = "baseline", vars = "", residual = b$est,
             lo = b$lo, hi = b$hi, p = b$p, pct_removed = 0),
  run_cf("rdw", "RDW", "primary"), run_cf("wbc", "WBC", "primary"),
  run_cf(c("rdw", "wbc"), "RDW + WBC", "primary"),
  run_cf(HAEM, "all four CBC inputs", "primary"),
  bind_rows(lapply(HAEM, function(v) run_cf(v, paste0("LOO: ", v), "haematology"))),
  bind_rows(lapply(CHEM, function(v) run_cf(v, paste0("PLACEBO: ", v), "placebo"))))
write.csv(cf, file.path(results_dir, "m_counterfactual.csv"), row.names = FALSE)
cat("\n===== COUNTERFACTUAL =====\n"); print(cf, row.names = FALSE)
cat("\nPrior run: RDW 77.7%, WBC 21.4%, both 99.1% (residual 0.02, P = .92),\n")
cat("all four 96.5%; every other input 8.0% or less. ALP has the largest\n")
cat("chemistry shift yet removes ~2% because its coefficient is 0.00188 --\n")
cat("displacement is shift x weight, not shift alone.\n")
gc_ <- function(tag) cf[cf$correction == tag, ][1, ]
put("CF_RDW_PCT",  gc_("RDW")$pct_removed, "%.1f")
put("CF_WBC_PCT",  gc_("WBC")$pct_removed, "%.1f")
put("CF_BOTH_PCT", gc_("RDW + WBC")$pct_removed, "%.1f")
put("CF_BOTH_RESID", gc_("RDW + WBC")$residual)
put("CF_BOTH_P", fmt_p(gc_("RDW + WBC")$p))
put("CF_MAX_OTHER", max(cf$pct_removed[cf$kind == "haematology" &
                        !grepl("rdw|wbc", cf$correction)]), "%.1f")
put("CF_PLACEBO_MAX", max(abs(cf$pct_removed[cf$kind == "placebo"])), "%.1f")
cat("\nPLACEBO (chemistry, no platform change at this boundary) -- the largest\n")
cat("of these is the counterfactual's noise floor, not a finding:\n")
print(cf[cf$kind == "placebo", c("correction", "residual", "pct_removed", "p")],
      row.names = FALSE)

## =====================================================================
## 5b. DERIVED MEASURES AT EVERY BOUNDARY, OBSERVED AND CORRECTED (2026-10-18)
##     For the letter's figure row D and table row group: the same adjacent-
##     cycle comparison as the analyte rows, but in ADULTS, where these
##     measures are used, before and after the RDW/WBC correction.
##       phenotypic age: years, age-adjusted, estimated on PhenoAge itself
##                       (same estimator as section 4)
##       RAR, NLR:       % of the earlier cycle's mean
##     RAR and NLR are recomputed here from the adult frame's own (trimmed)
##     values so observed and corrected use identical inputs. NLR contains
##     neither RDW nor WBC, so its corrected series equals its observed one -
##     that is the point of including it.
## =====================================================================
adu$pa_corr  <- project_pheno(adu, shift_back(adu, c("rdw", "wbc")))
adu$RAR_obs  <- adu$rdw / adu$albumin
adu$RAR_corr <- ifelse(adu$year >= 2013, adu$rdw - ADOL_STEP[["rdw"]], adu$rdw) / adu$albumin
adu$NLR_obs  <- ifelse(adu$lymph > 0, adu$neut / adu$lymph, NA_real_)
adu$NLR_corr <- adu$NLR_obs
## PLR = platelets / lymphocyte count, lymphocyte count = WBC x lymph% / 100.
## Correcting WBC changes each person's denominator, so the effect on the
## boundary step is empirical, not algebraic - same reason as RAR.
.wbc_corr    <- ifelse(adu$year >= 2013, adu$wbc - ADOL_STEP[["wbc"]], adu$wbc)
adu$PLR_obs  <- ifelse(adu$wbc * adu$lymph > 0, adu$plt / (adu$wbc * adu$lymph / 100), NA_real_)
adu$PLR_corr <- ifelse(.wbc_corr * adu$lymph > 0, adu$plt / (.wbc_corr * adu$lymph / 100), NA_real_)

adu$SII_obs  <- ifelse(adu$lymph > 0, adu$plt * adu$neut / adu$lymph, NA_real_)
adu$SII_corr <- adu$SII_obs    # WBC cancels; neither RDW nor WBC enters

DER <- list(pa  = c(obs = "pa",      corr = "pa_corr",  unit = "years"),
            RAR = c(obs = "RAR_obs", corr = "RAR_corr", unit = "%"),
            PLR = c(obs = "PLR_obs", corr = "PLR_corr", unit = "%"),
            SII = c(obs = "SII_obs", corr = "SII_corr", unit = "%"),
            NLR = c(obs = "NLR_obs", corr = "NLR_corr", unit = "%"))
yrs_a   <- sort(unique(adu$year))
pairs_a <- data.frame(from = head(yrs_a, -1), to = tail(yrs_a, -1))
drows <- list()
for (k in names(DER)) for (ver in c("obs", "corr")) for (i in seq_len(nrow(pairs_a))) {
  r <- step_est(adu, DER[[k]][[ver]], pairs_a$from[i], pairs_a$to[i], adjust_age = (k == "pa"))
  sc <- if (k == "pa") 1 else 100 / r$m0
  drows[[length(drows) + 1]] <- data.frame(
    measure = k, version = ifelse(ver == "obs", "observed", "corrected"),
    from = pairs_a$from[i], to = pairs_a$to[i],
    boundary = paste0(pairs_a$from[i], "-", pairs_a$to[i]), unit = DER[[k]][["unit"]],
    est = r$est * sc, lo = r$lo * sc, hi = r$hi * sc, p = r$p, n = r$n, row.names = NULL)
}
der <- bind_rows(drows)
write.csv(der, file.path(results_dir, "m_derived_boundary_steps.csv"), row.names = FALSE)

cat("\n===== DERIVED MEASURES, ADULTS: CHANGE AT EACH BOUNDARY, OBSERVED -> CORRECTED =====\n")
for (k in names(DER)) {
  cat(sprintf("\n%s (%s)\n", c(pa = "Phenotypic age", RAR = "RAR", PLR = "PLR", SII = "SII", NLR = "NLR")[[k]], DER[[k]][["unit"]]))
  o <- der[der$measure == k & der$version == "observed", ]
  c2 <- der[der$measure == k & der$version == "corrected", ]
  for (i in seq_len(nrow(o)))
    cat(sprintf("  %-10s %+7.2f -> %+7.2f%s\n", o$boundary[i], o$est[i], c2$est[i],
                if (o$boundary[i] == "2011-2013") "   <- analyzer boundary" else ""))
}
.g <- function(k, v) der$est[der$measure == k & der$version == v & der$boundary == "2011-2013"]
.gp <- function(k, v) der$p[der$measure == k & der$version == v & der$boundary == "2011-2013"]
put("DER_PA_OBS",  .g("pa", "observed"));  put("DER_PA_COR",  .g("pa", "corrected"))
put("DER_RAR_OBS", .g("RAR", "observed"), "%.1f"); put("DER_RAR_COR", .g("RAR", "corrected"), "%.1f")
put("DER_RAR_COR_P", fmt_p(.gp("RAR", "corrected")))
put("DER_NLR_OBS", .g("NLR", "observed"), "%.1f")
put("DER_PLR_OBS", .g("PLR", "observed"), "%.1f"); put("DER_PLR_COR", .g("PLR", "corrected"), "%.1f")
put("DER_PLR_COR_P", fmt_p(.gp("PLR", "corrected")))
put("DER_SII_OBS", .g("SII", "observed"), "%.1f"); put("DER_SII_OBS_P", fmt_p(.gp("SII", "observed")))
put("DER_NLR_OBS_P", fmt_p(.gp("NLR", "observed")))

## =====================================================================
## 6. Weighted quartile misclassification from pooling
##    The quantity that bears on the published literature, nearly all of which
##    cuts quartiles in the pooled sample. Weighted quantiles are computed
##    directly rather than via svyquantile, whose return shape varies by
##    survey package version.
## =====================================================================
wq <- function(x, w, probs) {
  ok <- !is.na(x) & !is.na(w) & w > 0; x <- x[ok]; w <- w[ok]
  o <- order(x); x <- x[o]; w <- w[o]; cw <- cumsum(w) / sum(w)
  sapply(probs, function(p) x[which(cw >= p)[1]])
}
qcut <- function(x, cuts) cut(x, c(-Inf, cuts, Inf), labels = 1:4, include.lowest = TRUE)
QMEAS <- c(RAR = "RAR", SIRI = "SIRI", SII = "SII", NLR = "NLR", PLR = "PLR")
qrows <- list()
qframe <- function(df, v) {
  x <- df[[v]]; w <- df$wt
  pooled <- qcut(x, wq(x, w, c(.25, .5, .75)))
  within <- rep(NA_integer_, length(x))
  for (y in unique(df$year)) {
    i <- df$year == y
    within[i] <- as.integer(qcut(x[i], wq(x[i], w[i], c(.25, .5, .75))))
  }
  ok <- !is.na(pooled) & !is.na(within) & !is.na(w) & w > 0
  pre <- ok & df$year < 2013; post <- ok & df$year >= 2013
  data.frame(
    measure = v,
    pct_moving = 100 * sum(w[ok] * (as.integer(pooled[ok]) != within[ok])) / sum(w[ok]),
    q4_pre  = 100 * sum(w[pre]  * (pooled[pre]  == "4")) / sum(w[pre]),
    q4_post = 100 * sum(w[post] * (pooled[post] == "4")) / sum(w[post]),
    row.names = NULL)
}
for (v in names(QMEAS)) if (sum(!is.na(ana[[v]])) > 5000)
  qrows[[length(qrows) + 1]] <- cbind(qframe(ana, v), population = "all ages")
## 2026-10-19: the letter's table is in adults throughout, so the ratios are
## also assigned quartiles within the adult phenotypic-age frame
for (v in names(QMEAS)) if (v %in% names(adu) && sum(!is.na(adu[[v]])) > 5000)
  qrows[[length(qrows) + 1]] <- cbind(qframe(adu, v), population = "adults")
qrows[[length(qrows) + 1]] <- cbind(qframe(adu %>% rename(PA = accel), "PA"), population = "adults")
qq <- bind_rows(qrows) %>% arrange(desc(pct_moving))
write.csv(qq, file.path(results_dir, "m_quartile_shift.csv"), row.names = FALSE)
cat("\n===== WEIGHTED QUARTILE MISCLASSIFICATION FROM POOLING =====\n")
print(qq, row.names = FALSE)
cat("\n25% in both of the last two columns would mean no distortion.\n")
## adults are the letter's population; all-ages tokens kept with a suffix for continuity
for (i in seq_len(nrow(qq))) put(paste0("Q_MOVE_", qq$measure[i], if (qq$population[i] == "all ages") "_ALLAGES" else ""),
                                qq$pct_moving[i], "%.1f")
put("Q4_PRE_RAR",  qq$q4_pre[qq$measure  == "RAR" & qq$population == "adults"], "%.1f")
put("Q4_POST_RAR", qq$q4_post[qq$measure == "RAR" & qq$population == "adults"], "%.1f")

## =====================================================================
## 7. Figure panel data
## =====================================================================
pan <- list()
for (v in c("rdw", "wbc")) for (fr in c("all", "ref_adol")) {
  dd <- if (fr == "all") ana else ana[ana$ref_adol, ]
  dd <- dd[!is.na(dd[[v]]), ]; dd$.y <- dd[[v]]
  m <- svyby(~.y, ~year, mk_design(dd), svymean, na.rm = TRUE, vartype = "ci")
  names(m) <- c("year", "mean", "lo", "hi")
  pan[[length(pan) + 1]] <- data.frame(panel = ifelse(v == "rdw", "A", "B"),
                                       series = v, frame = fr, m, row.names = NULL)
}
## Panel C. Each series' acceleration is its own pooled residual (series
## "none", "RDW", "RDW + WBC"), which re-centres on the pooled mean: removing
## the post-2013 inflation lowers the pooled mean and raises every earlier
## cycle by a constant. For the figure, the corrected series is therefore also
## expressed against the OBSERVED fit (series suffixed " (same reference)"),
## so the two coincide before 2013-2014 and differ only from the analyzer change.
.d0 <- adu[!is.na(adu$pa) & !is.na(adu$age), ]; .d0$.pa <- .d0$pa
.cf0 <- coef(svyglm(.pa ~ age, design = mk_design(.d0)))
for (tag in c("none", "RDW", "RDW + WBC")) {
  vars <- switch(tag, none = character(0), RDW = "rdw", `RDW + WBC` = c("rdw", "wbc"))
  pa <- if (!length(vars)) adu$pa else project_pheno(adu, shift_back(adu, vars))
  for (ref in c("own", "same")) {
    if (ref == "same" && !length(vars)) next
    dd <- adu
    dd$.y <- if (ref == "own") accel_of(pa, dd) else pa - (.cf0[[1]] + .cf0[[2]] * dd$age)
    dd <- dd[!is.na(dd$.y), ]
    m <- svyby(~.y, ~year, mk_design(dd), svymean, na.rm = TRUE, vartype = "ci")
    names(m) <- c("year", "mean", "lo", "hi")
    pan[[length(pan) + 1]] <- data.frame(panel = "C", series = if (ref == "own") tag else paste(tag, "(same reference)"),
                                         frame = "adults", m, row.names = NULL)
  }
}
pan[[length(pan) + 1]] <- data.frame(
  panel = "D",
  series = sub("^(LOO|PLACEBO): ", "", cf$correction[cf$kind %in% c("haematology", "placebo")]),
  frame = cf$kind[cf$kind %in% c("haematology", "placebo")], year = NA,
  mean = cf$pct_removed[cf$kind %in% c("haematology", "placebo")], lo = NA, hi = NA)
write.csv(bind_rows(pan), file.path(results_dir, "m_figure_panels.csv"), row.names = FALSE)

## =====================================================================
## 7b. Age-window sensitivity: re-run the primary step under the
##     cycle-specific top-code, so the effect of the window is a number.
## =====================================================================
## B1b in 01_build_data: a COMMON 20-79 ceiling instead of dropping only the
## top-coded category. Same measure, same fit, only the window differs.
alt <- ana %>% filter(age >= 20, age <= 79) %>% trim5(PA_INPUTS)
alt$pa <- project_pheno(alt); alt <- alt %>% filter(!is.na(pa))
ba <- step_est(alt, "pa", 2011, 2013, adjust_age = TRUE)
cat("\n===== AGE-WINDOW SENSITIVITY =====\n")
cat(sprintf("B1a, top-code dropped only (primary): %.2f (%.2f to %.2f), n = %d\n",
            b$est, b$lo, b$hi, b$n))
cat(sprintf("B1b, common 20-79 ceiling           : %.2f (%.2f to %.2f), n = %d\n",
            ba$est, ba$lo, ba$hi, ba$n))
cat("Both boundary cycles sit above 2007, so they share a ceiling either way;\n")
cat("a large difference here would mean the step is sensitive to the rest of\n")
cat("the pool, which the PhenoAge-based estimator is designed to prevent.\n")
put("PA_STEP_ALTWINDOW", ba$est)

## =====================================================================
## 8. Emit the placeholder table
## =====================================================================
ph <- tibble(token = names(PH), value = unlist(PH)) %>% arrange(token)
write.csv(ph, file.path(results_dir, "m_placeholders.csv"), row.names = FALSE)
cat("\n===== PLACEHOLDERS FOR THE BRIEF REPORT =====\n")
print(as.data.frame(ph), row.names = FALSE)
cat("\nFill every token in the Research Letter draft, then grep for the\n")
cat("opening bracket character to confirm none remain.\n")

cat("\nSaved to results/: m_analyte_steps.csv, m_all_boundary_steps.csv, m_derived_boundary_steps.csv,\n")
cat("  m_all_boundary_specificity.csv, m_correction_form_quantiles.csv,\n")
cat("  m_correction_form_call.csv, m_correction_factors.csv, m_index_steps.csv,\n")
cat("  m_pheno_steps.csv, m_counterfactual.csv, m_quartile_shift.csv,\n")
cat("  m_figure_panels.csv, m_placeholders.csv\n")

if (!is.null(.nosleep_handle) && !inherits(.nosleep_handle, "try-error"))
  try(NoSleepR::nosleep_off(.nosleep_handle), silent = TRUE)
message("Done: measurement analyses complete.")
