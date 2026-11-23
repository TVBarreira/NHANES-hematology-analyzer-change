###############################################################################
## M03_verification.R
##
## Run after M01_build_cbc_frame.R and M02_measurement_analyses.R. It reads
## data/m_cbc_frame.rds and results/m_correction_factors.csv, changes neither,
## and writes results/m_verification.csv. About 5-10 minutes; it downloads one
## table (CBC_L) for the 2021-2023 phlebotomy weight.
##
## Sections (labels as they appear in the "section" column of the output)
##  A. Age-group interaction. Tests whether the 2011-2012 to 2013-2014 change in
##     phenotypic age differs across age groups (20-39, 40-59, 60+), with an
##     age-group-by-cycle term in adults of those two cycles. This is the test
##     the manuscript reports (p = 0.47). For the record the file also gives
##     the age (continuous) x cycle term across all 11 cycles that M02 writes,
##     which is a different test.
##  B. Variance. Repeats the headline estimates with the survey design built on
##     the full sample and the subgroup analysed as a domain, alongside M02's
##     design built on the filtered rows. Point estimates and CIs agree.
##  C. 2021-2023 weight. NCHS recommends the phlebotomy weight WTPH2YR for
##     blood analytes in August 2021-August 2023. Repeats the 2017-2018 to
##     2021-2023 changes, the last point of Fig 1B and the quartile results
##     with that weight.
##  D. Quartiles after correction. Share of adults whose quartile differs when
##     cut in the pooled sample rather than within cycle, and share in the top
##     pooled quartile before and from 2013, as measured and after the RDW and
##     WBC correction.
###############################################################################

rm(list = ls())
.nosleep_handle <- NULL
if (requireNamespace("NoSleepR", quietly = TRUE))
  .nosleep_handle <- try(NoSleepR::nosleep_on(keep_display = TRUE), silent = TRUE)

library(dplyr); library(survey); library(BioAge); library(nhanesA)
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
options(survey.lonely.psu = "adjust")

## Project folder: set the environment variable NHANES_BIOAGE_DIR, or run from
## the project folder (it must contain data/ and results/ subfolders).
base_dir <- Sys.getenv("NHANES_BIOAGE_DIR", unset = getwd())
if (!dir.exists(base_dir)) stop("Project folder not found: ", base_dir)
data_dir <- file.path(base_dir, "data"); results_dir <- file.path(base_dir, "results")
if (!file.exists(file.path(data_dir, "m_cbc_frame.rds"))) stop("Run M01_build_cbc_frame.R and M02_measurement_analyses.R first")
ana <- readRDS(file.path(data_dir, "m_cbc_frame.rds"))
message("frame loaded: ", nrow(ana), " participants")

OUT <- list()
rec <- function(section, item, version, est, lo = NA, hi = NA, p = NA, n = NA)
  OUT[[length(OUT) + 1]] <<- data.frame(section = section, item = item, version = version,
                                        est = est, lo = lo, hi = hi, p = p, n = n)

## ---------------------------------------------------------------------------
## Rebuild the adult phenotypic-age frame exactly as M02 section 4 does
## ---------------------------------------------------------------------------
age_topcode <- function(year) ifelse(year <= 2005, 85, 80)
PA_INPUTS <- c("albumin", "creat", "glucose", "lymph", "mcv", "rdw", "alp", "wbc")
PHENO_BIOMARKERS <- c("albumin_gL", "lymph", "mcv", "glucose_mmol", "rdw", "creat_umol", "alp", "wbc")
trim5 <- function(df, cols) {
  df %>% group_by(cycle) %>%
    mutate(across(all_of(cols), ~ ifelse(
      . > mean(., na.rm = TRUE) + 5 * sd(., na.rm = TRUE) |
      . < mean(., na.rm = TRUE) - 5 * sd(., na.rm = TRUE), NA_real_, .))) %>%
    ungroup()
}
as_bioage <- function(df, rdw = df$rdw, wbc = df$wbc) data.frame(
  sampleID = paste0(df$year, "_", df$seqn), age = df$age, gender = df$sex,
  albumin_gL = df$albumin * 10, creat_umol = df$creat * 88.4017,
  glucose_mmol = df$glucose * 0.0555, lymph = df$lymph, mcv = df$mcv,
  rdw = rdw, alp = df$alp, wbc = wbc)
train_pheno <- phenoage_calc(NHANES3 %>% filter(age >= 20, age <= 84), PHENO_BIOMARKERS, fit = NULL)
project <- function(df, ...) {
  bd <- as_bioage(df, ...)
  out <- phenoage_calc(bd, PHENO_BIOMARKERS, fit = train_pheno$fit)$data
  as.numeric(out$phenoage[match(bd$sampleID, out$sampleID)])
}
adu <- ana %>% filter(age >= 20, age < age_topcode(year)) %>% trim5(PA_INPUTS)
adu$pa <- project(adu); adu <- adu %>% filter(!is.na(pa))
.cf <- read.csv(file.path(results_dir, "m_correction_factors.csv"))
STEP_RDW <- .cf$additive_step[.cf$analyte == "rdw"]; STEP_WBC <- .cf$additive_step[.cf$analyte == "wbc"]
stopifnot(length(STEP_RDW) == 1, length(STEP_WBC) == 1)
cat(sprintf("adolescent steps read from M02 output: RDW %.4f, WBC %.4f\n", STEP_RDW, STEP_WBC))
post13 <- adu$year >= 2013
adu$pa_corr  <- project(adu, rdw = ifelse(post13, adu$rdw - STEP_RDW, adu$rdw),
                             wbc = ifelse(post13, adu$wbc - STEP_WBC, adu$wbc))
adu$RAR_obs  <- adu$rdw / adu$albumin
adu$RAR_corr <- ifelse(post13, adu$rdw - STEP_RDW, adu$rdw) / adu$albumin
adu$agegrp <- cut(adu$age, c(19, 39, 59, Inf), labels = c("20-39", "40-59", "60+"))
message("adults with phenotypic age: ", nrow(adu), " (M02 reported 48,961)")
if (abs(nrow(adu) - 48961) > 5) warning("adult frame differs from the 48,961 reported by M02: ", nrow(adu))

## put the adult variables on the full frame so adults can be analysed as a domain
ana <- ana %>% left_join(adu %>% select(seqn, pa, pa_corr, RAR_obs, RAR_corr, agegrp), by = "seqn")
stopifnot(!any(grepl("\\.(x|y)$", names(ana))))

## ---------------------------------------------------------------------------
## Step estimators: M02's (design on the filtered rows) and the domain version
## ---------------------------------------------------------------------------
step_filtered <- function(df, yvar, y0, y1, cond, adjust_age = FALSE, wtvar = "wt") {
  cond <- !is.na(cond) & cond
  s <- df[df$year %in% c(y0, y1) & cond & !is.na(df[[yvar]]) & !is.na(df[[wtvar]]) & df[[wtvar]] > 0, ]
  s$.y <- s[[yvar]]; s$.post <- as.integer(s$year == y1); s$.w <- s[[wtvar]]
  d <- svydesign(id = ~psu, strata = ~stratum, weights = ~.w, nest = TRUE, data = s)
  m <- svyglm(if (adjust_age) .y ~ .post + age else .y ~ .post, design = d)
  co <- summary(m)$coefficients; tc <- qt(0.975, degf(d))
  c(est = co[".post", 1], lo = co[".post", 1] - tc * co[".post", 2],
    hi = co[".post", 1] + tc * co[".post", 2], p = co[".post", 4], n = nrow(s))
}
step_domain <- function(df, yvar, y0, y1, cond, adjust_age = FALSE, wtvar = "wt") {
  cond <- !is.na(cond) & cond
  s <- df[df$year %in% c(y0, y1) & !is.na(df[[wtvar]]) & df[[wtvar]] > 0, ]      # everyone in the two cycles
  s$.in <- cond[df$year %in% c(y0, y1) & !is.na(df[[wtvar]]) & df[[wtvar]] > 0] & !is.na(s[[yvar]])
  s$.y <- s[[yvar]]; s$.post <- as.integer(s$year == y1); s$.w <- s[[wtvar]]
  d  <- svydesign(id = ~psu, strata = ~stratum, weights = ~.w, nest = TRUE, data = s)
  ds <- subset(d, .in)
  m <- svyglm(if (adjust_age) .y ~ .post + age else .y ~ .post, design = ds)
  co <- summary(m)$coefficients; tc <- qt(0.975, degf(d))
  c(est = co[".post", 1], lo = co[".post", 1] - tc * co[".post", 2],
    hi = co[".post", 1] + tc * co[".post", 2], p = co[".post", 4], n = sum(s$.in))
}
both <- function(section, item, ...) {
  a <- step_filtered(...); b <- step_domain(...)
  rec(section, item, "as in M02 (design on filtered rows)", a[["est"]], a[["lo"]], a[["hi"]], a[["p"]], a[["n"]])
  rec(section, item, "domain analysis (full design, subset)", b[["est"]], b[["lo"]], b[["hi"]], b[["p"]], b[["n"]])
}

## ---------------------------------------------------------------------------
## B. Variance: the headline estimates both ways
## ---------------------------------------------------------------------------
adol <- ana$ref_adol; isadult <- !is.na(ana$pa)
both("B variance", "RDW, adolescents, 2011-2013", ana, "rdw", 2011, 2013, adol)
both("B variance", "WBC, adolescents, 2011-2013", ana, "wbc", 2011, 2013, adol)
both("B variance", "Phenotypic age, adults, 2011-2013", ana, "pa", 2011, 2013, isadult, TRUE)
both("B variance", "Phenotypic age after correction", ana, "pa_corr", 2011, 2013, isadult, TRUE)
for (v in c("rbc", "hct", "mch", "mchc", "mono", "hgb", "mcv", "lymph", "neut", "plt"))
  both("B variance", paste0(v, ", adolescents, 2011-2013"), ana, v, 2011, 2013, adol)

## ---------------------------------------------------------------------------
## A. The interaction test the Methods describe
## ---------------------------------------------------------------------------
s <- ana[ana$year %in% c(2011, 2013) & !is.na(ana$wt) & ana$wt > 0, ]
s$.post <- as.integer(s$year == 2013)
d2 <- subset(svydesign(id = ~psu, strata = ~stratum, weights = ~wt, nest = TRUE, data = s), !is.na(pa))
m_int <- svyglm(pa ~ .post * agegrp + age, design = d2)
p_int <- as.numeric(regTermTest(m_int, ~ .post:agegrp)$p)
rec("A interaction", "age-group x cycle, 2011-2012 vs 2013-2014, adults (test reported in the manuscript)",
    "domain analysis", NA, NA, NA, p_int, sum(!is.na(s$pa)))
## and what M02 actually computed, for the record
d_all <- svydesign(id = ~psu, strata = ~stratum, weights = ~wt, nest = TRUE, data = adu)
p_old <- as.numeric(regTermTest(svyglm(pa ~ age * factor(year), design = d_all), "age:factor(year)")$p)
rec("A interaction", "age (continuous) x cycle, all 11 cycles (as written by M02; a different test)",
    "as in M02", NA, NA, NA, p_old, nrow(adu))

## ---------------------------------------------------------------------------
## C. 2021-2023 phlebotomy weight
## ---------------------------------------------------------------------------
ph <- try(nhanes("CBC_L", translated = FALSE), silent = TRUE)
if (inherits(ph, "try-error") || is.null(ph) || !"WTPH2YR" %in% names(ph))
  ph <- try(nhanes("BIOPRO_L", translated = FALSE), silent = TRUE)
if (!inherits(ph, "try-error") && !is.null(ph) && "WTPH2YR" %in% names(ph)) {
  ana <- ana %>% left_join(data.frame(seqn = ph$SEQN, wtph2yr = as.numeric(ph$WTPH2YR)), by = "seqn")
  ncyc <- length(unique(ana$cycle))
  ana$wt_ph <- ifelse(ana$year == 2022, ana$wtph2yr / ncyc, ana$wt)
  cat("\n2021-2023 rows in the frame:", sum(ana$year == 2022),
      "| with a positive phlebotomy weight:", sum(ana$year == 2022 & !is.na(ana$wt_ph) & ana$wt_ph > 0), "\n")
  for (v in c("rdw", "wbc", "rbc", "hct", "mch", "mchc", "mono", "hgb", "mcv", "lymph", "neut", "plt")) {
    a <- step_filtered(ana, v, 2017, 2022, adol, FALSE, "wt")
    b <- step_filtered(ana, v, 2017, 2022, adol, FALSE, "wt_ph")
    rec("C phlebotomy weight", paste0(v, ", adolescents, 2017-2018 to 2021-2023"), "MEC weight (as in M02)", a[["est"]], a[["lo"]], a[["hi"]], a[["p"]], a[["n"]])
    rec("C phlebotomy weight", paste0(v, ", adolescents, 2017-2018 to 2021-2023"), "phlebotomy weight in 2021-2023", b[["est"]], b[["lo"]], b[["hi"]], b[["p"]], b[["n"]])
  }
  ## last point of Fig 1B: mean acceleration in 2021-2023 under each weight
  for (wv in c("wt", "wt_ph")) {
    dd <- ana[!is.na(ana$pa) & !is.na(ana[[wv]]) & ana[[wv]] > 0, ]; dd$.w <- dd[[wv]]
    dsg <- svydesign(id = ~psu, strata = ~stratum, weights = ~.w, nest = TRUE, data = dd)
    cf <- coef(svyglm(pa ~ age, design = dsg)); dd$.acc <- dd$pa - (cf[[1]] + cf[[2]] * dd$age)
    dsg <- svydesign(id = ~psu, strata = ~stratum, weights = ~.w, nest = TRUE, data = dd)
    mm <- svyby(~.acc, ~year, dsg, svymean, vartype = "ci")
    r <- mm[mm$year == 2022, ]
    rec("C phlebotomy weight", "mean phenotypic age acceleration, 2021-2023 (Fig 1B)",
        ifelse(wv == "wt", "MEC weight (as in M02)", "phlebotomy weight in 2021-2023"), r[[2]], r[[3]], r[[4]], NA, sum(dd$year == 2022))
  }
} else {
  message("WTPH2YR could not be downloaded; section C skipped. Check nhanesA / the connection and re-run.")
  ana$wt_ph <- ana$wt
}

## ---------------------------------------------------------------------------
## D. Quartile reassignment, measured and after correction; top-quartile shares
## ---------------------------------------------------------------------------
wq <- function(x, w, probs) { o <- order(x); x <- x[o]; w <- w[o]; cw <- cumsum(w) / sum(w)
  sapply(probs, function(p) x[which(cw >= p)[1]]) }
qcut <- function(x, cuts) cut(x, c(-Inf, cuts, Inf), labels = 1:4, include.lowest = TRUE)
quart <- function(x, yr, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0; x <- x[ok]; yr <- yr[ok]; w <- w[ok]
  pooled <- as.integer(qcut(x, wq(x, w, c(.25, .5, .75)))); within <- integer(length(x))
  for (y in unique(yr)) { i <- yr == y; within[i] <- as.integer(qcut(x[i], wq(x[i], w[i], c(.25, .5, .75)))) }
  pre <- yr < 2013
  c(moving = 100 * sum(w * (pooled != within)) / sum(w),
    q4_pre = 100 * sum(w[pre] * (pooled[pre] == 4)) / sum(w[pre]),
    q4_post = 100 * sum(w[!pre] * (pooled[!pre] == 4)) / sum(w[!pre]))
}
A <- ana[!is.na(ana$pa), ]
## acceleration = residual of phenotypic age on age in the pooled adult sample (as M02)
resid_age <- function(pa, df, w) {
  dd <- data.frame(.pa = pa, age = df$age, psu = df$psu, stratum = df$stratum, .w = w)
  d  <- svydesign(id = ~psu, strata = ~stratum, weights = ~.w, nest = TRUE, data = dd)
  cf <- coef(svyglm(.pa ~ age, design = d)); pa - (cf[[1]] + cf[[2]] * df$age)
}
SERIES <- list(c("Phenotypic age acceleration, measured", "pa", "resid"),
               c("Phenotypic age acceleration, after correction", "pa_corr", "resid"),
               c("RAR, measured", "RAR_obs", "raw"),
               c("RAR, after correction", "RAR_corr", "raw"))
for (wv in c("wt", "wt_ph")) {
  if (wv == "wt_ph" && isTRUE(all.equal(A$wt, A$wt_ph))) next
  lab <- ifelse(wv == "wt", "MEC weight (as in M02)", "phlebotomy weight in 2021-2023")
  ok <- !is.na(A[[wv]]) & A[[wv]] > 0; B <- A[ok, ]; w <- B[[wv]]
  for (k in SERIES) {
    x <- if (k[3] == "resid") resid_age(B[[k[2]]], B, w) else B[[k[2]]]
    q <- quart(x, B$year, w)
    rec("D quartiles", paste0(k[1], ": % of adults changing quartile"), lab, q[["moving"]], NA, NA, NA, sum(!is.na(x)))
    rec("D quartiles", paste0(k[1], ": % in top pooled quartile, before 2013"), lab, q[["q4_pre"]])
    rec("D quartiles", paste0(k[1], ": % in top pooled quartile, 2013 onward"), lab, q[["q4_post"]])
  }
}

res <- bind_rows(OUT)
write.csv(res, file.path(results_dir, "m_verification.csv"), row.names = FALSE)
cat("\n===== VERIFICATION =====\n")
print(res %>% mutate(across(c(est, lo, hi), ~ round(., 3)), p = signif(p, 3)), row.names = FALSE)
cat("\nSaved: results/m_verification.csv\n")

if (!is.null(.nosleep_handle) && !inherits(.nosleep_handle, "try-error"))
  try(NoSleepR::nosleep_off(.nosleep_handle), silent = TRUE)
message("Done: verification complete.")
