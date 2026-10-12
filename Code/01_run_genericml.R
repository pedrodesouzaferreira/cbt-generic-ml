# ==============================================================================
# Generic ML inference on heterogeneous treatment effects
# Application: Blattman, Jamison & Sheridan (2017), CBT for high-risk men in Liberia
# Method:      Chernozhukov, Demirer, Duflo & Fernandez-Val (2025), via the
#              GenericML R package (Welz, Alfons, Demirer & Chernozhukov)
#
# Stage 1 of 2: estimation.
#   Runs the six analyses in the paper and saves one GenericML object per
#   analysis to Output/gml_A*.rds. Tables and figures are produced by
#   Code/02_extract_results.R.
#
#   A1: CBT only vs. pure control        short-term antisocial behavior (ASB)
#   A2: CBT only vs. pure control        long-term ASB
#   A3: CBT only vs. pure control        short-term number of thefts
#   A4: CBT (any) vs. cash + control     short-term ASB      (factorial main effect)
#   A5: CBT (any) vs. cash + control     long-term ASB       (factorial main effect)
#   A6: CBT (any) vs. cash + control     short-term thefts   (factorial main effect)
#
# Usage (from the repository root):
#   Rscript Code/01_run_genericml.R
# Runtime: roughly 20-40 minutes on a laptop (100 splits x 6 analyses x 2 learners).
# ==============================================================================

suppressPackageStartupMessages({
  library(haven)
  library(GenericML)
})

# ------------------------------------------------------------------------------
# 0. Paths (relative to the repository root; override with env vars if needed)
# ------------------------------------------------------------------------------
# Locate the repository root whether the script is run via Rscript or sourced.
get_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f) == 1) return(normalizePath(file.path(dirname(f), "..")))
  if (file.exists("Code/01_run_genericml.R")) return(normalizePath("."))
  normalizePath("..")
}
ROOT <- get_root()

DATA_PATH <- Sys.getenv(
  "STYL_DATA",
  file.path(ROOT, "Data", "BLATTMAN", "20150503_data", "dtafiles", "STYL_Final.dta")
)
OUTPUT_PATH <- Sys.getenv("GML_OUTPUT", file.path(ROOT, "Output"))
dir.create(OUTPUT_PATH, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(DATA_PATH)) {
  stop(
    "Replication data not found at:\n  ", DATA_PATH,
    "\nDownload STYL_Final.dta from the AER / openICPSR replication package for\n",
    "Blattman, Jamison & Sheridan (2017) and place it there, or set the\n",
    "STYL_DATA environment variable to its location. The data are not\n",
    "redistributed in this repository."
  )
}

NUM_SPLITS <- as.integer(Sys.getenv("GML_SPLITS", "100"))
SEED       <- 2024
set.seed(SEED)

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------
cat("Loading data from", DATA_PATH, "\n")
df <- as.data.frame(read_dta(DATA_PATH))

# ------------------------------------------------------------------------------
# 2. Covariates
# ------------------------------------------------------------------------------
# All 57 baseline covariates from Blattman et al.'s balance table (their Table 1).
base_vars <- c(
  # Demographics
  "age_b", "livepartner_b", "mpartners_b", "hhunder15_b", "famseeoften_b",
  "muslim_b", "school_b", "schoolbasin_b",
  # Health & cognition
  "literacy_b", "mathscore_b", "health_resc_b", "disabled_b",
  "depression_b", "distress_b",
  # Conflict exposure
  "rel_commanders_b", "faction_b", "warexper_b",
  # Economic
  "profitsump99avg7d_b", "wealth_indexstd_b", "homeless_b",
  "slphungry7dx_b", "savstockp99_b", "loan50_b", "loan300_b",
  # Labor
  "illicit7da_zero_b", "agricul7da_zero_b", "nonagwage7da_zero_b",
  "allbiz7da_zero_b", "nonaghigh7da_zero_b",
  "agriculeveramt_b", "nonagbizeveramt_b", "nonaghigheveramt_b",
  # Substance use
  "drugssellever_b", "drinkboozeself_b", "druggrassself_b",
  "grassdailyuser_b", "harddrugsever_b", "harddrugsdailyuser_b",
  # Crime & violence
  "steals_b", "stealnb_nonviol_b", "stealnb_felony_b",
  "disputes_all_b", "asbhostil_b",
  # Personality / noncognitive
  "conscientious_b", "neurotic_b", "grit_b", "rewardresp_b",
  "locuscontr_b", "impulsive_b", "selfesteem_b",
  # Time & risk preferences
  "patient_game_real_b", "inconsistent_game_resc_b", "risk_game_resc_b",
  "timedecl_b", "riskdecl_b",
  # Executive function
  "cognitive_score_b", "ef_score_b"
)

# CLAN: an interpretable subset of 19 covariates used to characterise the
# most- and least-affected groups.
clan_vars <- c(
  "age_b", "school_b", "wealth_indexstd_b", "homeless_b",
  "depression_b", "distress_b", "impulsive_b", "conscientious_b",
  "grit_b", "selfesteem_b", "patient_game_real_b",
  "stealnb_felony_b", "stealnb_nonviol_b", "disputes_all_b",
  "drugssellever_b", "warexper_b", "rel_commanders_b",
  "locuscontr_b", "famseeoften_b"
)

missing_vars <- setdiff(c(base_vars, clan_vars), names(df))
if (length(missing_vars) > 0) {
  stop("Variables not found in the data: ", paste(missing_vars, collapse = ", "))
}

# ------------------------------------------------------------------------------
# 3. Helpers
# ------------------------------------------------------------------------------
# Median imputation of the (small) share of missing baseline values.
impute_median <- function(Z) {
  for (v in colnames(Z)) {
    miss <- is.na(Z[, v])
    if (any(miss)) Z[miss, v] <- median(Z[, v], na.rm = TRUE)
  }
  storage.mode(Z) <- "double"
  Z
}

run_analysis <- function(data, D_var, Y_var, label,
                         learners   = c("lasso", "random_forest"),
                         num_splits = NUM_SPLITS) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("Analysis ", label, ":  D = ", D_var, "   Y = ", Y_var, "\n", sep = "")

  s <- data[!is.na(data[[D_var]]) & !is.na(data[[Y_var]]), ]
  D <- as.numeric(s[[D_var]])
  Y <- as.numeric(s[[Y_var]])
  cat("N = ", nrow(s), "  (treated ", sum(D), ", control ", sum(1 - D), ")\n", sep = "")
  cat("Propensity score (known, constant): ", round(mean(D), 3), "\n", sep = "")
  cat("Running GenericML with ", num_splits, " splits ...\n", sep = "")

  GenericML(
    Z                        = impute_median(as.matrix(s[, base_vars])),
    D                        = D,
    Y                        = Y,
    learners_GenericML       = learners,
    learner_propensity_score = "constant",   # RCT: propensity = sample share treated
    num_splits               = num_splits,
    Z_CLAN                   = impute_median(as.matrix(s[, clan_vars])),
    HT                       = FALSE,
    quantile_cutoffs         = c(0.25, 0.5, 0.75),   # GATES quartiles
    significance_level       = 0.05,
    parallel                 = FALSE,
    seed                     = SEED
  )
}

# ------------------------------------------------------------------------------
# 4. Samples
# ------------------------------------------------------------------------------
# round 1 = endline ~5 weeks after therapy (short-term)
# round 5 = ~12-13 months after therapy (long-term)
df_r1 <- df[df$round == 1, ]
df_r5 <- df[df$round == 5, ]

# Narrow block: CBT-only arm vs. pure control.
narrow <- function(d) d[d$tpassonly == 1 | d$control == 1, ]

# Factorial block: any CBT (CBT-only or CBT+cash) vs. no CBT (cash-only or control).
# Cash assignment is orthogonal to CBT assignment by design, so this is the
# CBT main effect averaged over cash, using the full sample.
factorial <- function(d) {
  d$D_cbt <- as.numeric(d$tpassonly == 1 | d$tpcashass == 1)
  d
}

# ------------------------------------------------------------------------------
# 5. Run the six analyses and save the GenericML objects
# ------------------------------------------------------------------------------
analyses <- list(
  A1 = list(data = narrow(df_r1),    D = "tpassonly", Y = "fam_asb_st"),
  A2 = list(data = narrow(df_r5),    D = "tpassonly", Y = "fam_asb_lt"),
  A3 = list(data = narrow(df_r1),    D = "tpassonly", Y = "stealnb_stav"),
  A4 = list(data = factorial(df_r1), D = "D_cbt",     Y = "fam_asb_st"),
  A5 = list(data = factorial(df_r5), D = "D_cbt",     Y = "fam_asb_lt"),
  A6 = list(data = factorial(df_r1), D = "D_cbt",     Y = "stealnb_stav")
)

for (nm in names(analyses)) {
  a   <- analyses[[nm]]
  gml <- run_analysis(a$data, a$D, a$Y, nm)
  saveRDS(gml, file.path(OUTPUT_PATH, paste0("gml_", nm, ".rds")))

  blp <- get_BLP(gml, plot = FALSE)
  cat(sprintf("BLP  beta_1 (ATE) = %6.3f  [p = %.3f]   beta_2 (HET) = %6.3f  [p = %.3f]\n",
              blp$estimate["beta.1"], blp$p_value["beta.1"],
              blp$estimate["beta.2"], blp$p_value["beta.2"]))
  cat("Best learner by target:  BLP =", gml$best$BLP,
      "| GATES =", gml$best$GATES, "| CLAN =", gml$best$CLAN, "\n")
}

cat("\nDone. GenericML objects saved to", OUTPUT_PATH, "\n")
cat("Next: Rscript Code/02_extract_results.R\n")
