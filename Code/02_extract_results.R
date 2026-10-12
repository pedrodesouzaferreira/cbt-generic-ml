# ==============================================================================
# Generic ML inference on heterogeneous treatment effects
# Application: Blattman, Jamison & Sheridan (2017), CBT for high-risk men in Liberia
#
# Stage 2 of 2: reporting.
#   Reads Output/gml_A*.rds (from Code/01_run_genericml.R) and writes every
#   table and figure used in the paper:
#     Output/BLP_A*.csv, GATES_A*.csv, CLAN_A*.csv     tidy per-analysis results
#     Output/BLP_all.csv                               combined BLP table
#     Output/BLP_table.tex                             LaTeX table (paper Table 2)
#     Output/BLP_heterogeneity_summary.{pdf,png}       beta_2 across analyses
#     Output/GATES_A*.{pdf,png}, CLAN_A*.{pdf,png}     per-analysis figures
#     Output/GATES_combined_comparison.{pdf,png}       GATES for A1, A4, A5
#
# Usage (from the repository root):
#   Rscript Code/02_extract_results.R
# ==============================================================================

suppressPackageStartupMessages({
  library(GenericML)
  library(ggplot2)
})

get_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f) == 1) return(normalizePath(file.path(dirname(f), "..")))
  if (file.exists("Code/02_extract_results.R")) return(normalizePath("."))
  normalizePath("..")
}
ROOT        <- get_root()
OUTPUT_PATH <- Sys.getenv("GML_OUTPUT", file.path(ROOT, "Output"))

# ------------------------------------------------------------------------------
# Analysis metadata (labels match the paper)
# ------------------------------------------------------------------------------
analyses <- list(
  A1 = list(label = "A1: CBT only, ST ASB",    block = 1,
            D_label = "CBT only vs. control",         Y_label = "Antisocial behavior index (short-term)"),
  A2 = list(label = "A2: CBT only, LT ASB",    block = 1,
            D_label = "CBT only vs. control",         Y_label = "Antisocial behavior index (long-term)"),
  A3 = list(label = "A3: CBT only, ST Theft",  block = 1,
            D_label = "CBT only vs. control",         Y_label = "Number of thefts (short-term)"),
  A4 = list(label = "A4: CBT (any), ST ASB",   block = 2,
            D_label = "CBT (any) vs. cash + control", Y_label = "Antisocial behavior index (short-term)"),
  A5 = list(label = "A5: CBT (any), LT ASB",   block = 2,
            D_label = "CBT (any) vs. cash + control", Y_label = "Antisocial behavior index (long-term)"),
  A6 = list(label = "A6: CBT (any), ST Theft", block = 2,
            D_label = "CBT (any) vs. cash + control", Y_label = "Number of thefts (short-term)")
)

clan_var_labels <- c(
  age_b = "Age", school_b = "Years of schooling",
  wealth_indexstd_b = "Wealth index", homeless_b = "Homeless",
  depression_b = "Depression score", distress_b = "Distress score",
  impulsive_b = "Impulsivity", conscientious_b = "Conscientiousness",
  grit_b = "Grit", selfesteem_b = "Self-esteem",
  patient_game_real_b = "Patience (incentivized game)",
  stealnb_felony_b = "Felony thefts", stealnb_nonviol_b = "Non-violent thefts",
  disputes_all_b = "Disputes", drugssellever_b = "Ever sold drugs",
  warexper_b = "War experience", rel_commanders_b = "Ties to ex-commanders",
  locuscontr_b = "Locus of control", famseeoften_b = "Sees family often"
)

# Group ordering convention used throughout: GenericML sorts by the ML proxy
# for the CATE, so G1 has the LOWEST predicted CATE. Because CBT reduces
# antisocial behavior (negative effects), G1 is the group with the LARGEST
# predicted reduction and G4 the smallest.
group_labels <- c("G1 (largest effect)", "G2", "G3", "G4 (smallest effect)")

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
tidy_blp <- function(blp, label) data.frame(
  analysis = label,
  param    = names(blp$estimate),
  estimate = as.numeric(blp$estimate),
  ci_lower = as.numeric(blp$confidence_interval[, "lower"]),
  ci_upper = as.numeric(blp$confidence_interval[, "upper"]),
  p_value  = as.numeric(blp$p_value),
  learner  = blp$learner,
  stringsAsFactors = FALSE
)

tidy_gates <- function(gates, label) data.frame(
  analysis = label,
  group    = names(gates$estimate),
  estimate = as.numeric(gates$estimate),
  ci_lower = as.numeric(gates$confidence_interval[, "lower"]),
  ci_upper = as.numeric(gates$confidence_interval[, "upper"]),
  p_value  = as.numeric(gates$p_value),
  learner  = gates$learner,
  stringsAsFactors = FALSE
)

# CLAN for one covariate: median over splits of the group means and bands.
tidy_clan_var <- function(arr, variable, label) data.frame(
  analysis = label, variable = variable,
  group    = c(group_labels, "G4 - G1"),
  estimate = as.numeric(apply(arr[, "Estimate", , drop = FALSE], 1, median, na.rm = TRUE)),
  cb_lower = as.numeric(apply(arr[, "CB lower", , drop = FALSE], 1, median, na.rm = TRUE)),
  cb_upper = as.numeric(apply(arr[, "CB upper", , drop = FALSE], 1, median, na.rm = TRUE)),
  stringsAsFactors = FALSE
)

theme_paper <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "gray30"),
          strip.text    = element_text(face = "bold"))
}

save_both <- function(p, stem, width, height) {
  ggsave(file.path(OUTPUT_PATH, paste0(stem, ".pdf")), p, width = width, height = height)
  ggsave(file.path(OUTPUT_PATH, paste0(stem, ".png")), p, width = width, height = height, dpi = 150)
}

# GATES plot: four group effects with 90% CIs, the ATE (beta_1) as a reference
# band. Drawn by hand rather than with plot.GenericML() for a cleaner legend.
plot_gates_custom <- function(gates_tab, blp_tab, title, subtitle) {
  pts <- gates_tab[grepl("^gamma\\.[0-9]+$", gates_tab$group), ]
  pts$group_label <- factor(group_labels[as.integer(sub("gamma\\.", "", pts$group))],
                            levels = group_labels)
  ate <- blp_tab[blp_tab$param == "beta.1", ]

  ggplot(pts, aes(x = group_label, y = estimate, ymin = ci_lower, ymax = ci_upper)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = ate$ci_lower, ymax = ate$ci_upper,
             fill = "#1a6faf", alpha = 0.12) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.6) +
    geom_hline(yintercept = ate$estimate, color = "#1a6faf", linewidth = 0.9) +
    geom_errorbar(width = 0.18, linewidth = 0.9) +
    geom_point(size = 3.2, shape = 15) +
    annotate("text", x = 4.6, y = ate$estimate, label = "ATE", color = "#1a6faf",
             size = 3.2, hjust = 0, fontface = "italic") +
    annotate("text", x = 4.6, y = ate$ci_upper, label = "ATE\n90% CI", color = "#1a6faf",
             size = 2.8, hjust = 0, vjust = 1.1, fontface = "italic") +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.32))) +
    labs(title = title, subtitle = subtitle,
         x = "Impact group, sorted by ML-predicted CATE",
         y = "Group average treatment effect (90% CI)") +
    theme_paper() + theme(legend.position = "none")
}

# ------------------------------------------------------------------------------
# Main loop over analyses
# ------------------------------------------------------------------------------
all_blp <- list(); all_gates <- list(); all_clan <- list(); best <- list()

for (nm in names(analyses)) {
  meta <- analyses[[nm]]
  f <- file.path(OUTPUT_PATH, paste0("gml_", nm, ".rds"))
  if (!file.exists(f)) { cat("Skipping", nm, ": ", f, "not found\n"); next }
  cat("\n--- ", meta$label, " ---\n", sep = "")
  gml <- readRDS(f)
  best[[nm]] <- gml$best

  # BLP
  blp_tab <- tidy_blp(get_BLP(gml, plot = FALSE), meta$label)
  write.csv(blp_tab, file.path(OUTPUT_PATH, paste0("BLP_", nm, ".csv")), row.names = FALSE)
  all_blp[[nm]] <- blp_tab
  print(blp_tab[, c("param", "estimate", "ci_lower", "ci_upper", "p_value")], row.names = FALSE)

  # GATES
  gates_tab <- tidy_gates(get_GATES(gml, plot = FALSE), meta$label)
  write.csv(gates_tab, file.path(OUTPUT_PATH, paste0("GATES_", nm, ".csv")), row.names = FALSE)
  all_gates[[nm]] <- gates_tab

  save_both(
    plot_gates_custom(gates_tab, blp_tab,
                      paste0("GATES: ", meta$D_label), paste0("Outcome: ", meta$Y_label)),
    paste0("GATES_", nm), width = 7.5, height = 5
  )

  # CLAN (best CLAN learner), plus a G4 - G1 bar chart
  clan_raw <- gml$generic_targets[[gml$best$CLAN]]$CLAN
  if (length(clan_raw) > 0) {
    clan_tab <- do.call(rbind, lapply(names(clan_raw), function(v)
      tidy_clan_var(clan_raw[[v]], v, meta$label)))
    write.csv(clan_tab, file.path(OUTPUT_PATH, paste0("CLAN_", nm, ".csv")), row.names = FALSE)
    all_clan[[nm]] <- clan_tab

    d <- clan_tab[clan_tab$group == "G4 - G1", ]
    d$var_label <- ifelse(d$variable %in% names(clan_var_labels),
                          clan_var_labels[d$variable], d$variable)
    d$signif <- d$cb_lower > 0 | d$cb_upper < 0

    p_clan <- ggplot(d, aes(x = reorder(var_label, estimate), y = estimate,
                            ymin = cb_lower, ymax = cb_upper, fill = estimate > 0)) +
      geom_col(show.legend = FALSE, alpha = 0.9) +
      geom_errorbar(width = 0.25, linewidth = 0.4, color = "gray20") +
      scale_fill_manual(values = c("TRUE" = "#2166ac", "FALSE" = "#d73027")) +
      geom_hline(yintercept = 0, linewidth = 0.4) +
      coord_flip() +
      labs(title    = paste0("CLAN: ", meta$D_label),
           subtitle = paste0("Outcome: ", meta$Y_label,
                             "\nBaseline covariate mean, smallest-effect (G4) minus largest-effect (G1) group; 90% CI"),
           x = NULL, y = "Difference in baseline covariate mean (G4 - G1)") +
      theme_paper(12)
    save_both(p_clan, paste0("CLAN_", nm), width = 9, height = 7)
  }

  cat("Best learner:  BLP =", gml$best$BLP, "| GATES =", gml$best$GATES,
      "| CLAN =", gml$best$CLAN, "\n")
}

# ------------------------------------------------------------------------------
# Combined BLP table (CSV + LaTeX)
# ------------------------------------------------------------------------------
blp_all <- do.call(rbind, all_blp)
write.csv(blp_all, file.path(OUTPUT_PATH, "BLP_all.csv"), row.names = FALSE)
cat("\n=== Combined BLP results ===\n")
print(blp_all[, c("analysis", "param", "estimate", "ci_lower", "ci_upper", "p_value")], row.names = FALSE)

fmt_num <- function(x) ifelse(x < 0, sprintf("$%.3f$", x), sprintf("$\\phantom{-}%.3f$", x))
fmt_ci  <- function(lo, hi) if (hi < 0) sprintf("$[%.3f,{%.3f}]$", lo, hi) else sprintf("$[%.3f,\\ %.3f]$", lo, hi)
fmt_p   <- function(p) sprintf("%.3f", p)

rows_tex <- c(); prev_block <- 1
for (nm in names(analyses)) {
  meta <- analyses[[nm]]
  if (meta$block > prev_block) { rows_tex <- c(rows_tex, "\\midrule"); prev_block <- meta$block }
  if (!nm %in% names(all_blp)) {
    rows_tex <- c(rows_tex, sprintf("%s & \\multicolumn{6}{c}{[not run]} \\\\", meta$label)); next
  }
  b  <- all_blp[[nm]]
  b1 <- b[b$param == "beta.1", ]; b2 <- b[b$param == "beta.2", ]
  rows_tex <- c(rows_tex, sprintf(
    "%s & %s & %s & %s & %s & %s & %s \\\\", meta$label,
    fmt_num(b1$estimate), fmt_ci(b1$ci_lower, b1$ci_upper), fmt_p(b1$p_value),
    fmt_num(b2$estimate), fmt_ci(b2$ci_lower, b2$ci_upper), fmt_p(b2$p_value)))
}

pretty_learner <- function(x) c(lasso = "lasso", random_forest = "random forest")[x]
best_blp <- pretty_learner(vapply(best, function(b) b$BLP, character(1)))
learner_note <- if (length(unique(best_blp)) == 1) {
  paste0("Best learner (by BLP criterion): ", best_blp[1], ". ")
} else {
  paste0("Best learner by BLP criterion: ",
         paste(names(best), best_blp, sep = ": ", collapse = "; "), ". ")
}

writeLines(c(
  "\\resizebox{\\textwidth}{!}{",
  "\\begin{tabular}{lcccccc}",
  "\\toprule",
  "& \\multicolumn{3}{c}{$\\hat\\beta_1$ (ATE)} &",
  "  \\multicolumn{3}{c}{$\\hat\\beta_2$ (Heterogeneity)} \\\\",
  "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
  "Analysis & Estimate & 90\\% CI & $p$-value",
  "         & Estimate & 90\\% CI & $p$-value \\\\",
  "\\midrule",
  rows_tex,
  "\\bottomrule",
  "\\end{tabular}",
  "}",
  "\\begin{minipage}{\\textwidth}",
  paste0("\\footnotesize \\textit{Notes:} Median estimates and 90\\% confidence intervals over 100 random sample splits. ",
         "$\\hat\\beta_1$ estimates the ATE; $\\hat\\beta_2$ tests $H_0$: no heterogeneity. ",
         learner_note,
         "Outcome for A1, A2, A4, A5: standardized antisocial behavior index (lower = more prosocial). ",
         "A3, A6: number of thefts in the past two weeks. ",
         "A4--A6 use the full $2\\times2$ factorial design (any CBT vs.\\ cash + control)."),
  "\\end{minipage}"
), file.path(OUTPUT_PATH, "BLP_table.tex"))
cat("Wrote BLP_table.tex\n")

# ------------------------------------------------------------------------------
# Figure: beta_2 (heterogeneity) across the six analyses
# ------------------------------------------------------------------------------
b2 <- blp_all[blp_all$param == "beta.2", ]
to_short <- function(x) gsub(", ", "\n", sub(": ", "\n", x))   # "A1: CBT only, ST ASB" -> 3 lines
b2$short <- factor(to_short(b2$analysis), levels = to_short(sapply(analyses, `[[`, "label")))
b2$sig   <- b2$p_value < 0.10

p_b2 <- ggplot(b2, aes(x = short, y = estimate, ymin = ci_lower, ymax = ci_upper)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(width = 0.2, linewidth = 0.9) +
  geom_point(aes(shape = sig), size = 3.5, color = "#c0392b", show.legend = FALSE) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) +
  labs(title    = "BLP heterogeneity coefficient across analyses",
       subtitle = expression(paste(beta[2], " > 0 indicates the ML proxy predicts the CATE  |  filled point: p < 0.10  |  90% CI")),
       x = NULL, y = expression(hat(beta)[2])) +
  theme_paper()
save_both(p_b2, "BLP_heterogeneity_summary", width = 9, height = 5)

# ------------------------------------------------------------------------------
# Figure: GATES side by side for A1, A4, A5 (antisocial behavior)
# ------------------------------------------------------------------------------
cmp_keys <- c(A1 = "A1: CBT only\n(short-term, narrow)",
              A4 = "A4: any CBT\n(short-term, factorial)",
              A5 = "A5: any CBT\n(long-term, factorial)")
cmp <- do.call(rbind, lapply(names(cmp_keys), function(k) {
  if (!k %in% names(all_gates)) return(NULL)
  t <- all_gates[[k]]; t <- t[grepl("^gamma\\.[0-9]+$", t$group), ]
  t$g <- as.integer(sub("gamma\\.", "", t$group)); t$panel <- cmp_keys[k]; t
}))
cmp$panel <- factor(cmp$panel, levels = cmp_keys)

p_cmp <- ggplot(cmp, aes(x = g, y = estimate, ymin = ci_lower, ymax = ci_upper)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  geom_errorbar(width = 0.2, linewidth = 0.8) +
  geom_point(size = 3, color = "#1f4e79") +
  facet_wrap(~panel, ncol = 3) +
  scale_x_continuous(breaks = 1:4, labels = c("G1\n(largest)", "G2", "G3", "G4\n(smallest)")) +
  labs(title    = "GATES: sorted group average treatment effects on antisocial behavior",
       subtitle = "Groups are quartiles of the ML-predicted CATE; G1 = largest predicted reduction  |  90% CI",
       x = "Impact group", y = "Group average treatment effect (std. dev.)") +
  theme_paper()
save_both(p_cmp, "GATES_combined_comparison", width = 11, height = 5)

cat("\nAll tables and figures written to", OUTPUT_PATH, "\n")
