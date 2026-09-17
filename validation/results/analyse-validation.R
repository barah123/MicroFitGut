# Preliminary benchmark validation: the two bodies of evidence, kept separate.
# Body A: 106 checks against the course material's stated answers (software
#         verification, a known-correct reference exists).
# Body B: 57 claims from 10 published papers (reproducibility assessment, no
#         known-correct answer).
# They are NOT pooled. See README.md in this directory for why.

ROOT <- "/Users/phil/Desktop/Tools/MicroFitGut"
SC   <- file.path(ROOT, "microfitgut-plugin/skills/microfitgut/scripts")
for (f in sort(list.files(SC, pattern = "[.]R$", full.names = TRUE))) suppressWarnings(source(f))
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
OUT <- file.path(ROOT, "microfitgut-plugin/validation/results")
set.seed(42)

# Wilson score interval. Correct near the boundary, where Wald is not: at 10/10
# Wald gives a zero-width interval, which is nonsense.
wilson <- function(x, n, conf = 0.95) {
  if (n == 0) return(c(NA, NA))
  z <- stats::qnorm(1 - (1 - conf) / 2); ph <- x / n
  d <- 1 + z^2 / n
  c(lower = max(0, (ph + z^2/(2*n) - z*sqrt((ph*(1-ph) + z^2/(4*n))/n)) / d),
    upper = min(1, (ph + z^2/(2*n) + z*sqrt((ph*(1-ph) + z^2/(4*n))/n)) / d))
}
fmt_ci <- function(x, n) {
  ci <- wilson(x, n)
  sprintf("%d/%d = %.1f%% (95%% Wilson %.1f-%.1f%%)",
          x, n, 100*x/n, 100*ci[1], 100*ci[2])
}
# Clopper-Pearson: conservative, and the right default when the person reporting
# the validation also wrote the tool.
fmt_cp <- function(x, n) {
  b <- stats::binom.test(x, n)$conf.int
  sprintf("%d/%d = %.1f%% (95%% Clopper-Pearson %.1f-%.1f%%)",
          x, n, 100*x/n, 100*b[1], 100*b[2])
}

# ---------------------------------------------------------------- Body A
cw <- read.csv(file.path(ROOT,
  "output/regression-20260912-193250/tables/regression_results_full.csv"),
  stringsAsFactors = FALSE)
cw$set <- sub("^(ps[0-9]+|quiz[0-9]+).*$", "\\1", cw$source)

cat("=============== BODY A: course material regression ===============\n")
cat("Rows:", nrow(cw), "\n")
print(table(cw$status))

# Denominator. Three kinds of row are not checks and are named, not absorbed:
# INFO has no expected value, BLOCKED is a designed refusal, NODATA is a missing
# input file. A fourth kind IS a check but is not independent: the "ps10 repro"
# rows recompute the same six quantities under the course's own choices in order
# to attribute the discrepancy. Counting a discrepancy and its own explanation
# as separate data points inflates n symmetrically, so the ratio holds steady
# while the sample size becomes fiction.
attribution <- grepl("^ps10 repro", cw$source)
cat("\nrows that are attribution, not independent checks:", sum(attribution), "\n")
scored_A <- subset(cw, status %in% c("MATCH", "DIFFER") & !attribution)
cat("Row level:", fmt_cp(sum(scored_A$status == "MATCH"), nrow(scored_A)), "\n")
cat("  excluded and named separately: INFO", sum(cw$status == "INFO"),
    "| BLOCKED", sum(cw$status == "BLOCKED"), "| NODATA", sum(cw$status == "NODATA"), "\n")

# The honest unit is the exercise question, not the row. Several rows can come
# from one fitted model.
uq <- unique(scored_A[, c("source", "status")])
ux <- tapply(uq$status, uq$source, function(s) if (any(s == "DIFFER")) "DIFFER" else "MATCH")
cat("\nUNIT LEVEL (one unit per exercise question), the headline:\n  ",
    fmt_cp(sum(ux == "MATCH"), length(ux)), "\n")
cat("  questions that disagree:", paste(names(ux)[ux == "DIFFER"], collapse = ", "), "\n")

cat("\nBy problem set:\n")
tabA <- as.data.frame.matrix(table(cw$set, cw$status))
for (k in c("MATCH", "DIFFER", "BLOCKED", "NODATA", "INFO"))
  if (!k %in% names(tabA)) tabA[[k]] <- 0L
tabA <- tabA[, c("MATCH", "DIFFER", "BLOCKED", "NODATA", "INFO")]
tabA$scored <- tabA$MATCH + tabA$DIFFER
tabA$agree_pct <- ifelse(tabA$scored > 0, round(100 * tabA$MATCH / tabA$scored), NA)
tabA <- cbind(set = rownames(tabA), tabA); rownames(tabA) <- NULL
print(tabA)

# The six disagreeing questions trace to three causes.
causes <- data.frame(
  cause = c("DESeq2 size factors: 1 of 51 taxa present in every sample",
            "ps10 Q2 passes an object Q1 had overwritten with normalized counts",
            "Assumption check directs non-parametric; the question instructed ANOVA"),
  questions = c("ps10 Q2, Q3a, Q3b, Q3c", "ps10 Q2, Q3a, Q3b, Q3c", "ps7 Q6, Q8+"),
  set = c("ps10", "ps10", "ps7"), stringsAsFactors = FALSE)
cat("\nThe", sum(ux == "DIFFER"), "disagreeing questions trace to",
    nrow(causes), "causes:\n"); print(causes)

# ---------------------------------------------------------------- Body B
cl <- read.csv(file.path(OUT, "claims-scored.csv"), stringsAsFactors = FALSE)
cl$status <- factor(cl$status, levels = MFG_CLAIM_STATUSES)
cl$scorable   <- cl$status %in% MFG_CLAIM_SCORABLE

cat("\n=============== BODY B: ten published studies ===============\n")
cat("Claims:", nrow(cl), "across", length(unique(cl$study_id)), "studies",
    "(range", paste(range(table(cl$study_id)), collapse = "-"), "per study)\n\n")
print(table(cl$status))

cat("\n-- DESIGNATED HEADLINE CLAIM, one per study --\n")
cat("   Headline claims were designated during this analysis rather than in\n")
cat("   advance, so this proportion is descriptive of the designated set.\n")
pri <- subset(cl, primary == "TRUE")
stopifnot(nrow(pri) == 10, !anyDuplicated(pri$study_id))
print(table(pri$status))
cat("Reproduced:", fmt_ci(sum(pri$status == "reproduced"), nrow(pri)), "\n")

cat("\n-- SECONDARY: all scorable claims --\n")
cat("Reproduced:", fmt_ci(sum(cl$status == "reproduced"), sum(cl$scorable)), "\n")
cat("  denominator excludes", sum(!cl$scorable), "unscorable claims, each named:\n")
print(table(droplevels(cl$status[!cl$scorable])))

cat("\n-- RECONCILIATION AGAINST THE RUN LOGS --\n")
cat("claims in this table:", nrow(cl),
    "| passed through score_claims():", sum(cl$scored_by_function == "TRUE"),
    "| adjudicated in report prose only:", sum(cl$scored_by_function == "FALSE"), "\n")
print(table(cl$study_id, cl$scored_by_function))

cat("\n-- ADJUDICABILITY: could the claim be scored at all? --\n")
cat("Scorable:", fmt_ci(sum(cl$scorable), nrow(cl)), "\n")
adj <- cl %>% group_by(study_id, technology) %>%
  summarise(n = n(), scorable = sum(scorable), .groups = "drop") %>%
  mutate(pct = round(100*scorable/n))
print(as.data.frame(adj))

cat("\n-- 16S vs shotgun --\n")
tech <- cl %>% group_by(technology) %>%
  summarise(claims = n(), scorable = sum(scorable),
            reproduced = sum(status == "reproduced"), .groups = "drop")
print(as.data.frame(tech))
cat("Fisher on reproduced/not among scorable claims:\n")
ft <- fisher.test(table(cl$technology[cl$scorable],
                        cl$status[cl$scorable] == "reproduced"))
cat("  p =", signif(ft$p.value, 3), " OR =", signif(ft$estimate, 3), "\n")
cat("  With 5 studies per arm and claims clustered in studies, this is\n",
    " descriptive. It is not a test of a difference between technologies.\n")

# ---------------------------------------------------------------- robustness
cat("\n=============== ROBUSTNESS OF THE SCORING RULE ===============\n")
d <- subset(cl, !is.na(p_value) & status %in% c("reproduced", "not_reproduced"))
d$p <- as.numeric(d$p_value)
gap_lo <- max(d$p[d$claim_polarity == "difference" & d$status == "reproduced"])
gap_hi <- min(d$p[d$status == "not_reproduced"])
cat(sprintf("No claim has a p-value between %.3f and %.3f, so the 0.05 threshold\n",
            gap_lo, gap_hi))
cat("adjudicates no close call in this set.\n\n")
score_at <- function(alpha) {
  sig <- d$p < alpha
  ifelse(d$claim_polarity == "null", !sig, sig & d$direction_agrees == "TRUE")
}
ref <- score_at(0.05)
sens <- do.call(rbind, lapply(c(0.001, 0.01, 0.05, 0.10, 0.20), function(a)
  data.frame(alpha = a, reproduced = sum(score_at(a)), n = length(ref),
             pct = round(100 * mean(score_at(a))), changed = sum(score_at(a) != ref))))
print(sens)
cat("\nStable between alpha 0.01 and 0.05, where one claim moves. Outside that\n")
cat("range the scored outcome does move, so the threshold is a stated choice.\n")

# ---------------------------------------------------------------- tables out
dir.create(file.path(OUT, "tables"), showWarnings = FALSE, recursive = TRUE)
wr <- function(x, n) {
  write.csv(x, file.path(OUT, "tables", paste0(n, ".csv")), row.names = FALSE); n
}

T1 <- data.frame(
  body = c(rep("A: course material (verification)", 5),
           rep("B: published papers (reproducibility)", 3)),
  quantity = c("agreement per exercise question", "agreement per row",
               "checks with no expected value", "checks refused by a guard",
               "checks with no input data",
               "designated headline claim reproduced",
               "scorable claims reproduced", "claims adjudicable"),
  x = c(sum(ux == "MATCH"), sum(scored_A$status == "MATCH"),
        sum(cw$status == "INFO"), sum(cw$status == "BLOCKED"),
        sum(cw$status == "NODATA"),
        sum(pri$status == "reproduced"), sum(cl$status == "reproduced"),
        sum(cl$scorable)),
  n = c(length(ux), nrow(scored_A), NA, NA, NA,
        nrow(pri), sum(cl$scorable), nrow(cl)),
  stringsAsFactors = FALSE)
T1$pct <- ifelse(is.na(T1$n), NA, round(100 * T1$x / T1$n, 1))
T1$ci_95_clopper_pearson <- vapply(seq_len(nrow(T1)), function(i) {
  if (is.na(T1$n[i])) return(NA_character_)
  b <- stats::binom.test(T1$x[i], T1$n[i])$conf.int
  sprintf("%.1f-%.1f", 100 * b[1], 100 * b[2])
}, character(1))
wr(T1, "T1-headline"); print(T1)

wr(as.data.frame(adj), "T2-adjudicability-by-study")
wr(sens, "T3-threshold-sensitivity")
wr(tabA, "T4-coursework-by-set")
wr(causes, "T5-coursework-causes")
wr(cl %>% filter(!scorable) %>% count(status) %>% rename(claims = n),
   "T6-outside-the-denominator")
cat("\ntables written:", length(list.files(file.path(OUT, "tables"))), "\n")
