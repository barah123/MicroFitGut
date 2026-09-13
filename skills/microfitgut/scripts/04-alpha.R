# =============================================================================
# MicroFitGut — 04-alpha.R
#
# Within-sample diversity, and the test-selection logic that goes with it.
#
# Alpha diversity is where an analysis most often looks rigorous and is not. The
# three failure modes this file guards against:
#
#   1. Richness estimators (Chao1, ACE, Fisher) on unequal depth. They count how
#      many taxa were seen once or twice, so on unrarefied data they measure how
#      deeply each sample was sequenced. 03-normalize.R refuses this pairing.
#   2. A parametric test applied without checking its assumptions, or a
#      non-parametric one applied when a parametric test was fine and stronger.
#      check_assumptions() runs the sequence and reports the branch it implies.
#   3. Repeated measures treated as independent. choose_alpha_test() detects the
#      design and says so.
#
# Requires: mfg_require("alpha"); utils.R; 03-normalize.R
# =============================================================================

# ── Which measures are computable ────────────────────────────────────────────

#' Estimators that need integer read counts.
#'
#' Chao1, ACE and Fisher are estimated from how many taxa are seen exactly once
#' or twice, so they need counts. Relative-abundance profiles, as produced by
#' shotgun profilers, carry no such counts; Shannon, Simpson and inverse Simpson
#' are defined on proportions and remain valid.
COUNT_ONLY_MEASURES <- c("Chao1", "ACE", "Fisher")

available_alpha_measures <- function(is_relative) {
  all_m <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
  if (isTRUE(is_relative)) setdiff(all_m, COUNT_ONLY_MEASURES) else all_m
}

#' Alpha diversity for count or proportion data.
#'
#' phyloseq::estimate_richness refuses non-integer input outright, so proportion
#' data is routed through vegan directly for the indices that remain valid.
compute_alpha <- function(ps, measures = c("Observed", "Shannon", "Simpson", "Chao1"),
                          is_relative = NULL) {
  if (is.null(is_relative)) {
    is_relative <- identical(mfg_normalization(ps), "tss") ||
      looks_like_relative_abundance(as(phyloseq::otu_table(ps), "matrix"))
  }
  requested <- measures
  measures  <- intersect(measures, available_alpha_measures(is_relative))
  dropped   <- setdiff(requested, measures)
  if (length(dropped)) {
    message(sprintf("Dropping %s: these estimators need integer read counts.",
                    paste(dropped, collapse = ", ")))
    mfg_log("alpha", "measures_dropped",
            list(dropped = paste(dropped, collapse = ","),
                 reason = "singleton/doubleton-based estimators require counts"))
  }
  if (!length(measures)) measures <- "Shannon"

  mat <- as(phyloseq::otu_table(ps), "matrix")
  if (phyloseq::taxa_are_rows(ps)) mat <- t(mat)   # vegan wants samples as rows

  if (!isTRUE(is_relative)) {
    df <- suppressWarnings(phyloseq::estimate_richness(ps, measures = measures))
    rownames(df) <- phyloseq::sample_names(ps)
    out <- df[, intersect(c(measures, paste0("se.", measures)), colnames(df)), drop = FALSE]
  } else {
    o <- list()
    if ("Observed"   %in% measures) o$Observed   <- rowSums(mat > 0)
    if ("Shannon"    %in% measures) o$Shannon    <- vegan::diversity(mat, index = "shannon")
    if ("Simpson"    %in% measures) o$Simpson    <- vegan::diversity(mat, index = "simpson")
    if ("InvSimpson" %in% measures) o$InvSimpson <- vegan::diversity(mat, index = "invsimpson")
    out <- as.data.frame(o)
    rownames(out) <- rownames(mat)
  }
  mfg_log("alpha", "computed",
          list(measures = paste(colnames(out), collapse = ","),
               n_samples = nrow(out),
               normalization = mfg_normalization(ps)))
  out
}

#' Evenness measures from the microbiome package.
#'
#' Pielou, Simpson, Camargo, Evar and Bulla. Evenness is a different question
#' from richness — a community can be rich and dominated, or poor and even — so
#' reporting both is usually more informative than either alone.
compute_evenness <- function(ps, index = "all") {
  ev <- microbiome::evenness(ps, index = index)
  mfg_log("alpha", "evenness_computed",
          list(indices = paste(colnames(ev), collapse = ",")))
  ev
}

#' Faith's phylogenetic diversity.
#'
#' The total branch length spanning the taxa in a sample. No non-phylogenetic
#' index substitutes for it, which is why picante is a hard requirement rather
#' than optional. Refuses to run on a placeholder tree, because PD on a random
#' topology is a number with no meaning.
faith_pd <- function(ps, include_root = TRUE, require_real_tree = TRUE) {
  tree <- phyloseq::phy_tree(ps, errorIfNULL = FALSE)
  if (is.null(tree)) {
    stop("Faith's PD needs a phylogenetic tree and this object has none.",
         call. = FALSE)
  }
  real <- mfg_tree_is_real(ps)
  if (require_real_tree && !isTRUE(real)) {
    stop("Refusing to compute Faith's PD: the tree is ",
         if (is.na(real)) "of unrecorded provenance" else "marked as a placeholder",
         ".\nPD on a tree that was not estimated from the sequence data is a number ",
         "with no meaning. Confirm with mfg_mark_tree_real(ps, TRUE) if the tree is real.",
         call. = FALSE)
  }
  if (include_root && !ape::is.rooted(tree)) {
    stop("include_root = TRUE requires a rooted tree. Root it first ",
         "(phangorn::midpoint or ape::root), or pass include_root = FALSE.",
         call. = FALSE)
  }
  otu <- as.data.frame(t(as(phyloseq::otu_table(ps), "matrix")))
  if (!phyloseq::taxa_are_rows(ps)) otu <- as.data.frame(as(phyloseq::otu_table(ps), "matrix"))
  pd <- picante::pd(otu, tree, include.root = include_root)
  mfg_log("alpha", "faith_pd",
          list(include_root = include_root, n_samples = nrow(pd),
               tree_tips = length(tree$tip.label)))
  pd
}

#' Assemble alpha diversity alongside the metadata, ready for testing and plots.
#'
#' Returns one data frame: sample, every requested index, and every metadata
#' column. This is the object the tests and plots both take, so the numbers in a
#' figure and the numbers in a test can never diverge.
alpha_table <- function(ps, measures = c("Observed", "Shannon", "Simpson", "Chao1"),
                        include_pd = NULL, include_evenness = FALSE) {
  a <- compute_alpha(ps, measures)

  if (is.null(include_pd)) {
    include_pd <- isTRUE(mfg_tree_is_real(ps))
  }
  if (isTRUE(include_pd)) {
    pd <- faith_pd(ps)
    a <- cbind(a, PD = pd$PD[match(rownames(a), rownames(pd))],
               SR = pd$SR[match(rownames(a), rownames(pd))])
  }
  if (isTRUE(include_evenness)) {
    ev <- compute_evenness(ps)
    ev <- ev[match(rownames(a), rownames(ev)), , drop = FALSE]
    names(ev) <- paste0("evenness_", names(ev))
    a <- cbind(a, ev)
  }

  out <- data.frame(sample = rownames(a), a, stringsAsFactors = FALSE,
                    check.names = FALSE)
  if (!is.null(phyloseq::sample_data(ps, errorIfNULL = FALSE))) {
    md <- mfg_meta(ps)
    md$sample <- rownames(md)
    out <- merge(out, md, by = "sample", all.x = TRUE, sort = FALSE)
  }
  rownames(out) <- out$sample
  attr(out, "mfg_measures") <- setdiff(colnames(a), grep("^se\\.", colnames(a), value = TRUE))
  attr(out, "mfg_normalization") <- mfg_normalization(ps)
  out
}

# ── Assumption testing ───────────────────────────────────────────────────────

#' Run the assumption sequence and say which branch it implies.
#'
#' Shapiro-Wilk for normality, then Bartlett and Fligner-Killeen for equal
#' variances. Bartlett is more powerful but assumes normality itself; Fligner is
#' robust to non-normality, so both are reported and the robust one decides when
#' they disagree (Demo 7, ps7).
#'
#' Note what this does and does not license. Passing Shapiro-Wilk does not prove
#' normality — it fails to reject it, and with small n it has little power to
#' reject anything. The output is a defensible default, not a proof.
check_assumptions <- function(df, response, group, alpha = 0.05) {
  y <- df[[response]]
  g <- as.factor(df[[group]])
  ok <- stats::complete.cases(y, g)
  y <- y[ok]; g <- droplevels(g[ok])

  if (length(unique(g)) < 2) {
    stop("Group variable '", group, "' has fewer than 2 levels after dropping NAs.",
         call. = FALSE)
  }

  n_per <- table(g)
  sw <- if (length(y) >= 3 && length(y) <= 5000) stats::shapiro.test(y) else NULL
  # Per-group normality matters more than pooled for a group comparison.
  sw_by_group <- lapply(split(y, g), function(v)
    if (length(v) >= 3 && length(v) <= 5000) stats::shapiro.test(v) else NULL)
  sw_group_p <- vapply(sw_by_group, function(t) if (is.null(t)) NA_real_ else t$p.value,
                       numeric(1))

  bart <- try(stats::bartlett.test(y ~ g), silent = TRUE)
  flig <- try(stats::fligner.test(y ~ g), silent = TRUE)
  bart_p <- if (inherits(bart, "try-error")) NA_real_ else bart$p.value
  flig_p <- if (inherits(flig, "try-error")) NA_real_ else flig$p.value

  normal   <- !is.null(sw) && sw$p.value > alpha &&
              all(is.na(sw_group_p) | sw_group_p > alpha)
  # Fligner is the robust test, so it decides.
  equal_var <- !is.na(flig_p) && flig_p > alpha

  small_groups <- any(n_per < 5)
  branch <- if (normal && equal_var) "parametric"
            else if (normal && !equal_var) "parametric_welch"
            else "nonparametric"

  out <- list(
    response = response, group = group,
    n_total = length(y), n_per_group = n_per,
    shapiro_p = if (is.null(sw)) NA_real_ else sw$p.value,
    shapiro_p_by_group = sw_group_p,
    bartlett_p = bart_p, fligner_p = flig_p,
    normal = normal, equal_variance = equal_var,
    small_groups = small_groups,
    branch = branch,
    caveat = if (small_groups)
      paste("At least one group has fewer than 5 observations. Shapiro-Wilk has",
            "almost no power at this size, so 'normal' here means 'not rejected',",
            "not 'demonstrated'. Prefer the non-parametric branch.")
      else NA_character_
  )
  mfg_log("alpha", "assumptions_checked", list(
    response = response, group = group,
    shapiro_p = round(out$shapiro_p, 5), fligner_p = round(flig_p, 5),
    bartlett_p = round(bart_p, 5), branch = branch, small_groups = small_groups))
  class(out) <- c("mfg_assumptions", "list")
  out
}

#' @export
print.mfg_assumptions <- function(x, ...) {
  cat(sprintf("Assumptions for %s ~ %s   (n = %d; %s)\n", x$response, x$group,
              x$n_total, paste(sprintf("%s=%d", names(x$n_per_group),
                                       as.integer(x$n_per_group)), collapse = ", ")))
  cat(sprintf("  Shapiro-Wilk (pooled)   p = %s  %s\n", mfg_fmt_p(x$shapiro_p),
              if (isTRUE(x$shapiro_p > 0.05)) "-> normality not rejected" else "-> normality rejected"))
  cat(sprintf("  Shapiro-Wilk per group  %s\n",
              paste(sprintf("%s: %s", names(x$shapiro_p_by_group),
                            mfg_fmt_p(x$shapiro_p_by_group)), collapse = "  ")))
  cat(sprintf("  Bartlett (variance)     p = %s\n", mfg_fmt_p(x$bartlett_p)))
  cat(sprintf("  Fligner-Killeen (robust) p = %s  %s\n", mfg_fmt_p(x$fligner_p),
              if (isTRUE(x$equal_variance)) "-> equal variances not rejected" else "-> variances differ"))
  cat(sprintf("  => branch: %s\n", x$branch))
  if (!is.na(x$caveat)) cat("  ! ", x$caveat, "\n", sep = "")
  invisible(x)
}

# ── Test selection ───────────────────────────────────────────────────────────

#' Decide which group-comparison test to run.
#'
#' Returns the test name plus the reasoning, so the report can state why this
#' test and not another. Returns "friedman" only when a usable block variable
#' describes a complete block design; otherwise the blocked case falls back,
#' because running Kruskal-Wallis on repeated measures treats correlated samples
#' as independent draws and can materially mis-state significance.
#'
#' When repeated measures exist but the design is not a complete block, no
#' simple test is correct — the answer is a mixed model (07-models.R), and this
#' says so rather than picking the least-wrong simple test.
choose_alpha_test <- function(group, block = NULL, block_var = NULL, group_var = NULL,
                              assumptions = NULL) {
  n_levels <- length(unique(group[!is.na(group)]))
  blocked  <- !is.null(block) && !is.null(block_var) && nzchar(block_var %||% "") &&
              !identical(block_var, group_var)

  if (blocked) {
    if (is_complete_block_design(group, block)) {
      return(list(test = "friedman", n_levels = n_levels, blocked = TRUE,
        reason = sprintf(paste("'%s' forms a complete unreplicated block design",
          "(one sample per group x block cell), which is exactly what the Friedman",
          "test requires."), block_var)))
    }
    tab <- table(as.factor(group), as.factor(block))
    return(list(test = "mixed_model_required", n_levels = n_levels, blocked = TRUE,
      reason = sprintf(paste("'%s' indicates repeated measures but the design is",
        "not a complete block (cell counts range %d-%d). No simple test is valid:",
        "Kruskal-Wallis would treat correlated samples as independent and overstate",
        "significance, and Friedman requires exactly one observation per cell. Fit a",
        "mixed model with (1 | %s) instead — see 07-models.R and reference/07."),
        block_var, min(tab), max(tab), block_var)))
  }

  branch <- assumptions$branch %||% "nonparametric"

  if (n_levels == 2) {
    if (identical(branch, "parametric")) {
      return(list(test = "t_test", n_levels = 2, blocked = FALSE,
        reason = "Two independent groups; normality and equal variance were not rejected."))
    }
    if (identical(branch, "parametric_welch")) {
      return(list(test = "welch_t_test", n_levels = 2, blocked = FALSE,
        reason = "Two independent groups; normal but variances differ, so Welch's correction applies."))
    }
    return(list(test = "wilcoxon", n_levels = 2, blocked = FALSE,
      reason = "Two independent groups; normality was rejected or n is too small to assess it."))
  }

  if (identical(branch, "parametric")) {
    return(list(test = "anova", n_levels = n_levels, blocked = FALSE,
      reason = sprintf("%d independent groups; normality and equal variance were not rejected.", n_levels)))
  }
  list(test = "kruskal", n_levels = n_levels, blocked = FALSE,
       reason = sprintf(paste("%d independent groups; normality was rejected or n is",
         "too small to assess it, so the rank-based test applies."), n_levels))
}

# ── Running the test ─────────────────────────────────────────────────────────

#' Run the chosen test and return a tidy result with an effect size.
#'
#' reference/10 requires an effect size, not just a p-value, so one is always
#' computed: eta-squared for Kruskal-Wallis and ANOVA, rank-biserial r for
#' Wilcoxon, Cohen's d for t-tests.
run_alpha_test <- function(df, response, group_var, block_var = NULL,
                           assumptions = NULL, alpha = 0.05) {

  g <- as.factor(df[[group_var]])
  y <- df[[response]]
  b <- if (!is.null(block_var)) df[[block_var]] else NULL
  # complete.cases requires equal-length arguments, so the block is only included
  # when it exists rather than being padded with a scalar.
  ok <- if (is.null(b)) stats::complete.cases(y, g) else stats::complete.cases(y, g, b)
  y <- y[ok]; g <- droplevels(g[ok]); if (!is.null(b)) b <- b[ok]

  if (is.null(assumptions)) {
    assumptions <- check_assumptions(data.frame(y = y, g = g), "y", "g", alpha = alpha)
  }
  choice <- choose_alpha_test(g, b, block_var, group_var, assumptions)

  if (identical(choice$test, "mixed_model_required")) {
    mfg_log("alpha", "test_refused",
            list(response = response, group = group_var, reason = choice$reason))
    return(structure(list(response = response, group = group_var,
      test = "REFUSED", reason = choice$reason, assumptions = assumptions,
      p_value = NA_real_, statistic = NA_real_, effect_size = NA_real_),
      class = c("mfg_alpha_test", "list")))
  }

  n_per <- table(g)
  k     <- length(n_per)
  N     <- length(y)

  res <- switch(choice$test,
    t_test = {
      ft <- stats::t.test(y ~ g, var.equal = TRUE)
      d  <- mfg_cohens_d(y, g)
      list(fit = ft, statistic = unname(ft$statistic), p = ft$p.value,
           effect = d, effect_name = "Cohen's d",
           estimate = unname(diff(rev(ft$estimate))), ci = ft$conf.int)
    },
    welch_t_test = {
      ft <- stats::t.test(y ~ g, var.equal = FALSE)
      d  <- mfg_cohens_d(y, g)
      list(fit = ft, statistic = unname(ft$statistic), p = ft$p.value,
           effect = d, effect_name = "Cohen's d (Welch)",
           estimate = unname(diff(rev(ft$estimate))), ci = ft$conf.int)
    },
    wilcoxon = {
      ft <- stats::wilcox.test(y ~ g, exact = FALSE, conf.int = TRUE)
      list(fit = ft, statistic = unname(ft$statistic), p = ft$p.value,
           effect = mfg_rank_biserial(y, g), effect_name = "rank-biserial r",
           estimate = unname(ft$estimate), ci = ft$conf.int)
    },
    anova = {
      fit <- stats::aov(y ~ g)
      s   <- summary(fit)[[1]]
      ss_between <- s[["Sum Sq"]][1]; ss_total <- sum(s[["Sum Sq"]])
      list(fit = fit, statistic = s[["F value"]][1], p = s[["Pr(>F)"]][1],
           effect = ss_between / ss_total, effect_name = "eta-squared",
           estimate = NA_real_, ci = NULL)
    },
    kruskal = {
      ft <- stats::kruskal.test(y ~ g)
      # eta-squared for H: (H - k + 1) / (N - k), the standard conversion.
      eta <- (unname(ft$statistic) - k + 1) / (N - k)
      list(fit = ft, statistic = unname(ft$statistic), p = ft$p.value,
           effect = max(0, eta), effect_name = "eta-squared (from H)",
           estimate = NA_real_, ci = NULL)
    },
    friedman = {
      ft <- stats::friedman.test(y, g, as.factor(b))
      list(fit = ft, statistic = unname(ft$statistic), p = ft$p.value,
           effect = unname(ft$statistic) / (length(unique(b)) * (k - 1)),
           effect_name = "Kendall's W", estimate = NA_real_, ci = NULL)
    },
    stop("Unhandled test: ", choice$test, call. = FALSE)
  )

  out <- list(
    response = response, group = group_var, block = block_var,
    test = choice$test, test_reason = choice$reason,
    assumptions = assumptions,
    n_per_group = n_per,
    group_medians = vapply(split(y, g), stats::median, numeric(1)),
    group_means   = vapply(split(y, g), mean, numeric(1)),
    statistic = res$statistic, p_value = res$p,
    effect_size = res$effect, effect_name = res$effect_name,
    estimate = res$estimate, conf_int = res$ci,
    significant = !is.na(res$p) && res$p < alpha,
    fit = res$fit
  )
  mfg_log("alpha", "test_run", list(
    response = response, group = group_var, test = choice$test,
    statistic = round(res$statistic, 4), p_value = signif(res$p, 4),
    effect = sprintf("%s=%.4f", res$effect_name, res$effect),
    n = paste(sprintf("%s:%d", names(n_per), as.integer(n_per)), collapse = " ")))
  class(out) <- c("mfg_alpha_test", "list")
  out
}

mfg_cohens_d <- function(y, g) {
  s <- split(y, g)
  if (length(s) != 2) return(NA_real_)
  n1 <- length(s[[1]]); n2 <- length(s[[2]])
  sp <- sqrt(((n1 - 1) * stats::var(s[[1]]) + (n2 - 1) * stats::var(s[[2]])) / (n1 + n2 - 2))
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (mean(s[[2]]) - mean(s[[1]])) / sp
}

mfg_rank_biserial <- function(y, g) {
  s <- split(y, g)
  if (length(s) != 2) return(NA_real_)
  n1 <- length(s[[1]]); n2 <- length(s[[2]])
  U <- suppressWarnings(stats::wilcox.test(s[[1]], s[[2]], exact = FALSE)$statistic)
  1 - (2 * unname(U)) / (n1 * n2)
}

#' @export
print.mfg_alpha_test <- function(x, ...) {
  if (identical(x$test, "REFUSED")) {
    cat(sprintf("=== %s ~ %s : NO VALID SIMPLE TEST ===\n", x$response, x$group))
    cat(x$reason, "\n")
    return(invisible(x))
  }
  cat(sprintf("=== %s ~ %s ===\n", x$response, x$group))
  cat(sprintf("Test: %s\n  because %s\n", x$test, x$test_reason))
  cat(sprintf("n: %s\n", paste(sprintf("%s=%d", names(x$n_per_group),
                                       as.integer(x$n_per_group)), collapse = ", ")))
  cat(sprintf("Medians: %s\n", paste(sprintf("%s=%.3f", names(x$group_medians),
                                             x$group_medians), collapse = ", ")))
  cat(sprintf("Statistic = %.4f, p = %s\n", x$statistic, mfg_fmt_p(x$p_value)))
  cat(sprintf("%s = %.4f\n", x$effect_name, x$effect_size))
  if (!is.null(x$conf_int)) {
    cat(sprintf("95%% CI for the difference: [%.4f, %.4f]\n", x$conf_int[1], x$conf_int[2]))
  }
  cat(sprintf("=> %s at alpha = 0.05\n",
              if (x$significant) "significant" else "not significant"))
  invisible(x)
}

# ── Post-hoc ─────────────────────────────────────────────────────────────────

#' Dunn's test with FDR correction — the standard post-hoc after Kruskal-Wallis.
#'
#' Only run this when the omnibus test was significant. Running post-hoc
#' comparisons after a non-significant omnibus inflates the false positive rate
#' and is the most common multiple-testing error in this literature.
posthoc_dunn <- function(df, response, group_var, method = "bh",
                         omnibus_p = NULL, force = FALSE) {
  if (!is.null(omnibus_p) && omnibus_p >= 0.05 && !force) {
    message(sprintf(paste("Skipping post-hoc: the omnibus test was not significant",
      "(p = %s). Pairwise comparisons after a non-significant omnibus inflate the",
      "false positive rate. Pass force = TRUE to override, and say so in the report."),
      mfg_fmt_p(omnibus_p)))
    mfg_log("alpha", "posthoc_skipped",
            list(response = response, omnibus_p = signif(omnibus_p, 4),
                 reason = "omnibus not significant"))
    return(NULL)
  }
  f <- stats::as.formula(paste(response, "~", group_var))
  res <- FSA::dunnTest(f, data = df, method = method)
  tab <- res$res
  mfg_log("alpha", "posthoc_dunn",
          list(response = response, group = group_var, correction = method,
               n_comparisons = nrow(tab),
               n_significant = sum(tab$P.adj < 0.05, na.rm = TRUE)))
  tab
}

#' Nemenyi post-hoc — the alternative to Dunn, using the Tukey distribution.
posthoc_nemenyi <- function(df, response, group_var, omnibus_p = NULL, force = FALSE) {
  if (!is.null(omnibus_p) && omnibus_p >= 0.05 && !force) {
    message("Skipping post-hoc: omnibus test was not significant.")
    return(NULL)
  }
  res <- DescTools::NemenyiTest(x = df[[response]], g = as.factor(df[[group_var]]),
                               dist = "tukey")
  mfg_log("alpha", "posthoc_nemenyi", list(response = response, group = group_var))
  res
}

#' Tukey HSD after a significant ANOVA.
posthoc_tukey <- function(aov_fit, conf_level = 0.95, omnibus_p = NULL, force = FALSE) {
  if (!is.null(omnibus_p) && omnibus_p >= 0.05 && !force) {
    message("Skipping post-hoc: omnibus ANOVA was not significant.")
    return(NULL)
  }
  res <- stats::TukeyHSD(aov_fit, conf.level = conf_level)
  mfg_log("alpha", "posthoc_tukey", list(conf_level = conf_level))
  res
}

# ── Per-group summary ────────────────────────────────────────────────────────

#' Per-group descriptive table.
#'
#' reference/10 requires n per group after filtering in every summary, and this
#' is the table that satisfies it.
alpha_summary <- function(df, response, group_var) {
  s <- FSA::Summarize(stats::as.formula(paste(response, "~", group_var)), data = df)
  mfg_log("alpha", "summarised", list(response = response, group = group_var))
  s
}

#' Test every index against one grouping variable in one pass.
#'
#' Returns a tidy data frame, one row per index, with the test chosen, the
#' statistic, p, effect size and the reasoning. This is the table the report's
#' alpha section is built from.
alpha_test_all <- function(alpha_df, group_var, block_var = NULL,
                           measures = NULL, alpha = 0.05,
                           correct_across_indices = TRUE) {
  measures <- measures %||% attr(alpha_df, "mfg_measures") %||%
    intersect(c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson",
                "Fisher", "PD"), names(alpha_df))
  results <- list()
  rows <- list()
  for (m in measures) {
    r <- try(run_alpha_test(alpha_df, m, group_var, block_var, alpha = alpha), silent = TRUE)
    if (inherits(r, "try-error")) {
      message(sprintf("Skipping %s: %s", m, conditionMessage(attr(r, "condition"))))
      next
    }
    results[[m]] <- r
    rows[[length(rows) + 1]] <- data.frame(
      index = m, test = r$test, statistic = r$statistic, p_value = r$p_value,
      effect_name = r$effect_name %||% NA_character_, effect_size = r$effect_size,
      shapiro_p = r$assumptions$shapiro_p %||% NA_real_,
      fligner_p = r$assumptions$fligner_p %||% NA_real_,
      reason = r$test_reason, stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(list(table = data.frame(), results = results))
  tab <- do.call(rbind, rows)

  # Several indices tested against one variable is several tests. They are
  # correlated, so BH across them is conservative in the right direction rather
  # than exact, but reporting only raw p across six indices is not defensible.
  if (isTRUE(correct_across_indices) && nrow(tab) > 1) {
    tab$q_value <- stats::p.adjust(tab$p_value, method = "BH")
    mfg_log("alpha", "corrected_across_indices",
            list(n_indices = nrow(tab), method = "BH",
                 note = "indices are correlated, so this is conservative, not exact"))
  }
  tab$significant <- !is.na(tab$p_value) &
    (if ("q_value" %in% names(tab)) tab$q_value < alpha else tab$p_value < alpha)
  list(table = tab, results = results)
}
