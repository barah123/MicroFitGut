# =============================================================================
# MicroFitGut — 13-benchmark-concordance.R
#
# How closely does a reanalysis agree with a published result, and in what way?
#
# The rule this file is built on: do not reduce agreement to one number. A single
# "concordance score" hides the thing you want to know, which is *where* the
# agreement breaks. Alpha diversity can correlate at r = 0.98 per sample and
# still flip the group-level significance; a differential abundance list can
# recover 40% of the named taxa while the effect sizes rank-correlate at 0.9.
# Those are different findings and they license different conclusions.
#
# So each domain gets its own metrics:
#   alpha  per-sample correlation and Bland-Altman agreement, plus whether the
#          group-level test flips sign or significance
#   beta   Procrustes and Mantel between ordinations, plus the delta in
#          PERMANOVA R2
#   DA     Jaccard overlap, rank correlation of effect sizes over the union, and
#          separately the recovery rate of the paper's specific named taxa
#   profile Bray-Curtis and Aitchison distance between the two abundance tables
#
# Requires: mfg_require(c("beta")); utils.R; 12-benchmark-harmonize.R
# =============================================================================

# ── Alpha diversity concordance ──────────────────────────────────────────────

#' Agreement between published and recomputed per-sample alpha diversity.
#'
#' Three separate questions, reported separately:
#'   1. Do the per-sample values agree? (Pearson and Spearman)
#'   2. Is there a systematic offset or scale difference? (Bland-Altman bias and
#'      limits of agreement — a high correlation with a constant offset means the
#'      values are not interchangeable even though they track perfectly)
#'   3. Does the group-level conclusion survive? (sign and significance flip)
#'
#' The third is the one that matters for whether the paper's claim reproduces,
#' and it is the one a correlation coefficient cannot answer.
concordance_alpha <- function(published, reanalysis, group = NULL,
                              metric_name = "alpha", alpha_level = 0.05) {
  shared <- intersect(names(published), names(reanalysis))
  if (length(shared) < 3) {
    stop("Fewer than 3 samples in common (", length(shared), "). ",
         "Per-sample agreement cannot be assessed.", call. = FALSE)
  }
  p <- as.numeric(published[shared]); r <- as.numeric(reanalysis[shared])
  ok <- is.finite(p) & is.finite(r)
  p <- p[ok]; r <- r[ok]; shared <- shared[ok]

  pearson  <- suppressWarnings(stats::cor(p, r, method = "pearson"))
  spearman <- suppressWarnings(stats::cor(p, r, method = "spearman"))

  # Bland-Altman: difference against mean. Bias is the systematic offset; the
  # limits of agreement say how far an individual sample can differ.
  diff_ <- r - p
  mean_ <- (r + p) / 2
  bias <- mean(diff_)
  sd_diff <- stats::sd(diff_)
  loa <- c(bias - 1.96 * sd_diff, bias + 1.96 * sd_diff)
  # Proportional bias: does the difference grow with the magnitude?
  ba_slope <- if (length(p) >= 4) stats::coef(stats::lm(diff_ ~ mean_))[2] else NA_real_

  group_test <- NULL
  if (!is.null(group)) {
    g <- droplevels(as.factor(group[shared]))
    if (length(levels(g)) >= 2) {
      tp <- if (length(levels(g)) == 2)
        stats::wilcox.test(p ~ g, exact = FALSE) else stats::kruskal.test(p ~ g)
      tr <- if (length(levels(g)) == 2)
        stats::wilcox.test(r ~ g, exact = FALSE) else stats::kruskal.test(r ~ g)
      # Direction: which group is higher, by median.
      med_p <- vapply(split(p, g), stats::median, numeric(1))
      med_r <- vapply(split(r, g), stats::median, numeric(1))
      dir_p <- names(which.max(med_p)); dir_r <- names(which.max(med_r))

      sig_p <- tp$p.value < alpha_level; sig_r <- tr$p.value < alpha_level
      group_test <- list(
        test = if (length(levels(g)) == 2) "wilcoxon" else "kruskal",
        p_published = tp$p.value, p_reanalysis = tr$p.value,
        significant_published = sig_p, significant_reanalysis = sig_r,
        significance_flip = sig_p != sig_r,
        highest_group_published = dir_p, highest_group_reanalysis = dir_r,
        direction_flip = !identical(dir_p, dir_r),
        medians_published = med_p, medians_reanalysis = med_r)
    }
  }

  verdict <- if (!is.null(group_test) && group_test$significance_flip) {
      sprintf(paste("The group-level conclusion does NOT reproduce: published",
        "p = %s (%s), reanalysis p = %s (%s). Per-sample values correlate at",
        "r = %.3f, so this is a conclusion flip near the threshold rather than",
        "a measurement disagreement."),
        mfg_fmt_p(group_test$p_published),
        if (group_test$significant_published) "significant" else "not significant",
        mfg_fmt_p(group_test$p_reanalysis),
        if (group_test$significant_reanalysis) "significant" else "not significant",
        pearson)
    } else if (!is.null(group_test) && group_test$direction_flip) {
      sprintf(paste("Significance agrees but the DIRECTION flips: %s was highest",
        "in the published values, %s in the reanalysis. Check group coding and",
        "reference levels before attributing this to a pipeline difference."),
        group_test$highest_group_published, group_test$highest_group_reanalysis)
    } else if (pearson >= 0.95 && abs(bias) < 0.1 * stats::sd(p)) {
      sprintf("Strong agreement (r = %.3f) with negligible bias (%.3f).", pearson, bias)
    } else if (pearson >= 0.9) {
      sprintf(paste("Values track closely (r = %.3f) but with a bias of %.3f",
        "(limits of agreement %.3f to %.3f). They are correlated, not",
        "interchangeable — a systematic offset like this usually means a",
        "different normalization or rarefaction depth."),
        pearson, bias, loa[1], loa[2])
    } else {
      sprintf(paste("Weak per-sample agreement (r = %.3f). The two pipelines are",
        "not measuring the same quantity per sample; look upstream (denoiser,",
        "reference database, filtering) before comparing group-level results."),
        pearson)
    }

  out <- list(metric = metric_name, n = length(p),
    pearson = pearson, spearman = spearman,
    bias = bias, sd_difference = sd_diff, limits_of_agreement = loa,
    bland_altman_slope = ba_slope, group_test = group_test, verdict = verdict,
    values = data.frame(sample = shared, published = p, reanalysis = r,
                        difference = diff_, stringsAsFactors = FALSE))
  mfg_log("benchmark", "concordance_alpha", list(
    metric = metric_name, n = length(p),
    pearson = round(pearson, 4), spearman = round(spearman, 4),
    bias = round(bias, 4),
    significance_flip = group_test$significance_flip %||% NA,
    direction_flip = group_test$direction_flip %||% NA))
  class(out) <- c("mfg_concordance_alpha", "list")
  out
}

#' @export
print.mfg_concordance_alpha <- function(x, ...) {
  cat(sprintf("=== Alpha concordance: %s (n = %d) ===\n", x$metric, x$n))
  cat(sprintf("Pearson r = %.4f | Spearman rho = %.4f\n", x$pearson, x$spearman))
  cat(sprintf("Bland-Altman: bias = %.4f, limits of agreement %.4f to %.4f%s\n",
              x$bias, x$limits_of_agreement[1], x$limits_of_agreement[2],
              if (!is.na(x$bland_antman_slope %||% x$bland_altman_slope) &&
                  abs(x$bland_altman_slope) > 0.1)
                sprintf(" (proportional bias, slope = %.3f)", x$bland_altman_slope) else ""))
  if (!is.null(x$group_test)) {
    gt <- x$group_test
    cat(sprintf("Group test (%s): published p = %s, reanalysis p = %s\n",
                gt$test, mfg_fmt_p(gt$p_published), mfg_fmt_p(gt$p_reanalysis)))
    cat(sprintf("  significance flip: %s | direction flip: %s\n",
                gt$significance_flip, gt$direction_flip))
  }
  cat("\n"); cat(strwrap(x$verdict, width = 78), sep = "\n")
  invisible(x)
}

# ── Beta diversity concordance ───────────────────────────────────────────────

#' Agreement between two ordinations and two PERMANOVA results.
#'
#' Procrustes M2 measures how much the two ordinations differ after optimal
#' rotation, scaling and reflection — the right comparison, because an ordination
#' has no intrinsic axis orientation and comparing coordinates directly would
#' report a difference where there is none. Lower M2 is better; the permutation
#' p-value tests whether the configurations are more similar than chance.
#'
#' Mantel tests correlation between the two distance matrices themselves, which
#' is independent of the ordination step.
#'
#' The delta in PERMANOVA R2 is the effect-size comparison, and the significance
#' flip is the conclusion comparison.
concordance_beta <- function(dist_published, dist_reanalysis,
                             permanova_published = NULL, permanova_reanalysis = NULL,
                             ord_published = NULL, ord_reanalysis = NULL,
                             permutations = 999, seed = 42) {
  set.seed(seed)

  # Restrict both to shared samples in the same order.
  mp <- as.matrix(dist_published); mr <- as.matrix(dist_reanalysis)
  shared <- intersect(rownames(mp), rownames(mr))
  if (length(shared) < 4) {
    stop("Fewer than 4 samples in common; Procrustes and Mantel need more.",
         call. = FALSE)
  }
  dp <- stats::as.dist(mp[shared, shared]); dr <- stats::as.dist(mr[shared, shared])

  mantel <- vegan::mantel(dp, dr, method = "spearman", permutations = permutations)

  # Procrustes on PCoA configurations of each distance matrix, so the comparison
  # does not depend on which ordination method each side happened to use.
  k <- min(4, length(shared) - 1)
  cp <- stats::cmdscale(dp, k = k); cr <- stats::cmdscale(dr, k = k)
  pro <- vegan::procrustes(cp, cr, symmetric = TRUE)
  prot <- vegan::protest(cp, cr, permutations = permutations)

  perm_cmp <- NULL
  if (!is.null(permanova_published) && !is.null(permanova_reanalysis)) {
    tp <- permanova_published$results[[1]]
    tr <- permanova_reanalysis$results[[1]]
    sig_p <- !is.na(tp$p_value) && tp$p_value < 0.05
    sig_r <- !is.na(tr$p_value) && tr$p_value < 0.05
    perm_cmp <- list(
      term = tp$term,
      R2_published = tp$R2, R2_reanalysis = tr$R2,
      delta_R2 = tr$R2 - tp$R2,
      relative_change_R2 = (tr$R2 - tp$R2) / max(tp$R2, 1e-9),
      p_published = tp$p_value, p_reanalysis = tr$p_value,
      significant_published = sig_p, significant_reanalysis = sig_r,
      significance_flip = sig_p != sig_r)
  }

  verdict <- if (!is.null(perm_cmp) && perm_cmp$significance_flip) {
      sprintf(paste("The beta-diversity conclusion does NOT reproduce: PERMANOVA",
        "was %s in the paper (R2 = %.3f, p = %s) and %s here (R2 = %.3f, p = %s).",
        "The underlying distance matrices correlate at Mantel rho = %.3f, so the",
        "flip is in the test, not the distances."),
        if (perm_cmp$significant_published) "significant" else "not significant",
        perm_cmp$R2_published, mfg_fmt_p(perm_cmp$p_published),
        if (perm_cmp$significant_reanalysis) "significant" else "not significant",
        perm_cmp$R2_reanalysis, mfg_fmt_p(perm_cmp$p_reanalysis),
        mantel$statistic)
    } else if (prot$signif < 0.05 && mantel$statistic > 0.8) {
      sprintf(paste("Strong structural agreement: Mantel rho = %.3f, Procrustes",
        "M2 = %.4f (p = %.3f)%s."), mantel$statistic, prot$ss, prot$signif,
        if (!is.null(perm_cmp)) sprintf(", PERMANOVA R2 %.3f vs %.3f (delta %+.3f)",
          perm_cmp$R2_published, perm_cmp$R2_reanalysis, perm_cmp$delta_R2) else "")
    } else if (mantel$statistic > 0.5) {
      sprintf(paste("Moderate agreement (Mantel rho = %.3f, Procrustes M2 = %.4f).",
        "The two pipelines recover a related but not identical sample geometry."),
        mantel$statistic, prot$ss)
    } else {
      sprintf(paste("Poor structural agreement (Mantel rho = %.3f). The distance",
        "matrices themselves disagree, which points upstream of the ordination —",
        "different taxa retained, different normalization, or a different",
        "distance metric than the paper stated."), mantel$statistic)
    }

  out <- list(
    n_samples = length(shared),
    mantel_rho = mantel$statistic, mantel_p = mantel$signif,
    procrustes_m2 = prot$ss, procrustes_correlation = prot$scale,
    procrustes_p = prot$signif, procrustes = pro,
    permanova = perm_cmp, verdict = verdict, permutations = permutations)
  mfg_log("benchmark", "concordance_beta", list(
    n_samples = length(shared),
    mantel_rho = round(mantel$statistic, 4), mantel_p = signif(mantel$signif, 4),
    procrustes_m2 = round(prot$ss, 5), procrustes_p = signif(prot$signif, 4),
    delta_R2 = if (is.null(perm_cmp)) NA else round(perm_cmp$delta_R2, 4),
    significance_flip = perm_cmp$significance_flip %||% NA))
  class(out) <- c("mfg_concordance_beta", "list")
  out
}

#' @export
print.mfg_concordance_beta <- function(x, ...) {
  cat(sprintf("=== Beta concordance (n = %d samples) ===\n", x$n_samples))
  cat(sprintf("Mantel rho = %.4f (p = %s)\n", x$mantel_rho, mfg_fmt_p(x$mantel_p)))
  cat(sprintf("Procrustes M2 = %.5f, correlation = %.4f (p = %s)\n",
              x$procrustes_m2, x$procrustes_correlation, mfg_fmt_p(x$procrustes_p)))
  if (!is.null(x$permanova)) {
    p <- x$permanova
    cat(sprintf("PERMANOVA '%s': R2 %.4f -> %.4f (delta %+.4f, %+.1f%%), p %s -> %s%s\n",
                p$term, p$R2_published, p$R2_reanalysis, p$delta_R2,
                100 * p$relative_change_R2,
                mfg_fmt_p(p$p_published), mfg_fmt_p(p$p_reanalysis),
                if (p$significance_flip) "  [FLIP]" else ""))
  }
  cat("\n"); cat(strwrap(x$verdict, width = 78), sep = "\n")
  invisible(x)
}

# ── Differential abundance concordance ───────────────────────────────────────

#' Agreement between a published hit list and a reanalysis hit list.
#'
#' Reports two things that are routinely conflated:
#'   - set overlap (Jaccard) over the union of both lists
#'   - recovery rate of the paper's specific named taxa
#'
#' They differ whenever the two lists are different sizes. A reanalysis finding
#' 200 taxa including all 20 of the paper's has perfect recovery and a Jaccard of
#' 0.10; reporting only one of those numbers misrepresents the result.
#'
#' Effect-size rank correlation is computed over taxa present in both result
#' tables regardless of significance, which is the more robust comparison —
#' significance is a thresholded quantity and near-threshold taxa flip easily.
concordance_da <- function(published_taxa, reanalysis_taxa,
                           published_effects = NULL, reanalysis_effects = NULL,
                           published_tax = NULL, reanalysis_tax = NULL,
                           all_tested_reanalysis = NULL,
                           synonyms = MFG_TAXON_SYNONYMS) {

  # Recovery: match the paper's named taxa against the reanalysis hits.
  recovery <- match_taxa(published_taxa, reanalysis_taxa,
                         published_tax, reanalysis_tax, synonyms)

  # Was a missed taxon even testable? A taxon the reanalysis filtered out before
  # testing is a different kind of miss from one it tested and did not call.
  not_tested <- character()
  if (!is.null(all_tested_reanalysis)) {
    testable <- match_taxa(recovery$unmatched_published, all_tested_reanalysis,
                           published_tax, reanalysis_tax, synonyms)
    not_tested <- testable$unmatched_published
  }

  pub_h <- stats::na.omit(harmonize_taxon(published_taxa, synonyms))
  rea_h <- stats::na.omit(harmonize_taxon(reanalysis_taxa, synonyms))
  inter <- intersect(pub_h, rea_h)
  un    <- union(pub_h, rea_h)
  jaccard <- if (length(un)) length(inter) / length(un) else NA_real_

  effect_cor <- NA_real_; n_effect_shared <- 0L
  if (!is.null(published_effects) && !is.null(reanalysis_effects)) {
    pe <- published_effects; re <- reanalysis_effects
    names(pe) <- harmonize_taxon(names(pe), synonyms)
    names(re) <- harmonize_taxon(names(re), synonyms)
    pe <- pe[!is.na(names(pe))]; re <- re[!is.na(names(re))]
    sh <- intersect(names(pe), names(re))
    sh <- sh[is.finite(pe[sh]) & is.finite(re[sh])]
    n_effect_shared <- length(sh)
    if (n_effect_shared >= 4) {
      effect_cor <- suppressWarnings(stats::cor(pe[sh], re[sh], method = "spearman"))
      # Sign agreement matters more than magnitude: a taxon enriched in the
      # opposite group is a contradiction, not a smaller effect.
      sign_agree <- mean(sign(pe[sh]) == sign(re[sh]))
    } else sign_agree <- NA_real_
  } else sign_agree <- NA_real_

  verdict <- sprintf(paste("%d of %d published taxa recovered (%.0f%%);",
    "Jaccard over the union of both lists = %.3f (%d shared of %d total)."),
    recovery$n_matched, recovery$n_published, 100 * (recovery$recovery_rate %||% 0),
    jaccard, length(inter), length(un))
  if (length(not_tested)) {
    verdict <- paste(verdict, sprintf(paste("%d of the misses were never tested by",
      "the reanalysis — filtered out before testing — so they are a filtering",
      "difference, not a disagreement about significance."), length(not_tested)))
  }
  if (!is.na(effect_cor)) {
    verdict <- paste(verdict, sprintf(paste("Effect sizes over %d taxa present in",
      "both result tables rank-correlate at rho = %.3f with %.0f%% sign agreement.%s"),
      n_effect_shared, effect_cor, 100 * sign_agree,
      if (!is.na(effect_cor) && effect_cor > 0.7 && (recovery$recovery_rate %||% 1) < 0.6)
        paste(" Effects agree while the hit lists do not, which is the signature of",
              "a significance-threshold or power difference rather than a different",
              "biological answer.") else ""))
  }

  out <- list(
    recovery = recovery,
    n_published = recovery$n_published, n_reanalysis = length(reanalysis_taxa),
    n_recovered = recovery$n_matched,
    recovery_rate = recovery$recovery_rate,
    jaccard = jaccard, n_shared = length(inter), n_union = length(un),
    effect_spearman = effect_cor, n_effect_shared = n_effect_shared,
    sign_agreement = sign_agree,
    missed = recovery$unmatched_published,
    missed_because_not_tested = not_tested,
    novel_in_reanalysis = setdiff(rea_h, pub_h),
    verdict = verdict)
  mfg_log("benchmark", "concordance_da", list(
    n_published = out$n_published, n_reanalysis = out$n_reanalysis,
    n_recovered = out$n_recovered,
    recovery_rate = round(out$recovery_rate %||% NA, 3),
    jaccard = round(jaccard, 3),
    effect_spearman = if (is.na(effect_cor)) NA else round(effect_cor, 3),
    sign_agreement = if (is.na(sign_agree)) NA else round(sign_agree, 3),
    n_missed_not_tested = length(not_tested)))
  class(out) <- c("mfg_concordance_da", "list")
  out
}

#' @export
print.mfg_concordance_da <- function(x, ...) {
  cat("=== Differential abundance concordance ===\n")
  cat(sprintf("Published: %d taxa | Reanalysis: %d taxa | Recovered: %d\n",
              x$n_published, x$n_reanalysis, x$n_recovered))
  cat(sprintf("Recovery rate: %.1f%% | Jaccard: %.3f (%d shared of %d union)\n",
              100 * (x$recovery_rate %||% 0), x$jaccard, x$n_shared, x$n_union))
  if (!is.na(x$effect_spearman)) {
    cat(sprintf("Effect-size Spearman: %.3f over %d shared taxa | sign agreement %.0f%%\n",
                x$effect_spearman, x$n_effect_shared, 100 * x$sign_agreement))
  }
  if (length(x$missed)) {
    cat(sprintf("\nNot recovered (%d): %s\n", length(x$missed),
                paste(utils::head(x$missed, 12), collapse = ", ")))
  }
  if (length(x$missed_because_not_tested)) {
    cat(sprintf("  of which never tested by the reanalysis (%d): %s\n",
                length(x$missed_because_not_tested),
                paste(utils::head(x$missed_because_not_tested, 8), collapse = ", ")))
  }
  if (length(x$novel_in_reanalysis)) {
    cat(sprintf("\nNew in reanalysis (%d): %s\n", length(x$novel_in_reanalysis),
                paste(utils::head(x$novel_in_reanalysis, 12), collapse = ", ")))
  }
  cat("\n"); cat(strwrap(x$verdict, width = 78), sep = "\n")
  invisible(x)
}

# ── Abundance-table concordance ──────────────────────────────────────────────

#' Distance between a published abundance table and the reanalysis table.
#'
#' Only possible when the paper deposited its table, which is the case worth
#' starting a benchmark with, because it isolates the comparison engine from the
#' execution layer.
#'
#' Bray-Curtis per sample says how different the two profiles of the same sample
#' are. Aitchison (Euclidean on CLR) is the compositionally coherent counterpart.
#' Reporting both matters because they disagree in an informative way: a
#' difference confined to the abundant taxa shows in Bray-Curtis, one confined to
#' the rare tail shows in Aitchison.
concordance_profile <- function(alignment, pseudocount = 1e-6) {
  p <- alignment$published; r <- alignment$reanalysis

  bray <- vapply(seq_len(ncol(p)), function(i) {
    a <- p[, i]; b <- r[, i]
    sum(abs(a - b)) / max(sum(a + b), 1e-12)
  }, numeric(1))
  names(bray) <- colnames(p)

  clr <- function(x) { x <- x + pseudocount; log(x) - mean(log(x)) }
  aitch <- vapply(seq_len(ncol(p)), function(i)
    sqrt(sum((clr(p[, i]) - clr(r[, i]))^2)), numeric(1))
  names(aitch) <- colnames(p)

  # Per-taxon agreement: which taxa the two tables disagree most about.
  taxon_cor <- vapply(seq_len(nrow(p)), function(j) {
    a <- p[j, ]; b <- r[j, ]
    if (stats::sd(a) < 1e-12 || stats::sd(b) < 1e-12) return(NA_real_)
    suppressWarnings(stats::cor(a, b, method = "spearman"))
  }, numeric(1))
  names(taxon_cor) <- rownames(p)

  verdict <- if (mean(bray) < 0.1) {
      sprintf(paste("The abundance tables are near-identical (mean per-sample",
        "Bray-Curtis = %.4f). Any downstream disagreement comes from the",
        "statistical layer, not from the profiles."), mean(bray))
    } else if (mean(bray) < 0.3) {
      sprintf(paste("Profiles agree well (mean Bray-Curtis = %.3f, mean Aitchison",
        "= %.2f). Downstream differences are likely amplified by the statistics",
        "rather than caused by the profiles."), mean(bray), mean(aitch))
    } else {
      sprintf(paste("Profiles differ substantially (mean Bray-Curtis = %.3f). The",
        "reanalysis is not reconstructing the published table, so attributing",
        "downstream divergence to the statistical method would be wrong — the",
        "cause is upstream (denoiser, reference database, filtering)."), mean(bray))
    }

  out <- list(
    bray_per_sample = bray, aitchison_per_sample = aitch,
    mean_bray = mean(bray), median_bray = stats::median(bray),
    max_bray = max(bray), worst_sample = names(which.max(bray)),
    mean_aitchison = mean(aitch),
    taxon_spearman = taxon_cor,
    median_taxon_spearman = stats::median(taxon_cor, na.rm = TRUE),
    worst_taxa = names(sort(taxon_cor)[seq_len(min(10, sum(!is.na(taxon_cor))))]),
    n_samples = ncol(p), n_taxa = nrow(p),
    coverage_published = alignment$pct_published_abundance_retained,
    verdict = verdict)
  mfg_log("benchmark", "concordance_profile", list(
    n_samples = ncol(p), n_taxa = nrow(p),
    mean_bray = round(mean(bray), 4), max_bray = round(max(bray), 4),
    mean_aitchison = round(mean(aitch), 3),
    median_taxon_spearman = round(out$median_taxon_spearman, 3),
    coverage_pct = round(out$coverage_published, 1)))
  class(out) <- c("mfg_concordance_profile", "list")
  out
}

#' @export
print.mfg_concordance_profile <- function(x, ...) {
  cat(sprintf("=== Profile concordance (%d taxa x %d samples) ===\n",
              x$n_taxa, x$n_samples))
  cat(sprintf("Bray-Curtis per sample: mean %.4f, median %.4f, max %.4f (%s)\n",
              x$mean_bray, x$median_bray, x$max_bray, x$worst_sample))
  cat(sprintf("Aitchison per sample: mean %.3f\n", x$mean_aitchison))
  cat(sprintf("Per-taxon Spearman: median %.3f\n", x$median_taxon_spearman))
  cat(sprintf("Alignment covered %.1f%% of the published abundance\n", x$coverage_published))
  if (length(x$worst_taxa)) {
    cat(sprintf("Least concordant taxa: %s\n", paste(x$worst_taxa, collapse = ", ")))
  }
  cat("\n"); cat(strwrap(x$verdict, width = 78), sep = "\n")
  invisible(x)
}

# ── Overall verdict ──────────────────────────────────────────────────────────

#' Combine the domain concordances into one verdict, without averaging them.
#'
#' The verdict is categorical, not a score, and it is driven by conclusion flips
#' rather than by correlation magnitudes — because the question a benchmark
#' answers is "does the paper's claim hold", not "how similar are the numbers".
#'
#' reproduced          no conclusion flipped
#' partially_reproduced some claims held, others did not
#' diverged            the central claims did not hold
concordance_verdict <- function(alpha = NULL, beta = NULL, da = NULL,
                                profile = NULL,
                                da_recovery_threshold = 0.7) {
  flips <- character(); held <- character(); notes <- character()

  if (!is.null(alpha)) {
    for (a in if (inherits(alpha, "mfg_concordance_alpha")) list(alpha) else alpha) {
      gt <- a$group_test
      if (is.null(gt)) { notes <- c(notes, sprintf("%s: no group test to compare", a$metric)); next }
      if (isTRUE(gt$significance_flip)) {
        flips <- c(flips, sprintf("alpha (%s): significance flipped (p %s -> %s)",
          a$metric, mfg_fmt_p(gt$p_published), mfg_fmt_p(gt$p_reanalysis)))
      } else if (isTRUE(gt$direction_flip)) {
        flips <- c(flips, sprintf("alpha (%s): direction flipped", a$metric))
      } else {
        held <- c(held, sprintf("alpha (%s): conclusion held (r = %.3f)", a$metric, a$pearson))
      }
    }
  }

  if (!is.null(beta)) {
    p <- beta$permanova
    if (is.null(p)) {
      notes <- c(notes, sprintf("beta: structure agreement only (Mantel rho = %.3f)",
                                beta$mantel_rho))
    } else if (isTRUE(p$significance_flip)) {
      flips <- c(flips, sprintf("beta: PERMANOVA significance flipped (p %s -> %s)",
        mfg_fmt_p(p$p_published), mfg_fmt_p(p$p_reanalysis)))
    } else {
      held <- c(held, sprintf("beta: PERMANOVA conclusion held (R2 %.3f -> %.3f)",
                              p$R2_published, p$R2_reanalysis))
      if (abs(p$relative_change_R2) > 0.5) {
        notes <- c(notes, sprintf(paste("beta: R2 changed by %+.0f%% while staying",
          "significant — the effect size does not reproduce even though the",
          "conclusion does"), 100 * p$relative_change_R2))
      }
    }
  }

  if (!is.null(da)) {
    rr <- da$recovery_rate %||% 0
    if (rr >= da_recovery_threshold) {
      held <- c(held, sprintf("DA: %.0f%% of published taxa recovered", 100 * rr))
    } else {
      flips <- c(flips, sprintf("DA: only %.0f%% of published taxa recovered (%d of %d)",
        100 * rr, da$n_recovered, da$n_published))
    }
  }

  if (!is.null(profile) && profile$mean_bray > 0.3) {
    notes <- c(notes, sprintf(paste("The abundance profiles themselves differ",
      "(mean Bray-Curtis %.3f), so any downstream divergence must be attributed",
      "upstream, not to the statistics."), profile$mean_bray))
  }

  verdict <- if (!length(flips)) "reproduced"
             else if (length(held)) "partially_reproduced"
             else "diverged"

  out <- list(verdict = verdict, flips = flips, held = held, notes = notes,
              n_flips = length(flips), n_held = length(held))
  mfg_log("benchmark", "verdict", list(
    verdict = verdict, n_flips = length(flips), n_held = length(held),
    flips = paste(flips, collapse = " | ")))
  class(out) <- c("mfg_verdict", "list")
  out
}

#' @export
print.mfg_verdict <- function(x, ...) {
  cat("=== Benchmark verdict: ", toupper(gsub("_", " ", x$verdict)), " ===\n", sep = "")
  if (length(x$held)) {
    cat("\nClaims that held:\n")
    for (h in x$held) cat("  + ", h, "\n", sep = "")
  }
  if (length(x$flips)) {
    cat("\nClaims that did not hold:\n")
    for (f in x$flips) cat("  - ", f, "\n", sep = "")
  }
  if (length(x$notes)) {
    cat("\nQualifications:\n")
    for (n in x$notes) cat(strwrap(n, width = 76, prefix = "  * "), sep = "\n")
  }
  invisible(x)
}
