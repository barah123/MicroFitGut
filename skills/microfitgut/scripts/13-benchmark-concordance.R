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

  # Jaccard must use the SAME matcher as recovery. Comparing harmonised full
  # names with intersect() only ever finds exact string equality, so a pair that
  # match_taxa() resolved by genus truncation ("Bacteroides" vs "Bacteroides
  # ovatus et rel.") counts as recovered and simultaneously as not shared — two
  # numbers in one output that contradict each other, with the disagreement
  # growing exactly when the two naming schemes differ most.
  pub_h <- stats::na.omit(harmonize_taxon(published_taxa, synonyms))
  rea_h <- stats::na.omit(harmonize_taxon(reanalysis_taxa, synonyms))
  n_inter <- recovery$n_matched %||% 0L
  n_union <- length(pub_h) + length(rea_h) - n_inter
  inter <- if (!is.null(recovery$matches) && nrow(recovery$matches))
             recovery$matches$published else character(0)
  un <- seq_len(max(n_union, 0L))
  jaccard <- if (n_union > 0) n_inter / n_union else NA_real_

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

#' Mean relative abundance per taxon per group.
#'
#' The evidence a `dominance` claim is scored against.
group_mean_abundance <- function(ps, group_var) {
  rel <- phyloseq::transform_sample_counts(ps, function(x) {
    s <- sum(x); if (s > 0) x / s else x })
  mat <- as(phyloseq::otu_table(rel), "matrix")
  if (!phyloseq::taxa_are_rows(rel)) mat <- t(mat)
  g <- as.character(mfg_meta(ps)[[group_var]])
  out <- vapply(split(seq_along(g), g), function(idx)
    rowMeans(mat[, idx, drop = FALSE]), numeric(nrow(mat)))
  as.data.frame(out)
}

# ── Claim status vocabulary ──────────────────────────────────────────────────
#
# The pilot used four statuses and pooled the last two as "unadjudicated". That
# pooled roughly a third of all claims into one bucket whose members meant
# opposite things: a claim the evidence refuses to license is a finding about
# the paper, while a claim whose data was never deposited is a finding about
# data sharing, and a claim the analyst simply did not run is a protocol
# deviation. An endpoint built on the pooled category measures data availability
# rather than reproduction.

#' Statuses that carry evidence about the paper and enter the denominator.
MFG_CLAIM_SCORABLE <- c("reproduced", "not_reproduced", "not_licensed", "not_tested")

#' Statuses that do not. Each is reported as its own count.
MFG_CLAIM_UNSCORABLE <- c("schema_gap", "out_of_scope", "data_absent",
                          "contrast_mismatch", "not_attempted",
                          "pending_adjudication")

MFG_CLAIM_STATUSES <- c(MFG_CLAIM_SCORABLE, MFG_CLAIM_UNSCORABLE)

#' What each status means, for the report.
MFG_CLAIM_STATUS_MEANING <- c(
  reproduced          = "the reanalysis agrees with the paper",
  not_reproduced      = "the reanalysis does not support the claim as stated",
  not_licensed        = "the test performed cannot license the claim (evidence about the paper)",
  not_tested          = "the taxon was filtered out before testing (a filtering difference)",
  schema_gap          = "no claim type expresses this; evidence may exist (about the instrument)",
  out_of_scope        = "not a microbiome measurement; should not have been extracted",
  data_absent         = "the data needed was never deposited (about data sharing)",
  contrast_mismatch   = "the cohort or subgroup cannot be reconstructed from the deposit",
  not_attempted       = "the analyst did not run it (a protocol deviation)",
  pending_adjudication = "qualitative but humanly decidable; awaiting dual coding")

#' Score the paper's stated claims against the reanalysis artifacts.
#'
#' Eight outcomes in two groups. Which group a claim lands in decides whether it
#' enters the endpoint, so the split is the point.
#'
#' Scorable, evidence about the paper:
#'   reproduced      the reanalysis agrees
#'   not_reproduced  it does not support the claim as stated
#'   not_licensed    the test cannot license the claim (heterogeneous dispersion)
#'   not_tested      the taxon was filtered out before testing
#'
#' Unscorable, evidence about something else, reported separately:
#'   schema_gap, out_of_scope, data_absent, contrast_mismatch,
#'   not_attempted, pending_adjudication
#'
#' `contrast_mismatch` is the guard that matters most. A paper claims "ten taxa
#' changed after the diet switch"; a reanalysis of the baseline
#' between-population contrast returns a taxon count too. Scoring one against the
#' other is arithmetic that succeeds and means nothing.
#'
#' One decision rule governs every directional claim: significance first, then
#' direction. Scoring on sign alone holds with probability 0.5 under the null.
score_claims <- function(claims, evidence = list(),
                         synonyms = MFG_TAXON_SYNONYMS) {
  if (inherits(claims, "mfg_study_claim")) claims <- list(claims)
  if (!length(claims)) return(NULL)

  # A study usually makes claims about more than one contrast: a human arm and a
  # mouse arm, a disease comparison and a body-site comparison. Scoring each in
  # its own call produced one scored object per contrast, which then had to be
  # merged by hand before a verdict could be taken, and a hand merge is where a
  # claim quietly goes missing. Passing a NAMED list of evidence sets, each with
  # its own $contrast, scores them together and keeps one denominator.
  if (is.null(evidence$contrast) && length(evidence) &&
      all(vapply(evidence, function(e) is.list(e) && !is.null(e$contrast),
                 logical(1)))) {
    ev_contrasts <- vapply(evidence, function(e) e$contrast, character(1))
    if (anyDuplicated(ev_contrasts)) {
      stop("Two evidence sets describe the same contrast (",
           paste(unique(ev_contrasts[duplicated(ev_contrasts)]), collapse = ", "),
           "). One contrast means one set of results; merge them first.",
           call. = FALSE)
    }
    claim_contrasts <- vapply(claims, function(cl) cl$contrast %||% NA_character_,
                              character(1))
    scorable_claims <- !vapply(claims, function(cl) identical(cl$type, "unscored"),
                               logical(1))
    orphan <- setdiff(stats::na.omit(claim_contrasts[scorable_claims]), ev_contrasts)
    if (length(orphan)) {
      stop("No evidence set was supplied for contrast(s): ",
           paste(orphan, collapse = ", "),
           "\nA scorable claim with no matching evidence would silently become ",
           "contrast_mismatch, which reads as a fact about the deposit rather ",
           "than a gap in this run.", call. = FALSE)
    }
    parts <- lapply(seq_along(evidence), function(i) {
      # Unscored claims carry no contrast and must be scored exactly once.
      keep <- claim_contrasts %in% ev_contrasts[i] |
        (i == 1L & !scorable_claims)
      if (!any(keep)) return(NULL)
      score_claims(claims[keep], evidence = evidence[[i]], synonyms = synonyms)
    })
    out <- do.call(rbind, Filter(Negate(is.null), parts))
    out <- out[match(vapply(claims, function(cl) cl$id, character(1)), out$id), ,
               drop = FALSE]
    rownames(out) <- NULL
    # Each arm logged its own partial count. Log the merged total too, so the
    # run log carries the denominator the verdict was actually taken over.
    mfg_log("benchmark", "claims_scored_merged", c(
      list(n = nrow(out), scorable = sum(out$scorable),
           contrasts = paste(ev_contrasts, collapse = "|")),
      as.list(table(factor(out$status, levels = MFG_CLAIM_STATUSES)))))
    class(out) <- c("mfg_scored_claims", "data.frame")
    return(out)
  }

  ev_contrast <- evidence$contrast %||% NA_character_
  rows <- lapply(claims, function(cl) {
    status <- "data_absent"; observed <- NA_character_; detail <- ""

    if (identical(cl$type, "unscored")) {
      # "unscored" is not one thing. The caller says which, because the
      # distinction is an extraction decision and must be made before results
      # are seen, not inferred afterwards.
      status <- cl$reason %||% "pending_adjudication"
      if (!status %in% MFG_CLAIM_UNSCORABLE) {
        stop(sprintf("Claim '%s': reason '%s' is not one of: %s", cl$id, status,
                     paste(MFG_CLAIM_UNSCORABLE, collapse = ", ")), call. = FALSE)
      }
      detail <- MFG_CLAIM_STATUS_MEANING[[status]]
      return(data.frame(id = cl$id, type = cl$type,
                        contrast = cl$contrast %||% NA_character_,
                        status = status, observed = observed, detail = detail,
                        stringsAsFactors = FALSE))
    }

    if (!is.na(ev_contrast) && !identical(cl$contrast, ev_contrast)) {
      detail <- sprintf("claim is about '%s'; evidence describes '%s'",
                        cl$contrast, ev_contrast)
      return(data.frame(id = cl$id, type = cl$type, contrast = cl$contrast,
                        status = "contrast_mismatch", observed = NA_character_,
                        detail = detail, stringsAsFactors = FALSE))
    }

    if (identical(cl$type, "da_count")) {
      da <- evidence$da
      if (is.null(da)) { detail <- "no DA result supplied" }
      else {
        n_obs <- if ("significant" %in% names(da)) sum(da$significant %in% TRUE) else nrow(da)
        # No default. A silent +/-25% window (floor 2) decided endpoint outcomes
        # in the pilot: at n = 3 it accepted 1 through 5. The tolerance is an
        # analytical choice and belongs in the study card with its justification.
        if (identical(cl$comparator, "approx") && is.null(cl$tolerance)) {
          stop(sprintf(paste("Claim '%s' uses comparator 'approx' with no",
            "tolerance. Set it explicitly at extraction; it decides the outcome",
            "and must be recorded, not defaulted."), cl$id), call. = FALSE)
        }
        tol <- cl$tolerance %||% 0
        ok <- switch(cl$comparator,
          eq = n_obs == cl$n, lte = n_obs <= cl$n, gte = n_obs >= cl$n,
          approx = abs(n_obs - cl$n) <= tol, NA)
        status <- if (isTRUE(ok)) "reproduced" else "not_reproduced"
        observed <- as.character(n_obs)
        detail <- sprintf("paper %s %s, reanalysis %d%s", cl$comparator, cl$n, n_obs,
                          if (identical(cl$comparator, "approx"))
                            sprintf(" (tolerance %d)", tol) else "")
      }

    } else if (identical(cl$type, "da_direction")) {
      da <- evidence$da
      if (is.null(da) || is.null(evidence$positive_effect_group)) {
        detail <- "needs evidence$da and evidence$positive_effect_group"
      } else {
        m <- match_taxa(cl$taxon, as.character(da$taxon), synonyms = synonyms)
        if (!m$n_matched) {
          # Absent from the result table is not a contradiction. Distinguish a
          # taxon that was tested and not called from one that never reached
          # testing, which is a filtering difference. da_null already made this
          # distinction; da_direction did not, which biased the endpoint toward
          # "failed" whenever the reanalysis filtered more aggressively than the
          # paper did.
          tested <- evidence$all_tested
          if (!is.null(tested) &&
              match_taxa(cl$taxon, as.character(tested), synonyms = synonyms)$n_matched) {
            status <- "not_reproduced"
            detail <- "taxon was tested but not called by the reanalysis"
          } else {
            status <- "not_tested"
            detail <- paste("taxon absent from the result table;",
                            if (is.null(tested)) "pass evidence$all_tested to tell tested-and-not-called from never-tested"
                            else "it was filtered out before testing, which is a filtering difference")
          }
        } else {
          hit <- da[as.character(da$taxon) == m$matches$reanalysis[1], , drop = FALSE]
          eff <- hit$effect %||% hit$log2FoldChange %||% NA_real_
          sig <- isTRUE(hit$significant[1])
          obs_group <- if (is.na(eff[1])) NA_character_
                       else if (eff[1] > 0) evidence$positive_effect_group
                       else setdiff(evidence$groups %||% character(0),
                                    evidence$positive_effect_group)[1]
          observed <- obs_group %||% NA_character_
          # ONE rule for every directional claim: significance first, then
          # direction. Sign alone holds with probability 0.5 under the null, and
          # scoring alpha and da_direction differently meant identical evidence
          # scored held under one type and failed under the other.
          status <- if (!sig) "not_reproduced"
                    else if (identical(obs_group, cl$higher_in)) "reproduced"
                    else "not_reproduced"
          detail <- sprintf("paper says higher in %s; reanalysis %s, effect %.3f -> %s",
                            cl$higher_in,
                            if (sig) "significant" else "NOT significant",
                            eff[1], observed)
        }
      }

    } else if (identical(cl$type, "da_null")) {
      da <- evidence$da
      if (is.null(da)) { detail <- "no DA result supplied" }
      else {
        m <- match_taxa(cl$taxon, as.character(da$taxon), synonyms = synonyms)
        if (!m$n_matched) {
          # Absent from the result table is not the same as tested-and-null.
          status <- "not_tested"
          detail <- "taxon not present in the reanalysis result table, so it cannot be confirmed as tested"
        } else {
          hit <- da[as.character(da$taxon) == m$matches$reanalysis[1], , drop = FALSE]
          sig <- isTRUE(hit$significant[1])
          status <- if (!sig) "reproduced" else "not_reproduced"
          observed <- if (sig) "significant" else "not significant"
          detail <- sprintf("paper: tested, not significant; reanalysis: %s", observed)
        }
      }

    } else if (identical(cl$type, "dominance")) {
      ab <- evidence$abundance
      if (is.null(ab) || !cl$group %in% names(ab)) {
        detail <- "needs evidence$abundance with a column for this group"
      } else {
        v <- ab[[cl$group]]; names(v) <- rownames(ab)
        top <- names(sort(v, decreasing = TRUE))[1]
        m <- match_taxa(cl$taxon, names(v), synonyms = synonyms)
        # "Dominated by Bacteroides" where the data splits Bacteroides across
        # several features is a claim about the genus, not about one feature.
        gk <- harmonize_taxon(sub("\\s.*$", "", names(v)), synonyms)
        claim_gk <- harmonize_taxon(sub("\\s.*$", "", cl$taxon), synonyms)
        by_genus <- tapply(v, gk, sum)
        top_genus <- names(sort(by_genus, decreasing = TRUE))[1]
        status <- if (identical(top_genus, claim_gk)) "reproduced" else "not_reproduced"
        observed <- sprintf("%s (%.1f%%)", top_genus, 100 * max(by_genus))
        detail <- sprintf("paper: %s dominates in %s; reanalysis top genus: %s",
                          cl$taxon, cl$group, observed)
      }

    } else if (identical(cl$type, "alpha")) {
      gt <- evidence$alpha
      if (is.null(gt)) { detail <- "no alpha group test supplied" }
      else {
        sig <- isTRUE(gt$significant %||% (!is.na(gt$p) && gt$p < 0.05))
        obs_dir <- if (!sig) "none" else gt$direction %||% NA_character_
        observed <- obs_dir %||% NA_character_
        status <- if (identical(obs_dir, cl$direction)) "reproduced" else "not_reproduced"
        detail <- sprintf("paper: %s; reanalysis: %s", cl$direction, observed)
      }

    } else if (identical(cl$type, "beta")) {
      bt <- evidence$beta
      if (is.null(bt)) { detail <- "no beta test supplied" }
      else {
        sig <- isTRUE(bt$significant %||% (!is.na(bt$p) && bt$p < 0.05))
        status <- if (identical(sig, isTRUE(cl$differs))) "reproduced" else "not_reproduced"
        observed <- if (sig) "differs" else "does not differ"
        unqualified <- status
        # Heterogeneous dispersion means PERMANOVA cannot separate a shift in
        # location from a difference in spread, in EITHER direction: a
        # significant result is not a composition claim, and a null result is
        # not evidence of no difference. That is a finding about the paper, not
        # missing data, so it gets its own status rather than being pooled with
        # claims that were never testable.
        if (isTRUE(bt$dispersion_heterogeneous)) {
          status <- "not_licensed"
          detail <- if (isTRUE(cl$differs))
            "dispersion heterogeneous: a significant PERMANOVA does not license a composition claim"
          else
            "dispersion heterogeneous: a null PERMANOVA is not evidence of no difference"
        } else {
          detail <- sprintf("paper: %s; reanalysis: %s",
                            if (isTRUE(cl$differs)) "differs" else "does not differ", observed)
        }
      }
    }

    if (!status %in% MFG_CLAIM_STATUSES) {
      stop(sprintf("Claim '%s' has type '%s', which score_claims() does not handle. Add a branch for it rather than letting it fail open as unscorable.",
                   cl$id, cl$type), call. = FALSE)
    }
    data.frame(id = cl$id, type = cl$type, contrast = cl$contrast %||% NA_character_,
               status = status, observed = observed, detail = detail,
               stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, rows)
  # Statuses that carry evidence about the PAPER, and therefore enter the
  # endpoint denominator. Everything else is evidence about the instrument, the
  # deposit, or the analyst, and is reported separately with its own count.
  out$scorable <- out$status %in% MFG_CLAIM_SCORABLE
  mfg_log("benchmark", "claims_scored", c(
    list(n = nrow(out), scorable = sum(out$scorable)),
    as.list(table(factor(out$status, levels = MFG_CLAIM_STATUSES)))))
  class(out) <- c("mfg_scored_claims", "data.frame")
  out
}

#' @export
print.mfg_scored_claims <- function(x, ...) {
  cat("=== Claim scoring ===\n")
  for (i in seq_len(nrow(x))) {
    mark <- switch(x$status[i], reproduced = "+", not_reproduced = "-",
                   not_licensed = "!", not_tested = "~", "?")
    cat(sprintf("  %s [%-13s] %-22s %s\n", mark, x$status[i], x$id[i], x$detail[i]))
  }
  n_un <- sum(!(x$status %in% MFG_CLAIM_SCORABLE))
  if (n_un) {
    cat(sprintf("\n  %d of %d claim(s) are outside the endpoint denominator:\n",
                n_un, nrow(x)))
    tb <- table(x$status[!(x$status %in% MFG_CLAIM_SCORABLE)])
    for (k in names(tb)) {
      cat(sprintf("    %-20s %2d  %s\n", k, tb[[k]],
                  MFG_CLAIM_STATUS_MEANING[[k]] %||% ""))
    }
    if ("not_attempted" %in% names(tb)) {
      cat("\n  not_attempted is a protocol deviation, not a property of the paper.\n")
      cat("  It is under the analyst's control and must be reported by name.\n")
    }
  }
  invisible(x)
}

#' Combine the domain concordances into one verdict, without averaging them.
#'
#' The verdict is categorical, not a score, and it is driven by conclusion flips
#' rather than by correlation magnitudes — because the question a benchmark
#' answers is "does the paper's claim hold", not "how similar are the numbers".
#'
#' reproduced          no conclusion flipped
#' partially_reproduced some claims held, others did not
#' diverged            the central claims did not hold
concordance_verdict <- function(claims = NULL, alpha = NULL, beta = NULL,
                                da = NULL, profile = NULL,
                                primary_claim_id = NULL) {
  if (is.null(claims) || !nrow(claims)) {
    stop("concordance_verdict() scores CLAIMS. Pass the output of score_claims().\n",
         "Domain concordance objects (alpha, beta, da, profile) are reported as\n",
         "continuous secondaries and no longer decide the verdict: in the pilot a\n",
         "hard-coded recovery threshold of 0.7 supplied a study's only recorded\n",
         "failure, and it was not one of that paper's claims.", call. = FALSE)
  }
  for (nm in c("alpha", "beta", "da", "profile")) {
    if (!is.null(get(nm))) {
      warning(sprintf(paste("concordance_verdict() ignores '%s' when scoring.",
        "Report it as a prespecified continuous secondary instead."), nm),
        call. = FALSE)
    }
  }

  scorable <- claims[claims$status %in% MFG_CLAIM_SCORABLE, , drop = FALSE]
  excluded <- claims[!(claims$status %in% MFG_CLAIM_SCORABLE), , drop = FALSE]

  # ---- Primary: one pre-designated claim -----------------------------------
  # The unit of inference is the study, so the endpoint is one claim. Pooling
  # made the outcome depend on how finely the abstract was split: claim counts
  # ran 4 to 9 in the pilot, and a "did any claim fail" rule fires with
  # probability 1 - p^k.
  primary <- NULL
  if (!is.null(primary_claim_id)) {
    hit <- claims[claims$id == primary_claim_id, , drop = FALSE]
    if (!nrow(hit)) {
      stop(sprintf("primary_claim_id '%s' is not among the scored claims: %s",
                   primary_claim_id, paste(claims$id, collapse = ", ")), call. = FALSE)
    }
    # not_tested is scorable, but it is evidence about THIS reanalysis's
    # prevalence filter, not about the paper. Left in the H1 denominator it
    # dichotomizes to a non-reproduction, so a filter threshold under the
    # analyst's control would score as a failed claim. It is flagged for
    # exclusion instead, and named, never silently dropped.
    st <- hit$status[1]
    primary <- list(
      id = primary_claim_id, status = st, detail = hit$detail[1],
      adjudicable = st %in% MFG_CLAIM_SCORABLE,
      reproduced = identical(st, "reproduced"),
      h1_eligible = st %in% c("reproduced", "not_reproduced", "not_licensed"),
      exclusion_reason = if (identical(st, "not_tested"))
        paste("headline claim was filtered out before testing; the threshold is",
              "an analyst choice, so this is attrition to report, not a",
              "non-reproduction to count") else NULL)
  }

  # ---- Secondary: proportion of scorable claims reproduced -----------------
  n_scorable <- nrow(scorable)
  n_repro    <- sum(scorable$status == "reproduced")
  proportion <- if (n_scorable) n_repro / n_scorable else NA_real_

  # ---- Descriptive category, explicitly NOT the endpoint -------------------
  # Retained for cross-tabulation only. It is not ordinal: "incomplete" sits on
  # the adjudicability axis rather than the reproduction axis, so the four
  # levels cannot be ordered coherently.
  n_not <- sum(scorable$status %in% c("not_reproduced", "not_licensed"))
  category <- if (!n_scorable) "not_adjudicable"
              else if (n_not && n_repro) "partially_reproduced"
              else if (n_not) "diverged"
              else if (nrow(excluded)) "incomplete"
              else "reproduced"

  counts <- table(factor(claims$status, levels = MFG_CLAIM_STATUSES))
  out <- list(
    primary = primary,
    proportion_reproduced = proportion,
    n_scorable = n_scorable, n_reproduced = n_repro,
    n_excluded = nrow(excluded),
    status_counts = counts,
    category = category, category_is_descriptive_only = TRUE,
    scorable = scorable, excluded = excluded)
  mfg_log("benchmark", "verdict", c(
    list(primary_claim = primary_claim_id %||% NA_character_,
         primary_status = primary$status %||% NA_character_,
         proportion_reproduced = if (is.na(proportion)) NA else round(proportion, 4),
         n_scorable = n_scorable, category = category),
    as.list(counts)))
  class(out) <- c("mfg_verdict", "list")
  out
}

#' @export
print.mfg_verdict <- function(x, ...) {
  if (!is.null(x$primary)) {
    cat("=== PRIMARY ENDPOINT ===\n")
    cat(sprintf("  claim      %s\n", x$primary$id))
    cat(sprintf("  status     %s\n", x$primary$status))
    cat(sprintf("  H1 outcome %s\n",
        if (!isTRUE(x$primary$h1_eligible)) "EXCLUDED from the H1 denominator"
        else if (x$primary$reproduced) "reproduced" else "not reproduced"))
    if (!is.null(x$primary$exclusion_reason)) {
      cat(strwrap(x$primary$exclusion_reason, width = 74, prefix = "  "), sep = "\n")
    }
    cat(sprintf("  %s\n\n", x$primary$detail))
  } else {
    cat("=== No primary claim designated ===\n")
    cat("  Pass primary_claim_id. The endpoint is one pre-designated claim per\n")
    cat("  study, fixed at extraction (protocol section 7).\n\n")
  }
  cat(sprintf("Secondary: %d of %d scorable claims reproduced%s\n",
      x$n_reproduced, x$n_scorable,
      if (is.na(x$proportion_reproduced)) ""
      else sprintf(" (%.0f%%)", 100 * x$proportion_reproduced)))
  nz <- x$status_counts[x$status_counts > 0]
  cat("\nClaim statuses:\n")
  for (k in names(nz)) cat(sprintf("  %-20s %2d\n", k, nz[[k]]))
  if (x$n_excluded) {
    cat(sprintf("\n%d claim(s) outside the denominator. Each is reported with its\n", x$n_excluded))
    cat("own count; they mean different things and are never pooled.\n")
  }
  cat(sprintf("\nDescriptive category: %s (NOT the endpoint; not ordinal)\n", x$category))
  invisible(x)
}
