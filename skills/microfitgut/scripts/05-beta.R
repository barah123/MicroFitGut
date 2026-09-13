# =============================================================================
# MicroFitGut — 05-beta.R
#
# Between-sample community structure: distances, ordination, and the tests.
#
# The central trap this file is built around: a significant PERMANOVA does not
# mean the groups differ in location. adonis2 tests whether group membership
# explains variation in the distance matrix, and heterogeneous within-group
# dispersion produces a significant result with no difference in centroid. So
# run_permanova() computes the betadisper check alongside the test and refuses to
# report a location conclusion when dispersion is heterogeneous.
#
# The demos do not emphasise this. Reviewers do. See reference/05 and
# reference/11.
#
# Requires: mfg_require(c("beta")); utils.R; 03-normalize.R
# =============================================================================

MFG_BETA_ANALYSIS_KEY <- c(
  bray = "beta_bray", jaccard = "beta_jaccard", jsd = "beta_bray",
  unifrac = "beta_unifrac", wunifrac = "beta_unifrac",
  euclidean = "beta_euclidean", aitchison = "beta_aitchison"
)

# ── Distances ────────────────────────────────────────────────────────────────

#' Compute a distance matrix, refusing combinations that do not mean anything.
#'
#' Guards two things. Phylogenetic metrics need a real tree — UniFrac on a
#' placeholder topology returns a valid-looking matrix that encodes nothing.
#' And every metric has a normalization it requires; raw counts through
#' Bray-Curtis encode library size as community difference.
compute_distance <- function(ps, method = "bray", check_norm = TRUE) {
  method <- tolower(method)

  if (method %in% c("unifrac", "wunifrac")) {
    real <- mfg_tree_is_real(ps)
    if (!isTRUE(real)) {
      stop("Refusing to compute ", method, ": the tree is ",
           if (is.na(real)) "of unrecorded provenance" else "a placeholder",
           ".\nUniFrac on a tree that was not estimated from the sequence data ",
           "produces a matrix that looks valid and encodes nothing. Confirm with ",
           "mfg_mark_tree_real(ps, TRUE) if the tree is real, or use a ",
           "non-phylogenetic metric: ",
           paste(allowed_beta_distances(FALSE), collapse = ", "), call. = FALSE)
    }
  }

  if (isTRUE(check_norm)) {
    key <- MFG_BETA_ANALYSIS_KEY[[method]] %||% "beta_bray"
    mfg_check_normalization(ps, key)
  }

  d <- if (identical(method, "aitchison")) {
    # Aitchison distance is Euclidean distance on CLR values, by definition.
    if (!identical(mfg_normalization(ps), "clr")) {
      stop("Aitchison distance requires CLR-transformed values. Run clr_transform() first.",
           call. = FALSE)
    }
    stats::dist(mfg_otu_samples_as_rows(ps), method = "euclidean")
  } else {
    phyloseq::distance(ps, method = method)
  }

  mfg_log("beta", "distance_computed",
          list(method = method, normalization = mfg_normalization(ps),
               n_samples = attr(d, "Size"),
               mean_distance = round(mean(d), 4)))
  attr(d, "mfg_method") <- method
  attr(d, "mfg_normalization") <- mfg_normalization(ps)
  d
}

# ── Ordination ───────────────────────────────────────────────────────────────

#' Unconstrained ordination with the diagnostics that make it interpretable.
#'
#' PCoA reports variance explained per axis; NMDS reports stress, which is the
#' only thing that says whether the picture can be trusted. Kruskal's guidance:
#' stress < 0.05 excellent, < 0.10 good, < 0.20 usable, >= 0.20 means the
#' two-dimensional picture is misleading and should not be read as a map.
run_ordination <- function(ps, method = c("PCoA", "NMDS", "PCA", "RDA", "CCA"),
                           distance = "bray", dist_obj = NULL,
                           k = 2, trymax = 100, formula = NULL) {
  method <- match.arg(method)

  if (method %in% c("PCoA", "NMDS")) {
    d <- dist_obj %||% compute_distance(ps, distance)
  }

  if (method == "NMDS") {
    ord <- vegan::metaMDS(d, k = k, trymax = trymax, trace = 0, autotransform = FALSE)
    stress <- ord$stress
    verdict <- if (stress < 0.05) "excellent" else if (stress < 0.10) "good" else
               if (stress < 0.20) "usable, interpret with care" else
               "UNRELIABLE — do not read this plot as a map of community distance"
    mfg_log("beta", "nmds", list(distance = distance, k = k, stress = round(stress, 4),
                                 converged = ord$converged, verdict = verdict))
    if (stress >= 0.20) {
      warning(sprintf(paste("NMDS stress is %.3f (>= 0.20). A %d-dimensional",
        "solution cannot represent these distances; the plot will suggest structure",
        "that is an artefact of the projection. Increase k, or use PCoA, and say",
        "so in the report."), stress, k), call. = FALSE)
    }
    return(structure(list(method = "NMDS", ord = ord, distance = distance,
      stress = stress, stress_verdict = verdict, k = k,
      scores = as.data.frame(vegan::scores(ord, display = "sites")),
      axis_labels = paste0("NMDS", seq_len(k))),
      class = c("mfg_ordination", "list")))
  }

  if (method == "PCoA") {
    ord <- phyloseq::ordinate(ps, method = "PCoA", distance = distance)
    eig <- ord$values$Relative_eig
    pct <- round(100 * eig[seq_len(min(k, length(eig)))], 1)
    mfg_log("beta", "pcoa", list(distance = distance,
      axis1_pct = pct[1], axis2_pct = if (length(pct) > 1) pct[2] else NA,
      cumulative_pct = sum(pct)))
    return(structure(list(method = "PCoA", ord = ord, distance = distance,
      variance_explained = pct, k = k,
      scores = as.data.frame(ord$vectors[, seq_len(min(k, ncol(ord$vectors))), drop = FALSE]),
      axis_labels = sprintf("PCoA%d (%.1f%%)", seq_along(pct), pct)),
      class = c("mfg_ordination", "list")))
  }

  # PCA / RDA / CCA operate on the abundance matrix, not a distance matrix.
  mat <- mfg_otu_samples_as_rows(ps)
  if (method == "PCA") {
    mfg_check_normalization(ps, "ordination_pca")
    ord <- vegan::rda(mat)
  } else {
    if (is.null(formula)) {
      stop(method, " is a constrained ordination and needs a formula naming the ",
           "explanatory variables, e.g. formula = ~ region + wash.\n",
           "Without constraints use method = 'PCA' (unconstrained) instead.",
           call. = FALSE)
    }
    mfg_check_normalization(ps, if (method == "CCA") "ordination_cca" else "ordination_pca")
    meta <- mfg_meta(ps)
    f <- stats::as.formula(paste("mat", paste(deparse(formula), collapse = "")))
    ord <- if (method == "CCA") vegan::cca(f, data = meta) else vegan::rda(f, data = meta)
  }

  s <- summary(ord)
  total_inertia <- ord$tot.chi
  constrained   <- ord$CCA$tot.chi %||% NA_real_
  pct_constrained <- if (!is.na(constrained)) 100 * constrained / total_inertia else NA_real_
  eig <- vegan::eigenvals(ord)
  pct <- round(100 * as.numeric(eig[seq_len(min(k, length(eig)))]) / sum(eig), 1)

  mfg_log("beta", tolower(method), list(
    total_inertia = round(total_inertia, 4),
    constrained_inertia = if (is.na(constrained)) NA else round(constrained, 4),
    pct_constrained = if (is.na(pct_constrained)) NA else round(pct_constrained, 1),
    axis1_pct = pct[1]))

  structure(list(method = method, ord = ord, formula = formula,
    total_inertia = total_inertia, constrained_inertia = constrained,
    pct_constrained = pct_constrained, variance_explained = pct, k = k,
    scores = as.data.frame(vegan::scores(ord, display = "sites", choices = seq_len(k))),
    axis_labels = sprintf("%s%d (%.1f%%)", method, seq_along(pct), pct)),
    class = c("mfg_ordination", "list"))
}

#' @export
print.mfg_ordination <- function(x, ...) {
  cat(sprintf("=== %s ordination ===\n", x$method))
  if (!is.null(x$distance)) cat(sprintf("Distance: %s\n", x$distance))
  if (!is.null(x$stress)) {
    cat(sprintf("Stress: %.4f  -> %s\n", x$stress, x$stress_verdict))
  }
  if (!is.null(x$variance_explained)) {
    cat(sprintf("Variance explained: %s (cumulative %.1f%%)\n",
                paste(sprintf("axis%d %.1f%%", seq_along(x$variance_explained),
                              x$variance_explained), collapse = ", "),
                sum(x$variance_explained)))
  }
  if (!is.null(x$pct_constrained) && !is.na(x$pct_constrained)) {
    cat(sprintf("Constrained inertia: %.4f of %.4f total (%.1f%% explained by %s)\n",
                x$constrained_inertia, x$total_inertia, x$pct_constrained,
                paste(deparse(x$formula), collapse = "")))
    cat(sprintf("  %.1f%% remains unexplained.\n", 100 - x$pct_constrained))
  }
  invisible(x)
}

# ── Dispersion: the check that qualifies every PERMANOVA ─────────────────────

#' Test homogeneity of multivariate dispersion.
#'
#' betadisper reduces each sample to its distance from its group centroid, then
#' tests whether those distances differ by group. This matters because PERMANOVA
#' is sensitive to dispersion as well as location: groups with equal centroids
#' but different spread give a significant adonis2 result.
#'
#' A significant betadisper does not invalidate the PERMANOVA — it changes what
#' the PERMANOVA licenses you to say. "The communities differ" stays defensible;
#' "the community composition shifted between groups" does not, because the
#' difference may be in variability rather than central tendency.
check_dispersion <- function(dist_obj, group, permutations = 999, type = "median") {
  group <- droplevels(as.factor(group))
  bd <- vegan::betadisper(dist_obj, group, type = type)
  pt <- vegan::permutest(bd, permutations = permutations)
  an <- stats::anova(bd)
  p  <- pt$tab[1, "Pr(>F)"]

  centroid_dist <- tapply(bd$distances, bd$group, mean)
  out <- list(
    betadisper = bd, permutest = pt, anova = an,
    p_value = p,
    homogeneous = !is.na(p) && p >= 0.05,
    mean_distance_to_centroid = centroid_dist,
    spread_ratio = if (min(centroid_dist) > 0) max(centroid_dist) / min(centroid_dist) else NA_real_,
    interpretation = if (!is.na(p) && p < 0.05) paste(
      "Dispersion is heterogeneous. A significant PERMANOVA on these data cannot",
      "be reported as a difference in community composition (location), because",
      "unequal within-group variability alone produces that result. Report it as",
      "'the groups differ in their distance-matrix structure' and state which",
      "group is more variable.")
      else paste("Dispersion is homogeneous, so a significant PERMANOVA can be",
                 "interpreted as a difference in group centroids (location).")
  )
  mfg_log("beta", "dispersion_checked", list(
    p_value = signif(p, 4), homogeneous = out$homogeneous,
    spread_ratio = round(out$spread_ratio, 3),
    mean_dist = paste(sprintf("%s=%.3f", names(centroid_dist), centroid_dist), collapse = " ")))
  class(out) <- c("mfg_dispersion", "list")
  out
}

#' @export
print.mfg_dispersion <- function(x, ...) {
  cat("=== Multivariate dispersion (betadisper) ===\n")
  cat(sprintf("Mean distance to centroid: %s\n",
              paste(sprintf("%s=%.4f", names(x$mean_distance_to_centroid),
                            x$mean_distance_to_centroid), collapse = ", ")))
  cat(sprintf("Spread ratio (max/min): %.2f\n", x$spread_ratio))
  cat(sprintf("permutest p = %s -> %s\n", mfg_fmt_p(x$p_value),
              if (x$homogeneous) "homogeneous" else "HETEROGENEOUS"))
  cat(strwrap(x$interpretation, width = 78, prefix = "  "), sep = "\n")
  invisible(x)
}

# ── PERMANOVA ────────────────────────────────────────────────────────────────

#' PERMANOVA with the dispersion check attached and the conclusion qualified.
#'
#' `by = "margin"` is the default because it gives each predictor its own p-value
#' independent of the order terms were entered, which sequential testing does not
#' (Demo 7 notes this as preferable).
#'
#' `strata` blocks the permutations within a level, which is how repeated measures
#' are handled: permuting freely across subjects treats within-subject samples as
#' exchangeable when they are not. If the data has repeated measures and no
#' strata is given, this warns rather than quietly overstating significance.
run_permanova <- function(ps, formula, distance = "bray", dist_obj = NULL,
                          permutations = 999, by = "margin",
                          strata = NULL, check_disp = TRUE, seed = 42) {

  meta <- mfg_meta(ps)
  d <- dist_obj %||% compute_distance(ps, distance)

  rhs   <- attr(stats::terms(formula), "term.labels")
  vars  <- all.vars(formula)
  miss  <- setdiff(vars, names(meta))
  if (length(miss)) {
    stop("Formula variables not in metadata: ", paste(miss, collapse = ", "),
         "\nAvailable: ", paste(names(meta), collapse = ", "), call. = FALSE)
  }

  # Warn when repeated measures are present and permutations are unconstrained.
  rm_det <- mfg_detect_repeated_measures(meta)
  if (isTRUE(rm_det$repeated) && is.null(strata)) {
    warning(sprintf(paste("Repeated measures detected ('%s', up to %d samples per",
      "level) but permutations are unconstrained. This treats within-subject",
      "samples as exchangeable and inflates significance. Pass",
      "strata = '%s' to permute within subject, and see reference/07."),
      rm_det$subject_var, rm_det$max_per_subject, rm_det$subject_var), call. = FALSE)
    mfg_log("beta", "permanova_strata_warning",
            list(subject_var = rm_det$subject_var,
                 note = "unconstrained permutations with repeated measures"))
  }

  set.seed(seed)
  perm_design <- if (!is.null(strata)) {
    permute::how(nperm = permutations, blocks = as.factor(meta[[strata]]))
  } else permutations

  f <- stats::as.formula(paste("d", paste(deparse(formula), collapse = "")))
  fit <- vegan::adonis2(f, data = meta, permutations = perm_design, by = by)

  tab <- as.data.frame(fit)
  term_rows <- setdiff(rownames(tab), c("Residual", "Total"))
  results <- lapply(term_rows, function(tm) list(
    term = tm, R2 = tab[tm, "R2"], F = tab[tm, "F"],
    p_value = tab[tm, "Pr(>F)"], df = tab[tm, "Df"]))
  names(results) <- term_rows

  disp <- NULL
  if (isTRUE(check_disp) && length(term_rows) >= 1) {
    # Dispersion is only defined for a categorical grouping, so check the first
    # categorical term in the formula.
    cat_terms <- term_rows[vapply(term_rows, function(tm)
      tm %in% names(meta) && !is.numeric(meta[[tm]]), logical(1))]
    if (length(cat_terms)) {
      disp <- lapply(cat_terms, function(tm)
        check_dispersion(d, meta[[tm]], permutations = permutations))
      names(disp) <- cat_terms
    }
  }

  out <- list(
    formula = formula, distance = attr(d, "mfg_method") %||% distance,
    normalization = attr(d, "mfg_normalization") %||% mfg_normalization(ps),
    permutations = permutations, by = by, strata = strata,
    table = tab, results = results, dispersion = disp,
    n_samples = attr(d, "Size")
  )
  for (tm in term_rows) {
    mfg_log("beta", "permanova", list(
      term = tm, R2 = round(tab[tm, "R2"], 4), F = round(tab[tm, "F"], 4),
      p_value = signif(tab[tm, "Pr(>F)"], 4), distance = out$distance,
      permutations = permutations, by = by, strata = strata %||% "none"))
  }
  class(out) <- c("mfg_permanova", "list")
  out
}

#' @export
print.mfg_permanova <- function(x, ...) {
  cat(sprintf("=== PERMANOVA: %s ===\n", paste(deparse(x$formula), collapse = "")))
  cat(sprintf("Distance: %s on %s | n = %d | %d permutations | by = %s%s\n",
              x$distance, x$normalization, x$n_samples, x$permutations, x$by,
              if (!is.null(x$strata)) sprintf(" | blocked within %s", x$strata) else ""))
  cat("\n")
  print(x$table)
  cat("\n")
  for (tm in names(x$results)) {
    r <- x$results[[tm]]
    cat(sprintf("%s: R2 = %.4f (%.1f%% of variation), F = %.3f, p = %s\n",
                tm, r$R2, 100 * r$R2, r$F, mfg_fmt_p(r$p_value)))
  }
  if (!is.null(x$dispersion)) {
    for (tm in names(x$dispersion)) {
      dsp <- x$dispersion[[tm]]
      cat(sprintf("\n-- dispersion check for '%s': p = %s (%s)\n", tm,
                  mfg_fmt_p(dsp$p_value),
                  if (dsp$homogeneous) "homogeneous" else "HETEROGENEOUS"))
      r <- x$results[[tm]]
      if (!dsp$homogeneous && !is.na(r$p_value) && r$p_value < 0.05) {
        cat(strwrap(paste("!! PERMANOVA is significant but dispersion is",
          "heterogeneous. Do NOT report this as a shift in community composition.",
          "The significant result is consistent with the groups differing only in",
          "how variable they are. State the dispersion result alongside it."),
          width = 78, prefix = "   "), sep = "\n")
      }
    }
  }
  invisible(x)
}

#' One defensible sentence per term, for the results section.
#'
#' Deliberately refuses to write "community composition differed" when the
#' dispersion check does not license it.
permanova_sentence <- function(x) {
  out <- character()
  for (tm in names(x$results)) {
    r <- x$results[[tm]]
    dsp <- x$dispersion[[tm]]
    sig <- !is.na(r$p_value) && r$p_value < 0.05
    if (!sig) {
      out <- c(out, sprintf(paste("%s did not explain a significant share of",
        "variation in %s distance (PERMANOVA R2 = %.3f, F = %.2f, p = %s,",
        "%d permutations)."), tm, x$distance, r$R2, r$F, mfg_fmt_p(r$p_value),
        x$permutations))
      next
    }
    if (!is.null(dsp) && !dsp$homogeneous) {
      out <- c(out, sprintf(paste("%s was associated with %s distance",
        "(PERMANOVA R2 = %.3f, F = %.2f, p = %s, %d permutations), but",
        "within-group dispersion was heterogeneous (betadisper p = %s), so this",
        "cannot be attributed to a difference in community composition alone;",
        "the groups also differ in variability (mean distance to centroid %s)."),
        tm, x$distance, r$R2, r$F, mfg_fmt_p(r$p_value), x$permutations,
        mfg_fmt_p(dsp$p_value),
        paste(sprintf("%s = %.3f", names(dsp$mean_distance_to_centroid),
                      dsp$mean_distance_to_centroid), collapse = ", ")))
      next
    }
    out <- c(out, sprintf(paste("Community composition differed by %s",
      "(PERMANOVA on %s distance: R2 = %.3f, F = %.2f, p = %s, %d permutations),",
      "with homogeneous within-group dispersion (betadisper p = %s), so the",
      "effect is a difference in group centroids."),
      tm, x$distance, r$R2, r$F, mfg_fmt_p(r$p_value), x$permutations,
      if (is.null(dsp)) "not tested" else mfg_fmt_p(dsp$p_value)))
  }
  out
}

# ── The other two permutation tests ──────────────────────────────────────────

#' ANOSIM.
#'
#' R is bounded -1 to 1. Positive means between-group distances exceed
#' within-group distances, which is the expected direction for real structure. A
#' negative R means samples within a group are more dissimilar than samples from
#' different groups — the opposite of grouping having biological meaning (Demo 7).
#'
#' ANOSIM is more sensitive to dispersion differences than PERMANOVA, so it is
#' reported as a supporting test rather than a primary one.
run_anosim <- function(dist_obj, group, permutations = 999, seed = 42) {
  set.seed(seed)
  group <- droplevels(as.factor(group))
  fit <- vegan::anosim(dist_obj, group, permutations = permutations)
  strength <- if (is.na(fit$statistic)) "undefined"
    else if (fit$statistic < 0) "negative — within-group samples are MORE dissimilar than between-group, the opposite of meaningful grouping"
    else if (fit$statistic < 0.25) "very weak separation"
    else if (fit$statistic < 0.5)  "small to moderate separation"
    else if (fit$statistic < 0.75) "moderate to strong separation"
    else "strong separation"
  mfg_log("beta", "anosim", list(R = round(fit$statistic, 4),
    p_value = signif(fit$signif, 4), permutations = permutations, strength = strength))
  structure(list(fit = fit, R = fit$statistic, p_value = fit$signif,
                 permutations = permutations, strength = strength),
            class = c("mfg_anosim", "list"))
}

#' @export
print.mfg_anosim <- function(x, ...) {
  cat("=== ANOSIM ===\n")
  cat(sprintf("R = %.4f, p = %s (%d permutations)\n  %s\n",
              x$R, mfg_fmt_p(x$p_value), x$permutations, x$strength))
  invisible(x)
}

#' MRPP.
#'
#' delta summarises how dissimilar members of the same group are. A lower delta
#' means samples within groups are similar; the test asks whether the observed
#' delta is smaller than expected by chance (Demo 7, ps7).
run_mrpp <- function(dist_obj, group, permutations = 999, seed = 42) {
  set.seed(seed)
  group <- droplevels(as.factor(group))
  fit <- vegan::mrpp(dist_obj, group, permutations = permutations)
  mfg_log("beta", "mrpp", list(delta = round(fit$delta, 4),
    expected_delta = round(fit$E.delta, 4), A = round(fit$A, 5),
    p_value = signif(fit$Pvalue, 4), permutations = permutations))
  structure(list(fit = fit, delta = fit$delta, expected_delta = fit$E.delta,
                 A = fit$A, p_value = fit$Pvalue, permutations = permutations),
            class = c("mfg_mrpp", "list"))
}

#' @export
print.mfg_mrpp <- function(x, ...) {
  cat("=== MRPP ===\n")
  cat(sprintf("observed delta = %.4f, expected = %.4f, A = %.5f, p = %s (%d permutations)\n",
              x$delta, x$expected_delta, x$A, mfg_fmt_p(x$p_value), x$permutations))
  cat(sprintf("  %s\n", if (x$delta < x$expected_delta)
    "within-group samples are more similar than chance — consistent with real grouping"
    else "within-group samples are no more similar than chance"))
  invisible(x)
}

# ── Clustering ───────────────────────────────────────────────────────────────

#' Hierarchical clustering on a distance matrix.
#'
#' ward.D2 is the default because it is the correct Ward implementation for
#' squared distances and gives the compact clusters the microbiome literature
#' conventionally shows. Linkage choice changes the tree materially, so it is
#' recorded — average and complete are also common (Demo 6).
run_hclust <- function(dist_obj, method = "ward.D2", k = NULL) {
  hc <- stats::hclust(dist_obj, method = method)
  groups <- if (!is.null(k)) stats::cutree(hc, k = k) else NULL
  # The cophenetic correlation says how faithfully the tree represents the
  # original distances. Below ~0.7 the dendrogram is a poor summary.
  coph <- suppressWarnings(stats::cor(dist_obj, stats::cophenetic(hc)))
  mfg_log("beta", "hclust", list(method = method, k = k %||% NA,
    cophenetic_correlation = round(coph, 4),
    fidelity = if (coph >= 0.7) "tree represents distances well"
               else "tree is a poor summary of the distance matrix"))
  structure(list(hclust = hc, method = method, k = k, clusters = groups,
                 cophenetic_correlation = coph),
            class = c("mfg_hclust", "list"))
}

#' Non-hierarchical (k-means) clustering.
run_kmeans <- function(ps, k = 3, nstart = 25, seed = 42, normalize = "total") {
  mat <- mfg_otu_samples_as_rows(ps)
  if (!is.null(normalize)) mat <- vegan::decostand(mat, method = normalize)
  set.seed(seed)
  km <- stats::kmeans(mat, centers = k, nstart = nstart)
  mfg_log("beta", "kmeans", list(k = k, nstart = nstart, seed = seed,
    sizes = paste(km$size, collapse = ","),
    between_ss_pct = round(100 * km$betweenss / km$totss, 1)))
  structure(list(kmeans = km, k = k, clusters = km$cluster,
                 pct_variance_explained = 100 * km$betweenss / km$totss),
            class = c("mfg_kmeans", "list"))
}

#' Agreement between a clustering and a known grouping.
#'
#' Adjusted Rand index: 0 is chance agreement, 1 is identical. Useful for asking
#' whether unsupervised structure recovers the study design.
cluster_vs_group <- function(clusters, group) {
  tab <- table(cluster = clusters, group = group)
  n <- sum(tab)
  sum_comb <- function(x) sum(choose(x, 2))
  index    <- sum_comb(as.vector(tab))
  exp_idx  <- sum_comb(rowSums(tab)) * sum_comb(colSums(tab)) / choose(n, 2)
  max_idx  <- (sum_comb(rowSums(tab)) + sum_comb(colSums(tab))) / 2
  ari <- (index - exp_idx) / (max_idx - exp_idx)
  mfg_log("beta", "cluster_vs_group", list(adjusted_rand_index = round(ari, 4)))
  list(table = tab, adjusted_rand_index = ari)
}
