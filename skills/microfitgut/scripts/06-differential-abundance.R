# =============================================================================
# MicroFitGut — 06-differential-abundance.R
#
# Which taxa differ between groups. The stage with the widest spread between
# methods on the same data, and therefore the one where method choice has to be
# made before looking at results and stated plainly afterwards.
#
# Four methods, each with a different model of the data:
#   ALDEx2     Dirichlet-multinomial Monte Carlo, then CLR, then a rank or t test.
#              Treats compositionality head-on and reports dispersion.
#   ANCOM-BC2  log-linear model estimating true absolute abundance with bias
#              correction; detects structural zeros; handles random effects.
#   DESeq2     negative binomial with shrunken dispersion. Most power at small n,
#              least compositional awareness.
#   Wilcoxon / Kruskal-Wallis on CLR values. Simple, assumption-light, weakest.
#
# The decision table is reference/06. The honest summary: they disagree, the
# disagreement is informative, and reporting the overlap across two methods is
# more defensible than picking the one with the most hits.
#
# Requires: mfg_require("da"); utils.R; 03-normalize.R
# =============================================================================

# ── Method selection ─────────────────────────────────────────────────────────

#' Recommend a DA method from the design, with the reasoning.
#'
#' Does not pick silently — returns the recommendation, the alternatives, and why.
#' The agent states this in the report; the user can override.
choose_da_method <- function(ps, group_var, n_groups = NULL, subject_var = NULL,
                             covariates = NULL) {
  meta <- mfg_meta(ps)
  g <- meta[[group_var]]
  n_groups <- n_groups %||% length(unique(g[!is.na(g)]))
  n_per <- table(g)
  n_min <- min(n_per)
  n_taxa <- phyloseq::ntaxa(ps)
  zero_prop <- mean(mfg_otu_taxa_as_rows(ps) == 0)
  rm_det <- mfg_detect_repeated_measures(meta, subject_candidates = subject_var)
  repeated <- isTRUE(rm_det$repeated)

  notes <- character()
  if (repeated) {
    primary <- "ancombc2"
    reason <- sprintf(paste("Repeated measures on '%s' require a random effect.",
      "ANCOM-BC2 accepts rand_formula and is the only one of these four that",
      "models within-subject correlation directly. DESeq2, ALDEx2 and plain",
      "Wilcoxon all assume independent samples and would overstate significance."),
      rm_det$subject_var)
    alternatives <- c("glmmTMB negative-binomial GLMM per taxon (07-models.R)")
  } else if (!is.null(covariates) && length(covariates)) {
    primary <- "ancombc2"
    reason <- paste("Covariates must enter the model, not be ignored.",
      "ANCOM-BC2 and DESeq2 both accept a multi-term design; ALDEx2's glm mode",
      "does too but is slower. Plain Wilcoxon cannot adjust for anything.")
    alternatives <- c("deseq2", "aldex2 (glm mode)")
  } else if (n_min < 10) {
    primary <- "aldex2"
    reason <- sprintf(paste("Smallest group is n = %d. ALDEx2's Monte Carlo",
      "sampling propagates the uncertainty that small counts carry rather than",
      "assuming it away, and it reports effect size and dispersion so a",
      "borderline call can be judged. DESeq2 has more power here but its",
      "dispersion shrinkage is doing heavy lifting at this n."), n_min)
    alternatives <- c("deseq2 (more power, more assumption)", "ancombc2")
  } else if (zero_prop > 0.85) {
    primary <- "ancombc2"
    reason <- sprintf(paste("%.0f%% of the table is zeros. ANCOM-BC2 distinguishes",
      "structural zeros (the taxon is absent) from sampling zeros (present but",
      "unobserved), which is the distinction that matters at this sparsity and",
      "which the other three conflate."), 100 * zero_prop)
    alternatives <- c("aldex2", "zero-inflated models per taxon (07-models.R)")
  } else {
    primary <- "ancombc2"
    reason <- paste("Balanced design with adequate n and no repeated measures.",
      "ANCOM-BC2 is the default because bias correction addresses the",
      "compositional problem that makes raw counts misleading.")
    alternatives <- c("aldex2", "deseq2")
  }

  if (n_groups > 2) {
    notes <- c(notes, sprintf(paste("%d groups: ANCOM-BC2 gives a global test plus",
      "pairwise/Dunnett; DESeq2 needs an LRT against the reduced model for the",
      "overall effect, then Wald contrasts per pair; ALDEx2 uses kw/glm mode."),
      n_groups))
  }
  if (n_taxa > 500) {
    notes <- c(notes, sprintf(paste("%d taxa is a large multiple-testing burden.",
      "Agglomerating to Genus or Family, or tightening the prevalence filter,",
      "buys real power — but the choice must be made before testing."), n_taxa))
  }

  out <- list(primary = primary, reason = reason, alternatives = alternatives,
              notes = notes, n_groups = n_groups, n_per_group = n_per,
              zero_proportion = zero_prop, repeated_measures = repeated,
              subject_var = rm_det$subject_var, n_taxa = n_taxa)
  mfg_log("da", "method_chosen", list(primary = primary, n_groups = n_groups,
    smallest_group = n_min, zero_proportion = round(zero_prop, 3),
    repeated = repeated))
  class(out) <- c("mfg_da_choice", "list")
  out
}

#' @export
print.mfg_da_choice <- function(x, ...) {
  cat("=== Differential abundance method ===\n")
  cat(sprintf("Recommended: %s\n", x$primary))
  cat(strwrap(x$reason, width = 78, prefix = "  "), sep = "\n")
  cat(sprintf("\nAlternatives: %s\n", paste(x$alternatives, collapse = "; ")))
  cat(sprintf("Design: %d groups (%s), %d taxa, %.0f%% zeros%s\n",
    x$n_groups, paste(sprintf("%s=%d", names(x$n_per_group),
                              as.integer(x$n_per_group)), collapse = ", "),
    x$n_taxa, 100 * x$zero_proportion,
    if (x$repeated_measures) sprintf(", repeated on '%s'", x$subject_var) else ""))
  if (length(x$notes)) {
    cat("\nNotes:\n")
    for (n in x$notes) cat(strwrap(n, width = 76, prefix = "  - "), sep = "\n")
  }
  invisible(x)
}

# ── Wilcoxon / Kruskal-Wallis across taxa ────────────────────────────────────

#' Per-taxon rank test with FDR correction.
#'
#' The loop from Demo 7 and Demo 9, made explicit about what it does and does
#' not do. It has no model of compositionality, no covariate adjustment and no
#' dispersion estimate — which is exactly why it is the weakest of the four and
#' should be reported as a sanity check against a proper method, not on its own.
da_kruskal <- function(ps, group_var, p_adjust = "BH", alpha = 0.05,
                       min_prevalence = 0, check_norm = TRUE) {
  if (isTRUE(check_norm)) mfg_check_normalization(ps, "da_kruskal")

  meta <- mfg_meta(ps)
  g <- meta[[group_var]]
  if (is.null(g)) stop("'", group_var, "' is not in the metadata.", call. = FALSE)
  mat <- mfg_otu_samples_as_rows(ps)    # samples x taxa
  keep <- stats::complete.cases(g)
  mat <- mat[keep, , drop = FALSE]; g <- droplevels(as.factor(g[keep]))

  if (min_prevalence > 0) {
    prev <- colSums(mat > 0) / nrow(mat)
    mat  <- mat[, prev >= min_prevalence, drop = FALSE]
  }

  n_levels <- length(levels(g))
  test_name <- if (n_levels == 2) "wilcoxon" else "kruskal"

  rows <- lapply(colnames(mat), function(tx) {
    y <- mat[, tx]
    ft <- try(if (n_levels == 2) stats::wilcox.test(y ~ g, exact = FALSE)
              else stats::kruskal.test(y ~ g), silent = TRUE)
    if (inherits(ft, "try-error")) {
      return(data.frame(taxon = tx, statistic = NA_real_, p_value = NA_real_,
                        stringsAsFactors = FALSE))
    }
    data.frame(taxon = tx, statistic = unname(ft$statistic), p_value = ft$p.value,
               stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, rows)

  # Group medians and an effect size, because reference/10 requires one.
  med <- t(vapply(colnames(mat), function(tx)
    vapply(split(mat[, tx], g), stats::median, numeric(1)), numeric(n_levels)))
  colnames(med) <- paste0("median_", levels(g))
  res <- cbind(res, as.data.frame(med)[res$taxon, , drop = FALSE])
  if (n_levels == 2) {
    # A fold-change ratio is only meaningful on a non-negative scale. CLR and VST
    # values go negative, where log2 of a ratio is NaN, so the difference in
    # medians is reported instead of a ratio.
    if (all(mat >= 0, na.rm = TRUE)) {
      res$log2_median_ratio <- log2((med[res$taxon, 2] + 1e-9) / (med[res$taxon, 1] + 1e-9))
    } else {
      res$median_difference <- med[res$taxon, 2] - med[res$taxon, 1]
    }
    res$effect_size <- vapply(res$taxon, function(tx)
      mfg_rank_biserial(mat[, tx], g), numeric(1))
    res$effect_name <- "rank-biserial r"
  } else {
    N <- length(g); k <- n_levels
    res$effect_size <- pmax(0, (res$statistic - k + 1) / (N - k))
    res$effect_name <- "eta-squared (from H)"
  }
  res$prevalence <- colSums(mat > 0)[res$taxon] / nrow(mat)
  res$q_value <- stats::p.adjust(res$p_value, method = p_adjust)
  res <- res[order(res$q_value, res$p_value), ]
  res$significant <- !is.na(res$q_value) & res$q_value < alpha

  mfg_log("da", "kruskal", list(
    test = test_name, group = group_var, n_taxa_tested = nrow(res),
    p_adjust = p_adjust, n_significant_raw = sum(res$p_value < alpha, na.rm = TRUE),
    n_significant_adjusted = sum(res$significant),
    normalization = mfg_normalization(ps)))

  structure(list(method = test_name, group = group_var, results = res,
    n_tested = nrow(res), n_significant = sum(res$significant),
    p_adjust = p_adjust, alpha = alpha,
    caveat = paste("Rank tests on transformed abundances carry no model of",
      "compositionality, cannot adjust for covariates, and assume independent",
      "samples. Report alongside a compositional method, not instead of one.")),
    class = c("mfg_da_result", "list"))
}

# ── ALDEx2 ───────────────────────────────────────────────────────────────────

#' ALDEx2 differential abundance.
#'
#' Monte Carlo Dirichlet-multinomial sampling propagates the counting uncertainty
#' in sparse data, then CLR-transforms each instance, then tests. Because it
#' models sampling variation explicitly it is the method that behaves best at
#' small n.
#'
#' mc.samples: the authors recommend at least 128, and 1000 for a rigorous effect
#' size. 128 is the default here because it is what the course material uses;
#' raise it before reporting effect sizes as findings.
#'
#' Which p-value column to trust: wi.eBH (Wilcoxon, no normality assumption) for
#' two groups, glm.eBH for more than two. we.eBH assumes normality of CLR values.
#'
#' Covariates. Given `fix_formula`, the model-matrix path runs instead: the CLR
#' instances are fitted with aldex.glm() against model.matrix(fix_formula), so
#' the exposure is adjusted rather than compared marginally. Note that
#' aldex.glm() defaults to Holm; `p_adj_method` here defaults to BH and is
#' passed explicitly, because the default is a silent difference in what gets
#' called significant.
#'
#' ALDEx2 needs integer counts in either mode, because it Monte Carlo samples
#' from a Dirichlet-multinomial. Non-integer input is refused rather than
#' rounded: rounding a length-normalized rate such as HUMAnN RPK manufactures
#' counts that were never observed, and the whole point of the method is to
#' propagate genuine counting uncertainty.
#'
#' Expect the glm path to be much the more conservative of the two, and do not
#' read that as a failure. On a 300-feature simulation at n = 168 with twenty
#' planted four-fold signals, the marginal Wilcoxon path recovered all twenty at
#' a median p of 2e-08 while the glm path recovered none, at a median p of 0.15.
#' Dropping the covariates from the model matrix changed nothing, so the cost is
#' the per-instance GLM aggregation rather than adjustment, and raising
#' mc.samples from 64 to 256 barely moved it. Both paths held their nominal
#' false positive rate. When the two arms of a specification grid disagree this
#' is usually why, and it is a property of the method worth reporting rather
#' than a discrepancy worth hiding.
da_aldex2 <- function(ps, group_var, mc.samples = 128, test = NULL,
                      denom = "all", alpha = 0.1, effect = TRUE,
                      fix_formula = NULL, p_adj_method = "BH",
                      check_norm = TRUE, seed = 42) {
  if (isTRUE(check_norm)) mfg_check_normalization(ps, "da_aldex2")

  meta <- mfg_meta(ps)
  g <- meta[[group_var]]
  if (is.null(g)) stop("'", group_var, "' is not in the metadata.", call. = FALSE)
  keep <- stats::complete.cases(g)
  counts <- mfg_otu_taxa_as_rows(ps)[, keep, drop = FALSE]   # taxa x samples
  g <- droplevels(as.factor(g[keep]))
  n_levels <- length(levels(g))

  # Rows that are zero everywhere carry no information and slow the MC sampling.
  counts <- counts[rowSums(counts) > 0, , drop = FALSE]

  # Refuse non-integer input rather than coercing it. storage.mode() would
  # truncate silently, which turns a rate into a fabricated count.
  if (any(abs(counts - round(counts)) > 1e-8, na.rm = TRUE)) {
    stop("ALDEx2 needs integer counts; this matrix holds non-integer values.\n",
         "  It samples from a Dirichlet-multinomial, so there is nothing valid ",
         "to sample from.\n  Rounding or rescaling would invent counts that were ",
         "never observed. Use a\n  method that accepts proportions, such as ",
         "da_linda(), instead.", call. = FALSE)
  }
  storage.mode(counts) <- "integer"

  # ── Covariate-adjusted path ──────────────────────────────────────────────
  if (!is.null(fix_formula)) {
    fml <- if (is.character(fix_formula) && length(fix_formula) == 1) {
      stats::as.formula(if (grepl("^\\s*~", fix_formula)) fix_formula else paste("~", fix_formula))
    } else fix_formula
    mvars <- all.vars(fml)
    absent <- setdiff(mvars, names(meta))
    if (length(absent)) {
      stop("Not in the metadata: ", paste(absent, collapse = ", "), call. = FALSE)
    }
    md <- meta[keep, mvars, drop = FALSE]
    ok <- stats::complete.cases(md)
    if (any(!ok)) {
      warning(sprintf("%d sample(s) dropped for missing model variables.", sum(!ok)),
              call. = FALSE)
      md <- md[ok, , drop = FALSE]
      counts <- counts[, ok, drop = FALSE]
    }
    mm <- stats::model.matrix(fml, data = md)
    if (nrow(mm) != ncol(counts)) {
      stop("Model matrix has ", nrow(mm), " rows but ", ncol(counts),
           " samples remain. A covariate level was probably dropped.", call. = FALSE)
    }

    set.seed(seed)
    clr_obj <- ALDEx2::aldex.clr(counts, mm, mc.samples = mc.samples,
                                 denom = denom, verbose = FALSE)
    gl <- ALDEx2::aldex.glm(clr_obj, fdr.method = p_adj_method)

    # Pick the coefficient for the exposure. model.matrix() names a factor
    # column <var><level>, so match on prefix and say what was used.
    coefs <- unique(sub("[:.].*$", "", names(gl)))
    hit <- grep(paste0("^", group_var), names(gl), value = TRUE)
    if (!length(hit)) {
      stop("No aldex.glm coefficient for '", group_var, "'. Coefficients: ",
           paste(setdiff(coefs, "Intercept"), collapse = ", "), call. = FALSE)
    }
    term <- sub("[:].*$", "", hit[1])

    pick <- function(suffix) {
      col <- paste0(term, ":", suffix)
      if (col %in% names(gl)) gl[[col]] else rep(NA_real_, nrow(gl))
    }
    res <- data.frame(
      taxon    = rownames(gl),
      estimate = pick("Est"), std_error = pick("SE"), t_value = pick("t.val"),
      p_value  = pick("pval"), q_value = pick("pval.padj"),
      stringsAsFactors = FALSE)

    # aldex.glm.effect() returns a standardised effect only for binary terms.
    eff_col <- NA_character_
    if (isTRUE(effect)) {
      ge <- tryCatch(ALDEx2::aldex.glm.effect(clr_obj, verbose = FALSE),
                     error = function(e) NULL)
      if (!is.null(ge) && term %in% names(ge)) {
        e_df <- ge[[term]]
        res$effect      <- e_df[match(res$taxon, rownames(e_df)), "effect"]
        res$diff_btw    <- e_df[match(res$taxon, rownames(e_df)), "diff.btw"]
        res$overlap     <- e_df[match(res$taxon, rownames(e_df)), "overlap"]
        eff_col <- "effect"
      }
    }
    res$effect_size <- if (!is.na(eff_col)) res[[eff_col]] else res$estimate
    res$significant <- !is.na(res$q_value) & res$q_value < alpha
    res <- res[order(res$q_value, res$p_value), , drop = FALSE]

    mfg_log("da", "aldex2_glm", list(
      group = group_var, term = term, fix_formula = paste(deparse(fml), collapse = ""),
      mc.samples = mc.samples, denom = denom, p_adj_method = p_adj_method,
      alpha = alpha, n_samples = ncol(counts),
      n_taxa_tested = nrow(res), n_significant = sum(res$significant, na.rm = TRUE),
      seed = seed))

    return(structure(list(
      method = "aldex2", mode = "glm", group = group_var, term = term,
      fix_formula = paste(deparse(fml), collapse = ""),
      clr = clr_obj, results = res,
      n_tested = nrow(res), n_significant = sum(res$significant, na.rm = TRUE),
      q_column = paste0(term, ":pval.padj"), effect_column = eff_col,
      alpha = alpha, mc.samples = mc.samples, p_adj_method = p_adj_method,
      aldex_settings = sprintf(
        "glm mode, term = %s, mc.samples = %s, denom = %s, p_adj = %s",
        term, mc.samples, denom, p_adj_method),
      caveat = if (mc.samples < 1000 && isTRUE(effect)) paste(
        "mc.samples =", mc.samples, "- the ALDEx2 authors recommend 1000 for",
        "rigorous effect size estimates. Effect sizes here are indicative.") else NA_character_),
      class = c("mfg_da_result", "list")))
  }

  test <- test %||% if (n_levels == 2) "t" else "kw"
  if (n_levels > 2 && test %in% c("t")) {
    stop("ALDEx2 test = 't' compares two groups but '", group_var, "' has ",
         n_levels, " levels. Use test = 'kw' (Kruskal-Wallis and glm).",
         call. = FALSE)
  }

  set.seed(seed)
  res <- ALDEx2::aldex(counts, as.character(g), mc.samples = mc.samples,
                       test = test, effect = effect, denom = denom,
                       include.sample.summary = FALSE, verbose = FALSE)
  res <- as.data.frame(res)
  res$taxon <- rownames(res)

  q_col <- if (n_levels == 2) "wi.eBH" else "glm.eBH"
  p_col <- if (n_levels == 2) "wi.ep"  else "glm.ep"
  if (!q_col %in% names(res)) { q_col <- "we.eBH"; p_col <- "we.ep" }

  res$p_value <- res[[p_col]]
  res$q_value <- res[[q_col]]
  res$significant <- !is.na(res$q_value) & res$q_value < alpha
  res <- res[order(res$q_value, res$p_value), ]

  # ALDEx2 reports an effect size only for two-group comparisons. In kw/glm mode
  # (>2 groups) there is no single effect direction to report, so the column is
  # absent and downstream code must not assume it.
  eff_col <- intersect(c("effect", "diff.btw"), names(res))
  eff_col <- if (length(eff_col)) eff_col[1] else NA_character_

  mfg_log("da", "aldex2", list(
    group = group_var, test = test, mc.samples = mc.samples, denom = denom,
    q_column = q_col, alpha = alpha, n_taxa_tested = nrow(res),
    n_significant = sum(res$significant), seed = seed))

  structure(list(method = "aldex2", group = group_var, results = res,
    n_tested = nrow(res), n_significant = sum(res$significant),
    q_column = q_col, effect_column = eff_col,
    alpha = alpha, mc.samples = mc.samples, test = test, n_groups = n_levels,
    caveat = if (mc.samples < 1000 && effect) paste(
      "mc.samples =", mc.samples, "— the ALDEx2 authors recommend 1000 for",
      "rigorous effect size estimates. Effect sizes here are indicative.") else NA_character_),
    class = c("mfg_da_result", "list"))
}

# ── ANCOM-BC2 ────────────────────────────────────────────────────────────────

#' ANCOM-BC2 differential abundance.
#'
#' Estimates the true log absolute abundance of each taxon while correcting the
#' sampling-fraction bias that makes relative abundances misleading. Tests log
#' fold change with a Wald test.
#'
#' Two parameters do most of the work and must be reported:
#'   prv_cut   minimum prevalence. Taxa below it are dropped before testing, so
#'             this changes the multiple-testing denominator.
#'   struc_zero  when TRUE, a taxon absent from an entire group is declared a
#'             structural zero and handled separately rather than being tested as
#'             a fold change against zero. Turn this on for sparse data; see
#'             reference/11 on what it licenses you to claim.
da_ancombc2 <- function(ps, fix_formula, group = NULL, rand_formula = NULL,
                        tax_level = NULL, p_adj_method = "holm",
                        prv_cut = 0.10, lib_cut = 1000,
                        struc_zero = TRUE, neg_lb = FALSE,
                        alpha = 0.05, global = NULL, pairwise = FALSE,
                        dunnet = FALSE, check_norm = TRUE, verbose = FALSE) {
  if (isTRUE(check_norm)) mfg_check_normalization(ps, "da_ancombc2")

  if (is.character(fix_formula) && length(fix_formula) == 1) {
    fix_chr <- fix_formula
  } else {
    fix_chr <- paste(all.vars(fix_formula), collapse = " + ")
  }
  meta <- mfg_meta(ps)
  n_levels <- if (!is.null(group)) length(unique(meta[[group]][!is.na(meta[[group]])])) else NA
  global <- global %||% (!is.na(n_levels) && n_levels > 2)

  out <- ANCOMBC::ancombc2(
    data = ps, fix_formula = fix_chr, rand_formula = rand_formula,
    tax_level = tax_level, p_adj_method = p_adj_method,
    prv_cut = prv_cut, lib_cut = lib_cut, group = group,
    struc_zero = struc_zero, neg_lb = neg_lb,
    iter_control = list(tol = 1e-2, max_iter = 20, verbose = verbose),
    em_control = list(tol = 1e-5, max_iter = 100),
    alpha = alpha, global = global, pairwise = pairwise, dunnet = dunnet,
    verbose = verbose)

  primary <- out$res
  # Primary results carry one lfc_/q_ column set per model term. The intercept is
  # a baseline, not a comparison — a "significant" intercept only says the taxon's
  # abundance differs from zero in the reference group, so it is excluded from the
  # finding counts. diff_robust_* are the sensitivity-analysis-passing subset and
  # are reported separately rather than added to the same total.
  drop_intercept <- function(v) v[!grepl("\\(Intercept\\)", v)]
  q_cols    <- drop_intercept(grep("^q_", names(primary), value = TRUE))
  diff_cols <- drop_intercept(grep("^diff_", names(primary), value = TRUE))
  diff_main <- diff_cols[!grepl("^diff_robust_", diff_cols)]
  diff_rob  <- diff_cols[grepl("^diff_robust_", diff_cols)]
  n_sig     <- vapply(diff_main, function(c) sum(primary[[c]], na.rm = TRUE), integer(1))
  n_sig_robust <- vapply(diff_rob, function(c) sum(primary[[c]], na.rm = TRUE), integer(1))

  n_struc_zero <- if (!is.null(out$zero_ind)) {
    zi <- out$zero_ind
    sum(apply(zi[, -1, drop = FALSE], 1, any), na.rm = TRUE)
  } else NA_integer_

  mfg_log("da", "ancombc2", list(
    fix_formula = fix_chr, rand_formula = rand_formula %||% "none",
    tax_level = tax_level %||% "ASV/feature", p_adj_method = p_adj_method,
    prv_cut = prv_cut, lib_cut = lib_cut, struc_zero = struc_zero,
    global = global, alpha = alpha,
    n_taxa_tested = nrow(primary),
    n_significant = paste(sprintf("%s=%d", names(n_sig), n_sig), collapse = " "),
    n_significant_robust = paste(sprintf("%s=%d", names(n_sig_robust), n_sig_robust),
                                 collapse = " "),
    n_structural_zeros = n_struc_zero))

  structure(list(method = "ancombc2", group = group, fix_formula = fix_chr,
    rand_formula = rand_formula, out = out, results = primary,
    res_global = out$res_global, res_pair = out$res_pair, res_dunn = out$res_dunn,
    zero_ind = out$zero_ind, n_structural_zeros = n_struc_zero,
    n_tested = nrow(primary), n_significant = n_sig,
    n_significant_robust = n_sig_robust,
    q_columns = q_cols, diff_columns = diff_main,
    prv_cut = prv_cut, lib_cut = lib_cut, struc_zero = struc_zero,
    p_adj_method = p_adj_method, alpha = alpha),
    class = c("mfg_da_result", "list"))
}

# ── LinDA ────────────────────────────────────────────────────────────────────

#' LinDA differential abundance.
#'
#' Fits a linear model to CLR-transformed abundances and then corrects the
#' compositional bias in the coefficients, rather than trying to remove it
#' before fitting. Two consequences matter in practice. It takes a full model
#' formula, so covariates are handled natively. And it accepts proportions as
#' well as counts, which makes it the usable option when an assay cannot be
#' resolved to integers, where ALDEx2 and DESeq2 cannot run at all.
#'
#' feature_dat_type: "count" or "proportion". Left NULL it is inferred, but the
#' inference is a guess about the data and is logged as one. State it when you
#' know, because LinDA's zero handling differs between the two.
#'
#' Pseudo-count sensitivity. LinDA has no equivalent of ANCOM-BC2's built-in
#' sensitivity flag, yet a result that survives only one arbitrary zero
#' replacement is not a finding. With sensitivity = TRUE the model is refitted
#' across `sens_pseudo_cnt` with adaptive zero handling switched off, and a
#' feature is called robust only if it stays significant with the same sign in
#' every fit. The count of robust features is returned as `diff_robust_<term>`,
#' matching the ANCOM-BC2 field so that downstream code and printing treat the
#' two the same way.
da_linda <- function(ps, fix_formula, group = NULL,
                     feature_dat_type = NULL,
                     prv_cut = 0.10,
                     zero_handling = c("pseudo-count", "imputation"),
                     pseudo_cnt = 0.5, adaptive = TRUE,
                     is_winsor = TRUE, outlier_pct = 0.03, corr_cut = 0.1,
                     p_adj_method = "BH", alpha = 0.05,
                     sensitivity = TRUE, sens_pseudo_cnt = c(0.1, 0.5, 1),
                     n_cores = 1, check_norm = TRUE, verbose = FALSE) {

  if (!requireNamespace("MicrobiomeStat", quietly = TRUE)) {
    stop("LinDA needs the MicrobiomeStat package.\n",
         '  install.packages("MicrobiomeStat")', call. = FALSE)
  }
  if (isTRUE(check_norm)) mfg_check_normalization(ps, "da_linda")
  zero_handling <- match.arg(zero_handling)

  fix_chr <- if (is.character(fix_formula) && length(fix_formula) == 1) {
    sub("^\\s*~\\s*", "", fix_formula)
  } else {
    paste(all.vars(fix_formula), collapse = " + ")
  }
  model_vars <- all.vars(stats::as.formula(paste("~", fix_chr)))

  meta <- mfg_meta(ps)
  miss <- setdiff(model_vars, names(meta))
  if (length(miss)) {
    stop("Not in the metadata: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  # Complete cases only, and say how many were dropped. LinDA would otherwise
  # fail opaquely on an NA in a covariate.
  keep <- stats::complete.cases(meta[, model_vars, drop = FALSE])
  n_dropped <- sum(!keep)
  if (n_dropped) {
    warning(sprintf("%d sample(s) dropped for missing model variables.", n_dropped),
            call. = FALSE)
  }
  feat <- mfg_otu_taxa_as_rows(ps)[, keep, drop = FALSE]   # taxa x samples
  meta <- meta[keep, , drop = FALSE]

  feat <- feat[rowSums(feat) > 0, , drop = FALSE]

  if (is.null(feature_dat_type)) {
    feature_dat_type <- if (looks_like_relative_abundance(feat)) "proportion" else "count"
    inferred <- TRUE
    warning(sprintf(paste0(
      "feature.dat.type was inferred as '%s'.\n",
      "  This is not cosmetic: the count and proportion paths handle zeros ",
      "differently and\n  pseudo.cnt applies only to counts, so the same matrix ",
      "declared either way can give\n  materially different p-values. Pass ",
      "feature_dat_type explicitly."), feature_dat_type), call. = FALSE)
  } else {
    feature_dat_type <- match.arg(feature_dat_type, c("count", "proportion"))
    inferred <- FALSE
  }

  fit_once <- function(pc, adapt, zh) {
    MicrobiomeStat::linda(
      feature.dat = feat, meta.dat = meta,
      formula = paste("~", fix_chr),
      feature.dat.type = feature_dat_type,
      prev.filter = prv_cut, is.winsor = is_winsor, outlier.pct = outlier_pct,
      adaptive = adapt, zero.handling = zh, pseudo.cnt = pc,
      corr.cut = corr_cut, p.adj.method = p_adj_method, alpha = alpha,
      n.cores = n_cores, verbose = verbose)
  }

  out <- fit_once(pseudo_cnt, adaptive, zero_handling)

  # LinDA names each output element after a model coefficient. Pick the term
  # for `group` when given; otherwise take the first non-intercept term and say
  # which one was used rather than letting the choice go unrecorded.
  terms_all <- names(out$output)
  term <- if (!is.null(group)) {
    hit <- terms_all[startsWith(terms_all, group)]
    if (!length(hit)) {
      stop("No LinDA coefficient for group '", group, "'. Available: ",
           paste(terms_all, collapse = ", "), call. = FALSE)
    }
    hit[1]
  } else {
    terms_all[!grepl("\\(Intercept\\)", terms_all)][1]
  }

  res <- out$output[[term]]
  res$taxon <- rownames(res)
  # Column names the rest of the toolkit expects, so da_significant() and
  # da_compare() need no LinDA-specific branch.
  res$effect_size <- res$log2FoldChange
  res$q_value     <- res$padj
  res$p_value     <- res$pvalue
  res$significant <- !is.na(res$padj) & res$padj < alpha
  res <- res[order(res$q_value, res$p_value), , drop = FALSE]

  n_sig <- stats::setNames(sum(res$significant, na.rm = TRUE), term)

  # ── Pseudo-count sensitivity ───────────────────────────────────────────────
  robust_taxa <- NULL
  n_sig_robust <- NULL
  sens_detail  <- NULL
  if (isTRUE(sensitivity) && any(res$significant, na.rm = TRUE)) {
    grid <- sort(unique(c(pseudo_cnt, sens_pseudo_cnt)))
    per_fit <- lapply(grid, function(pc) {
      o <- tryCatch(fit_once(pc, FALSE, "pseudo-count"), error = function(e) NULL)
      if (is.null(o) || is.null(o$output[[term]])) return(NULL)
      d <- o$output[[term]]
      data.frame(taxon = rownames(d), sig = !is.na(d$padj) & d$padj < alpha,
                 sign = sign(d$log2FoldChange), stringsAsFactors = FALSE)
    })
    ok <- !vapply(per_fit, is.null, logical(1))
    if (any(ok)) {
      per_fit <- per_fit[ok]
      cand <- res$taxon[res$significant]
      robust_taxa <- Filter(function(tx) {
        rows <- lapply(per_fit, function(d) d[match(tx, d$taxon), , drop = FALSE])
        all(vapply(rows, function(r) isTRUE(r$sig), logical(1))) &&
          length(unique(vapply(rows, function(r) r$sign, numeric(1)))) == 1
      }, cand)
      res$passed_ss <- res$taxon %in% robust_taxa
      # Named to match ANCOM-BC2 so print.mfg_da_result and any downstream
      # reporting treat the two methods identically.
      n_sig_robust <- stats::setNames(length(robust_taxa),
                                      paste0("diff_robust_", term))
      sens_detail <- list(pseudo_counts = grid, n_fits = length(per_fit),
                          n_failed = sum(!ok))
    }
  }

  mfg_log("da", "linda", list(
    fix_formula = fix_chr, group = group %||% "none", term = term,
    feature_dat_type = feature_dat_type,
    feature_dat_type_inferred = inferred,
    prev_filter = prv_cut, zero_handling = zero_handling,
    pseudo_cnt = pseudo_cnt, adaptive = adaptive,
    is_winsor = is_winsor, outlier_pct = outlier_pct,
    p_adj_method = p_adj_method, alpha = alpha,
    n_samples = ncol(feat), n_samples_dropped = n_dropped,
    n_taxa_tested = nrow(res), n_significant = unname(n_sig),
    n_significant_robust = if (is.null(n_sig_robust)) NA_integer_ else unname(n_sig_robust),
    sensitivity_pseudo_counts = if (is.null(sens_detail)) "not run" else
      paste(sens_detail$pseudo_counts, collapse = ", ")))

  structure(list(
    method = "linda", group = group, fix_formula = fix_chr, term = term,
    out = out, results = res, all_terms = out$output,
    feature_dat_type = feature_dat_type,
    feature_dat_type_inferred = inferred,
    bias = out$bias,
    prv_cut_linda = prv_cut, zero_handling = zero_handling,
    pseudo_cnt = pseudo_cnt, adaptive = adaptive,
    p_adj_method = p_adj_method, alpha = alpha,
    n_tested = nrow(res), n_significant = n_sig,
    n_significant_robust = n_sig_robust,
    robust_taxa = robust_taxa, sensitivity = sens_detail,
    n_samples_dropped = n_dropped,
    linda_settings = sprintf(
      "input = %s%s, prev.filter = %s, zero handling = %s (pseudo.cnt %s, adaptive %s), p_adj = %s",
      feature_dat_type, if (inferred) " (inferred)" else "", prv_cut,
      zero_handling, pseudo_cnt, adaptive, p_adj_method),
    caveat = if (isTRUE(inferred)) paste(
      "feature.dat.type was inferred as", feature_dat_type, "rather than declared.",
      "The count and proportion paths handle zeros differently and pseudo.cnt",
      "applies only to counts, so the same matrix declared either way can give",
      "materially different p-values. Declare the type and re-run before",
      "reporting anything from this fit.") else NA_character_),
    class = c("mfg_da_result", "list"))
}

# ── DESeq2 ───────────────────────────────────────────────────────────────────

#' DESeq2 differential abundance.
#'
#' Negative binomial GLM with dispersion shrunk toward a fitted trend, which is
#' where its power at small n comes from. Requires raw counts: it estimates size
#' factors itself, and pre-normalized input breaks the mean-variance
#' relationship the shrinkage depends on.
#'
#' fitType = "local" is the default because microbiome dispersion-mean
#' relationships are often not well described by the parametric fit, and DESeq2
#' itself falls back to local in that case (Demo 10 uses local explicitly).
#'
#' alpha: DESeq2's default independent-filtering threshold is 0.1, not 0.05.
#' That matters — `results()` optimises filtering for the alpha you pass, so
#' passing 0.1 and then calling 0.05 significant is not the same as passing 0.05.
da_deseq2 <- function(ps, design, test = c("Wald", "LRT"), reduced = ~ 1,
                      fitType = "local", alpha = 0.1, contrast = NULL,
                      cooksCutoff = TRUE, shrink_lfc = TRUE,
                      check_norm = TRUE) {
  test <- match.arg(test)
  if (isTRUE(check_norm)) mfg_check_normalization(ps, "da_deseq2")

  if (is.character(design)) design <- stats::as.formula(design)
  dds <- phyloseq::phyloseq_to_deseq2(ps, design)

  # Microbiome data almost never has a taxon present in every sample, so the
  # default geometric mean is zero and size-factor estimation fails.
  #
  # When it does NOT fail it can still be degenerate, and that case is worth
  # reporting: DESeq2's default computes the reference geometric mean only over
  # rows with no zeros, so if a handful of taxa qualify, every sample's size
  # factor is estimated from those few. On a sparse table that is one taxon's
  # ratio masquerading as a median-of-ratios estimate, and it produces size
  # factors uncorrelated with the all-taxa estimate. The positive-count
  # geometric mean below uses every taxon and is the documented phyloseq
  # workaround (Demo 10, ps10 Q1).
  cm <- DESeq2::counts(dds)
  n_complete <- sum(rowSums(cm == 0) == 0)
  geo <- apply(cm, 1, gm_mean)
  dds <- DESeq2::estimateSizeFactors(dds, geoMeans = geo)
  sf <- DESeq2::sizeFactors(dds)

  if (n_complete < max(10, 0.05 * nrow(cm))) {
    message(sprintf(paste("Size factors: only %d of %d taxa are present in every",
      "sample, so DESeq2's default estimator would have had almost nothing to work",
      "from. Using the positive-count geometric mean over all taxa instead (median",
      "size factor %.3f). If you are comparing against an analysis that used the",
      "default, expect different results — and the default is the less reliable of",
      "the two here."), n_complete, nrow(cm), stats::median(sf)))
  }
  mfg_log("da", "deseq2_size_factors", list(
    n_taxa_complete = n_complete, n_taxa = nrow(cm),
    method = "positive-count geometric mean over all taxa",
    median_size_factor = round(stats::median(sf), 4),
    size_factor_range = sprintf("%.3f-%.3f", min(sf), max(sf)),
    default_estimator_degenerate = n_complete < max(10, 0.05 * nrow(cm))))

  dds <- if (identical(test, "LRT")) {
    DESeq2::DESeq(dds, test = "LRT", reduced = reduced, fitType = fitType, quiet = TRUE)
  } else {
    DESeq2::DESeq(dds, test = "Wald", fitType = fitType, quiet = TRUE)
  }

  res <- if (is.null(contrast)) {
    DESeq2::results(dds, alpha = alpha, cooksCutoff = cooksCutoff)
  } else {
    DESeq2::results(dds, contrast = contrast, alpha = alpha, cooksCutoff = cooksCutoff)
  }

  # Shrunken log fold changes are the ones to report: unshrunken LFCs for
  # low-count taxa are wildly inflated and make volcano plots misleading.
  res_shrunk <- NULL
  if (isTRUE(shrink_lfc) && identical(test, "Wald")) {
    res_shrunk <- try(DESeq2::lfcShrink(dds, res = res,
      coef = DESeq2::resultsNames(dds)[length(DESeq2::resultsNames(dds))],
      type = "apeglm", quiet = TRUE), silent = TRUE)
    if (inherits(res_shrunk, "try-error")) {
      res_shrunk <- try(DESeq2::lfcShrink(dds, res = res,
        coef = DESeq2::resultsNames(dds)[length(DESeq2::resultsNames(dds))],
        type = "normal", quiet = TRUE), silent = TRUE)
    }
    if (inherits(res_shrunk, "try-error")) res_shrunk <- NULL
  }

  df <- as.data.frame(res)
  df$taxon <- rownames(df)
  if (!is.null(res_shrunk)) {
    df$log2FoldChange_shrunk <- as.data.frame(res_shrunk)[df$taxon, "log2FoldChange"]
  }
  # Attach taxonomy so the result table is readable without a join.
  tt <- mfg_tax(ps)
  if (!is.null(tt)) {
    for (r in intersect(c("Phylum", "Family", "Genus", "Species"), names(tt))) {
      df[[r]] <- tt[df$taxon, r]
    }
  }
  df$significant <- !is.na(df$padj) & df$padj < alpha
  df <- df[order(df$padj, df$pvalue), ]

  # DESeq2's independent filtering sets padj to NA for low-count taxa, which is
  # not a missing value — it is a deliberate exclusion that increases power on
  # the rest. It must be reported, because "23 of 831 significant" means
  # something different from "23 of 83 tested".
  n_filtered <- sum(is.na(df$padj) & !is.na(df$pvalue))
  n_tested   <- sum(!is.na(df$padj))

  mfg_log("da", "deseq2", list(
    design = paste(deparse(design), collapse = ""), test = test,
    fitType = fitType, alpha = alpha,
    contrast = if (is.null(contrast)) "default" else paste(contrast, collapse = "/"),
    n_taxa_total = nrow(df), n_tested_after_filtering = n_tested,
    n_independent_filtered = n_filtered,
    n_significant = sum(df$significant),
    lfc_shrunk = !is.null(res_shrunk)))

  structure(list(method = "deseq2", dds = dds, res = res, res_shrunk = res_shrunk,
    results = df, design = design, test = test, fitType = fitType, alpha = alpha,
    contrast = contrast, n_tested = n_tested, n_total = nrow(df),
    n_taxa_complete = n_complete, size_factors = sf,
    n_independent_filtered = n_filtered, n_significant = sum(df$significant),
    caveat = sprintf(paste("%d taxa had padj set to NA by DESeq2's independent",
      "filtering (low mean count). The denominator for the FDR is %d, not %d —",
      "state both in the report."), n_filtered, n_tested, nrow(df))),
    class = c("mfg_da_result", "list"))
}

#' DESeq2 diagnostic plots.
#'
#' Not decoration. plotDispEsts says whether the dispersion model fits — if the
#' black points do not scatter around the red trend, the negative binomial
#' assumption is in trouble and the p-values are not trustworthy. The p-value
#' histogram should be flat with a spike near zero; other shapes mean something
#' is wrong with the test, not with the biology.
da_deseq2_diagnostics <- function(fit, save = TRUE, name_prefix = "deseq2") {
  stopifnot(inherits(fit, "mfg_da_result"), identical(fit$method, "deseq2"))
  paths <- character()

  if (isTRUE(save)) {
    p <- file.path(mfg_run_dir("figures"), paste0(name_prefix, "_MA.png"))
    grDevices::png(p, width = 7, height = 5, units = "in", res = 300)
    DESeq2::plotMA(fit$res, main = "MA plot: mean abundance vs log2 fold change")
    grDevices::dev.off(); paths <- c(paths, p)

    p <- file.path(mfg_run_dir("figures"), paste0(name_prefix, "_dispersion.png"))
    grDevices::png(p, width = 7, height = 5, units = "in", res = 300)
    DESeq2::plotDispEsts(fit$dds, main = "Dispersion estimates")
    grDevices::dev.off(); paths <- c(paths, p)

    p <- file.path(mfg_run_dir("figures"), paste0(name_prefix, "_pvalue_hist.png"))
    grDevices::png(p, width = 7, height = 5, units = "in", res = 300)
    graphics::hist(fit$results$pvalue, breaks = 50, col = "grey80",
                   main = "p-value distribution", xlab = "raw p-value")
    grDevices::dev.off(); paths <- c(paths, p)
  }

  # Quantify what the eye is supposed to check.
  pv <- fit$results$pvalue[!is.na(fit$results$pvalue)]
  # Under the null, p-values are uniform; the upper half should hold ~50%.
  upper_half_frac <- mean(pv > 0.5)
  hist_verdict <- if (upper_half_frac > 0.40 && upper_half_frac < 0.60)
      "p-value distribution looks well behaved (roughly uniform above 0.5)"
    else if (upper_half_frac >= 0.60)
      "p-values are conservative / over-dispersed — more large p-values than expected under the null; the test may be losing power"
    else
      "p-values are anti-conservative — too few large p-values; check for unmodelled structure (batch, repeated measures) inflating significance"

  disp <- DESeq2::dispersions(fit$dds)
  mfg_log("da", "deseq2_diagnostics", list(
    figures = paste(basename(paths), collapse = ","),
    pvalue_upper_half_fraction = round(upper_half_frac, 3),
    pvalue_verdict = hist_verdict,
    median_dispersion = round(stats::median(disp, na.rm = TRUE), 4)))

  list(figures = paths, pvalue_upper_half_fraction = upper_half_frac,
       pvalue_verdict = hist_verdict,
       median_dispersion = stats::median(disp, na.rm = TRUE))
}

# ── Printing and cross-method comparison ─────────────────────────────────────

#' @export
print.mfg_da_result <- function(x, ...) {
  cat(sprintf("=== Differential abundance: %s ===\n", x$method))
  if (!is.null(x$group))        cat(sprintf("Group: %s\n", x$group))
  if (!is.null(x$fix_formula))  cat(sprintf("Fixed: %s\n", x$fix_formula))
  if (!is.null(x$rand_formula)) cat(sprintf("Random: %s\n", x$rand_formula))
  if (!is.null(x$design))       cat(sprintf("Design: %s (%s test, fitType %s)\n",
                                            paste(deparse(x$design), collapse = ""),
                                            x$test, x$fitType))
  if (!is.null(x$prv_cut))      cat(sprintf("prv_cut = %s, lib_cut = %s, struc_zero = %s, p_adj = %s\n",
                                            x$prv_cut, x$lib_cut, x$struc_zero, x$p_adj_method))
  if (!is.null(x$linda_settings)) cat(x$linda_settings, "\n", sep = "")
  if (!is.null(x$aldex_settings)) cat(x$aldex_settings, "\n", sep = "")
  if (!is.null(x$n_taxa_complete)) {
    cat(sprintf("Size factors: positive-count geometric mean over all taxa (%d of %d taxa present in every sample; median SF %.3f)\n",
                x$n_taxa_complete, x$n_total, stats::median(x$size_factors)))
  }
  if (!is.null(x$q_column))     cat(sprintf("q-value column used: %s (mc.samples = %s)\n",
                                            x$q_column, x$mc.samples))
  if (!is.null(x$effect_column) && is.na(x$effect_column)) {
    cat(sprintf("No effect-size column: ALDEx2 reports one only for two-group comparisons (%d groups here).\n",
                x$n_groups %||% NA))
  }
  cat(sprintf("Tested %s taxa at alpha = %s\n", x$n_tested, x$alpha))
  if (!is.null(x$n_independent_filtered) && x$n_independent_filtered > 0) {
    cat(sprintf("  (%d of %d excluded by independent filtering)\n",
                x$n_independent_filtered, x$n_total))
  }
  if (length(x$n_significant) == 1) {
    cat(sprintf("Significant: %d\n", x$n_significant))
  } else {
    cat(sprintf("Significant: %s\n",
                paste(sprintf("%s=%d", names(x$n_significant), x$n_significant),
                      collapse = ", ")))
  }
  if (!is.null(x$n_significant_robust) && length(x$n_significant_robust)) {
    cat(sprintf("Passing sensitivity analysis (diff_robust): %s\n",
                paste(sprintf("%s=%d", sub("^diff_robust_", "",
                                           names(x$n_significant_robust)),
                              x$n_significant_robust), collapse = ", ")))
  }
  if (!is.null(x$n_structural_zeros) && !is.na(x$n_structural_zeros)) {
    cat(sprintf("Structural zeros detected: %d taxa\n", x$n_structural_zeros))
  }
  if (!is.null(x$caveat) && !is.na(x$caveat)) {
    cat("\n"); cat(strwrap(x$caveat, width = 78, prefix = "! "), sep = "\n")
  }
  invisible(x)
}

#' Significant taxa from any DA result, as one common shape.
da_significant <- function(x, n = Inf) {
  res <- x$results
  if (identical(x$method, "ancombc2")) {
    dcol <- x$diff_columns[1]
    if (length(dcol) == 0 || is.na(dcol) || is.null(dcol)) return(data.frame())
    sig <- res[which(res[[dcol]]), , drop = FALSE]
    term <- sub("^diff_", "", dcol)
    lfc <- grep(paste0("^lfc_", term, "$"), names(res), value = TRUE)
    if (!length(lfc)) lfc <- grep("^lfc_(?!\\(Intercept\\))", names(res),
                                  value = TRUE, perl = TRUE)[1]
    q   <- grep(paste0("^q_", term, "$"), names(res), value = TRUE)
    if (!length(q)) q <- x$q_columns[1]
    out <- data.frame(taxon = sig$taxon, effect = sig[[lfc]], q_value = sig[[q]],
                      stringsAsFactors = FALSE)
  } else if (identical(x$method, "deseq2")) {
    sig <- res[res$significant %in% TRUE, , drop = FALSE]
    eff <- if ("log2FoldChange_shrunk" %in% names(sig)) sig$log2FoldChange_shrunk
           else sig$log2FoldChange
    out <- data.frame(taxon = sig$taxon, effect = eff, q_value = sig$padj,
                      stringsAsFactors = FALSE)
  } else if (identical(x$method, "aldex2")) {
    sig <- res[res$significant %in% TRUE, , drop = FALSE]
    ec <- x$effect_column
    eff <- if (!is.na(ec) && ec %in% names(sig)) sig[[ec]] else
      rep(NA_real_, nrow(sig))   # kw/glm mode has no single effect direction
    out <- data.frame(taxon = sig$taxon, effect = eff,
                      q_value = sig$q_value, stringsAsFactors = FALSE)
  } else {
    sig <- res[res$significant %in% TRUE, , drop = FALSE]
    out <- data.frame(taxon = sig$taxon,
                      effect = sig$effect_size, q_value = sig$q_value,
                      stringsAsFactors = FALSE)
  }
  out <- out[order(out$q_value), , drop = FALSE]
  utils::head(out, n)
}

#' Compare hit lists across methods.
#'
#' Two methods agreeing on a taxon is a far stronger claim than one method
#' finding it. This reports the intersection, the per-method-only sets, and the
#' rank correlation of effect sizes over the union — which is the honest summary
#' of how much the method choice mattered.
da_compare <- function(...) {
  fits <- list(...)
  named <- !is.null(names(fits)) && all(nzchar(names(fits)))
  nm <- if (named) names(fits) else vapply(fits, function(f) f$method, character(1))
  names(fits) <- nm

  sig <- lapply(fits, function(f) da_significant(f)$taxon)
  all_taxa <- Reduce(union, sig)
  overlap <- Reduce(intersect, sig)

  # Namespace check. Zero overlap between methods is a real possibility, but it
  # is far more often a naming mismatch than a disagreement: ANCOM-BC2 with
  # tax_level = "Genus" returns genus names, while tax_glom() keeps the
  # representative ASV id, so the two label spaces are disjoint and every
  # comparison silently returns zero. Reporting "consensus = 0 taxa" in that case
  # is a clean-looking wrong answer, so it is detected rather than passed on.
  tested <- lapply(fits, function(f) f$results$taxon)
  pairwise_shared <- if (length(fits) >= 2) {
    min(vapply(utils::combn(seq_along(fits), 2, simplify = FALSE),
               function(ij) length(intersect(tested[[ij[1]]], tested[[ij[2]]])),
               integer(1)))
  } else Inf
  namespace_ok <- pairwise_shared > 0
  if (!namespace_ok) {
    warning(sprintf(paste("The methods share NO taxon labels at all, so every",
      "overlap here is zero for naming reasons rather than statistical ones.\n",
      "Label examples: %s\n",
      "This happens when one method was given a different taxonomic level or a",
      "differently-named object — e.g. ANCOM-BC2 with tax_level = \"Genus\"",
      "returns genus names while tax_glom() keeps the representative ASV id.\n",
      "Run both methods on the same object at the same rank, or map the labels",
      "through the taxonomy table, before comparing."),
      paste(sprintf("%s: %s", nm, vapply(tested, function(v)
        paste(utils::head(v, 2), collapse = ", "), character(1))), collapse = " | ")),
      call. = FALSE)
    mfg_log("da", "namespace_mismatch", list(
      methods = paste(nm, collapse = ","),
      shared_labels = 0,
      note = "comparison is meaningless until labels are reconciled"))
  }

  # Effect sizes over the union, for rank correlation.
  eff <- lapply(fits, function(f) {
    d <- da_significant(f, n = Inf)
    full <- switch(f$method,
      ancombc2 = { lfc <- grep("^lfc_", names(f$results), value = TRUE)[1]
                   stats::setNames(f$results[[lfc]], f$results$taxon) },
      deseq2   = { e <- if ("log2FoldChange_shrunk" %in% names(f$results))
                     f$results$log2FoldChange_shrunk else f$results$log2FoldChange
                   stats::setNames(e, f$results$taxon) },
      aldex2   = { ec <- f$effect_column
                   if (!is.na(ec) && ec %in% names(f$results))
                     stats::setNames(f$results[[ec]], f$results$taxon)
                   else stats::setNames(rep(NA_real_, nrow(f$results)),
                                        f$results$taxon) },
      stats::setNames(f$results$effect_size, f$results$taxon))
    full
  })

  pairs <- list()
  if (length(fits) >= 2) {
    cmb <- utils::combn(nm, 2, simplify = FALSE)
    for (p in cmb) {
      shared <- intersect(names(eff[[p[1]]]), names(eff[[p[2]]]))
      shared <- shared[!is.na(eff[[p[1]]][shared]) & !is.na(eff[[p[2]]][shared])]
      rho <- if (length(shared) >= 4)
        suppressWarnings(stats::cor(eff[[p[1]]][shared], eff[[p[2]]][shared],
                                    method = "spearman")) else NA_real_
      j <- length(intersect(sig[[p[1]]], sig[[p[2]]])) /
           max(1, length(union(sig[[p[1]]], sig[[p[2]]])))
      pairs[[paste(p, collapse = " vs ")]] <- list(
        spearman_effect = rho, jaccard_hits = j,
        n_shared_taxa = length(shared),
        n_both_significant = length(intersect(sig[[p[1]]], sig[[p[2]]])))
    }
  }

  out <- list(
    methods = nm,
    n_significant = vapply(sig, length, integer(1)),
    consensus = overlap,
    n_consensus = length(overlap),
    union = all_taxa,
    unique_to = lapply(nm, function(m) setdiff(sig[[m]], Reduce(union, sig[nm != m]))),
    pairs = pairs,
    namespace_ok = namespace_ok,
    n_shared_labels_min = pairwise_shared
  )
  names(out$unique_to) <- nm
  mfg_log("da", "methods_compared", list(
    methods = paste(nm, collapse = ","),
    n_significant = paste(sprintf("%s=%d", nm, out$n_significant), collapse = " "),
    n_consensus = length(overlap), n_union = length(all_taxa)))
  class(out) <- c("mfg_da_comparison", "list")
  out
}

#' @export
print.mfg_da_comparison <- function(x, ...) {
  cat("=== DA method comparison ===\n")
  if (!isTRUE(x$namespace_ok)) {
    cat(strwrap(paste("!! THE METHODS SHARE NO TAXON LABELS. Every overlap below",
      "is zero for naming reasons, not statistical ones. Reconcile the labels —",
      "same object, same rank — before reading anything from this table."),
      width = 78, prefix = "  "), sep = "\n")
    cat("\n")
  }
  for (m in x$methods) {
    cat(sprintf("  %-12s %d significant, %d unique to it\n", m,
                x$n_significant[[m]], length(x$unique_to[[m]])))
  }
  cat(sprintf("\nConsensus (significant in all %d methods): %d taxa\n",
              length(x$methods), x$n_consensus))
  if (x$n_consensus) {
    cat(sprintf("  %s\n", paste(utils::head(x$consensus, 15), collapse = ", ")))
  }
  cat(sprintf("Union: %d taxa\n", length(x$union)))
  if (length(x$pairs)) {
    cat("\nPairwise agreement:\n")
    for (p in names(x$pairs)) {
      pr <- x$pairs[[p]]
      cat(sprintf("  %-26s Jaccard(hits) = %.3f | Spearman(effect) = %s | both sig: %d\n",
                  p, pr$jaccard_hits,
                  if (is.na(pr$spearman_effect)) "NA" else sprintf("%.3f", pr$spearman_effect),
                  pr$n_both_significant))
    }
  }
  if (isTRUE(x$namespace_ok)) {
    cat(strwrap(paste("\nReport the consensus set as the primary finding and the",
      "method-specific sets as method-dependent. A taxon found by one method only",
      "is a hypothesis; a taxon found by all is a result."),
      width = 78), sep = "\n")
  }
  invisible(x)
}
