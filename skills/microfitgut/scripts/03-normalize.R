# =============================================================================
# MicroFitGut — 03-normalize.R
#
# The choice that most often silently changes conclusions.
#
# Sequencing depth is an artefact of the run, not biology, and it varies many
# fold across samples. Every comparative statistic therefore needs depth handled
# — but each method distorts the data differently, and a method that is right for
# one downstream test is wrong for another. Rarefied counts must not be fed to
# DESeq2; CLR values must not be fed to a count model; relative abundances must
# not be fed to Chao1.
#
# So this file does two things: it implements the transformations, and it
# enforces the pairing. mfg_check_normalization() refuses invalid combinations
# rather than letting them run.
#
# The decision table lives in reference/03-normalization.md.
#
# Requires: mfg_require(c("da")) for VST; utils.R
# =============================================================================

# ── What each downstream analysis accepts ────────────────────────────────────
#
# Keyed by analysis, listing the normalizations that are valid for it. This is
# the machine-readable form of the reference/03 decision table; if they ever
# disagree, the reference file is the authority and this is the bug.

MFG_VALID_NORMALIZATION <- list(
  # Richness estimators count how many taxa were seen once or twice, so they
  # need integer counts at equal depth. Nothing else is defensible.
  alpha_richness      = c("rarefied"),
  # Shannon/Simpson are defined on proportions and are far less depth-sensitive
  # than richness, so TSS is acceptable; rarefying is still the safer default.
  alpha_diversity     = c("rarefied", "tss"),
  faith_pd            = c("rarefied"),

  # Distances respond to depth. Bray-Curtis and Jaccard on raw counts encode
  # library size as if it were biology.
  beta_bray           = c("tss", "rarefied", "vst", "clr"),
  beta_jaccard        = c("tss", "rarefied", "presence"),
  beta_unifrac        = c("rarefied", "tss"),
  beta_aitchison      = c("clr"),
  beta_euclidean      = c("clr", "vst"),

  # Count models estimate their own size factors internally. Pre-normalizing
  # breaks the variance model they depend on.
  da_deseq2           = c("raw"),
  da_ancombc2         = c("raw"),
  da_aldex2           = c("raw"),
  # LinDA fits a linear model on CLR-transformed data and estimates the bias
  # itself, so unlike the count models above it accepts either counts or
  # proportions. It must be told which it was given. What it cannot take is
  # input that has already been CLR-transformed or rarefied: the first would
  # transform twice, the second discards the counts its winsorization needs.
  da_linda            = c("raw", "tss"),
  da_wilcoxon         = c("clr", "tss", "log10", "rarefied"),
  da_kruskal          = c("clr", "tss", "log10", "rarefied"),

  # Zero-inflated and hurdle models are count models with an offset for depth.
  models_count        = c("raw"),
  models_gaussian     = c("clr", "log10", "tss"),

  ordination_pca      = c("clr", "vst", "tss"),
  ordination_pcoa     = c("tss", "rarefied", "clr", "vst"),
  ordination_nmds     = c("tss", "rarefied"),
  ordination_cca      = c("tss", "raw"),

  composition_barplot = c("tss"),
  core_microbiome     = c("tss"),
  heatmap             = c("clr", "vst", "log10")
)

MFG_NORMALIZATION_LABELS <- c(
  raw       = "raw read counts (no normalization)",
  tss       = "total sum scaling (relative abundance)",
  rarefied  = "rarefied to even depth",
  clr       = "centred log-ratio transformed",
  vst       = "DESeq2 variance-stabilizing transformation",
  log10     = "log10(1 + x) of relative abundance",
  presence  = "presence/absence",
  deseq_rle = "DESeq2 median-of-ratios normalized counts"
)

#' What normalization is currently on this object?
#'
#' Tracked as an attribute, set by every function in this file. Returned as
#' "unknown" when an object arrives from outside the pipeline, which is itself
#' worth reporting: an untracked object cannot be validated.
mfg_normalization <- function(ps) {
  attr(ps, "mfg_normalization") %||% mfg_registry_get(ps, "normalization") %||% {
    mat <- as(phyloseq::otu_table(ps), "matrix")
    if (any(mat < 0, na.rm = TRUE)) "clr_or_vst_untracked"
    else if (looks_like_relative_abundance(mat)) "tss"
    else if (all(mat %% 1 == 0, na.rm = TRUE)) "raw_or_rarefied_untracked"
    else "unknown"
  }
}

mfg_set_normalization <- function(ps, value, detail = list()) {
  attr(ps, "mfg_normalization") <- value
  attr(ps, "mfg_normalization_detail") <- detail
  # Registry backstop: slot assignment (sample_data(ps)$x <- ...) drops
  # object-level attributes, and the caller never sees an intermediate to copy
  # from, so the state is recorded against a content fingerprint too.
  mfg_registry_set(ps, "normalization", value)
  mfg_registry_set(ps, "normalization_detail", detail)
  # Logged because this is an assertion about the data, not bookkeeping. Called
  # directly at intake it is the analyst's claim about an object built
  # elsewhere, and an unlogged claim is one the report cannot trace.
  mfg_log("normalize", "normalization_set",
          list(value = value,
               detail = if (length(detail))
                 paste(names(detail), unlist(detail), sep = "=", collapse = ",")
               else "none"))
  ps
}

#' Refuse an invalid normalization/analysis pairing.
#'
#' Called at the top of every analysis function. This is the single most useful
#' guard in MicroFitGut: it catches the error class that produces a clean-looking
#' result answering a question nobody asked.
mfg_check_normalization <- function(ps, analysis, strict = TRUE) {
  norm  <- mfg_normalization(ps)
  valid <- MFG_VALID_NORMALIZATION[[analysis]]
  if (is.null(valid)) {
    stop("Unknown analysis key '", analysis, "'. Known keys: ",
         paste(names(MFG_VALID_NORMALIZATION), collapse = ", "), call. = FALSE)
  }

  if (grepl("untracked|unknown", norm)) {
    msg <- sprintf(paste0(
      "This object's normalization is not tracked (inferred: %s).\n",
      "'%s' requires one of: %s.\n",
      "Run it through 03-normalize.R, or assert the current state with\n",
      "  ps <- mfg_set_normalization(ps, \"<one of the above>\")\n",
      "so the report can state what was done."), norm, analysis,
      paste(valid, collapse = ", "))
    if (strict) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
    return(invisible(FALSE))
  }

  if (!norm %in% valid) {
    msg <- sprintf(paste0(
      "Invalid normalization for this analysis.\n",
      "  object is: %s (%s)\n",
      "  '%s' accepts: %s\n",
      "  why: %s\n",
      "See reference/03-normalization.md."),
      norm, MFG_NORMALIZATION_LABELS[[norm]] %||% norm, analysis,
      paste(valid, collapse = ", "), mfg_normalization_reason(analysis, norm))
    if (strict) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
    return(invisible(FALSE))
  }

  mfg_log("normalize", "pairing_checked", list(analysis = analysis, normalization = norm))
  invisible(TRUE)
}

#' Why a particular pairing is refused, in one sentence for the error message.
mfg_normalization_reason <- function(analysis, norm) {
  if (grepl("^alpha_richness|^faith_pd", analysis) && norm != "rarefied") {
    return(paste("richness estimators count singletons and doubletons, so they",
                 "measure sequencing effort unless depth is equalized"))
  }
  if (analysis == "da_linda") {
    return(paste("LinDA CLR-transforms internally and estimates the bias term",
                 "itself; give it counts or proportions and declare which, but",
                 "never pre-transformed or rarefied input"))
  }
  if (grepl("^da_deseq2|^da_ancombc2|^da_aldex2|^models_count", analysis)) {
    return(paste("this method estimates its own size factors from raw counts;",
                 "pre-normalized input breaks its variance model and invalidates",
                 "the dispersion estimates"))
  }
  if (grepl("^beta|^ordination", analysis) && norm == "raw") {
    return(paste("distances computed on raw counts encode library size as if it",
                 "were community difference"))
  }
  if (analysis == "beta_aitchison" && norm != "clr") {
    return("Aitchison distance is by definition Euclidean distance on CLR values")
  }
  if (grepl("^models_gaussian", analysis) && norm == "raw") {
    return(paste("a Gaussian model on raw counts assumes constant variance,",
                 "which count data violates"))
  }
  "the method's assumptions do not hold for data in this state"
}

# ── The transformations ──────────────────────────────────────────────────────

#' Total sum scaling: proportions within each sample.
#'
#' The simplest depth correction and the right one for composition plots and
#' most distances. It does not address compositionality — proportions are still
#' constrained to sum to 1, so a rise in one taxon forces others down. For
#' anything inferential on individual taxa, prefer CLR or a count model.
tss_transform <- function(ps) {
  out <- phyloseq::transform_sample_counts(ps, function(x) if (sum(x) > 0) x / sum(x) else x)
  out <- mfg_carry_attrs(out, ps)
  mfg_log("normalize", "tss", list(n_samples = phyloseq::nsamples(out)))
  mfg_set_normalization(out, "tss")
}

# Kept under the app's name so existing code reads unchanged.
compositional <- tss_transform

#' Centred log-ratio transform.
#'
#' Moves compositional data into Euclidean space so ordinary statistics apply.
#' microbiome::transform("clr") adds a pseudocount internally for zeros. CLR is
#' relative to each sample's geometric mean, so a CLR value is a statement about
#' a taxon *relative to the rest of that sample*, not an absolute abundance —
#' state that in the report when interpreting coefficients.
clr_transform <- function(ps) {
  out <- mfg_carry_attrs(microbiome::transform(ps, "clr"), ps)
  mfg_log("normalize", "clr", list(n_taxa = phyloseq::ntaxa(out)))
  mfg_set_normalization(out, "clr")
}

#' log10(1 + x) of relative abundance (Demo 9, option 1).
#'
#' Cheaper than CLR and easier to explain, but it does not remove the
#' compositional constraint. Fine for heatmaps and visual scaling; weaker than
#' CLR for inference.
log10_transform <- function(ps, already_relative = FALSE) {
  base <- if (already_relative) ps else tss_transform(ps)
  out <- mfg_carry_attrs(microbiome::transform(base, "log10"), ps)
  mfg_log("normalize", "log10", list(note = "log10(1+x) applied; zeros present"))
  mfg_set_normalization(out, "log10")
}

#' Presence/absence.
presence_transform <- function(ps, detection = 0) {
  out <- mfg_carry_attrs(
    phyloseq::transform_sample_counts(ps, function(x) as.numeric(x > detection)), ps)
  mfg_log("normalize", "presence", list(detection = detection))
  mfg_set_normalization(out, "presence")
}

# ── Rarefaction ──────────────────────────────────────────────────────────────

#' Choose a rarefaction depth and show the cost of the choice.
#'
#' Rarefying to the minimum depth keeps every sample but throws away the most
#' reads; rarefying higher keeps more reads but drops samples. There is no
#' correct answer, only a documented one. Returns a table of the trade-off at
#' candidate depths so the choice can be made visibly and defended.
choose_rarefaction_depth <- function(ps, quantiles = c(0.05, 0.10, 0.25, 0.50),
                                     candidates = NULL, verbose = TRUE) {
  depths <- sort(phyloseq::sample_sums(ps))
  n <- length(depths)
  if (is.null(candidates)) {
    candidates <- unique(c(min(depths), stats::quantile(depths, quantiles, names = FALSE)))
    candidates <- unique(round(candidates))
  }
  rows <- lapply(candidates, function(d) {
    kept <- sum(depths >= d)
    data.frame(depth = d, samples_kept = kept, samples_dropped = n - kept,
               pct_samples_kept = round(100 * kept / n, 1),
               reads_retained = d * kept,
               pct_reads_retained = round(100 * d * kept / sum(depths), 1),
               stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  tab <- tab[order(tab$depth), ]
  if (verbose) {
    cat("Rarefaction depth trade-off:\n")
    print(tab, row.names = FALSE)
    cat("\nNo depth is correct. Pick one, record it, and check the rarefaction\n",
        "curves at that depth actually plateau (mfg_rarefaction_curve).\n", sep = "")
  }
  mfg_log("normalize", "rarefaction_depth_candidates",
          list(min_depth = min(depths), max_depth = max(depths),
               candidates = paste(candidates, collapse = ",")))
  tab
}

#' Rarefy to a fixed depth, with the seed recorded.
#'
#' Rarefaction subsamples at random, so an unseeded run is not reproducible and
#' two runs will disagree on richness. `rngseed` is required and logged.
#' `replace = FALSE` is sampling without replacement, which is what
#' "rarefaction" means; phyloseq's default of TRUE is faster but is not the
#' classical procedure.
rarefy_to_depth <- function(ps, depth = NULL, rngseed = 42, replace = FALSE,
                            justification = NULL, verbose = TRUE) {
  mat <- as(phyloseq::otu_table(ps), "matrix")
  if (looks_like_relative_abundance(mat)) {
    stop("Cannot rarefy relative abundances — there are no reads to subsample. ",
         "Rarefaction requires raw counts. See reference/03.", call. = FALSE)
  }
  if (!all(mat %% 1 == 0, na.rm = TRUE)) {
    stop("Cannot rarefy non-integer values. Rarefaction subsamples reads, so the ",
         "table must hold counts.", call. = FALSE)
  }
  if (is.null(depth)) {
    depth <- min(phyloseq::sample_sums(ps))
    if (verbose) message(sprintf("No depth given; using the minimum sample depth (%d). ",
                                 depth),
                         "This keeps every sample but discards the most reads.")
  }
  n_before <- phyloseq::nsamples(ps)
  # phyloseq::rarefy_even_depth saves and restores .Random.seed, which does not
  # exist yet in a fresh non-interactive session. Initialise it first so the
  # call does not fail on a bare Rscript run.
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) set.seed(rngseed)
  out <- mfg_carry_attrs(
    phyloseq::rarefy_even_depth(ps, sample.size = depth, replace = replace,
                                rngseed = rngseed, verbose = verbose), ps)
  dropped <- setdiff(phyloseq::sample_names(ps), phyloseq::sample_names(out))

  mfg_log("normalize", "rarefied", list(
    depth = depth, rngseed = rngseed, replace = replace,
    samples = sprintf("%d->%d", n_before, phyloseq::nsamples(out)),
    dropped_samples = paste(dropped, collapse = ";"),
    justification = justification %||% "NOT RECORDED"))

  if (length(dropped) && verbose) {
    message(sprintf("Rarefying to %d dropped %d sample(s): %s", depth, length(dropped),
                    paste(utils::head(dropped, 10), collapse = ", ")))
  }
  if (is.null(justification)) {
    warning("rarefy_to_depth called without a justification. reference/10 requires ",
            "the report to state the depth and why it was chosen.", call. = FALSE)
  }
  mfg_set_normalization(out, "rarefied",
                        detail = list(depth = depth, rngseed = rngseed,
                                      dropped = dropped))
}

#' Per-sample rarefaction curves, to check a depth actually plateaus.
#'
#' A curve still climbing at the chosen depth means richness is measuring
#' sequencing effort, not the community. Returns the data so the check is
#' recorded rather than eyeballed and forgotten.
#'
#' IMPORTANT: call this on the object at the depth you intend to use — usually a
#' rarefied one. On an unrarefied object every curve runs out to that sample's own
#' total depth, so each one necessarily ends flat and "has it plateaued" becomes
#' trivially true. That is not the question; the question is whether a curve is
#' still climbing *at the depth you chose*.
#'
#' Plateau is judged by the **terminal slope in taxa per 1000 reads**, not by the
#' slope relative to the start of the curve. The start of any rarefaction curve is
#' steep — the first reads are nearly all new taxa — so a ratio against the initial
#' slope declares a plateau almost regardless of saturation. The absolute terminal
#' slope answers the question that matters: how many new taxa would another 1000
#' reads reveal?
#'
#' `plateau_slope` defaults to 1 taxon per 1000 reads. That is a convention, not a
#' law; it is reported so a different threshold can be applied and defended.
mfg_rarefaction_curve <- function(ps, step = 500, max_samples = 60, seed = 42,
                                  plateau_slope = 1, tail_fraction = 0.2) {
  mat <- as(phyloseq::otu_table(ps), "matrix")
  if (phyloseq::taxa_are_rows(ps)) mat <- t(mat)
  set.seed(seed)
  if (nrow(mat) > max_samples) {
    idx <- sample(seq_len(nrow(mat)), max_samples)
    mat <- mat[idx, , drop = FALSE]
    message(sprintf("Curves computed on a random %d-sample subset for legibility.",
                    max_samples))
  }
  storage.mode(mat) <- "integer"
  curves <- vegan::rarecurve(mat, step = step, tidy = TRUE)
  names(curves) <- c("sample", "depth", "richness")

  # Terminal slope over the final `tail_fraction` of each curve, expressed as
  # taxa gained per 1000 additional reads. The secondary slope_ratio against the
  # initial slope is kept as a diagnostic but does NOT decide the verdict.
  flat <- do.call(rbind, lapply(split(curves, curves$sample), function(d) {
    d <- d[order(d$depth), ]
    n <- nrow(d)
    if (n < 4) {
      return(data.frame(sample = d$sample[1], terminal_slope_per_1k = NA_real_,
                        slope_ratio = NA_real_, final_richness = d$richness[n],
                        max_depth = d$depth[n], stringsAsFactors = FALSE))
    }
    i <- which.min(abs(d$depth - (1 - tail_fraction) * max(d$depth)))
    i <- min(i, n - 1)
    d_depth <- d$depth[n] - d$depth[i]
    term <- if (d_depth > 0) 1000 * (d$richness[n] - d$richness[i]) / d_depth else NA_real_

    head_i <- seq(1, max(2, floor(n / 5)))
    s_head <- diff(range(d$richness[head_i])) / max(1, diff(range(d$depth[head_i])))
    ratio  <- if (!is.na(term) && s_head > 0) (term / 1000) / s_head else NA_real_

    data.frame(sample = d$sample[1], terminal_slope_per_1k = term,
               slope_ratio = ratio, final_richness = d$richness[n],
               max_depth = d$depth[n], stringsAsFactors = FALSE)
  }))

  n_plateau <- sum(flat$terminal_slope_per_1k < plateau_slope, na.rm = TRUE)
  n_total   <- sum(!is.na(flat$terminal_slope_per_1k))
  med_slope <- stats::median(flat$terminal_slope_per_1k, na.rm = TRUE)
  frac      <- if (n_total > 0) n_plateau / n_total else NA_real_

  verdict <- if (is.na(frac)) "curves too short to judge"
    else if (frac > 0.8) sprintf(
      "most samples plateau (median terminal slope %.2f taxa per 1000 reads)", med_slope)
    else sprintf(paste("many samples have NOT plateaued — richness is depth-limited",
      "(median terminal slope %.2f taxa per 1000 reads; %d of %d below the",
      "%.2f threshold)"), med_slope, n_plateau, n_total, plateau_slope)

  depths <- range(flat$max_depth, na.rm = TRUE)
  if (diff(depths) > 0.01 * max(depths)) {
    warning(sprintf(paste("Curve end-depths range from %.0f to %.0f, so these curves",
      "were computed on an UNRAREFIED object: each one stops at its own sample's",
      "total depth and therefore ends flat by construction. Plateau is not",
      "meaningful here. Rarefy first, then call this on the rarefied object."),
      depths[1], depths[2]), call. = FALSE)
  }

  mfg_log("normalize", "rarefaction_curves", list(
    n_samples = length(unique(curves$sample)), step = step,
    plateau_threshold_per_1k = plateau_slope,
    median_terminal_slope = round(med_slope, 3),
    n_plateaued = n_plateau, n_judged = n_total,
    depth_range = sprintf("%.0f-%.0f", depths[1], depths[2]),
    verdict = verdict))

  list(curves = curves, plateau = flat, verdict = verdict,
       median_terminal_slope = med_slope, plateau_threshold = plateau_slope,
       n_plateaued = n_plateau, n_total = n_total)
}

# ── DESeq2-based normalization ───────────────────────────────────────────────

#' Median-of-ratios size factors, using the positive-count geometric mean.
#'
#' DESeq2's default size-factor estimation needs at least one taxon present in
#' every sample, which microbiome data rarely has — every row usually contains a
#' zero, so the default geometric mean is zero everywhere and estimation fails.
#' Restricting the product to positive counts is the standard workaround
#' (Demo 10, ps10).
deseq2_size_factors <- function(ps, design = ~ 1) {
  dds <- phyloseq::phyloseq_to_deseq2(ps, design)
  geo <- apply(DESeq2::counts(dds), 1, gm_mean)
  dds <- DESeq2::estimateSizeFactors(dds, geoMeans = geo)
  sf  <- DESeq2::sizeFactors(dds)
  mfg_log("normalize", "deseq2_size_factors",
          list(n = length(sf), range = sprintf("%.3f-%.3f", min(sf), max(sf)),
               method = "median-of-ratios with positive-count geometric means"))
  list(dds = dds, size_factors = sf)
}

#' Replace the count table with DESeq2-normalized counts.
#'
#' Useful for plots and for handing normalized counts to a method that expects
#' them, but note: do NOT then run DESeq2 on this object. DESeq2 must see raw
#' counts so it can model the mean-variance relationship itself.
deseq2_normalize <- function(ps, design = ~ 1, round_counts = TRUE) {
  fit <- deseq2_size_factors(ps, design)
  nc  <- DESeq2::counts(fit$dds, normalized = TRUE)
  if (round_counts) nc <- round(nc, 0)
  out <- mfg_carry_attrs(ps, ps)
  phyloseq::otu_table(out) <- phyloseq::otu_table(nc, taxa_are_rows = TRUE)
  mfg_log("normalize", "deseq2_normalized",
          list(rounded = round_counts,
               warning = "do not run DESeq2 on this object; it needs raw counts"))
  mfg_set_normalization(out, "deseq_rle")
}

#' Variance-stabilizing transformation, for ordination and heatmaps.
#'
#' VST removes the dependence of variance on the mean, which is what makes PCA
#' and Euclidean distance behave on count data. `blind = TRUE` keeps the
#' transformation independent of the design, which is what you want when the
#' transformed values feed an unsupervised method — otherwise the ordination is
#' partly told the answer.
vst_transform <- function(ps, design = ~ 1, blind = TRUE, fitType = "local") {
  dds <- phyloseq::phyloseq_to_deseq2(ps, design)
  geo <- apply(DESeq2::counts(dds), 1, gm_mean)
  dds <- DESeq2::estimateSizeFactors(dds, geoMeans = geo)
  vs  <- DESeq2::varianceStabilizingTransformation(dds, blind = blind, fitType = fitType)
  m   <- SummarizedExperiment::assay(vs)

  out <- mfg_carry_attrs(ps, ps)
  phyloseq::otu_table(out) <- phyloseq::otu_table(m, taxa_are_rows = TRUE)
  mfg_log("normalize", "vst",
          list(blind = blind, fitType = fitType,
               note = "negative values are expected and valid for VST output"))
  mfg_set_normalization(out, "vst", detail = list(blind = blind, fitType = fitType))
}

# ── Reporting ────────────────────────────────────────────────────────────────

#' One sentence naming the normalization, for the methods section.
mfg_normalization_sentence <- function(ps) {
  n <- mfg_normalization(ps)
  d <- attr(ps, "mfg_normalization_detail")
  base <- MFG_NORMALIZATION_LABELS[[n]] %||% n
  if (identical(n, "rarefied") && !is.null(d$depth)) {
    return(sprintf("Counts were rarefied without replacement to %s reads per sample (seed %s), which retained %d of %d samples.",
                   format(d$depth, big.mark = ","), d$rngseed,
                   phyloseq::nsamples(ps), phyloseq::nsamples(ps) + length(d$dropped %||% character())))
  }
  if (identical(n, "vst")) {
    return(sprintf("Counts were transformed with the DESeq2 variance-stabilizing transformation (blind = %s, fitType = '%s').",
                   d$blind, d$fitType))
  }
  sprintf("Abundances were expressed as %s.", base)
}
