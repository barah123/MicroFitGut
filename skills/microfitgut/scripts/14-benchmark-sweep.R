# =============================================================================
# MicroFitGut — 14-benchmark-sweep.R
#
# Turning a difference into an explanation.
#
# When a reanalysis disagrees with a paper, the useful output is not "we got 41
# of your 58 taxa" — it is "17 of the 58 were only significant under one DA
# method; under ANCOM-BC2 on the same counts they fall out." Getting there needs
# two things:
#
#   1. A sweep over the decisions the reanalysis had to guess. Every parameter
#      the paper did not state becomes an axis. One-at-a-time from the
#      replication baseline, not full factorial: 12-15 runs instead of 200, and
#      enough for marginal attribution.
#
#   2. A ranking of which axis moved the concordance metric most — computed in
#      code, so the attribution is measured rather than asserted.
#
# The escape valve matters as much as the ranking: when no tested axis explains
# the gap, "no tested factor accounts for this" is the correct answer, and the
# unexplained residual is reported rather than filled with a plausible story.
#
# Requires: utils.R; 03-normalize.R; the analysis scripts for whatever is swept
# =============================================================================

# ── Study card ───────────────────────────────────────────────────────────────
#
# What the paper said, what it did not say, and where each field came from. The
# "absent" values are the point: an unstated parameter is a benchmark finding,
# and it automatically becomes a sweep axis.

MFG_STUDY_CARD_FIELDS <- list(
  # --- Identification ---
  accessions        = list(type = "character", sweep = FALSE),
  n_samples         = list(type = "integer",   sweep = FALSE),
  design            = list(type = "character", sweep = FALSE),
  group_variable    = list(type = "character", sweep = FALSE),
  subject_variable  = list(type = "character", sweep = FALSE),

  # --- Upstream (recorded but not swept downstream) ---
  amplicon_region   = list(type = "character", sweep = FALSE,
                           options = c("V1-V2","V3-V4","V4","V4-V5","ITS1","shotgun","other")),
  denoiser          = list(type = "character", sweep = FALSE,
                           options = c("dada2","deblur","uparse","mothur_otu","metaphlan","kraken2","unknown")),
  reference_db      = list(type = "character", sweep = FALSE,
                           options = c("silva_138","silva_138_2","gg2_2024_09","gtdb_r220","rdp_19","unknown")),

  # --- Downstream: every one of these is a sweep axis when unstated ---
  prevalence_filter = list(type = "numeric", sweep = TRUE,
                           options = c(0, 0.05, 0.10, 0.20)),
  min_depth         = list(type = "integer", sweep = TRUE,
                           options = c(0, 1000, 5000, 10000)),
  normalization     = list(type = "character", sweep = TRUE,
                           options = c("rarefied","tss","clr","vst","raw")),
  rarefaction_depth = list(type = "integer", sweep = TRUE,
                           options = NULL),   # filled from the data's depth distribution
  taxonomic_level   = list(type = "character", sweep = TRUE,
                           options = c("ASV","Genus","Family","Phylum")),
  alpha_metric      = list(type = "character", sweep = FALSE),
  beta_distance     = list(type = "character", sweep = TRUE,
                           options = c("bray","jaccard","unifrac","wunifrac","jsd")),
  da_method         = list(type = "character", sweep = TRUE,
                           options = c("ancombc2","aldex2","deseq2","kruskal")),
  fdr_method        = list(type = "character", sweep = TRUE,
                           options = c("BH","holm","BY","none")),
  alpha_threshold   = list(type = "numeric", sweep = TRUE,
                           options = c(0.05, 0.1)),
  covariates        = list(type = "character", sweep = FALSE),
  random_effects    = list(type = "character", sweep = TRUE,
                           options = c("none","subject")),

  # --- The ground truth ---
  reported_alpha    = list(type = "list", sweep = FALSE),
  reported_permanova = list(type = "list", sweep = FALSE),
  reported_da_taxa  = list(type = "character", sweep = FALSE),
  reported_claims   = list(type = "character", sweep = FALSE)
)

#' Create a study card field with its provenance.
#'
#' confidence is the load-bearing part:
#'   "stated"   the paper says it — fixed, not swept
#'   "inferred" deduced from something else the paper says — swept, with the
#'              inference recorded
#'   "absent"   the paper does not say — swept across all options
#'
#' evidence should be a verbatim span of at most about 25 words, and location a
#' pointer such as "Methods 2.3" or "Supp. Table S1", so a reader can check the
#' extraction itself.
study_field <- function(value = NULL, confidence = c("stated", "inferred", "absent"),
                        evidence = NULL, location = NULL) {
  confidence <- match.arg(confidence)
  if (identical(confidence, "stated") && is.null(value)) {
    stop("confidence = 'stated' requires a value. Use 'absent' when the paper ",
         "does not report the parameter.", call. = FALSE)
  }
  if (identical(confidence, "stated") && is.null(evidence)) {
    warning("A 'stated' field without evidence cannot be checked. Quote the span ",
            "the value came from.", call. = FALSE)
  }
  list(value = value, confidence = confidence, evidence = evidence, location = location)
}

#' A single machine-checkable claim from the paper.
#'
#' `reported_claims` as free prose cannot drive a verdict: scoring it would mean
#' a model reading a sentence and deciding whether it held, which puts an
#' unreproducible judgement on the critical path of the primary endpoint. A claim
#' has to arrive already in a form a function can evaluate.
#'
#' `contrast` is load-bearing and not optional. A paper states several claims
#' about different comparisons, and scoring a claim against the wrong one
#' produces a confident wrong answer with no symptom — the arithmetic succeeds.
#' Scoring refuses unless the claim's contrast matches the analysis it is given.
#'
#' type            required fields
#' --------------  ----------------------------------------------------------
#' da_count        n, comparator ("eq","lte","gte","approx"); tolerance for approx
#' da_direction    taxon, higher_in
#' dominance       taxon, group
#' alpha           direction ("higher","lower","none"), groups
#' beta            differs (TRUE/FALSE)
#' unscored        text — explicitly not machine-checkable
#'
#' `unscored` is the escape valve, and it is the honest answer for a qualitative
#' claim. What it must never do is disappear: an unscored claim blocks a clean
#' "reproduced" verdict rather than being silently dropped.
study_claim <- function(id, type = c("da_count", "da_direction", "dominance",
                                     "alpha", "beta", "unscored"),
                        contrast = NULL, evidence = NULL, location = NULL, ...) {
  type <- match.arg(type)
  spec <- list(...)
  required <- switch(type,
    da_count     = c("n", "comparator"),
    da_direction = c("taxon", "higher_in"),
    dominance    = c("taxon", "group"),
    alpha        = c("direction"),
    beta         = c("differs"),
    unscored     = c("text"))
  missing <- setdiff(required, names(spec))
  if (length(missing)) {
    stop(sprintf("study_claim('%s', type = '%s') needs: %s", id, type,
                 paste(missing, collapse = ", ")), call. = FALSE)
  }
  if (!identical(type, "unscored") && is.null(contrast)) {
    stop(sprintf(paste("study_claim('%s') needs a contrast. A claim scored",
      "against the wrong comparison fails silently."), id), call. = FALSE)
  }
  if (is.null(evidence)) {
    warning(sprintf("Claim '%s' has no evidence span and cannot be checked.", id),
            call. = FALSE)
  }
  out <- c(list(id = id, type = type, contrast = contrast,
                evidence = evidence, location = location), spec)
  class(out) <- c("mfg_study_claim", "list")
  out
}

#' Build and validate a study card.
study_card <- function(..., study_id = NULL) {
  fields <- list(...)
  unknown <- setdiff(names(fields), names(MFG_STUDY_CARD_FIELDS))
  if (length(unknown)) {
    warning("Study card fields not in the schema (kept, but they will not be swept): ",
            paste(unknown, collapse = ", "), call. = FALSE)
  }
  # Anything the schema knows about and the card omits is 'absent', not missing —
  # which is what makes it a sweep axis rather than an oversight.
  for (f in names(MFG_STUDY_CARD_FIELDS)) {
    if (is.null(fields[[f]])) {
      fields[[f]] <- study_field(confidence = "absent")
    }
  }
  card <- list(study_id = study_id %||% "unnamed-study", fields = fields)

  conf <- vapply(fields[names(MFG_STUDY_CARD_FIELDS)],
                 function(f) f$confidence %||% "absent", character(1))
  card$n_stated   <- sum(conf == "stated")
  card$n_inferred <- sum(conf == "inferred")
  card$n_absent   <- sum(conf == "absent")

  mfg_log("benchmark", "study_card_built", list(
    study_id = card$study_id, n_stated = card$n_stated,
    n_inferred = card$n_inferred, n_absent = card$n_absent))
  class(card) <- c("mfg_study_card", "list")
  card
}

#' @export
print.mfg_study_card <- function(x, ...) {
  cat(sprintf("=== Study card: %s ===\n", x$study_id))
  cat(sprintf("%d fields stated, %d inferred, %d absent\n\n",
              x$n_stated, x$n_inferred, x$n_absent))
  for (nm in names(MFG_STUDY_CARD_FIELDS)) {
    f <- x$fields[[nm]]
    if (is.null(f)) next
    v <- if (is.null(f$value)) "-" else paste(format(unlist(f$value)), collapse = ", ")
    if (nchar(v) > 44) v <- paste0(substr(v, 1, 41), "...")
    cat(sprintf("  %-20s %-46s [%s]\n", nm, v, f$confidence %||% "absent"))
    if (!is.null(f$evidence)) {
      cat(sprintf("  %-20s   \"%s\"%s\n", "", f$evidence,
                  if (!is.null(f$location)) paste0(" — ", f$location) else ""))
    }
  }
  cat(strwrap(paste("\nEvery field marked 'inferred' or 'absent' becomes an axis in",
    "the sensitivity sweep. An unstated parameter is itself a finding about the",
    "paper's reporting and should appear in the benchmark report as one."),
    width = 78), sep = "\n")
  invisible(x)
}

# ── Sweep grid ───────────────────────────────────────────────────────────────

#' Build the one-at-a-time sweep grid from a study card.
#'
#' Baseline is the best reconstruction of the paper. Each further run changes
#' exactly one axis, so the effect of that axis is isolated. Full factorial over
#' five axes is 200+ runs; one-at-a-time is 12-15 and enough for marginal
#' attribution. Add a small factorial only for axes known to interact —
#' denoiser x DA method is the one that genuinely does.
build_sweep_grid <- function(card, ps = NULL, axes = NULL, max_runs = 20) {
  baseline <- list()
  for (nm in names(MFG_STUDY_CARD_FIELDS)) {
    f <- card$fields[[nm]]
    if (!is.null(f$value)) baseline[[nm]] <- f$value
  }
  # Sensible baseline for anything absent, so the replication run is runnable.
  # normalization is deliberately NOT fixed here: it depends on the DA method.
  defaults <- list(prevalence_filter = 0.10, min_depth = 0,
                   taxonomic_level = "ASV", beta_distance = "bray",
                   da_method = "ancombc2", fdr_method = "BH",
                   alpha_threshold = 0.05, random_effects = "none")
  for (nm in names(defaults)) {
    if (is.null(baseline[[nm]])) baseline[[nm]] <- defaults[[nm]]
  }
  # Every count-based DA method requires raw counts; the rank-based ones require
  # a transform and refuse raw. A baseline that pairs a DA method with a
  # normalization it does not accept is refused before it computes anything, so
  # the default is derived from the method rather than hardcoded. Hardcoding
  # "tss" alongside a default of ancombc2 made every DA benchmark unrunnable.
  da_norms <- function(m) MFG_VALID_NORMALIZATION[[paste0("da_", m)]]
  if (is.null(baseline$normalization)) {
    vn <- da_norms(baseline$da_method)
    baseline$normalization <- if (length(vn)) vn[[1]] else "tss"
  } else {
    vn <- da_norms(baseline$da_method)
    if (length(vn) && !baseline$normalization %in% vn) {
      warning(sprintf(paste("The card pairs normalization '%s' with da_method",
        "'%s', which accepts only: %s. Every DA run will be refused. Check the",
        "extraction, or state the DA method the paper actually used."),
        baseline$normalization, baseline$da_method,
        paste(vn, collapse = ", ")), call. = FALSE)
    }
  }

  # Which axes to sweep: those the paper left inferred or absent.
  sweepable <- names(MFG_STUDY_CARD_FIELDS)[vapply(names(MFG_STUDY_CARD_FIELDS),
    function(nm) isTRUE(MFG_STUDY_CARD_FIELDS[[nm]]$sweep), logical(1))]
  uncertain <- names(card$fields)[vapply(names(card$fields), function(nm) {
    f <- card$fields[[nm]]
    (f$confidence %||% "absent") %in% c("inferred", "absent")
  }, logical(1))]
  axes <- axes %||% intersect(sweepable, uncertain)

  # Axes ordered by how much they typically move a concordance metric, so that a
  # run cap costs levels within an axis rather than whole axes. Dropping
  # da_method entirely because prevalence_filter was earlier in the schema would
  # make the attribution useless.
  axis_priority <- c("da_method", "normalization", "random_effects",
                     "prevalence_filter", "rarefaction_depth", "fdr_method",
                     "taxonomic_level", "beta_distance", "min_depth",
                     "alpha_threshold")
  axes <- c(intersect(axis_priority, axes), setdiff(axes, axis_priority))

  # Candidate levels per axis, ordered so the first entry is the most informative
  # contrast against the baseline.
  levels_by_axis <- list()
  unrunnable <- list()
  paired <- list()
  for (ax in axes) {
    opts <- MFG_STUDY_CARD_FIELDS[[ax]]$options
    # Rarefaction depth options come from this dataset's depth distribution, not
    # from a fixed list.
    if (identical(ax, "rarefaction_depth") && !is.null(ps)) {
      d <- phyloseq::sample_sums(ps)
      opts <- unique(round(c(min(d), stats::quantile(d, c(0.10, 0.25), names = FALSE))))
    }
    # Varying normalization while holding the DA method fixed only makes sense
    # over normalizations that method accepts. The rest are refusals, not runs.
    if (identical(ax, "normalization")) {
      vn <- da_norms(baseline$da_method)
      if (length(vn)) {
        blocked <- setdiff(as.character(opts), as.character(vn))
        if (length(blocked)) unrunnable[[ax]] <- blocked
        opts <- opts[as.character(opts) %in% as.character(vn)]
      }
    }
    # A level this dataset cannot run is not a level. UniFrac against an object
    # with no tree, or a rank the taxonomy table does not carry, puts runs in the
    # grid that can only fail — and worse, inflates the denominator of any "share
    # of configurations in which the claim held" statistic computed from it.
    # Unrunnable is a different category from untested and is reported separately.
    if (!is.null(opts) && !is.null(ps)) {
      runnable <- switch(ax,
        beta_distance   = allowed_beta_distances(mfg_tree_is_real(ps)),
        taxonomic_level = c("ASV", phyloseq::rank_names(ps)),
        NULL)
      if (!is.null(runnable)) {
        blocked <- setdiff(as.character(opts), as.character(runnable))
        if (length(blocked)) unrunnable[[ax]] <- blocked
        opts <- opts[as.character(opts) %in% as.character(runnable)]
      }
    }
    if (is.null(opts)) next
    opts <- setdiff(opts, baseline[[ax]])
    if (length(opts)) levels_by_axis[[ax]] <- opts
  }
  if (length(unrunnable)) {
    warning(sprintf(paste("These axis levels cannot run on this dataset and were",
      "removed from the grid. Report them as unrunnable, not untested: %s"),
      paste(sprintf("%s (%s)", names(unrunnable),
                    vapply(unrunnable, paste, character(1), collapse = "/")),
            collapse = "; ")), call. = FALSE)
  }
  # The baseline itself is a run. If the paper's own stated configuration cannot
  # execute on this object, every concordance figure below is against nothing.
  if (!is.null(ps)) {
    bl_blocked <- character(0)
    if (!is.null(baseline$beta_distance) &&
        !baseline$beta_distance %in% allowed_beta_distances(mfg_tree_is_real(ps))) {
      bl_blocked <- c(bl_blocked, paste0("beta_distance=", baseline$beta_distance))
    }
    if (!is.null(baseline$taxonomic_level) &&
        !baseline$taxonomic_level %in% c("ASV", phyloseq::rank_names(ps))) {
      bl_blocked <- c(bl_blocked, paste0("taxonomic_level=", baseline$taxonomic_level))
    }
    if (length(bl_blocked)) {
      warning(sprintf(paste("The BASELINE cannot run on this object (%s). The",
        "baseline is the reconstruction of the paper, so no concordance can be",
        "computed until this is resolved."),
        paste(bl_blocked, collapse = ", ")), call. = FALSE)
    }
  }

  # Round-robin: one level from each axis per pass. Truncation then thins every
  # axis evenly instead of amputating the tail of the axis list.
  runs <- list(data.frame(run_id = "baseline", axis = "none", level = NA_character_,
                          stringsAsFactors = FALSE))
  params <- list(baseline = baseline)
  n_total_levels <- sum(vapply(levels_by_axis, length, integer(1)))
  max_depth <- if (length(levels_by_axis)) max(vapply(levels_by_axis, length, integer(1))) else 0L

  for (i in seq_len(max_depth)) {
    for (ax in names(levels_by_axis)) {
      opts <- levels_by_axis[[ax]]
      if (i > length(opts)) next
      if (length(runs) >= max_runs) break
      o <- opts[i]
      rid <- sprintf("%s=%s", ax, o)
      runs[[length(runs) + 1]] <- data.frame(run_id = rid, axis = ax,
                                             level = as.character(o),
                                             stringsAsFactors = FALSE)
      p <- baseline; p[[ax]] <- o
      # Swapping the DA method can invalidate the baseline normalization —
      # ancombc2 requires raw, kruskal refuses it. Carrying the baseline
      # normalization across would make every such run a refusal, so the
      # normalization moves with it. That is two changes in one run, which
      # weakens the isolation the design depends on, so it is recorded and
      # printed rather than done quietly.
      if (identical(ax, "da_method")) {
        vn <- da_norms(o)
        if (length(vn) && !p$normalization %in% vn) {
          paired[[rid]] <- sprintf("normalization %s -> %s (required by %s)",
                                   p$normalization, vn[[1]], o)
          p$normalization <- vn[[1]]
        }
      }
      params[[rid]] <- p
    }
    if (length(runs) >= max_runs) break
  }
  grid <- do.call(rbind, runs)

  n_scheduled <- nrow(grid) - 1L
  truncated <- max(0L, n_total_levels - n_scheduled)
  untested <- list()
  if (truncated > 0) {
    for (ax in names(levels_by_axis)) {
      tested <- grid$level[grid$axis == ax]
      missed <- setdiff(as.character(levels_by_axis[[ax]]), tested)
      if (length(missed)) untested[[ax]] <- missed
    }
    # A silent cap reads as "we covered everything" when it did not.
    warning(sprintf(paste("Sweep grid capped at %d runs; %d axis levels were NOT",
      "tested and must be listed as untested in the benchmark report: %s"),
      max_runs, truncated,
      paste(sprintf("%s (%s)", names(untested),
                    vapply(untested, paste, character(1), collapse = "/")),
            collapse = "; ")), call. = FALSE)
  }
  # Every axis still gets at least one level tested, which is what marginal
  # attribution requires.
  axes_tested <- setdiff(unique(grid$axis), "none")
  axes_dropped <- setdiff(names(levels_by_axis), axes_tested)
  if (length(axes_dropped)) {
    warning(sprintf(paste("These axes got no run at all and cannot be attributed:",
      "%s. Raise max_runs to at least %d to cover every axis once."),
      paste(axes_dropped, collapse = ", "), length(levels_by_axis) + 1L),
      call. = FALSE)
  }

  mfg_log("benchmark", "sweep_grid_built", list(
    n_runs = nrow(grid), axes = paste(axes, collapse = ","),
    n_axes = length(axes), truncated = truncated,
    unrunnable = if (length(unrunnable))
      paste(sprintf("%s=%s", names(unrunnable),
                    vapply(unrunnable, paste, character(1), collapse = "/")),
            collapse = " ") else "none",
    baseline = paste(sprintf("%s=%s", names(baseline),
                             vapply(baseline, function(v)
                               paste(format(unlist(v)), collapse = "/"), character(1))),
                     collapse = " ")))

  out <- list(grid = grid, params = params, baseline = baseline, axes = axes,
              n_truncated = truncated, untested = untested,
              unrunnable = unrunnable, paired = paired,
              axes_dropped = axes_dropped, n_total_levels = n_total_levels)
  class(out) <- c("mfg_sweep_grid", "list")
  out
}

#' @export
print.mfg_sweep_grid <- function(x, ...) {
  cat(sprintf("=== Sweep grid: %d runs over %d axes ===\n",
              nrow(x$grid), length(x$axes)))
  cat("\nBaseline (best reconstruction of the paper):\n")
  for (nm in names(x$baseline)) {
    cat(sprintf("  %-20s %s\n", nm,
                paste(format(unlist(x$baseline[[nm]])), collapse = ", ")))
  }
  cat("\nRuns:\n"); print(x$grid, row.names = FALSE)
  if (x$n_truncated > 0) {
    cat(sprintf("\n! %d of %d axis levels are UNTESTED (run cap):\n",
                x$n_truncated, x$n_total_levels))
    for (ax in names(x$untested)) {
      cat(sprintf("    %-20s %s\n", ax, paste(x$untested[[ax]], collapse = ", ")))
    }
    cat("  List these in the benchmark report — an untested axis cannot be ruled out.\n")
  }
  if (length(x$unrunnable)) {
    cat("\n! UNRUNNABLE on this dataset (removed from the grid, not a run cap):\n")
    for (ax in names(x$unrunnable)) {
      cat(sprintf("    %-20s %s\n", ax, paste(x$unrunnable[[ax]], collapse = ", ")))
    }
    cat("  Report these separately from untested levels. They were never candidates,\n")
    cat("  so they do not belong in the denominator of a fragility statistic.\n")
  }
  if (length(x$paired)) {
    cat("\n! Runs where a second parameter had to move as well:\n")
    for (r in names(x$paired)) cat(sprintf("    %-22s %s\n", r, x$paired[[r]]))
    cat("  These runs change two things, so their axis is not cleanly isolated.\n")
  }
  if (length(x$axes_dropped)) {
    cat(sprintf("\n!! %s got no run at all and cannot be attributed.\n",
                paste(x$axes_dropped, collapse = ", ")))
  }
  invisible(x)
}

# ── Attribution ──────────────────────────────────────────────────────────────

#' Rank the sweep axes by how much each moves the concordance metric.
#'
#' `results` is a data frame with one row per sweep run: run_id, axis, level, and
#' one or more concordance columns. The axis whose levels produce the widest
#' spread in a metric is the one that explains the divergence in that metric.
#'
#' This is computed, not asserted. The LLM's job afterwards is to name the cause
#' from the ranking plus the discrepancy rubric — not to guess the ranking.
attribute_divergence <- function(results, metrics = NULL, baseline_id = "baseline") {
  if (!"axis" %in% names(results)) {
    stop("`results` needs an 'axis' column naming which axis each run varied.",
         call. = FALSE)
  }
  metrics <- metrics %||% setdiff(names(results),
    c("run_id", "axis", "level", "error", "note"))
  metrics <- metrics[vapply(results[metrics], is.numeric, logical(1))]
  if (!length(metrics)) {
    stop("No numeric concordance columns found in `results`.", call. = FALSE)
  }

  base_row <- results[results$run_id == baseline_id, , drop = FALSE]
  if (!nrow(base_row)) {
    warning("No baseline run found (run_id == '", baseline_id,
            "'); deltas cannot be computed.", call. = FALSE)
  }

  rows <- list()
  for (m in metrics) {
    base_val <- if (nrow(base_row)) base_row[[m]][1] else NA_real_
    for (ax in setdiff(unique(results$axis), "none")) {
      sub <- results[results$axis == ax, , drop = FALSE]
      v <- sub[[m]][is.finite(sub[[m]])]
      if (!length(v)) next
      all_v <- c(v, base_val[is.finite(base_val)])
      rows[[length(rows) + 1]] <- data.frame(
        metric = m, axis = ax, n_levels = length(v),
        min = min(all_v), max = max(all_v),
        spread = max(all_v) - min(all_v),
        baseline = base_val,
        max_abs_delta = if (is.finite(base_val)) max(abs(v - base_val)) else NA_real_,
        direction = if (is.finite(base_val)) {
          d <- v[which.max(abs(v - base_val))] - base_val
          if (d > 0) "increases concordance" else "decreases concordance"
        } else NA_character_,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) stop("Nothing to attribute: no axis had finite metric values.",
                          call. = FALSE)
  tab <- do.call(rbind, rows)
  tab <- tab[order(tab$metric, -tab$spread), ]

  # Per-metric ranking and the dominant axis.
  dominant <- do.call(rbind, lapply(split(tab, tab$metric), function(d) {
    d <- d[order(-d$spread), ]
    data.frame(metric = d$metric[1], dominant_axis = d$axis[1],
               spread = d$spread[1],
               second_axis = if (nrow(d) > 1) d$axis[2] else NA_character_,
               second_spread = if (nrow(d) > 1) d$spread[2] else NA_real_,
               ratio = if (nrow(d) > 1 && d$spread[2] > 0) d$spread[1] / d$spread[2] else NA_real_,
               stringsAsFactors = FALSE)
  }))

  # Is the leading axis actually decisive, or is it a tie? A ratio near 1 means
  # the sweep does not distinguish the axes and the attribution is weak.
  dominant$confidence <- ifelse(is.na(dominant$ratio), "single axis tested",
    ifelse(dominant$ratio >= 2, "clear",
      ifelse(dominant$ratio >= 1.3, "moderate", "weak — axes are not distinguishable")))

  mfg_log("benchmark", "divergence_attributed", list(
    metrics = paste(metrics, collapse = ","),
    n_axes = length(unique(tab$axis)),
    dominant = paste(sprintf("%s:%s", dominant$metric, dominant$dominant_axis),
                     collapse = " "),
    confidence = paste(dominant$confidence, collapse = " | ")))

  out <- list(table = tab, dominant = dominant, metrics = metrics,
              axes = unique(tab$axis))
  class(out) <- c("mfg_attribution", "list")
  out
}

#' @export
print.mfg_attribution <- function(x, ...) {
  cat("=== Divergence attribution ===\n")
  cat("\nDominant axis per metric:\n")
  print(x$dominant, row.names = FALSE)
  cat("\nFull spread by axis:\n")
  print(x$table, row.names = FALSE)
  cat(strwrap(paste("\nSpread is the range a metric takes across that axis's tested",
    "levels, including the baseline. The axis with the widest spread is the one",
    "the divergence is most sensitive to. Where confidence reads 'weak', the",
    "sweep does not separate the axes and no single cause should be named."),
    width = 78), sep = "\n")
  invisible(x)
}

# ── Discrepancy rubric ───────────────────────────────────────────────────────
#
# Known divergence causes paired with what each looks like in the metrics. This
# is what lets the agent go from observing a difference to naming one, instead of
# reasoning it out fresh each time and producing something plausible.

MFG_DISCREPANCY_RUBRIC <- list(
  list(id = "upstream.sample_selection",
       stage = "sequencing/upstream",
       cause = "Different set of runs or samples included",
       signature = paste("Sample counts differ before any filtering; the shared-sample",
         "intersection is much smaller than either table. Profile concordance may be",
         "high on the shared samples while group-level results differ."),
       check = "Compare n per group in the study card against the reanalysis intake log."),

  list(id = "upstream.pretrimmed_reads",
       stage = "sequencing/upstream",
       cause = "Reads were already primer-trimmed or quality-trimmed differently",
       signature = paste("Read depths differ systematically by a roughly constant factor;",
         "richness shifts in one direction across all samples."),
       check = "Compare the depth distribution against any depths the paper reports."),

  list(id = "denoise.otu_vs_asv",
       stage = "denoising",
       cause = "97% OTU clustering versus exact ASVs",
       signature = paste("Feature counts differ several-fold (OTUs are fewer, ASVs more).",
         "Richness and Chao1 differ strongly; Shannon much less; beta-diversity",
         "structure is largely preserved (high Mantel rho) because the abundant",
         "community is the same."),
       check = "The denoiser field in the study card; feature count ratio."),

  list(id = "taxonomy.database_version",
       stage = "taxonomy",
       cause = "Different reference database release",
       signature = paste("The signature cause of 'disappearing' taxa. Many published",
         "names fail to match, and the unmatched list is dominated by renamed clades",
         "(Lactobacillus splits, phylum -ota renamings). Recovery rate is low while",
         "effect-size correlation on the matched subset stays high."),
       check = paste("The unresolvable/unmatched label list from match_taxa(). If those",
         "names appear in MFG_TAXON_SYNONYMS, this is the cause.")),

  list(id = "taxonomy.confidence_threshold",
       stage = "taxonomy",
       cause = "Different classifier confidence threshold",
       signature = paste("The share of reads unassigned at genus or species differs",
         "markedly; agglomerating to genus loses a different fraction on each side."),
       check = "reads_unresolved_pct from agglomerate_and_summarise() on both sides."),

  list(id = "normalization.rarefaction_depth",
       stage = "normalization",
       cause = "Different rarefaction depth, or rarefaction versus none",
       signature = paste("Observed richness and Chao1 differ with a systematic bias in",
         "Bland-Altman (offset scaling with depth), while Shannon and beta-diversity",
         "barely move. Sample counts may also differ if one depth dropped samples."),
       check = "Bland-Altman bias and proportional-bias slope for richness; the sweep axis rarefaction_depth."),

  list(id = "normalization.method",
       stage = "normalization",
       cause = "Compositional versus count-based normalization (CLR/TSS vs raw/DESeq2)",
       signature = paste("Differential abundance lists diverge substantially while alpha",
         "and beta agree. Effect-size sign agreement stays high but magnitudes are on",
         "different scales."),
       check = "The normalization sweep axis; whether the paper's method accepts what it was given."),

  list(id = "filtering.prevalence",
       stage = "filtering",
       cause = "Different prevalence or abundance filter",
       signature = paste("The multiple-testing denominator differs, so borderline taxa",
         "move in and out of significance. Missed taxa turn out never to have been",
         "tested rather than tested and not called."),
       check = "concordance_da()$missed_because_not_tested; the prevalence_filter sweep axis."),

  list(id = "statistics.da_method",
       stage = "statistics",
       cause = "Different differential abundance method",
       signature = paste("Recovery rate is moderate while effect-size rank correlation is",
         "high — the methods agree about the direction and size of effects and",
         "disagree about where to draw the line. LEfSe and plain Wilcoxon call more",
         "taxa than ANCOM-BC2 on the same counts."),
       check = "The da_method sweep axis; compare da_compare() consensus against each method's unique set."),

  list(id = "statistics.covariates",
       stage = "statistics",
       cause = "Covariates included in one analysis and not the other",
       signature = paste("PERMANOVA R2 for the main term differs substantially while the",
         "distance matrices agree (high Mantel rho). Some DA hits vanish when a",
         "covariate absorbs their variance."),
       check = "The covariates field in the study card against the reanalysis model formula."),

  list(id = "statistics.repeated_measures",
       stage = "statistics",
       cause = "Repeated measures treated as independent in one analysis",
       signature = paste("The independent-samples analysis reports many more significant",
         "results. The p-value histogram is anti-conservative. Adding a subject",
         "random effect removes most of the difference."),
       check = "mfg_detect_repeated_measures() on the metadata; the random_effects sweep axis."),

  list(id = "statistics.fdr_method",
       stage = "statistics",
       cause = "Different multiple-testing correction",
       signature = paste("Raw p-values agree closely; adjusted values differ. Holm is much",
         "more conservative than BH, so hit counts can differ severalfold with",
         "identical underlying statistics."),
       check = "The fdr_method sweep axis; compare raw p-value correlation against adjusted."),

  list(id = "metadata.recoded_groups",
       stage = "metadata",
       cause = "Groups recoded, merged, or a reference level changed",
       signature = paste("Effect signs flip wholesale while magnitudes are preserved. The",
         "direction_flip flag fires in concordance_alpha while correlation stays high."),
       check = "Group level names and order; concordance_alpha()$group_test$direction_flip."),

  list(id = "metadata.outlier_exclusion",
       stage = "metadata",
       cause = "Outliers excluded in one analysis without being reported",
       signature = paste("Sample counts differ by a handful with no stated rule. Removing",
         "the same few samples from the reanalysis reconciles the results."),
       check = "n per group against the paper's stated n; the QC exclusion log."),

  list(id = "stochastic.unseeded",
       stage = "stochasticity",
       cause = "Unseeded subsampling or too few permutations",
       signature = paste("Results differ slightly and irreproducibly between runs of the",
         "same configuration. p-values near the threshold flip; effect sizes barely",
         "move. Re-running the identical configuration gives a different answer."),
       check = paste("Re-run the baseline twice and compare. If the two baselines differ,",
         "the divergence is noise and no axis will explain it."))
)

#' Retrieve rubric entries relevant to a stage or a metric pattern.
discrepancy_rubric <- function(stage = NULL, id = NULL, search = NULL) {
  r <- MFG_DISCREPANCY_RUBRIC
  if (!is.null(id))    r <- Filter(function(e) e$id %in% id, r)
  if (!is.null(stage)) r <- Filter(function(e) grepl(stage, e$stage, ignore.case = TRUE), r)
  if (!is.null(search)) {
    r <- Filter(function(e) any(grepl(search, c(e$cause, e$signature, e$check),
                                      ignore.case = TRUE)), r)
  }
  r
}

#' The rubric as a printable table, for the report appendix.
rubric_table <- function(entries = MFG_DISCREPANCY_RUBRIC) {
  do.call(rbind, lapply(entries, function(e) data.frame(
    id = e$id, stage = e$stage, cause = e$cause,
    signature = e$signature, check = e$check, stringsAsFactors = FALSE)))
}

# ── The attribution record ───────────────────────────────────────────────────

#' Record an attribution with its evidence and its unexplained residual.
#'
#' `unexplained_residual` is the field that stops a plausible cause being
#' invented when none of the tested axes explain the gap. "No tested factor
#' accounts for this" is a valid and valuable answer, and the sweep is the
#' evidence for saying it.
attribution_record <- function(verdict, causes = list(),
                               unexplained_residual = NULL,
                               confidence = c("high", "medium", "low"),
                               attribution = NULL, sweep_grid = NULL) {
  confidence <- match.arg(confidence)
  valid_verdicts <- c("reproduced", "partially_reproduced", "diverged")
  if (!verdict %in% valid_verdicts) {
    stop("verdict must be one of: ", paste(valid_verdicts, collapse = ", "),
         call. = FALSE)
  }
  if (!identical(verdict, "reproduced") && is.null(unexplained_residual)) {
    stop("A non-reproducing verdict requires unexplained_residual — state what the ",
         "sweep does NOT account for, even if that is 'nothing, every gap is ",
         "attributed'. Leaving it empty is how a fabricated cause gets through.",
         call. = FALSE)
  }
  for (i in seq_along(causes)) {
    c_ <- causes[[i]]
    miss <- setdiff(c("rubric_id", "axis", "evidence"), names(c_))
    if (length(miss)) {
      stop("cause ", i, " is missing: ", paste(miss, collapse = ", "),
           "\nEvery cause must name a rubric entry, the sweep axis that showed it, ",
           "and the evidence (run ids and the metric movement).", call. = FALSE)
    }
    if (!c_$rubric_id %in% vapply(MFG_DISCREPANCY_RUBRIC, function(e) e$id, character(1))) {
      warning("cause ", i, " cites rubric_id '", c_$rubric_id,
              "' which is not in the rubric.", call. = FALSE)
    }
  }

  untested <- character()
  if (!is.null(sweep_grid) && sweep_grid$n_truncated > 0) {
    untested <- sprintf("%d axis levels were not tested because of the run cap",
                        sweep_grid$n_truncated)
  }

  out <- list(verdict = verdict, causes = causes,
              unexplained_residual = unexplained_residual,
              confidence = confidence, attribution = attribution,
              untested = untested)
  mfg_log("benchmark", "attribution_recorded", list(
    verdict = verdict, n_causes = length(causes), confidence = confidence,
    rubric_ids = paste(vapply(causes, function(c_) c_$rubric_id, character(1)),
                       collapse = ","),
    has_residual = !is.null(unexplained_residual)))
  class(out) <- c("mfg_attribution_record", "list")
  out
}

#' @export
print.mfg_attribution_record <- function(x, ...) {
  cat("=== Attribution ===\n")
  cat(sprintf("Verdict: %s (confidence: %s)\n", toupper(gsub("_", " ", x$verdict)),
              x$confidence))
  if (length(x$causes)) {
    cat("\nCauses:\n")
    for (c_ in x$causes) {
      entry <- discrepancy_rubric(id = c_$rubric_id)
      cat(sprintf("\n  [%s] %s\n", c_$rubric_id,
                  if (length(entry)) entry[[1]]$cause else "(not in rubric)"))
      cat(sprintf("    axis: %s\n", c_$axis))
      cat(strwrap(sprintf("evidence: %s", c_$evidence), width = 74,
                  prefix = "    "), sep = "\n")
      if (!is.null(c_$explanation)) {
        cat(strwrap(c_$explanation, width = 74, prefix = "    "), sep = "\n")
      }
    }
  } else {
    cat("\nNo causes attributed.\n")
  }
  if (!is.null(x$unexplained_residual)) {
    cat("\nUnexplained residual:\n")
    cat(strwrap(x$unexplained_residual, width = 76, prefix = "  "), sep = "\n")
  }
  if (length(x$untested)) {
    cat("\nCoverage limits:\n")
    for (u in x$untested) cat("  - ", u, "\n", sep = "")
  }
  invisible(x)
}
