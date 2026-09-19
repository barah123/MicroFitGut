# =============================================================================
# MicroFitGut — 02-qc.R
#
# Filtering, and an auditable record of everything it removed.
#
# The governing rule: no sample and no taxon leaves the dataset without an entry
# in the QC log saying which rule dropped it and why. reference/10 requires the
# report to state every exclusion, and this is where that record comes from.
#
# Order matters. Unwanted taxa first (they are not data), then prevalence, then
# singletons, then depth — because dropping taxa changes sample depths, and
# pruning on a depth computed before taxon removal drops the wrong samples.
#
# Requires: mfg_require("core") and utils.R
# =============================================================================

# ── The QC record ────────────────────────────────────────────────────────────

#' Start a QC record against a starting object.
qc_begin <- function(ps) {
  structure(list(
    start = list(taxa = phyloseq::ntaxa(ps), samples = phyloseq::nsamples(ps),
                 reads = sum(phyloseq::sample_sums(ps))),
    steps = list()
  ), class = c("mfg_qc_log", "list"))
}

#' Record one filtering step by comparing before and after.
qc_step <- function(log, step, before, after, rule, rationale) {
  dropped_taxa    <- setdiff(phyloseq::taxa_names(before),   phyloseq::taxa_names(after))
  dropped_samples <- setdiff(phyloseq::sample_names(before), phyloseq::sample_names(after))
  entry <- list(
    step = step, rule = rule, rationale = rationale,
    taxa_before = phyloseq::ntaxa(before), taxa_after = phyloseq::ntaxa(after),
    samples_before = phyloseq::nsamples(before), samples_after = phyloseq::nsamples(after),
    reads_before = sum(phyloseq::sample_sums(before)),
    reads_after  = sum(phyloseq::sample_sums(after)),
    dropped_taxa = dropped_taxa, dropped_samples = dropped_samples
  )
  log$steps[[length(log$steps) + 1]] <- entry
  mfg_log("qc", step, list(
    rule = rule, taxa = sprintf("%d->%d", entry$taxa_before, entry$taxa_after),
    samples = sprintf("%d->%d", entry$samples_before, entry$samples_after),
    n_dropped_samples = length(dropped_samples), n_dropped_taxa = length(dropped_taxa)))
  log
}

#' @export
print.mfg_qc_log <- function(x, ...) {
  cat("=== QC log ===\n")
  cat(sprintf("Start: %d taxa, %d samples, %s reads\n",
              x$start$taxa, x$start$samples, format(x$start$reads, big.mark = ",")))
  for (s in x$steps) {
    cat(sprintf("\n%-22s %s\n", s$step, s$rule))
    cat(sprintf("  taxa %d -> %d (-%d)   samples %d -> %d (-%d)\n",
                s$taxa_before, s$taxa_after, s$taxa_before - s$taxa_after,
                s$samples_before, s$samples_after, s$samples_before - s$samples_after))
    if (length(s$dropped_samples)) {
      cat(sprintf("  dropped samples: %s\n",
                  paste(utils::head(s$dropped_samples, 12), collapse = ", ")))
      if (length(s$dropped_samples) > 12)
        cat(sprintf("                   ... and %d more\n", length(s$dropped_samples) - 12))
    }
    cat(sprintf("  why: %s\n", s$rationale))
  }
  if (length(x$steps)) {
    last <- x$steps[[length(x$steps)]]
    cat(sprintf("\nEnd: %d taxa, %d samples, %s reads (%.1f%% of reads retained)\n",
                last$taxa_after, last$samples_after,
                format(last$reads_after, big.mark = ","),
                100 * last$reads_after / x$start$reads))
  }
  invisible(x)
}

#' Flatten the QC log to a data frame the report can print directly.
qc_table <- function(log) {
  if (!length(log$steps)) return(data.frame())
  do.call(rbind, lapply(log$steps, function(s) data.frame(
    step = s$step, rule = s$rule,
    taxa_before = s$taxa_before, taxa_after = s$taxa_after,
    samples_before = s$samples_before, samples_after = s$samples_after,
    reads_retained_pct = round(100 * s$reads_after / s$reads_before, 2),
    dropped_samples = paste(s$dropped_samples, collapse = ";"),
    rationale = s$rationale, stringsAsFactors = FALSE)))
}

#' Every sample ever dropped, with the rule that dropped it.
#'
#' reference/10 requires this in the report verbatim.
qc_exclusions <- function(log) {
  rows <- list()
  for (s in log$steps) {
    for (smp in s$dropped_samples) {
      rows[[length(rows) + 1]] <- data.frame(
        sample = smp, dropped_at = s$step, rule = s$rule, stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(data.frame(sample = character(), dropped_at = character(),
                                       rule = character()))
  do.call(rbind, rows)
}

# ── Individual filters ───────────────────────────────────────────────────────

#' Remove taxa that are not data.
#'
#' Mitochondria and chloroplasts are host or plant organellar 16S that amplifies
#' with bacterial primers; they are contamination in a bacterial census, not
#' low-abundance findings. Unassigned kingdom-level calls are usually chimeras or
#' off-target amplification. Leaving them in inflates richness and wastes
#' multiple-testing budget on sequences that cannot be interpreted.
remove_unwanted_taxa <- function(ps,
    drop_family = c("mitochondria", "Mitochondria"),
    drop_order  = c("Chloroplast", "chloroplast"),
    drop_kingdom_na = TRUE,
    drop_kingdom = c("Eukaryota", "Unassigned", "unassigned", "NA"),
    verbose = TRUE) {

  if (is.null(phyloseq::tax_table(ps, errorIfNULL = FALSE))) {
    if (verbose) message("No taxonomy table; nothing to filter by name.")
    return(ps)
  }
  tt <- mfg_tax(ps)
  keep <- rep(TRUE, nrow(tt))
  hits <- list()

  if ("Family" %in% names(tt) && length(drop_family)) {
    bad <- tt$Family %in% drop_family
    hits$Family <- sum(bad, na.rm = TRUE); keep <- keep & !(bad %in% TRUE)
  }
  if ("Order" %in% names(tt) && length(drop_order)) {
    bad <- tt$Order %in% drop_order
    hits$Order <- sum(bad, na.rm = TRUE); keep <- keep & !(bad %in% TRUE)
  }
  if ("Kingdom" %in% names(tt)) {
    bad <- tt$Kingdom %in% drop_kingdom
    if (drop_kingdom_na) bad <- bad | is.na(tt$Kingdom)
    hits$Kingdom <- sum(bad, na.rm = TRUE); keep <- keep & !(bad %in% TRUE)
  }
  if (verbose && sum(unlist(hits)) > 0) {
    message(sprintf("Removing %d unwanted taxa (%s)", sum(unlist(hits)),
                    paste(sprintf("%s:%d", names(hits), unlist(hits)), collapse = " ")))
  }
  mfg_carry_attrs(phyloseq::prune_taxa(keep, ps), ps)
}

#' Remove singletons: taxa with one read in total across the whole study.
#'
#' A taxon seen once cannot be distinguished from a sequencing error, and it
#' carries a full multiple-testing penalty for no possible signal. Threshold is
#' exposed because "singleton" means >1 read in the course material but some
#' pipelines use >2 or a per-sample rule.
remove_singletons <- function(ps, min_total_reads = 1) {
  mfg_carry_attrs(phyloseq::prune_taxa(phyloseq::taxa_sums(ps) > min_total_reads, ps), ps)
}

#' Prevalence filter: keep taxa present in at least `min_prevalence` of samples.
#'
#' `min_prevalence` below 1 is a fraction, 1 or above is a sample count. This is
#' the filter that most changes a differential abundance result, so the threshold
#' must be stated in the report and never tuned after seeing the p-values —
#' see reference/11.
filter_prevalence <- function(ps, min_prevalence = 0.10, detection = 0) {
  mat <- as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) mat <- t(mat)
  prev_n <- rowSums(mat > detection)
  n_samp <- ncol(mat)
  thresh <- if (min_prevalence < 1) ceiling(min_prevalence * n_samp) else min_prevalence
  mfg_carry_attrs(phyloseq::prune_taxa(prev_n >= thresh, ps), ps)
}

#' Drop HUMAnN bookkeeping rows and species-stratified rows.
#'
#' A HUMAnN pathway or gene family table is not a flat feature table. It carries
#' two sentinel rows, UNMAPPED and UNINTEGRATED, which count reads that did not
#' map or did not fall in a known pathway, and it repeats every community-level
#' feature once per contributing species, written FEATURE|g__Genus.s__species.
#'
#' Both have to go before testing, for different reasons. The sentinels are
#' quality measures rather than biology, and they are usually large enough to
#' dominate any compositional transform. The stratified rows are not independent
#' of the community rows they decompose, so leaving them in tests the same
#' signal many times over and inflates the multiple-testing denominator with
#' rows that cannot fail independently.
#'
#' Returns the pruned object with a summary attached as the "mfg_humann_filter"
#' attribute, so the counts can be reported rather than asserted.
filter_humann <- function(ps, drop_sentinels = TRUE, drop_stratified = TRUE,
                          sentinels = c("UNMAPPED", "UNINTEGRATED"),
                          verbose = TRUE) {
  feats <- phyloseq::taxa_names(ps)
  n0 <- length(feats)

  is_sentinel <- rep(FALSE, n0)
  if (isTRUE(drop_sentinels)) {
    pat <- paste0("^(", paste(sentinels, collapse = "|"), ")")
    is_sentinel <- grepl(pat, feats, ignore.case = TRUE)
  }
  # The pipe separates a feature from the species it was attributed to. It does
  # not otherwise occur in MetaCyc or UniRef identifiers.
  is_strat <- rep(FALSE, n0)
  if (isTRUE(drop_stratified)) is_strat <- grepl("\\|", feats, fixed = FALSE)

  keep <- !is_sentinel & !is_strat
  if (!any(keep)) {
    stop("Every feature was removed. Check that this is a HUMAnN table and that ",
         "the rows are features, not samples.", call. = FALSE)
  }

  summary_tbl <- data.frame(
    step = c("input", "sentinel rows removed", "stratified rows removed", "retained"),
    n    = c(n0, sum(is_sentinel), sum(is_strat & !is_sentinel), sum(keep)),
    stringsAsFactors = FALSE)

  out <- mfg_carry_attrs(phyloseq::prune_taxa(keep, ps), ps)
  attr(out, "mfg_humann_filter") <- summary_tbl

  mfg_log("qc", "filter_humann", list(
    n_input = n0,
    n_sentinel_removed = sum(is_sentinel),
    sentinels_found = paste(feats[is_sentinel], collapse = ", "),
    n_stratified_removed = sum(is_strat & !is_sentinel),
    n_retained = sum(keep)))

  if (isTRUE(verbose)) {
    message(sprintf(
      "filter_humann: %d features in, %d sentinel and %d stratified removed, %d retained.",
      n0, sum(is_sentinel), sum(is_strat & !is_sentinel), sum(keep)))
  }
  out
}

#' Abundance filter: keep taxa reaching a minimum total or mean relative share.
filter_abundance <- function(ps, min_total_reads = NULL, min_mean_relative = NULL) {
  keep <- rep(TRUE, phyloseq::ntaxa(ps))
  if (!is.null(min_total_reads)) {
    keep <- keep & (phyloseq::taxa_sums(ps) >= min_total_reads)
  }
  if (!is.null(min_mean_relative)) {
    rel <- phyloseq::transform_sample_counts(ps, function(x) if (sum(x) > 0) x / sum(x) else x)
    m <- as(phyloseq::otu_table(rel), "matrix")
    if (!phyloseq::taxa_are_rows(rel)) m <- t(m)
    keep <- keep & (rowMeans(m) >= min_mean_relative)
  }
  mfg_carry_attrs(phyloseq::prune_taxa(keep, ps), ps)
}

#' Drop samples below a sequencing depth.
#'
#' A shallow sample does not have a small amount of the community; it has a
#' biased sample of it, weighted toward abundant taxa. Its richness is an
#' artefact of depth. The threshold is a judgement call and must be justified in
#' the report, so it is required rather than defaulted silently — the course
#' material uses 5000 (Demo 8, Quiz 8) and 20000 (Demo 6) for different studies.
prune_low_depth <- function(ps, min_depth, justification = NULL) {
  if (missing(min_depth)) {
    stop("min_depth must be given explicitly. There is no universal threshold; ",
         "choose one from the depth distribution and the rarefaction curve, and ",
         "record why. See reference/02.", call. = FALSE)
  }
  if (is.null(justification)) {
    warning("prune_low_depth called without a justification. reference/10 ",
            "requires the report to state why the threshold was chosen.",
            call. = FALSE)
  }
  mfg_carry_attrs(phyloseq::prune_samples(phyloseq::sample_sums(ps) >= min_depth, ps), ps)
}

#' Drop taxa that are now empty after sample removal.
#'
#' Pruning samples can leave taxa with zero reads everywhere. They inflate the
#' multiple-testing count and appear in plots as flat lines, so they go.
drop_empty_taxa <- function(ps) {
  mfg_carry_attrs(phyloseq::prune_taxa(phyloseq::taxa_sums(ps) > 0, ps), ps)
}

# ── The standard pipeline ────────────────────────────────────────────────────

#' Run the QC sequence and return both the filtered object and the record.
#'
#' Set any threshold to NULL to skip that step. `min_depth` has no default on
#' purpose: choosing it is an analytical decision, not a default.
#'
#' Returns list(ps = <filtered>, log = <mfg_qc_log>).
run_qc <- function(ps,
                   humann            = NULL,
                   remove_unwanted   = TRUE,
                   prevalence        = 0.10,
                   singleton_min     = 1,
                   min_depth         = NULL,
                   depth_justification = NULL,
                   drop_empty        = TRUE,
                   verbose           = TRUE) {

  is_rel <- attr(ps, "mfg_is_relative") %||% mfg_registry_get(ps, "is_relative") %||%
    looks_like_relative_abundance(as(phyloseq::otu_table(ps), "matrix"))

  log <- qc_begin(ps)

  # A HUMAnN table carries sentinel rows and species-stratified rows that are
  # not community features. They go first: the sentinels are large enough to
  # distort anything abundance-based that follows, and the stratified rows would
  # otherwise be counted as independent features by every later step.
  if (is.null(humann)) {
    fn <- phyloseq::taxa_names(ps)
    humann <- any(grepl("^(UNMAPPED|UNINTEGRATED)", fn, ignore.case = TRUE)) ||
      any(grepl("\\|", fn))
  }
  if (isTRUE(humann)) {
    before <- ps
    ps <- filter_humann(ps, verbose = verbose)
    log <- qc_step(log, "humann_features", before, ps,
      rule = "drop UNMAPPED and UNINTEGRATED rows, and species-stratified rows (FEATURE|g__...s__...)",
      rationale = paste("The sentinels count reads that did not map or did not",
                        "fall in a known pathway; they are quality measures rather",
                        "than community features and are large enough to dominate",
                        "a compositional transform. The stratified rows decompose",
                        "the community rows they accompany, so testing both counts",
                        "the same signal repeatedly and pads the multiple-testing",
                        "denominator with rows that cannot fail independently."))
  }

  # The taxonomic filter needs a tax_table. A functional table has none, so
  # asking for it there is a category error rather than a choice.
  if (isTRUE(remove_unwanted) && is.null(phyloseq::access(ps, "tax_table"))) {
    if (verbose) message("Skipping unwanted-taxa removal: no tax_table (functional table).")
    mfg_log("qc", "unwanted_taxa_skipped",
            list(reason = "no tax_table; feature ids are pathways or gene families"))
    remove_unwanted <- FALSE
  }

  if (isTRUE(remove_unwanted)) {
    before <- ps
    ps <- remove_unwanted_taxa(ps, verbose = verbose)
    log <- qc_step(log, "unwanted_taxa", before, ps,
      rule = "drop Family=mitochondria, Order=Chloroplast, Kingdom unassigned/NA/Eukaryota",
      rationale = paste("Organellar and unassigned sequences are off-target",
                        "amplification, not community members. Retaining them",
                        "inflates richness and spends multiple-testing budget on",
                        "sequences that cannot be interpreted."))
  }

  if (!is.null(prevalence)) {
    before <- ps
    ps <- filter_prevalence(ps, min_prevalence = prevalence)
    log <- qc_step(log, "prevalence", before, ps,
      rule = sprintf("keep taxa present in >= %s of samples",
                     if (prevalence < 1) paste0(100 * prevalence, "%") else prevalence),
      rationale = paste("A taxon seen in a handful of samples cannot support a",
                        "group-level inference and costs a full multiple-testing",
                        "penalty. Threshold fixed before any test was run."))
  }

  if (!is.null(singleton_min) && !isTRUE(is_rel)) {
    before <- ps
    ps <- remove_singletons(ps, min_total_reads = singleton_min)
    log <- qc_step(log, "singletons", before, ps,
      rule = sprintf("drop taxa with <= %d total reads across all samples", singleton_min),
      rationale = paste("A taxon seen once study-wide is indistinguishable from",
                        "a sequencing error."))
  } else if (!is.null(singleton_min) && isTRUE(is_rel)) {
    if (verbose) message("Skipping singleton removal: values are relative abundances, not counts.")
    mfg_log("qc", "singletons_skipped",
            list(reason = "input is relative abundance; no read counts to threshold"))
  }

  if (!is.null(min_depth)) {
    if (isTRUE(is_rel)) {
      if (verbose) message("Skipping depth pruning: values are relative abundances.")
      mfg_log("qc", "depth_prune_skipped",
              list(reason = "input is relative abundance; sample sums carry no depth information"))
    } else {
      before <- ps
      ps <- prune_low_depth(ps, min_depth, justification = depth_justification)
      log <- qc_step(log, "low_depth_samples", before, ps,
        rule = sprintf("drop samples with < %s reads", format(min_depth, big.mark = ",")),
        rationale = depth_justification %||% paste(
          "Shallow samples give a depth-biased view of the community, so their",
          "richness is an artefact. NO JUSTIFICATION RECORDED — supply one."))
    }
  }

  if (isTRUE(drop_empty)) {
    before <- ps
    ps <- drop_empty_taxa(ps)
    if (phyloseq::ntaxa(before) != phyloseq::ntaxa(ps)) {
      log <- qc_step(log, "empty_taxa", before, ps,
        rule = "drop taxa with zero reads after sample filtering",
        rationale = "These became empty when samples were removed; they carry no information.")
    }
  }

  if (verbose) print(log)
  list(ps = ps, log = log)
}

#' Prevalence-vs-abundance diagnostic, for choosing a prevalence threshold.
#'
#' Plot this before fixing the threshold, not after — the shape of the cloud is
#' what justifies the cut, and picking it from the p-values is the pitfall in
#' reference/11.
qc_prevalence_table <- function(ps, rank = "Phylum") {
  mat <- as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) mat <- t(mat)
  df <- data.frame(
    taxon        = rownames(mat),
    prevalence   = rowSums(mat > 0),
    prevalence_frac = rowSums(mat > 0) / ncol(mat),
    total_abundance = rowSums(mat),
    mean_when_present = ifelse(rowSums(mat > 0) > 0,
                               rowSums(mat) / pmax(rowSums(mat > 0), 1), 0),
    stringsAsFactors = FALSE)
  r <- safe_rank(ps, rank)
  if (!is.null(r) && !is.null(phyloseq::tax_table(ps, errorIfNULL = FALSE))) {
    tt <- mfg_tax(ps)
    df[[r]] <- tt[df$taxon, r]
  }
  df[order(-df$prevalence), ]
}
