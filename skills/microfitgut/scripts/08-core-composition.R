# =============================================================================
# MicroFitGut — 08-core-composition.R
#
# Who is there, at what rank, and which taxa are shared across samples.
#
# Composition is where the two thresholds that define "core" have to be stated,
# because "core microbiome" means nothing without them: a detection threshold
# (how abundant before it counts as present) and a prevalence threshold (in what
# fraction of samples). Demo 6 uses detection = 0.2 and prevalence = 0.9, which
# is strict; much of the literature uses detection near zero and prevalence 0.5.
# Both are defensible, neither is default, and a core set quoted without them is
# not reproducible.
#
# Requires: mfg_require("core"); utils.R; 03-normalize.R
# =============================================================================

# ── Agglomeration ────────────────────────────────────────────────────────────

#' Agglomerate to a rank and summarise what happened.
#'
#' tax_glom drops taxa with NA at the target rank, which silently loses reads.
#' This reports how many and what share, because "we analysed at genus level"
#' often means "we discarded the 30% of reads that were unclassified at genus".
agglomerate_and_summarise <- function(ps, rank = "Genus", NArm = TRUE) {
  r <- safe_rank(ps, rank)
  if (is.null(r)) {
    stop("Rank '", rank, "' is not present. Available: ",
         paste(phyloseq::rank_names(ps), collapse = ", "), call. = FALSE)
  }
  reads_before <- sum(phyloseq::sample_sums(ps))
  taxa_before  <- phyloseq::ntaxa(ps)

  tt <- mfg_tax(ps)
  n_na <- sum(is.na(tt[[r]]) | tt[[r]] %in% c("", "unclassified", "Unclassified", "NA"))
  reads_na <- if (n_na > 0) {
    na_taxa <- rownames(tt)[is.na(tt[[r]]) |
      tt[[r]] %in% c("", "unclassified", "Unclassified", "NA")]
    sum(phyloseq::taxa_sums(ps)[na_taxa])
  } else 0

  out <- mfg_carry_attrs(phyloseq::tax_glom(ps, taxrank = r, NArm = NArm), ps)
  out <- mfg_set_normalization(out, mfg_normalization(ps))

  reads_after <- sum(phyloseq::sample_sums(out))
  pct_lost <- 100 * (reads_before - reads_after) / reads_before

  mfg_log("composition", "agglomerated", list(
    rank = r, NArm = NArm,
    taxa = sprintf("%d->%d", taxa_before, phyloseq::ntaxa(out)),
    unresolved_at_rank = n_na,
    reads_lost_pct = round(pct_lost, 2),
    reads_unresolved_pct = round(100 * reads_na / reads_before, 2)))

  if (pct_lost > 5) {
    warning(sprintf(paste("Agglomerating to %s discarded %.1f%% of reads (%d taxa",
      "unresolved at this rank). State this in the report — a genus-level result",
      "that silently drops a tenth of the data is not a genus-level result."),
      r, pct_lost, n_na), call. = FALSE)
  }

  structure(list(ps = out, rank = r,
    taxa_before = taxa_before, taxa_after = phyloseq::ntaxa(out),
    n_unresolved = n_na, reads_lost_pct = pct_lost,
    reads_unresolved_pct = 100 * reads_na / reads_before),
    class = c("mfg_agglomeration", "list"))
}

#' @export
print.mfg_agglomeration <- function(x, ...) {
  cat(sprintf("=== Agglomerated to %s ===\n", x$rank))
  cat(sprintf("Taxa: %d -> %d\n", x$taxa_before, x$taxa_after))
  cat(sprintf("Unresolved at %s: %d taxa (%.1f%% of reads)\n", x$rank,
              x$n_unresolved, x$reads_unresolved_pct))
  cat(sprintf("Reads lost: %.2f%%\n", x$reads_lost_pct))
  invisible(x)
}

# ── Relative abundance summaries ─────────────────────────────────────────────

#' Mean relative abundance per taxon, overall and by group.
#'
#' Reports both mean and prevalence, because they answer different questions and
#' are routinely confused: a taxon can be highly abundant in three samples and
#' absent elsewhere (high mean, low prevalence) or present everywhere at 0.1%
#' (low mean, high prevalence).
taxa_summary <- function(ps, rank = NULL, group_var = NULL, top_n = NULL) {
  work <- if (!is.null(rank)) agglomerate_and_summarise(ps, rank)$ps else ps
  rel  <- if (identical(mfg_normalization(work), "tss")) work else tss_transform(work)
  mat  <- mfg_otu_taxa_as_rows(rel)

  df <- data.frame(
    taxon = rownames(mat),
    mean_relative = rowMeans(mat),
    median_relative = apply(mat, 1, stats::median),
    max_relative = apply(mat, 1, max),
    prevalence = rowSums(mat > 0) / ncol(mat),
    mean_when_present = apply(mat, 1, function(r) if (any(r > 0)) mean(r[r > 0]) else 0),
    stringsAsFactors = FALSE)

  tt <- mfg_tax(work)
  if (!is.null(tt)) {
    for (r in intersect(c("Phylum", "Class", "Order", "Family", "Genus", "Species",
                          "Pathways"), names(tt))) {
      df[[r]] <- tt[df$taxon, r]
    }
  }

  if (!is.null(group_var)) {
    meta <- mfg_meta(rel)
    g <- as.factor(meta[[group_var]])
    for (lv in levels(g)) {
      cols <- rownames(meta)[which(g == lv)]
      cols <- intersect(cols, colnames(mat))
      df[[paste0("mean_", lv)]] <- rowMeans(mat[, cols, drop = FALSE])
      df[[paste0("prev_", lv)]] <- rowSums(mat[, cols, drop = FALSE] > 0) / length(cols)
    }
  }

  df <- df[order(-df$mean_relative), ]
  if (!is.null(top_n)) df <- utils::head(df, top_n)
  mfg_log("composition", "taxa_summarised",
          list(rank = rank %||% "as-is", n_taxa = nrow(df),
               group = group_var %||% "none",
               top_taxon = df$taxon[1],
               top_mean_pct = round(100 * df$mean_relative[1], 2)))
  df
}

#' Collapse everything below the top N into "Other", for readable bar plots.
#'
#' A stacked bar with 140 taxa communicates nothing. This makes the collapse
#' explicit and records how much abundance went into "Other" — which the plot
#' itself cannot show.
collapse_to_top <- function(ps, rank = "Phylum", top_n = 10, other_label = "Other") {
  work <- agglomerate_and_summarise(ps, rank)$ps
  rel  <- if (identical(mfg_normalization(work), "tss")) work else tss_transform(work)
  mat  <- mfg_otu_taxa_as_rows(rel)
  tt   <- mfg_tax(rel)

  means <- sort(rowMeans(mat), decreasing = TRUE)
  keep  <- names(means)[seq_len(min(top_n, length(means)))]
  other <- setdiff(rownames(mat), keep)

  new_mat <- mat[keep, , drop = FALSE]
  if (length(other)) {
    new_mat <- rbind(new_mat, colSums(mat[other, , drop = FALSE]))
    rownames(new_mat)[nrow(new_mat)] <- other_label
  }
  labels <- c(tt[keep, rank], if (length(other)) other_label else NULL)
  rownames(new_mat) <- make.unique(as.character(labels))

  pct_other <- if (length(other)) 100 * mean(colSums(mat[other, , drop = FALSE])) else 0
  mfg_log("composition", "collapsed_to_top", list(
    rank = rank, top_n = top_n, n_collapsed = length(other),
    mean_pct_in_other = round(pct_other, 2)))

  list(matrix = new_mat, kept = keep, collapsed = other,
       mean_pct_other = pct_other, rank = rank)
}

# ── Core microbiome ──────────────────────────────────────────────────────────

#' Core taxa at stated detection and prevalence thresholds.
#'
#' Both thresholds are required arguments in spirit — they have defaults matching
#' Demo 6, but the return value carries them and the report must quote them,
#' because "the core microbiome" is not a property of the data alone.
#'
#' detection is on the relative-abundance scale, so 0.2 means 20% of a sample.
#' That is a strict threshold; 0.001 (0.1%) with prevalence 0.5 is the commoner
#' convention and finds a much larger core.
core_taxa <- function(ps, detection = 0.2, prevalence = 0.9, rank = NULL) {
  work <- if (!is.null(rank)) agglomerate_and_summarise(ps, rank)$ps else ps
  rel  <- if (identical(mfg_normalization(work), "tss")) work else tss_transform(work)

  ps_core <- microbiome::core(rel, detection = detection, prevalence = prevalence)
  taxa <- phyloseq::taxa_names(ps_core)

  tt <- mfg_tax(rel)
  tax_df <- if (!is.null(tt) && length(taxa)) {
    d <- tt[taxa, , drop = FALSE]; d$taxon <- rownames(d); d
  } else data.frame(taxon = taxa, stringsAsFactors = FALSE)

  # What share of each sample the core accounts for. A "core" holding 5% of the
  # community is a different claim from one holding 80%.
  mat <- mfg_otu_taxa_as_rows(rel)
  core_share <- if (length(taxa)) colSums(mat[taxa, , drop = FALSE]) else rep(0, ncol(mat))

  mfg_log("composition", "core_taxa", list(
    detection = detection, prevalence = prevalence,
    rank = rank %||% "as-is", n_core = length(taxa),
    n_total = phyloseq::ntaxa(rel),
    mean_core_share_pct = round(100 * mean(core_share), 1)))

  structure(list(ps_core = ps_core, taxa = taxa, taxonomy = tax_df,
    detection = detection, prevalence = prevalence, rank = rank,
    n_core = length(taxa), n_total = phyloseq::ntaxa(rel),
    core_share_per_sample = core_share,
    mean_core_share = mean(core_share),
    threshold_sentence = sprintf(paste("Core taxa were defined as those present",
      "at a relative abundance of at least %.3g in at least %.0f%% of samples;",
      "%d of %d taxa met this, accounting for a mean %.1f%% of each sample."),
      detection, 100 * prevalence, length(taxa), phyloseq::ntaxa(rel),
      100 * mean(core_share))),
    class = c("mfg_core", "list"))
}

#' @export
print.mfg_core <- function(x, ...) {
  cat("=== Core microbiome ===\n")
  cat(strwrap(x$threshold_sentence, width = 78), sep = "\n")
  if (x$n_core) {
    cat("\nCore taxa:\n")
    show_cols <- intersect(c("taxon", "Phylum", "Family", "Genus", "Species"),
                           names(x$taxonomy))
    print(utils::head(x$taxonomy[, show_cols, drop = FALSE], 25), row.names = FALSE)
  } else {
    cat("\nNo taxa met these thresholds. Either the community has no shared core",
        "at this strictness, or the thresholds are too strict for this data —",
        "try core_threshold_scan() to see how the count varies.\n")
  }
  invisible(x)
}

#' How the core size varies with both thresholds.
#'
#' The honest way to present a core: not one number at one arbitrary cut, but the
#' surface. If the core is 3 taxa at prevalence 0.9 and 47 at 0.5, that is the
#' finding, and quoting either alone is a choice the reader cannot see.
core_threshold_scan <- function(ps, detections = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2),
                                prevalences = c(0.5, 0.7, 0.8, 0.9, 0.95, 1.0),
                                rank = NULL) {
  work <- if (!is.null(rank)) agglomerate_and_summarise(ps, rank)$ps else ps
  rel  <- if (identical(mfg_normalization(work), "tss")) work else tss_transform(work)
  mat  <- mfg_otu_taxa_as_rows(rel)
  n    <- ncol(mat)

  rows <- list()
  for (d in detections) for (p in prevalences) {
    n_core <- sum(rowSums(mat > d) >= ceiling(p * n))
    rows[[length(rows) + 1]] <- data.frame(detection = d, prevalence = p,
                                           n_core = n_core, stringsAsFactors = FALSE)
  }
  tab <- do.call(rbind, rows)
  mfg_log("composition", "core_threshold_scan", list(
    n_combinations = nrow(tab), rank = rank %||% "as-is",
    range_n_core = sprintf("%d-%d", min(tab$n_core), max(tab$n_core))))
  tab
}

#' Core taxa counted per phylum.
#'
#' Demo 6's core-per-phylum view: which lineages contribute the shared backbone.
core_by_phylum <- function(ps, detection = 0.001, prevalence = 0.5,
                           rank = "Phylum") {
  cr <- core_taxa(ps, detection = detection, prevalence = prevalence)
  if (!cr$n_core) {
    return(data.frame(level = character(), n_core = integer(),
                      n_total = integer(), pct_core = numeric()))
  }
  tt_all  <- mfg_tax(ps)
  r <- safe_rank(ps, rank)
  if (is.null(r)) stop("Rank '", rank, "' not present.", call. = FALSE)

  core_by <- table(tt_all[cr$taxa, r], useNA = "ifany")
  all_by  <- table(tt_all[[r]], useNA = "ifany")
  lv <- names(all_by)
  out <- data.frame(
    level = lv,
    n_core = as.integer(core_by[lv]),
    n_total = as.integer(all_by[lv]),
    stringsAsFactors = FALSE)
  out$n_core[is.na(out$n_core)] <- 0L
  out$pct_core <- round(100 * out$n_core / out$n_total, 1)
  out <- out[order(-out$n_core, -out$n_total), ]
  mfg_log("composition", "core_by_rank",
          list(rank = r, detection = detection, prevalence = prevalence,
               n_levels_with_core = sum(out$n_core > 0)))
  out
}

# ── Prevalence/abundance structure ───────────────────────────────────────────

#' Rank-abundance (dominance) profile.
#'
#' How evenly the community is distributed. Reports the share held by the single
#' most abundant taxon and by the top five, which is the summary that makes
#' "dominated" quantitative instead of impressionistic.
dominance_summary <- function(ps, group_var = NULL) {
  rel <- if (identical(mfg_normalization(ps), "tss")) ps else tss_transform(ps)
  mat <- mfg_otu_taxa_as_rows(rel)

  per_sample <- apply(mat, 2, function(x) {
    s <- sort(x, decreasing = TRUE)
    c(top1 = s[1], top5 = sum(s[seq_len(min(5, length(s)))]),
      n_to_50pct = which(cumsum(s) >= 0.5)[1])
  })
  df <- as.data.frame(t(per_sample))
  names(df) <- c("top1_share", "top5_share", "n_taxa_to_50pct")
  df$sample <- rownames(df)

  if (!is.null(group_var)) {
    meta <- mfg_meta(rel); meta$sample <- rownames(meta)
    df <- merge(df, meta[, c("sample", group_var)], by = "sample", all.x = TRUE)
  }
  mfg_log("composition", "dominance",
          list(mean_top1_pct = round(100 * mean(df$top1_share), 1),
               mean_top5_pct = round(100 * mean(df$top5_share), 1),
               median_taxa_to_50pct = stats::median(df$n_taxa_to_50pct)))
  df
}
