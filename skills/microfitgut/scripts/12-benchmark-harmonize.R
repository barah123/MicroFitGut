# =============================================================================
# MicroFitGut — 12-benchmark-harmonize.R
#
# Making a published result and a fresh reanalysis comparable before comparing
# them.
#
# This is the step benchmark systems usually skip, and skipping it manufactures
# disagreement. A taxon that "disappears" between the paper and the replication
# is most often the same organism under a different name: SILVA 138 splits
# Lactobacillus into Lacticaseibacillus, Limosilactobacillus and others;
# Eubacterium rectale becomes Agathobacter rectalis in GTDB; Bacteroidetes was
# renamed Bacteroidota. Counting those as misses inflates the divergence and then
# the attribution step invents a cause for it.
#
# So: resolve names first, log every label that could not be resolved, and
# compare at the deepest rank where both sides have confident assignments.
#
# Requires: utils.R
# =============================================================================

# ── Synonym table ────────────────────────────────────────────────────────────
#
# Curated renamings that actually bite in 16S microbiome comparisons. This is
# deliberately a small, auditable table rather than a live NCBI/GTDB lookup:
# a benchmark has to be reproducible, and a name resolution that changes when a
# remote database updates makes last month's benchmark unrepeatable.
#
# Extend it per project in reference/13 and record what was added.

MFG_TAXON_SYNONYMS <- list(
  # --- Phylum renamings (2021 ICNP validation of phylum names) ---
  "Bacteroidetes"      = "Bacteroidota",
  "Firmicutes"         = "Bacillota",
  "Proteobacteria"     = "Pseudomonadota",
  "Actinobacteria"     = "Actinomycetota",
  "Fusobacteria"       = "Fusobacteriota",
  "Verrucomicrobia"    = "Verrucomicrobiota",
  "Acidobacteria"      = "Acidobacteriota",
  "Cyanobacteria"      = "Cyanobacteriota",
  "Spirochaetes"       = "Spirochaetota",
  "Tenericutes"        = "Mycoplasmatota",
  "Deinococcus-Thermus" = "Deinococcota",
  "Chloroflexi"        = "Chloroflexota",
  "Planctomycetes"     = "Planctomycetota",
  "Synergistetes"      = "Synergistota",
  "Epsilonbacteraeota" = "Campylobacterota",
  "Campylobacteria"    = "Campylobacterota",

  # --- Lactobacillus split (Zheng et al. 2020) ---
  "Lactobacillus casei"        = "Lacticaseibacillus casei",
  "Lactobacillus paracasei"    = "Lacticaseibacillus paracasei",
  "Lactobacillus rhamnosus"    = "Lacticaseibacillus rhamnosus",
  "Lactobacillus plantarum"    = "Lactiplantibacillus plantarum",
  "Lactobacillus reuteri"      = "Limosilactobacillus reuteri",
  "Lactobacillus fermentum"    = "Limosilactobacillus fermentum",
  "Lactobacillus brevis"       = "Levilactobacillus brevis",
  "Lactobacillus salivarius"   = "Ligilactobacillus salivarius",
  "Lactobacillus sakei"        = "Latilactobacillus sakei",
  "Lactobacillus curvatus"     = "Latilactobacillus curvatus",

  # --- Gut genera commonly reassigned ---
  "Eubacterium rectale"        = "Agathobacter rectalis",
  "Ruminococcus gnavus"        = "Mediterraneibacter gnavus",
  "Ruminococcus torques"       = "Mediterraneibacter torques",
  "Clostridium difficile"      = "Clostridioides difficile",
  "Clostridium clostridioforme" = "Enterocloster clostridioformis",
  "Bacteroides dorei"          = "Phocaeicola dorei",
  "Bacteroides vulgatus"       = "Phocaeicola vulgatus",
  "Bacteroides massiliensis"   = "Phocaeicola massiliensis",
  "Propionibacterium acnes"    = "Cutibacterium acnes",
  "Propionibacterium"          = "Cutibacterium",
  "Peptoclostridium difficile" = "Clostridioides difficile",
  "Blautia producta"           = "Blautia producta",
  "Escherichia"                = "Escherichia-Shigella",
  "Shigella"                   = "Escherichia-Shigella",

  # --- Skin genera ---
  "Staphylococcus epidermidis" = "Staphylococcus epidermidis",
  "Malassezia furfur"          = "Malassezia furfur"
)

# Labels that mean "no assignment" rather than naming a taxon. Treating these as
# names creates a spurious shared taxon called "unclassified" that matches
# everything.
MFG_UNRESOLVED_LABELS <- c("", "NA", "na", "unclassified", "Unclassified",
                           "unidentified", "Unidentified", "unknown", "Unknown",
                           "uncultured", "Uncultured", "metagenome",
                           "unassigned", "Unassigned", "no_rank", "None",
                           "__", "g__", "s__", "f__", "o__", "c__", "p__", "d__", "k__")

#' Normalise a taxon label to a comparable form.
#'
#' Strips rank prefixes, collapses whitespace and underscores, drops bracketed
#' provisional-name markers such as [Eubacterium], and maps through the synonym
#' table. Returns NA for anything that means "unassigned", so those never match.
harmonize_taxon <- function(x, synonyms = MFG_TAXON_SYNONYMS,
                            case_insensitive = TRUE) {
  x <- as.character(x)
  x <- sub("^[dkpcofgst]__", "", trimws(x))
  x <- gsub("[\\[\\]]", "", x, perl = TRUE)
  x <- gsub("_", " ", x)
  x <- gsub("\\s+", " ", trimws(x))
  # Trailing qualifiers that are not part of the name.
  x <- sub("\\s+(group|clade|complex|sensu stricto\\s*\\d*)$", "", x,
           ignore.case = TRUE)
  x[x %in% MFG_UNRESOLVED_LABELS] <- NA_character_
  x[!nzchar(x %||% "")] <- NA_character_

  if (length(synonyms)) {
    keys <- names(synonyms)
    if (case_insensitive) {
      idx <- match(tolower(x), tolower(keys))
    } else {
      idx <- match(x, keys)
    }
    hit <- !is.na(idx)
    x[hit] <- unlist(synonyms)[idx[hit]]
  }
  x
}

#' Harmonise a taxonomy table rank by rank.
harmonize_taxonomy <- function(tax_df, ranks = NULL, synonyms = MFG_TAXON_SYNONYMS) {
  ranks <- ranks %||% intersect(TAX_RANKS, names(tax_df))
  out <- tax_df
  for (r in ranks) out[[r]] <- harmonize_taxon(tax_df[[r]], synonyms)
  n_changed <- sum(vapply(ranks, function(r)
    sum(!identical(as.character(tax_df[[r]]), as.character(out[[r]]))), integer(1)))
  mfg_log("benchmark", "taxonomy_harmonized",
          list(ranks = paste(ranks, collapse = ","),
               n_cells_changed = sum(vapply(ranks, function(r)
                 sum(as.character(tax_df[[r]]) != as.character(out[[r]]), na.rm = TRUE),
                 integer(1)))))
  out
}

# ── Matching two taxon lists ─────────────────────────────────────────────────

#' Match a published taxon list against a reanalysis taxon list.
#'
#' Tries progressively looser matching and records which rule matched each pair,
#' because "matched after genus truncation" is a weaker claim than "matched
#' exactly" and the benchmark should say which it relied on.
#'
#' Returns the matches, the unmatched on each side, and a resolution log. Every
#' unresolved label is reported rather than dropped, since an unresolvable name
#' is itself a benchmark finding about the paper's reporting.
match_taxa <- function(published, reanalysis,
                       published_tax = NULL, reanalysis_tax = NULL,
                       synonyms = MFG_TAXON_SYNONYMS) {

  pub_h <- harmonize_taxon(published, synonyms)
  rea_h <- harmonize_taxon(reanalysis, synonyms)

  pub_df <- data.frame(published = published, published_h = pub_h,
                       stringsAsFactors = FALSE)
  rea_df <- data.frame(reanalysis = reanalysis, reanalysis_h = rea_h,
                       stringsAsFactors = FALSE)

  matches <- list()
  used_rea <- character()

  # Rule 1: exact match on harmonised names.
  for (i in seq_len(nrow(pub_df))) {
    if (is.na(pub_df$published_h[i])) next
    j <- which(rea_df$reanalysis_h == pub_df$published_h[i] &
               !rea_df$reanalysis %in% used_rea)
    if (length(j)) {
      matches[[length(matches) + 1]] <- data.frame(
        published = pub_df$published[i], reanalysis = rea_df$reanalysis[j[1]],
        matched_as = pub_df$published_h[i], rule = "exact (harmonised)",
        confidence = "high", stringsAsFactors = FALSE)
      used_rea <- c(used_rea, rea_df$reanalysis[j[1]])
    }
  }

  # Rule 2: genus-level match — the published name is a species and the
  # reanalysis resolves only to genus, or the reverse. This is a real match at a
  # shallower rank, not a failure, but it is a weaker claim.
  matched_pub <- vapply(matches, function(m) m$published, character(1))
  # Truncate FIRST, then harmonise the genus token — not the reverse. A synonym
  # that rewrites a bare genus into a compound form (Escherichia ->
  # Escherichia-Shigella) does not fire on a multi-word name that contains that
  # genus, so harmonising before truncating leaves the two sides in different
  # namespaces and misses a real match. That silently deflates recovery rate,
  # which is worse than not harmonising at all.
  genus_key <- function(x) harmonize_taxon(sub("\\s.*$", "", x), synonyms)
  pub_genus_all <- genus_key(pub_df$published)
  rea_genus_all <- genus_key(rea_df$reanalysis)
  for (i in seq_len(nrow(pub_df))) {
    if (pub_df$published[i] %in% matched_pub || is.na(pub_df$published_h[i])) next
    pub_genus <- pub_genus_all[i]
    if (is.na(pub_genus)) next
    j <- which(rea_genus_all == pub_genus &
               !rea_df$reanalysis %in% used_rea & !is.na(rea_genus_all))
    if (length(j)) {
      matches[[length(matches) + 1]] <- data.frame(
        published = pub_df$published[i], reanalysis = rea_df$reanalysis[j[1]],
        matched_as = pub_genus, rule = "genus-level truncation",
        confidence = "medium", stringsAsFactors = FALSE)
      used_rea <- c(used_rea, rea_df$reanalysis[j[1]])
    }
  }

  # Rule 3: match through the taxonomy tables, when the reanalysis names are ASV
  # ids and the published names are taxon names.
  if (!is.null(reanalysis_tax)) {
    matched_pub <- vapply(matches, function(m) m$published, character(1))
    rt <- harmonize_taxonomy(reanalysis_tax, synonyms = synonyms)
    ranks <- intersect(rev(TAX_RANKS), names(rt))
    for (i in seq_len(nrow(pub_df))) {
      if (pub_df$published[i] %in% matched_pub || is.na(pub_df$published_h[i])) next
      target <- pub_df$published_h[i]
      hit <- NULL; hit_rank <- NA_character_
      for (r in ranks) {
        cand <- rownames(rt)[which(rt[[r]] == target)]
        cand <- setdiff(intersect(cand, reanalysis), used_rea)
        if (length(cand)) { hit <- cand[1]; hit_rank <- r; break }
      }
      if (!is.null(hit)) {
        matches[[length(matches) + 1]] <- data.frame(
          published = pub_df$published[i], reanalysis = hit,
          matched_as = target, rule = sprintf("via taxonomy table at %s", hit_rank),
          confidence = if (hit_rank %in% c("Species", "Genus")) "high" else "medium",
          stringsAsFactors = FALSE)
        used_rea <- c(used_rea, hit)
      }
    }
  }

  match_df <- if (length(matches)) do.call(rbind, matches) else
    data.frame(published = character(), reanalysis = character(),
               matched_as = character(), rule = character(), confidence = character())

  unmatched_published <- setdiff(published, match_df$published)
  unmatched_reanalysis <- setdiff(reanalysis, match_df$reanalysis)
  unresolvable <- published[is.na(pub_h)]

  out <- list(
    matches = match_df,
    n_matched = nrow(match_df),
    unmatched_published = unmatched_published,
    unmatched_reanalysis = unmatched_reanalysis,
    unresolvable_published_labels = unresolvable,
    n_published = length(published), n_reanalysis = length(reanalysis),
    recovery_rate = if (length(published)) nrow(match_df) / length(published) else NA_real_,
    by_rule = table(match_df$rule),
    by_confidence = table(match_df$confidence)
  )
  mfg_log("benchmark", "taxa_matched", list(
    n_published = length(published), n_reanalysis = length(reanalysis),
    n_matched = nrow(match_df),
    recovery_rate = round(out$recovery_rate %||% NA, 3),
    n_unmatched_published = length(unmatched_published),
    n_unresolvable_labels = length(unresolvable),
    rules = paste(sprintf("%s=%d", names(out$by_rule), as.integer(out$by_rule)),
                  collapse = " ")))
  class(out) <- c("mfg_taxon_match", "list")
  out
}

#' @export
print.mfg_taxon_match <- function(x, ...) {
  cat("=== Taxon harmonisation ===\n")
  cat(sprintf("Published: %d taxa | Reanalysis: %d taxa | Matched: %d (%.1f%% of published recovered)\n",
              x$n_published, x$n_reanalysis, x$n_matched, 100 * (x$recovery_rate %||% 0)))
  if (length(x$by_rule)) {
    cat("\nMatched by rule:\n")
    for (r in names(x$by_rule)) cat(sprintf("  %-28s %d\n", r, x$by_rule[[r]]))
  }
  if (length(x$unmatched_published)) {
    cat(sprintf("\nPublished taxa with no counterpart (%d):\n", length(x$unmatched_published)))
    cat("  ", paste(utils::head(x$unmatched_published, 15), collapse = ", "), "\n", sep = "")
    if (length(x$unmatched_published) > 15)
      cat(sprintf("  ... and %d more\n", length(x$unmatched_published) - 15))
  }
  if (length(x$unresolvable_published_labels)) {
    cat(sprintf("\n! %d published labels could not be resolved to a taxon name (%s).\n",
                length(x$unresolvable_published_labels),
                paste(utils::head(x$unresolvable_published_labels, 5), collapse = ", ")))
    cat(strwrap(paste("These are unresolvable as reported. That is a finding about",
      "the paper's reporting, not a replication failure, and must be stated",
      "separately from genuine misses."), width = 78, prefix = "  "), sep = "\n")
  }
  invisible(x)
}

# ── Aligning two abundance tables ────────────────────────────────────────────

#' Align a published abundance table with a reanalysis table.
#'
#' Returns both matrices restricted to shared samples and matched taxa, in the
#' same order, ready for the distance and correlation metrics in 13-benchmark.
#' Reports what was lost on each side, because a concordance computed on the 12
#' taxa two tables happen to share is not a concordance between the tables.
align_abundance_tables <- function(published_mat, reanalysis_mat,
                                   published_tax = NULL, reanalysis_tax = NULL,
                                   synonyms = MFG_TAXON_SYNONYMS,
                                   relative = TRUE) {

  shared_samples <- intersect(colnames(published_mat), colnames(reanalysis_mat))
  if (!length(shared_samples)) {
    stop("No sample IDs in common between the two tables.\n",
         "Published: ", paste(utils::head(colnames(published_mat), 3), collapse = ", "),
         "\nReanalysis: ", paste(utils::head(colnames(reanalysis_mat), 3), collapse = ", "),
         "\nSample naming must be reconciled before any concordance metric means anything.",
         call. = FALSE)
    }

  tm <- match_taxa(rownames(published_mat), rownames(reanalysis_mat),
                   published_tax, reanalysis_tax, synonyms)
  if (!tm$n_matched) {
    stop("No taxa could be matched between the two tables. Check that both use ",
         "taxon names (not ASV ids on one side only) and that a taxonomy table ",
         "was supplied for the reanalysis.", call. = FALSE)
  }

  p <- published_mat[tm$matches$published, shared_samples, drop = FALSE]
  r <- reanalysis_mat[tm$matches$reanalysis, shared_samples, drop = FALSE]
  rownames(p) <- rownames(r) <- tm$matches$matched_as

  if (isTRUE(relative)) {
    p <- sweep(p, 2, pmax(colSums(p), 1e-12), "/")
    r <- sweep(r, 2, pmax(colSums(r), 1e-12), "/")
  }

  out <- list(
    published = p, reanalysis = r,
    shared_samples = shared_samples,
    taxon_match = tm,
    n_samples_published_only = length(setdiff(colnames(published_mat), shared_samples)),
    n_samples_reanalysis_only = length(setdiff(colnames(reanalysis_mat), shared_samples)),
    pct_published_abundance_retained = 100 * mean(
      colSums(published_mat[tm$matches$published, shared_samples, drop = FALSE]) /
      pmax(colSums(published_mat[, shared_samples, drop = FALSE]), 1e-12)),
    pct_reanalysis_abundance_retained = 100 * mean(
      colSums(reanalysis_mat[tm$matches$reanalysis, shared_samples, drop = FALSE]) /
      pmax(colSums(reanalysis_mat[, shared_samples, drop = FALSE]), 1e-12)),
    relative = relative
  )
  mfg_log("benchmark", "tables_aligned", list(
    shared_samples = length(shared_samples),
    matched_taxa = tm$n_matched,
    pct_published_abundance_retained = round(out$pct_published_abundance_retained, 1),
    pct_reanalysis_abundance_retained = round(out$pct_reanalysis_abundance_retained, 1),
    samples_published_only = out$n_samples_published_only,
    samples_reanalysis_only = out$n_samples_reanalysis_only))

  if (out$pct_published_abundance_retained < 70) {
    warning(sprintf(paste("Matched taxa account for only %.0f%% of the published",
      "table's abundance. Concordance metrics computed on this alignment describe",
      "a minority of the published community, and must be reported with that",
      "coverage figure."), out$pct_published_abundance_retained), call. = FALSE)
  }
  class(out) <- c("mfg_alignment", "list")
  out
}

#' @export
print.mfg_alignment <- function(x, ...) {
  cat("=== Table alignment ===\n")
  cat(sprintf("Shared samples: %d (published-only %d, reanalysis-only %d)\n",
              length(x$shared_samples), x$n_samples_published_only,
              x$n_samples_reanalysis_only))
  cat(sprintf("Matched taxa: %d\n", x$taxon_match$n_matched))
  cat(sprintf("Abundance coverage: %.1f%% of the published table, %.1f%% of the reanalysis\n",
              x$pct_published_abundance_retained, x$pct_reanalysis_abundance_retained))
  cat(sprintf("Values: %s\n", if (x$relative) "relative abundance" else "as supplied"))
  invisible(x)
}
