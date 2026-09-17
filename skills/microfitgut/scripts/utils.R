# =============================================================================
# MicroFitGut — utils.R
#
# Shared helpers. Most are lifted from CanisLupus 2.0 (app_original.R) where the
# comments explain a decision the code alone would not; those comments are kept
# because they are the reason the function exists.
#
# Nothing here computes a statistic. These are the guards that stop a later
# stage from producing a number that looks valid and encodes nothing.
# =============================================================================

TAX_RANKS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a

# ── Run logging ──────────────────────────────────────────────────────────────
#
# Every stage appends to one log. 11-report.R reads it back, and the verifier
# subagent checks report claims against it. A number that is not in here did not
# come from a run.

MFG_LOG <- new.env(parent = emptyenv())
MFG_LOG$entries <- list()
MFG_LOG$run_id  <- NULL

#' Start a run and fix its identifier.
#'
#' The run ID goes into every output filename and into the report header, so a
#' figure can always be traced back to the log that produced it.
mfg_start_run <- function(label = "analysis", outdir = "output") {
  stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  MFG_LOG$run_id  <- paste0(label, "-", stamp)
  MFG_LOG$entries <- list()
  MFG_LOG$outdir  <- outdir
  dir.create(file.path(outdir, MFG_LOG$run_id), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, MFG_LOG$run_id, "figures"), showWarnings = FALSE)
  dir.create(file.path(outdir, MFG_LOG$run_id, "tables"),  showWarnings = FALSE)
  mfg_log("run", "started", list(run_id = MFG_LOG$run_id, outdir = outdir))
  invisible(MFG_LOG$run_id)
}

mfg_run_id <- function() {
  if (is.null(MFG_LOG$run_id)) mfg_start_run()
  MFG_LOG$run_id
}

mfg_run_dir <- function(sub = NULL) {
  d <- file.path(MFG_LOG$outdir %||% "output", mfg_run_id())
  if (!is.null(sub)) d <- file.path(d, sub)
  d
}

#' Append one structured record to the run log.
mfg_log <- function(stage, event, detail = list()) {
  MFG_LOG$entries[[length(MFG_LOG$entries) + 1]] <- list(
    stage = stage, event = event,
    time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    detail = detail
  )
  invisible(TRUE)
}

#' Write the log to disk as JSON-ish text and a flat CSV.
mfg_write_log <- function() {
  d <- mfg_run_dir()
  flat <- do.call(rbind, lapply(MFG_LOG$entries, function(e) {
    data.frame(stage = e$stage, event = e$event, time = e$time,
               detail = paste(names(e$detail), unlist(lapply(e$detail, function(x)
                 paste(format(x), collapse = "/"))), sep = "=", collapse = "; "),
               stringsAsFactors = FALSE)
  }))
  utils::write.csv(flat, file.path(d, "run_log.csv"), row.names = FALSE)
  invisible(flat)
}

mfg_get_log <- function() MFG_LOG$entries

# ── Saving outputs ───────────────────────────────────────────────────────────

#' Save a ggplot at the house dimensions and log it.
#'
#' Every figure goes through here so nothing reaches the report without a log
#' entry naming the file and the stage that made it.
mfg_save_plot <- function(p, name, stage = "plot",
                          width = 7, height = 5, dpi = 300, device = "png") {
  path <- file.path(mfg_run_dir("figures"), paste0(name, ".", device))
  if (inherits(p, "ggplot") || inherits(p, "patchwork")) {
    ggplot2::ggsave(path, p, width = width, height = height, dpi = dpi)
  } else {
    grDevices::png(path, width = width, height = height, units = "in", res = dpi)
    on.exit(grDevices::dev.off(), add = TRUE)
    # A gtable or grob (pheatmap, heat trees, base-graphics recordings) has to be
    # drawn, not printed — print() on a gtable writes its structure to the console
    # and leaves the device blank.
    if (inherits(p, c("gtable", "grob", "gTree"))) {
      grid::grid.newpage(); grid::grid.draw(p)
    } else if (inherits(p, "pheatmap")) {
      grid::grid.newpage(); grid::grid.draw(p$gtable)
    } else if (is.function(p)) {
      p()
    } else {
      print(p)
    }
  }
  mfg_log(stage, "figure_saved",
          list(file = basename(path), width = width, height = height, dpi = dpi))
  invisible(path)
}

#' Save a table and log it.
mfg_save_table <- function(df, name, stage = "table") {
  path <- file.path(mfg_run_dir("tables"), paste0(name, ".csv"))
  utils::write.csv(df, path, row.names = FALSE)
  mfg_log(stage, "table_saved",
          list(file = basename(path), rows = nrow(df), cols = ncol(df)))
  invisible(path)
}

# ── phyloseq guards ──────────────────────────────────────────────────────────

#' Carry MicroFitGut's bookkeeping attributes across a phyloseq operation.
#'
#' phyloseq's own functions (prune_taxa, rarefy_even_depth, transform_sample_counts,
#' tax_glom, subset_samples) rebuild the object and drop any attributes we attached.
#' Without this, rarefying loses the tree-provenance marker — so Faith's PD and
#' UniFrac get refused at exactly the point in the pipeline where they become
#' valid. Every function that returns a derived object passes through here.
MFG_CARRIED_ATTRS <- c("mfg_tree_is_real", "mfg_is_relative", "mfg_is_functional",
                       "mfg_normalization", "mfg_normalization_detail")

mfg_carry_attrs <- function(to, from) {
  for (a in MFG_CARRIED_ATTRS) {
    v <- attr(from, a)
    if (!is.null(v) && is.null(attr(to, a))) attr(to, a) <- v
  }
  to
}

#' Subset samples or taxa without losing MicroFitGut's tracked state.
#'
#' Every QC function in 02-qc.R routes through mfg_carry_attrs(), but a direct
#' phyloseq::prune_samples() or subset_samples() call does not: the returned
#' object is rebuilt and the normalization marker and tree provenance go with it.
#' The next analysis then reports "normalization is not tracked (inferred:
#' unknown)" and refuses — on an object that WAS correctly normalized a line
#' earlier. Found when subsetting a TSS-transformed object down to root samples.
#'
#' Use these in place of the phyloseq calls whenever the object is mid-pipeline.
mfg_prune_samples <- function(keep, ps) {
  mfg_carry_attrs(phyloseq::prune_samples(keep, ps), ps)
}

mfg_prune_taxa <- function(keep, ps) {
  mfg_carry_attrs(phyloseq::prune_taxa(keep, ps), ps)
}

#' Write metadata back onto an object without losing tracked state.
#'
#' `sample_data(ps) <- ...` rebuilds the object and drops MicroFitGut's
#' attributes. mfg_set_normalization() records a registry backstop keyed on a
#' content fingerprint for exactly this reason, but the fingerprint changes when
#' the object has been subset — so subset-then-assign loses the state and the
#' next analysis refuses an object that was correctly normalized two lines
#' earlier. Found adding a derived column to a subset of a TSS-transformed table.
mfg_set_meta <- function(ps, md) {
  out <- ps
  phyloseq::sample_data(out) <- phyloseq::sample_data(as.data.frame(md))
  mfg_carry_attrs(out, ps)
}

#' Add or replace one metadata column, carrying state.
mfg_add_meta <- function(ps, name, value) {
  md <- mfg_meta(ps)
  md[[name]] <- value
  mfg_set_meta(ps, md)
}

#' Metadata as a genuine base data.frame.
#'
#' as.data.frame() on a phyloseq sample_data returns an object still carrying the
#' sample_data class, and merge() silently produces all-NA columns against it.
#' Every stage goes through here so that failure mode cannot recur.
mfg_meta <- function(ps, add_sample_col = FALSE) {
  sd <- phyloseq::sample_data(ps, errorIfNULL = FALSE)
  if (is.null(sd)) return(NULL)
  md <- as(sd, "data.frame")
  md <- as.data.frame(md, stringsAsFactors = FALSE)
  class(md) <- "data.frame"
  if (add_sample_col) md[[".sample"]] <- rownames(md)
  md
}

#' Taxonomy as a genuine base data.frame, for the same reason.
mfg_tax <- function(ps) {
  tt <- phyloseq::tax_table(ps, errorIfNULL = FALSE)
  if (is.null(tt)) return(NULL)
  td <- as.data.frame(as(tt, "matrix"), stringsAsFactors = FALSE)
  class(td) <- "data.frame"
  td
}

#' Abundance matrix with samples as rows, which is what vegan expects.
mfg_otu_samples_as_rows <- function(ps) {
  m <- as(phyloseq::otu_table(ps), "matrix")
  if (phyloseq::taxa_are_rows(ps)) m <- t(m)
  m
}

#' Abundance matrix with taxa as rows.
mfg_otu_taxa_as_rows <- function(ps) {
  m <- as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) m <- t(m)
  m
}

#' Safe tax_glom: returns the object unchanged when the rank is absent.
safe_tax_glom <- function(ps, rank) {
  if (rank %in% phyloseq::rank_names(ps)) phyloseq::tax_glom(ps, rank) else ps
}

#' Resolve a requested taxonomic rank against what the dataset actually holds.
#'
#' Rank availability varies by source: an amplicon taxonomy often stops at
#' Genus, while shotgun profiles usually resolve to Species. Returns NULL when
#' the rank is absent so callers can explain, rather than pushing an
#' unresolvable column name into phyloseq.
safe_rank <- function(ps, rank) {
  ranks <- phyloseq::rank_names(ps)
  if (is.null(rank) || !length(rank) || is.na(rank[1]) || !nzchar(rank[1])) {
    return(if (length(ranks)) ranks[length(ranks)] else NULL)
  }
  if (rank[1] %in% ranks) rank[1] else NULL
}

# ── State that survives phyloseq object surgery ──────────────────────────────
#
# Attributes attached to a phyloseq object are dropped not only by phyloseq's
# functions but also by slot assignment — `sample_data(ps)$x <- ...` rebuilds the
# object. mfg_carry_attrs() cannot help there, because the caller never sees an
# intermediate to copy from.
#
# Two storage strategies that do survive:
#   - the tree-provenance marker lives on the phylo object itself, which travels
#     inside the slot through prune, rarefy, tax_glom, subset and slot assignment
#   - everything else falls back to a session registry keyed on a cheap
#     fingerprint of the abundance table

MFG_REGISTRY <- new.env(parent = emptyenv())

#' Cheap content fingerprint, stable across metadata edits.
#'
#' Deliberately ignores sample_data, so editing metadata does not lose the
#' record. Changes when taxa, samples or counts change — which is when the
#' recorded state genuinely may no longer apply.
mfg_fingerprint <- function(ps) {
  otu <- phyloseq::otu_table(ps, errorIfNULL = FALSE)
  if (is.null(otu)) return(NULL)
  m <- as(otu, "matrix")
  paste(phyloseq::ntaxa(ps), phyloseq::nsamples(ps),
        format(sum(m), digits = 15),
        substr(paste(sort(phyloseq::taxa_names(ps)), collapse = ""), 1, 64),
        sep = "|")
}

mfg_registry_set <- function(ps, key, value) {
  fp <- mfg_fingerprint(ps)
  if (is.null(fp)) return(invisible(NULL))
  slot <- MFG_REGISTRY[[fp]] %||% list()
  slot[[key]] <- value
  assign(fp, slot, envir = MFG_REGISTRY)
  invisible(value)
}

mfg_registry_get <- function(ps, key) {
  fp <- mfg_fingerprint(ps)
  if (is.null(fp)) return(NULL)
  slot <- MFG_REGISTRY[[fp]]
  if (is.null(slot)) return(NULL)
  slot[[key]]
}

#' Is the phylogeny in this object real, or a placeholder?
#'
#' A random tree merged in to satisfy phyloseq's slot is not a phylogeny.
#' UniFrac computed on one returns numbers that look valid but encode nothing,
#' so every phylogenetic metric must check here first. Demo 3 and Demo 6 both
#' attach `rtree()` for teaching purposes, which is exactly the case this
#' catches.
#'
#' Returns NA when provenance is unrecorded — absent evidence we cannot prove a
#' tree is real, so the caller is told to confirm rather than handed a silent
#' TRUE.
mfg_tree_is_real <- function(ps, attr_name = "mfg_tree_is_real") {
  tree <- phyloseq::phy_tree(ps, errorIfNULL = FALSE)
  if (is.null(tree)) return(FALSE)
  # Primary store: on the phylo object, which survives slot assignment.
  flag <- attr(tree, "mfg_tree_is_real")
  if (!is.null(flag)) return(isTRUE(flag))
  # Object-level attribute, for objects marked before this became the primary.
  flag <- attr(ps, attr_name)
  if (!is.null(flag)) return(isTRUE(flag))
  flag <- mfg_registry_get(ps, "tree_is_real")
  if (!is.null(flag)) return(isTRUE(flag))
  NA
}

#' Mark a tree as genuinely estimated from the sequence data.
#'
#' Writes the marker onto the phylo object so it survives `sample_data(ps) <- `,
#' `prune_taxa`, `rarefy_even_depth`, `tax_glom` and `subset_samples`. Also
#' records it in the session registry as a backstop.
mfg_mark_tree_real <- function(ps, real = TRUE) {
  tree <- phyloseq::phy_tree(ps, errorIfNULL = FALSE)
  if (is.null(tree)) {
    warning("No tree in this object; nothing to mark.", call. = FALSE)
    return(ps)
  }
  attr(tree, "mfg_tree_is_real") <- isTRUE(real)
  phyloseq::phy_tree(ps) <- tree
  attr(ps, "mfg_tree_is_real") <- isTRUE(real)
  mfg_registry_set(ps, "tree_is_real", isTRUE(real))
  mfg_log("intake", "tree_provenance_marked",
          list(real = isTRUE(real), n_tips = length(tree$tip.label)))
  ps
}

#' Distance metrics that are meaningful for this object.
#'
#' UniFrac is phylogenetic, so it is only offered with a real tree. Without one
#' those options are withheld instead of silently returning a wrong number.
allowed_beta_distances <- function(tree_is_real) {
  all_metrics <- c("bray", "jaccard", "jsd", "euclidean", "unifrac", "wunifrac")
  if (isTRUE(tree_is_real)) all_metrics
  else setdiff(all_metrics, c("unifrac", "wunifrac"))
}

#' Resolve the distance metric to compute with, or NULL to refuse.
resolve_beta_distance <- function(selected, tree_is_real) {
  if (is.null(selected) || !nzchar(selected)) return("bray")
  if (selected %in% c("unifrac", "wunifrac") && !isTRUE(tree_is_real)) return(NULL)
  selected
}

# ── Design detection ─────────────────────────────────────────────────────────

#' Is this a complete, unreplicated block design?
#'
#' TRUE when every combination of group and block holds exactly one sample,
#' which is what a Friedman test requires. Repeated-measures data that fails
#' this check (unbalanced, or several samples per cell) must fall back to an
#' unblocked test.
is_complete_block_design <- function(group, block) {
  if (is.null(group) || is.null(block)) return(FALSE)
  if (length(group) != length(block)) return(FALSE)
  if (!length(group)) return(FALSE)
  tab <- table(as.factor(group), as.factor(block))
  length(tab) > 0 && all(tab == 1)
}

#' Does this metadata contain repeated measures on the same unit?
#'
#' Returns the subject variable and how many observations each subject has.
#' A design that is repeated-measures and gets an independent-samples test is
#' the single most common error the statistician subagent catches, so this is
#' checked explicitly rather than left to inspection.
mfg_detect_repeated_measures <- function(meta, subject_candidates = NULL) {
  # Record this BEFORE the argument is reassigned below: R's missing() stops
  # reporting the original state once the formal has been written to.
  # A caller who names a variable means it: honour the choice, but say so when
  # its shape is unusual rather than accepting it silently.
  explicit <- !is.null(subject_candidates)
  if (is.null(subject_candidates)) {
    # Clustering in microbiome studies is not always called "subject". Ecology
    # clusters by site, plot, block, colony, nest, cage, tank or litter; lab work
    # clusters by batch, run or plate. Missing these means samples that are not
    # independent are silently treated as if they were. The shape test below
    # still rejects the ones that are grouping factors, so widening the candidate
    # list costs nothing and makes the rejection visible instead of absent.
    subject_candidates <- grep(paste0("patient|subject|indiv|host|animal|mouse|",
                               "donor|participant|site|plot|block|colony|nest|",
                               "cage|tank|litter|family|pair|batch|run|plate|",
                               "id$|_id$"),
                               names(meta), ignore.case = TRUE, value = TRUE)
  }
  subject_candidates <- intersect(subject_candidates, names(meta))
  out <- list(repeated = FALSE, subject_var = NULL, max_per_subject = 1L,
              candidates = subject_candidates, rejected = character(0))
  for (v in subject_candidates) {
    tab <- table(meta[[v]])
    # A variable that is unique per row is a sample name, not a subject.
    if (!(length(tab) >= 2 && max(tab) >= 2 && length(tab) < nrow(meta))) next
    # A subject variable partitions samples into MANY SMALL clusters; a grouping
    # factor partitions them into FEW LARGE ones. If one "subject" contributed
    # more samples than there are subjects, it is a factor — sample type, site,
    # treatment — and fitting (1 | that) models the group effect as noise.
    # Caught on a 4-level `host` column over 72 samples that named the sample
    # type, not an organism that was sampled repeatedly.
    if (length(tab) <= max(tab)) {
      if (!explicit) {
        out$rejected <- c(out$rejected, sprintf(
          "%s (%d levels, up to %d samples each - looks like a grouping factor)",
          v, length(tab), max(tab)))
        next
      }
      warning(sprintf(paste("'%s' was given as the subject variable but has only",
        "%d levels with up to %d samples each, which is the shape of a grouping",
        "factor rather than a subject identifier. Using it as instructed; check",
        "that it is not the study's group variable."), v, length(tab), max(tab)),
        call. = FALSE)
    }
    out$repeated        <- TRUE
    out$subject_var     <- v
    out$max_per_subject <- as.integer(max(tab))
    out$n_subjects      <- length(tab)
    break
  }
  # Say what was considered and set aside, so a real subject variable that looks
  # factor-shaped is not silently ignored.
  if (!out$repeated && length(out$rejected)) {
    out$note <- paste0("Considered and rejected as subject identifiers: ",
                       paste(out$rejected, collapse = "; "),
                       ". Pass subject_var explicitly to override.")
  }
  out
}

# ── Shotgun / taxonomic-profile input ────────────────────────────────────────
#
# Shotgun profilers (MetaPhlAn, Kraken2/Bracken and similar) emit a single table
# of lineage strings by sample, rather than the separate abundance and taxonomy
# tables an amplicon pipeline produces. These turn such a table into the same
# two matrices every other stage already works with.

#' Guess the separator used in lineage strings.
#'
#' MetaPhlAn uses "|", QIIME-style and several Kraken exports use ";".
detect_lineage_separator <- function(lineages) {
  if (!length(lineages)) return("|")
  n_pipe <- sum(grepl("|", lineages, fixed = TRUE))
  n_semi <- sum(grepl(";", lineages, fixed = TRUE))
  if (n_pipe >= n_semi && n_pipe > 0) "|" else if (n_semi > 0) ";" else "|"
}

#' Strip rank prefixes such as "k__" or "s__" and tidy separators.
clean_taxon_name <- function(x) {
  x <- sub("^[dkpcofgst]__", "", trimws(x))
  x <- gsub("_", " ", x)
  x[!nzchar(x)] <- NA_character_
  x
}

#' Keep only the deepest (leaf) lineages in a hierarchical profile.
#'
#' MetaPhlAn-style tables repeat each clade at every rank, so the same reads
#' appear in a kingdom row, a phylum row, and so on. Summing those rows counts
#' the same organism many times over. A row is a leaf when no other row extends
#' it, which selects one non-overlapping set regardless of how deep individual
#' lineages go.
leaf_lineages <- function(lineages, sep = "|") {
  if (!length(lineages)) return(logical(0))
  parts <- strsplit(lineages, sep, fixed = TRUE)
  ancestors <- new.env(hash = TRUE, parent = emptyenv())
  for (pp in parts) {
    if (length(pp) < 2) next
    for (i in seq_len(length(pp) - 1)) {
      assign(paste(pp[seq_len(i)], collapse = sep), TRUE, envir = ancestors)
    }
  }
  !vapply(lineages, function(l) exists(l, envir = ancestors, inherits = FALSE),
          logical(1), USE.NAMES = FALSE)
}

#' Are these values already relative abundances rather than read counts?
#'
#' Profilers such as MetaPhlAn report percentages, so depth-based QC, singleton
#' removal, rarefaction and the count-based richness estimators do not apply.
looks_like_relative_abundance <- function(mat) {
  if (!length(mat)) return(FALSE)
  if (any(mat < 0, na.rm = TRUE)) return(FALSE)
  totals <- colSums(mat, na.rm = TRUE)
  totals <- totals[totals > 0]
  if (!length(totals)) return(FALSE)

  # Every column landing on 100 (or on 1) is decisive on its own: read counts
  # would not do that across a whole study.
  on_100 <- all(abs(totals - 100) < 0.5)
  on_1   <- all(abs(totals - 1) < 0.01)
  if (on_100 || on_1) return(TRUE)

  # Otherwise the totals alone are ambiguous, because MetaPhlAn leaves an
  # unclassified fraction out and a real profile often sums to the nineties.
  # Fractional values then settle it: read counts are whole numbers.
  fractional <- any(mat %% 1 != 0, na.rm = TRUE)
  if (!fractional) return(FALSE)

  as_proportion <- all(totals > 0.5 & totals <= 1.01)
  as_percentage <- all(totals > 50 & totals <= 100.5)
  as_proportion || as_percentage
}

#' Read a profiler export, skipping any banner lines above the header.
#'
#' MetaPhlAn writes a version banner such as "#mpa_v30_CHOCOPhlAn_201901" on the
#' first line, which would otherwise be read as the header and push the real
#' header down into the data. The header is the first line that actually
#' contains the separator; anything above it is a banner.
read_profile_file <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (!length(lines)) stop("The taxonomic profile is empty.", call. = FALSE)

  probe  <- seq_len(min(20, length(lines)))
  n_tab  <- vapply(lines[probe],
                   function(l) lengths(regmatches(l, gregexpr("\t", l))), integer(1))
  n_com  <- vapply(lines[probe],
                   function(l) lengths(regmatches(l, gregexpr(",", l))), integer(1))
  sep    <- if (max(n_tab) >= max(n_com)) "\t" else ","
  counts <- if (identical(sep, "\t")) n_tab else n_com

  header_idx <- which(counts > 0)[1]
  if (is.na(header_idx)) {
    stop("Could not find a header row with more than one column in this file.",
         call. = FALSE)
  }

  df <- utils::read.delim(
    text = paste(lines[seq(header_idx, length(lines))], collapse = "\n"),
    sep = sep, check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
  # Some exports mark the header itself with a leading '#'.
  names(df)[1] <- sub("^#\\s*", "", names(df)[1])
  df
}

#' Parse a shotgun taxonomic profile into abundance and taxonomy matrices.
#'
#' Accepts a data frame whose first column (or row names) holds lineage strings
#' and whose remaining columns are samples. Returns the same pieces an amplicon
#' upload provides, so every downstream stage is unchanged.
parse_taxonomic_profile <- function(df, lineage_col = NULL) {
  if (is.null(df) || !nrow(df)) stop("The taxonomic profile is empty.", call. = FALSE)

  if (is.null(lineage_col)) {
    first_is_lineage <- !is.numeric(df[[1]])
    lineages <- if (first_is_lineage) as.character(df[[1]]) else rownames(df)
    abund    <- if (first_is_lineage) df[, -1, drop = FALSE] else df
  } else {
    lineages <- as.character(df[[lineage_col]])
    abund    <- df[, setdiff(names(df), lineage_col), drop = FALSE]
  }

  if (is.null(lineages) || !length(lineages) || all(is.na(lineages))) {
    stop("Could not find lineage strings. Expected a first column such as ",
         "'clade_name' holding entries like 'k__Bacteria|p__Firmicutes'.",
         call. = FALSE)
  }

  numeric_cols <- vapply(abund, is.numeric, logical(1))
  if (!any(numeric_cols)) {
    stop("No numeric sample columns found in the taxonomic profile.", call. = FALSE)
  }
  # Drop profiler bookkeeping columns such as NCBI taxid, which are numeric but
  # are not samples.
  drop_names <- grepl("tax(onomy)?_?id|ncbi|clade_taxid", names(abund), ignore.case = TRUE)
  abund <- abund[, numeric_cols & !drop_names, drop = FALSE]
  if (!ncol(abund)) stop("No sample columns remain in the taxonomic profile.", call. = FALSE)

  sep  <- detect_lineage_separator(lineages)
  keep <- leaf_lineages(lineages, sep = sep)
  n_collapsed <- sum(!keep)
  lineages <- lineages[keep]
  abund    <- abund[keep, , drop = FALSE]
  if (!nrow(abund)) stop("No taxa remained after collapsing the lineage hierarchy.", call. = FALSE)

  parts   <- strsplit(lineages, sep, fixed = TRUE)
  depth   <- max(lengths(parts))
  n_ranks <- max(length(TAX_RANKS), depth)
  tax <- t(vapply(parts, function(pp) {
    out <- rep(NA_character_, n_ranks)
    pp  <- clean_taxon_name(pp)
    if (length(pp)) out[seq_along(pp)] <- pp
    out
  }, character(n_ranks)))
  colnames(tax) <- if (n_ranks <= length(TAX_RANKS)) {
    TAX_RANKS[seq_len(n_ranks)]
  } else {
    c(TAX_RANKS, paste0("Rank", seq_len(n_ranks - length(TAX_RANKS))))
  }

  # Name each taxon by its deepest resolved rank, kept unique for phyloseq.
  leaf_name <- apply(tax, 1, function(r) {
    r <- r[!is.na(r)]
    if (length(r)) r[length(r)] else NA_character_
  })
  leaf_name[is.na(leaf_name)] <- "Unclassified"
  ids <- make.unique(as.character(leaf_name), sep = "_")

  otu <- as.matrix(abund)
  mode(otu) <- "numeric"
  otu[is.na(otu)] <- 0
  rownames(otu) <- ids
  rownames(tax) <- ids

  list(otu = otu, tax = tax,
       is_relative      = looks_like_relative_abundance(otu),
       n_ranks          = n_ranks,
       separator        = sep,
       n_rows_collapsed = n_collapsed)
}

# ── Input coercion ───────────────────────────────────────────────────────────

input_num <- function(x, default) {
  if (is.null(x) || !length(x) || !is.finite(suppressWarnings(as.numeric(x)[1]))) default
  else as.numeric(x)[1]
}

input_chr <- function(x, default) {
  if (is.null(x) || !length(x) || is.na(x[1]) || !nzchar(x[1])) default else as.character(x)[1]
}

#' Geometric mean over positive values only.
#'
#' DESeq2's default size-factor estimation fails when every taxon has at least
#' one zero, which is the norm for microbiome data. Restricting the product to
#' positive counts is the standard phyloseq workaround (Demo 10, ps10).
gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

#' Format a p-value for prose without ever rounding it to zero.
mfg_fmt_p <- function(p, digits = 3) {
  ifelse(is.na(p), "NA",
    ifelse(p < 1e-16, "< 1e-16",
      ifelse(p < 0.001, format(p, scientific = TRUE, digits = 2),
             formatC(p, format = "f", digits = digits))))
}
