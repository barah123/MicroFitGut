# =============================================================================
# MicroFitGut — 01-intake.R
#
# Getting data into a phyloseq object, and refusing to proceed when it is not
# what it claims to be. Every check here exists because skipping it produces an
# analysis that runs cleanly and answers the wrong question.
#
# Accepted inputs:
#   - a saved phyloseq object            .RDS
#   - an OTU/ASV + taxonomy + metadata triplet   .csv / .txt / .tsv
#   - a shotgun taxonomic profile        MetaPhlAn, Kraken2/Bracken
#   - a functional profile               PICRUSt2, HUMAnN pathway abundance
#   - a BIOM table                       via phyloseq::import_biom
#
# Requires: mfg_require(c("intake"))  and utils.R
# =============================================================================

# ── Which dataset is the one to analyse ──────────────────────────────────────
#
# A data directory rarely holds one dataset. It holds a full study object, two
# teaching subsets, a triplet that is a third subset, and last month's rerun.
# Every one of them loads without complaint and every one gives a different
# answer, so the choice between them is an analytical decision and gets recorded
# as one.

MFG_DATASET_PATTERNS <- c(
  "\\.rds$", "\\.biom$", "\\.qza$", "\\.tsv$", "\\.txt$", "\\.csv$"
)

#' Inventory the candidate inputs in a directory before choosing one.
#'
#' Returns every file that could be a microbiome input, with the size and
#' modification time that usually distinguish a full study from a subset. This
#' does not choose; it makes the choice visible so that it has to be made.
mfg_input_inventory <- function(dir = "data", recursive = TRUE) {
  if (!dir.exists(dir)) stop("Not a directory: ", dir, call. = FALSE)

  files <- list.files(dir, recursive = recursive, full.names = TRUE)
  files <- files[file.exists(files) & !dir.exists(files)]
  keep  <- Reduce(`|`, lapply(MFG_DATASET_PATTERNS, function(p)
    grepl(p, files, ignore.case = TRUE)), init = rep(FALSE, length(files)))
  files <- files[keep]

  if (!length(files)) {
    mfg_log("intake", "inputs_inventoried", list(dir = dir, n_files = 0))
    message("No candidate microbiome inputs found under ", dir, ".")
    return(invisible(data.frame()))
  }

  info <- file.info(files)
  inv <- data.frame(
    file     = basename(files),
    path     = files,
    format   = tolower(tools::file_ext(files)),
    bytes    = as.numeric(info$size),
    modified = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
  inv <- inv[order(-inv$bytes), ]
  rownames(inv) <- NULL

  # A self-contained object is a candidate dataset on its own. A .csv is usually
  # one third of a triplet, so it is listed but not counted as a rival dataset.
  standalone <- sum(inv$format %in% c("rds", "biom", "qza"))
  MFG_LOG$sources$inventory  <- inv
  MFG_LOG$sources$standalone <- standalone

  mfg_log("intake", "inputs_inventoried",
          list(dir = dir, n_files = nrow(inv), n_standalone = standalone))

  if (standalone > 1) {
    message("\n", standalone, " self-contained datasets found under ", dir,
            ". They will give different answers.\n",
            "Call mfg_declare_authoritative(<path>, reason = \"...\") ",
            "before loading, or the report will record the choice as ",
            "undocumented.\n")
  }
  inv
}

#' Record which input is the authoritative one, and why.
#'
#' The reason is required and is quoted in the report. "It was the largest file"
#' is a reason; an undocumented pick between four files is not.
mfg_declare_authoritative <- function(path, reason) {
  if (missing(reason) || !is.character(reason) || !nzchar(trimws(reason))) {
    stop("mfg_declare_authoritative() needs a reason. Name what makes this file ",
         "the right one: the study it covers, the version, the documentation ",
         "that says so.", call. = FALSE)
  }
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)

  MFG_LOG$sources$authoritative <- list(
    path = normalizePath(path, mustWork = FALSE),
    file = basename(path),
    reason = trimws(reason)
  )
  mfg_log("intake", "authoritative_source_declared",
          list(file = basename(path), reason = trimws(reason)))
  invisible(MFG_LOG$sources$authoritative)
}

#' What was declared, or NULL. Read by assemble_summary() and the report.
mfg_authoritative_source <- function() MFG_LOG$sources$authoritative

# ── Reading a delimited table ────────────────────────────────────────────────

#' Read a table whose first column holds row names, whatever the delimiter.
#'
#' The course material mixes .csv and tab-delimited .txt freely, sometimes with
#' the same content, so the delimiter is detected rather than assumed.
mfg_read_table <- function(path, row_names = 1) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  # Detect the delimiter from the HEADER, not from line 1. QIIME and BIOM
  # exports begin "# Constructed from biom file", a comment containing neither
  # a tab nor a comma. Reading that line to choose the separator picks comma for
  # a tab-separated file, and the result is a data frame with ZERO columns and
  # no error: every downstream step then fails somewhere unrelated.
  head_lines <- readLines(path, n = 50, warn = FALSE)
  head_lines <- sub("^\ufeff", "", head_lines)          # strip a UTF-8 BOM
  cand <- head_lines[nzchar(trimws(head_lines))]
  # The header is the first non-blank line that is not a pure comment, or the
  # last comment line when the header itself is commented ("#OTU ID\t...").
  is_comment <- grepl("^#", cand)
  hdr <- if (any(!is_comment)) {
    commented_header <- which(is_comment & grepl("[\t,]", cand))
    if (length(commented_header)) cand[max(commented_header)] else cand[which(!is_comment)[1]]
  } else cand[1]
  ntab <- lengths(regmatches(hdr, gregexpr("\t", hdr)))
  ncom <- lengths(regmatches(hdr, gregexpr(",",  hdr)))
  sep <- if (ntab > ncom) "\t" else ","
  # A commented header must be kept, so only skip comment lines above it.
  skip <- if (grepl("^#", hdr)) which(cand == hdr)[1] - 1L else 0L
  df <- utils::read.delim(path, sep = sep, row.names = row_names, header = TRUE,
                          check.names = FALSE, stringsAsFactors = FALSE,
                          comment.char = "", skip = skip,
                          na.strings = c("NA", "", "NaN"))
  names(df) <- sub("^\ufeff", "", names(df))
  if (ncol(df) == 0) {
    stop(sprintf(paste("Read %s and got %d rows but no columns. The delimiter",
      "was detected as %s from the header line. Check the file's actual",
      "separator and comment lines."), basename(path), nrow(df),
      if (identical(sep, "\t")) "tab" else "comma"), call. = FALSE)
  }
  df
}

# ── Orientation ──────────────────────────────────────────────────────────────

#' Decide whether taxa are rows or columns in an abundance table.
#'
#' Getting this backwards is silent and catastrophic: every per-sample statistic
#' becomes a per-taxon statistic and nothing errors. Decided by matching against
#' the names in the taxonomy and metadata, not by guessing from the shape, since
#' a study can easily have more samples than taxa or the reverse.
mfg_detect_orientation <- function(abund, tax_ids = NULL, sample_ids = NULL) {
  rn <- rownames(abund); cn <- colnames(abund)
  score_rows_are_taxa <- 0L
  score_cols_are_taxa <- 0L

  if (!is.null(tax_ids)) {
    score_rows_are_taxa <- score_rows_are_taxa + sum(rn %in% tax_ids)
    score_cols_are_taxa <- score_cols_are_taxa + sum(cn %in% tax_ids)
  }
  if (!is.null(sample_ids)) {
    score_rows_are_taxa <- score_rows_are_taxa + sum(cn %in% sample_ids)
    score_cols_are_taxa <- score_cols_are_taxa + sum(rn %in% sample_ids)
  }

  if (score_rows_are_taxa == score_cols_are_taxa) {
    return(list(taxa_are_rows = NA, evidence = "ambiguous",
                rows_score = score_rows_are_taxa, cols_score = score_cols_are_taxa))
  }
  list(taxa_are_rows = score_rows_are_taxa > score_cols_are_taxa,
       evidence = "matched against taxonomy and metadata names",
       rows_score = score_rows_are_taxa, cols_score = score_cols_are_taxa)
}

# ── Building a phyloseq object ────────────────────────────────────────────────

#' Build a phyloseq object from separate tables.
#'
#' Follows Demo 3 but adds the orientation check and the sample-ID intersection
#' report, and refuses a random tree by default. `tree` accepts a path to a
#' Newick file or a phylo object; `tree_is_real` records whether it was actually
#' estimated from the sequence data, which every phylogenetic metric checks.
build_phyloseq <- function(abund, tax = NULL, meta = NULL,
                           tree = NULL, seqs = NULL,
                           taxa_are_rows = NULL, tree_is_real = NULL) {

  # Record the paths before reading turns them into data frames and the file
  # they came from is no longer recoverable from the object.
  for (p in list(abund, tax, meta, tree, seqs)) {
    if (is.character(p) && length(p) == 1 && file.exists(p)) mfg_record_input(p, "raw")
  }

  if (is.character(abund)) abund <- mfg_read_table(abund)
  if (is.character(tax) && length(tax) == 1)  tax  <- mfg_read_table(tax)
  if (is.character(meta) && length(meta) == 1) meta <- mfg_read_table(meta)

  abund_m <- as.matrix(abund)
  mode(abund_m) <- "numeric"

  tax_ids    <- if (!is.null(tax))  rownames(tax)  else NULL
  sample_ids <- if (!is.null(meta)) rownames(meta) else NULL

  if (is.null(taxa_are_rows)) {
    det <- mfg_detect_orientation(abund_m, tax_ids, sample_ids)
    if (is.na(det$taxa_are_rows)) {
      stop("Cannot determine table orientation: row and column names match the ",
           "taxonomy/metadata equally well (scores ", det$rows_score, " vs ",
           det$cols_score, ").\nPass taxa_are_rows = TRUE or FALSE explicitly.",
           call. = FALSE)
    }
    taxa_are_rows <- det$taxa_are_rows
    mfg_log("intake", "orientation_detected",
            list(taxa_are_rows = taxa_are_rows, evidence = det$evidence))
  }

  parts <- list(phyloseq::otu_table(abund_m, taxa_are_rows = taxa_are_rows))
  if (!is.null(tax))  parts <- c(parts, list(phyloseq::tax_table(as.matrix(tax))))
  if (!is.null(meta)) parts <- c(parts, list(phyloseq::sample_data(as.data.frame(meta))))

  if (!is.null(tree)) {
    if (is.character(tree)) tree <- ape::read.tree(tree)
    parts <- c(parts, list(phyloseq::phy_tree(tree)))
  }
  if (!is.null(seqs)) {
    if (is.character(seqs) && length(seqs) == 1 && file.exists(seqs)) {
      seqs <- Biostrings::readDNAStringSet(seqs)
    }
    parts <- c(parts, list(phyloseq::refseq(seqs)))
  }

  ps <- do.call(phyloseq::phyloseq, parts)

  if (!is.null(tree)) {
    if (is.null(tree_is_real)) {
      warning("A tree was supplied but not marked as real or placeholder. ",
              "Phylogenetic metrics (UniFrac, Faith's PD) will be withheld ",
              "until you call mfg_mark_tree_real(ps, TRUE). See reference/01.",
              call. = FALSE)
    } else {
      ps <- mfg_mark_tree_real(ps, tree_is_real)
    }
  }

  mfg_log("intake", "phyloseq_built",
          list(taxa = phyloseq::ntaxa(ps), samples = phyloseq::nsamples(ps),
               has_tree = !is.null(tree), has_seqs = !is.null(seqs),
               tree_is_real = tree_is_real %||% NA))
  ps
}

#' Load a shotgun taxonomic profile as a phyloseq object.
#'
#' Wraps read_profile_file + parse_taxonomic_profile from utils.R, attaches
#' metadata, and records that the values are relative abundances when they are —
#' which blocks rarefaction and the count-based richness estimators downstream.
build_phyloseq_from_profile <- function(profile_path, meta = NULL, lineage_col = NULL) {
  # A taxonomic profile has already been through classification and summarisation,
  # so it is processed input however it is labelled upstream.
  mfg_record_input(profile_path, "processed")
  if (is.character(meta) && length(meta) == 1 && file.exists(meta)) {
    mfg_record_input(meta, "raw")
  }
  df  <- read_profile_file(profile_path)
  pp  <- parse_taxonomic_profile(df, lineage_col = lineage_col)

  parts <- list(
    phyloseq::otu_table(pp$otu, taxa_are_rows = TRUE),
    phyloseq::tax_table(pp$tax)
  )
  if (!is.null(meta)) {
    if (is.character(meta) && length(meta) == 1) meta <- mfg_read_table(meta)
    parts <- c(parts, list(phyloseq::sample_data(as.data.frame(meta))))
  }
  ps <- do.call(phyloseq::phyloseq, parts)
  attr(ps, "mfg_is_relative") <- pp$is_relative
  mfg_registry_set(ps, "is_relative", pp$is_relative)

  mfg_log("intake", "profile_parsed",
          list(file = basename(profile_path), separator = pp$separator,
               taxa = phyloseq::ntaxa(ps), samples = phyloseq::nsamples(ps),
               rows_collapsed = pp$n_rows_collapsed, is_relative = pp$is_relative))
  ps
}
#' Build a phyloseq object from a SummarizedExperiment.
#'
#' curatedMetagenomicData's curatedMetagenomicData() and returnSamples() return
#' a (Tree)SummarizedExperiment, while the rest of this toolkit works on
#' phyloseq. This is the bridge, and it does three things beyond reshaping.
#'
#' It carries colData across as sample data, dropping any list column, because
#' phyloseq's sample_data() requires a plain data frame and will fail opaquely
#' on curated metadata that holds one.
#'
#' It detects whether the assay holds counts or proportions and records the
#' answer in the normalization registry, rather than leaving a later method to
#' assume. For functional assays this is the decision that determines which
#' differential abundance methods can run at all: ALDEx2 and DESeq2 need integer
#' counts, LinDA does not.
#'
#' It marks functional assays, whose rows are pathways or gene families rather
#' than taxa, so downstream code does not try to read a lineage out of them.
build_phyloseq_from_se <- function(se, assay_name = NULL, meta = NULL,
                                   normalization = NULL, functional = NULL) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Needs SummarizedExperiment.\n",
         '  BiocManager::install("SummarizedExperiment")', call. = FALSE)
  }
  assays_available <- SummarizedExperiment::assayNames(se)
  an <- assay_name %||% (if (length(assays_available)) assays_available[1] else 1L)
  mat <- as.matrix(SummarizedExperiment::assay(se, an))   # features x samples

  cd <- as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)
  # A list column is legal in a DataFrame and fatal in sample_data().
  is_list_col <- vapply(cd, function(x) is.list(x) && !is.data.frame(x), logical(1))
  if (any(is_list_col)) {
    dropped <- names(cd)[is_list_col]
    cd <- cd[, !is_list_col, drop = FALSE]
    message("Dropped list column(s) from colData: ", paste(dropped, collapse = ", "))
  }
  if (!is.null(meta)) {
    if (is.character(meta) && length(meta) == 1) {
      mfg_record_input(meta, "raw")
      meta <- utils::read.csv(meta, header = TRUE, row.names = 1, na.strings = "NA")
    }
    common <- intersect(rownames(cd), rownames(meta))
    if (!length(common)) stop("Supplied metadata shares no sample ids with the object.",
                              call. = FALSE)
    cd <- cbind(cd[common, , drop = FALSE],
                meta[common, setdiff(names(meta), names(cd)), drop = FALSE])
    mat <- mat[, common, drop = FALSE]
  }
  if (!identical(colnames(mat), rownames(cd))) {
    common <- intersect(colnames(mat), rownames(cd))
    if (!length(common)) stop("Assay columns and colData rows share no ids.", call. = FALSE)
    mat <- mat[, common, drop = FALSE]; cd <- cd[common, , drop = FALSE]
  }

  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(mat, taxa_are_rows = TRUE),
    phyloseq::sample_data(cd))

  # Rows that carry a HUMAnN-style lineage separator, or the unmapped and
  # unintegrated sentinels, are functional rather than taxonomic.
  if (is.null(functional)) {
    functional <- any(grepl("^(UNMAPPED|UNINTEGRATED)", rownames(mat))) ||
      any(grepl("\\|", rownames(mat))) ||
      any(grepl("PWY|^UniRef", rownames(mat)))
  }
  if (isTRUE(functional)) {
    attr(ps, "mfg_is_functional") <- TRUE
    # Side effect only: mfg_registry_set() returns the value, not the object.
    mfg_registry_set(ps, "is_functional", TRUE)
  }

  detected <- if (looks_like_relative_abundance(mat)) "tss" else "raw"
  norm <- normalization %||% detected
  ps <- mfg_set_normalization(ps, norm)

  integerish <- all(abs(mat - round(mat)) < 1e-8, na.rm = TRUE)
  mfg_log("intake", "from_summarized_experiment", list(
    assay = as.character(an),
    assays_available = paste(assays_available, collapse = ", "),
    n_features = nrow(mat), n_samples = ncol(mat),
    functional = isTRUE(functional),
    normalization_detected = detected,
    normalization_set = norm,
    integer_valued = integerish,
    sample_sum_median = stats::median(colSums(mat), na.rm = TRUE)))

  message(sprintf(
    "SummarizedExperiment -> phyloseq: %d features x %d samples, assay '%s'.",
    nrow(mat), ncol(mat), as.character(an)))
  message(sprintf(
    "  values look like %s and are %sinteger-valued; normalization set to '%s'.",
    if (detected == "tss") "proportions" else "counts",
    if (integerish) "" else "not ", norm))
  if (!integerish) {
    message("  note: ALDEx2 and DESeq2 require integer counts and cannot run on this assay.")
  }
  ps
}


#' Load a functional profile (PICRUSt2 / HUMAnN) as a phyloseq object.
#'
#' The pathway table takes the place of the OTU table and the pathway
#' description table takes the place of the taxonomy, giving a `Pathways` rank
#' (Demo 6). Diversity metrics are meaningful on this only with care — see
#' reference/08 — but composition, DA and ordination all work unchanged.
build_phyloseq_functional <- function(abund_path, pathway_tax_path, meta = NULL,
                                      skip = 1) {
  mfg_record_input(abund_path, "processed")
  mfg_record_input(pathway_tax_path, "raw")
  if (is.character(meta) && length(meta) == 1 && file.exists(meta)) {
    mfg_record_input(meta, "raw")
  }
  otu  <- utils::read.delim(abund_path, skip = skip, row.names = 1,
                            check.names = FALSE, stringsAsFactors = FALSE)
  taxa <- utils::read.delim(pathway_tax_path, row.names = 1,
                            check.names = FALSE, stringsAsFactors = FALSE)

  parts <- list(
    phyloseq::otu_table(as.matrix(otu), taxa_are_rows = TRUE),
    phyloseq::tax_table(as.matrix(taxa))
  )
  if (!is.null(meta)) {
    if (is.character(meta) && length(meta) == 1) {
      meta <- utils::read.csv(meta, header = TRUE, row.names = 1, na.strings = "NA")
    }
    parts <- c(parts, list(phyloseq::sample_data(as.data.frame(meta))))
  }
  ps <- do.call(phyloseq::phyloseq, parts)
  attr(ps, "mfg_is_functional") <- TRUE
  mfg_registry_set(ps, "is_functional", TRUE)

  mfg_log("intake", "functional_profile_loaded",
          list(file = basename(abund_path), features = phyloseq::ntaxa(ps),
               samples = phyloseq::nsamples(ps),
               ranks = paste(phyloseq::rank_names(ps), collapse = ",")))
  ps
}

#' Load whatever was handed over, dispatching on file type.
mfg_load <- function(path, meta = NULL, ...) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("rds")) {
    ps <- readRDS(path)
    if (!methods::is(ps, "phyloseq")) {
      stop(path, " holds a ", class(ps)[1], ", not a phyloseq object.", call. = FALSE)
    }
    # A saved phyloseq object has been built, and usually filtered, by whoever
    # saved it. Calling it raw would misdescribe it.
    mfg_record_input(path, "processed", taxa = phyloseq::ntaxa(ps),
                     samples = phyloseq::nsamples(ps))
    mfg_log("intake", "rds_loaded",
            list(file = basename(path), taxa = phyloseq::ntaxa(ps),
                 samples = phyloseq::nsamples(ps)))
    return(ps)
  }
  if (ext %in% c("biom")) {
    ps <- phyloseq::import_biom(path, ...)
    mfg_record_input(path, "processed", taxa = phyloseq::ntaxa(ps),
                     samples = phyloseq::nsamples(ps))
    mfg_log("intake", "biom_loaded",
            list(file = basename(path), taxa = phyloseq::ntaxa(ps),
                 samples = phyloseq::nsamples(ps)))
    return(ps)
  }
  stop("mfg_load handles .RDS and .biom. For a table triplet use ",
       "build_phyloseq(); for a shotgun profile use build_phyloseq_from_profile().",
       call. = FALSE)
}

# ── Validation ───────────────────────────────────────────────────────────────

#' The mandatory pre-analysis check.
#' Read an HDF5 BIOM that phyloseq's import_biom() cannot open.
#'
#' `phyloseq::import_biom()` assumes a BIOM carries observation metadata. When
#' the `/observation/metadata` group exists but holds no datasets, which is what
#' QIIME 2 writes for a feature table exported without taxonomy, it fails with
#'
#'     length of 'dimnames' [2] not equal to array extent
#'
#' That message names neither the file nor the missing taxonomy, so the usual
#' response is to suspect a corrupt download. The table itself is fine. This
#' reader builds the object from the sparse matrix directly and attaches
#' taxonomy and sample data only when they are actually present, so a table
#' without taxonomy loads as a table without taxonomy rather than failing.
#'
#' Requires rhdf5. Returns a phyloseq object with whatever the file contains.
mfg_read_biom_hdf5 <- function(path) {
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("Reading an HDF5 BIOM needs rhdf5:\n",
         "  BiocManager::install(\"rhdf5\")", call. = FALSE)
  }
  h5 <- function(x) rhdf5::h5read(path, x)
  ls <- rhdf5::h5ls(path)
  obs  <- as.character(h5("/observation/ids"))
  samp <- as.character(h5("/sample/ids"))
  # BIOM stores CSR by observation: indptr walks rows, indices give columns.
  data    <- as.numeric(h5("/observation/matrix/data"))
  indices <- as.integer(h5("/observation/matrix/indices"))
  indptr  <- as.integer(h5("/observation/matrix/indptr"))
  m <- matrix(0, nrow = length(obs), ncol = length(samp),
              dimnames = list(obs, samp))
  for (i in seq_along(obs)) {
    lo <- indptr[i] + 1L; hi <- indptr[i + 1L]
    if (hi >= lo) m[i, indices[lo:hi] + 1L] <- data[lo:hi]
  }
  parts <- list(phyloseq::otu_table(m, taxa_are_rows = TRUE))
  has_tax  <- any(ls$group == "/observation/metadata" & ls$otype == "H5I_DATASET")
  has_meta <- any(ls$group == "/sample/metadata"      & ls$otype == "H5I_DATASET")
  if (has_tax) {
    tx <- h5("/observation/metadata/taxonomy")
    tm <- if (is.matrix(tx)) t(tx) else as.matrix(tx)
    rownames(tm) <- obs
    parts <- c(parts, list(phyloseq::tax_table(tm)))
  }
  if (has_meta) {
    nm <- ls$name[ls$group == "/sample/metadata" & ls$otype == "H5I_DATASET"]
    md <- as.data.frame(lapply(setNames(nm, nm),
                               function(k) as.character(h5(paste0("/sample/metadata/", k)))),
                        stringsAsFactors = FALSE)
    rownames(md) <- samp
    parts <- c(parts, list(phyloseq::sample_data(md)))
  }
  ps <- do.call(phyloseq::phyloseq, parts)
  mfg_log("intake", "biom_hdf5_read", list(
    file = basename(path), taxa = length(obs), samples = length(samp),
    taxonomy = has_tax, sample_metadata = has_meta,
    generated_by = tryCatch(as.character(rhdf5::h5readAttributes(path, "/")[["generated-by"]])[1],
                            error = function(e) NA_character_)))
  if (!has_tax) {
    warning(sprintf(paste("%s carries no taxonomy. Feature IDs are all that",
      "identify a row, so nothing can be agglomerated or matched by name.",
      "Assign taxonomy upstream, or restrict the analysis to diversity and",
      "ordination."), basename(path)), call. = FALSE)
  }
  if (!has_meta) {
    warning(sprintf(paste("%s carries no sample metadata. Any grouping must",
      "come from another source."), basename(path)), call. = FALSE)
  }
  ps
}

#' Split a single-column lineage into proper rank columns.
#'
#' `phyloseq::import_biom()` on a QIIME-style BIOM yields one column, usually
#' named Rank1, holding the whole lineage as `k__Bacteria;p__Firmicutes;...`.
#' Everything downstream that agglomerates, filters or matches taxa by rank sees
#' a taxonomy table with no ranks in it and silently has nothing to work with —
#' `tax_glom(ps, "Genus")` errors, and taxon matching against a published genus
#' name finds nothing. This is the commonest shape of deposited 16S data, so it
#' is handled rather than left to the caller.
#'
#' Accepts `;` or `|` separators and `k__`/`d__` prefixes, and pads short
#' lineages with NA rather than recycling, so a lineage that stops at Family does
#' not acquire a fabricated Genus.
mfg_split_lineage <- function(ps, sep = NULL, ranks = TAX_RANKS) {
  tt <- phyloseq::tax_table(ps, errorIfNULL = FALSE)
  if (is.null(tt)) stop("No taxonomy table to split.", call. = FALSE)
  m <- as(tt, "matrix")
  if (ncol(m) > 1) {
    # Already split, but import_biom names the columns Rank1..Rank7 rather than
    # Kingdom..Species. Same practical failure as an unsplit lineage: there is no
    # "Genus" column, so tax_glom(ps, "Genus") errors and taxon matching by rank
    # finds nothing. Rename when the count matches the rank schema; leave alone
    # when it does not, since guessing an alignment would be worse than failing.
    if (all(grepl("^Rank[0-9]+$", colnames(m))) && ncol(m) == length(ranks)) {
      old <- colnames(m)
      colnames(m) <- ranks
      m[] <- sub("^[dkpcofgst]__", "", trimws(m))
      m[!nzchar(m) | m %in% c("NA", "unidentified")] <- NA_character_
      phyloseq::tax_table(ps) <- phyloseq::tax_table(m)
      mfg_log("intake", "ranks_renamed", list(
        from = paste(old, collapse = ","), to = paste(ranks, collapse = ","),
        resolved_to_genus = sum(!is.na(m[, "Genus"]))))
      message("Renamed ", paste(old, collapse = "/"), " to ",
              paste(ranks, collapse = "/"), ".")
      return(ps)
    }
    message("Taxonomy already has ", ncol(m), " columns; nothing to split.")
    return(ps)
  }
  lin <- as.character(m[, 1])
  if (is.null(sep)) {
    sep <- if (mean(grepl(";", lin, fixed = TRUE)) >= mean(grepl("|", lin, fixed = TRUE))) ";" else "|"
  }
  parts <- strsplit(lin, sep, fixed = TRUE)
  n <- length(ranks)
  out <- t(vapply(parts, function(p) {
    p <- trimws(sub("^[dkpcofgst]__", "", trimws(p)))
    p[!nzchar(p)] <- NA_character_
    length(p) <- n                      # pads with NA, never recycles
    p
  }, character(n)))
  colnames(out) <- ranks
  rownames(out) <- rownames(m)
  phyloseq::tax_table(ps) <- phyloseq::tax_table(out)
  mfg_log("intake", "lineage_split", list(
    sep = sep, ranks = paste(ranks, collapse = ","),
    n_taxa = nrow(out),
    resolved_to_genus = sum(!is.na(out[, "Genus"])),
    resolved_to_species = sum(!is.na(out[, "Species"]))))
  ps
}

#'
#' Returns a structured report and, unless `strict = FALSE`, stops on anything
#' that would invalidate a downstream result. Nothing in MicroFitGut runs before
#' this passes. Each check maps to a failure mode in reference/01 and
#' reference/11.
validate_inputs <- function(ps, group_var = NULL, subject_var = NULL,
                            strict = TRUE, min_depth_warn = 1000) {

  problems <- character()
  warnings_ <- character()
  rep <- list()

  # --- 1. Slots present -----------------------------------------------------
  rep$has_otu  <- !is.null(phyloseq::otu_table(ps, errorIfNULL = FALSE))
  rep$has_tax  <- !is.null(phyloseq::tax_table(ps, errorIfNULL = FALSE))
  rep$has_meta <- !is.null(phyloseq::sample_data(ps, errorIfNULL = FALSE))
  rep$has_tree <- !is.null(phyloseq::phy_tree(ps, errorIfNULL = FALSE))
  rep$has_seqs <- !is.null(phyloseq::refseq(ps, errorIfNULL = FALSE))
  if (!rep$has_otu) problems <- c(problems, "No abundance table (otu_table) present.")

  # --- 2. Dimensions --------------------------------------------------------
  rep$n_taxa    <- phyloseq::ntaxa(ps)
  rep$n_samples <- phyloseq::nsamples(ps)
  rep$ranks     <- if (rep$has_tax) phyloseq::rank_names(ps) else character()
  if (rep$n_samples < 3) {
    problems <- c(problems, sprintf("Only %d sample(s). No group comparison is possible.",
                                    rep$n_samples))
  }

  # --- 3. Sample IDs match across tables ------------------------------------
  # phyloseq intersects silently on construction, so a mismatch shows up as
  # missing samples rather than an error. Reported explicitly here.
  if (rep$has_meta) {
    otu_samples  <- phyloseq::sample_names(ps)
    meta_samples <- rownames(mfg_meta(ps))
    rep$samples_in_both   <- length(intersect(otu_samples, meta_samples))
    rep$samples_otu_only  <- setdiff(otu_samples, meta_samples)
    rep$samples_meta_only <- setdiff(meta_samples, otu_samples)
    if (length(rep$samples_otu_only) || length(rep$samples_meta_only)) {
      problems <- c(problems, sprintf(
        "Sample IDs do not match: %d in abundance table only, %d in metadata only.",
        length(rep$samples_otu_only), length(rep$samples_meta_only)))
    }
  }

  # --- 4. Taxa IDs match ----------------------------------------------------
  if (rep$has_tax) {
    otu_taxa <- phyloseq::taxa_names(ps)
    tax_taxa <- rownames(mfg_tax(ps))
    rep$taxa_otu_only <- setdiff(otu_taxa, tax_taxa)
    rep$taxa_tax_only <- setdiff(tax_taxa, otu_taxa)
    if (length(rep$taxa_otu_only) || length(rep$taxa_tax_only)) {
      problems <- c(problems, sprintf(
        "Taxa IDs do not match: %d in abundance table only, %d in taxonomy only.",
        length(rep$taxa_otu_only), length(rep$taxa_tax_only)))
    }
  }

  # --- 5. Are these counts or proportions? ---------------------------------
  mat <- as(phyloseq::otu_table(ps), "matrix")
  rep$is_relative <- attr(ps, "mfg_is_relative") %||%
    mfg_registry_get(ps, "is_relative") %||% looks_like_relative_abundance(mat)
  rep$all_integer <- all(mat %% 1 == 0, na.rm = TRUE)
  rep$has_negative <- any(mat < 0, na.rm = TRUE)
  if (rep$has_negative) {
    warnings_ <- c(warnings_, paste(
      "The abundance table contains negative values, so it has already been",
      "transformed (CLR or similar). Counts-based methods (DESeq2, ANCOM-BC2,",
      "rarefaction, Chao1/ACE/Fisher) are invalid on this input."))
  }

  # --- 6. Read depth --------------------------------------------------------
  depths <- phyloseq::sample_sums(ps)
  rep$depth_min    <- min(depths)
  rep$depth_max    <- max(depths)
  rep$depth_median <- stats::median(depths)
  rep$depth_mean   <- mean(depths)
  rep$depth_fold_range <- if (rep$depth_min > 0) rep$depth_max / rep$depth_min else Inf
  rep$samples_below_warn <- names(depths)[depths < min_depth_warn]
  rep$zero_depth_samples <- names(depths)[depths == 0]

  if (length(rep$zero_depth_samples)) {
    problems <- c(problems, sprintf("%d sample(s) have zero total reads: %s",
      length(rep$zero_depth_samples),
      paste(utils::head(rep$zero_depth_samples, 5), collapse = ", ")))
  }
  if (!rep$is_relative && length(rep$samples_below_warn)) {
    warnings_ <- c(warnings_, sprintf(
      "%d sample(s) below %d reads. Decide explicitly whether to drop them (02-qc.R) — %s",
      length(rep$samples_below_warn), min_depth_warn,
      paste(utils::head(rep$samples_below_warn, 5), collapse = ", ")))
  }
  if (!rep$is_relative && is.finite(rep$depth_fold_range) && rep$depth_fold_range > 10) {
    warnings_ <- c(warnings_, sprintf(
      "Depth varies %.0f-fold across samples (%.0f to %.0f). Normalization is not optional here; see reference/03.",
      rep$depth_fold_range, rep$depth_min, rep$depth_max))
  }

  # --- 7. Sparsity ----------------------------------------------------------
  rep$zero_proportion <- mean(mat == 0, na.rm = TRUE)
  rep$taxa_all_zero   <- sum(phyloseq::taxa_sums(ps) == 0)
  if (rep$zero_proportion > 0.80) {
    warnings_ <- c(warnings_, sprintf(
      "%.1f%% of the table is zeros. Consider zero-inflated models (reference/07) and read reference/11 on structural zeros.",
      100 * rep$zero_proportion))
  }
  if (rep$taxa_all_zero > 0) {
    warnings_ <- c(warnings_, sprintf(
      "%d taxa have zero reads in every sample and contribute nothing but multiple-testing burden.",
      rep$taxa_all_zero))
  }

  # --- 8. Tree provenance ---------------------------------------------------
  if (rep$has_tree) {
    rep$tree_is_real <- mfg_tree_is_real(ps)
    rep$tree_tips    <- length(phyloseq::phy_tree(ps)$tip.label)
    rep$tree_rooted  <- ape::is.rooted(phyloseq::phy_tree(ps))
    if (is.na(rep$tree_is_real)) {
      warnings_ <- c(warnings_, paste(
        "A tree is present but its provenance is unrecorded. UniFrac and Faith's PD",
        "on a placeholder tree return numbers that look valid and encode nothing.",
        "Confirm with mfg_mark_tree_real(ps, TRUE/FALSE) before using them."))
    }
    if (!isTRUE(rep$tree_rooted)) {
      warnings_ <- c(warnings_, paste(
        "The tree is unrooted. Faith's PD with include.root = TRUE and unweighted",
        "UniFrac both depend on rooting; root it (ape::midpoint or phangorn) first."))
    }
  }

  # --- 9. Metadata variables ------------------------------------------------
  if (rep$has_meta) {
    md <- mfg_meta(ps)
    rep$variables <- names(md)
    rep$variable_types <- vapply(md, function(x) class(x)[1], character(1))
    rep$variable_levels <- lapply(md, function(x) {
      if (is.numeric(x)) NULL else sort(unique(as.character(x[!is.na(x)])))
    })
    rep$variable_n_missing <- vapply(md, function(x) sum(is.na(x)), integer(1))

    # A character column that looks categorical is not a factor yet, so its
    # reference level is alphabetical rather than chosen. Every model
    # coefficient is relative to that level, so it must be deliberate.
    chr_cats <- names(md)[vapply(md, function(x)
      is.character(x) && length(unique(x[!is.na(x)])) %in% 2:20, logical(1))]
    rep$unset_factors <- chr_cats
    if (length(chr_cats)) {
      warnings_ <- c(warnings_, sprintf(
        "These are character columns, not factors, so their reference level is alphabetical: %s. Set it deliberately with factor(..., levels = ) — every coefficient is relative to it.",
        paste(chr_cats, collapse = ", ")))
    }

    if (!is.null(group_var)) {
      if (!group_var %in% names(md)) {
        problems <- c(problems, sprintf("Group variable '%s' is not in the metadata. Available: %s",
                                        group_var, paste(names(md), collapse = ", ")))
      } else {
        g <- md[[group_var]]
        rep$group_var <- group_var
        rep$group_n   <- table(g, useNA = "ifany")   # display: missing stays visible
        rep$group_n_missing <- sum(is.na(g))
        # Balance and minimum-size checks run over REAL levels only. Counting NA
        # as a level reports a balanced design as lopsided (52/52/6 reads as
        # 8.7:1) and, worse, buries the thing that actually matters: samples with
        # no group assignment are dropped silently by almost every test.
        gn <- table(g[!is.na(g)])
        rep$group_n_observed <- gn
        if (rep$group_n_missing > 0) {
          warnings_ <- c(warnings_, sprintf(paste(
            "%d of %d samples have no value for '%s'. Most tests drop them without",
            "saying so, which changes n and the multiple-testing denominator.",
            "Exclude them explicitly and log the rule, or recode them."),
            rep$group_n_missing, length(g), group_var))
        }
        if (length(gn) && any(gn < 3)) {
          warnings_ <- c(warnings_, sprintf(
            "Group sizes are %s. Groups under ~3 support no meaningful test; see reference/12 on power.",
            paste(sprintf("%s=%d", names(gn), as.integer(gn)), collapse = ", ")))
        }
        if (length(gn) >= 2) {
          bal <- max(gn) / min(gn)
          rep$group_imbalance <- bal
          if (bal > 3) {
            warnings_ <- c(warnings_, sprintf(
              "Group sizes are unbalanced %.1f:1. PERMANOVA and ANOSIM are both sensitive to this; see reference/11.", bal))
          }
        }
      }
    }

    # --- 10. Repeated measures -------------------------------------------
    rm_det <- mfg_detect_repeated_measures(md, subject_candidates = subject_var)
    rep$repeated_measures <- rm_det
    if (isTRUE(rm_det$repeated)) {
      warnings_ <- c(warnings_, sprintf(
        "Repeated measures detected: '%s' has up to %d samples per level across %d levels. Independent-samples tests (Wilcoxon, Kruskal-Wallis, plain PERMANOVA) treat these as independent draws and overstate significance. Use a mixed model or a blocked/strata design; see reference/07.",
        rm_det$subject_var, rm_det$max_per_subject, rm_det$n_subjects))
    }
  }

  rep$problems <- problems
  rep$warnings <- warnings_
  rep$passed   <- length(problems) == 0

  mfg_log("intake", "validated", list(
    passed = rep$passed, n_problems = length(problems), n_warnings = length(warnings_),
    taxa = rep$n_taxa, samples = rep$n_samples,
    depth_min = rep$depth_min, depth_max = rep$depth_max,
    zero_proportion = round(rep$zero_proportion, 4),
    is_relative = rep$is_relative))

  if (strict && !rep$passed) {
    stop("Input validation failed:\n  - ", paste(problems, collapse = "\n  - "),
         "\n\nFix these before any analysis. Pass strict = FALSE only to inspect.",
         call. = FALSE)
  }
  class(rep) <- c("mfg_validation", "list")
  rep
}

#' @export
print.mfg_validation <- function(x, ...) {
  cat("=== MicroFitGut input validation ===\n")
  cat(sprintf("%d taxa x %d samples", x$n_taxa, x$n_samples))
  if (length(x$ranks)) cat(sprintf("   ranks: %s", paste(x$ranks, collapse = ", ")))
  cat("\n")
  cat(sprintf("Slots: otu=%s tax=%s meta=%s tree=%s seqs=%s\n",
              x$has_otu, x$has_tax, x$has_meta, x$has_tree, x$has_seqs))
  cat(sprintf("Values: %s, %s\n",
              if (isTRUE(x$is_relative)) "relative abundances" else "read counts",
              if (isTRUE(x$all_integer)) "all integer" else "non-integer present"))
  cat(sprintf("Depth: min %.0f | median %.0f | max %.0f  (%.1f-fold range)\n",
              x$depth_min, x$depth_median, x$depth_max, x$depth_fold_range))
  cat(sprintf("Sparsity: %.1f%% zeros | %d taxa all-zero\n",
              100 * x$zero_proportion, x$taxa_all_zero))
  if (!is.null(x$group_n)) {
    cat(sprintf("Group '%s': %s\n", x$group_var,
                paste(sprintf("%s=%d", names(x$group_n), as.integer(x$group_n)),
                      collapse = ", ")))
  }
  if (isTRUE(x$repeated_measures$repeated)) {
    cat(sprintf("Repeated measures: '%s' (%d levels, up to %d each)\n",
                x$repeated_measures$subject_var, x$repeated_measures$n_subjects,
                x$repeated_measures$max_per_subject))
  }
  if (length(x$problems)) {
    cat("\nPROBLEMS (block the analysis):\n")
    for (p in x$problems) cat("  - ", p, "\n", sep = "")
  }
  if (length(x$warnings)) {
    cat("\nWARNINGS (must be stated in the report):\n")
    for (w in x$warnings) cat("  - ", w, "\n", sep = "")
  }
  if (x$passed && !length(x$warnings)) cat("\nAll checks passed with no warnings.\n")
  invisible(x)
}

#' Rarefaction-style depth summary table for the report.
mfg_depth_table <- function(ps) {
  d <- phyloseq::sample_sums(ps)
  data.frame(sample = names(d), depth = as.numeric(d),
             observed_taxa = as.numeric(phyloseq::estimate_richness(
               ps, measures = "Observed")[["Observed"]]),
             stringsAsFactors = FALSE)[order(d), ]
}
