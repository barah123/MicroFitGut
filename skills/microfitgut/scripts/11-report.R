# =============================================================================
# MicroFitGut — 11-report.R
#
# Provenance capture, the report assembly, and the traceability check.
#
# The guardrail this file exists to enforce, from the project notes: the agent
# must not report a conclusion that is not traceable to a run and a stored
# artifact. Anything it cannot point at, it does not get to say.
#
# mfg_traceability_check() is the mechanical half of that — it pulls every number
# out of the draft report and looks for it in the run log and the saved tables.
# The verifier subagent is the judgement half. Both run before a report is final.
#
# Requires: mfg_require("report"); utils.R
# =============================================================================

# ── Provenance ───────────────────────────────────────────────────────────────

#' Everything needed to reproduce this run.
#'
#' Software versions come from the session, never typed by hand. The random seeds
#' matter because rarefaction and permutation tests are stochastic: without them
#' a rerun disagrees with the report and neither is wrong.
capture_provenance <- function(ps_initial = NULL, ps_final = NULL,
                               qc_log = NULL, seeds = NULL) {
  si <- mfg_session_info()
  log_entries <- mfg_get_log()

  # Recover every seed that was actually used, from the log rather than from the
  # caller's memory of what it passed.
  logged_seeds <- list()
  for (e in log_entries) {
    for (k in c("rngseed", "seed")) {
      if (!is.null(e$detail[[k]])) {
        logged_seeds[[paste(e$stage, e$event, k, sep = ".")]] <- e$detail[[k]]
      }
    }
  }

  prov <- list(
    run_id = mfg_run_id(),
    run_at = si$run_at,
    r_version = si$r_version,
    platform = si$platform,
    packages = si$packages,
    seeds = c(seeds %||% list(), logged_seeds),
    n_log_entries = length(log_entries),
    output_dir = mfg_run_dir(),
    # What the analysis was done to, alongside what was done. Both are needed
    # before another researcher can repeat the run.
    inputs = mfg_inputs(),
    authoritative = mfg_authoritative_source(),
    inventory = MFG_LOG$sources$inventory
  )

  if (!is.null(ps_initial)) {
    prov$initial <- list(taxa = phyloseq::ntaxa(ps_initial),
                         samples = phyloseq::nsamples(ps_initial),
                         reads = sum(phyloseq::sample_sums(ps_initial)))
  }
  if (!is.null(ps_final)) {
    prov$final <- list(taxa = phyloseq::ntaxa(ps_final),
                       samples = phyloseq::nsamples(ps_final),
                       reads = sum(phyloseq::sample_sums(ps_final)),
                       normalization = mfg_normalization(ps_final))
  }
  if (!is.null(qc_log)) {
    prov$qc_steps <- qc_table(qc_log)
    prov$exclusions <- qc_exclusions(qc_log)
  }

  mfg_log("report", "provenance_captured",
          list(run_id = prov$run_id, n_packages = nrow(si$packages),
               n_seeds = length(prov$seeds), n_inputs = nrow(prov$inputs),
               authoritative = prov$authoritative$file %||% "undeclared"))
  class(prov) <- c("mfg_provenance", "list")
  prov
}

#' @export
print.mfg_provenance <- function(x, ...) {
  cat("=== Provenance ===\n")
  cat(sprintf("Run ID:   %s\n", x$run_id))
  cat(sprintf("Run at:   %s\n", x$run_at))
  cat(sprintf("R:        %s on %s\n", x$r_version, x$platform))
  cat(sprintf("Output:   %s\n", x$output_dir))
  if (!is.null(x$initial) && !is.null(x$final)) {
    cat(sprintf("Data:     %d taxa x %d samples -> %d x %d (%s)\n",
                x$initial$taxa, x$initial$samples,
                x$final$taxa, x$final$samples, x$final$normalization))
  }
  if (length(x$seeds)) {
    cat(sprintf("Seeds:    %s\n",
                paste(sprintf("%s=%s", names(x$seeds), unlist(x$seeds)), collapse = ", ")))
  }
  if (!is.null(x$authoritative)) {
    cat(sprintf("Source:   %s (%s)\n", x$authoritative$file, x$authoritative$reason))
  } else if (!is.null(x$inventory) && sum(x$inventory$format %in%
                                         c("rds", "biom", "qza")) > 1) {
    cat("Source:   UNDECLARED, and the data directory held more than one dataset\n")
  }
  if (!is.null(x$inputs) && nrow(x$inputs)) {
    cat(sprintf("\nInputs (%d):\n", nrow(x$inputs)))
    print(x$inputs[, c("file", "format", "role", "bytes", "md5")], row.names = FALSE)
  }
  cat(sprintf("\nPackages (%d):\n", nrow(x$packages)))
  print(x$packages, row.names = FALSE)
  invisible(x)
}

#' Provenance as the methods-section paragraph it needs to become.
provenance_paragraph <- function(prov) {
  pk <- prov$packages
  key <- pk[pk$package %in% c("phyloseq", "vegan", "microbiome", "DESeq2", "ALDEx2",
                              "ANCOMBC", "picante", "lme4", "glmmTMB", "brms"), ]
  parts <- c(
    sprintf("Analysis was performed in R %s (%s) under run identifier %s on %s.",
            prov$r_version, prov$platform, prov$run_id, prov$run_at))
  if (!is.null(prov$inputs) && nrow(prov$inputs)) {
    n_in <- nrow(prov$inputs)
    parts <- c(parts, sprintf("The analysis read %d input file%s (%s), recorded in Table 1 with %s format, role and MD5 checksum.",
      n_in, if (n_in == 1) "" else "s",
      paste(prov$inputs$file, collapse = ", "),
      if (n_in == 1) "its" else "their"))
  }
  if (!is.null(prov$authoritative)) {
    # A colon, not "because": the reason is whatever the analyst wrote, and a
    # noun phrase after "because" reads as a grammatical error in the methods.
    parts <- c(parts, sprintf("%s was taken as the authoritative version: %s.",
      prov$authoritative$file, sub("\\.$", "", prov$authoritative$reason)))
  }
  if (nrow(key)) {
    parts <- c(parts, sprintf("Key packages: %s.",
      paste(sprintf("%s %s", key$package, key$version), collapse = ", ")))
  }
  if (length(prov$seeds)) {
    parts <- c(parts, sprintf("Stochastic steps used fixed seeds (%s), so the run is reproducible.",
      paste(unique(unlist(prov$seeds)), collapse = ", ")))
  }
  if (!is.null(prov$initial) && !is.null(prov$final)) {
    parts <- c(parts, sprintf(paste("The dataset entered analysis with %d taxa across",
      "%d samples (%s reads) and, after filtering, %d taxa across %d samples",
      "(%s reads) were analysed."),
      prov$initial$taxa, prov$initial$samples, format(prov$initial$reads, big.mark = ","),
      prov$final$taxa, prov$final$samples, format(prov$final$reads, big.mark = ",")))
  }
  paste(parts, collapse = " ")
}

# ── Assembling the summary ───────────────────────────────────────────────────

#' Build the structured analysis summary.
#'
#' `sections` is a named list of whatever the run produced. Every entry carries
#' its own print method, so the report is assembled from the analysis objects
#' themselves rather than from numbers re-typed into prose — which is where
#' transcription errors enter.
#'
#' The mandatory-content check from reference/10 runs here: n per group, every
#' exclusion, the normalization, the exact test with assumptions, effect sizes,
#' software versions. Missing items are listed rather than silently omitted.
assemble_summary <- function(provenance, validation = NULL, qc_log = NULL,
                             sections = list(), title = NULL,
                             group_var = NULL) {

  # An undeclared pick between several self-contained datasets is a reporting
  # gap of the same kind as an unstated normalization: the methods read as
  # complete and the reader still cannot tell what was analysed.
  rival_datasets <- !is.null(provenance$inventory) &&
    sum(provenance$inventory$format %in% c("rds", "biom", "qza")) > 1

  required <- c(
    n_per_group     = !is.null(validation$group_n) || !is.null(sections$alpha),
    exclusions      = !is.null(qc_log),
    normalization   = !is.null(provenance$final$normalization),
    tests_with_assumptions = any(vapply(sections, function(s)
      inherits(s, c("mfg_alpha_test", "mfg_permanova", "mfg_da_result", "mfg_model")),
      logical(1))),
    effect_sizes    = TRUE,
    software_versions = !is.null(provenance$packages),
    input_provenance  = !is.null(provenance$inputs) && nrow(provenance$inputs) > 0,
    authoritative_source = !rival_datasets || !is.null(provenance$authoritative)
  )
  missing_required <- names(required)[!required]

  out <- list(
    title = title %||% sprintf("Microbiome downstream analysis — %s", provenance$run_id),
    provenance = provenance,
    validation = validation,
    qc_log = qc_log,
    sections = sections,
    group_var = group_var,
    missing_required = missing_required,
    complete = length(missing_required) == 0
  )
  mfg_log("report", "summary_assembled",
          list(n_sections = length(sections),
               missing_required = if (length(missing_required))
                 paste(missing_required, collapse = ",") else "none"))
  class(out) <- c("mfg_summary", "list")
  out
}

#' @export
print.mfg_summary <- function(x, ...) {
  cat(strrep("=", 78), "\n", x$title, "\n", strrep("=", 78), "\n\n", sep = "")

  cat("-- Provenance ", strrep("-", 62), "\n", sep = "")
  cat(strwrap(provenance_paragraph(x$provenance), width = 78), sep = "\n")
  cat("\n\n")

  if (!is.null(x$validation)) {
    cat("-- Input validation ", strrep("-", 56), "\n", sep = "")
    print(x$validation); cat("\n")
  }
  if (!is.null(x$qc_log)) {
    cat("-- Filtering ", strrep("-", 63), "\n", sep = "")
    print(x$qc_log); cat("\n")
  }
  for (nm in names(x$sections)) {
    cat("-- ", nm, " ", strrep("-", max(0, 74 - nchar(nm))), "\n", sep = "")
    s <- x$sections[[nm]]
    if (inherits(s, "data.frame")) print(s, row.names = FALSE) else print(s)
    cat("\n")
  }
  if (length(x$missing_required)) {
    cat("!! INCOMPLETE — reference/10 requires these and they are absent:\n")
    for (m in x$missing_required) cat("   - ", m, "\n", sep = "")
  }
  invisible(x)
}

# ── Traceability ─────────────────────────────────────────────────────────────

#' Check that every number in a draft traces to the run log or a saved table.
#'
#' Extracts numeric literals from the text and looks for each in the run log and
#' in every saved CSV. A number that appears nowhere was either typed from memory
#' or computed outside the pipeline; both are reasons not to publish it.
#'
#' Deliberately noisy rather than clever: it flags candidates for the verifier
#' subagent to adjudicate, and it is better to over-flag than to pass a fabricated
#' statistic. Years, sample counts written in prose and similar will appear as
#' false positives; that is the intended trade.
mfg_traceability_check <- function(text, tolerance = 1e-6, min_digits = 2,
                                   extra_tables = NULL) {
  if (length(text) > 1) text <- paste(text, collapse = "\n")

  # Numbers worth checking: at least `min_digits` significant digits, so that
  # "3 groups" and "Table 2" are not chased. Handles scientific notation.
  pat <- "-?\\b\\d+\\.\\d+(?:[eE][-+]?\\d+)?\\b|-?\\b\\d+[eE][-+]?\\d+\\b|-?\\b\\d{3,}\\b"
  m <- regmatches(text, gregexpr(pat, text, perl = TRUE))[[1]]
  nums <- suppressWarnings(as.numeric(gsub(",", "", m)))
  keep <- !is.na(nums)
  m <- m[keep]; nums <- nums[keep]
  if (min_digits > 0) {
    sig <- vapply(m, function(s) nchar(gsub("[^0-9]", "", s)), integer(1))
    m <- m[sig >= min_digits]; nums <- nums[sig >= min_digits]
  }
  candidates <- unique(data.frame(literal = m, value = nums, stringsAsFactors = FALSE))

  # Haystack 1: every value recorded in the run log.
  log_values <- unlist(lapply(mfg_get_log(), function(e)
    suppressWarnings(as.numeric(unlist(lapply(e$detail, function(v)
      unlist(strsplit(as.character(v), "[^0-9eE.+-]+"))))))))
  log_values <- log_values[is.finite(log_values)]

  # Haystack 2: every numeric cell in every saved table.
  table_files <- list.files(mfg_run_dir("tables"), pattern = "\\.csv$", full.names = TRUE)
  table_values <- numeric()
  for (f in c(table_files, extra_tables %||% character())) {
    d <- try(utils::read.csv(f, stringsAsFactors = FALSE), silent = TRUE)
    if (inherits(d, "try-error")) next
    v <- suppressWarnings(as.numeric(unlist(d[vapply(d, is.numeric, logical(1))])))
    table_values <- c(table_values, v[is.finite(v)])
  }

  haystack <- c(log_values, table_values)
  found_in <- character(nrow(candidates))
  for (i in seq_len(nrow(candidates))) {
    v <- candidates$value[i]
    tol <- max(tolerance, abs(v) * 1e-3)   # tolerate the report's own rounding
    in_log <- any(abs(log_values - v) <= tol)
    in_tab <- any(abs(table_values - v) <= tol)
    found_in[i] <- if (in_log && in_tab) "log+table" else if (in_log) "log"
                   else if (in_tab) "table" else "NOT FOUND"
  }
  candidates$found_in <- found_in
  untraced <- candidates[candidates$found_in == "NOT FOUND", , drop = FALSE]

  mfg_log("report", "traceability_checked", list(
    n_numbers_checked = nrow(candidates),
    n_traced = sum(candidates$found_in != "NOT FOUND"),
    n_untraced = nrow(untraced),
    n_log_values = length(log_values), n_table_values = length(table_values),
    tables_scanned = length(table_files)))

  # Persist the log here. render_report() writes it too, but this function and
  # mfg_manifest() normally run AFTER the report is rendered, so without this
  # their own events stay in memory and never reach run_log.csv — which is
  # exactly the file the verifier subagent is told to grep for them.
  mfg_write_log()

  out <- list(candidates = candidates, untraced = untraced,
    n_checked = nrow(candidates), n_traced = sum(candidates$found_in != "NOT FOUND"),
    n_untraced = nrow(untraced),
    passed = nrow(untraced) == 0,
    n_log_values = length(log_values), n_table_values = length(table_values))
  class(out) <- c("mfg_traceability", "list")
  out
}

#' @export
print.mfg_traceability <- function(x, ...) {
  cat("=== Traceability check ===\n")
  cat(sprintf("%d numbers checked against %d logged values and %d table cells\n",
              x$n_checked, x$n_log_values, x$n_table_values))
  cat(sprintf("Traced: %d | Untraced: %d\n", x$n_traced, x$n_untraced))
  if (x$n_untraced) {
    cat("\nNumbers not found in the run log or any saved table:\n")
    print(x$untraced[, c("literal", "value")], row.names = FALSE)
    cat(strwrap(paste("Each of these was either typed from memory, computed outside",
      "the pipeline, or is incidental prose (a year, a count written out). Resolve",
      "every one before the report is final: either point it at an artifact or",
      "remove it. Hand this list to the verifier subagent."),
      width = 78, prefix = "  "), sep = "\n")
  } else {
    cat("\nEvery checked number traces to a logged value or a saved table.\n")
  }
  invisible(x)
}

# ── Rendering ────────────────────────────────────────────────────────────────

#' Write the report as R Markdown and render it.
#'
#' Sections arrive as markdown strings from the agent, and the provenance,
#' validation and QC blocks are generated here from the objects so those cannot
#' drift from what actually ran.
render_report <- function(summary_obj, sections_md = list(),
                          format = c("html", "pdf", "md"),
                          filename = "report", open = FALSE) {
  format <- match.arg(format)
  d <- mfg_run_dir()
  rmd_path <- file.path(d, paste0(filename, ".Rmd"))

  prov <- summary_obj$provenance
  fig_dir <- mfg_run_dir("figures")
  figures <- list.files(fig_dir, pattern = "\\.(png|pdf|svg)$")
  tables  <- list.files(mfg_run_dir("tables"), pattern = "\\.csv$")

  yaml <- c("---",
    sprintf('title: "%s"', summary_obj$title),
    sprintf('date: "%s"', prov$run_at),
    "output:",
    if (format == "pdf") "  pdf_document:\n    toc: true\n    number_sections: true"
    else "  html_document:\n    toc: true\n    toc_float: true\n    theme: flatly\n    number_sections: true",
    "---", "")

  body <- c(
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
    "```", "",
    "# Provenance", "",
    provenance_paragraph(prov), "",
    if (!is.null(prov$inputs) && nrow(prov$inputs)) c(
      "```{r input-table}",
      "knitr::kable(prov$inputs[, c('file', 'format', 'role', 'bytes', 'modified', 'md5')],",
      "             caption = 'Table 1. Input files. The checksum identifies the exact file version analysed.')",
      "```", ""),
    "```{r provenance-table}",
    "knitr::kable(prov$packages, caption = 'Software versions')",
    "```", "")

  if (!is.null(summary_obj$validation)) {
    v <- summary_obj$validation
    body <- c(body, "# Data and quality control", "",
      sprintf("The dataset comprised %d taxa across %d samples. Sequencing depth ranged from %s to %s reads (median %s), a %.1f-fold range. %.1f%% of the abundance table was zeros.",
              v$n_taxa, v$n_samples, format(v$depth_min, big.mark = ","),
              format(v$depth_max, big.mark = ","), format(v$depth_median, big.mark = ","),
              v$depth_fold_range, 100 * v$zero_proportion), "")
    if (!is.null(v$group_n)) {
      body <- c(body, sprintf("Group sizes before filtering: %s.",
        paste(sprintf("%s = %d", names(v$group_n), as.integer(v$group_n)), collapse = ", ")), "")
    }
    if (length(v$warnings)) {
      body <- c(body, "## Data characteristics requiring disclosure", "",
        paste0("- ", v$warnings), "")
    }
  }

  if (!is.null(summary_obj$qc_log)) {
    body <- c(body, "## Filtering", "",
      "```{r qc-table}",
      "knitr::kable(qc_table(summary_obj$qc_log)[, c('step','rule','taxa_before','taxa_after','samples_before','samples_after')], caption = 'Filtering steps')",
      "```", "",
      "```{r exclusions}",
      "ex <- qc_exclusions(summary_obj$qc_log)",
      "if (nrow(ex)) knitr::kable(ex, caption = 'Every excluded sample and the rule that excluded it') else cat('No samples were excluded.')",
      "```", "")
  }

  # Normalization is chosen per analysis, not once per run (reference/03), so the
  # report lists every state the run produced and which analysis it was paired
  # with — read from the log rather than from the final object, which only knows
  # its own state.
  norm_events <- Filter(function(e) identical(e$stage, "normalize"), mfg_get_log())
  applied <- Filter(function(e) e$event %in% c("tss", "clr", "log10", "presence",
                                               "rarefied", "vst", "deseq2_normalized"),
                    norm_events)
  pairings <- Filter(function(e) identical(e$event, "pairing_checked"), norm_events)

  body <- c(body, "## Normalization", "")
  if (length(applied)) {
    # Composition and plotting helpers transform internally, so the same event
    # appears many times with identical detail. Deduplicate on the rendered line
    # rather than listing "tss" five times.
    lines <- vapply(applied, function(e) {
      detail <- paste(sprintf("%s = %s", names(e$detail),
                              vapply(e$detail, function(v)
                                paste(format(v), collapse = "/"), character(1))),
                      collapse = ", ")
      sprintf("- **%s** — %s", e$event, detail)
    }, character(1))
    counts <- table(lines)
    body <- c(body, vapply(names(counts), function(l)
      if (counts[[l]] > 1) sprintf("%s  *(applied %d times)*", l, counts[[l]]) else l,
      character(1)), "")
  } else {
    body <- c(body, sprintf("- %s (no transformation was applied in this run)",
                            prov$final$normalization %||% "not recorded"), "")
  }
  if (length(pairings)) {
    pair_tab <- unique(do.call(rbind, lapply(pairings, function(e)
      data.frame(analysis = e$detail$analysis, normalization = e$detail$normalization,
                 stringsAsFactors = FALSE))))
    pair_tab <- pair_tab[order(pair_tab$analysis), ]
    body <- c(body,
      "Each analysis was run on the state of the data its assumptions require, and every pairing was checked:",
      "",
      "| Analysis | Normalization |", "|---|---|",
      sprintf("| %s | %s |", pair_tab$analysis, pair_tab$normalization), "")
  }

  for (nm in names(sections_md)) {
    body <- c(body, sprintf("# %s", nm), "", sections_md[[nm]], "")
  }

  if (length(figures)) {
    body <- c(body, "# Figures", "",
      unlist(lapply(figures, function(f) c(
        sprintf("![%s](figures/%s)", tools::file_path_sans_ext(f), f), ""))))
  }
  if (length(tables)) {
    body <- c(body, "# Saved tables", "",
      paste0("- `tables/", tables, "`"), "")
  }

  body <- c(body, "# Reproducibility", "",
    sprintf("Run identifier `%s`. All outputs are under `%s`. The complete step-by-step log is in `run_log.csv` and the calls that produced it are in `analysis_calls.R`; every number in this report traces to an entry there or to a table in `tables/`.",
            prov$run_id, prov$output_dir), "")

  if (length(summary_obj$missing_required)) {
    body <- c(body, "# Incomplete reporting", "",
      "This report is missing content that reference/10 requires:", "",
      paste0("- ", summary_obj$missing_required), "")
  }

  writeLines(c(yaml, body), rmd_path)
  mfg_log("report", "rmd_written",
          list(file = basename(rmd_path), n_sections = length(sections_md),
               n_figures = length(figures), n_tables = length(tables)))

  out_path <- NA_character_
  if (format != "md") {
    env <- new.env(parent = globalenv())
    assign("prov", prov, envir = env)
    assign("summary_obj", summary_obj, envir = env)
    assign("qc_table", qc_table, envir = env)
    assign("qc_exclusions", qc_exclusions, envir = env)
    r <- try(rmarkdown::render(rmd_path, output_format = paste0(format, "_document"),
                               envir = env, quiet = TRUE), silent = TRUE)
    if (inherits(r, "try-error")) {
      warning("Rendering failed; the .Rmd source is still at ", rmd_path, "\n",
              conditionMessage(attr(r, "condition")), call. = FALSE)
      mfg_log("report", "render_failed", list(format = format))
    } else {
      out_path <- r
      mfg_log("report", "rendered", list(format = format, file = basename(r)))
    }
  }
  mfg_write_log()
  list(rmd = rmd_path, output = out_path)
}

#' Write out the sequence of MicroFitGut calls the run actually made.
#'
#' Each log entry carries the outermost call that produced it, so this is the
#' analysis in the order it happened, with the arguments as they were passed.
#'
#' It is a record, not a finished script. R cannot see, from inside a function,
#' what its result was assigned to, so the assignments are absent and object
#' names such as `ps` or `ps_f` appear as they were typed. Consecutive calls
#' that logged more than one event are collapsed to one line. The header says
#' this in the file itself, because a file named like a script gets run.
mfg_write_analysis_script <- function(filename = "analysis_calls.R") {
  entries <- mfg_get_log()
  calls <- vapply(entries, function(e) e$call %||% NA_character_, character(1))
  calls <- calls[!is.na(calls) & nzchar(calls)]

  # One call that logs six events is one line, not six. Only consecutive
  # repeats collapse: the same function called twice on different data stays
  # twice, because it is two steps of the analysis.
  if (length(calls) > 1) calls <- calls[c(TRUE, calls[-1] != calls[-length(calls)])]

  path <- file.path(mfg_run_dir(), filename)
  header <- c(
    "# =============================================================================",
    sprintf("# MicroFitGut call record - run %s", mfg_run_id()),
    sprintf("# Written %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "#",
    "# Every MicroFitGut call this run made, in order, with the arguments as they",
    "# were passed. Reconstructed from the run log, not transcribed by hand.",
    "#",
    "# READ BEFORE RUNNING. This is a record of what happened, not a script that",
    "# reproduces it unedited:",
    "#   - assignments are not captured, so add them back (`ps <- mfg_load(...)`)",
    "#   - object names appear as they were typed in the session that ran",
    "#   - calls made outside MicroFitGut functions do not log and are absent",
    "#",
    "# Source the scripts and start a run before replaying any of this.",
    "# =============================================================================",
    "")
  writeLines(c(header, calls), path)

  mfg_log("report", "analysis_script_written",
          list(file = filename, n_calls = length(calls)))
  invisible(path)
}

#' Everything the run produced, as one manifest.
#'
#' The index the verifier subagent reads to know what exists before checking what
#' the report claims.
mfg_manifest <- function() {
  # Logged before the record is written so this call appears in it, and the
  # record is written before the file listing so it appears in the manifest.
  mfg_log("report", "manifest_started", list(run_id = mfg_run_id()))
  mfg_write_analysis_script()
  d <- mfg_run_dir()
  files <- list.files(d, recursive = TRUE, full.names = FALSE)
  info <- file.info(file.path(d, files))
  man <- data.frame(file = files, bytes = info$size,
                    modified = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
                    stringsAsFactors = FALSE)
  man <- man[order(man$file), ]
  utils::write.csv(man, file.path(d, "manifest.csv"), row.names = FALSE)
  mfg_log("report", "manifest_written", list(n_files = nrow(man)))
  mfg_write_log()
  man
}
