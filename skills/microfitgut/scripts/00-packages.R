# =============================================================================
# MicroFitGut — 00-packages.R
#
# Loads the package set and asserts versions. Fails loudly when something is
# missing rather than letting an analysis improvise a substitute, because a
# silent fallback (say, Kruskal-Wallis standing in for ANCOM-BC2) changes the
# statistical claim without changing the narrative around it.
#
# Usage:
#   source("scripts/00-packages.R")
#   mfg_require("beta")           # load one stage's packages
#   mfg_require(c("alpha","da"))  # or several
#   mfg_session_info()            # provenance block for the report
# =============================================================================

# ── The registry ─────────────────────────────────────────────────────────────
#
# Grouped by analysis stage so a run loads only what it needs. `min` is the
# lowest version the reference scripts were written against; `repo` tells the
# user where to get it. `note` is what breaks if it is absent.

MFG_PACKAGES <- list(

  core = list(
    list(pkg = "phyloseq",   min = "1.40.0", repo = "bioc",
         note = "the container every stage reads and writes"),
    list(pkg = "microbiome", min = "1.18.0", repo = "bioc",
         note = "transform(), core(), meta(), summarize_phyloseq()"),
    list(pkg = "vegan",      min = "2.6.0",  repo = "cran",
         note = "distances, adonis2, betadisper, anosim, mrpp, rarecurve"),
    list(pkg = "ape",        min = "5.6.0",  repo = "cran",
         note = "tree handling for UniFrac and plot_tree"),
    list(pkg = "dplyr",      min = "1.1.0",  repo = "cran",
         note = "result-table filtering and ordering"),
    list(pkg = "ggplot2",    min = "3.4.0",  repo = "cran",
         note = "every figure")
  ),

  intake = list(
    list(pkg = "Biostrings", min = "2.64.0", repo = "bioc",
         note = "reading representative sequences into refseq()")
  ),

  alpha = list(
    list(pkg = "picante",    min = "1.8.2",  repo = "cran",
         note = "pd() — Faith's phylogenetic diversity. No substitute exists"),
    list(pkg = "FSA",        min = "0.9.0",  repo = "cran",
         note = "dunnTest() post-hoc and Summarize() per-group tables"),
    list(pkg = "ggpubr",     min = "0.5.0",  repo = "cran",
         note = "stat_compare_means() significance annotation, ggboxplot"),
    list(pkg = "DescTools",  min = "0.99.40", repo = "cran",
         note = "NemenyiTest() — the alternative non-parametric post-hoc"),
    list(pkg = "multcomp",   min = "1.4.20", repo = "cran",
         note = "general linear hypotheses for parametric post-hoc")
  ),

  beta = list(
    list(pkg = "dendextend", min = "1.16.0", repo = "cran",
         note = "coloured dendrograms from hclust"),
    list(pkg = "cowplot",    min = "1.1.0",  repo = "cran",
         note = "plot_grid() panel assembly, panel_border()")
  ),

  da = list(
    list(pkg = "DESeq2",     min = "1.36.0", repo = "bioc",
         note = "negative-binomial DA, VST, plotMA, plotDispEsts"),
    list(pkg = "ALDEx2",     min = "1.28.0", repo = "bioc",
         note = "Dirichlet-multinomial Monte Carlo CLR differential abundance"),
    list(pkg = "ANCOMBC",    min = "2.0.0",  repo = "bioc",
         note = "ancombc2() — bias correction and structural zeros"),
    list(pkg = "MicrobiomeStat", min = "1.1", repo = "cran",
         note = "linda() — CLR-based linear model; accepts counts or proportions"),
    list(pkg = "Maaslin2",   min = "1.10.0", repo = "bioc",
         note = "Maaslin2() — per-feature general linear models with covariates")
  ),

  models = list(
    list(pkg = "lme4",       min = "1.1.30", repo = "cran",
         note = "lmer() and glmer() mixed models"),
    list(pkg = "lmerTest",   min = "3.1.0",  repo = "cran",
         note = "Satterthwaite p-values for lmer, which lme4 withholds"),
    list(pkg = "glmmTMB",    min = "1.1.5",  repo = "cran",
         note = "negative-binomial GLMM (nbinom2) for overdispersed counts"),
    list(pkg = "car",        min = "3.1.0",  repo = "cran",
         note = "Anova() Type-II/III tests for glmer output"),
    list(pkg = "pscl",       min = "1.5.5",  repo = "cran",
         note = "zeroinfl() and hurdle() — ZIP, ZINB, ZHP, ZHNB"),
    list(pkg = "lmtest",     min = "0.9.40", repo = "cran",
         note = "lrtest() for nested zero-model comparison")
  ),

  bayes = list(
    list(pkg = "brms",       min = "2.18.0", repo = "cran",
         note = "Bayesian joint models over multiple responses via mvbind()"),
    list(pkg = "loo",        min = "2.5.0",  repo = "cran",
         note = "LOO-CV model comparison — elpd_diff and se_diff")
  ),

  exploratory = list(
    list(pkg = "metacoder",  min = "0.3.5",  repo = "cran",
         note = "heat trees"),
    list(pkg = "igraph",     min = "1.3.0",  repo = "cran",
         note = "correlation networks and sample networks"),
    list(pkg = "ggtree",     min = "3.4.0",  repo = "bioc",
         note = "publication phylogenies"),
    list(pkg = "cluster",    min = "2.1.0",  repo = "cran",
         note = "k-means / PAM non-hierarchical clustering"),
    list(pkg = "phangorn",   min = "2.10.0", repo = "cran",
         note = "NJ trees, dist.ml, midpoint rooting"),
    list(pkg = "pheatmap",   min = "1.0.12", repo = "cran",
         note = "CLR heatmaps")
  ),

  plots = list(
    list(pkg = "patchwork",  min = "1.1.0",  repo = "cran",
         note = "figure composition"),
    list(pkg = "viridis",    min = "0.6.0",  repo = "cran",
         note = "continuous colour scales"),
    list(pkg = "ggrepel",    min = "0.9.0",  repo = "cran",
         note = "non-overlapping labels on ordinations and volcanoes"),
    list(pkg = "scales",     min = "1.2.0",  repo = "cran",
         note = "axis formatting")
  ),

  report = list(
    list(pkg = "knitr",      min = "1.40",   repo = "cran",
         note = "kable() tables and report rendering"),
    list(pkg = "rmarkdown",  min = "2.18",   repo = "cran",
         note = "render() the analysis report")
  )
)

# ── Install hints ────────────────────────────────────────────────────────────

mfg_install_hint <- function(pkg, repo) {
  switch(repo,
    bioc = sprintf('BiocManager::install("%s")', pkg),
    cran = sprintf('install.packages("%s")', pkg),
    github = sprintf('remotes::install_github("%s")', pkg),
    sprintf('install.packages("%s")', pkg)
  )
}

#' Load one or more stage groups, asserting presence and minimum version.
#'
#' Stops on the first unmet requirement with the exact install command and what
#' the analysis loses without it. A missing package is a blocked analysis, not a
#' reason to pick a different method.
mfg_require <- function(groups = "core", quiet = TRUE) {
  groups <- unique(c("core", groups))
  unknown <- setdiff(groups, names(MFG_PACKAGES))
  if (length(unknown)) {
    stop("Unknown package group(s): ", paste(unknown, collapse = ", "),
         "\nAvailable: ", paste(names(MFG_PACKAGES), collapse = ", "),
         call. = FALSE)
  }

  specs <- unlist(MFG_PACKAGES[groups], recursive = FALSE, use.names = FALSE)
  missing_specs <- list()
  stale_specs   <- list()

  for (s in specs) {
    if (!requireNamespace(s$pkg, quietly = TRUE)) {
      missing_specs[[length(missing_specs) + 1]] <- s
      next
    }
    if (!is.null(s$min) && utils::packageVersion(s$pkg) < s$min) {
      stale_specs[[length(stale_specs) + 1]] <- s
    }
  }

  if (length(missing_specs) || length(stale_specs)) {
    lines <- c("", "MicroFitGut cannot run this stage: package requirements are unmet.", "")
    if (length(missing_specs)) {
      lines <- c(lines, "MISSING:")
      for (s in missing_specs) {
        lines <- c(lines,
          sprintf("  %-12s needed for: %s", s$pkg, s$note),
          sprintf("  %-12s install with: %s", "", mfg_install_hint(s$pkg, s$repo)))
      }
      lines <- c(lines, "")
    }
    if (length(stale_specs)) {
      lines <- c(lines, "TOO OLD:")
      for (s in stale_specs) {
        lines <- c(lines, sprintf("  %-12s have %s, need >= %s  (%s)",
                                  s$pkg, utils::packageVersion(s$pkg), s$min, s$note))
      }
      lines <- c(lines, "")
    }
    lines <- c(lines,
      "Install the packages above and re-run. Do not substitute a different",
      "method for a missing one — that changes the statistical claim.", "")
    stop(paste(lines, collapse = "\n"), call. = FALSE)
  }

  for (s in specs) {
    suppressPackageStartupMessages(
      library(s$pkg, character.only = TRUE, quietly = quiet, warn.conflicts = FALSE)
    )
  }
  invisible(vapply(specs, function(s) s$pkg, character(1)))
}

#' Report which groups are satisfied without stopping.
#'
#' For the pre-flight check at the top of a run, so the agent can tell the user
#' up front that (say) Bayesian joint models are unavailable, instead of
#' discovering it nine stages in.
mfg_check_packages <- function(groups = names(MFG_PACKAGES)) {
  rows <- list()
  for (g in groups) {
    for (s in MFG_PACKAGES[[g]]) {
      have <- requireNamespace(s$pkg, quietly = TRUE)
      ver  <- if (have) as.character(utils::packageVersion(s$pkg)) else NA_character_
      ok   <- have && (is.null(s$min) || utils::packageVersion(s$pkg) >= s$min)
      rows[[length(rows) + 1]] <- data.frame(
        group = g, package = s$pkg, installed = have, version = ver,
        min_required = s$min %||% NA_character_, ok = ok,
        install_with = if (ok) NA_character_ else mfg_install_hint(s$pkg, s$repo),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out[!duplicated(out$package), ]
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Versions of everything actually loaded, for the report's provenance block.
#'
#' 10-reporting-standards.md requires software versions in every summary. This
#' is where that comes from — read off the session, never typed by hand.
mfg_session_info <- function() {
  loaded <- rev(names(utils::sessionInfo()$otherPkgs))
  pkgs <- vapply(loaded, function(p) as.character(utils::packageVersion(p)), character(1))
  list(
    r_version   = paste(R.version$major, R.version$minor, sep = "."),
    platform    = R.version$platform,
    run_at      = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    packages    = data.frame(package = loaded, version = unname(pkgs),
                             stringsAsFactors = FALSE)
  )
}

#' Print the pre-flight table in a form that reads well in a terminal.
mfg_print_package_status <- function(groups = names(MFG_PACKAGES)) {
  st <- mfg_check_packages(groups)
  bad <- st[!st$ok, , drop = FALSE]
  cat(sprintf("MicroFitGut package check: %d of %d requirements satisfied\n",
              sum(st$ok), nrow(st)))
  if (nrow(bad)) {
    cat("\nUnavailable — these analyses are blocked until installed:\n")
    for (i in seq_len(nrow(bad))) {
      cat(sprintf("  [%s] %s  -> %s\n", bad$group[i], bad$package[i], bad$install_with[i]))
    }
  } else {
    cat("All stages available.\n")
  }
  invisible(st)
}
