# =============================================================================
# MicroFitGut — 10-plots.R
#
# One palette, one theme, applied everywhere, so every figure in a run reads as
# part of the same publication set. Carried from CanisLupus 2.0.
#
# The colour rule, which matters more in figures than in a logo: in an ordination
# or a stacked bar, colour IS the encoding, so a palette that collapses under
# colour-vision deficiency destroys the data. Okabe-Ito is the standard
# colourblind-safe qualitative set and is the default here. viridis for
# continuous. Never a raw rainbow — its lightness is non-monotonic, so it invents
# banding that is not in the data.
#
# Every plot function returns a ggplot object rather than drawing, so the caller
# composes and saves through mfg_save_plot(), which logs the file.
#
# Requires: mfg_require(c("plots")); utils.R
# =============================================================================

# ── Palette ──────────────────────────────────────────────────────────────────

# Okabe-Ito extended with further distinguishable hues for taxa-heavy plots.
# The first eight are the Okabe-Ito set proper; beyond that, distinguishability
# degrades no matter what, which is why collapse_to_top() exists.
canis_categorical <- c(
  "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9",
  "#D55E00", "#F0E442", "#5D3A9B", "#117733", "#882255",
  "#44AA99", "#DDCC77", "#AA4499", "#88CCEE", "#999933",
  "#661100", "#6699CC", "#332288", "#AA7744", "#BBBBBB"
)

canis_ui <- list(
  primary   = "#0B6E4F",
  secondary = "#14746F",
  accent    = "#E9C46A",
  ink       = "#1B2A33",
  muted     = "#5A6C77",
  surface   = "#FFFFFF",
  canvas    = "#F5F7F9",
  border    = "#DFE5EA",
  danger    = "#B00020"
)

# MicroFitGut mark colours: teal spiral, amber dots. The amber rather than coral
# variant, because it opens a real luminance gap (ratio ~1.6:1) and sits on the
# blue-yellow axis, which survives every common colour-vision deficiency and
# separates in greyscale for black-and-white print.
mfg_brand <- list(teal = "#1D9E75", amber = "#E8A33D")

#' Categorical colours, recycled smoothly when a plot needs more than the set.
canis_colors <- function(n) {
  if (is.null(n) || is.na(n) || n < 1) return(canis_categorical[1])
  if (n <= length(canis_categorical)) return(canis_categorical[seq_len(n)])
  grDevices::colorRampPalette(canis_categorical)(n)
}

scale_fill_canis <- function(...) {
  ggplot2::discrete_scale("fill", palette = function(n) canis_colors(n), ...)
}
scale_color_canis <- function(...) {
  ggplot2::discrete_scale("colour", palette = function(n) canis_colors(n), ...)
}
scale_colour_canis <- scale_color_canis

#' Continuous scales: viridis, which is perceptually uniform and CVD-safe.
scale_fill_canis_c  <- function(...) viridis::scale_fill_viridis(option = "D", ...)
scale_color_canis_c <- function(...) viridis::scale_color_viridis(option = "D", ...)

#' Diverging scale for anything centred on zero (log fold change, CLR, correlation).
#'
#' A diverging scale must be symmetric about its midpoint or it misreads the
#' sign of the data, so the limits are forced symmetric.
scale_fill_canis_diverging <- function(limits = NULL, ...) {
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                               midpoint = 0, limits = limits, ...)
}

#' Shared plot theme: light, uncluttered, readable at publication size.
theme_canis <- function(base_size = 13) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      text             = ggplot2::element_text(colour = canis_ui$ink),
      plot.title       = ggplot2::element_text(face = "bold", size = base_size * 1.1,
                                               colour = canis_ui$ink,
                                               margin = ggplot2::margin(b = 8)),
      plot.subtitle    = ggplot2::element_text(colour = canis_ui$muted, size = base_size * 0.9),
      plot.caption     = ggplot2::element_text(colour = canis_ui$muted, size = base_size * 0.72,
                                               hjust = 0),
      axis.title       = ggplot2::element_text(colour = canis_ui$muted, size = base_size * 0.9),
      axis.text        = ggplot2::element_text(colour = canis_ui$muted, size = base_size * 0.82),
      panel.grid.major = ggplot2::element_line(colour = canis_ui$border, linewidth = 0.35),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = canis_ui$surface, colour = NA),
      plot.background  = ggplot2::element_rect(fill = canis_ui$surface, colour = NA),
      strip.text       = ggplot2::element_text(face = "bold", colour = canis_ui$ink,
                                               size = base_size * 0.88),
      legend.position  = "right",
      legend.title     = ggplot2::element_text(colour = canis_ui$muted, size = base_size * 0.85),
      legend.text      = ggplot2::element_text(size = base_size * 0.82),
      plot.margin      = ggplot2::margin(12, 12, 12, 12)
    )
}

#' Hide per-sample axis labels once there are too many to read.
sample_axis_theme <- function(n_samples, max_labels = 25) {
  if (n_samples > max_labels) {
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())
  } else {
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8))
  }
}

# Kept for backwards compatibility with any external code referencing it.
wolf_pal <- canis_categorical

#' Wrap caption text to the width a figure of this size can actually show.
#'
#' A caption is set at roughly base_size * 0.72 pt, which fits about 13 characters
#' per inch. Wrapping at a fixed 110 overflows a single- or double-column figure
#' and the sentence gets clipped at the page edge — which is worse than no caption,
#' because the reader cannot tell it was truncated.
mfg_wrap_caption <- function(x, width_in = 7.2, chars_per_inch = 13) {
  if (is.null(x) || !length(x) || all(is.na(x))) return(NULL)
  paste(strwrap(paste(x, collapse = " "),
                width = max(40, floor(width_in * chars_per_inch))), collapse = "\n")
}

# ── Figure dimensions ────────────────────────────────────────────────────────
#
# Journal figure widths. Passing these rather than eyeballing means text in a
# figure lands at the size it will actually print at, which is the difference
# between a readable axis and a resubmission.

MFG_FIGURE_SIZES <- list(
  single_column = list(width = 3.5, height = 3.0, dpi = 300,
                       note = "~89 mm, most journals' single column"),
  onehalf_column = list(width = 5.5, height = 4.0, dpi = 300,
                       note = "~140 mm"),
  double_column = list(width = 7.2, height = 5.0, dpi = 300,
                       note = "~183 mm, full width"),
  full_page     = list(width = 7.2, height = 9.0, dpi = 300,
                       note = "full page portrait"),
  slide         = list(width = 10.0, height = 5.6, dpi = 200,
                       note = "16:9 presentation")
)

#' Save at a named journal size.
mfg_save_figure <- function(p, name, size = "double_column", stage = "plot",
                            device = "png", height = NULL, width = NULL) {
  s <- MFG_FIGURE_SIZES[[size]]
  if (is.null(s)) {
    stop("Unknown figure size '", size, "'. Available: ",
         paste(names(MFG_FIGURE_SIZES), collapse = ", "), call. = FALSE)
  }
  mfg_save_plot(p, name, stage = stage, width = width %||% s$width,
                height = height %||% s$height, dpi = s$dpi, device = device)
}

# ── Alpha diversity ──────────────────────────────────────────────────────────

#' Alpha diversity boxplot with points and optional significance annotation.
#'
#' Individual points over the boxplot by default, because a boxplot of n = 5 that
#' looks like a boxplot of n = 500 misrepresents the evidence. The reader should
#' be able to count the samples.
#'
#' `test_result` accepts the object from run_alpha_test() so the annotation shows
#' the test that was actually run, rather than stat_compare_means picking its own
#' default — which is how a figure ends up disagreeing with the results text.
plot_alpha <- function(alpha_df, measure, group_var, facet_var = NULL,
                       test_result = NULL, show_points = TRUE,
                       comparisons = NULL, label = "p.format",
                       caption_width_in = 7.2, compact = FALSE) {
  df <- alpha_df[stats::complete.cases(alpha_df[[measure]], alpha_df[[group_var]]), ]
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[group_var]], y = .data[[measure]])) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = .data[[group_var]]),
                          outlier.shape = if (show_points) NA else 16,
                          outlier.size = 1, alpha = 0.75, width = 0.6)
  if (show_points) {
    p <- p + ggplot2::geom_jitter(width = 0.18, height = 0, size = 1.1,
                                  shape = 21, fill = "white", alpha = 0.7)
  }
  p <- p + scale_fill_canis() +
    ggplot2::labs(x = NULL, y = measure,
                  title = sprintf("%s diversity by %s", measure, group_var)) +
    theme_canis() + ggplot2::theme(legend.position = "none")

  if (!is.null(facet_var)) p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_var)))

  # n per group under each box — reference/10 requires n to be visible.
  n_lab <- as.data.frame(table(df[[group_var]]))
  names(n_lab) <- c(group_var, "n")
  p <- p + ggplot2::geom_text(data = n_lab, inherit.aes = FALSE,
      ggplot2::aes(x = .data[[group_var]], y = -Inf, label = paste0("n=", .data$n)),
      vjust = -0.6, size = 2.6, colour = canis_ui$muted)

  if (!is.null(test_result) && !identical(test_result$test, "REFUSED")) {
    # A panel gives each plot roughly total_width/ncol inches, so the full
    # subtitle and the reasoning caption overlap their neighbours. The compact
    # form keeps what identifies the test and drops what repeats across panels;
    # plot_alpha_panel() prints the reasoning once underneath instead.
    sub <- if (compact) {
      sprintf("%s, p = %s", test_result$test, mfg_fmt_p(test_result$p_value))
    } else {
      sprintf("%s: statistic = %.3f, p = %s, %s = %.3f",
              test_result$test, test_result$statistic,
              mfg_fmt_p(test_result$p_value),
              test_result$effect_name, test_result$effect_size)
    }
    p <- p + ggplot2::labs(subtitle = sub,
                           caption = if (compact) NULL else
                             mfg_wrap_caption(test_result$test_reason,
                                              caption_width_in))
  } else if (!is.null(comparisons)) {
    p <- p + ggpubr::stat_compare_means(comparisons = comparisons, label = label,
                                        hide.ns = FALSE)
  }
  p
}

#' Panel of every index at once.
#'
#' Each panel gets a compact subtitle; the test-selection reasoning — identical
#' across indices, because the same design drives it — is printed once beneath
#' the panel rather than repeated and clipped in every facet.
plot_alpha_panel <- function(alpha_df, measures, group_var, results = NULL,
                             ncol = 2, total_width_in = 7.2, titles = TRUE) {
  plots <- lapply(measures, function(m) {
    p <- plot_alpha(alpha_df, m, group_var, compact = TRUE,
                    test_result = if (!is.null(results)) results[[m]] else NULL)
    if (!titles) p <- p + ggplot2::labs(title = NULL)
    p
  })
  out <- patchwork::wrap_plots(plots, ncol = ncol)

  reason <- NULL
  if (!is.null(results)) {
    rs <- Filter(Negate(is.null), lapply(results[measures], function(r) r$test_reason))
    if (length(rs)) reason <- unique(unlist(rs))[1]
  }
  if (!is.null(reason)) {
    out <- out + patchwork::plot_annotation(
      caption = mfg_wrap_caption(reason, total_width_in),
      theme = ggplot2::theme(
        plot.caption = ggplot2::element_text(colour = canis_ui$muted,
                                             size = 9, hjust = 0)))
  }
  out
}

# ── Ordination ───────────────────────────────────────────────────────────────

#' Ordination scatter with the diagnostics in the subtitle.
#'
#' Stress for NMDS, variance explained for PCoA/PCA, constrained fraction for
#' CCA/RDA. These belong on the figure, not only in the text: an ordination
#' without them cannot be judged, and a two-axis plot of a 17%-explained
#' ordination looks identical to one explaining 80%.
#'
#' Ellipses are 95% confidence ellipses of the group centroid by default, not
#' data ellipses — they describe where the mean is, which is what a group
#' comparison is about.
plot_ordination_std <- function(ord, ps, color_var = NULL, shape_var = NULL,
                                ellipse = TRUE, ellipse_level = 0.95,
                                ellipse_type = "t", label_samples = FALSE,
                                permanova = NULL, caption_width_in = 7.2) {
  scores <- ord$scores
  names(scores)[1:2] <- c("Axis1", "Axis2")
  meta <- mfg_meta(ps)
  df <- cbind(scores[, 1:2, drop = FALSE],
              meta[rownames(scores), , drop = FALSE])
  df$.sample <- rownames(scores)

  aes_args <- list(x = quote(.data$Axis1), y = quote(.data$Axis2))
  if (!is.null(color_var)) aes_args$colour <- bquote(.data[[.(color_var)]])
  if (!is.null(shape_var)) aes_args$shape  <- bquote(.data[[.(shape_var)]])

  p <- ggplot2::ggplot(df, do.call(ggplot2::aes, aes_args)) +
    ggplot2::geom_point(size = 2.4, alpha = 0.85)

  if (isTRUE(ellipse) && !is.null(color_var)) {
    # An ellipse needs enough points per group to mean anything.
    ok <- names(which(table(df[[color_var]]) >= 4))
    if (length(ok)) {
      p <- p + ggplot2::stat_ellipse(data = df[df[[color_var]] %in% ok, ],
                                     level = ellipse_level, type = ellipse_type,
                                     linewidth = 0.5)
    }
  }
  if (isTRUE(label_samples)) {
    p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = .data$.sample),
                                      size = 2.2, max.overlaps = 20)
  }

  ax <- ord$axis_labels %||% c("Axis 1", "Axis 2")
  sub <- switch(ord$method,
    NMDS = sprintf("NMDS on %s distance | stress = %.3f (%s)",
                   ord$distance, ord$stress, ord$stress_verdict),
    PCoA = sprintf("PCoA on %s distance | axes 1-2 explain %.1f%% of variation",
                   ord$distance, sum(ord$variance_explained)),
    sprintf("%s | axes 1-2 explain %.1f%%%s", ord$method,
            sum(ord$variance_explained),
            if (!is.null(ord$pct_constrained) && !is.na(ord$pct_constrained))
              sprintf("; %.1f%% of total inertia constrained", ord$pct_constrained) else ""))

  cap <- NULL
  if (!is.null(permanova)) {
    cap <- paste(permanova_sentence(permanova), collapse = " ")
  }

  p + scale_color_canis() +
    ggplot2::labs(x = ax[1], y = ax[2],
                  title = sprintf("%s ordination", ord$method),
                  subtitle = sub,
                  caption = mfg_wrap_caption(cap, caption_width_in),
                  colour = color_var, shape = shape_var) +
    theme_canis() +
    ggplot2::coord_fixed(ratio = 1)
}

# ── Composition ──────────────────────────────────────────────────────────────

#' Stacked relative abundance bar plot.
#'
#' Always on relative abundance — a stacked bar of raw counts shows sequencing
#' depth, not composition. The "Other" share is reported in the caption, because
#' a plot showing the top 10 of 140 taxa is not showing the community.
plot_stacked_bar <- function(ps, rank = "Phylum", top_n = 10, group_var = NULL,
                             sample_order_var = NULL, facet_var = NULL,
                             by_group_mean = FALSE) {
  cl <- collapse_to_top(ps, rank = rank, top_n = top_n)
  mat <- cl$matrix
  long <- data.frame(
    taxon = rep(rownames(mat), times = ncol(mat)),
    sample = rep(colnames(mat), each = nrow(mat)),
    abundance = as.vector(mat), stringsAsFactors = FALSE)
  meta <- mfg_meta(ps)
  long <- cbind(long, meta[long$sample, , drop = FALSE])

  # "Other" belongs last in the stack, not sorted among the named taxa.
  lv <- c(setdiff(rownames(mat), "Other"), intersect("Other", rownames(mat)))
  long$taxon <- factor(long$taxon, levels = rev(lv))

  if (isTRUE(by_group_mean) && !is.null(group_var)) {
    agg <- stats::aggregate(abundance ~ taxon + get(group_var), data = long, FUN = mean)
    names(agg)[2] <- group_var
    p <- ggplot2::ggplot(agg, ggplot2::aes(x = .data[[group_var]], y = .data$abundance,
                                           fill = .data$taxon))
    xlab <- group_var
  } else {
    if (!is.null(sample_order_var)) {
      ord <- order(long[[sample_order_var]])
      long$sample <- factor(long$sample, levels = unique(long$sample[ord]))
    }
    p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$sample, y = .data$abundance,
                                            fill = .data$taxon))
    xlab <- NULL
  }

  p <- p + ggplot2::geom_bar(stat = "identity", width = 0.95, colour = NA) +
    scale_fill_canis() +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                expand = c(0, 0)) +
    ggplot2::labs(x = xlab, y = "Relative abundance", fill = rank,
                  title = sprintf("%s composition (top %d)", rank, top_n),
                  caption = sprintf("%d further %s collapsed into 'Other' (mean %.1f%% per sample)",
                                    length(cl$collapsed), tolower(rank), cl$mean_pct_other)) +
    theme_canis()

  if (is.null(group_var) || !isTRUE(by_group_mean)) {
    p <- p + sample_axis_theme(ncol(mat))
  }
  if (!is.null(facet_var)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_var)),
                                 scales = "free_x")
  }
  p
}

#' Prevalence versus abundance, for choosing a filter threshold visibly.
plot_prevalence <- function(ps, rank = "Phylum", threshold = NULL) {
  df <- qc_prevalence_table(ps, rank = rank)
  r <- safe_rank(ps, rank)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$total_abundance,
                                        y = .data$prevalence_frac)) +
    ggplot2::geom_point(ggplot2::aes(colour = if (!is.null(r)) .data[[r]] else NULL),
                        alpha = 0.7, size = 1.8) +
    ggplot2::scale_x_log10(labels = scales::comma) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_color_canis() +
    ggplot2::labs(x = "Total abundance (log scale)", y = "Prevalence",
                  colour = r, title = "Prevalence versus abundance",
                  subtitle = "Choose the prevalence filter from this shape, before testing") +
    theme_canis()
  if (!is.null(threshold)) {
    p <- p + ggplot2::geom_hline(yintercept = threshold, linetype = "dashed",
                                 colour = canis_ui$danger) +
      ggplot2::annotate("text", x = Inf, y = threshold, hjust = 1.05, vjust = -0.6,
                        label = sprintf("filter at %.0f%%", 100 * threshold),
                        size = 3, colour = canis_ui$danger)
  }
  p
}

#' Rarefaction curves with the chosen depth marked.
plot_rarefaction <- function(curve_obj, depth = NULL, group_map = NULL) {
  df <- curve_obj$curves
  if (!is.null(group_map)) df$group <- group_map[df$sample]
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$depth, y = .data$richness,
                                        group = .data$sample)) +
    ggplot2::geom_line(ggplot2::aes(colour = if (!is.null(group_map)) .data$group else NULL),
                       alpha = 0.6, linewidth = 0.4) +
    scale_color_canis() +
    ggplot2::scale_x_continuous(labels = scales::comma) +
    ggplot2::labs(x = "Sequencing depth", y = "Observed taxa",
                  colour = if (!is.null(group_map)) "Group" else NULL,
                  title = "Rarefaction curves",
                  subtitle = sprintf("%d of %d samples plateau (terminal slope < %.2g taxa per 1000 reads; median %.2f)",
                                     curve_obj$n_plateaued, curve_obj$n_total,
                                     curve_obj$plateau_threshold %||% 1,
                                     curve_obj$median_terminal_slope %||% NA_real_)) +
    theme_canis()
  if (!is.null(depth)) {
    p <- p + ggplot2::geom_vline(xintercept = depth, linetype = "dashed",
                                 colour = canis_ui$danger) +
      ggplot2::annotate("text", x = depth, y = Inf, vjust = 1.6, hjust = -0.08,
                        label = sprintf("depth = %s", format(depth, big.mark = ",")),
                        size = 3, colour = canis_ui$danger)
  }
  p
}

# ── Differential abundance ───────────────────────────────────────────────────

#' Volcano plot.
#'
#' Effect size on x, significance on y. The deliberate choices: the y axis is
#' -log10 of the ADJUSTED p-value, not the raw one, so the visual threshold is
#' the one the claim rests on; and both thresholds are drawn, so a reader can see
#' what "significant" meant here.
plot_volcano <- function(da_result, effect_col = NULL, q_col = NULL,
                         q_cutoff = NULL, effect_cutoff = 1,
                         label_top = 10, label_col = NULL) {
  res <- da_result$results
  q_cutoff <- q_cutoff %||% da_result$alpha %||% 0.05

  if (is.null(effect_col)) {
    effect_col <- switch(da_result$method,
      deseq2 = if ("log2FoldChange_shrunk" %in% names(res)) "log2FoldChange_shrunk"
               else "log2FoldChange",
      aldex2 = da_result$effect_column %||% NA_character_,
      ancombc2 = grep("^lfc_(?!\\(Intercept\\))", names(res), value = TRUE, perl = TRUE)[1],
      "effect_size")
  }
  if (is.null(q_col)) {
    q_col <- switch(da_result$method,
      deseq2 = "padj",
      ancombc2 = grep("^q_(?!\\(Intercept\\))", names(res), value = TRUE, perl = TRUE)[1],
      "q_value")
  }
  if (length(effect_col) != 1 || is.na(effect_col) || !effect_col %in% names(res)) {
    stop("No effect-size column is available for this result", 
         if (identical(da_result$method, "aldex2"))
           " — ALDEx2 reports one only for two-group comparisons. Use plot_da_effects() with a q-value ordering, or run the pairwise contrasts."
         else ". Pass effect_col explicitly.", call. = FALSE)
  }
  if (length(q_col) != 1 || is.na(q_col) || !q_col %in% names(res)) {
    stop("Could not find q-value column '", q_col, "' in the result table. ",
         "Pass q_col explicitly.", call. = FALSE)
  }

  df <- data.frame(taxon = res$taxon, effect = res[[effect_col]],
                   q = res[[q_col]], stringsAsFactors = FALSE)
  df <- df[stats::complete.cases(df$effect, df$q), ]
  df$neglog10q <- -log10(pmax(df$q, .Machine$double.xmin))
  df$status <- ifelse(df$q >= q_cutoff, "not significant",
                ifelse(abs(df$effect) < effect_cutoff, "significant, small effect",
                       "significant"))
  if (!is.null(label_col) && label_col %in% names(res)) {
    df$label <- res[[label_col]][match(df$taxon, res$taxon)]
  } else df$label <- df$taxon

  lim <- max(abs(df$effect), na.rm = TRUE)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$effect, y = .data$neglog10q)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$status), size = 1.8, alpha = 0.8) +
    ggplot2::geom_hline(yintercept = -log10(q_cutoff), linetype = "dashed",
                        colour = canis_ui$muted, linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = c(-effect_cutoff, effect_cutoff),
                        linetype = "dotted", colour = canis_ui$muted, linewidth = 0.4) +
    ggplot2::scale_colour_manual(values = c(
      "not significant" = "#BBBBBB",
      "significant, small effect" = canis_categorical[2],
      "significant" = canis_categorical[6])) +
    ggplot2::xlim(-lim, lim) +
    ggplot2::labs(x = effect_col, y = sprintf("-log10(%s)", q_col),
      colour = NULL,
      title = sprintf("Differential abundance: %s", da_result$method),
      subtitle = sprintf("%d of %d taxa significant at %s < %.3g",
                         sum(df$q < q_cutoff, na.rm = TRUE), nrow(df), q_col, q_cutoff),
      caption = sprintf("Dashed line: %s = %.3g. Dotted lines: |effect| = %g. y axis uses the adjusted p-value.",
                        q_col, q_cutoff, effect_cutoff)) +
    theme_canis()

  if (label_top > 0) {
    top <- utils::head(df[order(df$q), ], label_top)
    p <- p + ggrepel::geom_text_repel(data = top,
      ggplot2::aes(label = .data$label), size = 2.4, max.overlaps = 20,
      colour = canis_ui$ink, min.segment.length = 0)
  }
  p
}

#' Effect sizes of the significant taxa, sorted — usually more informative than a volcano.
plot_da_effects <- function(da_result, top_n = 25, label_col = NULL) {
  sig <- da_significant(da_result, n = top_n)
  if (!nrow(sig)) {
    return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0,
      label = "No significantly differentially abundant taxa") +
      theme_canis() + ggplot2::theme(axis.text = ggplot2::element_blank(),
                                     panel.grid = ggplot2::element_blank()) +
      ggplot2::labs(x = NULL, y = NULL))
  }
  if (!is.null(label_col) && label_col %in% names(da_result$results)) {
    sig$label <- da_result$results[[label_col]][match(sig$taxon, da_result$results$taxon)]
    # Several ASVs commonly share a genus, and an unclassified rank gives NA.
    # Both would collapse rows onto one axis position, so the taxon id
    # disambiguates and NA falls back to it.
    sig$label[is.na(sig$label) | !nzchar(sig$label)] <- sig$taxon[is.na(sig$label) |
                                                                 !nzchar(sig$label)]
    dup <- sig$label %in% sig$label[duplicated(sig$label)]
    sig$label[dup] <- paste0(sig$label[dup], " (", sig$taxon[dup], ")")
  } else sig$label <- sig$taxon
  sig$label <- factor(sig$label, levels = unique(sig$label[order(sig$effect)]))

  ggplot2::ggplot(sig, ggplot2::aes(x = .data$effect, y = .data$label,
                                    fill = .data$effect > 0)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_vline(xintercept = 0, colour = canis_ui$ink, linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = canis_categorical[1],
                                          `FALSE` = canis_categorical[6]),
                               guide = "none") +
    ggplot2::labs(x = "Effect size", y = NULL,
      title = sprintf("Differentially abundant taxa (%s)", da_result$method),
      subtitle = sprintf("Top %d of %d significant, ordered by effect",
                         nrow(sig), da_result$n_significant[1] %||% nrow(sig))) +
    theme_canis()
}

# ── Model output ─────────────────────────────────────────────────────────────

#' Coefficient plot with confidence intervals.
#'
#' Effect sizes with uncertainty, which is what reference/10 asks for and what a
#' table of p-values does not give.
plot_model_coefficients <- function(model, exclude_intercept = TRUE, level = 0.95) {
  fit <- if (inherits(model, "mfg_model")) model$fit else model
  ci <- try(stats::confint(fit, level = level), silent = TRUE)
  co <- if (inherits(model, "mfg_model")) model$coefficients else summary(fit)$coefficients
  est <- co[, 1]
  df <- data.frame(term = rownames(co), estimate = est, stringsAsFactors = FALSE)
  if (!inherits(ci, "try-error")) {
    ci <- as.data.frame(ci)
    df$lower <- ci[df$term, 1]; df$upper <- ci[df$term, 2]
  } else {
    se <- co[, 2]
    z <- stats::qnorm(1 - (1 - level) / 2)
    df$lower <- est - z * se; df$upper <- est + z * se
  }
  if (exclude_intercept) df <- df[!grepl("Intercept", df$term), , drop = FALSE]
  df$crosses_zero <- df$lower <= 0 & df$upper >= 0
  df$term <- factor(df$term, levels = df$term[order(df$estimate)])

  ggplot2::ggplot(df, ggplot2::aes(x = .data$estimate, y = .data$term)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = canis_ui$muted) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$lower, xmax = .data$upper,
                                         colour = .data$crosses_zero),
                            height = 0.2, linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$crosses_zero), size = 2.4) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = canis_categorical[1],
                                            `TRUE` = "#BBBBBB"), guide = "none") +
    ggplot2::labs(x = sprintf("Estimate (%.0f%% CI)", 100 * level), y = NULL,
      title = "Model coefficients",
      subtitle = paste(deparse(model$formula %||% stats::formula(fit)), collapse = ""),
      caption = "Grey intervals cross zero. Coefficients are relative to the reference level of each factor.") +
    theme_canis()
}

#' Count and zero distribution, the Demo 11 / ps11 diagnostic pair.
plot_zero_distribution <- function(counts, taxon_name = "taxon") {
  d <- as.data.frame(table(counts), stringsAsFactors = FALSE)
  names(d) <- c("count", "frequency")
  d$count <- as.numeric(d$count)
  p1 <- ggplot2::ggplot(d, ggplot2::aes(x = .data$count, y = .data$frequency)) +
    ggplot2::geom_col(fill = canis_categorical[1], width = 0.8) +
    ggplot2::labs(x = "Read count", y = "Number of samples",
                  title = sprintf("%s: count distribution", taxon_name),
                  subtitle = sprintf("%d of %d samples have zero reads (%.0f%%)",
                    sum(counts == 0), length(counts), 100 * mean(counts == 0))) +
    theme_canis()
  p2 <- ggplot2::ggplot(data.frame(i = seq_along(counts), v = sort(counts)),
                        ggplot2::aes(x = .data$i, y = .data$v)) +
    ggplot2::geom_point(size = 0.9, colour = canis_categorical[2]) +
    ggplot2::labs(x = "Sample (sorted)", y = "Observed read count",
                  title = "Sorted counts") +
    theme_canis()
  patchwork::wrap_plots(p1, p2, ncol = 2)
}

# ── Dispersion ───────────────────────────────────────────────────────────────

#' Distance-to-centroid by group — the visual form of the betadisper check.
plot_dispersion <- function(disp, group_name = "group", caption_width_in = 7.2) {
  bd <- disp$betadisper
  df <- data.frame(distance = bd$distances, group = bd$group, stringsAsFactors = FALSE)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$group, y = .data$distance,
                                   fill = .data$group)) +
    ggplot2::geom_boxplot(alpha = 0.75, width = 0.6, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.18, size = 1, shape = 21, fill = "white", alpha = 0.7) +
    scale_fill_canis() +
    ggplot2::labs(x = NULL, y = "Distance to group centroid",
      title = "Multivariate dispersion",
      subtitle = sprintf("permutest p = %s — %s", mfg_fmt_p(disp$p_value),
                         if (disp$homogeneous) "homogeneous" else "HETEROGENEOUS"),
      caption = mfg_wrap_caption(disp$interpretation, caption_width_in)) +
    theme_canis() + ggplot2::theme(legend.position = "none")
}
