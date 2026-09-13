# =============================================================================
# MicroFitGut — 09-exploratory.R
#
# The looking-at-the-data stage: heat trees, CLR heatmaps, dendrograms,
# correlation networks, phylogenies, sample networks.
#
# One rule governs this whole file, and it is the reason exploratory analysis has
# its own module rather than being mixed into the testing stages: anything found
# here is a hypothesis, not a result. Choosing which comparison to test after
# seeing a heatmap is circular, and the p-value from that test does not mean what
# it says. Exploration comes first, the analysis plan is fixed after it, and the
# report says which is which.
#
# Correlation networks carry a second, sharper problem: correlations between
# relative abundances are spurious by construction, because proportions must sum
# to one. See the note on cooccurrence_network().
#
# Requires: mfg_require(c("exploratory","plots")); utils.R; 03-normalize.R
# =============================================================================

# ── Heat tree ────────────────────────────────────────────────────────────────

#' Metacoder heat tree: taxonomy as a tree with abundance mapped to size/colour.
#'
#' Good at showing where in the taxonomy the signal sits, which a bar plot cannot.
#' Needs a taxonomy with several populated ranks to be worth anything.
plot_heat_tree <- function(ps, rank_limit = "Genus", top_n = 100,
                           node_label_max = 30, seed = 42) {
  rel <- if (identical(mfg_normalization(ps), "tss")) ps else tss_transform(ps)

  # Heat trees become unreadable past a few dozen nodes, so the set is trimmed
  # explicitly rather than letting the plot degrade silently.
  keep <- names(sort(phyloseq::taxa_sums(rel), decreasing = TRUE))
  keep <- keep[seq_len(min(top_n, length(keep)))]
  rel  <- phyloseq::prune_taxa(keep, rel)

  tt <- mfg_tax(rel)
  ranks <- intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                     names(tt))
  if (!is.null(rank_limit) && rank_limit %in% ranks) {
    ranks <- ranks[seq_len(which(ranks == rank_limit))]
  }
  if (length(ranks) < 2) {
    stop("A heat tree needs at least two populated taxonomic ranks; this object has: ",
         paste(ranks, collapse = ", "), call. = FALSE)
  }

  lineage <- apply(tt[, ranks, drop = FALSE], 1, function(r) {
    r <- r[!is.na(r) & nzchar(r)]
    paste(paste0(substr(tolower(ranks[seq_along(r)]), 1, 1), "__", r), collapse = ";")
  })
  df <- data.frame(taxon_id = rownames(tt), lineage = lineage,
                   stringsAsFactors = FALSE)
  mat <- mfg_otu_taxa_as_rows(rel)
  df <- cbind(df, as.data.frame(mat[df$taxon_id, , drop = FALSE]))

  obj <- metacoder::parse_tax_data(df, class_cols = "lineage", class_sep = ";",
    class_regex = "^([a-z])__(.*)$", class_key = c(tax_rank = "taxon_rank",
                                                   tax_name = "taxon_name"))
  obj$data$tax_abund <- metacoder::calc_taxon_abund(obj, "tax_data",
    cols = phyloseq::sample_names(rel))

  set.seed(seed)
  p <- metacoder::heat_tree(obj,
    node_label = metacoder::taxon_names(obj),
    node_size = obj$data$tax_abund$total %||% rowSums(
      obj$data$tax_abund[, -1, drop = FALSE]),
    node_color = obj$data$tax_abund$total %||% rowSums(
      obj$data$tax_abund[, -1, drop = FALSE]),
    node_size_axis_label = "Abundance",
    node_color_axis_label = "Abundance",
    layout = "davidson-harel", initial_layout = "reingold-tilford")

  mfg_log("exploratory", "heat_tree",
          list(ranks = paste(ranks, collapse = ","), n_taxa_shown = length(keep),
               seed = seed,
               note = "exploratory only; do not select comparisons from this figure"))
  p
}

# ── Heatmap ──────────────────────────────────────────────────────────────────

#' CLR heatmap of the most variable taxa.
#'
#' CLR rather than raw or relative abundance because a heatmap of proportions is
#' dominated by the few abundant taxa and shows nothing else. Selecting by
#' variance rather than by mean is what makes the picture informative — but note
#' that selecting the most variable taxa and then testing those same taxa is
#' circular; this is a display, not a screen.
plot_clr_heatmap <- function(ps, top_n = 40, rank = NULL,
                            annotation_vars = NULL,
                            cluster_rows = TRUE, cluster_cols = TRUE,
                            show_colnames = NULL) {
  work <- if (!is.null(rank)) agglomerate_and_summarise(ps, rank)$ps else ps
  cl   <- if (identical(mfg_normalization(work), "clr")) work else clr_transform(work)
  mat  <- mfg_otu_taxa_as_rows(cl)

  v <- apply(mat, 1, stats::var)
  keep <- names(sort(v, decreasing = TRUE))[seq_len(min(top_n, length(v)))]
  mat <- mat[keep, , drop = FALSE]

  # Label rows by the deepest resolved name rather than an ASV id, where possible.
  tt <- mfg_tax(cl)
  if (!is.null(tt)) {
    lbl <- apply(tt[keep, , drop = FALSE], 1, function(r) {
      r <- r[!is.na(r) & nzchar(r) & !r %in% c("unclassified", "Unclassified")]
      if (length(r)) utils::tail(r, 1) else NA_character_
    })
    lbl[is.na(lbl)] <- keep[is.na(lbl)]
    rownames(mat) <- make.unique(paste0(lbl, " (", keep, ")"))
  }

  ann <- NULL
  if (!is.null(annotation_vars)) {
    meta <- mfg_meta(cl)
    av <- intersect(annotation_vars, names(meta))
    if (length(av)) {
      ann <- meta[colnames(mat), av, drop = FALSE]
      ann[] <- lapply(ann, function(x) if (is.numeric(x)) x else as.factor(x))
    }
  }

  show_colnames <- show_colnames %||% (ncol(mat) <= 40)
  p <- pheatmap::pheatmap(mat, scale = "none",
    cluster_rows = cluster_rows, cluster_cols = cluster_cols,
    annotation_col = ann, show_colnames = show_colnames,
    fontsize_row = 7, fontsize_col = 7, silent = TRUE,
    main = sprintf("Top %d most variable taxa (CLR)", nrow(mat)))

  mfg_log("exploratory", "clr_heatmap",
          list(top_n = nrow(mat), rank = rank %||% "as-is",
               selected_by = "variance across samples",
               annotations = paste(annotation_vars, collapse = ","),
               note = "display only; selecting by variance then testing the same taxa is circular"))
  p
}

# ── Dendrogram ───────────────────────────────────────────────────────────────

#' Coloured dendrogram of samples, with a known grouping overlaid.
#'
#' The useful reading is whether the clustering recovers the design. Returns the
#' adjusted Rand index alongside, so "the groups separate" becomes a number.
plot_dendrogram <- function(ps, distance = "bray", method = "ward.D2",
                            group_var = NULL, k = NULL) {
  d  <- compute_distance(ps, distance)
  hc <- run_hclust(d, method = method, k = k)
  dend <- stats::as.dendrogram(hc$hclust)

  ari <- NA_real_
  if (!is.null(group_var)) {
    meta <- mfg_meta(ps)
    g <- as.factor(meta[[group_var]])
    names(g) <- rownames(meta)
    lab_order <- labels(dend)
    cols <- canis_colors(length(levels(g)))[as.integer(g[lab_order])]
    dend <- dendextend::set(dend, "labels_colors", cols)
    dend <- dendextend::set(dend, "labels_cex", 0.5)
    if (!is.null(hc$clusters)) {
      ari <- cluster_vs_group(hc$clusters, g[names(hc$clusters)])$adjusted_rand_index
    }
  }
  mfg_log("exploratory", "dendrogram",
          list(distance = distance, linkage = method, k = k %||% NA,
               cophenetic_correlation = round(hc$cophenetic_correlation, 3),
               adjusted_rand_index = if (is.na(ari)) NA else round(ari, 3)))
  list(dendrogram = dend, hclust = hc, adjusted_rand_index = ari,
       cophenetic_correlation = hc$cophenetic_correlation)
}

# ── Correlation network ──────────────────────────────────────────────────────

#' Taxon co-occurrence network from Spearman correlations on CLR values.
#'
#' The reason CLR is not optional here: correlations between relative abundances
#' are spurious by construction. Proportions sum to one, so if one taxon rises
#' the others must fall, and a naive correlation matrix on proportions produces
#' strong negative correlations that reflect arithmetic rather than ecology.
#' CLR removes the constraint. This is the single most common error in microbiome
#' network papers, and it is why this function refuses non-CLR input.
#'
#' Even on CLR values, treat the edges as hypotheses: SpiecEasi or SPARCC are the
#' methods designed for this problem, and a Spearman network is a first look.
cooccurrence_network <- function(ps, min_abs_rho = 0.5, p_cutoff = 0.05,
                                p_adjust = "BH", top_n = 60, rank = NULL,
                                require_clr = TRUE) {
  work <- if (!is.null(rank)) agglomerate_and_summarise(ps, rank)$ps else ps

  if (isTRUE(require_clr) && !identical(mfg_normalization(work), "clr")) {
    if (identical(mfg_normalization(work), "tss")) {
      stop("Refusing to build a correlation network on relative abundances.\n",
           "Proportions are constrained to sum to 1, so correlations between them ",
           "are spurious by construction — you will get strong negative edges that ",
           "are arithmetic, not ecology.\nCLR-transform first: ",
           "cooccurrence_network(clr_transform(ps), ...)", call. = FALSE)
    }
    work <- clr_transform(work)
    message("CLR-transformed before correlating, because proportions cannot be correlated directly.")
  }

  mat <- mfg_otu_taxa_as_rows(work)
  keep <- names(sort(apply(mat, 1, stats::var), decreasing = TRUE))
  keep <- keep[seq_len(min(top_n, length(keep)))]
  mat <- mat[keep, , drop = FALSE]

  ct <- stats::cor(t(mat), method = "spearman")
  # p-values for every pair, then corrected: a 60-taxon network is 1770 tests.
  n <- ncol(mat)
  tstat <- ct * sqrt((n - 2) / pmax(1 - ct^2, 1e-12))
  pv <- 2 * stats::pt(-abs(tstat), df = n - 2)
  diag(pv) <- 1
  upper <- upper.tri(pv)
  q <- rep(NA_real_, length(pv)); dim(q) <- dim(pv)
  q[upper] <- stats::p.adjust(pv[upper], method = p_adjust)
  q[lower.tri(q)] <- t(q)[lower.tri(q)]

  adj <- (abs(ct) >= min_abs_rho) & (q <= p_cutoff)
  diag(adj) <- FALSE
  n_tests <- sum(upper)

  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  igraph::E(g)$rho <- ct[igraph::as_edgelist(g, names = FALSE)]
  igraph::E(g)$sign <- ifelse(igraph::E(g)$rho > 0, "positive", "negative")
  igraph::V(g)$degree <- igraph::degree(g)

  edges <- data.frame(
    from = igraph::as_edgelist(g)[, 1], to = igraph::as_edgelist(g)[, 2],
    rho = igraph::E(g)$rho, sign = igraph::E(g)$sign, stringsAsFactors = FALSE)
  edges$q_value <- q[cbind(match(edges$from, rownames(q)), match(edges$to, colnames(q)))]
  edges <- edges[order(-abs(edges$rho)), ]

  mfg_log("exploratory", "cooccurrence_network", list(
    n_taxa = length(keep), n_pairs_tested = n_tests,
    min_abs_rho = min_abs_rho, p_adjust = p_adjust, p_cutoff = p_cutoff,
    n_edges = nrow(edges),
    n_positive = sum(edges$sign == "positive"),
    n_negative = sum(edges$sign == "negative"),
    normalization = mfg_normalization(work),
    caveat = "Spearman on CLR is a first look; use SpiecEasi/SPARCC for inference"))

  structure(list(graph = g, edges = edges, correlation = ct, q_value = q,
    n_taxa = length(keep), n_pairs_tested = n_tests, n_edges = nrow(edges),
    min_abs_rho = min_abs_rho, p_adjust = p_adjust, p_cutoff = p_cutoff,
    caveat = paste("Edges are Spearman correlations on CLR values with", p_adjust,
      "correction across", n_tests, "pairs. They are hypotheses about",
      "co-occurrence, not inferred interactions — use SpiecEasi or SPARCC if the",
      "network itself is the finding.")),
    class = c("mfg_network", "list"))
}

#' @export
print.mfg_network <- function(x, ...) {
  cat("=== Co-occurrence network ===\n")
  cat(sprintf("%d taxa, %d pairs tested, %d edges retained (|rho| >= %.2f, %s q <= %.3f)\n",
              x$n_taxa, x$n_pairs_tested, x$n_edges, x$min_abs_rho, x$p_adjust, x$p_cutoff))
  if (x$n_edges) {
    cat(sprintf("  %d positive, %d negative\n",
                sum(x$edges$sign == "positive"), sum(x$edges$sign == "negative")))
    cat("\nStrongest edges:\n")
    print(utils::head(x$edges, 10), row.names = FALSE)
  }
  cat("\n"); cat(strwrap(x$caveat, width = 78, prefix = "! "), sep = "\n")
  invisible(x)
}

#' Sample-similarity network, as in Demo 6.
#'
#' Nodes are samples, edges join samples closer than max.dist. Lowering max.dist
#' reduces the number of connected nodes; the threshold is arbitrary and must be
#' stated, because the picture changes completely with it.
sample_network <- function(ps, distance = "bray", max_dist = 0.6,
                          color_var = NULL, shape_var = NULL) {
  g <- phyloseq::make_network(ps, distance = distance, max.dist = max_dist)
  n_nodes <- igraph::gorder(g)
  mfg_log("exploratory", "sample_network",
          list(distance = distance, max_dist = max_dist,
               n_samples = phyloseq::nsamples(ps), n_connected = n_nodes,
               n_edges = igraph::gsize(g),
               note = "max.dist is arbitrary and changes the picture; state it"))
  list(graph = g, n_connected = n_nodes, n_total = phyloseq::nsamples(ps),
       max_dist = max_dist, distance = distance)
}

# ── Phylogeny ────────────────────────────────────────────────────────────────

#' Plot the phylogeny with abundance or taxonomy annotation.
#'
#' Refuses a placeholder tree, for the same reason UniFrac does: a picture of a
#' random topology presented as a phylogeny is a false claim about evolutionary
#' relationship, and it looks exactly like a real one.
plot_phylogeny <- function(ps, color_rank = "Phylum", layout = "rectangular",
                           show_tip_labels = NULL, require_real_tree = TRUE) {
  tree <- phyloseq::phy_tree(ps, errorIfNULL = FALSE)
  if (is.null(tree)) stop("No phylogenetic tree in this object.", call. = FALSE)
  real <- mfg_tree_is_real(ps)
  if (isTRUE(require_real_tree) && !isTRUE(real)) {
    stop("Refusing to plot the phylogeny: the tree is ",
         if (is.na(real)) "of unrecorded provenance" else "a placeholder",
         ".\nA figure of a random topology labelled as a phylogeny is a false claim ",
         "about evolutionary relationship. Confirm with mfg_mark_tree_real(ps, TRUE).",
         call. = FALSE)
  }
  show_tip_labels <- show_tip_labels %||% (length(tree$tip.label) <= 60)

  r <- safe_rank(ps, color_rank)
  p <- ggtree::ggtree(tree, layout = layout)
  if (!is.null(r)) {
    tt <- mfg_tax(ps)
    dat <- data.frame(label = rownames(tt), group = tt[[r]], stringsAsFactors = FALSE)
    p <- p %<+% dat + ggtree::geom_tippoint(ggplot2::aes(color = group), size = 2) +
      ggplot2::labs(color = r) + scale_color_canis()
  }
  if (show_tip_labels) p <- p + ggtree::geom_tiplab(size = 2)
  mfg_log("exploratory", "phylogeny_plotted",
          list(n_tips = length(tree$tip.label), layout = layout,
               colored_by = r %||% "none", rooted = ape::is.rooted(tree)))
  p
}

#' Build a tree from representative sequences.
#'
#' Neighbour-joining on ML distances, midpoint-rooted — the Demo 8 procedure. NJ
#' is fast and adequate for UniFrac; it is not a phylogenetic inference you would
#' publish as a phylogeny. Trees built here are marked real, because they were
#' estimated from the sequence data.
build_tree_nj <- function(seqs, model = "JC69", root = TRUE) {
  if (is.character(seqs) && length(seqs) == 1 && file.exists(seqs)) {
    seqs <- phangorn::read.phyDat(seqs, format = "fasta")
  }
  if (inherits(seqs, "DNAStringSet")) {
    aln <- DECIPHER::AlignSeqs(seqs, verbose = FALSE)
    seqs <- phangorn::as.phyDat(as(aln, "matrix"), type = "DNA")
  }
  dm  <- phangorn::dist.ml(seqs, model = model)
  tre <- ape::nj(dm)
  if (isTRUE(root)) tre <- phangorn::midpoint(tre)
  mfg_log("exploratory", "tree_built",
          list(method = "neighbour-joining", distance_model = model,
               n_tips = length(tre$tip.label), rooted = ape::is.rooted(tre),
               note = "adequate for UniFrac; not a publication phylogeny"))
  tre
}
