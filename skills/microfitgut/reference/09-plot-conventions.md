# Plot conventions

The house style, carried from CanisLupus 2.0. One palette and one theme applied
everywhere, so every figure in a run reads as part of the same set.

## Colour

**Categorical: Okabe-Ito.** `scale_fill_canis()` / `scale_color_canis()`.

```
#0072B2 blue     #E69F00 orange   #009E73 green    #CC79A7 pink
#56B4E9 sky      #D55E00 vermilion #F0E442 yellow  #5D3A9B purple
```

extended with further distinguishable hues for taxa-heavy plots.

**Continuous: viridis.** `scale_fill_canis_c()` / `scale_color_canis_c()`.
Perceptually uniform and CVD-safe.

**Diverging, centred on zero: symmetric blue–white–red.**
`scale_fill_canis_diverging()`. For log fold change, CLR values, correlations.
A diverging scale must be symmetric about its midpoint or it misreads the sign of
the data.

**Never a raw rainbow.** Its lightness is non-monotonic, so it invents banding
that is not in the data and makes ordered values look categorical.

### Why this matters more in a figure than in a logo

In a logo, colour carries no meaning — a spiral is a spiral whatever hue it is.
In an ordination or a stacked bar, **colour is the encoding**, so a palette that
collapses under colour-vision deficiency destroys the data rather than merely
looking different.

Two specific checks:

- A true green and a true red at similar lightness is the classic failure pair.
  Teal and orange, or blue and orange, survive deuteranopia and protanopia.
- **Greyscale separation is a separate requirement.** Two hues can be
  distinguishable in colour and have near-identical relative luminance, collapsing
  to the same tone in black-and-white print. If the target journal prints figures
  in greyscale, check luminance, not just hue — a 1.6:1 luminance ratio separates,
  1.15:1 does not.

Beyond about eight categories, distinguishability degrades no matter what
palette is used. That is why `collapse_to_top()` exists: collapse to a readable
number and report what went into `Other`.

## Dimensions

`MFG_FIGURE_SIZES`, passed via `mfg_save_figure(p, name, size)`:

| Name | Inches | Use |
|---|---|---|
| `single_column` | 3.5 × 3.0 | ~89 mm, most journals' single column |
| `onehalf_column` | 5.5 × 4.0 | ~140 mm |
| `double_column` | 7.2 × 5.0 | ~183 mm, full width (default) |
| `full_page` | 7.2 × 9.0 | full page portrait, for panels |
| `slide` | 10.0 × 5.6 | 16:9 presentation, 200 dpi |

300 dpi for print, 200 for slides.

**Size the figure before you write the text into it.** Text in a figure scales
with the saved dimensions, so a plot designed at 10 inches and printed at 3.5 has
unreadable axes. This is the commonest reason a figure comes back at revision.

`mfg_wrap_caption()` wraps caption text to what the chosen width can actually
show (~13 characters per inch at the caption's point size). A caption wrapped too
wide gets clipped at the page edge, which is worse than no caption because the
reader cannot tell it was truncated.

## What must appear on each figure type

| Figure | Must show |
|---|---|
| Alpha boxplot | **n per group**, individual points, the test actually run with its p and effect size |
| Ordination | variance explained per axis (PCoA/PCA) or **stress** (NMDS); the PERMANOVA result with its dispersion qualification |
| Stacked bar | relative abundance (never counts), the count collapsed into `Other` and its mean share |
| Volcano | −log10 of the **adjusted** p-value, both thresholds drawn, the significant/tested denominator |
| Prevalence plot | the filter threshold, if one was applied |
| Rarefaction curves | the chosen depth marked, and how many samples plateau |
| Dispersion plot | the betadisper p-value and its interpretation |
| Heatmap | what the values are (CLR/VST/log10) and how rows were selected |
| Coefficient plot | confidence intervals, and which intervals cross zero |
| Phylogeny | only if the tree is real (`01`) |

### Individual points over boxplots

Default `show_points = TRUE`. A boxplot of n = 5 looks identical to a boxplot of
n = 500, which misrepresents the evidence. The reader should be able to count the
samples. `plot_alpha()` also prints `n=` under each box.

### Significance annotation

`plot_alpha(ad, measure, group_var, test_result = <run_alpha_test output>)` is
preferred over bare `stat_compare_means()`, because it annotates **the test that
was actually run** rather than letting ggpubr pick its own default. That is how a
figure ends up disagreeing with the results text.

When using `stat_compare_means()` directly, pass the same method and comparisons
the results section reports.

### Ellipses

95% confidence ellipses of the group **centroid** by default (`type = "t"`), not
data ellipses — they describe where the mean is, which is what a group comparison
is about. Drawn only for groups with at least 4 points; below that an ellipse is
a shape, not an interval.

## Sample axis labels

`sample_axis_theme(n_samples, max_labels = 25)` hides per-sample labels beyond 25.
Unreadable labels are worse than none — they occupy space and suggest
information the reader cannot extract.

## Saving

Everything goes through `mfg_save_plot()` or `mfg_save_figure()`, which write to
`output/<run_id>/figures/` and **log the file**. A figure that is not logged
cannot be traced by `mfg_traceability_check()` or by the verifier subagent, and
under the rule in `10` it does not get discussed in the report.

Plot functions return ggplot objects rather than drawing, so the caller composes
panels (`patchwork::wrap_plots`, `cowplot::plot_grid`) before saving.

## File naming

`<stage>_<content>.png`: `alpha_shannon.png`, `ordination_pcoa.png`,
`da_volcano_region.png`, `composition_phylum.png`. Descriptive, lowercase,
underscore-separated — so a directory listing tells you what the run produced
without opening anything.
