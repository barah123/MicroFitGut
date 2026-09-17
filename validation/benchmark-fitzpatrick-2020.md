# Benchmark: Fitzpatrick & Schneider (2020), a parasitic plant and its host

Study **A1** of the validation matrix, intended as the negative control, and
the first study to exercise the phylogenetic distance path on a real tree.

| | |
|---|---|
| **Paper** | Fitzpatrick CR, Schneider AC. Unique bacterial assembly, composition, and interactions in a parasitic plant and its host. *J Exp Bot.* 2020;71(6):2198–2209. |
| **Data** | Dryad [10.5061/dryad.7wm37pvnk](https://doi.org/10.5061/dryad.7wm37pvnk), ASV table, RDP taxonomy, **phylogeny**, and the authors' analysis code |
| **Design** | *Orobanche hederae* holoparasite and its *Hedera* host, 4 sites, leaf/root/soil, infected and uninfected patches |
| **Verdict** | **Incomplete**: 2 claims held, 0 failed, 3 unadjudicated |

---

## 1. The matrix assignment was wrong, and that is useful

A1 was placed in the matrix as a **negative control**: a clean cross-sectional
design where no guard should fire. The validation plan flagged that designs were
inferred from titles and file listings and warned they might not hold.

They did not. A1 is a **paired infected/uninfected design across four sites with
multiple tissue types per plant**: nine groups in total:

| | leaf | root | infected root | soil |
|---|---|---|---|---|
| **Parasite (P)** | 12 | 12 |, |, |
| **Infected host (I)** | 12 | 12 | 12 | 8 |
| **Uninfected host (U)** | 12 | 12 |, | 7 |

Plus four technical controls (`mock`, `PAO1`, `water`, `Undetermined`), dropped
and logged.

So A1 is not a negative control. **The validation set still lacks one**, and that
matters: without a clean study where nothing fires, there is no evidence
distinguishing "the guards catch real problems" from "the guards fire on
everything". Recorded as a gap to fill.

What A1 *is* good for turns out to be more valuable: it has a real phylogeny
covering all 21,865 ASVs, which no earlier study in the series did.

---

## 2. Reconstruction

Minimum depth 1,000, prevalence 0.02: **21,865 → 9,488 taxa**, 99 samples.
Depth range 28,701–261,684 (9.1-fold). Sparsity 94.0%.

The tree was marked genuine (`mfg_mark_tree_real`), `PPM_bac.tre` was built by
the authors from these ASVs and its tips match the table exactly, which makes
UniFrac available.

---

## 3. Results

### 3.1 Parasite roots are less diverse: strongly confirmed

| Index | Parasite root | Host root | Wilcoxon |
|---|---|---|---|
| Observed ASVs | **392.5** | 734.0 | **p = 1.9e-09** |
| Shannon | **3.26** | 5.10 | **p = 1.9e-08** |

The paper's central claim. Median richness by group shows the pattern across the
whole design:

```
I-IR 721.0   I-L  56.5   I-R 762.0   I-S 1505.5
P-L   62.0   P-R 392.5
U-L   40.0   U-R 783.5   U-S 1263.0
```

Soil is richest, leaves poorest, roots intermediate, and the parasite root sits
roughly half way between host roots and leaves, which is itself the shape of the
paper's homogenisation argument.

### 3.2 The distance metric decides whether the claim is licensed

Parasite root versus host root, same samples, same test:

| Distance | R² | p | betadisper | Composition claim |
|---|---|---|---|---|
| Bray-Curtis | 0.2585 | 0.001 | **p = 0.018, heterogeneous** | **not licensed** |
| **UniFrac** | 0.1060 | 0.001 | p = 0.663, homogeneous | **licensed** |

This is the most methodologically interesting result in the series so far. The
two metrics agree that the groups differ. They disagree about whether that
difference can be read as a shift in composition, because dispersion is
heterogeneous under one and homogeneous under the other.

Bray-Curtis responds to which exact ASVs are present; parasite roots are far less
diverse, so their community occupies a smaller, more variable region of that
space. UniFrac collapses the ASVs onto shared phylogeny, and on that scale the
two groups are equally dispersed.

**The choice of distance metric determines not only the effect size but whether
the conclusion is permitted at all.** Neither the paper nor the tool's own
reference documentation makes that point, and it belongs in `reference/05`.

The claim is scored against UniFrac, where dispersion is homogeneous.

### 3.3 The parasite root resembles its own leaf, not the host root

Mean Bray-Curtis from parasite root:

| To | Distance |
|---|---|
| **its own leaf (P-L)** | **0.613** |
| infected host root (I-R) | 0.918 |
| uninfected host root (U-R) | 0.926 |
| soil (I-S) | 0.947 |

Two of the paper's claims fall directly out of this. The parasite's root and leaf
communities are far more alike than root communities are to each other, the
"increased homogenization between shoot and root tissues" claim, supported with a
0.3 gap. And the parasite is closer to host root than to soil, supporting
"congruency with *Hedera* root bacteria ... but not the surrounding soil".

Both are recorded as `unscored`, because no claim type expresses *X resembles Y
more than Z*, the same schema gap A3 raised.

---

## 4. Verdict

```
=== Benchmark verdict: INCOMPLETE ===

Claims that held:
  + parasite-less-diverse : lower, p = 1.9e-09
  + parasite-distinct     : UniFrac p = 0.001, dispersion homogeneous

Claims NOT adjudicated:
  ? shoot-root-homogenised  : supported by the distance matrix; no claim type
  ? congruent-root-not-soil : supported by the distance matrix; no claim type
  ? fewer-coassociations    : network property, not computed
```

Nothing failed. Both scorable claims held, and two of the three unadjudicated
ones are supported by evidence the schema cannot yet express.

---

## 5. Defects this run exposed

### #14: the clustering-variable candidate list was too narrow

`mfg_detect_repeated_measures()` looked only for names matching
`patient|subject|indiv|host|animal|mouse|id$`. A1 clusters by **site**: four
sites, 24–25 samples each, and `site` was never considered, so the report said
nothing about it at all.

Microbiome studies cluster by many names: site, plot, block, colony, nest, cage,
tank, litter, family, pair in ecology; batch, run, plate in the lab. The list now
covers these. The shape test from defect #13 still rejects `site` here, four
levels holding 25 samples each is a blocking factor, not a subject, but the
rejection is now **reported** rather than absent:

> Considered and rejected as subject identifiers: site (4 levels, up to 25
> samples each, looks like a grouping factor).

Silence and "checked, and here is why not" are different things, and only the
second is usable.

### #15: assigning metadata destroys the tracked normalization

Subsetting a TSS-transformed object and adding a derived column produced:

> This object's normalization is not tracked (inferred: unknown)

on an object normalized two lines earlier. Every distance computation was
refused.

Traced precisely:

| Step | Tracked state |
|---|---|
| `tss_transform(ps)` | `tss` |
| `prune_samples(...)` | `tss` (wrapper carries it) |
| **`sample_data(sub) <- sd`** | **`unknown`** |

`mfg_set_normalization()` already documents this hazard and writes a registry
backstop keyed on a content fingerprint, but subsetting *changes the content*,
so the fingerprint no longer matches and the backstop cannot fire. The guard
against the failure was itself defeated by the subset.

Added `mfg_set_meta()` and `mfg_add_meta()`, which perform the assignment and
carry the attributes, plus `mfg_prune_samples()` / `mfg_prune_taxa()` for the
subsetting step.

This one is worth noting for a reason beyond the fix: the failure is *loud and
wrong-looking*. It reports that normalization is untracked, which invites the
user to assert a state, and asserting the wrong one would silently produce an
invalid analysis. A guard that fails toward a plausible user error is more
dangerous than one that fails toward a stop.

---

## 6. Limitations

1. **Site clustering not modeled.** Four sites, ~25 samples each. Permuting
   within site would be the correct treatment and would likely reduce R².
2. **The paper's own analysis code ships with the deposit** (`PPM_Analysis_Code_Dryad.R`)
   and was not used. Comparing this reconstruction against it would separate
   "different choices" from "different implementation" and is the obvious next step.
3. **Network claims not tested.** "Fewer co-associations" needs co-occurrence
   network construction, which was not performed.
4. **Taxonomy stops at Genus.** RDP assignment gives six ranks; no species level.
5. **No sweep.** Single reconstruction at one prevalence threshold.

---

## 7. Reproducing this

```r
sp <- readRDS("A1/spe_table.rds")      # 103 x 21865, ASV sequences as colnames
tx <- readRDS("A1/RDP_taxa.rds")       # 21865 x 6, Kingdom..Genus
tr <- ape::read.tree("A1/PPM_bac.tre") # tips match the ASVs exactly

md <- do.call(rbind, strsplit(rownames(sp), "-", fixed = TRUE))
md <- data.frame(site = md[,1], patch = md[,2], tissue = md[,3], rep = md[,4],
                 row.names = rownames(sp), stringsAsFactors = FALSE)

ps <- phyloseq(otu_table(sp, taxa_are_rows = FALSE),
               tax_table(as.matrix(tx)), sample_data(md))
ps <- merge_phyloseq(ps, phy_tree(tr))
ps <- mfg_prune_samples(md$patch %in% c("I","P","U"), ps)   # drops 4 controls
ps <- mfg_mark_tree_real(ps, TRUE)
ps <- mfg_set_normalization(filter_prevalence(prune_low_depth(ps, 1000), 0.02), "raw")

pst   <- tss_transform(ps)
roots <- mfg_prune_samples(mfg_meta(pst)$tissue %in% c("R","IR"), pst)
roots <- mfg_add_meta(roots, "kind",
           ifelse(mfg_meta(roots)$patch == "P", "parasite", "host"))

for (dm in c("bray", "unifrac")) {
  d <- compute_distance(roots, dm)
  print(run_permanova(roots, ~ kind, dist_obj = d))
  print(check_dispersion(d, factor(mfg_meta(roots)$kind)))
}
```

Seeds fixed at 42. Every value above is written to `run_log.csv`.
