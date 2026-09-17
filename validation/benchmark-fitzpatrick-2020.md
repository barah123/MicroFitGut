# Benchmark: Fitzpatrick & Schneider (2020), a parasitic plant and its host

A full benchmark-mode run against a published study, in the same form as the
Pérez-Losada benchmark: reconstruct the analysis, score the paper's stated
claims, and report which survive. Study **A1** of the validation matrix,
intended as the negative control, and the first study to exercise the
phylogenetic distance path on a real tree.

| | |
|---|---|
| **Paper** | Fitzpatrick CR, Schneider AC. Unique bacterial assembly, composition, and interactions in a parasitic plant and its host. *J Exp Bot.* 2020;71(6):2198-2209. [10.1093/jxb/erz572](https://doi.org/10.1093/jxb/erz572) |
| **Raw data** | Dryad [10.5061/dryad.7wm37pvnk](https://doi.org/10.5061/dryad.7wm37pvnk): ASV table, RDP taxonomy, **phylogeny**, and the authors' analysis code |
| **Analyzed** | 21,865 ASVs × 103 samples, filtered to 9,488 taxa × 99 samples |
| **Design** | *Orobanche hederae* holoparasite and its *Hedera* host, 4 sites, leaf/root/soil, infected and uninfected patches |
| **Hypothesis** | A holoparasitic plant assembles a bacterial community distinct from its host's, less diverse and more homogeneous between shoot and root |
| **Verdict** | **Incomplete.** 2 of 2 scorable claims reproduced, 3 outside the denominator. No sweep |

---

## 1. Study card

**8 fields stated · 3 inferred · 12 absent.**

The paper deposits its own analysis code, which is unusual and valuable, and
which makes several nominally absent fields recoverable. They were not recovered
for this run; see the limitations.

### The matrix assignment was wrong, and that is useful

A1 was placed in the matrix as a **negative control**: a clean cross-sectional
design where no guard should fire. The validation plan flagged that designs were
inferred from titles and file listings, and warned they might not hold.

They did not. A1 is a **paired infected and uninfected design across four sites
with multiple tissue types per plant**, nine groups in total:

| | leaf | root | infected root | soil |
|---|---|---|---|---|
| **Parasite (P)** | 12 | 12 | | |
| **Infected host (I)** | 12 | 12 | 12 | 8 |
| **Uninfected host (U)** | 12 | 12 | | 7 |

Plus four technical controls (`mock`, `PAO1`, `water`, `Undetermined`), dropped
and logged.

So A1 is not a negative control. **The validation set still lacks one**, and that
matters: without a clean study where nothing fires, there is no evidence
distinguishing "the guards catch real problems" from "the guards fire on
everything". This is recorded as a gap in the set, and it is the reason the set
cannot support any claim about guard specificity.

What A1 *is* good for turns out to be more valuable: it has a real phylogeny
covering all 21,865 ASVs, which no earlier study in the series did.

### The paper's five claims

| id | Claim | Evidence in the paper |
|---|---|---|
| `parasite-less-diverse` | Parasite roots are less diverse than host roots | stated as reduced bacterial diversity in the parasite |
| `parasite-distinct` | Parasite and host root communities differ | stated as unique composition |
| `shoot-root-homogenised` | Parasite shoot and root communities are more alike | "increased homogenization between shoot and root tissues" |
| `congruent-root-not-soil` | The parasite resembles host root, not soil | "congruency with *Hedera* root bacteria ... but not the surrounding soil" |
| `fewer-coassociations` | The parasite has fewer bacterial co-associations | network result |

---

## 2. Reconstruction

Minimum depth 1,000, prevalence 0.02: **21,865 taxa reduced to 9,488**, 99
samples retained. Depth range 28,701 to 261,684, a 9.1-fold spread. Sparsity
94.0%.

The tree was marked genuine with `mfg_mark_tree_real()`. `PPM_bac.tre` was built
by the authors from these ASVs and its tips match the table exactly, which is
what makes UniFrac available here and unavailable in most of the series.

---

## 3. Results

### 3.1 Parasite roots are less diverse

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
roughly halfway between host roots and leaves, which is itself the shape of the
paper's homogenization argument.

### 3.2 The distance metric decides whether the claim is licensed

Parasite root against host root, same samples, same test:

| Distance | R² | p | betadisper | Composition claim |
|---|---|---|---|---|
| Bray-Curtis | 0.2585 | 0.001 | **p = 0.018, heterogeneous** | **not licensed** |
| **UniFrac** | 0.1060 | 0.001 | p = 0.663, homogeneous | **licensed** |

This is the most methodologically interesting result in the series. The two
metrics agree that the groups differ. They disagree about whether that difference
can be read as a shift in composition, because dispersion is heterogeneous under
one and homogeneous under the other.

Bray-Curtis responds to which exact ASVs are present. Parasite roots are far less
diverse, so their community occupies a smaller and more variable region of that
space. UniFrac collapses the ASVs onto shared phylogeny, and on that scale the
two groups are equally dispersed.

**The choice of distance metric determines not only the effect size but whether
the conclusion is permitted at all.** Neither the paper nor the tool's own
reference documentation made that point before this run; it has since been added
to `reference/05`.

The claim is scored against UniFrac, where dispersion is homogeneous. Scoring it
against Bray-Curtis would have returned `not_licensed`, and the choice between
the two was made before the dispersion tests were seen, not after.

### 3.3 The parasite root resembles its own leaf, not the host root

Mean Bray-Curtis from parasite root:

| To | Distance |
|---|---|
| **its own leaf (P-L)** | **0.613** |
| infected host root (I-R) | 0.918 |
| uninfected host root (U-R) | 0.926 |
| soil (I-S) | 0.947 |

Two of the paper's claims fall directly out of this. The parasite's root and leaf
communities are far more alike than root communities are to each other, which is
the homogenization claim, supported with a 0.3 gap. And the parasite is closer to
host root than to soil, supporting the congruency claim.

Both are recorded as `schema_gap`, because no claim type expresses *X resembles Y
more than Z*. This is the same gap A3 raised, and it is now the most common
unscorable reason in the series.

---

## 4. Verdict

```
=== Benchmark verdict: INCOMPLETE ===

Reproduced (2):
  + parasite-less-diverse : lower in parasite, p = 1.9e-09
  + parasite-distinct     : UniFrac p = 0.001, dispersion homogeneous

Not reproduced (0)
Not licensed  (0)

Outside the denominator (3):
  ? shoot-root-homogenised  : schema_gap    - supported by the distance
                              matrix (0.613 vs 0.918); no claim type
  ? congruent-root-not-soil : schema_gap    - supported by the distance
                              matrix; no claim type
  ? fewer-coassociations    : not_attempted - network not computed

Scorable: 2 of 5.  Reproduced: 2 of 2.
```

Nothing was contradicted. Both scorable claims reproduced, and two of the three
outside the denominator are supported by evidence the schema cannot yet express.

### What this means

Five of the paper's claims were extracted and two could be scored. That ratio is
the study's most useful result, and it is about the instrument rather than the
paper: **three claims went unscored and only one of those was a data problem.**
Two are cases where the evidence exists, points the paper's way, and has no slot
to go in.

That is a concrete specification for the next version of the claim schema. A
`similarity` claim type taking three groups and asserting that the first is
closer to the second than to the third would have scored both, using a distance
matrix the reanalysis already computed.

### Limitations

1. **Site clustering not modeled.** Four sites with roughly 25 samples each.
   Permuting within site would be the correct treatment and would likely reduce
   R².
2. **The paper's own analysis code ships with the deposit**
   (`PPM_Analysis_Code_Dryad.R`) and was not used. Comparing this reconstruction
   against it would separate "different choices" from "different implementation",
   and it is the obvious next step for this study specifically, since almost no
   other deposit in the set makes it possible.
3. **Network claims not tested.** Fewer co-associations needs co-occurrence
   network construction, which was not performed.
4. **Taxonomy stops at Genus.** RDP assignment gives six ranks, no species level.
5. **No sweep.** Single reconstruction at one prevalence threshold. Given that
   section 3.2 shows the distance metric changing whether a claim is licensed, a
   sweep over the metric would be more informative here than anywhere else in the
   series.

---

## 5. Defects this run exposed

### #14: the clustering-variable candidate list was too narrow

`mfg_detect_repeated_measures()` looked only for names matching
`patient|subject|indiv|host|animal|mouse|id$`. A1 clusters by **site**: four
sites with 24 to 25 samples each, and `site` was never considered, so the report
said nothing about it at all.

Microbiome studies cluster by many names: site, plot, block, colony, nest, cage,
tank, litter, family and pair in ecology; batch, run and plate in the lab. The
list now covers these. The shape test from defect #13 still rejects `site` here,
because four levels holding 25 samples each is a blocking factor rather than a
subject, but the rejection is now **reported** rather than absent:

> Considered and rejected as subject identifiers: site (4 levels, up to 25
> samples each, looks like a grouping factor).

Silence and "checked, and here is why not" are different things, and only the
second is usable.

**Direction of bias.** This defect is **anti-conservative**. An undetected
cluster means the naive unstratified test runs, and significance is inflated,
which makes a claim of difference *more* likely to be scored as reproduced. It is
one of two defects in the series that run this way, the other being #16 in S5.
The series-level statement that every defect biased toward non-reproduction does
not hold for it.

### #15: assigning metadata destroys the tracked normalization

Subsetting a TSS-transformed object and adding a derived column produced:

> This object's normalization is not tracked (inferred: unknown)

on an object normalized two lines earlier. Every distance computation was
refused.

Traced precisely:

| Step | Tracked state |
|---|---|
| `tss_transform(ps)` | `tss` |
| `prune_samples(...)` | `tss`, the wrapper carries it |
| **`sample_data(sub) <- sd`** | **`unknown`** |

`mfg_set_normalization()` already documents this hazard and writes a registry
backstop keyed on a content fingerprint, but subsetting *changes the content*, so
the fingerprint no longer matches and the backstop cannot fire. The guard against
the failure was itself defeated by the subset.

Added `mfg_set_meta()` and `mfg_add_meta()`, which perform the assignment and
carry the attributes, plus `mfg_prune_samples()` and `mfg_prune_taxa()` for the
subsetting step.

This one is worth noting for a reason beyond the fix: the failure is *loud and
wrong-looking*. It reports that normalization is untracked, which invites the
user to assert a state, and asserting the wrong one would silently produce an
invalid analysis. A guard that fails toward a plausible user error is more
dangerous than one that fails toward a stop.

---

## 6. Reproducing this

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

# Both metrics, because the dispersion check disagrees between them.
for (dm in c("bray", "unifrac")) {
  d <- compute_distance(roots, dm)
  print(run_permanova(roots, ~ kind, dist_obj = d))
  print(check_dispersion(d, factor(mfg_meta(roots)$kind)))
}
```

Seeds fixed at 42. Every value above is written to `run_log.csv`.
