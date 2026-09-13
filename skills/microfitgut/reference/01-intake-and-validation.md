# Intake and validation

Every check here exists because skipping it produces an analysis that runs
cleanly and answers the wrong question.

## Accepted inputs

| Input | How | Notes |
|---|---|---|
| Saved phyloseq object | `mfg_load("x.RDS")` | Fastest path. Normalization state is **untracked** on arrival — assert it with `mfg_set_normalization()` before any analysis, or the pairing guard will stop you. |
| ASV/OTU + taxonomy + metadata triplet | `build_phyloseq(abund, tax, meta, tree, seqs)` | Delimiter is detected per file; `.csv` and tab-delimited `.txt` mix freely in practice. |
| Shotgun taxonomic profile | `build_phyloseq_from_profile(path, meta)` | MetaPhlAn, Kraken2/Bracken. Handles the version banner line and collapses the repeated-clade hierarchy. |
| Functional profile | `build_phyloseq_functional(abund, pathway_tax, meta)` | PICRUSt2, HUMAnN. Pathways occupy the taxonomy slot, giving a `Pathways` rank. |
| BIOM table | `mfg_load("x.biom")` | Via `phyloseq::import_biom`. |

## The mandatory checks

`validate_inputs(ps, group_var, subject_var)` runs all of these and **stops** on
anything in the "blocks" column. Nothing else runs until it passes.

### 1. Sample IDs match across all three tables — blocks

phyloseq silently intersects on construction. A mismatch does not error; it
produces an object with fewer samples than you think you have. The validator
reports the counts on each side.

If they do not match, the cause is almost always one of: row names read as a
data column (`row.names = 1` omitted), a trailing-whitespace difference, or
`_S1_L001`-style suffixes present in one file and not the other.

### 2. Taxa IDs match — blocks

Same mechanism, same silence.

### 3. Table orientation — blocks if ambiguous

Getting taxa-are-rows backwards is silent and catastrophic: every per-sample
statistic becomes a per-taxon statistic and nothing errors.

`mfg_detect_orientation()` decides by **matching row and column names against the
taxonomy and metadata**, not by comparing the two dimensions. Shape is not
evidence: a study can easily have more samples than taxa.

When the evidence is tied, it refuses and asks. Pass `taxa_are_rows` explicitly.

### 4. Counts or proportions — changes what is legal downstream

`looks_like_relative_abundance()` decides. Every column summing to exactly 100 or
1 settles it. Otherwise the test is fractional values: read counts are whole
numbers. MetaPhlAn leaves an unclassified fraction out, so a real profile often
sums to the nineties, which is why totals alone are not decisive.

Consequences of proportions:

- Chao1, ACE and Fisher are **unavailable** — they are estimated from taxa seen
  exactly once or twice, and there are no counts.
- Rarefaction is **impossible** — there are no reads to subsample.
- Singleton removal and depth pruning are **meaningless** — sample sums carry no
  depth information.

Negative values mean the table has already been CLR- or VST-transformed. Every
count-based method is then invalid on it.

### 5. Read depth distribution — warns

Reported as min / median / max and the fold range. A fold range above 10 means
normalization is not optional. Samples below 1,000 reads are listed individually
so the drop decision is made visibly rather than by a default.

Zero-depth samples **block**: they break `log()` offsets, distances and every
per-sample proportion.

### 6. Sparsity — warns

Zero proportion across the whole table. Above 80%, zero-inflated models
(`07`) become relevant and the structural-versus-sampling-zero distinction in
`11` starts to matter. Taxa that are zero in every sample are counted; they
contribute nothing but multiple-testing burden.

### 7. Tree provenance — warns, and gates phylogenetic metrics

**This is the check the course material does not have, and it matters.**

Demo 3 and Demo 6 both attach `rtree()` — a random tree — to satisfy phyloseq's
slot for teaching purposes. UniFrac and Faith's PD computed on a random topology
return a matrix and a vector that look entirely valid and encode nothing about
evolutionary relationship.

So a tree's provenance must be recorded:

```r
ps <- mfg_mark_tree_real(ps, TRUE)   # estimated from the sequence data
ps <- mfg_mark_tree_real(ps, FALSE)  # placeholder
```

Until it is, `faith_pd()`, `compute_distance(method = "unifrac")` and
`plot_phylogeny()` all refuse. Unmarked is treated as unknown, not as real —
absent provenance we cannot prove a tree is genuine, so the caller is asked.

Also checked: is the tree rooted? Faith's PD with `include_root = TRUE` and
unweighted UniFrac both depend on rooting. Root with `phangorn::midpoint()`.

### 8. Metadata variable types — warns

Character columns that look categorical are **not factors**, so their reference
level is alphabetical rather than chosen. Every model coefficient and every
DESeq2/ANCOM-BC2 log fold change is relative to that level. Set it deliberately:

```r
sample_data(ps)$status <- factor(sample_data(ps)$status,
                                 levels = c("healthy", "disease"))
```

Getting this wrong does not error; it flips the sign of every reported effect.

### 9. Group sizes and balance — warns

n per level, plus the imbalance ratio. Groups under about 3 support no meaningful
test. An imbalance beyond 3:1 matters for PERMANOVA and ANOSIM specifically — see
`11`.

### 10. Repeated measures — warns loudly

`mfg_detect_repeated_measures()` looks for a variable with several samples per
level and fewer levels than rows. This is the check that changes the whole
analysis plan, because independent-samples tests are then invalid.

On the bundled skin dataset it fires: `patient` has up to 5 samples across 128
participants. The course demos run Kruskal-Wallis and unstratified PERMANOVA on
exactly this data, which treats correlated samples as independent draws and
overstates significance. `04` and `05` will refuse or warn; `07` has the right
models.

## Reading the validation object

```r
v <- validate_inputs(ps, group_var = "region", subject_var = "patient")
print(v)
```

`v$problems` blocks. `v$warnings` must be **stated in the report** —
`10-reporting-standards.md` requires it, and `assemble_summary()` carries them
into the output automatically.

Run with `strict = FALSE` only to inspect a broken dataset, never to proceed.
