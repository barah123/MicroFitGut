# Benchmark: Reed et al. (2019), microbiota along the squirrel gut

Study **A2** of the validation matrix. It tests the QIIME 2 BIOM intake path, and
it is the study where two of the paper's own claims turn out to constrain each
other.

| | |
|---|---|
| **Paper** | Reed A, Pigage JC, Pigage HK, Glickman C, Bono JM. Comparative analysis of microbiota along the length of the gastrointestinal tract of two tree squirrel species (*Sciurus aberti* and *S. niger*) living in sympatry. *Ecol Evol.* 2019;9(23):13344–13358. |
| **Data** | Dryad [10.5061/dryad.931zcrjfn](https://doi.org/10.5061/dryad.931zcrjfn), one HDF5 BIOM and a FASTA of representative sequences |
| **Design** | Two squirrel species, 4 animals each, up to 8 gut regions per animal, 57 samples |
| **Verdict** | **Partially reproduced.** Two claims held, one failed, three unadjudicated |

---

## 1. What the deposit contains, and does not

The BIOM is structurally sound: 2,455 features by 57 samples, 11,969 non-zero
entries, written by QIIME 2 in December 2017. It carries **no taxonomy and no
sample metadata**. Both HDF5 metadata groups exist and are empty, the same shape
as A3's missing sample data.

The sample identifiers carry the design: `A1C`, `F3SI2` and so on encode species
(A for Abert's, F for fox), animal number, and gut region. That is enough to
rebuild the grouping without an external lookup, which is a better position than
A3, where the key had to come from NCBI.

Taxonomy is a different matter. The filename says `silva`, so taxonomy was
assigned at some point, but it is not in the deposit. Representative sequences
are, so taxonomy could be reassigned, but that is upstream work MicroFitGut does
not perform. **Every taxon-level claim in this paper is therefore untestable from
what was deposited**, including the fiber-degradation finding that motivates the
paper's discussion.

This is the third distinct shape of the same problem in ten studies. A3 was
missing sample metadata, recoverable from another repository. A2 is missing
taxonomy, recoverable only by rerunning classification. S5 mixed two file formats
in one folder.

---

## 2. Results

Quality control kept all 57 samples. Depth ranges 14,886 to 153,760, sparsity
91.4%.

### 2.1 The upper gut is less diverse, decisively

| Index | upper GI | lower GI | Wilcoxon |
|---|---|---|---|
| Observed features | 85.0 | 272.5 | **6.2e-11** |
| Shannon | 1.8 | 4.1 | **6.2e-09** |

Upper means stomach and the three small-intestine sections; lower means caecum,
the two large-intestine sections, and feces. A three-fold difference in richness.
The claim holds.

### 2.2 Species differ in the lower gut but not the upper

| Tract | PERMANOVA R² | p | betadisper |
|---|---|---|---|
| upper | 0.0878 | 0.071 | p = 0.426, homogeneous |
| lower | 0.1901 | **0.001** | **p = 0.001, heterogeneous** |

The upper-gut claim, that there is "relatively little differentiation between
the species", **holds**: no significant difference, and dispersion is homogeneous,
so the null is interpretable rather than masked.

The lower-gut claim is significant but **not licensed**. Dispersion differs
significantly between the species, so PERMANOVA cannot separate a shift in
composition from a difference in spread.

### 2.3 The stability claim explains why

The paper also states that "the Abert's squirrel lower GI community was more
stable in composition." Tested directly:

| | mean within-species Bray-Curtis | distance to centroid |
|---|---|---|
| Abert's | **0.606** | **0.426** |
| Fox | 0.783 | 0.558 |

Wilcoxon on pairwise distances p = 2.4e-09; betadisper permutest p = 0.001. The
stability claim **holds, and holds strongly**.

**The two claims constrain each other.** Heterogeneous dispersion is not an
artifact here, it is the paper's own finding stated in the language of a
different test. Abert's squirrels being more uniform is exactly what makes
dispersion unequal, and unequal dispersion is exactly what prevents the
significant PERMANOVA from being read as a pure composition shift.

The paper reports both claims without noting the interaction. Neither is wrong.
What the reanalysis establishes is that the second one qualifies the first, so
"clearly delineated" needs restating as something like: the species differ in
their lower-gut communities, and part of that difference is that one species
varies more than the other.

### 2.4 The diversity claim fails

The paper states that "overall microbial diversity was higher in fox squirrels."

| Index | Abert's | Fox | p |
|---|---|---|---|
| Observed features | 218.0 | 244.5 | 0.148 |
| Shannon | 3.5 | 3.4 | 0.806 |

Richness is numerically higher in fox squirrels, matching the direction, but not
significantly. Shannon shows no difference at all, and if anything runs the other
way. Scored as failed: the claim as written is a difference, and there is none
here.

---

## 3. Verdict

```
Claims that held:
  + upper-less-diverse     : upper 85.0 vs lower 272.5 features, p = 6.2e-11
  + upper-no-species-diff  : p = 0.071, dispersion homogeneous

Claims that did not hold:
  - fox-higher-diversity   : p = 0.148 richness, p = 0.806 Shannon

Claims NOT adjudicated:
  ? lower-species-diff     : significant at p = 0.001 but dispersion
                             heterogeneous, so not licensed
  ? aberts-more-stable     : holds on direct test (0.426 vs 0.558,
                             permutest p = 0.001); no claim type encodes
                             a dispersion hypothesis
  ? fiber-degraders        : untestable, no taxonomy in the deposit
```

---

## 4. Defects this run exposed

### #17, an opaque failure on a BIOM without taxonomy

`phyloseq::import_biom()` failed with:

```
length of 'dimnames' [2] not equal to array extent
```

The message names neither the file nor the cause. The natural reading is a
corrupt download. The table is fine; it simply has no taxonomy, and
`import_biom()` assumes there is some.

Added `mfg_read_biom_hdf5()`, which reads the sparse matrix directly and attaches
taxonomy and sample data only when datasets actually exist under those groups. A
table without taxonomy now loads as a table without taxonomy, with a warning
naming the consequence:

> carries no taxonomy. Feature IDs are all that identify a row, so nothing can be
> agglomerated or matched by name.

### #18, a single-component phyloseq is not a phyloseq

`phyloseq::phyloseq()` given one component returns **that component**, not a
phyloseq object. So a BIOM with no taxonomy and no metadata came back as a bare
`otu_table`, and the next call failed with "no slot of name otu_table for this
object of class otu_table".

This one was in code written during this study, which is a fair illustration of
how the trap works: the object looks right, prints plausibly, and fails two steps
later with a message pointing at the wrong place. `mfg_set_meta()` now constructs
rather than assigns when it receives a bare component.

### #19, stratified permutation cannot test a between-subject term

The first PERMANOVA run returned **p = 1.000 in both tracts**, for a difference
the lower-gut data clearly contains.

The cause is that species is constant within animal: `A1` is always an Abert's
squirrel. Permuting within animal therefore reproduces the observed species
assignment in every permutation, so the p-value is 1 by construction.

This is a serious silent failure. It produces the most convincing possible null
result, p exactly 1.000, from a design error rather than from the data. An analyst
following the tool's own earlier advice, which was to pass `strata` whenever
repeated measures are detected, would land in it.

`run_permanova()` now refuses:

> Cannot test species with permutations restricted to 'animal': species does not
> vary within a single animal, so every permutation reproduces the observed
> assignment and the p-value is 1 by construction. This is a between-animal
> comparison, not a within-animal one.

The repeated-measures warning was also made conditional, since it had been
recommending exactly what the new guard rejects. When the tested term does not
vary within subject it now says so and gives the effective sample size:

> the tested term does not vary within animal, so this is a between-animal
> comparison and strata cannot help. Treat each animal as one observation: the
> effective n is 8, not 32.

That second number matters. With 8 animals rather than 32 samples, the lower-gut
result rests on four animals per species.

---

## 5. A boundary case in the subject detector

`animal` has exactly 8 levels with up to 8 samples each. The shape test from
defect #13 rejects a candidate when one level holds at least as many samples as
there are levels, so this design sits exactly on the boundary and was rejected by
automatic detection. Passing `subject_var = "animal"` explicitly worked and
produced the warning added for that case.

The heuristic cannot resolve a tie of this shape, and the explicit override is
the correct escape. Worth recording because 8 by 8 designs are not rare in animal
studies.

---

## 6. Limitations

1. **No taxonomy.** Every taxon-level claim is untestable from the deposit.
2. **Effective n is 8, not 57.** Four animals per species. The lower-gut
   difference rests on that.
3. **Regions treated as independent within animal.** A proper analysis would
   model animal as a random effect or test region as a within-animal term.
4. **Upper and lower boundaries are my assignment**, not the paper's: stomach and
   small intestine against caecum, large intestine and feces. The paper may draw
   the line differently.
5. **No sweep.**

---

## 7. Reproducing this

```r
# phyloseq::import_biom() fails on this file; the taxonomy group is empty
ps <- mfg_read_biom_hdf5("A2/table_no_mitochondria_no_chloroplast_silva.biom")

sn <- sample_names(ps)
md <- data.frame(
  species = ifelse(substr(sn, 1, 1) == "A", "Aberts", "Fox"),
  animal  = substr(sn, 1, 2),
  region  = sub("^[AF][0-9]", "", sn),
  row.names = sn, stringsAsFactors = FALSE)
md$tract <- ifelse(md$region %in% c("S","SI1","SI2","SI3"), "upper", "lower")
ps <- mfg_set_normalization(mfg_set_meta(ps, md), "raw")

# strata = "animal" is refused here: species does not vary within animal
pst <- tss_transform(ps)
lower <- mfg_prune_samples(sn[md$tract == "lower"], pst)
d <- compute_distance(lower, "bray")
run_permanova(lower, ~ species, dist_obj = d, check_disp = FALSE)
check_dispersion(d, factor(mfg_meta(lower)$species))
```

Seeds fixed at 42. Every value above is written to `run_log.csv`.
