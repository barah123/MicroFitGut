# Benchmark: Reed et al. (2019), microbiota along the squirrel gut

A full benchmark-mode run against a published study, in the same form as the
Pérez-Losada benchmark: reconstruct the analysis, score the paper's stated
claims, and report which survive. Study **A2** of the validation matrix. It
tests the QIIME 2 BIOM intake path, and it is the study where two of the paper's
own claims turn out to constrain each other.

| | |
|---|---|
| **Paper** | Reed A, Pigage JC, Pigage HK, Glickman C, Bono JM. Comparative analysis of microbiota along the length of the gastrointestinal tract of two tree squirrel species (*Sciurus aberti* and *S. niger*) living in sympatry. *Ecol Evol.* 2019;9(23):13344-13358. [10.1002/ece3.5789](https://doi.org/10.1002/ece3.5789) |
| **Raw data** | Dryad [10.5061/dryad.931zcrjfn](https://doi.org/10.5061/dryad.931zcrjfn): one HDF5 BIOM and a FASTA of representative sequences |
| **Analyzed** | 2,455 features × 57 samples. No taxonomy and no sample metadata in the deposit |
| **Design** | Two squirrel species, 4 animals each, up to 8 gut regions per animal. Effective n is 8 animals, not 57 samples |
| **Hypothesis** | Gut community composition varies along the gastrointestinal tract, and the two sympatric squirrel species differ in their gut microbiota in ways reflecting diet |
| **Verdict** | **Partially reproduced.** 2 of 4 scorable claims reproduced, 1 not reproduced, 1 not licensed, 2 outside the denominator. No sweep |

---

## 1. Study card

**7 fields stated · 2 inferred · 14 absent.**

### What the deposit contains, and does not

The BIOM is structurally sound: 2,455 features by 57 samples, 11,969 non-zero
entries, written by QIIME 2 in December 2017. It carries **no taxonomy and no
sample metadata**. Both HDF5 metadata groups exist and are empty, the same shape
as A3's missing sample data.

The sample identifiers carry the design. `A1C` and `F3SI2` encode species, A for
Abert's and F for fox, then animal number, then gut region. That is enough to
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

### The paper's six claims

| id | Claim | Evidence in the paper |
|---|---|---|
| `upper-less-diverse` | The upper GI tract is less diverse than the lower | stated as a diversity gradient along the tract |
| `upper-no-species-diff` | The species do not differ in the upper tract | "relatively little differentiation between the species" |
| `lower-species-diff` | The species differ in the lower tract | "the community composition of the lower GI tract was clearly delineated" |
| `aberts-more-stable` | Abert's lower GI community is more stable | "the Abert's squirrel lower GI community was more stable in composition" |
| `fox-higher-diversity` | Fox squirrels have higher overall diversity | "overall microbial diversity was higher in fox squirrels" |
| `fiber-degraders` | Abert's is enriched for fiber degraders | discussion of plant fiber degradation |

---

## 2. Reconstruction

Quality control kept all 57 samples. Depth ranges 14,886 to 153,760, sparsity
91.4%.

Upper tract means stomach and the three small-intestine sections. Lower means
caecum, the two large-intestine sections, and feces. **That boundary is my
assignment, not the paper's**, and it is a reconstruction choice that could
reasonably be drawn elsewhere.

`animal` has exactly 8 levels with up to 8 samples each, which sits exactly on
the boundary of the shape test from defect #13 and was rejected by automatic
detection. Passing `subject_var = "animal"` explicitly worked and produced the
warning added for that case. The heuristic cannot resolve a tie of this shape,
and the explicit override is the correct escape. Worth recording, because 8 by 8
designs are not rare in animal studies.

---

## 3. Results

### 3.1 The upper gut is less diverse, decisively

| Index | upper GI | lower GI | Wilcoxon |
|---|---|---|---|
| Observed features | 85.0 | 272.5 | **6.2e-11** |
| Shannon | 1.8 | 4.1 | **6.2e-09** |

A threefold difference in richness. The claim reproduces.

### 3.2 Species differ in the lower gut but not the upper

| Tract | PERMANOVA R² | p | betadisper |
|---|---|---|---|
| upper | 0.0878 | 0.071 | p = 0.426, homogeneous |
| lower | 0.1901 | **0.001** | **p = 0.001, heterogeneous** |

The upper-gut claim **reproduces**: no significant difference, and dispersion is
homogeneous, so the null is interpretable rather than a failure to detect a shift
masked by unequal spread.

The lower-gut claim is significant but **not licensed**. Dispersion differs
significantly between the species, so PERMANOVA cannot separate a shift in
composition from a difference in spread.

### 3.3 The stability claim explains why

The paper also states that the Abert's squirrel lower GI community was more
stable in composition. Tested directly:

| | mean within-species Bray-Curtis | distance to centroid |
|---|---|---|
| Abert's | **0.606** | **0.426** |
| Fox | 0.783 | 0.558 |

Wilcoxon on pairwise distances p = 2.4e-09; betadisper permutest p = 0.001. The
stability claim holds, and holds strongly, but it is recorded as `schema_gap`
because no claim type encodes a dispersion hypothesis.

**The two claims constrain each other.** Heterogeneous dispersion is not an
artifact here. It is the paper's own finding stated in the language of a
different test. Abert's squirrels being more uniform is exactly what makes
dispersion unequal, and unequal dispersion is exactly what prevents the
significant PERMANOVA from being read as a pure composition shift.

The paper reports both claims without noting the interaction. Neither is wrong.
What the reanalysis establishes is that the second one qualifies the first, so
"clearly delineated" needs restating as something like: the species differ in
their lower-gut communities, and part of that difference is that one species
varies more than the other.

### 3.4 The diversity claim does not reproduce

| Index | Abert's | Fox | p |
|---|---|---|---|
| Observed features | 218.0 | 244.5 | 0.148 |
| Shannon | 3.5 | 3.4 | 0.806 |

Richness is numerically higher in fox squirrels, matching the direction, but not
significantly. Shannon shows no difference at all and if anything runs the other
way. The claim as written asserts a difference, and there is none here.

---

## 4. Verdict

```
=== Benchmark verdict: PARTIALLY REPRODUCED ===

Reproduced (2):
  + upper-less-diverse     : upper 85.0 vs lower 272.5 features, p = 6.2e-11
  + upper-no-species-diff  : p = 0.071, dispersion homogeneous

Not reproduced (1):
  - fox-higher-diversity   : p = 0.148 richness, p = 0.806 Shannon

Not licensed (1):
  ! lower-species-diff     : significant at p = 0.001, but betadisper
                             p = 0.001, so PERMANOVA does not license a
                             composition claim

Outside the denominator (2):
  ? aberts-more-stable     : schema_gap  - holds on direct test
                             (0.426 vs 0.558, permutest p = 0.001); no
                             claim type encodes a dispersion hypothesis
  ? fiber-degraders        : data_absent - no taxonomy in the deposit

Scorable: 4 of 6.  Reproduced: 2 of 4.
```

**The `not_licensed` status is new since this report was first written.** The
lower-gut claim was previously pooled with the unadjudicated claims, which put a
finding about the paper in the same bucket as a missing metadata column. It now
carries its own status and sits inside the denominator, because the test ran and
returned a result; what is missing is the licence to read that result as a
composition claim.

### What this means

Four of the paper's six claims could be scored, and two reproduced. The two that
did not fail in quite different ways, and the difference is the point of this
study.

`fox-higher-diversity` is a failure to detect. The direction agrees on richness
and reverses on Shannon, neither significantly. It is not a contradiction.

`lower-species-diff` is not a failure at all. The effect is there, at p = 0.001
with R² = 0.19. What the reanalysis withholds is the interpretation, and it
withholds it on the strength of the paper's *own* other claim. This is the
cleanest example in the series of a paper containing the qualification to its own
result without connecting the two.

### Limitations

1. **No taxonomy.** Every taxon-level claim is untestable from the deposit.
2. **Effective n is 8, not 57.** Four animals per species. The lower-gut
   difference rests on that, and no interval reported here reflects it.
3. **Regions treated as independent within animal.** A proper analysis would
   model animal as a random effect or test region as a within-animal term.
4. **The upper and lower boundary is my assignment**, not the paper's. The paper
   may draw the line differently, and the diversity contrast in 3.1 depends on
   where it is drawn.
5. **No sweep.** The boundary choice in point 4 is the obvious axis to vary.

---

## 5. Defects this run exposed

### #17: an opaque failure on a BIOM without taxonomy

`phyloseq::import_biom()` failed with:

```
length of 'dimnames' [2] not equal to array extent
```

The message names neither the file nor the cause. The natural reading is a
corrupt download. The table is fine. It simply has no taxonomy, and
`import_biom()` assumes there is some.

Added `mfg_read_biom_hdf5()`, which reads the sparse matrix directly and attaches
taxonomy and sample data only when datasets actually exist under those groups. A
table without taxonomy now loads as a table without taxonomy, with a warning
naming the consequence:

> carries no taxonomy. Feature IDs are all that identify a row, so nothing can be
> agglomerated or matched by name.

### #18: a single-component phyloseq is not a phyloseq

`phyloseq::phyloseq()` given one component returns **that component**, not a
phyloseq object. So a BIOM with no taxonomy and no metadata came back as a bare
`otu_table`, and the next call failed with "no slot of name otu_table for this
object of class otu_table".

This one was in code written during this study, which is a fair illustration of
how the trap works: the object looks right, prints plausibly, and fails two steps
later with a message pointing at the wrong place. `mfg_set_meta()` now constructs
rather than assigns when it receives a bare component.

### #19: stratified permutation cannot test a between-subject term

The first PERMANOVA run returned **p = 1.000 in both tracts**, for a difference
the lower-gut data clearly contains.

The cause is that species is constant within animal. `A1` is always an Abert's
squirrel, so permuting within animal reproduces the observed species assignment
in every permutation, and the p-value is 1 by construction.

This is a serious silent failure. It produces the most convincing possible null
result, p exactly 1.000, from a design error rather than from the data. An
analyst following the tool's own earlier advice, which was to pass `strata`
whenever repeated measures are detected, would land in it.

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

This is the one demonstrable correctness result in the series: a case where the
tool caught a p-value produced by construction rather than by the data. It is
worth distinguishing from the general claim that the guards work, which the
validation set cannot establish, since it contains no clean negative control.

---

## 6. Reproducing this

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
pst   <- tss_transform(ps)
lower <- mfg_prune_samples(sn[md$tract == "lower"], pst)
d <- compute_distance(lower, "bray")
run_permanova(lower, ~ species, dist_obj = d, check_disp = FALSE)
check_dispersion(d, factor(mfg_meta(lower)$species))
```

Seeds fixed at 42. Every value above is written to `run_log.csv`.
