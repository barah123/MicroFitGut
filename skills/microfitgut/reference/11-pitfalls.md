# Pitfalls

The expert-knowledge file. Read once fully, then revisit the relevant entries per
stage. Each of these produces an analysis that runs cleanly and is wrong.

---

## 1. Compositionality

**The problem.** Sequencing measures *relative* abundance. Total microbial load is
not recoverable from the data. A table of proportions is constrained to sum to
one, so an increase in one taxon forces the others down arithmetically.

**What it breaks.**

- A taxon can appear to "decrease" when its absolute abundance is unchanged and
  something else rose.
- Correlations between proportions are spurious by construction — they carry a
  built-in negative bias. This is why `cooccurrence_network()` **refuses**
  relative-abundance input.
- Any test treating proportions as independent measurements per taxon is testing
  something other than what it claims.

**What to do.** CLR for descriptive and ordination work; ANCOM-BC2 or ALDEx2 for
inference, both of which model the problem rather than ignoring it. And state
plainly that conclusions are about relative abundance, because absolute abundance
was not measured.

**What CLR does not fix.** A CLR value is relative to the sample's geometric mean,
so a positive CLR coefficient still admits two readings: the taxon rose, or
everything else fell. Say which the design can distinguish — usually neither,
without spike-ins or qPCR total load.

---

## 2. The rarefaction dispute

**Both sides are right about different things, and the literature reads as though
they are not.**

McMurdie & Holmes (2014) showed rarefying is statistically inadmissible for
differential abundance: it discards data and adds variance, and a model-based
method using all the counts is strictly better.

That critique is correct **for differential abundance**, which is why every DA
method here takes raw counts. It does not extend to richness estimation, where
Chao1 and Observed are estimated from singletons and doubletons and equal
sampling effort is the only way the estimator means anything.

**What to do.** Rarefy for richness. Do not rarefy for DA. State which you did
and for what — the pairing guard in `03` enforces it, and a reviewer from either
camp will accept a stated, consistent choice more readily than an unstated one.

---

## 3. Repeated measures treated as independent

**The most consequential error in this literature, and the easiest to miss** —
because nothing errors, the analysis looks normal, and the results look stronger.

Several samples per subject are correlated. An independent-samples test treats
them as independent draws, underestimates the standard errors, and produces false
positives. Demo 12 shows a predictor that is "highly significant" under `lm` and
has no effect once `(1 | patient)` is added.

**Where it hides.** Longitudinal designs are obvious. Less obvious:

- multiple body sites per participant (the bundled skin dataset: 128 patients,
  up to 5 samples each)
- technical replicates pooled with biological ones
- litter, cage or tank effects in animal studies — cagemates share a microbiome,
  so the cage is the experimental unit, not the animal
- family or household members

**Detection.** `mfg_detect_repeated_measures()` flags a variable with several
samples per level and fewer levels than rows. It runs inside
`validate_inputs()`.

**What to do.**

| Analysis | Fix |
|---|---|
| alpha diversity | `fit_lmm(Shannon ~ group + (1 \| subject))`. `run_alpha_test()` refuses the simple test. |
| beta diversity | `run_permanova(..., strata = "subject")` — permute within subject. |
| differential abundance | ANCOM-BC2 with `rand_formula = "(1 \| subject)"`, or per-taxon GLMM. |
| any regression | a subject random effect (`07`). |

**The diagnostic if you suspect it after the fact.** A DESeq2 p-value histogram
with too few large p-values is anti-conservative, and unmodelled repeated measures
is the commonest cause. `da_deseq2_diagnostics()` computes and reads this.

---

## 4. The betadisper confound

**A significant PERMANOVA does not mean the groups differ in composition.**

`adonis2` responds to both location (centroid) and dispersion (spread). Groups
with identical centroids and different spread give a significant result.

On the bundled skin data: PERMANOVA R² = 0.30, p = 0.005 — and betadisper
p = 0.005, with mean distance-to-centroid 0.272 / 0.325 / 0.468 across the three
regions. The significant PERMANOVA **cannot** be reported as a composition shift.

**What to do.** Always run `check_dispersion()` — `run_permanova()` does it
automatically — and let `permanova_sentence()` write the conclusion, because it
refuses to say "composition differed" when dispersion is heterogeneous.

Heterogeneous dispersion is itself a finding. A group that is more variable than
another is a real biological statement (dysbiosis is often *more variable*, not
differently centred — the Anna Karenina principle).

ANOSIM is *more* dispersion-sensitive than PERMANOVA, so it does not resolve the
ambiguity; it deepens it.

---

## 5. Unequal group sizes

PERMANOVA and ANOSIM are both sensitive to imbalance, and ANOSIM more so. Beyond
about 3:1 the permutation null is distorted and the test can reject for reasons
unrelated to group difference.

Combined with heterogeneous dispersion — which imbalance often accompanies — the
result is uninterpretable.

**What to do.** Report the imbalance (`validate_inputs()` computes the ratio).
Consider subsampling the larger group as a sensitivity check, reporting both. Do
not report only the analysis that gave the answer you preferred.

---

## 6. Multiple testing across taxa

Testing 600 taxa at α = 0.05 yields ~30 false positives by construction.

**Correction is necessary and not sufficient.** Two further points:

- **The denominator matters as much as the method.** Filtering 600 taxa to 143
  before testing legitimately increases power — but the report must say 143 were
  tested, not 600. See `10`.
- **Correction happens once, across the family of tests you are drawing
  conclusions from.** Correcting within each of three DA methods separately and
  then reporting the union of the three hit lists has no valid error rate.
  Report the consensus (`da_compare()`), or pre-specify one primary method.

Also: **post-hoc tests after a non-significant omnibus** inflate the error rate.
`posthoc_dunn()` refuses unless forced.

---

## 7. Choosing a threshold after seeing the results

A garden-of-forking-paths problem, and it is easy to do accidentally.

- trying prevalence 0.05, 0.10 and 0.20 and reporting the one with most hits
- trying three DA methods and reporting the one that found your taxon of interest
- trying two rarefaction depths and keeping the significant one

Any of these invalidates the reported FDR. The p-value assumes the analysis was
specified before the data were seen.

**What to do.** Fix filters and the primary method before testing. If you do
explore several, report all of them as a sensitivity analysis — which is
informative and honest — rather than presenting one as the analysis. The
benchmark sweep machinery (`14`, `15`) is exactly this done deliberately.

---

## 8. Structural zeros

A zero means either "absent" or "present but not sequenced deeply enough to see".
These license different claims, and conflating them is how a shallow-sequenced
sample becomes a biological finding.

ANCOM-BC2's `struc_zero = TRUE` declares a taxon absent from an entire group a
structural zero and handles it separately. That is the right machinery — but a
structural zero is a **strong claim**, and at low depth it is frequently
undersampling.

**What to do.** Check the depth of the samples where the taxon is absent. If they
are the shallow ones, the structural-zero call is a depth artefact. Report
structural zeros separately from fold changes, and state the depth range.

---

## 9. Taxonomic name instability

**The commonest cause of "the paper's taxa disappeared".**

Names change between database releases, and the same organism appears under
different names in two analyses:

- **Phylum renamings (2021 ICNP validation):** Bacteroidetes → Bacteroidota,
  Firmicutes → Bacillota, Proteobacteria → Pseudomonadota, Actinobacteria →
  Actinomycetota, and so on.
- **The *Lactobacillus* split (Zheng et al. 2020):** one genus became 25.
  *L. casei* → *Lacticaseibacillus casei*, *L. plantarum* →
  *Lactiplantibacillus plantarum*, *L. reuteri* → *Limosilactobacillus reuteri*.
- **Individual reassignments:** *Eubacterium rectale* → *Agathobacter rectalis*,
  *Propionibacterium acnes* → *Cutibacterium acnes*, *Clostridium difficile* →
  *Clostridioides difficile*, *Bacteroides vulgatus* → *Phocaeicola vulgatus*.
- **Merged labels:** SILVA reports `Escherichia-Shigella` as one genus.

**What to do.** Always state the reference database **and release** — SILVA 138
versus 138.2 versus GTDB r220 matters enormously. In benchmark mode,
`harmonize_taxon()` maps through `MFG_TAXON_SYNONYMS` before comparing, and every
unresolvable label is reported rather than counted as a miss (`12`, `15`).

Note also: the same ASV sequence can receive different genus assignments from
different databases at the same confidence threshold. The sequence is the stable
identity; the name is not.

---

## 10. Placeholder phylogenies

Demo 3 and Demo 6 both attach `rtree()` — a random tree — to satisfy phyloseq's
slot for teaching purposes. UniFrac and Faith's PD on a random topology return a
matrix and a vector that look **entirely valid** and encode nothing.

This is not a hypothetical: the teaching pattern is widely copied.

**What to do.** `mfg_mark_tree_real(ps, TRUE/FALSE)`. Until marked, `faith_pd()`,
`compute_distance("unifrac")` and `plot_phylogeny()` all refuse. Unmarked is
treated as unknown, not as real.

Also check rooting: unweighted UniFrac and PD with `include_root = TRUE` both
depend on it.

---

## 11. Correlation networks on proportions

Correlating relative abundances produces strong spurious negative edges, because
proportions must sum to one. A network built this way describes arithmetic.

`cooccurrence_network()` refuses TSS input and requires CLR. Even on CLR values,
treat edges as hypotheses — SpiecEasi and SPARCC are the methods designed for
this problem, and a Spearman network with 1,770 pairwise tests needs correction
(which the function applies and reports).

---

## 12. Table orientation

Taxa-are-rows backwards is silent: every per-sample statistic becomes a per-taxon
statistic and nothing errors.

`mfg_detect_orientation()` decides by matching names against the taxonomy and
metadata, **not** by comparing dimensions — shape is not evidence, since a study
can have more samples than taxa. When the evidence ties, it refuses and asks.

---

## 13. Alphabetical reference levels

A character column is not a factor, so its reference level is alphabetical rather
than chosen. Every model coefficient, every DESeq2 log2 fold change and every
ANCOM-BC2 lfc is relative to it.

Getting this wrong does not error. It **flips the sign of every reported effect**,
and the analysis still looks fine.

```r
sample_data(ps)$status <- factor(sample_data(ps)$status,
                                 levels = c("healthy", "disease"))
```

`validate_inputs()` lists the character columns that look categorical. Before
interpreting any direction, confirm the reference with
`DESeq2::resultsNames(dds)` or the ANCOM-BC2 column names.

---

## 14. Batch effects

Sequencing run, extraction kit, plate position and collection date all leave
signals that can exceed the biological effect.

Filtering does not remove them. **A batch variable belongs in the model**, not in
the QC step — as a covariate (`06`, `07`) or a PERMANOVA term.

The failure mode to watch for: batch confounded with the variable of interest. If
all cases were run on plate 1 and all controls on plate 2, the design cannot
separate them and no statistical method will. That is a design finding, and the
honest report says the effect is not estimable.

---

## 15. Pseudo-replication

Related to (3) but distinct: the experimental unit is not always the sample.

- Cagemates share a microbiome through coprophagy. **The cage is the unit**, not
  the mouse. Ten mice in two cages is n = 2, not n = 10.
- Several aliquots of one stool sample are technical replicates, not biological.
- Multiple time points on one subject are one subject.

Treating samples as the unit when the unit is coarser inflates n, sometimes by an
order of magnitude, and every p-value with it.

---

## 16. Independent filtering read as missing data

DESeq2 sets `padj` to NA for low-mean-count taxa. That is a deliberate exclusion
that increases power on the remainder, not a failure to compute.

Dropping those rows and reporting "31 of 41 significant" is correct. Reporting
"31 of 49" without noting that 8 were not tested is not. `da_deseq2()` reports
both numbers and carries the caveat.

---

## 17. Over-interpreting alpha diversity

"Diversity was reduced in disease" is one of the most-reported and least-specific
findings in the field.

- Shannon confounds richness and evenness. Report both, or say which one moved.
- Richness is strongly depth-dependent (`03`), so an unrarefied richness
  difference may be a depth difference.
- A null alpha result with n = 10 per group is uninformative, not evidence of no
  effect (`12`).
- Alpha diversity is a summary statistic with no mechanism attached. A composition
  or DA result says more about biology than a Shannon difference does.

---

## 18. Stress and variance explained left off ordinations

An NMDS at stress 0.25 shows structure that is an artefact of projection. A PCoA
of axes explaining 17% looks identical to one explaining 80%. Without these
numbers the figure cannot be assessed, and readers will assess it anyway.

`plot_ordination_std()` puts them in the subtitle. See `09`.

---

## 19. Unseeded stochastic steps

Rarefaction subsamples at random. Permutation tests sample at random. Without a
fixed seed, a rerun disagrees with the report and neither run is wrong — and the
disagreement is indistinguishable from a real methodological difference, which
makes it poison for benchmarking.

Every stochastic function here takes and logs a seed, and
`capture_provenance()` recovers them all from the run log.

**The benchmark diagnostic:** re-run the identical baseline configuration twice.
If the two baselines differ, the divergence is noise and no sweep axis will
explain it (`15`, `stochastic.unseeded`).

---

## 20. Functional prediction read as measurement

PICRUSt2 infers function from 16S using reference genomes. It is a **prediction
conditional on the reference database**, not a measurement of gene content.
Confidence is low for poorly represented lineages and for strain-level functional
variation, which is precisely where interesting differences often live.

Report predicted pathways as predicted. Shotgun metagenomics measures gene
content; 16S plus PICRUSt2 does not.
