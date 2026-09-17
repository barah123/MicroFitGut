# Benchmark: Lucas et al. (2017), the ant-built home

A full benchmark-mode run against a published study, in the same form as the
Pérez-Losada benchmark: reconstruct the analysis, score the paper's stated
claims, and report which survive. Study **A3** of the validation matrix, the
minimal-input test, and the study that showed "data available" and "data
reusable" are different tests.

| | |
|---|---|
| **Paper** | Lucas J, Bill B, Stevenson B, Kaspari M. The microbiome of the ant-built home: the microbial communities of a tropical arboreal ant and its nest. *Ecosphere.* 2017;8(2):e01639. [10.1002/ecs2.1639](https://doi.org/10.1002/ecs2.1639) |
| **Raw data** | Dryad [10.5061/dryad.ph2c5](https://doi.org/10.5061/dryad.ph2c5): one HDF5 BIOM, 14.4 MB. **Metadata not in the deposit**, recovered from NCBI BioSample SAMN04576300 to SAMN04576371 |
| **Analyzed** | 48,137 taxa × 79 samples, filtered to 11,278 taxa × 72 samples |
| **Design** | Four sample types from 20 *Azteca trigona* colonies, Panama: ant 17, refuse 17, nest 20, soil 18 |
| **Hypothesis** | The ant microbiome is distinct from its environment, and the nests ants build carry microbiota reflecting the ants rather than the surrounding soil |
| **Verdict** | **Incomplete.** 1 of 2 scorable claims reproduced, 1 not licensed, 2 outside the denominator. No sweep |

---

## 1. Study card

**6 fields stated · 2 inferred · 15 absent.**

The paper does not state a depth floor or a prevalence filter, so both are absent
card fields and would become sweep axes.

### The deposit is available but not reusable

The Dryad record contains exactly one file: `all_otu_table.biom`, 48,137 taxa by
79 samples. It has a DOI, it downloads, it parses, and it contains a valid OTU
table. **It would pass any automated data-availability check.**

It is also, on its own, impossible to analyze. The sample identifiers are bare
numeric barcodes, `100044`, `100181`, `100205`, and every claim in the paper is a
comparison between ants, refuse, nests and soil. Nothing in the file says which
sample is which.

Inspecting the HDF5 structure confirms this is not a parser failure:

```
/sample/metadata          <- group exists
  (no datasets)           <- and is empty
```

The paper's data availability statement points somewhere else entirely:

> "All microbial data have been uploaded and are available at NCBI's BioSamples
> databank (accession nos. SAMN04576300-SAMN04576371)."

Those records carry the key:

```
BioSample: SAMN04576300; Sample name: 100008
  /Replicate="ant 1"
  /host="Azteca trigona ant"
  /geographic location="Panama: Gigante"
```

`Sample name` matches the BIOM's identifiers. Fetching all 72 records rebuilds
the design and the study becomes analyzable. Seven of the 79 BIOM samples have no
BioSample record at all and were dropped.

**This is the finding worth carrying into the thesis.** The data passed every
mechanical test of availability and was still unusable until a cross-repository
join recovered the grouping. A data-availability statement that resolves is not
the same as a deposit that can be reanalyzed, and no automated check in current
use distinguishes them.

### The paper's four claims

| id | Claim | Evidence in the paper |
|---|---|---|
| `types-differ` | The four sample types have distinct communities | stated as distinct microbial communities |
| `lactobacillus-ants` | *Lactobacillus* dominates the ant community | identified as the dominant ant-associated genus |
| `nest-like-soil` | Nest microbiota resemble soil, not ants | contrary to the paper's own prediction |
| `dimorphic-colonies` | Ant communities are dimorphic across colonies | stated qualitatively |

---

## 2. Reconstruction

Quality control at minimum depth 1,000 and prevalence 0.05: **48,137 taxa reduced
to 11,278**, 72 samples retained.

Intake flagged, before the metadata was joined: depth varying **72,449-fold**,
since one sample had a single read, 95.2% sparsity, and 10,077 taxa with zero
reads in every sample.

---

## 3. Results

### 3.1 The four sample types differ, but the claim is not licensed

| Test | Result |
|---|---|
| PERMANOVA (4 types, Bray-Curtis) | R² = 0.2391, F = 7.12, **p = 0.001** |
| **betadisper** | **p = 0.027, dispersion heterogeneous** |

Sample type explains 24% of variation in community composition, which is a large
effect. But dispersion differs significantly between types, so PERMANOVA cannot
separate a shift in composition from a difference in spread. The paper's claim
that the ant microbiome is distinct is **not contradicted**. It is not licensed
by this test.

Given what the data shows, heterogeneous dispersion is expected rather than
surprising. The paper's own thesis is that ant microbiomes vary dramatically
across colonies while soil does not. Unequal spread between groups is part of the
finding, not an artifact, which is exactly why PERMANOVA alone cannot carry the
claim here. The same pattern appears in A2, where the paper's stability claim and
its composition claim constrain each other in the same way.

### 3.2 Mean between-type Bray-Curtis

| | ant | refuse | nest | soil |
|---|---|---|---|---|
| **ant** | 0.675 | 0.810 | 0.857 | **0.878** |
| **refuse** | 0.810 | 0.662 | 0.810 | 0.829 |
| **nest** | 0.857 | 0.810 | 0.739 | **0.831** |
| **soil** | 0.878 | 0.829 | 0.831 | 0.782 |

Ant and soil are the **most dissimilar pair in the matrix** at 0.878, consistent
with the paper's headline.

The nest claim is directly supported: **nest to soil 0.831 is less than nest to
ant 0.857**. Nests are built from ant exudates and chewed plant fiber, and the
paper predicted their microbiota would reflect the ants. They resemble the soil
instead. The paper reports this as a result contrary to its own prediction, and
the reanalysis agrees.

### 3.3 *Lactobacillus* drives the variation

| Type | mean | sd | range |
|---|---|---|---|
| **ant** | **33.46%** | 23.58 | 0.56 to 78.12 |
| refuse | 2.30% | 1.52 | 0.50 to 6.11 |
| nest | 1.31% | 1.55 | 0.26 to 7.45 |
| soil | 1.18% | 0.89 | 0.11 to 3.38 |

*Lactobacillus* is **25-fold more abundant in ants** than in any other sample
type, Wilcoxon ant against soil p = 2.2e-07, and its range within ants spans
almost the entire possible scale, from 0.6% to 78% of the community. That is the
dimorphic-across-colonies pattern the paper describes, and it is the single
clearest result in this study.

---

## 4. Verdict

```
=== Benchmark verdict: INCOMPLETE ===

Reproduced (1):
  + lactobacillus-ants : higher in ant, 33.5% vs 1.2% in soil, p = 2.2e-07

Not reproduced (0)

Not licensed (1):
  ! types-differ       : significant at p = 0.001, but betadisper p = 0.027,
                         so PERMANOVA does not license a composition claim

Outside the denominator (2):
  ? nest-like-soil     : schema_gap - the distance matrix supports it
                         (nest-soil 0.831 < nest-ant 0.857), but no claim
                         type encodes "X is more similar to Y than to Z"
  ? dimorphic-colonies : schema_gap - qualitative

Scorable: 2 of 4.  Reproduced: 1 of 2.
```

Nothing was contradicted. The verdict is incomplete because three of four claims
could not be settled cleanly: one because the evidence does not license it, two
because they are comparative statements the schema cannot express.

### What this means

Every one of the paper's four claims is supported by the data in the direction it
states. Only one could be scored as reproduced, and the reasons for the other
three are all about the instrument or the test rather than about the paper.

That is worth stating plainly, because a bare verdict of `incomplete` reads
worse than the evidence warrants. A reader of this study's row in a summary table
would see 1 of 2, and the underlying position is that the paper's findings hold
and the tool could not certify three of them.

### A claim type the schema is missing

`nest-like-soil` is a **relative similarity** claim: *X resembles Y more than Z*.
The distance matrix answers it unambiguously here, but no claim type encodes it,
so it could only be recorded as `schema_gap`.

This is the same gap `da_null` filled for S1, in a different domain. Comparative
distance claims are common in community ecology; "the gut community resembled the
diet more than the environment" is a standard form. A `beta_closer` claim type
taking `subject`, `closer_to` and `than_to` would settle it from the same matrix
that is already computed. Together with A1's two claims of the same shape, this
is three of the seven `schema_gap` claims in the series, and it is the single
highest-value extension to the schema.

Not implemented here. Recorded as the next schema extension.

### Limitations

1. **Colony structure not modeled.** BioSample's `Replicate` field encodes the
   colony number, `ant 1`, `nest 1`, so the design is genuinely paired within 20
   colonies. This reanalysis treats samples as independent, exactly as the
   PERMANOVA above does. Extracting colony and permuting within it would be the
   correct next step and would likely reduce the apparent effect.
2. **Seven samples unmatched.** 79 in the BIOM, 72 in BioSample.
3. **No sweep.** Single reconstruction.
4. **QC thresholds are the tool's, not the paper's.** The paper does not state a
   depth floor or prevalence filter, so both would become sweep axes.

---

## 5. Defects this run exposed

### #12: `Rank1..Rank7` taxonomy has no usable ranks

`import_biom()` returned seven taxonomy columns named `Rank1` through `Rank7`.
Functionally this is the same failure as S4's single-column lineage: there is no
`Genus` column, so `tax_glom(ps, "Genus")` errors and matching a published genus
name finds nothing.

`mfg_split_lineage()` now renames generic `RankN` columns to the standard schema
when the count matches, strips rank prefixes, and normalizes empty strings and
`"unidentified"` to `NA`. On this dataset that resolved **30,644 of 48,137 taxa
to genus** and made `tax_glom` work, yielding 1,676 genera.

### #13: the repeated-measures detector fired on a grouping factor

With BioSample metadata joined, `validate_inputs()` reported:

> Repeated measures detected: `host` has up to 20 samples per level across 4
> levels.

`host` names the **sample type**, `Azteca trigona ant`, soil and so on. It is the
study's grouping variable. Fitting `(1 | host)` would have modeled the effect
under investigation as noise, which is the exact inversion of what the guard
exists to prevent.

The detector's only structural test was "not unique per row", which any
categorical column passes. A subject variable partitions samples into **many
small** clusters; a grouping factor into **few large** ones. A candidate is now
rejected when one level holds more samples than there are levels, rejections are
reported rather than silent, and a variable named explicitly by the caller is
still honored, with a warning about its shape.

Regression checked against three known cases: `dietswap` with 38 subjects and a
maximum of 6 samples each, the skin cohort with 128 patients and a maximum of 5,
and S1, correctly reporting no repeated measures.

---

## 6. Reproducing this

```r
b  <- phyloseq::import_biom("A3/all_otu_table.biom")
b  <- mfg_split_lineage(b)                 # Rank1..Rank7 -> Kingdom..Species

# Metadata is NOT in the deposit: rebuild it from NCBI BioSample.
md <- read.csv("A3/biosample_metadata.csv", stringsAsFactors = FALSE)
md$type <- sub("\\s*\\d+$", "", md$replicate)     # "ant 1" -> "ant"
rownames(md) <- md$SampleID
ov <- intersect(sample_names(b), md$SampleID)     # 72 of 79
ps <- prune_samples(ov, b)
sample_data(ps) <- sample_data(md[ov, , drop = FALSE])

ps <- filter_prevalence(prune_low_depth(ps, 1000), 0.05)
ps <- mfg_set_normalization(ps, "raw")
validate_inputs(ps, group_var = "type")

d <- compute_distance(tss_transform(ps), "bray")
run_permanova(tss_transform(ps), ~ type, dist_obj = d)
check_dispersion(d, mfg_meta(ps)$type)
```

Seeds fixed at 42. Every value above is written to `run_log.csv`.
