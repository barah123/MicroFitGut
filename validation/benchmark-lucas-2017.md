# Benchmark: Lucas et al. (2017), the ant-built home

Study **A3** of the validation matrix — the minimal-input test, and the study
that showed "data available" and "data reusable" are different tests.

| | |
|---|---|
| **Paper** | Lucas J, Bill B, Stevenson B, Kaspari M. The microbiome of the ant-built home: the microbial communities of a tropical arboreal ant and its nest. *Ecosphere.* 2017;8(2):e01639. |
| **Data** | Dryad [10.5061/dryad.ph2c5](https://doi.org/10.5061/dryad.ph2c5) — one HDF5 BIOM, 14.4 MB |
| **Metadata** | **Not in the deposit.** Recovered from NCBI BioSample SAMN04576300–371 |
| **Design** | Four sample types from 20 *Azteca trigona* colonies, Panama: ant, refuse, nest, soil |
| **Verdict** | **Incomplete** — 1 claim held, 3 unadjudicated, 0 failed |

---

## 1. The deposit is available but not reusable

The Dryad record contains exactly one file: `all_otu_table.biom`, 48,137 taxa ×
79 samples. It has a DOI, it downloads, it parses, and it contains a valid OTU
table. It would pass any automated data-availability check.

It is also, on its own, **impossible to analyse**. The sample identifiers are
bare numeric barcodes — `100044`, `100181`, `100205` — and every claim in the
paper is a comparison between ants, refuse, nests and soil. Nothing in the file
says which sample is which.

Inspecting the HDF5 structure confirms this is not a parser failure:

```
/sample/metadata          <- group exists
  (no datasets)           <- and is empty
```

The paper's data availability statement points somewhere else entirely:

> "All microbial data have been uploaded and are available at NCBI's BioSamples
> databank (accession nos. SAMN04576300–SAMN04576371)."

Those records carry the key:

```
BioSample: SAMN04576300; Sample name: 100008
  /Replicate="ant 1"
  /host="Azteca trigona ant"
  /geographic location="Panama: Gigante"
```

`Sample name` matches the BIOM's identifiers. Fetching all 72 records rebuilds
the design — **ant 17, refuse 17, nest 20, soil 18** — and the study becomes
analysable.

**This is the finding worth carrying into the thesis.** The data passed every
mechanical test of availability and was still unusable until a cross-repository
join recovered the grouping. Seven of the 79 BIOM samples have no BioSample
record at all and were dropped.

---

## 2. Reconstruction

Quality control at minimum depth 1,000 and prevalence 0.05: **48,137 → 11,278
taxa**, 72 samples retained.

Intake flagged, before the metadata was joined: depth varying **72,449-fold**
(one sample had a single read), 95.2% sparsity, and 10,077 taxa with zero reads
in every sample.

---

## 3. Results

### 3.1 The four sample types differ — but the claim is not licensed

| Test | Result |
|---|---|
| PERMANOVA (4 types, Bray-Curtis) | R² = 0.2391, F = 7.12, **p = 0.001** |
| **betadisper** | **p = 0.027 — dispersion heterogeneous** |

Sample type explains 24% of variation in community composition, which is a large
effect. But dispersion differs significantly between types, so PERMANOVA cannot
separate a shift in composition from a difference in spread. The paper's claim
that the ant microbiome is "distinct" is **not contradicted** — it is not
licensed by this test.

Given what the data shows, heterogeneous dispersion is expected rather than
surprising: the paper's own thesis is that ant microbiomes vary *dramatically*
across colonies while soil does not. Unequal spread between groups is part of
the finding, not an artefact — which is exactly why PERMANOVA alone cannot carry
the claim here.

### 3.2 Mean between-type Bray-Curtis

|  | ant | refuse | nest | soil |
|---|---|---|---|---|
| **ant** | 0.675 | 0.810 | 0.857 | **0.878** |
| **refuse** | 0.810 | 0.662 | 0.810 | 0.829 |
| **nest** | 0.857 | 0.810 | 0.739 | **0.831** |
| **soil** | 0.878 | 0.829 | 0.831 | 0.782 |

Ant and soil are the **most dissimilar pair in the matrix** (0.878), consistent
with the paper's headline.

The nest claim is directly supported: **nest–soil 0.831 < nest–ant 0.857**. Nests
are built from ant exudates and chewed plant fibre, and the paper predicted their
microbiota would reflect the ants. They resemble the soil instead.

### 3.3 *Lactobacillus* drives the variation — strongly supported

| Type | mean | sd | range |
|---|---|---|---|
| **ant** | **33.46%** | 23.58 | 0.56 – 78.12 |
| refuse | 2.30% | 1.52 | 0.50 – 6.11 |
| nest | 1.31% | 1.55 | 0.26 – 7.45 |
| soil | 1.18% | 0.89 | 0.11 – 3.38 |

*Lactobacillus* is **25-fold more abundant in ants** than in any other sample
type (Wilcoxon ant vs soil p = 2.2e-07), and its range within ants spans almost
the entire possible scale — from 0.6% to 78% of the community. That is the
"dimorphic across colonies" pattern the paper describes, and it is the single
clearest result in this study.

---

## 4. Verdict

```
=== Benchmark verdict: INCOMPLETE ===

Claims that held:
  + lactobacillus-ants : higher in ant, 33.5% vs 1.2% in soil, p = 2.2e-07

Claims NOT adjudicated:
  ? types-differ       : dispersion heterogeneous — PERMANOVA does not license
                         a composition claim
  ? nest-like-soil     : qualitative; the distance matrix supports it
                         (nest-soil 0.831 < nest-ant 0.857) but no claim type
                         encodes "X is more similar to Y than to Z"
  ? dimorphic-colonies : qualitative
```

Nothing failed. The verdict is `incomplete` rather than `reproduced` because
three of four claims could not be settled by the scoring machinery — one because
the evidence does not license it, two because they are comparative statements the
schema cannot express.

---

## 5. Defects this run exposed

### #12 — `Rank1..Rank7` taxonomy has no usable ranks

`import_biom()` returned seven taxonomy columns named `Rank1` through `Rank7`.
Functionally this is the same failure as S4's single-column lineage: there is no
`Genus` column, so `tax_glom(ps, "Genus")` errors and matching a published genus
name finds nothing.

`mfg_split_lineage()` now renames generic `RankN` columns to the standard schema
when the count matches, strips rank prefixes, and normalises empty strings and
`"unidentified"` to `NA`. On this dataset that resolved **30,644 of 48,137 taxa
to genus** and made `tax_glom` work, yielding 1,676 genera.

### #13 — the repeated-measures detector fired on a grouping factor

With BioSample metadata joined, `validate_inputs()` reported:

> Repeated measures detected: `host` has up to 20 samples per level across 4
> levels.

`host` names the **sample type** (`Azteca trigona ant`, soil, and so on). It is
the study's grouping variable. Fitting `(1 | host)` would have modelled the
effect under investigation as noise — the exact inversion of what the guard
exists to prevent.

The detector's only structural test was "not unique per row", which any
categorical column passes. A subject variable partitions samples into **many
small** clusters; a grouping factor into **few large** ones. A candidate is now
rejected when one level holds more samples than there are levels, rejections are
reported rather than silent, and a variable named explicitly by the caller is
still honoured with a warning about its shape.

Regression-checked against three known cases: `dietswap` (38 subjects, max 6),
the skin cohort (128 patients, max 5), and S1 (correctly no repeated measures).

---

## 6. A claim type the schema is missing

`nest-like-soil` is a **relative similarity** claim: *X resembles Y more than Z*.
The distance matrix answers it unambiguously here, but no claim type encodes it,
so it could only be recorded as `unscored`.

This is the same gap `da_null` filled for S1, in a different domain. Comparative
distance claims are common in community ecology — "the gut community resembled
the diet more than the environment" is a standard form. A `beta_closer` claim
type taking `subject`, `closer_to` and `than_to` would settle it from the same
matrix that is already computed.

Not implemented here; recorded as the next schema extension.

---

## 7. Limitations

1. **Colony structure not modelled.** BioSample's `Replicate` field encodes the
   colony number (`ant 1`, `nest 1`), so the design is genuinely paired within
   20 colonies. This reanalysis treats samples as independent, exactly as the
   PERMANOVA above does. Extracting colony and permuting within it would be the
   correct next step and would likely reduce the apparent effect.
2. **Seven samples unmatched.** 79 in the BIOM, 72 in BioSample.
3. **No sweep.** Single reconstruction.
4. **QC thresholds are the tool's, not the paper's.** The paper does not state a
   depth floor or prevalence filter, so both are `absent` card fields and would
   become sweep axes.

---

## 8. Reproducing this

```r
b  <- phyloseq::import_biom("A3/all_otu_table.biom")
b  <- mfg_split_lineage(b)                 # Rank1..Rank7 -> Kingdom..Species

# metadata is NOT in the deposit — rebuild it from NCBI BioSample
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
