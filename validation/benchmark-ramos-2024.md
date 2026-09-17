# Benchmark: Ramos et al. (2024), bisphenol S and the murine gut microbiome

A full benchmark-mode run against a published study, in the same form as the
Pérez-Losada benchmark: reconstruct the analysis, score the paper's stated
claims, and report which survive. Study **S5** of the validation matrix. It
exercises MicroFitGut's own MetaPhlAn profile parser, and it is the only study in
the series whose headline claim is a null result.

| | |
|---|---|
| **Paper** | Ramos C, et al. Integrated Metagenomic and Metabolomic Analysis of In Vitro Murine Gut Microbiome. *Metabolites.* 2024;14(12):713. [10.3390/metabo14120713](https://doi.org/10.3390/metabo14120713) |
| **Raw data** | Zenodo [10.5281/zenodo.13917959](https://doi.org/10.5281/zenodo.13917959): 20 MetaPhlAn profiles and 80 HUMAnN files |
| **Analyzed** | The 20 MetaPhlAn profiles, merged, giving 50 taxa × 20 samples |
| **Design** | In vitro anaerobic cultures, 5 donor animals, 2 conditions (BPS, vehicle), 2 timepoints. Balanced 5 × 2 × 2 |
| **Hypothesis** | A supraphysiologic dose of bisphenol S does not substantially distort the gut microbial community of cultured murine communities |
| **Verdict** | **Incomplete.** 1 of 2 scorable claims reproduced, 1 not reproduced, 2 outside the denominator. No sweep |

---

## 1. Study card

**7 fields stated · 3 inferred · 13 absent.**

S5 tests the functional profile path and MicroFitGut's own
`build_phyloseq_from_profile()` parser against real MetaPhlAn output, rather than
against the pre-assembled objects curatedMetagenomicData supplies.

The paper's central claim being a null one changes what the absent fields cost.
For a positive claim an unstated filter threatens a false negative. For a null
claim it threatens the opposite: an underpowered analysis will agree with a null
claim for the wrong reason, so the design's power matters more than usual and is
treated below as a limitation rather than a footnote.

### The paper's four claims

| id | Claim | Evidence in the paper |
|---|---|---|
| `no-distortion` | BPS does not distort the community | "did not overtly distort the metagenomic or metabolomic profiles of exposed cultures compared to controls" |
| `lactobacillus-bps` | *Lactobacillus* is higher under BPS | identified by a discriminant model |
| `inter-animal-variation` | Inter-animal variation dominates | stated in the discussion |
| `metabolite-profile` | The metabolite profile is unchanged | metabolomics result |

---

## 2. Reconstruction

The Zenodo record holds **per-sample** MetaPhlAn profiles rather than one merged
table, which is how MetaPhlAn is normally run but not how downstream tools expect
to receive it. Twenty samples, coded by donor animal (A to E), condition (S or V)
and timepoint (1 or 2).

`build_phyloseq_from_profile()` accepted a single per-sample profile directly,
returning 18 taxa for one sample. Merging the twenty into one table and parsing
that gave **50 taxa across 20 samples**, with the hierarchy correctly collapsed:
`leaf_lineages()` reduced 160 lineage rows to the 50 non-overlapping leaves,
which is the behavior that prevents the same organism being counted once per
rank.

MetaPhlAn output is relative abundance, so every sample sums to 100 and the
normalization state is `tss`. Count-based differential abundance is unavailable,
as in S4.

### A format inconsistency in the deposit

Nineteen of the twenty files are MetaPhlAn community profiles in the standard
lineage format. One, `46.AS2_...`, is in a different format entirely: rank
columns (`kingdom`, `phylum`, through `abundance`) rather than pipe-separated
lineage strings. Sample AS2 exists **only** in that format, so it cannot be
skipped without losing a sample.

Merging the deposit therefore requires detecting and converting two formats. A
naive read of the directory fails, which is what happened on the first attempt
here. This is a data deposition issue rather than a tool defect, but it is the
kind of friction that makes reuse expensive, and it sits alongside A3's missing
metadata as a second instance of the same theme.

---

## 3. Results

### 3.1 The null claim reproduces, and is licensed

| Test | Result |
|---|---|
| PERMANOVA, permuted within donor | R² = 0.0213, F = 0.392, **p = 0.493** |
| betadisper | p = 0.472, dispersion homogeneous |

`permanova_sentence()` wrote the conclusion:

> cond did not explain a significant share of variation in bray distance
> (PERMANOVA R² = 0.021, F = 0.39, p = 0.493, 999 permutations).

Treatment explains about 2% of variation in community composition, and that share
is not distinguishable from zero. The claim reproduces.

Dispersion being homogeneous matters here in the opposite direction to the usual
case. A null PERMANOVA under heterogeneous dispersion would be weak evidence,
since the test could be failing to detect a location shift masked by unequal
spread. With homogeneous dispersion the null result is interpretable, and the
claim is licensed rather than merely unrefuted.

### 3.2 Inter-animal variation dominates

| Term | R² | p |
|---|---|---|
| Donor animal | **0.6756** | 0.001 |
| BPS treatment | 0.0213 | 0.493 |

Which animal a culture came from explains **68%** of variation in composition.
The treatment explains 2%. That is a thirtyfold difference, and it supports the
paper's closing point that inter-animal variation needs accounting for before an
exposure effect could be detected at this sample size.

The claim is recorded as `schema_gap` because it is qualitative and no claim type
expresses "X varies more than Y", though the evidence is unambiguous.

### 3.3 *Lactobacillus*: direction agrees, evidence does not

| | BPS | vehicle |
|---|---|---|
| Mean *Lactobacillus* | 5.247% | 4.231% |

Higher under BPS exposure, matching the paper's direction, but Wilcoxon
**p = 0.955**.

The paper's claim came from a discriminant model rather than a univariate test,
and a discriminant model can weight a feature that no single test would flag. The
direction agrees. The univariate evidence for it is absent, and p = 0.955 is as
close to no information as a p-value gets.

---

## 4. Verdict

```
=== Benchmark verdict: INCOMPLETE ===

Reproduced (1):
  + no-distortion     : PERMANOVA p = 0.493, dispersion homogeneous

Not reproduced (1):
  - lactobacillus-bps : direction agrees, p = 0.955, no univariate evidence

Not licensed (0)

Outside the denominator (2):
  ? inter-animal-variation : schema_gap   - qualitative; donor R2 = 0.676
                             is consistent with it
  ? metabolite-profile     : out_of_scope - metabolomics

Scorable: 2 of 4.  Reproduced: 1 of 2.
```

**This verdict changed on rescoring.** The report originally read two claims
held, with `lactobacillus-bps` scored on direction alone. At p = 0.955 the
direction of the effect is not evidence of anything, and under the
significance-first rule the claim is `not_reproduced`.

That change matters more here than in the other studies, because it is the
clearest example in the series of what the old rule was doing wrong. A p-value of
0.955 means the observed difference is almost exactly what noise would produce.
Recording it as support for the paper's claim, on the basis that 5.247% exceeds
4.231%, was scoring the sign of a coin flip.

### What this means

The paper's headline claim, that BPS does not distort these communities,
reproduces cleanly and is licensed by a homogeneous dispersion check. That is the
strongest result available for a null claim.

The *Lactobacillus* claim does not reproduce univariately. This is not a
contradiction of the paper, which derived it from a discriminant model, and a
discriminant model legitimately uses information a single rank-sum test discards.
What the reanalysis establishes is that the claim does not survive the simplest
test one would apply to it, which is worth knowing about a feature named in an
abstract.

### Limitations

1. **Power, and it cuts against the headline claim.** Five animals and twenty
   cultures, with donor explaining 68% of variance. A treatment effect smaller
   than the donor effect would not be detectable at this size. The null result is
   consistent with no effect and also with an effect this design cannot see. The
   paper says as much itself, and the reanalysis cannot separate the two.
2. **HUMAnN files not analyzed.** The deposit contains 80 files of gene family
   and pathway abundances. Parsing them requires `build_phyloseq_functional()`
   with a separate pathway taxonomy file, which the deposit does not include in
   the expected form. The functional path is therefore only partly exercised.
3. **In vitro cultures, not animals.** The design measures what BPS does to
   cultured communities, not to a host.
4. **Metabolomic claims outside scope.** Two of the paper's four abstract claims
   concern metabolites, which MicroFitGut does not analyze.
5. **No sweep.** Single reconstruction.

---

## 5. Defects this run exposed

### #16: the subject detector missed an abbreviated column name

The metadata column naming the donor animal is `subj`. The candidate pattern
required the full word `subject`, so `subj` was never considered and
`validate_inputs()` reported no clustering at all, in a design where one animal
contributes four of twenty samples.

The consequence is not hypothetical. Without the detection, a PERMANOVA on this
design would run unstratified, treating four cultures from one animal as four
independent observations. Given that donor explains 68% of the variance, that
would be a serious error.

**Direction of bias.** This defect is **anti-conservative**: an undetected
cluster means the naive test runs and significance is inflated. In this study
that happened to work against the paper, because the paper's claim here is a null
one and inflated significance would have contradicted it. In a study claiming a
difference the same defect would have made the claim *more* likely to be scored
as reproduced. The direction is a property of the claim, not of the defect, and
the series-level statement that every defect biased toward non-reproduction does
not hold for this one.

The pattern now matches on stems rather than whole words: `subj` in place of
`subject`, `particip` in place of `participant`. Real metadata abbreviates, and a
detector that recognizes only the unabbreviated form will miss a large share of
deposited studies.

Regression checked against five cases: `dietswap` (subject), the skin cohort
(patient), S5 (subj), A3 (host correctly rejected as a grouping factor), and S1
(correctly no clustering).

---

## 6. Reproducing this

```r
D  <- "S5/metaphlan/MetaPhlan"
fs <- grep("^\\._", list.files(D, pattern = "\\.tsv$"), value = TRUE, invert = TRUE)

# Nineteen files are lineage format; one (AS2) is rank columns. Both are needed.
PRE <- c("k__","p__","c__","o__","f__","g__","s__","t__")
read_one <- function(p) {
  if (grepl("^kingdom\t", readLines(p, n = 1))) {
    d  <- read.delim(p, check.names = FALSE, stringsAsFactors = FALSE)
    ab <- d[[ncol(d)]]; rk <- d[, seq_len(ncol(d) - 1), drop = FALSE]
    lin <- apply(rk, 1, function(r) {
      r <- r[nzchar(r) & !is.na(r)]
      paste(paste0(PRE[seq_along(r)], r), collapse = "|")
    })
    setNames(ab[nzchar(lin)], lin[nzchar(lin)])
  } else {
    d <- read.delim(p, comment.char = "#", header = FALSE,
                    col.names = c("lineage","abund"), stringsAsFactors = FALSE)
    setNames(d$abund, d$lineage)
  }
}
# merge, write with a clade_name column, then:
ps <- build_phyloseq_from_profile("merged.tsv")
ps <- mfg_set_normalization(ps, "tss")     # MetaPhlAn output is relative abundance

dist <- compute_distance(ps, "bray")
run_permanova(ps, ~ cond, dist_obj = dist, strata = "subj")
check_dispersion(dist, factor(mfg_meta(ps)$cond))
run_permanova(ps, ~ subj, dist_obj = dist, check_disp = FALSE)   # donor effect
```

Seeds fixed at 42. Every value above is written to `run_log.csv`.
