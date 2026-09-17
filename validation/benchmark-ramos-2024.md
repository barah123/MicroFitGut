# Benchmark: Ramos et al. (2024), bisphenol S and the murine gut microbiome

Study **S5** of the validation matrix. It exercises MicroFitGut's own MetaPhlAn
profile parser, and it is the first study in the series whose headline claim is a
null result.

| | |
|---|---|
| **Paper** | Ramos C, et al. Integrated Metagenomic and Metabolomic Analysis of In Vitro Murine Gut Microbiome. *Metabolites.* 2024;14(12):713. |
| **Data** | Zenodo [10.5281/zenodo.13917959](https://doi.org/10.5281/zenodo.13917959), 20 MetaPhlAn profiles and 80 HUMAnN files |
| **Design** | In vitro anaerobic cultures, 5 donor animals, 2 conditions (BPS, vehicle), 2 timepoints |
| **Verdict** | **Incomplete**. Two claims held, none failed, two unadjudicated |

---

## 1. Role in the matrix

S5 tests the functional profile path and MicroFitGut's own
`build_phyloseq_from_profile()` parser against real MetaPhlAn output, rather than
against the pre-assembled objects curatedMetagenomicData supplies. It is also the
only study in the set whose central claim is that an exposure did **not** change
the microbiome.

---

## 2. What the deposit contains

The Zenodo record holds per-sample MetaPhlAn profiles rather than one merged
table, which is how MetaPhlAn is normally run but not how downstream tools expect
to receive it. Twenty samples, coded by donor animal (A to E), condition (S or V)
and timepoint (1 or 2), giving a balanced 5 by 2 by 2 design.

`build_phyloseq_from_profile()` accepted a single per-sample profile directly,
returning 18 taxa for one sample. Merging the twenty into one table and parsing
that gave **50 taxa across 20 samples**, with the hierarchy correctly collapsed:
`leaf_lineages()` reduced 160 lineage rows to the 50 non-overlapping leaves,
which is the behavior that prevents the same organism being counted once per rank.

### A format inconsistency in the deposit

Nineteen of the twenty files are MetaPhlAn community profiles in the standard
lineage format. One, `46.AS2_...`, is in a different format entirely: rank columns
(`kingdom`, `phylum`, ... `abundance`) rather than pipe-separated lineage strings.
Sample AS2 exists **only** in that format, so it cannot be skipped without losing
a sample.

Merging the deposit therefore requires detecting and converting two formats. A
naive read of the directory fails, which is what happened on the first attempt
here. This is a data deposition issue rather than a tool defect, but it is the
kind of friction that makes reuse expensive, and it is worth recording alongside
A3's missing metadata as a second instance of the same theme.

---

## 3. Results

The MetaPhlAn output is relative abundance, so every sample sums to 100 and the
normalization state is `tss`. Count based differential abundance is unavailable,
as in S4.

### 3.1 The null claim holds

The paper states that a supraphysiologic BPS dose "did not overtly distort the
metagenomic or metabolomic profiles of exposed cultures compared to controls."

| Test | Result |
|---|---|
| PERMANOVA, permuted within donor | R² = 0.0213, F = 0.392, **p = 0.493** |
| betadisper | p = 0.472, dispersion homogeneous |

`permanova_sentence()` wrote the conclusion:

> cond did not explain a significant share of variation in bray distance
> (PERMANOVA R² = 0.021, F = 0.39, p = 0.493, 999 permutations).

Treatment explains about 2% of variation in community composition, and that share
is not distinguishable from zero. The claim holds.

Dispersion being homogeneous matters here in the opposite direction to the usual
case. A null PERMANOVA on heterogeneous dispersion would be weak evidence, since
the test could be failing to detect a location shift masked by unequal spread.
With homogeneous dispersion the null result is interpretable.

### 3.2 Inter-animal variation dominates

| Term | R² | p |
|---|---|---|
| Donor animal | **0.6756** | 0.001 |
| BPS treatment | 0.0213 | 0.493 |

Which animal a culture came from explains **68%** of variation in composition.
The treatment explains 2%. This is a thirty-fold difference, and it supports the
paper's closing point that inter-animal variation needs accounting for before an
exposure effect could be detected at this sample size.

Recorded as unadjudicated because the claim is qualitative, though the evidence
is unambiguous.

### 3.3 Lactobacillus direction holds, significance does not

| | BPS | vehicle |
|---|---|---|
| Mean *Lactobacillus* | 5.247% | 4.231% |

Higher under BPS exposure, matching the paper's direction, but Wilcoxon
p = 0.955. The paper's claim came from a discriminant model rather than a
univariate test, and a discriminant model can weight a feature that no single
test would flag. Direction agrees; the univariate evidence for it is absent.

Scored as held, because the claim encoded is directional. The lack of univariate
significance is recorded as a qualification.

---

## 4. Verdict

```
=== Benchmark verdict: INCOMPLETE ===

Claims that held:
  + no-distortion     : PERMANOVA p = 0.493, dispersion homogeneous
  + lactobacillus-bps : higher under BPS, 5.25% vs 4.23%

Claims NOT adjudicated:
  ? inter-animal-variation : qualitative, though donor R2 = 0.676 supports it
  ? metabolite-profile     : metabolomics, outside MicroFitGut's scope
```

Nothing failed.

---

## 5. Defect this run exposed

### #16, the subject detector missed an abbreviated column name

The metadata column naming the donor animal is `subj`. The candidate pattern
required the full word `subject`, so `subj` was never considered and
`validate_inputs()` reported no clustering at all in a design where one animal
contributes four of twenty samples.

The consequence is not hypothetical. Without the detection, a PERMANOVA on this
design would be run unstratified, treating four cultures from one animal as four
independent observations. Given that donor explains 68% of the variance, that
would be a serious error, and it would inflate significance in a study whose
finding is a null result.

The pattern now matches on stems rather than whole words: `subj` in place of
`subject`, `particip` in place of `participant`. Real metadata abbreviates, and a
detector that only recognizes the unabbreviated form will miss a large share of
deposited studies.

Regression checked against five cases: `dietswap` (subject), the skin cohort
(patient), S5 (subj), A3 (host correctly rejected as a grouping factor), and S1
(correctly no clustering).

---

## 6. Limitations

1. **HUMAnN files not analyzed.** The deposit contains 80 files of gene family
   and pathway abundances. Parsing them requires
   `build_phyloseq_functional()` with a separate pathway taxonomy file, which the
   deposit does not include in the expected form. The functional path therefore
   remains only partly exercised.
2. **Small n.** Five animals, twenty cultures. A treatment effect smaller than
   the donor effect would not be detectable at this size, which is the paper's
   own conclusion.
3. **In vitro cultures, not animals.** The design measures what BPS does to
   cultured communities, not to a host.
4. **Metabolomic claims outside scope.** Two of the paper's four abstract claims
   concern metabolites, which MicroFitGut does not analyze.
5. **No sweep.** Single reconstruction.

---

## 7. Reproducing this

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
