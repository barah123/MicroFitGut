# Benchmark: Berer et al. (2017), gut microbiota in MS-discordant twins

A full benchmark-mode run against a published study, in the same form as the
Pérez-Losada benchmark: reconstruct the analysis, score the paper's stated
claims, and report which survive. Study **A4** of the validation matrix, the
last of the ten. It tests the loose-table assembly path, and it is the one
design in the set where stratified permutation is genuinely valid.

| | |
|---|---|
| **Paper** | Berer K, Gerdes LA, Cekanaviciute E, et al. Gut microbiota from multiple sclerosis patients enables spontaneous autoimmune encephalomyelitis in mice. *PNAS.* 2017;114(40):10719-10724. [10.1073/pnas.1711233114](https://doi.org/10.1073/pnas.1711233114) |
| **Raw data** | UC San Francisco Dash [10.7272/Q6RX997G](https://doi.org/10.7272/Q6RX997G): an OTU table and a metadata file, both tab-separated text |
| **Analyzed** | 8,856 OTUs × 115 samples, 68 human and 47 mouse |
| **Design** | 34 monozygotic twin pairs discordant for MS, plus germ-free mice colonized from a subset of donors |
| **Hypothesis** | Gut microbiota from multiple sclerosis patients differ from those of their healthy twins, and transferring them to germ-free mice raises the incidence of spontaneous autoimmune encephalomyelitis |
| **Verdict** | **Incomplete.** 3 of 3 scorable claims reproduced, 3 outside the denominator. No sweep |

---

## 1. Study card

**9 fields stated · 2 inferred · 12 absent.**

A4 was chosen for the **loose table** path: an OTU table and a metadata file as
separate text files, rather than a packaged object. That shape covers 73 of the
136 papers in the 16S corpus, so if it is painful the main study inherits the
pain. It is also the only study in the set with a genuinely paired design.

### The paper's six claims

| id | Claim | Evidence in the paper |
|---|---|---|
| `humans-no-major-diff` | MS and healthy twins do not differ overall | "no major differences in the overall microbial profiles" |
| `akkermansia-untreated-ms` | *Akkermansia* is increased in untreated MS twins | "a significant increase in some taxa such as *Akkermansia* in untreated MS twins" |
| `mice-differ` | Colonized mice differ by donor disease state | stated in the results |
| `sutterella-mice` | *Sutterella* differs between mouse groups | described as "an organism shown to induce a protective immunoregulatory profile in vitro" |
| `autoimmunity-incidence` | MS-colonized mice show higher autoimmunity incidence | the paper's central biological finding |
| `il10` | IL-10 production is reduced | immunological assay |

---

## 2. Reconstruction

The deposit is two tab-separated text files. Both carry a UTF-8 byte-order mark,
and the metadata uses classic Mac carriage-return line endings, which `wc -l`
reports as zero lines. R's `read.delim` handles the carriage returns, so the
metadata read correctly at 115 rows by 7 columns.

The OTU table did not. It produced **8,857 rows and zero columns, without an
error**. Defect #20 below.

Once read: 8,856 OTUs by 115 samples, 68 human and 47 mouse, 57 control and 58
MS, across 34 twin pairs. The taxonomy arrived as a semicolon-separated lineage
in a trailing `taxonomy` column, which `mfg_split_lineage()` resolved to the
seven standard ranks.

### Why stratification is valid here

In A2, species never varied within an animal, so permuting within animal could
not test it and the tool refused. A4 is the opposite case, and it shows the guard
is not simply blocking stratification:

```
twin pairs among humans: 34 | disease varies within pair: TRUE
```

The pairs are **discordant** by design: each monozygotic pair contains one twin
with MS and one without. Disease state therefore varies within the stratum, the
permutation is meaningful, and `run_permanova(..., strata = "Twin_code")` runs.

This is what a paired design is for. Twins share genetics and much of their early
environment, so permuting within pair removes that shared variation from the null
rather than leaving it in the residual.

---

## 3. Results

### 3.1 The human null claim reproduces

| Test | Result |
|---|---|
| PERMANOVA, permuted within twin pair | R² = 0.0153, **p = 0.210** |
| betadisper | p = 0.134, dispersion homogeneous |

Disease state explains 1.5% of variation in composition, and that is not
distinguishable from zero. Dispersion is homogeneous, so the null is
interpretable rather than a failure to detect a shift hidden by unequal spread.

### 3.2 The *Akkermansia* claim cannot be tested from this deposit

The paper reports a significant increase in some taxa such as *Akkermansia* in
**untreated** MS twins.

The word untreated is load bearing. MS is commonly treated with immunomodulatory
drugs that alter the gut microbiota, so the claim is about a subset of MS twins,
not all of them.

The deposited metadata has seven columns: `Human_mouse`, `Human_subject`,
`Weeks`, `Disease_state`, `Weeks_disease_state`, `Twin_code`, `Donor`. **There is
no treatment column.** The untreated subset cannot be identified.

Testing all MS twins against all controls gives *Akkermansia* at 0.152% against
0.129%, higher in MS, matching the direction, at p = 0.869. That is not a test of
the paper's claim. The contrast guard scored it `contrast_mismatch` and named the
mismatch:

> claim is about 'untreated_MS_vs_control'; evidence describes 'human_disease_state'

Reporting p = 0.869 as a failure to reproduce would have been wrong, and it is the
kind of wrong that looks entirely reasonable on the page.

### 3.3 Both mouse claims reproduce

| Claim | Result |
|---|---|
| Mouse profiles differ by donor disease state | R² = 0.0571, **p = 0.001**, betadisper p = 0.806, homogeneous |
| *Sutterella* differs | MS-donor 0.467%, control-donor 0.633%, **p = 1.3e-05** |

The colonized mice separate by the disease state of their human donor, and
because dispersion is homogeneous that is a licensed composition claim rather
than a dispersion artifact.

*Sutterella* is higher in mice colonized from healthy twins. The paper describes
it as an organism shown to induce a protective immunoregulatory profile in vitro,
so more of it in the healthy-donor mice is the direction the paper's argument
requires. The effect is small in absolute terms, about 0.17 percentage points,
but consistent: p = 1.3e-05 across 19 OTUs.

Note that the mouse samples come from only 4 donor twin pairs, so this rests on a
narrow base regardless of the sample count.

---

## 4. Verdict

```
=== Benchmark verdict: INCOMPLETE ===

Reproduced (3):
  + humans-no-major-diff : R2 = 0.0153, p = 0.210, dispersion homogeneous
  + mice-differ          : R2 = 0.0571, p = 0.001, dispersion homogeneous
  + sutterella-mice      : higher in control-donor mice, p = 1.3e-05

Not reproduced (0)
Not licensed  (0)

Outside the denominator (3):
  ? akkermansia-untreated-ms : contrast_mismatch - the deposit has no
                               treatment column, so the untreated subset
                               cannot be identified
  ? autoimmunity-incidence   : out_of_scope - mouse disease phenotype, not
                               microbiome data
  ? il10                     : out_of_scope - immunological assay

Scorable: 3 of 6.  Reproduced: 3 of 3.
```

Nothing was contradicted. Every claim that could be tested reproduced.

### What this means

This is the cleanest paired design in the series, and it is the only study where
the repeated-measures machinery could be used as intended rather than as a
refusal. The human null claim reproduces under a permutation scheme that respects
the twin structure, which is a stronger result than the same null obtained by
ignoring it.

Three claims sit outside the denominator for two quite different reasons, and the
distinction matters more here than the count does. Two are simply not microbiome
measurements: a disease incidence and an immunological assay. The tool has
nothing to say about either, and recording them as anything other than out of
scope would be misleading.

The third is the interesting one. *Akkermansia* in untreated MS twins is a claim
the deposit cannot support, not because the taxon is missing but because the
covariate that defines the subgroup was never deposited. A reanalysis that
ignored the word "untreated" would have produced p = 0.869 and read as a
contradiction of the paper. **That is the failure mode this guard exists for**,
and this is the clearest instance of it in the series.

### Limitations

1. **Treatment status is absent**, so the paper's principal taxon claim is not
   testable from the deposit.
2. **The mouse arm rests on 4 donor pairs.** 47 mouse samples, but they derive
   from at most four human donors per arm.
3. **Repeated sampling of mice not modeled.** `Weeks` records longitudinal
   sampling; this analysis treats mouse samples as independent within donor,
   which overstates precision.
4. **Phenotype and immunology claims are outside scope.** The paper's central
   biological finding, higher autoimmunity incidence in MS-colonized mice, is a
   disease-incidence measurement, not a microbiome measurement. The benchmark
   therefore cannot speak to the paper's main argument at all, only to the
   microbiome observations that support it.
5. **No sweep.**

---

## 5. Defects this run exposed

### #20: the delimiter was detected from a comment line

`mfg_read_table()` chose the delimiter by counting tabs and commas in the **first
line of the file**. The OTU table's first line is the standard QIIME export
header:

```
# Constructed from biom file
```

That line contains neither a tab nor a comma, so the detector fell through to
comma. Reading a tab-separated file as comma-separated puts every row into a
single field, that field becomes the row names, and the result is a data frame
with **8,857 rows and zero columns**.

No error is raised. The object exists, prints, and has the right number of rows.
The failure surfaces later, somewhere unrelated, as a dimension mismatch.

`mfg_read_table()` now finds the header rather than assuming it is line 1. It
reads the first 50 lines, strips any byte-order mark, and selects the first
non-blank line that is not a pure comment, or the last commented line that
contains a delimiter when the header itself is commented, which is the
`#OTU ID<TAB>...` convention. It also refuses to return a zero-column frame,
since that result is never correct.

Regression checked on a plain CSV and a plain TSV, both unchanged.

This matters beyond one file. "# Constructed from biom file" is what QIIME writes
on every table exported to text, so any deposit in that form would have hit it.

---

## 6. Reproducing this

```r
# The OTU table begins "# Constructed from biom file"; the delimiter must be
# detected from the header, not line 1.
ot <- mfg_read_table("A4/Max_Planck_Twin_OTU_table.txt")
md <- mfg_read_table("A4/Max_Planck_Twin_metadata.txt")   # CR line endings

lin <- ot$taxonomy; ot$taxonomy <- NULL
tax <- matrix(lin, ncol = 1, dimnames = list(rownames(ot), "Rank1"))

ov <- intersect(colnames(ot), rownames(md))
ps <- phyloseq(otu_table(as.matrix(ot[, ov]), taxa_are_rows = TRUE),
               sample_data(md[ov, ]), tax_table(tax))
ps <- mfg_set_normalization(mfg_split_lineage(ps), "raw")

hum <- mfg_prune_samples(sample_names(ps)[mfg_meta(ps)$Human_mouse == "Human"], ps)
hum <- filter_prevalence(mfg_set_normalization(hum, "raw"), 0.05)
d   <- compute_distance(tss_transform(hum), "bray")

# strata is valid here: the twin pairs are MS-discordant
run_permanova(tss_transform(hum), ~ Disease_state, dist_obj = d,
              strata = "Twin_code", check_disp = FALSE)
check_dispersion(d, factor(mfg_meta(hum)$Disease_state))
```

Seeds fixed at 42. Every value above is written to `run_log.csv`.
