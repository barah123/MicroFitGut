# Benchmark: Kim et al. (2020), gut microbiome in pulmonary arterial hypertension

A full benchmark-mode run against a published study, in the same form as the
Pérez-Losada benchmark: reconstruct the analysis, score the paper's stated
claims, and report which survive. Study **S4** of the validation matrix, the
BIOM intake path, and the first study to exercise a deposited file that is
already normalized.

| | |
|---|---|
| **Paper** | Kim S, Rigatto K, Gazzana MB, et al. Altered Gut Microbiome Profile in Patients With Pulmonary Arterial Hypertension. *Hypertension.* 2020;75(4):1063-1071. [10.1161/HYPERTENSIONAHA.119.14294](https://doi.org/10.1161/HYPERTENSIONAHA.119.14294) |
| **Raw data** | Dryad [10.5061/dryad.stqjq2c03](https://doi.org/10.5061/dryad.stqjq2c03): taxonomic BIOM, KEGG-orthology BIOM, QIIME metadata |
| **Analyzed** | The taxonomic BIOM, 249 taxa × 30 samples after removing one duplicated metadata row |
| **Design** | Case-control: PAH patients against reference subjects, shotgun metagenomics. Paper states 18 PAH and 13 reference |
| **Hypothesis** | The gut microbiome of patients with pulmonary arterial hypertension differs from that of reference subjects, with butyrate- and propionate-producing bacteria depleted in disease |
| **Verdict** | **Partially reproduced.** 3 of 7 claims reproduced, 4 not reproduced, 2 outside the denominator. No sweep |

---

## 1. Study card

**8 fields stated · 2 inferred · 13 absent.**

S4 replaced an earlier candidate that proved to be a false positive: the corpus
classifier had read "Bracken" in a filename as the Kraken2/Bracken tool when it
was in fact the author's surname. This study was chosen instead because it ships
**two** BIOM files, a taxonomic table and a KEGG-orthology table, so it tests the
BIOM parser and functional data in one study.

The most consequential absent field is the normalization of the deposit itself,
which the paper does not state and which turns out to decide what can be run at
all. See section 2.

### The paper's nine claims

| id | Claim | Evidence in the paper |
|---|---|---|
| `butyrivibrio-control` | *Butyrivibrio* increased in reference | named among groups "increased in reference cohort" |
| `akkermansia-control` | *Akkermansia* increased in reference | named among the same group |
| `bacteroides-control` | *Bacteroides* increased in reference | named among the same group |
| `lachnospiraceae-control` | Lachnospiraceae increased in reference | named among the same group |
| `coprococcus-control` | *Coprococcus* increased in reference | named among the same group |
| `eubacterium-control` | *Eubacterium* increased in reference | named among the same group |
| `beta-differs` | Communities differ between PAH and reference | stated as an altered gut microbiome profile |
| `random-forest` | A classifier separates the groups | 83% predictive accuracy |
| `virome` | Enterococcal and Lactococcal phage differ | stated in the results |

---

## 2. Reconstruction

Three findings before any analysis ran, all of which determine how the study can
be scored.

**The BIOM holds relative abundances, not counts.** Every sample sums to exactly
100. `validate_inputs()` reported `Values: relative abundances, non-integer
present` and a depth range of 100 to 100, contradicting the `raw` state that had
been asserted. Every count-based differential abundance method is therefore
correctly refused:

```
da_ancombc2(...) -> REFUSED: Invalid normalization for this analysis.
                    object is: tss (total sum scaling (relative abundance))
```

This is the normalization guard working exactly as designed, and it confines the
reanalysis to rank-based tests. It also means **the paper's own methods are not
reproducible from its deposit**, which is a finding about the deposit rather than
about the paper's conclusions.

**The metadata has a duplicated row.** `Kim1` appears twice with identical
content, at lines 2 and 32 of the deposited metadata. Excluded explicitly and
logged.

**The deposit is one sample short of the paper.** The paper reports 18 PAH and 13
reference subjects, n = 31. After removing the duplicate, the deposit holds 18
PAH and **12** reference, n = 30, and the BIOM carries 30 samples. One reference
subject is absent from the deposited data with no note. Recorded as a
discrepancy, not a reanalysis failure.

---

## 3. Results

Reanalysis used Wilcoxon rank-sum tests on the deposited relative abundances,
the only family the data supports.

### 3.1 Named taxa: two reproduce, four do not

The paper states that six butyrate- and propionate-producing groups were
"increased in reference cohort".

| Taxon | Reference | PAH | Higher in | p | Status |
|---|---|---|---|---|---|
| *Butyrivibrio* | 4.39% | 0.00% | reference | **5.7e-05** | `reproduced` |
| *Akkermansia* | 0.99% | 0.87% | reference | **0.009** | `reproduced` |
| *Bacteroides* | 23.41% | 18.71% | reference | 0.059 | `not_reproduced` |
| Lachnospiraceae | 15.25% | 11.54% | reference | 0.518 | `not_reproduced` |
| *Coprococcus* | 2.03% | 2.80% | **PAH** | 0.465 | `not_reproduced` |
| *Eubacterium* | 17.77% | 18.08% | **PAH** | 0.662 | `not_reproduced` |

**Only *Butyrivibrio* and *Akkermansia* reach significance.**

Two of the six reverse direction. Neither reversal is significant, so the honest
reading is that *Coprococcus* and *Eubacterium* show **no difference** in this
reanalysis rather than an opposite one. The paper's claim as written is that they
were increased in the reference cohort, and that does not hold.

Two more, *Bacteroides* at p = 0.059 and Lachnospiraceae at p = 0.518, point the
right way without reaching significance. These were originally scored as holding
on direction alone. Under the significance-first rule they do not: the paper
asserted a difference, and the reanalysis does not detect one. The Lachnospiraceae
case makes the point plainly, since p = 0.518 carries no information about
direction at all.

### 3.2 Beta diversity: the claim holds, and is licensed

| Test | Result |
|---|---|
| PERMANOVA (Bray-Curtis) | R² = 0.0673, F = 2.019, **p = 0.004** |
| **betadisper** | **p = 0.091, dispersion homogeneous** |

This is the instructive contrast with A0, where a significant PERMANOVA was
accompanied by heterogeneous dispersion and the composition claim was therefore
not licensed. Here dispersion is homogeneous, so the significant PERMANOVA
**does** support "the microbial communities differ", and the claim is scored
`reproduced` rather than `not_licensed`.

The effect is small. R² = 0.067 means group membership explains under 7% of
variation in community composition.

### 3.3 Alpha diversity: not a claim, but worth recording

Shannon diversity is significantly **lower** in PAH, median 2.572 against 3.057,
Welch t-test p = 0.004. The paper makes no alpha-diversity claim in its abstract,
so this is a reanalysis observation rather than a claim test, and it is not
scored.

---

## 4. Verdict

```
=== Benchmark verdict: PARTIALLY REPRODUCED ===

Reproduced (3):
  + butyrivibrio-control    : higher in reference, p = 5.7e-05
  + akkermansia-control     : higher in reference, p = 0.009
  + beta-differs            : PERMANOVA p = 0.004, dispersion homogeneous

Not reproduced (4):
  - bacteroides-control     : direction agrees, p = 0.059, not significant
  - lachnospiraceae-control : direction agrees, p = 0.518, not significant
  - coprococcus-control     : direction reversed, p = 0.465, no difference
  - eubacterium-control     : direction reversed, p = 0.662, no difference

Not licensed (0)

Outside the denominator (2):
  ? random-forest : not_attempted  - classifier not re-fitted
  ? virome        : data_absent    - no viral data deposited

Scorable: 7 of 9.  Reproduced: 3 of 7.
```

**This verdict changed on rescoring.** The report originally read 5 held and 2
failed, with *Bacteroides* and Lachnospiraceae counted as holding on direction.
Under the significance-first rule both are `not_reproduced`, which moves the
study from 5 of 7 to 3 of 7.

The two items outside the denominator mean different things and are counted
separately. The random-forest claim **could** be tested but was not, which is a
protocol deviation under the analyst's control. The virome claim **cannot** be
tested, because the deposit contains no viral data, so a claim the paper makes in
its results is unverifiable from what was shared.

### What this means

The paper's central composition claim reproduces and is licensed, which is the
strongest form of agreement available here. Its two strongest taxon claims also
reproduce, and *Butyrivibrio* does so emphatically, at 4.39% in reference against
0.00% in PAH.

The weaker four do not. None of them is contradicted: two show no difference in
the claimed direction, and two show no difference in the opposite direction. What
this reanalysis establishes is that four of the six named groups do not carry
detectable univariate evidence in the deposited data. The paper's list reads as
six equally supported findings, and the deposit supports two.

### Limitations

1. **Deposited data is already normalized.** Relative abundances only, so
   ANCOM-BC2, ALDEx2 and DESeq2 are all unavailable and the reanalysis rests on
   rank-based tests.
2. **One reference subject missing** from the deposit relative to the paper's
   stated n = 31.
3. **Functional (KEGG) claims not tested.** The KO BIOM was downloaded, but the
   paper's pathway claims covering arginine, proline and ornithine synthesis,
   TMA and TMAO, and purine metabolism require KO-to-pathway mapping that was not
   performed.
4. **Virome claims untestable.** No viral data in the deposit.
5. **No sweep.** Single reconstruction. The prevalence filter and taxonomic level
   were not varied, and with four claims resting on p-values between 0.059 and
   0.662 a sweep would be more informative here than in most of the series.
6. **Small n.** 18 against 12 after exclusions. The marginal *Bacteroides* result
   at p = 0.059 would plausibly move either way with a few more samples, and that
   is precisely why it is reported as not reproduced rather than as refuted.

---

## 5. Defects this run exposed

### #11: QIIME-style BIOM taxonomy arrives in a single column

`phyloseq::import_biom()` on this file returned a taxonomy table with **one
column named `Rank1`**, holding the entire lineage as a string:

```
k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Leuconostocaceae;g__Weissella;s__Weissella_confusa
```

Everything downstream that works by rank then has nothing to work with.
`tax_glom(ps, "Genus")` errors, matching a published genus name against the
taxonomy finds nothing, and `validate_inputs()` reports `ranks: Rank1`. Nothing
crashes loudly. The analysis simply cannot proceed, and the reason is not obvious
from the error.

This is the commonest shape of deposited 16S and shotgun BIOM data, so it was
handled rather than left to the caller. Added `mfg_split_lineage()`:

- Detects `;` or `|` as the separator rather than assuming one.
- Strips `k__` and `d__` rank prefixes.
- **Pads short lineages with NA rather than recycling**, so a lineage that stops
  at Family does not acquire a fabricated Genus.
- Logs the separator used and how many taxa resolved to genus and species.

On this dataset it resolved **249 of 249 taxa to genus and species**. A2 and A3
both ship QIIME-style BIOMs and would have hit the same wall.

---

## 6. Reproducing this

```r
b  <- phyloseq::import_biom("S4/taxon_biom_table_PH.biom")
md <- read.csv("S4/Qiime_MG_Metadata_PH_Data_comma_separated_test.txt",
               check.names = FALSE, stringsAsFactors = FALSE)
names(md)[1] <- sub("^#", "", names(md)[1])
md <- md[!duplicated(md$SampleID), ]          # removes the duplicated Kim1 row
rownames(md) <- md$SampleID

ps <- phyloseq(otu_table(b), tax_table(b),
               sample_data(md[sample_names(b), , drop = FALSE]))
ps <- mfg_split_lineage(ps)                   # Rank1 -> Kingdom..Species
ps <- mfg_set_normalization(ps, "tss")        # the deposit IS relative abundance
validate_inputs(ps, group_var = "Group")

d <- compute_distance(ps, "bray")
run_permanova(ps, ~ Group, dist_obj = d, check_disp = TRUE)
check_dispersion(d, mfg_meta(ps)$Group)
```

Seeds fixed at 42. Every value above is written to `run_log.csv`.
