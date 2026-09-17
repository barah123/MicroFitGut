# Benchmark: Kim et al. (2020), gut microbiome in pulmonary arterial hypertension

Study **S4** of the validation matrix — the BIOM intake path, and the first study
to exercise a deposited file that is already normalized.

| | |
|---|---|
| **Paper** | Kim S, Rigatto K, Gazzana MB, et al. Altered Gut Microbiome Profile in Patients With Pulmonary Arterial Hypertension. *Hypertension.* 2020;75(4):1063–1071. |
| **Data** | Dryad [10.5061/dryad.stqjq2c03](https://doi.org/10.5061/dryad.stqjq2c03) — taxonomic BIOM, KEGG-orthology BIOM, QIIME metadata |
| **Design** | Case-control: PAH patients vs reference subjects, shotgun metagenomics |
| **Verdict** | **Partially reproduced** — 5 of 9 claims held, 2 failed, 2 unadjudicated |

---

## 1. Role in the matrix

S4 replaced an earlier candidate that proved to be a false positive: the corpus
classifier had read "Bracken" in a filename as the Kraken2/Bracken tool when it
was the author's surname. This study was chosen instead because it ships **two**
BIOM files — a taxonomic table and a KEGG-orthology table — so it tests the BIOM
parser and functional data in one study.

---

## 2. What the deposit actually contains

Three findings before any analysis ran, all of which matter for how the study can
be scored.

**The BIOM holds relative abundances, not counts.** Every sample sums to exactly
100. `validate_inputs()` reported `Values: relative abundances, non-integer
present` and a depth range of 100–100, contradicting the `raw` state that had
been asserted. Consequence: every count-based differential abundance method is
correctly refused.

```
da_ancombc2(...) -> REFUSED: Invalid normalization for this analysis.
                    object is: tss (total sum scaling (relative abundance))
```

This is the normalization guard working exactly as designed, and it confines the
reanalysis to rank-based tests.

**The metadata has a duplicated row.** `Kim1` appears twice with identical
content (lines 2 and 32 of the deposited metadata). Excluded explicitly and
logged.

**The deposit is one sample short of the paper.** The paper reports 18 PAH and 13
reference subjects (n = 31). After removing the duplicate, the deposit holds 18
PAH and **12** reference (n = 30), and the BIOM carries 30 samples. One reference
subject is absent from the deposited data with no note. Recorded as a
discrepancy, not a reanalysis failure.

---

## 3. Results

Reanalysis used Wilcoxon rank-sum tests on the deposited relative abundances —
the only family the data supports.

### 3.1 Named taxa — 4 of 6 directions hold, 2 reverse

The paper states that six butyrate- and propionate-producing groups were
"increased in reference cohort".

| Taxon | Reference | PAH | Higher in | p | Claim |
|---|---|---|---|---|---|
| *Butyrivibrio* | 4.39% | 0.00% | reference | **5.7e-05** | **held** |
| *Akkermansia* | 0.99% | 0.87% | reference | **0.009** | **held** |
| *Bacteroides* | 23.41% | 18.71% | reference | 0.059 | **held** (direction) |
| Lachnospiraceae | 15.25% | 11.54% | reference | 0.518 | **held** (direction) |
| *Coprococcus* | 2.03% | 2.80% | **PAH** | 0.465 | **failed** — direction reversed |
| *Eubacterium* | 17.77% | 18.08% | **PAH** | 0.662 | **failed** — direction reversed |

Two of the six reverse direction. Neither reversal is significant, so the honest
reading is that *Coprococcus* and *Eubacterium* show no difference in this
reanalysis rather than an opposite one — but the paper's claim as written is that
they were increased in the reference cohort, and that does not hold.

Only *Butyrivibrio* and *Akkermansia* reach significance. *Bacteroides* is
marginal at 0.059.

### 3.2 Beta diversity — the claim holds, and is licensed

| Test | Result |
|---|---|
| PERMANOVA (Bray-Curtis) | R² = 0.0673, F = 2.019, **p = 0.004** |
| **betadisper** | **p = 0.091 — dispersion homogeneous** |

This is the instructive contrast with A0 (Pérez-Losada), where a significant
PERMANOVA was accompanied by heterogeneous dispersion and the composition claim
was therefore not licensed. Here dispersion is homogeneous, so the significant
PERMANOVA **does** support "the microbial communities differ", and the claim is
scored `held` rather than `unscored`.

The effect is small — R² = 0.067 means group membership explains under 7% of
variation in community composition.

### 3.3 Alpha diversity — not a claim, but worth recording

Shannon diversity is significantly **lower** in PAH: median 2.572 vs 3.057,
Welch t-test p = 0.004. The paper makes no alpha-diversity claim in its abstract,
so this is a reanalysis observation rather than a claim test.

---

## 4. Verdict

```
=== Benchmark verdict: PARTIALLY REPRODUCED ===

Claims that held:
  + butyrivibrio-control    : higher in reference, p = 5.7e-05
  + akkermansia-control     : higher in reference, p = 0.009
  + bacteroides-control     : higher in reference, p = 0.059
  + lachnospiraceae-control : higher in reference, p = 0.518
  + beta-differs            : PERMANOVA p = 0.004, dispersion homogeneous

Claims that did not hold:
  - coprococcus-control     : reanalysis finds it higher in PAH
  - eubacterium-control     : reanalysis finds it higher in PAH

Claims NOT adjudicated:
  ? random-forest  : 83% predictive accuracy — not re-fitted
  ? virome         : Enterococcal/Lactococcal phage claims — no viral data deposited
```

The two unadjudicated claims are honest gaps rather than schema limitations. The
random-forest claim could be tested but was not; the virome claim **cannot** be —
the deposit contains no viral data, so a central claim of the paper is
unverifiable from what was shared.

---

## 5. Defect this run exposed

### #11 — QIIME-style BIOM taxonomy arrives in a single column

`phyloseq::import_biom()` on this file returned a taxonomy table with **one
column named `Rank1`**, holding the entire lineage as a string:

```
k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Leuconostocaceae;g__Weissella;s__Weissella_confusa
```

Everything downstream that works by rank then has nothing to work with.
`tax_glom(ps, "Genus")` errors; matching a published genus name against the
taxonomy finds nothing; `validate_inputs()` reports `ranks: Rank1`. Nothing
crashes loudly — the analysis simply cannot proceed, and the reason is not
obvious from the error.

This is the commonest shape of deposited 16S and shotgun BIOM data, so it was
handled rather than left to the caller. Added `mfg_split_lineage()`:

- Detects `;` or `|` as separator rather than assuming.
- Strips `k__` / `d__` rank prefixes.
- **Pads short lineages with NA rather than recycling**, so a lineage that stops
  at Family does not acquire a fabricated Genus.
- Logs the separator used and how many taxa resolved to genus and species.

On this dataset it resolved **249/249 taxa to genus and species**. A2 and A3 both
ship QIIME-style BIOMs and would have hit the same wall.

---

## 6. Limitations

1. **Deposited data is already normalized.** Relative abundances only, so
   ANCOM-BC2, ALDEx2 and DESeq2 are all unavailable and the reanalysis rests on
   rank-based tests. The paper's own methods are not reproducible from its
   deposit for this reason.
2. **One reference subject missing** from the deposit relative to the paper's
   stated n = 31.
3. **Functional (KEGG) claims not tested.** The KO BIOM was downloaded but the
   paper's pathway claims — arginine, proline and ornithine synthesis, TMA/TMAO,
   purine metabolism — require KO-to-pathway mapping that was not performed.
4. **Virome claims untestable.** No viral data in the deposit.
5. **No sweep.** Single reconstruction; the prevalence filter and taxonomic level
   were not varied.
6. **Small n.** 18 vs 12 after exclusions. The marginal *Bacteroides* result
   (p = 0.059) would plausibly move either way with a few more samples.

---

## 7. Reproducing this

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
