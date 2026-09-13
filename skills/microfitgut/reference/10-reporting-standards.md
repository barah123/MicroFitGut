# Reporting standards

What every analysis summary must state, regardless of what was run.

## The overriding rule

**A conclusion that cannot be traced to a run identifier and a stored artifact
does not get reported.** If a number is not in `run_log.csv` or in a file under
`tables/`, it does not go in the report.

Two checks enforce this:

1. `mfg_traceability_check(report_text)` — mechanical. Extracts every numeric
   literal from the draft and looks for it in the run log and every saved CSV.
   Deliberately noisy: it over-flags rather than passing a fabricated statistic.
2. The **`verifier` subagent** — judgement. Adjudicates the flagged list, checks
   that the methods text describes what the code actually did, and checks every
   citation.

Both run before a report is final. The traceability check produces the candidate
list; the verifier decides.

## Mandatory content

`assemble_summary()` checks for these and lists what is missing rather than
silently omitting it.

### 1. n per group, after filtering

Not the number collected — the number analysed. Both, if they differ.

A result reported against the enrolled n when the analysed n is smaller overstates
the evidence. State them separately: "of 344 samples collected, 301 passed
filtering (BE 121, BT 122, FA 101 before filtering; 109, 104, 88 after)".

### 2. Every sample dropped, and why

The **list**, not the count. `qc_exclusions(log)` produces it. A reader who cannot
see which samples left cannot judge whether the exclusions were selective.

Include the rule that dropped each one, and the justification for each threshold.

### 3. The normalization, per analysis

Not one sentence for the run — one per stage, because a run legitimately holds
several states of the same data (`03`). `mfg_normalization_sentence(ps)` writes
each, including the rarefaction depth, the seed, and how many samples the
rarefaction dropped.

### 4. The exact test, with its assumptions

Name the test, say why it was chosen, and give the assumption-check results that
justified it — Shapiro-Wilk, Bartlett/Fligner, the design detection. `run_alpha_test()`
carries `test_reason` for exactly this.

"We used non-parametric tests" is not sufficient. "Shapiro-Wilk rejected
normality (p = 8.1e-06) and Fligner-Killeen rejected equal variances (p = 0.011),
so Kruskal-Wallis was used with Dunn's post-hoc and BH correction" is.

### 5. Effect sizes, not only p-values

Every test in this library computes one (`04`, `05`, `06`). A p-value says
whether an effect is distinguishable from zero; only the effect size says whether
it matters.

Include confidence or credible intervals wherever the method provides them.

### 6. Multiple-testing correction and its denominator

The method (BH, holm, BY) **and how many tests it corrected across**. The
denominator is as important as the method:

- DESeq2's independent filtering sets `padj` to NA for low-count taxa. "31 of 41
  tested" and "31 of 49 in the table" are different claims.
- ANCOM-BC2's `prv_cut` and `lib_cut` remove taxa and samples before testing.
- A prevalence filter applied in QC changes it too.

### 7. Software versions

Read from the session, never typed by hand. `capture_provenance()` collects them
and `provenance_paragraph()` writes the sentence. Include R itself, the platform,
and every analysis package version.

### 8. Random seeds

Rarefaction and permutation tests are stochastic. Without the seed a rerun
disagrees with the report and neither is wrong. `capture_provenance()` recovers
every seed from the run log rather than from the caller's memory of what it
passed.

### 9. The run identifier and output location

So a reader — or you in six months — can find the artifacts behind every number.

## Data characteristics that must be disclosed

`validate_inputs()` produces these as warnings, and they carry into the report
automatically. They are not caveats to bury; each one changes how a result should
be read:

- depth range and fold variation
- zero proportion
- group imbalance ratio
- **repeated measures**, and how they were handled
- tree provenance, if phylogenetic metrics were used
- heterogeneous dispersion, if PERMANOVA was run
- an anti-conservative p-value distribution, if DESeq2 was run

## Language

| Do not write | Write |
|---|---|
| "no difference" | "no statistically significant difference" |
| "trend toward significance" | the p-value and the effect size, and let the reader judge |
| "X increased Y" (observational data) | "X was associated with higher Y" |
| "community composition differed" (heterogeneous dispersion) | "the groups differed in their distance-matrix structure; dispersion was also heterogeneous" |
| "the core microbiome comprised N taxa" | "at a detection threshold of D in at least P% of samples, N taxa met the criterion" |
| "significant after correction" | the q-value and the correction method |
| "we normalized the data" | the specific method, and which analysis it fed |

Causal language requires a design that supports it. Microbiome association
studies almost never do.

## Reporting a null result

"No significant difference" is weak evidence of no difference unless the study had
power to detect one. Before writing it, check `12-power-and-sample-size.md` and
state the effect size the design could have detected.

An underpowered null is not a finding of no effect; it is a finding of no
information.

## Reporting disagreement between methods

When several DA methods were run (`06`), report the **consensus** as the primary
finding and the method-specific sets as method-dependent. Do not pick the method
with the most hits and present it as the result.

A taxon found by one method only is a hypothesis. A taxon found by all is a
result. Saying which is which is more useful than a longer list.

## The report structure

`render_report()` produces this skeleton, with the provenance, validation and QC
sections generated from the objects so they cannot drift from what ran:

```
1. Provenance          run id, versions, seeds, data dimensions
2. Data and QC         depth, sparsity, group sizes, warnings,
                       filtering table, exclusion list, normalization
3. <analysis sections> supplied by the agent as markdown
4. Figures             every logged figure
5. Saved tables        every logged table
6. Reproducibility     run id, output path, pointer to run_log.csv
7. Incomplete reporting  anything mandatory that is missing
```

Section 7 appears only when something is missing, and its presence in a finished
report is a failure, not a disclosure.
