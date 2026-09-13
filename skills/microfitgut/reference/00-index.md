# Reference index

One line per file saying when to open it. Route here first, then open only the
files the task needs. Opening everything wastes context and buries the decision
you actually have to make.

| File | Open it when |
|---|---|
| `01-intake-and-validation.md` | Loading any dataset. What the accepted inputs are and the checks that must pass before anything runs. |
| `02-qc-and-filtering.md` | Deciding what to remove and justifying the thresholds. Also: what must be logged when samples are dropped. |
| `03-normalization.md` | **Always, before any comparative statistic.** The decision table mapping normalization to the analyses it is valid for. The single highest-value file here. |
| `04-alpha-diversity.md` | Computing or testing within-sample diversity. Which metrics, the assumption-testing sequence, which test, which post-hoc. |
| `05-beta-diversity.md` | Distances, ordination, PERMANOVA/ANOSIM/MRPP. Includes the dispersion check that qualifies every PERMANOVA conclusion. |
| `06-differential-abundance.md` | Asking which taxa differ. Decision table mapping study design to method, and what each method's p-value means. |
| `07-regression-and-longitudinal.md` | Repeated measures, covariates, zero-inflated counts, multiple responses. LM → GLM → LMM → GLMM → zero models → Bayesian joint. |
| `08-core-and-composition.md` | Core microbiome, agglomeration, relative abundance summaries, dominance. |
| `09-plot-conventions.md` | Making any figure. Palette, dimensions, what must appear on each figure type. |
| `10-reporting-standards.md` | Writing any summary. What every report must state regardless of what was run. |
| `11-pitfalls.md` | **Read once fully, then revisit per stage.** The expert-knowledge file: compositionality, rarefaction disputes, pseudo-replication, name instability, the dispersion confound. |
| `12-power-and-sample-size.md` | Before collecting data, or when asked whether a null result means anything. |
| `13-benchmark-study-card.md` | Benchmark mode. Extracting what a paper says and — more importantly — what it does not. |
| `14-benchmark-concordance.md` | Benchmark mode. Which agreement metrics to compute per domain and how to read them. |
| `15-benchmark-discrepancy-taxonomy.md` | Benchmark mode. Known divergence causes paired with their diagnostic signatures. |

## The shortest useful path

For a standard two-group comparison on 16S counts:

1. `01` → load and validate. Stop if validation fails.
2. `11` → skim the pitfalls that apply to this design (sparsity, repeated measures, group imbalance).
3. `02` → set and justify filters.
4. `03` → **pick the normalization per downstream analysis, not once for the run.** Different stages need different states of the same data.
5. `04`, `05`, `06` → the analyses themselves.
6. `09` → figures.
7. `10` → the report, then the traceability check, then the `verifier` subagent.

For a design with repeated measures, `07` replaces most of `04`–`06`, because
the simple tests in those files are not valid and the scripts will refuse them.

## The rule that overrides everything here

A conclusion that cannot be traced to a run identifier and a stored artifact does
not get reported. If a number is not in `run_log.csv` or a file in `tables/`, it
does not go in the report. `mfg_traceability_check()` is the mechanical check;
the `verifier` subagent is the judgement.
