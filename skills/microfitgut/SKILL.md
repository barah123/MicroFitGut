---
name: microfitgut
description: >-
  Downstream microbiome data analysis agent for 16S amplicon and shotgun
  metagenomics. Use when analysing a phyloseq object (.RDS), an ASV/OTU table
  with taxonomy and metadata, a MetaPhlAn or Kraken2/Bracken taxonomic profile,
  or a PICRUSt2/HUMAnN functional profile — covering QC and filtering,
  normalization, alpha diversity (Shannon, Chao1, Faith's PD), beta diversity
  (Bray-Curtis, UniFrac, PCoA, NMDS, PERMANOVA), differential abundance
  (ANCOM-BC2, ALDEx2, DESeq2), mixed and zero-inflated models for longitudinal
  designs, core microbiome, publication figures, and a reproducible report. Also
  use for benchmark mode: reanalysing a published dataset, measuring concordance
  with the paper's reported results, and attributing any divergence to a specific
  analytical decision. Triggers on: analyse microbiome data, 16S analysis,
  metagenomics, phyloseq, alpha diversity, beta diversity, PERMANOVA,
  differential abundance, ANCOM-BC, ALDEx2, microbiome report, reproduce this
  paper's microbiome analysis, benchmark microbiome pipeline.
---

# MicroFitGut

Downstream microbiome analysis. The data is already sequenced, denoised and
classified; this takes it from a feature table to plots, statistics and a
reproducible report.

**The model orchestrates and interprets. The R functions in `scripts/` compute.**
Do not write ad-hoc analysis code when a script function exists — every function
you bypass is a chance to produce something subtly wrong, and an unlogged
computation cannot be traced or reported.

## Two modes

| Mode | When | Path |
|---|---|---|
| **analyse** | A dataset and a question. | Stages 1–11 below. |
| **benchmark** | A dataset whose published result is known; measure whether it reproduces and explain any difference. | Stages 1–11, then B1–B5. |

Benchmark mode depends on analyse mode being right, so the analysis stages run
unchanged and the benchmark layer sits on top.

## Before anything

```r
# Resolve the script directory. Works whether MicroFitGut is installed as a
# plugin, symlinked into a project, or run from a checkout — the location is
# found, never assumed, so a wrong path fails here with a clear message rather
# than nine stages later.
S <- local({
  root <- Sys.getenv("CLAUDE_PLUGIN_ROOT", "")
  cands <- c(
    if (nzchar(root)) file.path(root, "skills", "microfitgut", "scripts"),
    ".claude/skills/microfitgut/scripts",
    "skills/microfitgut/scripts",
    "microfitgut-plugin/skills/microfitgut/scripts"
  )
  hit <- cands[file.exists(file.path(cands, "00-packages.R"))]
  if (!length(hit)) {
    stop("Cannot locate the MicroFitGut scripts. Looked in:\n  ",
         paste(cands, collapse = "\n  "),
         "\nSet MICROFITGUT_SCRIPTS or pass the path explicitly.", call. = FALSE)
  }
  hit[1]
})
source(file.path(S, "00-packages.R"))
mfg_print_package_status()        # tell the user up front if a stage is blocked
mfg_require(c("intake", "alpha", "beta", "da"))   # load what this run needs
for (f in c("utils.R", "01-intake.R", "02-qc.R", "03-normalize.R", "04-alpha.R",
            "05-beta.R", "06-differential-abundance.R", "07-models.R",
            "08-core-composition.R", "09-exploratory.R", "10-plots.R",
            "11-report.R")) source(file.path(S, f))

mfg_start_run("<short-label>", outdir = "output")
```

`mfg_start_run()` fixes the run identifier and creates
`output/<run_id>/{figures,tables}`. Everything downstream writes there and logs.

Read `reference/00-index.md` and route to the files the task needs. Do not open
all of them.

---

## The eleven stages

### 1. Intake — `reference/01`

**Inventory before loading.** A data directory rarely holds one dataset.

```r
inv <- mfg_input_inventory("data"); print(inv)
```

If it reports more than one self-contained dataset, say which one you are
analysing and why, in the response, and record it:

```r
mfg_declare_authoritative("data/study.RDS",
  reason = "the full study: 344 samples with tree and 12 metadata variables; ps9/ps10/ps12 are teaching subsets")
```

The reason is required, it is quoted in the methods, and `assemble_summary()`
lists `authoritative_source` as missing content if several datasets were present
and none was declared. Two files of identical size with different names are the
case this exists for.

Then load. Every loader records the file's path, format, size, modification time
and MD5 checksum, so the report can name the exact file version analysed.

```r
ps <- mfg_load("data/study.RDS")
# or build_phyloseq(abund, tax, meta, tree, seqs, tree_is_real = TRUE)
# or build_phyloseq_from_profile("metaphlan.tsv", meta)       # shotgun
# or build_phyloseq_functional("pathways.tsv", "pwtax.txt", meta)  # PICRUSt2
```

Then **always**:

```r
v <- validate_inputs(ps, group_var = "<group>", subject_var = "<subject>")
print(v)
```

It stops on anything that would invalidate a result. `v$warnings` must reach the
report — they are not optional caveats.

**Two things to settle here and not later:**

- **Tree provenance.** `mfg_mark_tree_real(ps, TRUE/FALSE)`. Until marked,
  UniFrac, Faith's PD and phylogeny plots all refuse. A random tree merged in to
  fill phyloseq's slot gives numbers that look valid and encode nothing.
- **Factor levels.** A character column's reference level is alphabetical, and
  every coefficient and fold change is relative to it. Set it deliberately.

If the object came from outside the pipeline, assert its normalization:
`ps <- mfg_set_normalization(ps, "raw")`.

### 2. Plan — state it before running it

Say, in the response: the design, the question, the tests you will run and why,
and the thresholds you will use. Fix filters and the primary DA method **now**.
Choosing them after seeing p-values invalidates the FDR (`reference/11`, item 7).

If the design has repeated measures — `validate_inputs()` will have said so —
the plan changes: `reference/07` replaces most of the simple tests, and the
scripts in stages 4–6 will refuse them.

### 3. QC — `reference/02`

```r
qc <- run_qc(ps, remove_unwanted = TRUE, prevalence = 0.10,
             min_depth = 5000,
             depth_justification = "<why this threshold, from the curves>")
ps_f <- qc$ps
```

`min_depth` has no default. Choose it from `choose_rarefaction_depth()` and the
rarefaction curves, and record why.

### 4. Normalize — `reference/03`, the highest-stakes file

**Per analysis, not once per run.** A typical run holds three states at once:

```r
ps_rare <- rarefy_to_depth(ps_f, depth = 5000, rngseed = 42, justification = "...")
ps_tss  <- tss_transform(ps_f)
ps_raw  <- mfg_set_normalization(ps_f, "raw")
```

`mfg_check_normalization()` runs inside every analysis function and **refuses**
invalid pairings. If it stops you, it is right — read the reason, do not work
around it.

### 5. Alpha diversity — `reference/04`

```r
ad <- alpha_table(ps_rare, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
a  <- check_assumptions(ad, "Shannon", "<group>")
r  <- alpha_test_all(ad, "<group>")
d  <- posthoc_dunn(ad, "Shannon", "<group>", omnibus_p = r$results$Shannon$p_value)
```

### 6. Beta diversity — `reference/05`

```r
dist <- compute_distance(ps_tss, "bray")
ord  <- run_ordination(ps_tss, "PCoA", dist_obj = dist)
pm   <- run_permanova(ps_tss, ~ <group>, dist_obj = dist,
                      permutations = 999, strata = "<subject or NULL>")
```

`run_permanova()` runs the betadisper check automatically. **Use
`permanova_sentence(pm)` to write the conclusion** — it refuses to say
"composition differed" when dispersion is heterogeneous, which is the qualification
most often dropped.

### 7. Differential abundance — `reference/06`

```r
choice <- choose_da_method(ps_raw, "<group>"); print(choice)
ab <- da_ancombc2(ps_raw, fix_formula = "<group>", group = "<group>",
                  rand_formula = NULL, prv_cut = 0.10, struc_zero = TRUE)
al <- da_aldex2(ps_raw, "<group>", mc.samples = 128)
de <- da_deseq2(ps_raw, ~ <group>); da_deseq2_diagnostics(de)
print(da_compare(ancombc2 = ab, aldex2 = al, deseq2 = de))
```

Run at least two methods. Report the **consensus** as the primary finding and the
method-specific sets as method-dependent.

### 8. Models — `reference/07`

For repeated measures, covariates, several responses, or excess zeros. Climb only
as far as the design requires; a complex model that does not converge is worse
than a simple one that does.

### 9. Composition and exploration — `reference/08`, `reference/09`

`core_taxa()` with **both thresholds stated**, `core_threshold_scan()` to show the
surface, `taxa_summary()`, `dominance_summary()`.

Exploratory output (heat trees, heatmaps, networks) is **hypothesis, not result**.
Say which is which in the report.

### 10. Figures — `reference/09`

```r
p <- plot_alpha(ad, "Shannon", "<group>", test_result = r$results$Shannon)
mfg_save_figure(p, "alpha_shannon", size = "double_column")
```

Every figure through `mfg_save_figure()` or `mfg_save_plot()`, which log the
file. An unlogged figure cannot be traced, and under the rule below it does not
get discussed.

### 11. Report — `reference/10`

```r
prov <- capture_provenance(ps_initial = ps, ps_final = ps_f, qc_log = qc$log)
sm   <- assemble_summary(prov, validation = v, qc_log = qc$log,
                         sections = list(...), group_var = "<group>")
out  <- render_report(sm, sections_md = list("Alpha diversity" = "...", ...))
mfg_manifest()
```

`assemble_summary()` lists any mandatory content that is missing. Section 7 of the
rendered report is "Incomplete reporting" — its presence in a finished report is a
failure, not a disclosure.

`mfg_manifest()` also writes `analysis_calls.R`: every MicroFitGut call the run
made, in order, with the arguments as they were passed, reconstructed from the
log rather than transcribed. It is a record, not a runnable script — assignments
are not captured, because R cannot see from inside a function what its result was
bound to. Read it before the report is final. A call you expected and cannot find
there is a step that did not log, and a step that did not log cannot be reported.

---

## Benchmark mode

### B1. Study card — `reference/13`

Extract what the paper says, what it does not, and where each value came from.
`confidence = "absent"` is a first-class value: an unstated parameter is a
benchmark finding and automatically becomes a sweep axis.

**Present the card to the user and get confirmation before running the sweep.**
Extraction is the one step no artifact can validate.

### B2. Replication run

Run stages 1–11 at the card's baseline — the best reconstruction of the paper.

### B3. Harmonise and measure — `reference/14`

```r
al <- align_abundance_tables(published_mat, reanalysis_mat, ptax, rtax)
ca <- concordance_alpha(pub_shannon, rea_shannon, group)
cb <- concordance_beta(dist_pub, dist_rea, pm_pub, pm_rea)
cd <- concordance_da(pub_taxa, rea_taxa, pub_eff, rea_eff,
                     all_tested_reanalysis = rownames(de$results))
verdict <- concordance_verdict(alpha = ca, beta = cb, da = cd)
```

Harmonise before comparing. Skipping it **manufactures** disagreement.

### B4. Sweep — `reference/13`

```r
sg <- build_sweep_grid(card, ps, max_runs = 20); print(sg)
```

Run each configuration, collect the concordance metrics into one data frame with
`run_id`, `axis`, `level` and the metric columns.

### B5. Attribute — `reference/15`

```r
att <- attribute_divergence(results)          # computed ranking, not asserted
rec <- attribution_record(verdict$verdict, causes = list(...),
                          unexplained_residual = "...",
                          confidence = "medium", attribution = att, sweep_grid = sg)
```

Every cause needs a `rubric_id`, the `axis` that demonstrated it, and `evidence`
naming run ids and metric movement. Where the attribution confidence reads
"weak", **name no single cause**.

---

## Subagents

Both ship with this plugin. Installed, they are addressed as
`microfitgut:statistician` and `microfitgut:verifier`; in a local checkout with
the agents copied to `~/.claude/agents/` the bare names work. Use whichever form
the agent listing shows — if neither resolves, say so rather than skipping the
step, because both exist to catch what the scripts cannot.

### `microfitgut:statistician` — before the report is written

Delegate whenever:

- the design is not covered by the reference files (nested random effects,
  crossover, a method not in the demos)
- several plausible model specifications exist and the choice changes the answer
- the statistical section is about to be finalised

It audits design-to-method match, multiple-testing correction, whether numbers
trace to output, and effect-size reporting. It does **not** run new analyses, so
give it the scripts and the output files, and carry its findings into the report
rather than quietly acting on them.

### `microfitgut:verifier` — before the report is final

Always, in analyse mode and benchmark mode both. Run the mechanical check first:

```r
tc <- mfg_traceability_check(report_text); print(tc)
```

then hand the flagged list, `run_log.csv`, `manifest.csv` and the draft to the
verifier. It checks every number against the artifacts, that the methods text
describes what the code actually did, and every citation.

---

## Standing rules

1. **Never report a conclusion that is not traceable to a run identifier and a
   stored artifact.** If a number is not in `run_log.csv` or a file under
   `tables/`, it does not go in the report.
2. **State assumptions and intermediate numbers, not just the final plot.** n per
   group after filtering, every sample dropped and why, the normalization, the
   test with its assumption checks, effect sizes with intervals.
3. **Never drop a sample or a taxon without logging it.** `run_qc()` does this;
   do not filter around it.
4. **A guard that refuses is right.** The normalization pairing check, the tree
   provenance gate, the repeated-measures refusal and the REML comparison block
   exist because each catches an error that otherwise produces a clean, wrong
   result. Read the reason and change the analysis, not the guard.
5. **Fix thresholds before testing.** Choosing a filter or a method after seeing
   p-values invalidates the FDR.
6. **No silent caps.** If coverage was bounded — top-N taxa, a truncated sweep, a
   subsampled figure — say what was left out.
7. **Uncertain and uncovered decisions stop and ask**, after completing everything
   that does not depend on the answer. Consult `microfitgut:statistician` first,
   then put
   the question to the user with the options and their consequences.
8. **Report failures faithfully.** If a model did not converge, a method errored,
   or a stage was skipped, say so with the output. A run with three of five
   analyses complete is reported as three of five.
