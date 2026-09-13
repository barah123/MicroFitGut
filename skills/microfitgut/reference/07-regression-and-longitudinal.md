# Regression, mixed models, and longitudinal designs

Where repeated measures are handled properly rather than warned about.

## The governing distinction

A fixed-effects model (`lm`, `glm`) assumes every observation is an independent
draw. With several samples per subject that is false: within-subject samples are
correlated, the standard errors are underestimated, and the result is false
positives.

A mixed model absorbs subject-level variation into a random effect. It often
removes an "effect" that was really between-subject variation all along.

Demo 12 shows this directly: `Moraxella_log` is a strong predictor of cytokine
level under `lm`, and has **no apparent effect** once `(1 | patient)` is added.
The lm result was not a weaker version of the truth; it was an artefact of
treating 32 observations per patient as 32 independent samples.

`fit_lm()` and `fit_glm()` detect repeated measures and warn. They do not refuse —
sometimes a fixed-effects model is the right comparison to show — but the warning
must reach the report.

## The ladder

| Model | Response | Handles | Function |
|---|---|---|---|
| LM | continuous | fixed effects only | `fit_lm()` |
| GLM | continuous, binary, count | non-Gaussian families | `fit_glm()` |
| LMM | continuous | + subject random effects | `fit_lmm()` |
| GLMM | any | + random effects, overdispersion | `fit_glmm()` |
| ZIP / ZINB | count | excess zeros | `fit_zip()`, `fit_zinb()` |
| ZHP / ZHNB | count | excess zeros, hurdle story | `fit_zhp()`, `fit_zhnb()` |
| Bayesian joint | **several at once** | multiple responses + random effects | `fit_brms_joint()` |

Climb only as far as the design requires. A more complex model that does not
converge is worse than a simpler one that does.

## Offsets

```r
df <- add_offset(df, depth_col = "total_reads", offset_name = "Offset")
```

Modelling raw counts without an offset asks "does this taxon have more reads in
group A", which is partly a question about sequencing depth. The offset
`log(total reads)` turns it into a question about **rate**, which is the
biological question.

Offsets belong in the **count component only**, not the zero component (Demo 11).

## Random effects

`random = "(1 | patient)"` — random intercepts. Each patient has their own
baseline level of the response. This handles within-subject correlation and is
what most longitudinal microbiome designs need.

`random = "(1 + time | patient)"` — random intercepts **and slopes**. Each patient
also has their own rate of change over time.

Random slopes cost degrees of freedom and frequently fail on microbiome-sized
data. **A singular fit is the model telling you there is not enough information
to estimate them** — a variance at or near zero, or a perfect correlation between
intercept and slope. Drop the slopes and refit.

`fit_lmm()` reports the **ICC**: the share of total variance sitting between
subjects. A high ICC is the quantitative reason the mixed model was necessary. An
ICC of zero with random intercepts only means the grouping explains no variation
beyond the fixed effects, and a plain `lm` would give the same answer — which is
worth reporting rather than presenting a mixed model with a null variance
component.

## p-values

`lme4::lmer` withholds p-values by design (the denominator degrees of freedom are
not well defined). Two ways round it, both used here:

- `lmerTest::lmer` — Satterthwaite approximation. `fit_lmm()` uses this, so
  `anova(fit)` gives F tests with p-values.
- `car::Anova(fit, type = "II")` — Wald chi-square tests. `fit_glmm()` attaches
  this, because `glmer` and `glmmTMB` give no p-values in their summaries
  (Demo 12 does exactly this).

## Model comparison — the rules that are easy to get wrong

`compare_models(..., what = "fixed" | "random")` enforces them.

### Comparing fixed-effect structures requires ML, not REML

REML likelihoods computed under different fixed effects are **not comparable**.
Refit with `REML = FALSE` and compare those. `compare_models(what = "fixed")`
**refuses** to compare REML fits rather than returning a misleading LRT.

### LRT is unreliable for comparing random-effect structures

The null hypothesis sits on the boundary of the parameter space (a variance of
zero), so the p-value is conservative and can point the wrong way. Demo 12 notes
this twice.

Prefer AIC/BIC and parsimony, and **prefer the simpler structure when they
disagree**. Demo 12's AIM IV is exactly this case: the LRT favours random slopes
(p = 0.041) while AIC and BIC favour intercepts only, and the conclusion is that
intercepts only is the better model.

### Delta AIC below 2

Two models within 2 AIC units are not meaningfully distinguishable. Prefer the
simpler or more interpretable one; do not claim the lower-AIC model is better.
`compare_models()` and `compare_zero_models()` both say so when it applies.

### Keeping a non-significant predictor

Demo 12's answer, which generalises: it depends on the goal. For **prediction**,
keep it — it may explain variance in combination with others. For
**interpretation**, drop it — that simplifies the model and avoids overfitting.
State which goal you are serving.

## Zero-inflated and hurdle models

Two families telling two different stories about the zeros.

**Zero-inflated** — zeros come from two processes mixed together: structural zeros
(the taxon is genuinely absent) and sampling zeros (present but not sequenced
deeply enough to see). The count component can itself emit zeros.

**Zero-hurdle** — one process decides presence/absence, a second models abundance
given presence. All zeros are structural; the count part is truncated at zero.

**This is a biological claim, not a statistical convenience.** Hurdle says every
zero means absent. Zero-inflated says some zeros are undersampling. For microbiome
data at realistic depth, zero-inflated is usually the more honest story.

### The formula

```r
f <- zero_model_formula("nReadsLv", count_predictors = "pregnancy",
                        zero_predictors = "pregnancy", offset_name = "Offset")
# nReadsLv ~ pregnancy + offset(Offset) | pregnancy
```

Before `|` is the count process; after it is the probability of a structural zero.
They can carry different predictors.

### Diagnose before choosing

```r
zero_diagnostics(counts, "Lactobacillus vaginalis")
```

On `Lvaginalis.csv`: 73% zeros (660 of 900), a Poisson with mean 2.66 predicts 63
zeros, so there are 597 **excess** zeros; variance/mean = 30.5. Both
zero-inflated and overdispersed → ZINB or ZHNB.

Without this check the temptation is to fit a zero-inflated model whenever zeros
look numerous. Overdispersion alone, without excess zeros, needs a plain negative
binomial and a zero component is not justified.

### Compare them correctly

```r
compare_zero_models(f, data)
```

LRT applies **only to nested pairs** — ZIP within ZINB, ZHP within ZHNB. ZIP and
ZHP are not nested in each other, so AIC compares all four while the LRT is
restricted to the two valid pairs (Demo 11).

On `Lvaginalis.csv` the result reproduces Demo 11 exactly: ZINB 2457.6,
ZHNB 2467.7, ZIP 4576.4, ZHP 4577.1. ZINB wins, and both LRTs strongly favour the
negative binomial over the Poisson.

### Interpreting a two-component model

From Demo 11, and worth quoting because the agent will not infer it reliably:

| Count model | Zero model | Reading |
|---|---|---|
| significant | significant | The predictor affects both occurrence and abundance — a dual role. |
| significant | not significant | Affects abundance but not occurrence: meaningful **only if** the taxon is present. |
| not significant | significant | Affects only whether zeros occur, not abundance when present. |
| not significant | not significant | No effect on either. |

Signs are relative to the reference level. In Demo 11, pregnancy has a positive
count coefficient (0.81) and a negative zero-inflation coefficient (−1.36):
pregnancy raises expected abundance when the taxon is present **and** reduces the
chance it is structurally absent.

## Bayesian joint models

LMM and GLMM fit one response at a time. When two outcomes are measured on the
same subjects — health status and a cytokine, say — modelling them jointly borrows
strength across them and estimates their residual correlation. `mvbind()` binds
the responses and each can have its own family (Demo 13).

```r
fit <- fit_brms_joint(bf(mvbind(status, cyt) ~ Moraxella + time + (1 | patient)),
                      data = df, family = NULL, chains = 4, iter = 2000)
```

**Check convergence before reading anything.** `fit_brms_joint()` warns when
max Rhat ≥ 1.01, the minimum effective-sample-size ratio ≤ 0.1, or there are any
divergent transitions. Those mean the chains have not explored the posterior and
the intervals are not trustworthy. Raise `iter`, raise `adapt_delta`, or simplify.

### Interpretation — the rule from Demo 12/13

```r
brms_hypothesis(fit, "Moraxella > 0")
```

**If the 95% credible interval excludes 0, the model is confident the effect is
nonzero and the null is rejected in Bayesian terms.** brms marks these with a
star. `Evid.Ratio` is the posterior odds in favour of the hypothesis.

This is not a p-value and should not be described as one. It is a statement about
the posterior, conditional on the model and priors.

### LOO-CV model comparison

```r
compare_loo(full = fit, reduced = fit_r)
```

**The rule: `elpd_diff / se_diff > 2` is strong evidence for the model with higher
ELPD.** The first row is always the best model, so its `elpd_diff` is 0.

Demo 13: elpd_diff = 166.9, se_diff = 16.3, ratio ≈ 10.2 → very strong evidence
for the full model. quiz13: elpd_diff = −2.3, se_diff = 3.1, ratio = 0.74 → **not**
strong evidence, so the models are not distinguishable and the simpler one is
preferred. `compare_loo()` computes the ratio and states which case applies.

## Per-taxon modelling

```r
df <- taxon_frame(ps, taxon = "ASV1")
```

One row per sample with the taxon's count, sample depth, a log offset, and every
metadata column — the shape all the fitters take.

Running a model per taxon across hundreds of taxa is a multiple-testing problem
and needs correction across taxa, which is what the methods in `06` do properly.
Per-taxon modelling here is for a handful of pre-specified taxa, not a screen.
