# Power and sample size

From Xia, Sun & Chen (2018) Ch. 5. Open this before collecting data, or when
asked whether a null result means anything.

## Why it belongs in a downstream-analysis reference

Two questions arrive constantly and both need it:

1. "We found no significant difference — does that mean there is no effect?"
2. "Is n = 12 per group enough?"

Neither can be answered from the p-value.

## The shape of the problem

Microbiome power calculation is harder than the textbook case for three reasons:

- **The outcome is not one number.** Power for alpha diversity, for a PERMANOVA,
  and for per-taxon differential abundance are three different calculations with
  three different answers.
- **Overdispersion.** Counts vary far more than Poisson, so effective sample size
  is well below nominal.
- **Multiple testing.** Per-taxon power must be computed at the corrected
  threshold, not at 0.05. Testing 200 taxa at BH-0.05 means individual taxa face
  an effective threshold far below 0.05.

Consequence: a design adequately powered for alpha diversity is often badly
underpowered for differential abundance of individual taxa. Those are not the
same study.

## Alpha diversity

The standard two-sample calculation applies, on the diversity index rather than
on counts.

```r
# Parametric, for a normally distributed index such as Shannon
power.t.test(n = 20, delta = 0.5, sd = 1.0, sig.level = 0.05)
power.t.test(power = 0.8, delta = 0.5, sd = 1.0, sig.level = 0.05)  # solve for n
```

`delta` and `sd` must come from pilot data or a comparable published study, in the
units of the index. A Shannon difference of 0.5 with SD 1.0 is a moderate effect.

**Correction for the non-parametric test.** Wilcoxon has roughly 95% of the
efficiency of the t-test under normality and can exceed it under skew. Inflating
the parametric n by ~15% is a conservative rule of thumb when the
non-parametric branch will be used.

**Repeated measures change this.** With several samples per subject, effective n
is closer to the number of **subjects** than the number of samples, scaled by the
intra-class correlation:

```
n_effective ≈ n_samples / (1 + (m - 1) * ICC)
```

where `m` is samples per subject. With m = 5 and ICC = 0.5, 100 samples from 20
subjects carry the information of about 33 independent samples. `fit_lmm()`
reports the ICC, so this can be computed after the fact.

## Beta diversity / PERMANOVA

No closed form. The practical approach is simulation, or reasoning from R²:

- PERMANOVA power depends on R², group sizes, and the number of permutations.
- With 999 permutations the smallest achievable p-value is 0.001, so tiny n cannot
  produce small p-values however large the effect.
- **The floor:** with n = 3 per group and two groups there are only 10 distinct
  permutations of the group labels, so the minimum attainable p-value is 0.1. The
  test cannot reject at 0.05 regardless of effect size. This is worth checking
  before running it.

Rough guide from the literature: detecting R² ≈ 0.05 (a small but real community
effect) needs roughly 20–30 per group; R² ≈ 0.15 needs roughly 10–15.

## Differential abundance

The hardest case, and the one where intuition is most wrong.

Power per taxon depends on its mean abundance, its dispersion, the fold change,
sequencing depth, **and** the multiple-testing burden. A taxon at 0.01% relative
abundance in a 5,000-read sample is expected to appear 0.5 times; no method can
detect a two-fold change in it at any n.

Practical implications:

- **Prevalence filtering is a power intervention**, not just tidying. Removing
  taxa that cannot be detected raises power on those that can (`02`).
- **Agglomerating to genus or family raises power** by reducing both the number of
  tests and the sparsity per feature. 143 genera is a much better-powered analysis
  than 619 ASVs, at the cost of resolution.
- **Sequencing depth and sample size are not interchangeable.** Depth helps rare
  taxa; sample size helps everything. For a fixed budget, more samples at moderate
  depth generally beats fewer at high depth for group comparison.

`DESeq2` has no built-in power function for this; the usual routes are simulation
from a fitted negative binomial, or published power curves for the chosen method.

## Reporting a null result

This is the most frequent use of this file.

**An underpowered null is not a finding of no effect; it is a finding of no
information.** The report must say which it is.

Instead of "there was no significant difference in Shannon diversity", write:

> Shannon diversity did not differ significantly between groups (Wilcoxon
> p = 0.42, rank-biserial r = 0.11). With n = 12 and 14 and an observed SD of
> 0.9, this design had 80% power to detect a difference of 1.05 Shannon units;
> the observed difference was 0.18. A difference smaller than about one Shannon
> unit cannot be excluded.

That is a bounded negative result. The bare version is not.

The same applies to differential abundance: "no taxa were significant after FDR
correction" should be accompanied by the abundance and fold change the design
could have detected.

## Post-hoc power

Computing power from the observed effect size is circular — it is a monotone
function of the p-value and adds nothing. Do not do it.

The useful version is the **minimum detectable effect**: given the realised n and
variance, what effect size would have reached significance? That is a genuine
bound on the negative result, and it is what the wording above reports.

## When a design cannot answer the question

Some designs are not underpowered, they are unidentifiable, and no n fixes them:

- batch fully confounded with the variable of interest (`11`, item 14)
- one sample per group
- pseudo-replication where the true n is 2 (`11`, item 15)

The honest report says the effect is not estimable from this design. That is more
useful than a p-value computed anyway.
