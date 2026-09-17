# Validation

Ten published studies, five 16S amplicon and five shotgun metagenomic, run
through MicroFitGut's benchmark mode. The set is a purposive coverage matrix
chosen to exercise every input path the tool claims to support. **It is not a
sample, and no rate in these reports describes the literature.**

`validation-set-manifest.csv` and `.xlsx` list all ten with paper DOI, data
repository, accession, format and design. Neither the papers nor the datasets are
redistributed here; each study is identified by DOI so it can be retrieved from
source.

## The reports predate a revision of the scoring scheme

Every `benchmark-*.md` file records the verdict produced by the scheme in force
when it was run. That scheme was audited after the tenth study finished, and the
endpoint was redefined. **The reports have not been rescored**, so their verdict
blocks use vocabulary the current code no longer emits.

| Reports say | Current code says |
|---|---|
| `held` | `reproduced` |
| `failed` | `not_reproduced` |
| `unscored` / `not_evaluable` pooled as "unadjudicated" | one of `not_licensed`, `not_tested`, `schema_gap`, `out_of_scope`, `data_absent`, `contrast_mismatch`, `not_attempted`, `pending_adjudication` |
| a four-level study verdict as the outcome | one pre-designated claim per study as the outcome; the category is descriptive only |

The revision and its reasons are in `Expert notes/claim-survival-protocol.md`
§14. Four changes matter for reading these reports:

1. **Directional claims now require significance, not just sign.** Several claims
   recorded as `held` were scored on the direction of the effect alone, at
   p-values as high as 0.955. Under the current rule they would be
   `not_reproduced`.
2. **A claim the evidence cannot license is no longer pooled with missing data.**
   The dispersion downgrades in A0, A2 and A3 are findings about those papers,
   and they now carry their own status.
3. **Concordance measures no longer decide a verdict.** A0's only recorded
   failure was a differential-abundance recovery rate falling below a hard-coded
   0.7, not one of its three claims.
4. **A taxon filtered out before testing is `not_tested`, not a failure.**

## What the pilot did and did not establish

**Established.** The tool runs across every input path in the matrix: phyloseq
objects, QIIME 2 BIOM with and without taxonomy, HDF5 BIOM with empty metadata
groups, loose delimited tables, MetaPhlAn profiles both merged and per-sample,
and curatedMetagenomicData objects. Twenty defects were found and fixed, all
before any corpus study was scored, and **every one of them biased toward
reporting non-reproduction**, which is the direction that would have flattered
the tool.

**Not established: a failure rate.** The pilot contains **no statistically
significant contradiction of any paper**. The four items recorded as failures are
artifacts of scoring convention: one is the 0.7 recovery threshold, one is the
significance-gating asymmetry between claim types, and two are direction
reversals at p = 0.465 and p = 0.662 that the reports themselves describe as
showing no difference rather than an opposite one.

**Not established: parameters for a power calculation.** Between-study ICC of the
old composite endpoint was 0.191, but the ICC of adjudicability alone was 0.199
and of reproduction conditional on adjudicability approximately zero, on ten
studies with seven at the boundary. Those numbers describe data deposition, not
analytical fragility.

## Known gap

No run artifacts survive from these ten studies: there is no `output/` directory,
`run_log.csv` or `tables/*.csv` for any of them, so no number in these reports can
be checked against the run that produced it. The project's first standing rule
requires exactly that check. Retaining one complete run directory per study,
along with the encoded claim objects, is a prerequisite for the main study.
