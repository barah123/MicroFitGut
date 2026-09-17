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
before any corpus study was scored.

**Correction on the direction of those defects.** An earlier version of this file
said every one of them biased toward reporting non-reproduction. That claim does
not hold for all twenty, and the reports themselves show why. Defects #14 and #16
are both failures to detect a clustering variable, and an undetected cluster
means the naive unstratified test runs, which inflates significance. That is the
anti-conservative direction, and it makes a claim of difference *more* likely to
be scored as reproduced. The Ramos report says so in as many words. It happened
to work the other way in that study only because the claim there was a null one.

The defensible statement is narrower: every defect was found and fixed before any
corpus study was scored, and the direction of bias must be recorded per defect
rather than asserted in aggregate. That per-defect register does not yet exist
and is listed as outstanding work.

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

## What can and cannot be traced

**Correction.** An earlier version of this file stated that no run artifacts
survived. That was wrong, and it understated the evidence. What actually exists:

| | |
|---|---|
| Run directories under `output/` for the ten studies | 48 |
| Carrying a `run_log.csv` | 39, holding 293 logged rows |
| Carrying `tables/*.csv` | **0** |
| Carrying `figures/` | **0** |

So a value that was written to a log can be traced to the run that produced it.
A value that lived only in a result table cannot, because no result table was
saved. Under the project's first standing rule that gap still has to be closed
before the main study, but it is narrower than previously stated.

Two further traceability gaps the logs expose:

1. **The claim counts in these reports do not all match the runs.** Summing the
   final per-study `claims_scored` events gives **53 claims**; the reports
   describe **57**. The four extra are in A2 and A4, where claims were settled in
   the report prose without being passed through `score_claims()`. They are real
   claims the papers made, but they carry no logged adjudication.
2. **The encoded claim objects were not saved.** The `type`, `comparator`,
   `tolerance` and `contrast` used for each claim exist only as report prose.
   A0's `ten-genera` verdict, for instance, rests on a tolerance of 3 against a
   claim of 10, a window of 7 to 13, which is recoverable only from the text.

## Known gap

Two verdicts changed after the scoring code was changed in response to seeing a
result, and both should be read with that in mind:

- **S1 Vogtmann.** Two verdict runs 40 seconds apart, `held=6` then `held=8`,
  either side of adding the `da_null` claim type. The change encodes a kind of
  claim the schema could not previously express, which is a defensible reason,
  but the sequence was score, inspect, change the scorer, rescore.
- **A0 Pérez-Losada.** Two verdict runs, recovery 20% then 30%, either side of a
  taxon-separator fix.

Disclosed, both are ordinary instrument development. Undisclosed, either is
something an adversarial reader would find in the logs.
