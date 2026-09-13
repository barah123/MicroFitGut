# =============================================================================
# MicroFitGut — 07-models.R
#
# Regression, mixed models, zero-inflated counts, and Bayesian joint models.
#
# This is where repeated measures are handled properly rather than warned about.
# The governing distinction, from Demo 12: a fixed-effects model (lm, glm) assumes
# every observation is an independent draw. With several samples per subject that
# is false, standard errors are underestimated, and the result is false positives.
# A mixed model absorbs subject-level variation into a random effect, which often
# removes an "effect" that was really between-subject variation all along.
#
# Model comparison rules, also from Demo 12 and 13:
#   - Comparing fixed-effect structures with LRT or AIC requires ML, not REML
#     (REML = FALSE), because REML likelihoods with different fixed effects are
#     not comparable.
#   - LRT is not reliable for comparing random-effect structures; prefer AIC/BIC
#     and parsimony there.
#   - For Bayesian models use LOO-CV: elpd_diff / se_diff > 2 is strong evidence.
#
# Requires: mfg_require(c("models")) and, for joint models, mfg_require("bayes")
# =============================================================================


# ── Safe accessors ───────────────────────────────────────────────────────────

#' Number of observations, for model classes that do not implement nobs().
#'
#' pscl's zeroinfl and hurdle objects have no nobs method, so stats::nobs()
#' errors on them. Falls back through the pieces that are always present.
mfg_nobs <- function(fit) {
  n <- suppressWarnings(try(stats::nobs(fit), silent = TRUE))
  if (!inherits(n, "try-error") && is.numeric(n) && length(n) == 1) return(as.integer(n))
  if (!is.null(fit$n)) return(as.integer(fit$n))
  r <- suppressWarnings(try(stats::residuals(fit), silent = TRUE))
  if (!inherits(r, "try-error")) return(length(r))
  mf <- suppressWarnings(try(stats::model.frame(fit), silent = TRUE))
  if (!inherits(mf, "try-error")) return(nrow(mf))
  NA_integer_
}

# ── Offsets ──────────────────────────────────────────────────────────────────

#' Add a log library-size offset column.
#'
#' Modelling raw counts without an offset asks "does this taxon have more reads
#' in group A", which is partly a question about sequencing depth. The offset
#' log(total reads) turns it into a question about rate, which is the biological
#' question. Offsets belong in the count component only (Demo 11).
add_offset <- function(df, depth_col, offset_name = "Offset") {
  if (!depth_col %in% names(df)) {
    stop("'", depth_col, "' is not a column in the data. The offset needs total ",
         "reads per sample.", call. = FALSE)
  }
  if (any(df[[depth_col]] <= 0, na.rm = TRUE)) {
    stop("Depth column '", depth_col, "' contains zero or negative values; ",
         "log() of those is undefined. Drop zero-depth samples first (02-qc.R).",
         call. = FALSE)
  }
  df[[offset_name]] <- log(df[[depth_col]])
  mfg_log("models", "offset_added",
          list(from = depth_col, name = offset_name,
               range = sprintf("%.2f-%.2f", min(df[[offset_name]]), max(df[[offset_name]]))))
  df
}

#' Build a per-taxon modelling frame from a phyloseq object.
#'
#' Returns one row per sample with the taxon's count, the sample depth, a log
#' offset, and every metadata column — the shape all the model fitters below take.
taxon_frame <- function(ps, taxon, depth_col = NULL) {
  mat <- mfg_otu_taxa_as_rows(ps)
  if (!taxon %in% rownames(mat)) {
    stop("Taxon '", taxon, "' not found. Available include: ",
         paste(utils::head(rownames(mat), 5), collapse = ", "), call. = FALSE)
  }
  meta <- mfg_meta(ps)
  df <- data.frame(sample = colnames(mat), count = as.numeric(mat[taxon, ]),
                   stringsAsFactors = FALSE)
  df$depth <- as.numeric(phyloseq::sample_sums(ps))[match(df$sample, phyloseq::sample_names(ps))]
  if (!is.null(depth_col) && depth_col %in% names(meta)) {
    df$depth <- meta[[depth_col]][match(df$sample, rownames(meta))]
  }
  meta$sample <- rownames(meta)
  df <- merge(df, meta, by = "sample", all.x = TRUE, sort = FALSE)
  df <- add_offset(df, "depth")
  rownames(df) <- df$sample
  attr(df, "mfg_taxon") <- taxon
  df
}

# ── Fixed-effects models ─────────────────────────────────────────────────────

#' Linear model, with the independence assumption checked rather than assumed.
#'
#' Warns when the data has repeated measures, because that is the condition under
#' which lm's standard errors are wrong and its p-values are too small.
fit_lm <- function(formula, data, subject_var = NULL, check_independence = TRUE) {
  fit <- stats::lm(formula, data = data)
  warn <- NA_character_
  if (isTRUE(check_independence)) {
    rm_det <- mfg_detect_repeated_measures(data, subject_candidates = subject_var)
    if (isTRUE(rm_det$repeated)) {
      warn <- sprintf(paste("'%s' has up to %d observations per level, so these",
        "observations are not independent. lm() underestimates the standard errors",
        "here and will report effects that a mixed model does not. Fit",
        "fit_lmm(..., random = '(1 | %s)') and compare."),
        rm_det$subject_var, rm_det$max_per_subject, rm_det$subject_var)
      warning(warn, call. = FALSE)
    }
  }
  s <- summary(fit)
  av <- stats::anova(fit)
  mfg_log("models", "lm", list(
    formula = paste(deparse(formula), collapse = ""), n = stats::nobs(fit),
    r_squared = round(s$r.squared, 4), adj_r_squared = round(s$adj.r.squared, 4),
    aic = round(stats::AIC(fit), 2),
    independence_warning = !is.na(warn)))
  structure(list(fit = fit, type = "lm", formula = formula, summary = s, anova = av,
    coefficients = s$coefficients, r_squared = s$r.squared,
    aic = stats::AIC(fit), bic = stats::BIC(fit), n = stats::nobs(fit),
    independence_warning = warn),
    class = c("mfg_model", "list"))
}

#' Generalized linear model.
#'
#' family must match the response: gaussian for a continuous outcome, binomial
#' for a two-level factor, poisson or negative binomial for counts. A binomial
#' GLM needs the response as a factor, not a character column (Demo 12).
fit_glm <- function(formula, data, family = stats::gaussian(),
                    subject_var = NULL, check_independence = TRUE) {
  fit <- stats::glm(formula, data = data, family = family)
  warn <- NA_character_
  if (isTRUE(check_independence)) {
    rm_det <- mfg_detect_repeated_measures(data, subject_candidates = subject_var)
    if (isTRUE(rm_det$repeated)) {
      warn <- sprintf(paste("'%s' indicates repeated measures; glm() treats these",
        "as independent. Use fit_glmm() with (1 | %s)."),
        rm_det$subject_var, rm_det$subject_var)
      warning(warn, call. = FALSE)
    }
  }
  s <- summary(fit)
  av <- stats::anova(fit, test = "Chisq")
  # Overdispersion: for poisson/binomial, residual deviance well above residual df
  # means the variance assumption fails and the standard errors are too small.
  disp <- if (family$family %in% c("poisson", "binomial")) {
    s$deviance / s$df.residual
  } else NA_real_
  overdispersed <- !is.na(disp) && disp > 1.5
  if (overdispersed) {
    warning(sprintf(paste("Residual deviance / df = %.2f, well above 1. The %s",
      "family's variance assumption is violated (overdispersion), so these standard",
      "errors are too small and the p-values too optimistic. Use a negative",
      "binomial (glmmTMB nbinom2) or quasi-family instead."), disp, family$family),
      call. = FALSE)
  }
  mfg_log("models", "glm", list(
    formula = paste(deparse(formula), collapse = ""), family = family$family,
    n = stats::nobs(fit), aic = round(stats::AIC(fit), 2),
    dispersion_ratio = if (is.na(disp)) NA else round(disp, 3),
    overdispersed = overdispersed))
  structure(list(fit = fit, type = "glm", formula = formula, family = family$family,
    summary = s, anova = av, coefficients = s$coefficients,
    aic = stats::AIC(fit), bic = stats::BIC(fit), n = stats::nobs(fit),
    dispersion_ratio = disp, overdispersed = overdispersed,
    independence_warning = warn),
    class = c("mfg_model", "list"))
}

# ── Mixed models ─────────────────────────────────────────────────────────────

#' Linear mixed model.
#'
#' `random` is the random-effects term as a string: "(1 | patient)" for random
#' intercepts, "(1 + time | patient)" for random intercepts and slopes.
#'
#' Random intercepts say each subject has its own baseline. Random slopes say the
#' effect of a predictor also varies by subject. Slopes cost degrees of freedom
#' and frequently fail to converge on microbiome-sized data; a singular fit is the
#' model telling you there is not enough information to estimate them (Demo 12).
#'
#' REML = TRUE is right for reporting a single model; REML = FALSE is required
#' when comparing different fixed-effect structures.
fit_lmm <- function(formula, data, random = NULL, REML = TRUE) {
  full <- if (is.null(random)) formula else
    stats::as.formula(paste(paste(deparse(formula), collapse = ""), "+", random))
  fit <- lmerTest::lmer(full, data = data, REML = REML)

  sing <- lme4::isSingular(fit, tol = 1e-4)
  conv_warn <- length(fit@optinfo$conv$lme4$messages) > 0
  if (sing) {
    has_slopes <- !is.null(random) && grepl("\\|", random) &&
                  !grepl("^\\(\\s*1\\s*\\|", trimws(random))
    warning(paste("Singular fit: a variance component was estimated at or near",
      "zero, or two are perfectly correlated — the random structure asks for more",
      "than this data identifies.",
      if (has_slopes)
        "Drop the random slopes and refit with random intercepts only."
      else paste("With random intercepts only, this means between-subject variance",
        "is effectively zero: the grouping explains no variation beyond the fixed",
        "effects, and a plain lm() would give the same answer. Report that rather",
        "than the mixed model's zero variance component.")),
      call. = FALSE)
  }
  s <- summary(fit)
  an <- stats::anova(fit)      # lmerTest gives Satterthwaite F tests and p-values
  vc <- as.data.frame(lme4::VarCorr(fit))

  # ICC: how much of the total variance sits between subjects. A high ICC is the
  # quantitative reason the mixed model was necessary.
  icc <- NA_real_
  if (nrow(vc) >= 2) {
    resid_var <- vc$vcov[is.na(vc$var1) | vc$grp == "Residual"]
    resid_var <- if (length(resid_var)) utils::tail(resid_var, 1) else NA_real_
    group_var <- sum(vc$vcov[vc$grp != "Residual" & is.na(vc$var2)], na.rm = TRUE)
    if (!is.na(resid_var) && (group_var + resid_var) > 0) {
      icc <- group_var / (group_var + resid_var)
    }
  }

  mfg_log("models", "lmm", list(
    formula = paste(deparse(full), collapse = ""), REML = REML,
    n = stats::nobs(fit), aic = round(stats::AIC(fit), 2),
    singular = sing, converged = !conv_warn,
    icc = if (is.na(icc)) NA else round(icc, 4)))
  structure(list(fit = fit, type = "lmm", formula = full, random = random,
    REML = REML, summary = s, anova = an, varcorr = vc,
    coefficients = s$coefficients, icc = icc,
    aic = stats::AIC(fit), bic = stats::BIC(fit), n = stats::nobs(fit),
    singular = sing, converged = !conv_warn),
    class = c("mfg_model", "list"))
}

#' Generalized linear mixed model.
#'
#' `engine = "glmmTMB"` is the default for counts because it offers nbinom1 and
#' nbinom2 for overdispersion, plus zero-inflation, which lme4 does not.
#' `engine = "lme4"` uses glmer.
#'
#' glmer output has no p-values by design; car::Anova() supplies Wald tests
#' (Demo 12 does exactly this).
fit_glmm <- function(formula, data, random = NULL, family = "nbinom2",
                     engine = c("glmmTMB", "lme4"), offset_col = NULL,
                     ziformula = ~ 0) {
  engine <- match.arg(engine)
  full <- if (is.null(random)) formula else
    stats::as.formula(paste(paste(deparse(formula), collapse = ""), "+", random))

  if (identical(engine, "glmmTMB")) {
    fit <- glmmTMB::glmmTMB(full, data = data, family = family,
                            ziformula = ziformula,
                            offset = if (!is.null(offset_col)) data[[offset_col]] else NULL)
    s <- summary(fit)
    # glmmTMB reports convergence = 0 even when the Hessian is not
    # positive-definite, which means at least one variance is unidentified and the
    # standard errors are not usable. Both must hold for a fit to be trusted.
    pd_hess <- isTRUE(fit$sdr$pdHess)
    conv <- isTRUE(fit$fit$convergence == 0) && pd_hess && is.finite(stats::AIC(fit))
    if (!pd_hess) {
      warning(paste("Non-positive-definite Hessian: at least one variance or",
        "dispersion parameter is not identified by this data, so the standard",
        "errors and AIC are unreliable. Simplify the random structure or the",
        "zero-inflation formula before reporting anything from this fit."),
        call. = FALSE)
    }
    an <- try(car::Anova(fit, type = "II"), silent = TRUE)
    coefs <- s$coefficients$cond
    disp_note <- if (grepl("^nbinom", family))
      "negative binomial: variance exceeds the mean, which is what count microbiome data does"
      else family
    mfg_log("models", "glmm", list(
      engine = engine, formula = paste(deparse(full), collapse = ""),
      family = family, ziformula = paste(deparse(ziformula), collapse = ""),
      n = mfg_nobs(fit), aic = round(stats::AIC(fit), 2), converged = conv,
      positive_definite_hessian = pd_hess))
    return(structure(list(fit = fit, type = "glmm", engine = engine, formula = full,
      random = random, family = family, summary = s,
      anova = if (inherits(an, "try-error")) NULL else an,
      coefficients = coefs, aic = stats::AIC(fit), bic = stats::BIC(fit),
      n = mfg_nobs(fit), converged = conv, pd_hessian = pd_hess,
      dispersion_note = disp_note),
      class = c("mfg_model", "list")))
  }

  fam <- if (is.character(family)) get(family, mode = "function")() else family
  fit <- lme4::glmer(full, data = data, family = fam,
                     offset = if (!is.null(offset_col)) data[[offset_col]] else NULL)
  s <- summary(fit)
  # glmer withholds p-values; car::Anova supplies Wald chi-square tests.
  an <- try(car::Anova(fit, type = "II"), silent = TRUE)
  sing <- lme4::isSingular(fit, tol = 1e-4)
  mfg_log("models", "glmm", list(
    engine = engine, formula = paste(deparse(full), collapse = ""),
    family = fam$family, n = stats::nobs(fit), aic = round(stats::AIC(fit), 2),
    singular = sing,
    note = "glmer reports no p-values; car::Anova Wald tests attached"))
  structure(list(fit = fit, type = "glmm", engine = engine, formula = full,
    random = random, family = fam$family, summary = s,
    anova = if (inherits(an, "try-error")) NULL else an,
    coefficients = s$coefficients, aic = stats::AIC(fit), bic = stats::BIC(fit),
    n = stats::nobs(fit), singular = sing),
    class = c("mfg_model", "list"))
}

# ── Zero-inflated and hurdle models ──────────────────────────────────────────
#
# Two families for two different stories about the zeros.
#
# Zero-INFLATED: zeros come from two processes mixed together — structural zeros
#   (the taxon is genuinely absent) and sampling zeros (present but not sequenced
#   deeply enough to see). The count component can itself emit zeros.
# Zero-HURDLE: one process decides presence/absence, a second models abundance
#   given presence. All zeros are structural; the count part is truncated at zero.
#
# Which to prefer is a biological claim, not a statistical convenience: hurdle
# says "every zero means absent", zero-inflated says "some zeros are undersampling".
# For microbiome data the zero-inflated story is usually the more honest one.

#' Formula for a two-component model.
#'
#' The part before | models the count process; the part after models the
#' probability of a structural zero. They can carry different predictors, and
#' the offset belongs only in the count part.
zero_model_formula <- function(response, count_predictors, zero_predictors = NULL,
                               offset_name = "Offset") {
  zero_predictors <- zero_predictors %||% count_predictors
  cnt <- paste(count_predictors, collapse = " + ")
  if (!is.null(offset_name) && nzchar(offset_name)) {
    cnt <- paste0(cnt, " + offset(", offset_name, ")")
  }
  stats::as.formula(sprintf("%s ~ %s | %s", response, cnt,
                            paste(zero_predictors, collapse = " + ")))
}

#' Zero-inflated Poisson.
#'
#' Assumes variance equals the mean in the count component. Microbiome counts
#' almost always violate that, which is why ZINB usually wins the comparison —
#' but fitting both and comparing is the point (Demo 11).
fit_zip <- function(formula, data, link = "logit") {
  fit <- pscl::zeroinfl(formula, data = data, dist = "poisson", link = link)
  mfg_log("models", "zip", list(formula = paste(deparse(formula), collapse = ""),
    n = mfg_nobs(fit), aic = round(stats::AIC(fit), 2),
    loglik = round(as.numeric(stats::logLik(fit)), 2)))
  structure(list(fit = fit, type = "zip", formula = formula, summary = summary(fit),
    aic = stats::AIC(fit), bic = stats::BIC(fit),
    loglik = as.numeric(stats::logLik(fit)), n = mfg_nobs(fit)),
    class = c("mfg_zero_model", "mfg_model", "list"))
}

#' Zero-inflated negative binomial.
fit_zinb <- function(formula, data, link = "logit") {
  fit <- pscl::zeroinfl(formula, data = data, dist = "negbin", link = link)
  mfg_log("models", "zinb", list(formula = paste(deparse(formula), collapse = ""),
    n = mfg_nobs(fit), aic = round(stats::AIC(fit), 2),
    theta = round(fit$theta, 4),
    loglik = round(as.numeric(stats::logLik(fit)), 2)))
  structure(list(fit = fit, type = "zinb", formula = formula, summary = summary(fit),
    aic = stats::AIC(fit), bic = stats::BIC(fit), theta = fit$theta,
    loglik = as.numeric(stats::logLik(fit)), n = mfg_nobs(fit)),
    class = c("mfg_zero_model", "mfg_model", "list"))
}

#' Zero-hurdle Poisson.
fit_zhp <- function(formula, data, link = "logit") {
  fit <- pscl::hurdle(formula, data = data, dist = "poisson", link = link)
  mfg_log("models", "zhp", list(formula = paste(deparse(formula), collapse = ""),
    n = mfg_nobs(fit), aic = round(stats::AIC(fit), 2),
    loglik = round(as.numeric(stats::logLik(fit)), 2)))
  structure(list(fit = fit, type = "zhp", formula = formula, summary = summary(fit),
    aic = stats::AIC(fit), bic = stats::BIC(fit),
    loglik = as.numeric(stats::logLik(fit)), n = mfg_nobs(fit)),
    class = c("mfg_zero_model", "mfg_model", "list"))
}

#' Zero-hurdle negative binomial.
fit_zhnb <- function(formula, data, link = "logit") {
  fit <- pscl::hurdle(formula, data = data, dist = "negbin", link = link)
  mfg_log("models", "zhnb", list(formula = paste(deparse(formula), collapse = ""),
    n = mfg_nobs(fit), aic = round(stats::AIC(fit), 2),
    loglik = round(as.numeric(stats::logLik(fit)), 2)))
  structure(list(fit = fit, type = "zhnb", formula = formula, summary = summary(fit),
    aic = stats::AIC(fit), bic = stats::BIC(fit),
    loglik = as.numeric(stats::logLik(fit)), n = mfg_nobs(fit)),
    class = c("mfg_zero_model", "mfg_model", "list"))
}

#' Fit all four zero models and compare them correctly.
#'
#' LRT applies only to nested pairs — ZIP within ZINB, ZHP within ZHNB. ZIP and
#' ZHP are not nested in each other, so AIC compares all four while the LRT is
#' restricted to the two valid pairs (Demo 11).
compare_zero_models <- function(formula, data, link = "logit", verbose = TRUE) {
  fits <- list(
    ZIP  = try(fit_zip(formula, data, link), silent = TRUE),
    ZINB = try(fit_zinb(formula, data, link), silent = TRUE),
    ZHP  = try(fit_zhp(formula, data, link), silent = TRUE),
    ZHNB = try(fit_zhnb(formula, data, link), silent = TRUE)
  )
  failed <- names(fits)[vapply(fits, function(f) inherits(f, "try-error"), logical(1))]
  fits <- fits[!names(fits) %in% failed]
  if (!length(fits)) stop("All four zero models failed to fit.", call. = FALSE)

  aic <- data.frame(
    model = names(fits),
    aic = vapply(fits, function(f) f$aic, numeric(1)),
    bic = vapply(fits, function(f) f$bic, numeric(1)),
    loglik = vapply(fits, function(f) f$loglik, numeric(1)),
    stringsAsFactors = FALSE)
  aic$delta_aic <- aic$aic - min(aic$aic)
  aic <- aic[order(aic$aic), ]
  best <- aic$model[1]

  lrt <- list()
  if (all(c("ZIP", "ZINB") %in% names(fits))) {
    lrt$`ZIP vs ZINB` <- lmtest::lrtest(fits$ZIP$fit, fits$ZINB$fit)
  }
  if (all(c("ZHP", "ZHNB") %in% names(fits))) {
    lrt$`ZHP vs ZHNB` <- lmtest::lrtest(fits$ZHP$fit, fits$ZHNB$fit)
  }

  interp <- sprintf(paste("%s has the lowest AIC (%.1f).%s The negative binomial",
    "variants model overdispersion; the Poisson variants assume variance equals",
    "the mean. Hurdle models treat every zero as a true absence, zero-inflated",
    "models allow some zeros to be undersampling — that difference is a",
    "biological claim, so prefer the zero-inflated form unless absence is",
    "genuinely certain."),
    best, min(aic$aic),
    if (nrow(aic) > 1 && aic$delta_aic[2] < 2)
      sprintf(" But %s is within 2 AIC units (delta = %.2f), so the two are not meaningfully distinguishable — prefer the simpler or more interpretable one.",
              aic$model[2], aic$delta_aic[2]) else "")

  mfg_log("models", "zero_models_compared", list(
    fitted = paste(names(fits), collapse = ","),
    failed = if (length(failed)) paste(failed, collapse = ",") else "none",
    best_by_aic = best,
    aic = paste(sprintf("%s=%.1f", aic$model, aic$aic), collapse = " "),
    delta_second = if (nrow(aic) > 1) round(aic$delta_aic[2], 2) else NA))

  out <- list(fits = fits, aic_table = aic, lrt = lrt, best = best,
              failed = failed, interpretation = interp)
  if (verbose) print.mfg_zero_comparison(out)
  class(out) <- c("mfg_zero_comparison", "list")
  out
}

#' @export
print.mfg_zero_comparison <- function(x, ...) {
  cat("=== Zero-model comparison ===\n")
  print(x$aic_table, row.names = FALSE)
  if (length(x$failed)) cat(sprintf("\nFailed to fit: %s\n", paste(x$failed, collapse = ", ")))
  for (nm in names(x$lrt)) {
    cat(sprintf("\nLRT %s (nested):\n", nm)); print(x$lrt[[nm]])
  }
  cat("\n"); cat(strwrap(x$interpretation, width = 78), sep = "\n")
  invisible(x)
}

#' Zero and count distribution diagnostic.
#'
#' Run this before choosing a zero model. The point is to see whether the excess
#' of zeros is large enough to need a two-component model at all — the Demo 11
#' and ps11 diagnostic.
zero_diagnostics <- function(counts, taxon_name = "taxon") {
  n <- length(counts)
  n_zero <- sum(counts == 0)
  # What a Poisson with this mean would predict for P(0).
  lambda <- mean(counts)
  expected_zero_poisson <- n * stats::dpois(0, lambda)
  excess <- n_zero - expected_zero_poisson
  out <- list(
    taxon = taxon_name, n = n, n_zero = n_zero,
    prop_zero = n_zero / n, mean = lambda, variance = stats::var(counts),
    var_mean_ratio = stats::var(counts) / max(lambda, 1e-9),
    expected_zero_poisson = expected_zero_poisson,
    excess_zeros = excess,
    overdispersed = stats::var(counts) / max(lambda, 1e-9) > 1.5,
    zero_inflated = excess > 0.05 * n
  )
  out$verdict <- paste(
    sprintf("%.0f%% zeros (%d of %d).", 100 * out$prop_zero, n_zero, n),
    sprintf("A Poisson with mean %.2f predicts %.0f zeros; there are %.0f excess.",
            lambda, expected_zero_poisson, excess),
    sprintf("Variance/mean = %.1f.", out$var_mean_ratio),
    if (out$zero_inflated && out$overdispersed)
      "Both zero-inflated and overdispersed: ZINB or ZHNB."
    else if (out$zero_inflated) "Zero-inflated but not strongly overdispersed: ZIP may suffice, still compare against ZINB."
    else if (out$overdispersed) "Overdispersed but zeros are not in excess: a plain negative binomial is enough; a zero component is not justified."
    else "Neither strongly zero-inflated nor overdispersed: a Poisson GLM is defensible.")
  mfg_log("models", "zero_diagnostics", list(
    taxon = taxon_name, prop_zero = round(out$prop_zero, 3),
    var_mean_ratio = round(out$var_mean_ratio, 2),
    zero_inflated = out$zero_inflated, overdispersed = out$overdispersed))
  class(out) <- c("mfg_zero_diagnostics", "list")
  out
}

#' @export
print.mfg_zero_diagnostics <- function(x, ...) {
  cat(sprintf("=== Zero/count distribution: %s ===\n", x$taxon))
  cat(strwrap(x$verdict, width = 78), sep = "\n")
  invisible(x)
}

# ── Model comparison ─────────────────────────────────────────────────────────

#' Compare nested models, enforcing the ML requirement.
#'
#' Comparing fixed-effect structures on REML-fitted models is invalid, because
#' REML likelihoods computed under different fixed effects are not comparable.
#' This refuses to do it rather than returning a misleading LRT.
compare_models <- function(..., what = c("fixed", "random")) {
  what <- match.arg(what)
  fits <- list(...)
  objs <- lapply(fits, function(f) if (inherits(f, "mfg_model")) f$fit else f)
  nm <- names(fits)
  if (is.null(nm) || !all(nzchar(nm))) nm <- paste0("model", seq_along(fits))

  if (identical(what, "fixed")) {
    reml <- vapply(fits, function(f)
      inherits(f, "mfg_model") && identical(f$type, "lmm") && isTRUE(f$REML),
      logical(1))
    if (any(reml)) {
      stop("Refusing to compare fixed-effect structures on REML fits (",
           paste(nm[reml], collapse = ", "), ").\n",
           "REML likelihoods under different fixed effects are not comparable. ",
           "Refit with fit_lmm(..., REML = FALSE) and compare those.", call. = FALSE)
    }
  }

  # unname() matters: do.call on a named list passes the models as named
  # arguments, and anova.merMod does not accept them that way — it silently
  # fails and the LRT comes back NULL.
  an <- try(do.call(stats::anova, unname(objs)), silent = TRUE)
  if (inherits(an, "try-error")) {
    warning("Likelihood ratio test unavailable: ",
            conditionMessage(attr(an, "condition")),
            "\nAIC/BIC comparison is still reported.", call. = FALSE)
  }
  tab <- data.frame(model = nm,
    df = vapply(objs, function(o) tryCatch(attr(stats::logLik(o), "df"),
                                           error = function(e) NA_real_), numeric(1)),
    aic = vapply(objs, stats::AIC, numeric(1)),
    bic = vapply(objs, stats::BIC, numeric(1)),
    loglik = vapply(objs, function(o) as.numeric(stats::logLik(o)), numeric(1)),
    stringsAsFactors = FALSE)
  tab$delta_aic <- tab$aic - min(tab$aic)
  best <- tab$model[which.min(tab$aic)]

  caveat <- if (identical(what, "random")) paste(
    "Comparing random-effect structures: the likelihood ratio test is not",
    "reliable here because the null hypothesis sits on the boundary of the",
    "parameter space (a variance of zero), so the reported p-value is",
    "conservative and can point the wrong way. Prefer AIC/BIC and parsimony —",
    "and prefer the simpler structure when they disagree.") else NA_character_

  mfg_log("models", "compared", list(
    what = what, models = paste(nm, collapse = ","),
    aic = paste(sprintf("%s=%.1f", tab$model, tab$aic), collapse = " "),
    best_by_aic = best,
    lrt_available = !inherits(an, "try-error")))

  out <- list(what = what, table = tab,
              lrt = if (inherits(an, "try-error")) NULL else an,
              best = best, caveat = caveat)
  class(out) <- c("mfg_model_comparison", "list")
  out
}

#' @export
print.mfg_model_comparison <- function(x, ...) {
  cat(sprintf("=== Model comparison (%s effects) ===\n", x$what))
  print(x$table, row.names = FALSE)
  if (!is.null(x$lrt)) { cat("\nLikelihood ratio test:\n"); print(x$lrt) }
  cat(sprintf("\nLowest AIC: %s\n", x$best))
  if (nrow(x$table) > 1) {
    srt <- x$table[order(x$table$aic), ]
    if (srt$delta_aic[2] < 2) {
      cat(sprintf("  %s is within 2 AIC units (delta = %.2f) — not meaningfully distinguishable; prefer the simpler model.\n",
                  srt$model[2], srt$delta_aic[2]))
    }
  }
  if (!is.na(x$caveat)) { cat("\n"); cat(strwrap(x$caveat, width = 78, prefix = "! "), sep = "\n") }
  invisible(x)
}

# ── Bayesian joint models ────────────────────────────────────────────────────

#' Bayesian joint model over several responses.
#'
#' The reason to reach for this: LMM and GLMM fit one response at a time. When
#' two outcomes are measured on the same subjects — say health status and a
#' cytokine — modelling them jointly borrows strength across them and estimates
#' their residual correlation. mvbind() binds the responses; each can have its
#' own family (Demo 13).
#'
#' Slow. Set chains/iter deliberately and check convergence before reading
#' anything: Rhat above 1.01 or low effective sample size means the chains have
#' not explored the posterior and the intervals are not trustworthy.
fit_brms_joint <- function(formula, data, family = NULL, chains = 4, iter = 2000,
                           cores = getOption("mc.cores", 2), seed = 42,
                           control = list(adapt_delta = 0.95), ...) {
  fit <- brms::brm(formula, data = data, family = family, chains = chains,
                   iter = iter, cores = cores, seed = seed, control = control,
                   refresh = 0, ...)
  s <- summary(fit)
  rh <- brms::rhat(fit)
  ness <- brms::neff_ratio(fit)
  max_rhat <- max(rh, na.rm = TRUE)
  min_neff <- min(ness, na.rm = TRUE)
  ndiv <- sum(brms::nuts_params(fit, pars = "divergent__")$Value)

  converged <- max_rhat < 1.01 && min_neff > 0.1 && ndiv == 0
  if (!converged) {
    warning(sprintf(paste("Convergence is not clean: max Rhat = %.4f, min",
      "effective-sample-size ratio = %.3f, %d divergent transitions. Do not report",
      "these intervals yet — raise iter, raise adapt_delta, or simplify the model."),
      max_rhat, min_neff, ndiv), call. = FALSE)
  }

  mfg_log("models", "brms_joint", list(
    formula = paste(utils::capture.output(print(formula))[1], collapse = ""),
    chains = chains, iter = iter, seed = seed,
    max_rhat = round(max_rhat, 4), min_neff_ratio = round(min_neff, 3),
    divergent_transitions = ndiv, converged = converged))

  structure(list(fit = fit, type = "brms", formula = formula, summary = s,
    max_rhat = max_rhat, min_neff_ratio = min_neff,
    divergent_transitions = ndiv, converged = converged,
    n = stats::nobs(fit)),
    class = c("mfg_brms", "mfg_model", "list"))
}

#' Bayesian hypothesis test on one predictor.
#'
#' The interpretation rule, lifted from Demo 12/13 because the agent will not
#' infer it reliably: if the 95% credible interval excludes 0, the model is
#' confident the effect is nonzero and the null is rejected in Bayesian terms.
#' brms marks these with a star. Evid.Ratio is the posterior odds in favour of
#' the hypothesis.
brms_hypothesis <- function(model, hypothesis) {
  fit <- if (inherits(model, "mfg_brms")) model$fit else model
  h <- brms::hypothesis(fit, hypothesis)
  tab <- h$hypothesis
  credible <- !is.na(tab$CI.Lower) & !is.na(tab$CI.Upper) &
              ((tab$CI.Lower > 0 & tab$CI.Upper > 0) |
               (tab$CI.Lower < 0 & tab$CI.Upper < 0))
  tab$credible_effect <- credible
  tab$interpretation <- ifelse(credible,
    sprintf("The 95%% credible interval [%.3f, %.3f] excludes 0, so the model is confident this effect is nonzero (direction: %s).",
            tab$CI.Lower, tab$CI.Upper, ifelse(tab$Estimate > 0, "positive", "negative")),
    sprintf("The 95%% credible interval [%.3f, %.3f] includes 0, so the data do not establish a nonzero effect.",
            tab$CI.Lower, tab$CI.Upper))
  mfg_log("models", "brms_hypothesis", list(
    hypothesis = paste(hypothesis, collapse = "; "),
    estimate = round(tab$Estimate, 4),
    ci = sprintf("[%.3f, %.3f]", tab$CI.Lower, tab$CI.Upper),
    evidence_ratio = round(tab$Evid.Ratio, 3), credible = credible))
  structure(list(hypothesis = h, table = tab), class = c("mfg_brms_hypothesis", "list"))
}

#' @export
print.mfg_brms_hypothesis <- function(x, ...) {
  cat("=== Bayesian hypothesis test ===\n")
  print(x$table[, intersect(c("Hypothesis", "Estimate", "Est.Error", "CI.Lower",
                              "CI.Upper", "Evid.Ratio", "Post.Prob", "Star"),
                            names(x$table))], row.names = FALSE)
  for (i in seq_len(nrow(x$table))) {
    cat("\n"); cat(strwrap(x$table$interpretation[i], width = 78, prefix = "  "), sep = "\n")
  }
  invisible(x)
}

#' LOO-CV comparison of Bayesian models.
#'
#' The decision rule, from Demo 13: elpd_diff / se_diff > 2 is strong evidence
#' for the model with higher ELPD. The first row is always the best model, so its
#' elpd_diff is 0. Below 2 the models are not meaningfully distinguishable and
#' the simpler one should be preferred.
compare_loo <- function(..., model_names = NULL) {
  models <- list(...)
  fits <- lapply(models, function(m) if (inherits(m, "mfg_brms")) m$fit else m)
  nm <- model_names %||% names(models) %||% paste0("model", seq_along(models))
  names(fits) <- nm

  if (length(fits) < 2) {
    stop("compare_loo needs at least two models.", call. = FALSE)
  }
  fits <- lapply(fits, function(f) brms::add_criterion(f, "loo"))
  # Build the call from the list rather than naming positions: passing a literal
  # NULL for an absent third model makes loo_compare reject the whole call with
  # "Only model objects can be passed to '...'". do.call also lifts the
  # three-model ceiling.
  cmp <- do.call(brms::loo_compare,
                 c(unname(fits), list(model_names = nm[seq_along(fits)])))
  tab <- as.data.frame(cmp)

  # Row 2 onward carry the difference from the best model.
  ratio <- if (nrow(tab) > 1 && tab$se_diff[2] > 0) abs(tab$elpd_diff[2]) / tab$se_diff[2] else NA_real_
  verdict <- if (is.na(ratio)) "Cannot evaluate: standard error of the difference is zero or missing."
    else if (ratio > 2) sprintf(paste("elpd_diff / se_diff = %.1f (> 2): strong",
      "evidence that %s predicts better than %s."), ratio, rownames(tab)[1], rownames(tab)[2])
    else sprintf(paste("elpd_diff / se_diff = %.1f (<= 2): the models are NOT",
      "meaningfully distinguishable on predictive accuracy. Prefer the simpler",
      "one — do not claim the higher-ELPD model is better."), ratio)

  mfg_log("models", "loo_compared", list(
    models = paste(nm, collapse = ","), best = rownames(tab)[1],
    elpd_diff = if (nrow(tab) > 1) round(tab$elpd_diff[2], 2) else NA,
    se_diff = if (nrow(tab) > 1) round(tab$se_diff[2], 2) else NA,
    ratio = if (is.na(ratio)) NA else round(ratio, 2), verdict = verdict))

  out <- list(comparison = cmp, table = tab, ratio = ratio, verdict = verdict,
              best = rownames(tab)[1])
  class(out) <- c("mfg_loo_comparison", "list")
  out
}

#' @export
print.mfg_loo_comparison <- function(x, ...) {
  cat("=== LOO-CV model comparison ===\n")
  print(x$comparison)
  cat("\n"); cat(strwrap(x$verdict, width = 78), sep = "\n")
  invisible(x)
}

# ── Printing ─────────────────────────────────────────────────────────────────

#' @export
print.mfg_brms <- function(x, ...) {
  # brms models carry no single coefficient matrix (each response has its own),
  # so the generic mfg_model printer would show an empty table. Defer to brms's
  # own summary and add the convergence line, which is the part that decides
  # whether anything here can be reported at all.
  cat("=== Bayesian joint model ===\n")
  cat(sprintf("n = %d | max Rhat = %.4f | min effective-sample-size ratio = %.3f | divergent transitions = %d\n",
              x$n, x$max_rhat, x$min_neff_ratio, x$divergent_transitions))
  cat(sprintf("=> %s\n\n", if (x$converged) "converged — intervals are usable"
              else "NOT converged — do not report these intervals yet"))
  print(x$summary)
  invisible(x)
}

#' @export
print.mfg_model <- function(x, ...) {
  cat(sprintf("=== %s: %s ===\n", toupper(x$type),
              paste(deparse(x$formula), collapse = "")))
  cat(sprintf("n = %d", x$n))
  if (!is.null(x$family)) cat(sprintf(" | family = %s", x$family))
  if (!is.null(x$REML))   cat(sprintf(" | REML = %s", x$REML))
  cat(sprintf(" | AIC = %.2f | BIC = %.2f\n", x$aic, x$bic))
  if (!is.null(x$r_squared)) cat(sprintf("R-squared = %.4f\n", x$r_squared))
  if (!is.null(x$icc) && !is.na(x$icc)) {
    cat(sprintf("ICC = %.4f — %.1f%% of variance is between subjects%s\n", x$icc,
                100 * x$icc,
                if (x$icc > 0.1) ", which is why the random effect is needed" else ""))
  }
  if (isTRUE(x$singular)) cat("! Singular fit — random structure too complex for the data\n")
  if (!is.null(x$pd_hessian) && !isTRUE(x$pd_hessian)) {
    cat("! Non-positive-definite Hessian — standard errors and AIC are unreliable\n")
  }
  if (isTRUE(x$overdispersed)) cat(sprintf("! Overdispersed (deviance/df = %.2f)\n", x$dispersion_ratio))
  if (!is.null(x$independence_warning) && !is.na(x$independence_warning)) {
    cat("\n"); cat(strwrap(x$independence_warning, width = 78, prefix = "! "), sep = "\n")
  }
  cat("\nCoefficients:\n"); print(round(as.data.frame(x$coefficients), 5))
  if (!is.null(x$anova)) { cat("\nEffects:\n"); print(x$anova) }
  invisible(x)
}
