
## helper to omogenize names
transnames <- function(original, ref) {
  unlist(lapply(original, function(x) {
    i <- names(ref)[sapply(ref, function(y) any(y %in% trimws(x)))]
    ifelse(length(i) > 0, i, x)
  }))
}

## Genuine Type-III likelihood-ratio test for `clm` (ordinal::clm) objects. Neither
## car::Anova() (Anova.clm just relabels the object and delegates to Anova.default, which
## is a Wald test based on vcov()) nor ordinal::clm's own single-model anova() (explicitly
## labeled "Wald chi-square tests" in its own output heading) provide a refit-based
## per-term LR test for this model class -- unlike glm/multinom, whose car::Anova() methods
## do refit reduced models internally. This refits the model once per term (dropping that
## term, holding all others -- i.e. Type III), and compares deviances directly, producing a
## table shaped like car::Anova(type=3, test="Chisq")'s output (Df, Chisq, Pr(>Chisq)) so it
## can be dropped into eta2.default()/eta2_partial.default() via their `anova_table` override.
.clm_anova_lr <- function(object) {

  full_terms <- attr(stats::terms(object), "term.labels")
  if (length(full_terms) == 0)
    stop("model has no terms to test")

  full_ll <- as.numeric(stats::logLik(object))
  full_df <- length(stats::coef(object))

  rows <- lapply(full_terms, function(term) {
    # data=object$model (rather than relying on update()'s default re-evaluation of the
    # original call in parent.frame()) keeps this self-contained: the original `data`
    # argument may reference a variable that isn't in scope wherever this helper is called
    # from, but the fitted object's own stored model frame always is.
    reduced <- stats::update(object, stats::as.formula(paste("~ . -", term)), data = object$model)
    df <- full_df - length(stats::coef(reduced))
    chisq <- 2 * (full_ll - as.numeric(stats::logLik(reduced)))
    c(Df = df, Chisq = chisq)
  })

  tab <- as.data.frame(do.call(rbind, rows))
  rownames(tab) <- full_terms
  tab[["Pr(>Chisq)"]] <- stats::pchisq(tab$Chisq, tab$Df, lower.tail = FALSE)
  tab
}
