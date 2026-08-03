
## helper to omogenize names
transnames <- function(original, ref) {
  unlist(lapply(original, function(x) {
    i <- names(ref)[sapply(ref, function(y) any(y %in% trimws(x)))]
    ifelse(length(i) > 0, i, x)
  }))
}

## Genuine Type-III likelihood-ratio test for `clm` (ordinal::clm) objects. Neither
## car::Anova() nor ordinal::clm::anova() provides a refit-based per-term LR test
## for this model class. The reduced fits below use the original model-matrix
## columns so that dropping a main effect does not re-encode an interaction.
.clm_anova_lr <- function(object) {
  full_terms <- attr(stats::terms(object), "term.labels")
  if (length(full_terms) == 0) {
    stop("model has no terms to test")
  }

  model_frame <- stats::model.frame(object)
  terms_object <- stats::delete.response(stats::terms(object))
  design <- stats::model.matrix(terms_object, data = model_frame)
  assignment <- attr(design, "assign")
  predictor_columns <- which(assignment != 0L)
  predictor_names <- paste0(".clm_x", seq_along(predictor_columns))
  response_name <- names(model_frame)[1]

  design_data <- model_frame[1]
  if (length(predictor_columns) > 0) {
    predictors <- as.data.frame(design[, predictor_columns, drop = FALSE])
    names(predictors) <- predictor_names
    design_data[predictor_names] <- predictors
  }

  full_ll <- as.numeric(stats::logLik(object))
  full_df <- length(stats::coef(object))
  rows <- lapply(seq_along(full_terms), function(term_index) {
    keep <- assignment[predictor_columns] != term_index
    reduced_formula <- stats::reformulate(
      predictor_names[keep],
      response = response_name,
      intercept = attr(terms_object, "intercept") == 1L
    )
    reduced <- stats::update(object, formula = reduced_formula, data = design_data)
    df <- full_df - length(stats::coef(reduced))
    chisq <- 2 * (full_ll - as.numeric(stats::logLik(reduced)))
    c(Df = df, Chisq = chisq)
  })

  tab <- as.data.frame(do.call(rbind, rows))
  rownames(tab) <- full_terms
  tab[["Pr(>Chisq)"]] <- stats::pchisq(tab$Chisq, tab$Df, lower.tail = FALSE)
  tab
}
