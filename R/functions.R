
## helper to omogenize names
transnames <- function(original, ref) {
  unlist(lapply(original, function(x) {
    i <- names(ref)[sapply(ref, function(y) any(y %in% trimws(x)))]
    ifelse(length(i) > 0, i, x)
  }))
}
