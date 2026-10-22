#' @export
venn_abers_one_hot <- function(x) {
  # Represent existing class decisions; never choose a class from probabilities.
  if (is.data.frame(x)) {
    if (!".pred_class" %in% names(x) || !is.factor(x[[".pred_class"]])) {
      stop("`x` must contain a factor `.pred_class` column.", call. = FALSE)
    }
    x <- x[[".pred_class"]]
  }
  if (is.logical(x) && is.null(dim(x))) {
    classes <- c("FALSE", "TRUE")
    index <- as.integer(x) + 1L
  } else if (is.factor(x) && is.null(dim(x))) {
    classes <- levels(x)
    if (!length(classes)) stop("`x` must be a factor with at least one level.", call. = FALSE)
    index <- as.integer(x)
  } else {
    stop("`x` must be a factor, a logical vector, or a data frame with a factor `.pred_class` column.",
         call. = FALSE)
  }
  result <- matrix(0, length(index), length(classes), dimnames = list(NULL, classes))
  known <- !is.na(index)
  result[cbind(which(known), index[known])] <- 1
  result[!known, ] <- NA_real_
  result
}
