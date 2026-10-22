.valid_probability <- function(x) {
  is.vector(x) && is.numeric(x) && all(is.finite(x) & x >= 0 & x <= 1)
}

.validate_binary_input <- function(truth, probability) {
  if (!is.vector(truth) || !(is.logical(truth) || is.numeric(truth)) ||
      anyNA(truth) || any(!is.finite(truth)) || any(!truth %in% c(0, 1))) {
    stop("`truth` must be a logical or numeric vector containing only 0 and 1.", call. = FALSE)
  }
  if (!length(truth)) {
    stop("`truth` must be non-empty.", call. = FALSE)
  }
  if (!.valid_probability(probability)) {
    stop("`probability` must be a numeric vector of finite values between 0 and 1.", call. = FALSE)
  }
  if (length(probability) != length(truth)) {
    stop("`probability` must have the same length as `truth`.", call. = FALSE)
  }
}

.validate_choice <- function(value, choices, argument) {
  if (!is.character(value) || !is.vector(value) || length(value) != 1L ||
      is.na(value) || !value %in% choices) {
    stop("`", argument, "` must be one of: ", paste(choices, collapse = ", "), ".", call. = FALSE)
  }
}

# One finite whole number in [min, max].
.is_count <- function(x, min = 1, max = .Machine$integer.max) {
  is.numeric(x) && is.vector(x) && length(x) == 1L && is.finite(x) &&
    x == floor(x) && x >= min && x <= max
}

# One finite number strictly between 0 and 1.
.is_proportion <- function(x) {
  is.numeric(x) && is.vector(x) && length(x) == 1L && is.finite(x) && x > 0 && x < 1
}

.validate_precision <- function(precision) {
  if (!is.null(precision) && !.is_count(precision, -Inf, Inf)) {
    stop("`precision` must be NULL or one finite numeric whole number.", call. = FALSE)
  }
}

.validate_dots <- function(...) {
  if (...length()) stop("`...` must be empty; unsupported argument: ",
                        paste(names(list(...)), collapse = ", "), ".", call. = FALSE)
}

.validate_names <- function(x, role) {
  if (!is.character(x) || !is.vector(x) || !length(x)) {
    stop("The ", role, " names must be a nonempty character vector.", call. = FALSE)
  }
  bad <- is.na(x) | !nzchar(x) | duplicated(x)
  if (any(bad)) stop("Invalid ", role, " name ", sQuote(x[which(bad)[1L]]),
                     ": names must be unique, nonmissing and nonempty; rename explicitly.", call. = FALSE)
}

.validate_truth <- function(y, outcome, classification, mode) {
  if (mode == "regression") {
    .validate_regression_vector(y, paste0(outcome, " regression truth"))
    return(NULL)
  }
  .validate_class_levels(y, outcome, classification)
}

.validate_class_levels <- function(y, outcome, classification) {
  binary <- classification == "binary"
  if (!is.factor(y) || anyNA(y) || nlevels(y) < 2L || (binary && nlevels(y) != 2L)) {
    stop("Outcome ", sQuote(outcome), " must be a nonmissing factor with ",
         if (binary) "exactly two" else "at least two", " levels.", call. = FALSE)
  }
  classes <- levels(y)
  .validate_names(classes, "class")
  bad <- classes == "class" | startsWith(classes, ".pred_")
  if (any(bad)) stop("Class ", sQuote(classes[which(bad)[1L]]),
                     " is reserved for native predictions; rename it explicitly.", call. = FALSE)
  classes
}
