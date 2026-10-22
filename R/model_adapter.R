# An adapter is an ordinary unfitted parsnip spec. Workflows keep preprocessing,
# splits and calibration; every local fit calls `fit(x, y)` on its own prepared
# predictors, and `predict(state, new_data)` sees the same preparation.
#' @export
model_adapter <- function(mode, fit, predict) {
  .exact_argument_names()
  if (missing(mode) || missing(fit) || missing(predict)) {
    stop("Supply `mode`, `fit` and `predict`.", call. = FALSE)
  }
  .validate_choice(mode, c("classification", "regression"), "mode")
  if (!is.function(fit)) stop("`fit` must be a function(x, y).", call. = FALSE)
  if (!is.function(predict)) stop("`predict` must be a function(object, new_data).", call. = FALSE)
  .integration_packages()
  .register_model_adapter()
  callbacks <- structure(list(mode = mode, fit = fit, predict = predict), class = "vennabers_adapter_callbacks")
  spec <- parsnip::new_model_spec("vennabers_adapter", args = list(), eng_args = NULL,
                                  mode = mode, method = NULL, engine = "callbacks")
  parsnip::set_engine(spec, "callbacks", parameters = callbacks)
}

.fit_model_adapter <- function(x, y, callbacks) {
  if (is.data.frame(y) && ncol(y) == 1L) y <- y[[1L]]
  state <- tryCatch(callbacks$fit(as.data.frame(x), y),
                    error = function(e) stop("model_adapter() `fit` callback: ", conditionMessage(e), call. = FALSE))
  structure(list(state = state, predict = callbacks$predict, mode = callbacks$mode,
                 class_levels = if (callbacks$mode == "classification") levels(y)),
            class = "vennabers_adapter_fit")
}

.predict_fit.vennabers_adapter_fit <- function(object, new_data) {
  new_data <- as.data.frame(new_data)
  n <- nrow(new_data)
  # Empty predictions keep their shape without calling the callback.
  if (!n) return(if (object$mode == "regression") tibble::tibble(.pred = double()) else
    tibble::as_tibble(stats::setNames(rep(list(double()), length(object$class_levels)), paste0(".pred_", object$class_levels))))
  p <- tryCatch(object$predict(object$state, new_data),
                error = function(e) stop("model_adapter() `predict` callback: ", conditionMessage(e), call. = FALSE))
  if (object$mode == "regression") {
    if (!is.numeric(p) || !is.null(dim(p)) || length(p) != n || any(!is.finite(p))) {
      stop("model_adapter() `predict` must return a numeric vector with one finite value per input row; ",
           "received ", .adapter_shape(p), " for ", n, " input rows.", call. = FALSE)
    }
    return(tibble::tibble(.pred = as.double(p)))
  }
  classes <- object$class_levels
  expected <- paste(sQuote(classes), collapse = ", ")
  if (!is.matrix(p) || !is.numeric(p) || is.null(colnames(p))) {
    stop("model_adapter() `predict` must return a numeric matrix with columns named by the classes ",
         expected, "; received ", .adapter_shape(p), ".", call. = FALSE)
  }
  found <- colnames(p)
  problem <- if (any(duplicated(found))) c("duplicate column", found[duplicated(found)][1L]) else
    if (any(!found %in% classes)) c("column for unknown class", found[!found %in% classes][1L]) else
      if (any(!classes %in% found)) c("no column for class", classes[!classes %in% found][1L])
  if (!is.null(problem)) {
    stop("model_adapter() `predict` returned ", problem[1L], " ", sQuote(problem[2L]),
         "; columns must be named by the classes ", expected, ".", call. = FALSE)
  }
  if (nrow(p) != n) stop("model_adapter() `predict` returned ", nrow(p), " rows for ", n, " input rows.", call. = FALSE)
  # Columns follow the stored class order by name, never by position.
  p <- p[, classes, drop = FALSE]
  if (!.valid_probability(as.vector(p)) || any(abs(rowSums(p) - 1) > 1e-8)) {
    stop("model_adapter() `predict` must return finite probabilities in [0, 1] that sum to 1 in every row.", call. = FALSE)
  }
  probabilities <- as.data.frame(matrix(as.double(p), n, length(classes)))
  names(probabilities) <- paste0(".pred_", classes)
  tibble::as_tibble(probabilities)
}

.adapter_shape <- function(x) {
  if (is.null(dim(x))) paste0("a ", class(x)[1L], " of length ", length(x)) else
    paste0("a ", class(x)[1L], " with dimensions ", paste(dim(x), collapse = " x "))
}
