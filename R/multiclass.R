.class_probability <- function(probability, classes, argument = "probability") {
  if (!is.matrix(probability) || !is.numeric(probability) ||
      any(!is.finite(probability) | probability < 0 | probability > 1)) {
    stop("`", argument, "` must be a numeric matrix of finite values between 0 and 1.", call. = FALSE)
  }
  columns <- colnames(probability)
  if (is.null(columns) || anyNA(columns) || any(!nzchar(columns)) ||
      anyDuplicated(columns) || !setequal(columns, classes)) {
    stop("`", argument, "` column names must be non-empty, unique, and match the classes.", call. = FALSE)
  }
  probability <- probability[, classes, drop = FALSE]
  # rowSums uses extended precision on some platforms. A one-bit normalization
  # difference can cross a calibration knot; use binary64 additions everywhere.
  totals <- Reduce(`+`, as.data.frame(probability))
  if (any(abs(totals - 1) > 1e-8)) {
    stop("`", argument, "` row sums must be within 1e-8 of one.", call. = FALSE)
  }
  probability / totals
}

.validate_multiclass_input <- function(truth, probability) {
  classes <- levels(truth)
  if (!is.factor(truth) || !is.null(dim(truth)) || !length(truth) || anyNA(truth) ||
      length(classes) < 2L || anyNA(classes) || any(!nzchar(classes))) {
    stop("`truth` must be a non-empty factor without missing values and with at least two non-empty levels.", call. = FALSE)
  }
  probability <- .class_probability(probability, classes)
  if (nrow(probability) != length(truth)) {
    stop("`probability` must have one row per `truth` value.", call. = FALSE)
  }
  probability
}

.multiclass_eligibility <- function(truth, method, fold = NULL) {
  classes <- levels(truth)
  absent <- classes[tabulate(as.integer(truth), nbins = length(classes)) == 0L]
  if (length(absent)) {
    stop("Every class in `truth` must be observed in the calibration set; missing: ",
         paste(sQuote(absent), collapse = ", "), ".", call. = FALSE)
  }
  pairs <- if (method == "ovo") utils::combn(classes, 2L, simplify = FALSE) else NULL
  if (!is.null(fold)) for (pair in pairs) {
    absent <- setdiff(unique(fold), fold[truth %in% pair])
    if (length(absent)) {
      stop("Empty calibration for pair [", paste(pair, collapse = ", "),
           "] in fold [", absent[1L], "].", call. = FALSE)
    }
  }
  pairs
}

.pair_probability <- function(probability, pair) {
  total <- rowSums(probability[, pair, drop = FALSE])
  result <- rep(.5, nrow(probability))
  positive <- total > 0
  result[positive] <- probability[positive, pair[2L]] / total[positive]
  result
}

.subproblem_probability <- function(probability, problem) {
  if (length(problem) == 1L) probability[, problem] else .pair_probability(probability, problem)
}

.fit_multiclass <- function(truth, probability, precision, method, fold = NULL) {
  probability <- .validate_multiclass_input(truth, probability)
  classes <- levels(truth)
  pairs <- .multiclass_eligibility(truth, method, fold)
  problems <- if (method == "ovo") pairs else as.list(classes)
  calibrators <- lapply(problems, function(problem) {
    rows <- if (method == "ovo") truth %in% problem else rep(TRUE, length(truth))
    event <- truth[rows] == problem[length(problem)]
    scores <- .subproblem_probability(probability[rows, problem, drop = FALSE], problem)
    if (is.null(fold)) venn_abers(event, scores, precision) else
      cross_venn_abers(event, scores, fold[rows], precision)
  })
  result <- list(classes = classes, method = method, pairs = pairs, calibrators = calibrators,
                 observations = length(truth))
  if (!is.null(fold)) result$folds <- unique(fold)
  structure(result, class = if (is.null(fold)) "venn_abers" else "cross_venn_abers")
}

.predict_multiclass <- function(object, newdata, loss, type) {
  .validate_choice(type, c("response", "class", "details"), "type")
  cross <- inherits(object, "cross_venn_abers")
  if (cross) {
    ids <- names(newdata)
    if (!is.list(newdata) || is.data.frame(newdata) || !is.null(dim(newdata)) ||
        is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) ||
        !setequal(ids, object$folds)) {
      stop("`newdata` must be a list named uniquely by the fitted fold IDs.", call. = FALSE)
    }
    newdata <- lapply(object$folds, function(id) {
      .class_probability(newdata[[id]], object$classes, paste0("newdata fold ", id))
    })
    names(newdata) <- object$folds
    if (length(unique(vapply(newdata, nrow, 0L))) != 1L) {
      stop("`newdata` fold matrices must have the same number of rows.", call. = FALSE)
    }
  } else {
    newdata <- list(.class_probability(newdata, object$classes, "newdata"))
  }
  n <- nrow(newdata[[1L]])
  ovr <- identical(object$method, "ovr")
  problems <- if (ovr) as.list(object$classes) else object$pairs
  details <- vector("list", length(problems))
  points <- matrix(0, n, length(problems))
  for (i in seq_along(problems)) {
    problem <- problems[[i]]
    scores <- lapply(newdata, .subproblem_probability, problem = problem)
    scores <- if (cross) do.call(cbind, scores) else scores[[1L]]
    binary <- stats::predict(object$calibrators[[i]], scores, loss = loss,
                             type = if (cross && type == "details") "details" else "response")
    point <- if (cross && type == "details") binary$prediction$probability else binary$probability
    event <- problem[length(problem)]
    points[, i] <- point
    if (type == "details") {
      parts <- if (cross) binary$folds else data.frame(row = seq_len(n), binary)
      identity <- if (ovr) data.frame(event = rep(event, nrow(parts))) else
        data.frame(non_event = rep(problem[1L], nrow(parts)), event = rep(event, nrow(parts)))
      details[[i]] <- cbind(parts["row"], identity, parts[setdiff(names(parts), "row")])
    }
  }
  prediction <- if (ovr) points / rowSums(points) else .couple_ovo(points, problems, object$classes)
  dimnames(prediction) <- dimnames(newdata[[1L]])
  if (type == "class") return(.class_decision(prediction, object$classes))
  if (type == "response") return(prediction)
  if (ovr) list(prediction = prediction, classes = do.call(rbind, details)) else
    list(prediction = prediction, pairs = do.call(rbind, details))
}

# Both shared class matrices and separately trained binary pairs use this merger.
.couple_ovo <- function(points, pairs, classes) {
  prediction <- matrix(2 - length(classes), nrow(points), length(classes),
                       dimnames = list(NULL, classes))
  for (i in seq_along(pairs)) {
    pair <- pairs[[i]]
    prediction[, pair[1L]] <- prediction[, pair[1L]] + 1 / (1 - points[, i])
    prediction[, pair[2L]] <- prediction[, pair[2L]] + 1 / points[, i]
  }
  prediction <- 1 / prediction
  prediction / rowSums(prediction)
}
