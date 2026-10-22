# Shared candidate eligibility, primary scores, reports, and decisions.

.tuning_score <- function(metric, predictions, truth, outcome, classes, event_level,
                          classification) {
  regression <- is.null(classes)
  expected <- if (regression) ".pred" else paste0(".pred_", classes)
  if (!is.data.frame(predictions) || nrow(predictions) != length(truth) ||
      !identical(sort(names(predictions)), sort(expected)) ||
      !all(vapply(predictions, function(x) is.numeric(x) && length(x) == length(truth) &&
                    all(is.finite(x)) && (regression || all(x >= 0 & x <= 1)), logical(1))) ||
      (!regression && any(abs(rowSums(predictions[expected]) - 1) > 1e-8))) {
    stop("Candidate predictions must contain every assessment row and finite ",
         if (regression) "numeric predictions." else "class probabilities.", call. = FALSE)
  }
  frame <- data.frame(truth, predictions[expected], check.names = FALSE)
  names(frame)[1L] <- outcome
  if (regression) {
    if (!inherits(metric, "numeric_metric")) stop("Regression tuning requires a numeric metric.", call. = FALSE)
    result <- rlang::inject(metric(frame, truth = !!rlang::sym(outcome), estimate = !!rlang::sym(".pred")))
  } else if (inherits(metric, "class_metric")) {
    frame$.pred_class <- .class_decision(predictions[expected], classes,
      if (classification == "binary") .event_class(classes, event_level))
    result <- rlang::inject(metric(frame, truth = !!rlang::sym(outcome),
                                   estimate = !!rlang::sym(".pred_class"), event_level = event_level))
  } else if (inherits(metric, "prob_metric")) {
    columns <- if (length(classes) == 2L) .event_class(expected, event_level) else expected
    result <- rlang::inject(metric(frame, truth = !!rlang::sym(outcome),
                                   !!!rlang::syms(columns), event_level = event_level))
  } else stop("Only probability or class metrics are supported for classification tuning.", call. = FALSE)
  if (nrow(result) != 1L || !is.numeric(result$.estimate) ||
      length(result$.estimate) != 1L || !is.finite(result$.estimate)) {
    stop("The selection metric must return one finite estimate.", call. = FALSE)
  }
  result$.estimate[[1L]]
}

.same_candidate <- function(predictions, values) {
  if (!is.data.frame(predictions) || !all(names(values) %in% names(predictions))) return(rep(FALSE, nrow(predictions)))
  Reduce(`&`, Map(function(name) as.character(predictions[[name]]) == as.character(values[[name]]), names(values)))
}

# Keys align with candidate rows: native candidate_id for Grid, parameter values for Bayes.
.score_tuning_candidates <- function(search, result, data, plans, outcome, event_level, classification,
                                    candidates, keys, counter, initial_fits, budget) {
  classes <- levels(data[[outcome]])
  prediction_columns <- if (is.null(classes)) ".pred" else paste0(".pred_", classes)
  candidates$status <- rep("excluded", nrow(candidates))
  candidates$estimate <- rep(NA_real_, nrow(candidates))
  candidates$reason <- rep(NA_character_, nrow(candidates))
  scores <- reports <- decisions <- vector("list", nrow(candidates))
  result_ids <- if ("id2" %in% names(result)) paste(result$id, result$id2, sep = "/") else result$id
  for (i in seq_len(nrow(candidates))) {
    scores[[i]] <- stats::setNames(rep(NA_real_, length(plans)), names(plans))
    reports[[i]] <- decisions[[i]] <- vector("list", length(plans))
    names(reports[[i]]) <- names(decisions[[i]]) <- names(plans)
    errors <- if (!identical(result_ids, names(plans))) "Native resample IDs changed." else character()
    candidate_values <- keys[i, , drop = FALSE]
    for (j in seq_along(plans)) {
      id <- names(plans)[j]
      split <- plans[[j]]
      checked <- tryCatch({
        predictions <- result$.predictions[[j]]
        if (!is.data.frame(predictions)) {
          stop("Native predictions are missing for this resample.", call. = FALSE)
        }
        predictions <- predictions[.same_candidate(predictions, candidate_values), , drop = FALSE]
        rows <- predictions$.row
        if (is.null(rows) || anyNA(rows) || anyDuplicated(rows) ||
            !identical(sort(as.integer(rows)), sort(split$assessment)) ||
            !identical(predictions[[outcome]], data[[outcome]][rows])) {
          stop("Native predictions omit, duplicate, or mismatch assessment rows (",
               nrow(predictions), " predictions for ", length(split$assessment), " assessment rows).", call. = FALSE)
        }
        score_values <- predictions[match(split$assessment, rows), prediction_columns, drop = FALSE]
        primary <- .tuning_score(search$primary_metric, score_values,
          data[[outcome]][split$assessment], outcome, classes, event_level, classification)
        extra <- lapply(search$report_metrics, function(metric) tryCatch(
          .tuning_score(metric, score_values, data[[outcome]][split$assessment],
                        outcome, classes, event_level, classification), error = conditionMessage))
        extract <- result$.extracts[[j]]
        extract <- extract[.same_candidate(extract, candidate_values), , drop = FALSE]
        details <- if (nrow(extract)) extract$.extracts[[1L]] else NULL
        list(primary = primary, extra = extra, details = details)
      }, error = identity)
      if (inherits(checked, "error")) errors <- c(errors, paste0(id, ": ", conditionMessage(checked))) else {
        scores[[i]][[id]] <- checked$primary
        reports[[i]][[id]] <- checked$extra
        decisions[[i]][id] <- list(checked$details$decision)
      }
    }
    if (length(errors)) candidates$reason[i] <- paste(errors, collapse = " | ") else {
      candidates$status[i] <- "eligible"
      candidates$estimate[i] <- mean(scores[[i]])
    }
  }
  eligible <- which(candidates$status == "eligible")
  list(candidates = candidates, eligible = eligible, scores = stats::setNames(scores, candidates$id),
       report_metrics = stats::setNames(reports, candidates$id),
       decisions = stats::setNames(decisions, candidates$id),
       cost = list(candidate_evaluations = nrow(candidates) * length(plans),
                   eligible_candidates = length(eligible), excluded_candidates = nrow(candidates) - length(eligible),
                   observed_base_fits = counter$n - initial_fits), budget = budget)
}
