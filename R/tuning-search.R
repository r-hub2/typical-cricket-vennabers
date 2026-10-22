# Native Grid and Bayes search orchestration.

.tuning_resamples <- function(resamples, data, label) {
  plan <- if (is.function(resamples)) resamples(data) else resamples
  if (!inherits(plan, "rset") || !nrow(plan) ||
      !all(c("id", "splits") %in% names(plan))) {
    stop(label, " resamples must be a nonempty rsample rset.", call. = FALSE)
  }
  ids <- if ("id2" %in% names(plan)) paste(plan$id, plan$id2, sep = "/") else as.character(plan$id)
  .validate_names(ids, paste(label, "resample ID"))
  splits <- lapply(seq_len(nrow(plan)), function(i) .split_positions(plan$splits[[i]], data, ids[i]))
  list(rset = plan, positions = stats::setNames(splits, ids))
}

# Public only as the parsnip proxy engine entry point.
#' @export
fit_vennabers_proxy <- function(x, y, outcome, parameters, ...) {
  candidate_values <- list(...)
  if (!is.null(parameters$candidates)) {
    index <- match(candidate_values$candidate_id, parameters$candidates$id)
    if (length(index) != 1L || is.na(index)) stop("Unknown grid candidate ID.", call. = FALSE)
    candidate_values <- as.list(parameters$candidates[index,
      setdiff(names(parameters$candidates), "id"), drop = FALSE])
  }
  if (is.data.frame(y) && ncol(y) == 1L) y <- y[[1L]]
  record <- if (is.null(parameters$fit_log)) function(values, error) NULL else parameters$fit_log$record
  fitted <- withCallingHandlers({
    candidate <- parameters$candidate_plan(candidate_values)
    if (is.null(candidate$parameters)) {
      data <- as.data.frame(x)
      data[[outcome]] <- y
      .fit_base_model(candidate$workflow, data, parameters$counter)
    } else {
      data <- .calibrated_fit_data(x, y, candidate$workflow, outcome, candidate$parameters)
      seed <- candidate$parameters$requested$seed
      if (!is.null(seed)) withr::local_seed(seed)
      .fit_calibrated_components(candidate$workflow, data, outcome,
        candidate$parameters, candidate$tuning_plan)
    }
  }, error = function(e) record(candidate_values, conditionMessage(e)))
  record(candidate_values, NA_character_)
  structure(list(fitted = fitted, class_levels = levels(y),
    event_class = if (identical(parameters$classification, "binary"))
      .event_class(levels(y), parameters$event_level)),
    class = "vennabers_grid_fit")
}

# A candidate is a plain base fit or a nested calibrated fit.
.predict_fit.vennabers_grid_fit <- function(object, new_data) {
  if (inherits(object$fitted, "vennabers_fit")) {
    if (!is.null(object$fitted$schema)) .schema_encoded_data(new_data, object$fitted$schema)
    return(.predict_fit(object$fitted, new_data))
  }
  classes <- object$class_levels
  if (is.null(classes)) return(stats::predict(object$fitted, new_data, type = "numeric"))
  stats::predict(object$fitted, new_data, type = "prob")[paste0(".pred_", classes)]
}

.bayes_proxy_workflow <- function(data, outcome, event_level, classification,
                                  candidate_plan, counter, ids, fit_log) {
  .register_grid_proxy()
  mode <- if (is.numeric(data[[outcome]])) "regression" else "classification"
  spec <- parsnip::new_model_spec("vennabers_grid", args = list(), eng_args = NULL,
                                  mode = mode, method = NULL, engine = "proxy")
  spec <- parsnip::set_engine(spec, "proxy", outcome = outcome,
    parameters = list(candidate_plan = candidate_plan,
      candidates = NULL, event_level = event_level, classification = classification,
      counter = counter, fit_log = fit_log))
  for (id in unname(ids)) spec$eng_args[[id]] <- rlang::new_quosure(tune::tune(id), env = emptyenv())
  carrier <- rlang::inject(workflows::add_variables(
    workflows::workflow(),
    outcomes = .integration_select_outcome(!!outcome, !!classification, !!mode),
    predictors = tidyselect::everything(),
    blueprint = hardhat::default_xy_blueprint(intercept = FALSE, allow_novel_levels = TRUE)))
  workflows::add_model(carrier, spec)
}

# The stand-in ID lets one native tune_grid() run the complete candidate grid.
.grid_proxy_workflow <- function(data, outcome, event_level, classification,
                                 candidates, candidate_plan, counter) {
  .register_grid_proxy()
  mode <- if (is.numeric(data[[outcome]])) "regression" else "classification"
  spec <- parsnip::new_model_spec("vennabers_grid", args = list(), eng_args = NULL,
                                  mode = mode, method = NULL, engine = "proxy")
  spec <- parsnip::set_engine(spec, "proxy", candidate_id = tune::tune("candidate_id"),
                              outcome = outcome,
                              parameters = list(candidate_plan = candidate_plan,
                                candidates = candidates, event_level = event_level,
                                classification = classification,
                                counter = counter))
  carrier <- rlang::inject(workflows::add_variables(
    workflows::workflow(),
    outcomes = .integration_select_outcome(!!outcome, !!classification, !!mode),
    predictors = tidyselect::everything(),
    blueprint = hardhat::default_xy_blueprint(intercept = FALSE, allow_novel_levels = TRUE)))
  workflows::add_model(carrier, spec)
}

.grid_budget <- function(search, plans, candidate_plan, classes) {
  candidate <- candidate_plan(as.list(search$grid[1L, , drop = FALSE]))
  paths <- 1
  model_candidates <- 0L
  model_resamples <- 0L
  inner_paths <- 1
  inner <- 1L
  pair_count <- 1L
  if (!is.null(candidate$parameters)) {
    settings <- candidate$parameters$requested
    if (identical(settings$model_training, "pairwise")) pair_count <- choose(length(classes), 2L)
    folds <- if (is.function(settings$splits)) NA_integer_ else
      if (settings$method == "cross") settings$v else 1L
    # Each VA/full path holds J inner B/P pairs, and each B runs its own search.
    post <- settings$postprocessing
    if (!is.null(post) && post$strategy == "ensemble") inner <- if (is.null(post$splits)) as.integer(post$v) else NA_integer_
    sources <- folds + as.integer(identical(settings$prediction_model, "full"))
    paths <- pair_count * sources * inner
    model <- candidate$tuning_plan$model
    if (!is.null(model)) {
      model_candidates <- if (identical(model$method, "bayes")) model$initial + model$iter else
        if (!is.null(model$grid_size)) model$grid_size else nrow(model$grid)
      model_resamples <- if (is.function(model$resamples)) NA_integer_ else nrow(model$resamples)
      if (candidate$tuning_plan$strategy == "joint") inner_paths <- sources * inner
    }
  }
  unit <- paths * (1 + model_candidates * model_resamples * inner_paths)
  list(candidates = nrow(search$grid), resamples = length(plans), pair_count = pair_count,
       gross_base_fit_estimate = nrow(search$grid) * length(plans) * unit,
       paths = paths, inner_pairs = inner,
       model_candidates = model_candidates,
       model_resamples = model_resamples, inner_paths = inner_paths,
       reused_base_fits = 0L)
}

.grid_select <- function(search, data, outcome, event_level, classification,
                         candidate_plan, label, counter) {
  # yardstick caches column selections in tune's shared progress environment.
  # Nested binary searches must not reuse or overwrite multiclass selections.
  rlang::local_bindings(known_selections = NULL,
    .env = utils::getFromNamespace("tune_env", "tune"))
  if (!is.null(search$grid_size)) {
    params <- dials::finalize(search$param_info, data)
    generated <- if (is.null(search$seed)) dials::grid_space_filling(params, size = search$grid_size) else
      withr::with_seed(search$seed, dials::grid_space_filling(params, size = search$grid_size))
    search$grid <- .checked_grid(generated, search$ids, label)
    search$candidate_ids <- sprintf("%s_%03d", label, seq_len(nrow(search$grid)))
  }
  plan <- .tuning_resamples(search$resamples, data, label)
  plans <- plan$positions
  budget <- .grid_budget(search, plans, candidate_plan, classes = levels(data[[outcome]]))
  initial_fits <- counter$n
  candidates <- data.frame(id = search$candidate_ids)
  candidates_with_values <- cbind(candidates["id"], search$grid)
  native_metrics <- if (inherits(search$primary_metric, "class_metric"))
    yardstick::metric_set(search$primary_metric, yardstick::brier_class) else
    yardstick::metric_set(search$primary_metric)
  control <- if (is.null(search$control)) tune::control_grid(save_pred = TRUE) else search$control
  control$allow_par <- FALSE
  control$event_level <- event_level
  control$extract <- function(w) {
    proxy <- workflows::extract_fit_engine(w)
    state <- proxy$fitted
    if (inherits(state, "vennabers_fit")) list(decision = state$tuning$models) else
      list(decision = NULL)
  }
  result <- tryCatch(suppressMessages(suppressWarnings(tune::tune_grid(
    .grid_proxy_workflow(data, outcome, event_level, classification,
                         candidates_with_values, candidate_plan, counter),
    resamples = plan$rset,
    grid = data.frame(candidate_id = candidates$id),
    param_info = dials::parameters(dials::new_qual_param(
      type = "character", values = candidates$id,
      label = c(candidate_id = "Candidate"))),
    metrics = native_metrics,
    control = control))), error = identity)
  if (inherits(result, "error")) stop("Native ", label, " grid failed: ",
                                      conditionMessage(result), call. = FALSE)
  for (option in intersect(c("m", "epsilon"), names(search$ids))) {
    for (j in seq_len(nrow(result))) {
      predictions <- result$.predictions[[j]]
      if (is.data.frame(predictions)) {
        column <- if (identical(outcome, option)) paste0("va_", option) else option
        predictions[[column]] <- search$grid[[search$ids[[option]]]][
          match(predictions$candidate_id, candidates$id)]
        result$.predictions[[j]] <- predictions
      }
    }
  }
  scored <- .score_tuning_candidates(search, result, data, plans, outcome, event_level,
    classification, candidates, data.frame(candidate_id = candidates$id), counter, initial_fits, budget)
  candidates <- scored$candidates
  eligible <- scored$eligible
  if (!length(eligible)) stop("No eligible ", label, " grid candidate: ",
    paste(c(stats::na.omit(candidates$reason), unique(tune::collect_notes(result)$note)),
      collapse = " | "), call. = FALSE)
  direction <- attr(search$primary_metric, "direction")
  if (!direction %in% c("minimize", "maximize")) stop("Selection metric needs a direction.", call. = FALSE)
  best <- eligible[if (direction == "minimize") which.min(candidates$estimate[eligible]) else
                     which.max(candidates$estimate[eligible])]
  list(values = as.list(search$grid[best, , drop = FALSE]),
       diagnostics = list(chosen = candidates$id[best], candidates = candidates, grid = search$grid,
         tie_break = "first candidate in grid order",
         scores = scored$scores,
         report_metrics = scored$report_metrics,
         decisions = scored$decisions,
         native = result,
         resamples = plans, primary_metric = search$primary_metric,
         direction = direction,
         budget = budget,
         cost = scored$cost,
         seed = NULL))
}

# Parameter values of the iterations tune_bayes added after `after`.
.bayes_new_values <- function(result, ids, after) {
  metrics <- tryCatch(tune::collect_metrics(result, summarize = FALSE), error = identity)
  if (inherits(metrics, "error") || !all(c(unname(ids), ".iter") %in% names(metrics))) return(data.frame())
  unique(as.data.frame(metrics[metrics$.iter > after, unname(ids), drop = FALSE], stringsAsFactors = FALSE))
}

# tune_bayes keeps no trace of a proposal that failed on every resample, so the
# proxy reports each fit it was asked for: the values and the error, if any.
.bayes_fit_log <- function(ids) {
  log <- new.env()
  log$fits <- list()
  log$record <- function(values, error) {
    log$fits[[length(log$fits) + 1L]] <- list(
      values = as.data.frame(values[unname(ids)], stringsAsFactors = FALSE), error = error)
  }
  log
}

.note_rows <- function(step = character(), resample = NA_character_, location = NA_character_, note = character()) {
  if (!length(note)) step <- resample <- location <- character()
  data.frame(step = step, resample = resample, location = location, note = note, stringsAsFactors = FALSE)
}

# Native notes of excluded candidates would be lost with the filtered GP input.
.bayes_notes <- function(result, step, rows = seq_len(nrow(result))) {
  ids <- if ("id2" %in% names(result)) paste(result$id, result$id2, sep = "/") else result$id
  notes <- lapply(rows, function(i) {
    x <- result$.notes[[i]]
    if (!is.data.frame(x) || !nrow(x) || !"note" %in% names(x)) return(NULL)
    .note_rows(step, ids[[i]], if ("location" %in% names(x)) x$location else NA_character_, x$note)
  })
  do.call(rbind, c(list(.note_rows()), notes))
}

# Keep parameter columns distinct from candidate metadata.
.bayes_candidate_columns <- function(values) {
  names(values) <- make.unique(c("id", "status", "estimate", "reason", names(values)), sep = "_")[-seq_len(4L)]
  values
}

# A proposal dropped by tune_bayes has no native scores; it is excluded as a whole.
.excluded_candidates <- function(values, reason, plans) {
  parameters <- values
  values <- .bayes_candidate_columns(values)
  values$id <- rep(NA_character_, nrow(values))
  values$status <- rep("excluded", nrow(values))
  values$estimate <- rep(NA_real_, nrow(values))
  values$reason <- rep(reason, nrow(values))
  missing <- stats::setNames(rep(NA_real_, length(plans)), names(plans))
  empty <- stats::setNames(vector("list", length(plans)), names(plans))
  list(candidates = values, values = parameters, scores = rep(list(missing), nrow(values)),
       report_metrics = rep(list(empty), nrow(values)), decisions = rep(list(empty), nrow(values)))
}

.append_history <- function(history, scored) {
  history$candidates <- rbind(history$candidates, scored$candidates)
  history$values <- rbind(history$values, scored$values)
  for (part in c("scores", "report_metrics", "decisions")) history[[part]] <- c(history[[part]], scored[[part]])
  history
}

.candidate_key <- function(values) {
  apply(values, 1L, function(row) paste(names(values), as.character(row), sep = "=", collapse = "\r"))
}

.stable_bayes_ids <- function(scored, ids, label, known = character()) {
  keys <- .candidate_key(scored$values[unname(ids)])
  for (key in keys) if (!key %in% names(known)) known[key] <- sprintf("%s_%03d", label, length(known) + 1L)
  current <- unname(known[keys])
  scored$candidates$id <- current
  names(scored$scores) <- names(scored$report_metrics) <- names(scored$decisions) <- current
  list(scored = scored, known = known)
}

.bayes_best <- function(estimates, direction) {
  if (direction == "minimize") min(estimates) else max(estimates)
}

# The one record of a Bayes search: every attempt with a stable ID, the
# complete eligible results that feed the GP, per-step diagnostics, and the
# stall count behind no_improve.
.bayes_ledger <- function(scored, initial, notes, ids, label, direction) {
  identified <- .stable_bayes_ids(scored, ids, label)
  scored <- identified$scored
  values <- scored$values[scored$eligible, unname(ids), drop = FALSE]
  list(history = scored, known = identified$known, ids = ids, label = label, direction = direction,
       values = values, native = tune::filter_parameters(initial, parameters = values),
       steps = list(initial = scored$candidates), notes = notes,
       best = .bayes_best(scored$candidates$estimate[scored$eligible], direction), stale = 0L)
}

# One continuation step as one outcome: "failed" (tune_bayes raised),
# "repeated" (no new candidate), "excluded" (tune_bayes dropped a proposal that
# failed on every resample), or "evaluated" (scored rows of the new iteration).
# Messages are kept as notes when tune_bayes left no native notes for the step.
.bayes_step <- function(name, ledger, run, fit_log, score, plans) {
  fit_log$fits <- list()
  messages <- character()
  keep <- function(restart) function(condition) {
    messages <<- c(messages, trimws(conditionMessage(condition)))
    invokeRestart(restart)
  }
  result <- tryCatch(withCallingHandlers(run(ledger$native),
    message = keep("muffleMessage"), warning = keep("muffleWarning")), error = identity)
  logged <- .note_rows(name, note = unique(messages))
  if (inherits(result, "error")) return(list(kind = "failed", error = conditionMessage(result), notes = logged))
  after <- if (".iter" %in% names(ledger$native)) max(ledger$native$.iter) else 0L
  rows <- which(result$.iter > after)
  attempted <- lapply(fit_log$fits, `[[`, "values")
  values <- if (length(rows)) .bayes_new_values(result, ledger$ids, after) else
    if (length(attempted)) unique(do.call(rbind, attempted)) else ledger$values[0L, , drop = FALSE]
  notes <- if (length(rows)) .bayes_notes(result, name, rows) else logged
  if (!nrow(values) || all(.candidate_key(values) %in% names(ledger$known))) {
    return(list(kind = "repeated", values = values, notes = notes))
  }
  if (!length(rows)) {
    errors <- unique(stats::na.omit(vapply(fit_log$fits, `[[`, "", "error")))
    return(list(kind = "excluded", notes = notes, scored = .excluded_candidates(values,
      paste(c("All resamples failed", errors), collapse = " | "), plans)))
  }
  list(kind = "evaluated", notes = notes, result = result,
       scored = score(result[rows, , drop = FALSE], values))
}

.bayes_record <- function(ledger, name, step) {
  ledger$notes <- rbind(ledger$notes, step$notes)
  if (step$kind %in% c("failed", "repeated")) {
    ledger$steps[[name]] <- if (step$kind == "failed") list(error = step$error) else list(repeated = step$values)
    ledger$stale <- ledger$stale + 1L
    return(ledger)
  }
  identified <- .stable_bayes_ids(step$scored, ledger$ids, ledger$label, ledger$known)
  scored <- identified$scored
  ledger$known <- identified$known
  ledger$history <- .append_history(ledger$history, scored)
  ledger$steps[[name]] <- scored$candidates
  eligible <- which(scored$candidates$status == "eligible")
  if (!length(eligible)) {
    ledger$stale <- ledger$stale + 1L
    return(ledger)
  }
  ledger$values <- rbind(ledger$values, scored$values[eligible, unname(ledger$ids), drop = FALSE])
  ledger$native <- tune::filter_parameters(step$result, parameters = ledger$values)
  estimate <- .bayes_best(scored$candidates$estimate[eligible], ledger$direction)
  improved <- if (ledger$direction == "minimize") estimate < ledger$best else estimate > ledger$best
  ledger$stale <- if (improved) 0L else ledger$stale + 1L
  if (improved) ledger$best <- estimate
  ledger
}

.bayes_select <- function(search, data, outcome, event_level, classification,
                          candidate_plan, label, counter) {
  # Keep the same metric-cache boundary as Grid, including Bayes initialization.
  rlang::local_bindings(known_selections = NULL,
    .env = utils::getFromNamespace("tune_env", "tune"))
  plan <- .tuning_resamples(search$resamples, data, label)
  plans <- plan$positions
  initial_fits <- counter$n
  started <- proc.time()[3L]
  control <- search$control
  control$time_limit <- NA_real_
  control$event_level <- event_level
  control$extract <- function(w) {
    proxy <- workflows::extract_fit_engine(w)
    state <- proxy$fitted
    if (inherits(state, "vennabers_fit")) list(decision = state$tuning$models) else list(decision = NULL)
  }
  fit_log <- .bayes_fit_log(search$ids)
  proxy <- .bayes_proxy_workflow(data, outcome, event_level, classification, candidate_plan,
    counter, search$ids, fit_log)
  params <- dials::finalize(search$param_info, data)
  initial_grid <- withr::with_seed(control$seed, dials::grid_space_filling(params, size = search$initial))
  budget <- .grid_budget(list(grid = initial_grid), plans, candidate_plan, classes = levels(data[[outcome]]))
  budget$candidates <- nrow(initial_grid) + search$iter
  budget$gross_base_fit_estimate <- budget$gross_base_fit_estimate *
    budget$candidates / nrow(initial_grid)
  native_metrics <- if (inherits(search$primary_metric, "class_metric"))
    yardstick::metric_set(search$primary_metric, yardstick::brier_class) else yardstick::metric_set(search$primary_metric)
  grid_control <- tune::control_grid(save_pred = TRUE, allow_par = FALSE,
    event_level = event_level, extract = control$extract, parallel_over = NULL)
  initial <- tryCatch(suppressMessages(suppressWarnings(tune::tune_grid(proxy, resamples = plan$rset,
    grid = initial_grid, param_info = params, metrics = native_metrics, control = grid_control))), error = identity)
  if (inherits(initial, "error")) stop("Native ", label, " Bayes initialization failed: ", conditionMessage(initial), call. = FALSE)
  score <- function(result, values) {
    candidates <- .bayes_candidate_columns(values)
    candidates$id <- sprintf("%s_%03d", label, seq_len(nrow(candidates)))
    scored <- .score_tuning_candidates(search, result, data, plans, outcome, event_level, classification,
      candidates, values, counter, initial_fits, budget)
    scored$values <- values
    scored
  }
  scored <- score(initial, as.data.frame(initial_grid))
  required <- nrow(params) + 1L
  if (!nrow(scored$candidates)) {
    stop("Bayes needs at least ", required, " complete initial ", label,
         " candidates; found 0.", call. = FALSE)
  }
  notes <- .bayes_notes(initial, "initial")
  if (length(scored$eligible) < required) {
    # As in Grid: excluded candidates and native notes explain the shortfall.
    reasons <- unique(c(stats::na.omit(scored$candidates$reason), notes$note))
    stop("Bayes needs at least ", required, " complete initial ", label, " candidates; found ",
         length(scored$eligible), if (length(reasons)) paste0(": ", paste(reasons, collapse = " | ")) else ".",
         call. = FALSE)
  }
  direction <- attr(search$primary_metric, "direction")
  ledger <- .bayes_ledger(scored, initial, notes, search$ids, label, direction)
  metric_names <- unique(tune::collect_metrics(initial)$.metric)
  metric_name <- if (inherits(search$primary_metric, "class_metric"))
    setdiff(metric_names, "brier_class")[[1L]] else metric_names[[1L]]
  metric_functions <- list(search$primary_metric)
  names(metric_functions) <- metric_name
  if (inherits(search$primary_metric, "class_metric")) metric_functions$brier_class <- yardstick::brier_class
  native_metrics <- do.call(yardstick::metric_set, metric_functions)
  reason <- "budget"
  for (i in seq_len(search$iter)) {
    if (!is.na(search$control$time_limit) && proc.time()[3L] - started >= search$control$time_limit) {
      reason <- "time_limit"
      break
    }
    step_control <- control
    step_control$seed <- control$seed + i
    run <- function(valid) tune::tune_bayes(proxy, resamples = plan$rset, initial = valid, iter = 1L,
      param_info = params, metrics = native_metrics, control = step_control)
    name <- paste0("step_", i)
    step <- .bayes_step(name, ledger, run, fit_log, score, plans)
    ledger <- .bayes_record(ledger, name, step)
    if (step$kind == "repeated") {
      reason <- "repeated_candidate"
      break
    }
    if (ledger$stale >= search$control$no_improve) {
      reason <- "no_improve"
      break
    }
  }
  history <- ledger$history
  valid <- ledger$native
  eligible <- which(history$candidates$status == "eligible")
  values <- history$values[eligible, unname(search$ids), drop = FALSE]
  pick <- if (direction == "minimize") which.min(history$candidates$estimate[eligible]) else
    which.max(history$candidates$estimate[eligible])
  best <- eligible[pick]
  attempts <- nrow(initial_grid) + length(ledger$steps) - 1L
  history$cost$candidate_evaluations <- attempts * length(plans)
  history$cost$eligible_candidates <- length(eligible)
  history$cost$excluded_candidates <- nrow(history$candidates) - length(eligible)
  history$cost$observed_base_fits <- counter$n - initial_fits
  for (option in intersect(c("m", "epsilon"), names(search$ids))) {
    column <- if (identical(outcome, option)) paste0("va_", option) else option
    valid$.predictions <- lapply(valid$.predictions, function(predictions) {
      predictions[[column]] <- predictions[[search$ids[[option]]]]
      predictions
    })
  }
  list(values = as.list(values[pick, , drop = FALSE]), diagnostics = c(list(chosen = history$candidates$id[best],
    selected = as.list(values[pick, , drop = FALSE]), candidates = history$candidates, grid = NULL,
    tie_break = "first candidate in attempt order", scores = history$scores,
    report_metrics = history$report_metrics, decisions = history$decisions, native = valid,
    resamples = plans, primary_metric = search$primary_metric, direction = direction,
    seed = control$seed, bayes_seed = control$seed, steps = ledger$steps, notes = ledger$notes, stop_reason = reason,
    attempts = attempts, no_improve = ledger$stale, elapsed_seconds = proc.time()[3L] - started),
    history[c("budget", "cost")]))
}
