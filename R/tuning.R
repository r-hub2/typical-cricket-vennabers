.tune_id <- function(x) {
  if (!rlang::is_call(x, "tune")) return(NULL)
  if (length(x) != 2L || !is.character(x[[2L]]) || length(x[[2L]]) != 1L ||
      is.na(x[[2L]]) || !nzchar(x[[2L]])) {
    stop("Calibration tune() placeholders need one explicit nonempty ID.", call. = FALSE)
  }
  x[[2L]]
}

.fit_base_model <- function(workflow, data, counter) {
  if (!is.null(counter)) counter$n <- counter$n + 1L
  generics::fit(workflow, data)
}

.tuning_names <- function(x, allowed, label) {
  if (!is.list(x) || is.null(names(x)) || anyNA(names(x)) ||
      any(!nzchar(names(x))) || anyDuplicated(names(x)) ||
      !all(names(x) %in% allowed)) {
    stop("`", label, "` must be a named list using only: ",
         paste(allowed, collapse = ", "), ".", call. = FALSE)
  }
}

.validate_tuning <- function(tuning, workflow, mode, model_training,
                             loss, precision, seed, m, epsilon, postprocessing = NULL) {
  model_ids <- generics::tune_args(workflow)$id
  calibration_ids <- if (mode == "regression") c(m = .tune_id(m), epsilon = .tune_id(epsilon)) else
    c(loss = .tune_id(loss), precision = .tune_id(precision))
  if (is.null(tuning)) {
    if (length(calibration_ids)) stop("Calibration tune() placeholders require `tuning$calibration`.", call. = FALSE)
    return(NULL)
  }
  .tuning_names(tuning, c("strategy", "model", "calibration", "selection_prop", "selection_split",
                          "joint_resamples"), "tuning")
  .validate_choice(tuning$strategy, c("joint", "two_stage"), "tuning$strategy")
  if (anyNA(model_ids) || any(!nzchar(model_ids)) || anyDuplicated(model_ids) ||
      anyDuplicated(calibration_ids) || any(model_ids %in% calibration_ids)) {
    stop("Model/recipe and calibration tune() IDs must be nonempty and unique.", call. = FALSE)
  }
  if (xor(length(model_ids) > 0L, !is.null(tuning$model))) {
    stop("`tuning$model` is required exactly when the base workflow has tune() placeholders.", call. = FALSE)
  }
  if (xor(length(calibration_ids) > 0L, !is.null(tuning$calibration))) {
    stop("`tuning$calibration` is required exactly when a calibration option uses tune().", call. = FALSE)
  }
  if (!length(model_ids) && !length(calibration_ids)) {
    stop("`tuning` needs at least one model/recipe or calibration tune() placeholder.", call. = FALSE)
  }
  pair_joint <- identical(model_training, "pairwise") && identical(tuning$strategy, "joint") && length(model_ids)
  pair_grid_joint <- .pairwise_grid_joint(tuning, identical(model_training, "pairwise"), postprocessing)
  if (!is.null(tuning$joint_resamples) &&
      (!pair_joint || length(calibration_ids) ||
       !(is.function(tuning$joint_resamples) || inherits(tuning$joint_resamples, "rset")))) {
    stop("`tuning$joint_resamples` requires model-only pairwise joint tuning and a global rset or provider function.", call. = FALSE)
  }
  if (pair_joint && !length(calibration_ids) && is.null(tuning$joint_resamples)) {
    stop("Pairwise joint model tuning without calibration options requires `tuning$joint_resamples` on global data.", call. = FALSE)
  }
  if (length(calibration_ids)) {
    if (is.null(tuning$selection_prop) == is.null(tuning$selection_split)) {
      stop("Supply exactly one of `selection_prop` and `selection_split` for calibration selection.", call. = FALSE)
    }
    if (!is.null(tuning$selection_prop) &&
        !.is_proportion(tuning$selection_prop)) {
      stop("`selection_prop` must be one number strictly between 0 and 1.", call. = FALSE)
    }
    if (!is.null(tuning$selection_split) &&
        !inherits(tuning$selection_split, "rsplit") && !is.function(tuning$selection_split)) {
      stop("`selection_split` must be an rsplit or function of local data.", call. = FALSE)
    }
  } else if (!is.null(tuning$selection_prop) || !is.null(tuning$selection_split)) {
    stop("A selection S/F split requires automatic calibration selection.", call. = FALSE)
  }
  if (length(model_ids)) tuning$model <- .validate_search(tuning$model, model_ids, workflow, seed,
    "model", require_resamples = !pair_grid_joint)
  if (length(calibration_ids)) tuning$calibration <- .validate_search(
    tuning$calibration, calibration_ids, workflow, seed, "calibration")
  if (pair_grid_joint && !is.null(tuning$model$resamples)) {
    stop("Pairwise joint model grids use global calibration resamples or `joint_resamples`; omit `tuning$model$resamples`.", call. = FALSE)
  }
  if (pair_grid_joint && length(calibration_ids) &&
      !identical(tuning$model$primary_metric, tuning$calibration$primary_metric)) {
    stop("Pairwise joint model and calibration grids need the same primary metric.", call. = FALSE)
  }
  if (pair_grid_joint && length(calibration_ids) && !is.null(tuning$model$control) &&
      !is.null(tuning$calibration$control) &&
      !identical(tuning$model$control, tuning$calibration$control)) {
    stop("Pairwise joint model and calibration grids need one consistent control.", call. = FALSE)
  }
  tuning
}

.validate_search <- function(search, ids, workflow, seed, label, require_resamples = TRUE) {
  if (!is.list(search) || is.null(search$method)) {
    stop("`tuning", "$", label, "$method` is required.", call. = FALSE)
  }
  if (identical(search$method, "grid")) {
    return(.validate_grid_search(search, ids, workflow, seed, label, require_resamples))
  }
  .validate_bayes_search(search, ids, workflow, label, require_resamples)
}

.validate_bayes_search <- function(search, ids, workflow, label, require_resamples = TRUE) {
  .tuning_names(search, c("method", "initial", "iter", "param_info", "resamples", "primary_metric",
                          "report_metrics", "control"), paste0("tuning$", label))
  .validate_choice(search$method, "bayes", paste0("tuning$", label, "$method"))
  for (name in c("initial", "iter")) {
    value <- search[[name]]
    if (!.is_count(value)) {
      stop("`tuning$", label, "$", name, "` must be one positive whole number.", call. = FALSE)
    }
    search[[name]] <- as.integer(value)
  }
  if ((require_resamples && is.null(search$resamples)) ||
      (!is.null(search$resamples) &&
       !(is.function(search$resamples) || inherits(search$resamples, "rset")))) {
    stop("`tuning$", label, "$resamples` must be an rset or function of local data.", call. = FALSE)
  }
  if (!inherits(search$param_info, "parameters") || !setequal(search$param_info$id, unname(ids))) {
    stop("Bayes `tuning$", label, "$param_info` must describe exactly the tune() IDs.", call. = FALSE)
  }
  if (search$initial <= nrow(search$param_info)) {
    stop("Bayes `tuning$", label, "$initial` must exceed the number of tuned parameters.", call. = FALSE)
  }
  if (is.null(search$control)) search$control <- tune::control_bayes()
  control <- search$control
  if (!inherits(control, "control_bayes") || !isTRUE(control$save_pred) ||
      !identical(control$uncertain, Inf) || !.is_count(control$seed, 0, Inf) ||
      !is.numeric(control$no_improve) || length(control$no_improve) != 1L || control$no_improve < 1 ||
      (!is.na(control$time_limit) && (!is.numeric(control$time_limit) || length(control$time_limit) != 1L ||
                                      !is.finite(control$time_limit) || control$time_limit <= 0))) {
    stop("Bayes control needs save_pred = TRUE, uncertain = Inf, a whole seed, positive no_improve, and a positive or NA time_limit.", call. = FALSE)
  }
  if (!is.null(control$extract)) stop("Custom `control_bayes(extract = ...)` conflicts with selection diagnostics.", call. = FALSE)
  search$control$allow_par <- FALSE
  if (is.null(search$primary_metric)) search$primary_metric <-
    if (workflows::extract_spec_parsnip(workflow)$mode == "regression") yardstick::rmse else yardstick::mn_log_loss
  if (!inherits(search$primary_metric, "metric")) {
    stop("`tuning$", label, "$primary_metric` must be one yardstick metric.", call. = FALSE)
  }
  reports <- search$report_metrics
  if (inherits(reports, "metric_set")) reports <- attr(reports, "metrics")
  if (inherits(reports, "metric")) reports <- list(reports)
  if (is.null(reports)) reports <- list()
  if (!is.list(reports) || !all(vapply(reports, inherits, logical(1), "metric"))) {
    stop("`tuning$", label, "$report_metrics` must contain yardstick metrics.", call. = FALSE)
  }
  if (length(reports) && is.null(names(reports))) names(reports) <- paste0("report_", seq_along(reports))
  search$report_metrics <- reports
  search$ids <- ids
  search
}

.validate_grid_search <- function(search, ids, workflow, seed, label, require_resamples = TRUE) {
  .tuning_names(search, c("method", "grid", "param_info", "resamples", "primary_metric",
                          "report_metrics", "control"), paste0("tuning$", label))
  .validate_choice(search$method, "grid", paste0("tuning$", label, "$method"))
  if (is.null(search$grid)) stop("`tuning$", label, "$grid` is required.", call. = FALSE)
  if ((require_resamples && is.null(search$resamples)) ||
      (!is.null(search$resamples) &&
       !(is.function(search$resamples) || inherits(search$resamples, "rset")))) {
    stop("`tuning$", label, "$resamples` must be an rset or function of local data.", call. = FALSE)
  }
  if (!is.null(search$control) &&
      (!inherits(search$control, "control_grid") || !isTRUE(search$control$save_pred))) {
    stop("`tuning$", label, "$control` must be control_grid(save_pred = TRUE).", call. = FALSE)
  }
  if (!is.null(search$control$extract)) {
    stop("Custom `control_grid(extract = ...)` conflicts with selection diagnostics.", call. = FALSE)
  }
  if (is.null(search$primary_metric)) search$primary_metric <-
    if (workflows::extract_spec_parsnip(workflow)$mode == "regression") yardstick::rmse else yardstick::mn_log_loss
  if (!inherits(search$primary_metric, "metric")) {
    stop("`tuning$", label, "$primary_metric` must be one yardstick metric.", call. = FALSE)
  }
  reports <- search$report_metrics
  if (inherits(reports, "metric_set")) reports <- attr(reports, "metrics")
  if (inherits(reports, "metric")) reports <- list(reports)
  if (is.null(reports)) reports <- list()
  if (!is.list(reports) || !all(vapply(reports, inherits, logical(1), "metric"))) {
    stop("`tuning$", label, "$report_metrics` must contain yardstick metrics.", call. = FALSE)
  }
  if (length(reports) && is.null(names(reports))) {
    names(reports) <- paste0("report_", seq_along(reports))
  }
  search$report_metrics <- reports
  grid <- search$grid
  if (.is_count(grid)) {
    params <- search$param_info
    if (is.null(params) && label == "model") params <- tune::extract_parameter_set_dials(workflow)
    if (!inherits(params, "parameters") || !setequal(params$id, unname(ids))) {
      stop("An integer `tuning$", label, "$grid` needs dials `param_info`.", call. = FALSE)
    }
    search$param_info <- params
    search$grid_size <- as.integer(grid)
    search$grid <- NULL
    search$seed <- seed
    search$ids <- ids
    return(search)
  }
  search$grid <- .checked_grid(grid, ids, label)
  search$ids <- ids
  search$candidate_ids <- sprintf("%s_%03d", label, seq_len(nrow(search$grid)))
  search
}

.checked_grid <- function(grid, ids, label) {
  if (!is.data.frame(grid) || !nrow(grid) || anyNA(names(grid)) ||
      anyDuplicated(names(grid)) || !setequal(names(grid), unname(ids))) {
    stop("`tuning$", label, "$grid` needs rows with exactly these ID columns: ",
         paste(unname(ids), collapse = ", "), ".", call. = FALSE)
  }
  grid <- as.data.frame(grid[unname(ids)], stringsAsFactors = FALSE)
  if (any(vapply(grid, function(x) anyNA(x), logical(1)))) {
    stop("`tuning$", label, "$grid` cannot contain missing values.", call. = FALSE)
  }
  if (label == "calibration") {
    for (option in intersect(c("m", "epsilon"), names(ids))) {
      values <- grid[[ids[[option]]]]
      # Each row is one candidate; a column that is not a plain numeric vector fails whole.
      for (value in if (is.numeric(values) && is.vector(values)) values else list(values)) {
        .validate_moderation(if (option == "m") value, if (option == "epsilon") value,
                             label = "Calibration grid: ")
      }
    }
    if ("loss" %in% names(ids)) {
      value <- as.character(grid[[ids[["loss"]]]])
      if (!all(value %in% c("log", "brier"))) stop("Calibration loss grid allows only log and brier.", call. = FALSE)
      grid[[ids[["loss"]]]] <- value
    }
    if ("precision" %in% names(ids)) {
      value <- as.character(grid[[ids[["precision"]]]])
      if (!all(value %in% c("none", "0", "2", "-2"))) {
        stop("Calibration precision grid allows only none, 0, 2 and -2.", call. = FALSE)
      }
      grid[[ids[["precision"]]]] <- value
    }
  }
  grid
}

# A resolved plan contains only the model search (if any), fixed pair choices,
# diagnostics and the shared fit counter. Nested candidates never repeat global
# calibration or pairwise joint selection.
.model_tuning_plan <- function(settings, counter, pairs = NULL, pair_grid_id = NULL) {
  list(model = if (is.null(pairs)) settings$tuning$model else NULL, strategy = settings$tuning$strategy,
    pairs = pairs, pair_grid_id = pair_grid_id, counter = counter,
    config = settings$tuning)
}

.calibrated_candidate <- function(workflow, settings, schema, counter,
                                   pairs = NULL, pair_grid_id = NULL) {
  if (!is.function(settings$splits)) settings$splits <- NULL
  list(workflow = workflow, parameters = list(requested = settings, schema = schema),
    tuning_plan = .model_tuning_plan(settings, counter, pairs, pair_grid_id))
}

.choose_model <- function(workflow, data, domain, settings, tuning_plan, path_id = NULL) {
  if (!is.null(tuning_plan$pairs)) {
    scope <- if (identical(path_id, "Full1")) "full" else "source"
    values <- tuning_plan$pairs[[domain$pair_id]][[scope]]
    return(list(workflow = tune::finalize_workflow(workflow, values),
      diagnostics = list(chosen = tuning_plan$pair_grid_id, selected = values,
                         pair_id = domain$pair_id, fold_id = path_id)))
  }
  search <- tuning_plan$model
  if (is.null(search)) return(list(workflow = workflow, diagnostics = NULL))
  outcome <- domain$outcome
  select <- if (identical(search$method, "bayes")) .bayes_select else .grid_select
  candidate <- if (tuning_plan$strategy == "two_stage") {
    # The first two-stage step scores raw engine/recipe output: no Tailor, no VA.
    # Stripped once per search: finding an absent postprocessor costs a caught error.
    unprocessed <- .without_postprocessor(workflow)
    function(values) list(workflow = tune::finalize_workflow(unprocessed, values))
  } else {
    function(values) {
      nested <- settings
      nested$tuning <- NULL
      .calibrated_candidate(tune::finalize_workflow(workflow, values), nested, domain$schema, tuning_plan$counter)
    }
  }
  result <- select(search, data, outcome, settings$event_level, settings$classification,
                   candidate, "model", tuning_plan$counter)
  result$diagnostics$selected <- result$values
  result$diagnostics$seed <- settings$seed
  list(workflow = tune::finalize_workflow(workflow, result$values),
       diagnostics = result$diagnostics)
}

.without_postprocessor <- function(workflow) {
  if (is.null(tryCatch(workflows::extract_postprocessor(workflow), error = function(e) NULL))) return(workflow)
  workflows::remove_tailor(workflow)
}

.calibration_values <- function(settings, ids, values) {
  for (option in intersect(c("m", "epsilon"), names(ids))) settings[[option]] <- values[[ids[[option]]]]
  if ("loss" %in% names(ids)) settings$loss <- values[[ids[["loss"]]]]
  if ("precision" %in% names(ids)) {
    precision <- values[[ids[["precision"]]]]
    settings$precision <- if (identical(precision, "none")) NULL else as.numeric(precision)
  }
  settings
}

.selection_roles <- function(data, outcome, settings, config) {
  split <- if (!is.null(config$selection_prop)) {
    if (settings$mode == "regression" || identical(settings$strata, FALSE)) rsample::initial_split(data, prop = 1 - config$selection_prop) else
      rlang::inject(rsample::initial_split(data, prop = 1 - config$selection_prop,
                                           strata = !!outcome, pool = 0))
  } else if (is.function(config$selection_split)) config$selection_split(data) else config$selection_split
  .split_positions(split, data, "selection S/F")
}

.choose_calibration <- function(data, workflow, outcome, parameters, counter) {
  settings <- parameters$requested
  config <- settings$tuning
  roles <- .selection_roles(data, outcome, settings, config)
  selection <- data[roles$assessment, , drop = FALSE]
  search <- config$calibration
  select <- if (identical(search$method, "bayes")) .bayes_select else .grid_select
  result <- select(search, selection, outcome, settings$event_level, settings$classification,
    function(values) {
      nested <- .calibration_values(settings, search$ids, values)
      nested$tuning$calibration <- NULL
      .calibrated_candidate(workflow, nested, parameters$schema, counter)
    }, "calibration", counter)
  settings <- .calibration_values(settings, search$ids, result$values)
  settings$tuning$calibration <- NULL
  result$diagnostics$selected <- result$values
  result$diagnostics$seed <- settings$seed
  result$diagnostics$roles <- roles
  result$diagnostics$selection_rows <- roles$assessment
  result$diagnostics$final_rows <- roles$analysis
  list(data = data[roles$analysis, , drop = FALSE], settings = settings,
       diagnostics = result$diagnostics)
}

#' Inspect classification or regression grid selection after fitting
#' @param object A fitted calibrated workflow.
#' @return Search configuration, candidate scores, selected values, plans, and data roles.
#' @export
venn_abers_tuning <- function(object) {
  state <- if (inherits(object, "vennabers_fit")) object else workflows::extract_fit_engine(object)
  if (!inherits(state, "vennabers_fit") || is.null(state$tuning)) {
    stop("`object` must be a fitted tuned calibrated workflow.", call. = FALSE)
  }
  state$tuning
}
