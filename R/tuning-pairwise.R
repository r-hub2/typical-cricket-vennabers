# Global pairwise joint selection builds candidates for the common search path.

# With postprocessing every B selects locally, as in pairwise Bayes; the global
# Cartesian pair grid would choose outside B and share one choice across B sets.
.pairwise_grid_joint <- function(tuning, pairwise, postprocessing = NULL) {
  isTRUE(pairwise) && is.null(postprocessing) && identical(tuning$strategy, "joint") &&
    !is.null(tuning$model) && identical(tuning$model$method, "grid") &&
    (is.null(tuning$calibration) || identical(tuning$calibration$method, "grid"))
}

.assess_pairwise_joint <- function(data, workflow, outcome, parameters, counter) {
  settings <- parameters$requested
  search <- list(grid = data.frame(pairwise = 1L), candidate_ids = "joint_0001",
    ids = c(pairwise = "pairwise"), resamples = settings$tuning$joint_resamples,
    primary_metric = settings$tuning$model$primary_metric,
    report_metrics = settings$tuning$model$report_metrics)
  result <- .grid_select(search, data, outcome, settings$event_level, settings$classification,
    function(values) {
      .calibrated_candidate(workflow, settings, parameters$schema, counter)
    }, "joint", counter)
  result$diagnostics$pair_ids <- paste0("pair_", seq_len(choose(nlevels(data[[outcome]]), 2L)))
  result$diagnostics$seed <- settings$seed
  list(data = data, settings = settings, diagnostics = result$diagnostics)
}

.choose_pairwise_joint <- function(data, workflow, outcome, parameters, counter) {
  settings <- parameters$requested
  config <- settings$tuning
  roles <- NULL
  selection <- data
  if (!is.null(config$calibration)) {
    roles <- .selection_roles(data, outcome, settings, config)
    selection <- data[roles$assessment, , drop = FALSE]
  }
  pairs <- utils::combn(levels(data[[outcome]]), 2L, simplify = FALSE)
  pair_ids <- paste0("pair_", seq_along(pairs))
  model <- config$model
  grids <- stats::setNames(lapply(seq_along(pairs), function(i) {
    if (is.null(model$grid_size)) return(model$grid)
    local <- selection[selection[[outcome]] %in% pairs[[i]], , drop = FALSE]
    local[[outcome]] <- factor(local[[outcome]], levels = pairs[[i]])
    params <- dials::finalize(model$param_info, local)
    generated <- if (is.null(settings$seed)) dials::grid_space_filling(params, size = model$grid_size) else
      withr::with_seed(settings$seed, dials::grid_space_filling(params, size = model$grid_size))
    .checked_grid(generated, model$ids, "model")
  }), pair_ids)
  dimensions <- list()
  for (id in pair_ids) {
    dimensions[[paste0("..", id, "_source")]] <- seq_len(nrow(grids[[id]]))
    if (identical(settings$prediction_model, "full"))
      dimensions[[paste0("..", id, "_full")]] <- seq_len(nrow(grids[[id]]))
  }
  calibration <- config$calibration
  if (!is.null(calibration)) {
    if (!is.null(calibration$grid_size)) {
      params <- dials::finalize(calibration$param_info, selection)
      generated <- if (is.null(settings$seed)) dials::grid_space_filling(params, size = calibration$grid_size) else
        withr::with_seed(settings$seed, dials::grid_space_filling(params, size = calibration$grid_size))
      calibration$grid <- .checked_grid(generated, calibration$ids, "calibration")
    }
    dimensions[["..calibration"]] <- seq_len(nrow(calibration$grid))
  }
  if (!is.null(calibration) && any(names(calibration$grid) %in% c(names(dimensions), "..candidate_id"))) {
    stop("Calibration IDs conflict with reserved pairwise joint grid columns.", call. = FALSE)
  }
  count <- prod(vapply(dimensions, length, 0L))
  # Coupling and scoring still enumerate the Cartesian grid, even with cached fits.
  if (!is.finite(count) || count > 4096L) {
    stop("Pairwise joint Cartesian grid exceeds 4096 combinations; reduce the grids.", call. = FALSE)
  }
  combinations <- expand.grid(dimensions, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  combinations[["..candidate_id"]] <- sprintf("joint_%04d", seq_len(nrow(combinations)))
  pair_values <- function(values) stats::setNames(lapply(pair_ids, function(id) {
    source <- as.list(grids[[id]][values[[paste0("..", id, "_source")]], , drop = FALSE])
    full <- if (identical(settings$prediction_model, "full"))
      as.list(grids[[id]][values[[paste0("..", id, "_full")]], , drop = FALSE]) else source
    list(source = source, full = full)
  }), pair_ids)
  search <- if (is.null(calibration)) model else calibration
  if (!is.null(calibration)) {
    search$control <- if (is.null(calibration$control)) model$control else calibration$control
    search$report_metrics <- c(
      if (length(calibration$report_metrics)) stats::setNames(calibration$report_metrics,
        paste0("calibration_", names(calibration$report_metrics))) else list(),
      if (length(model$report_metrics)) stats::setNames(model$report_metrics,
        paste0("model_", names(model$report_metrics))) else list())
  }
  search$grid <- combinations
  search$grid_size <- NULL
  search$candidate_ids <- combinations[["..candidate_id"]]
  search$resamples <- if (is.null(calibration)) config$joint_resamples else calibration$resamples
  cache <- new.env(parent = emptyenv())
  cache$seed <- if (is.null(settings$seed)) sample.int(.Machine$integer.max, 1L) else settings$seed
  cache$reused <- 0L
  # Native results retain the candidate closure; release temporary models on exit.
  on.exit({ cache$data <- cache$entries <- cache$plans <- NULL }, add = TRUE)
  result <- .grid_select(search, selection, outcome, settings$event_level, settings$classification,
    function(values) {
      nested <- settings
      if (!is.null(calibration)) nested <- .calibration_values(nested, calibration$ids,
        as.list(calibration$grid[values[["..calibration"]], , drop = FALSE]))
      nested$tuning$calibration <- NULL
      candidate <- .calibrated_candidate(workflow, nested, parameters$schema, counter,
        pairs = pair_values(values), pair_grid_id = values[["..candidate_id"]])
      candidate$tuning_plan$pair_cache <- cache
      candidate$tuning_plan$pair_indices <- values
      candidate
    }, "joint", counter)
  result$diagnostics$budget$reused_base_fits <- cache$reused
  values <- result$values
  if (!is.null(calibration)) {
    settings <- .calibration_values(settings, calibration$ids,
      as.list(calibration$grid[values[["..calibration"]], , drop = FALSE]))
    settings$tuning$calibration <- NULL
  }
  selected_pairs <- pair_values(values)
  result$diagnostics$selected <- list(pairs = selected_pairs,
    calibration = if (!is.null(calibration)) list(loss = settings$loss, precision = settings$precision) else NULL)
  result$diagnostics$seed <- settings$seed
  result$diagnostics$pair_ids <- pair_ids
  if (!is.null(roles)) {
    result$diagnostics$roles <- roles
    result$diagnostics$selection_rows <- roles$assessment
    result$diagnostics$final_rows <- roles$analysis
  }
  list(data = if (is.null(roles)) data else data[roles$analysis, , drop = FALSE],
       settings = settings, pairs = selected_pairs, pair_grid_id = result$diagnostics$chosen,
       diagnostics = result$diagnostics)
}

# The native grid runs serially by resample. Retain only the active data's pair
# fits; they read no calibration settings, so every candidate calibrates and
# couples its own. Cached models memoise their last prediction, which the
# candidates sharing them repeat on the same assessment rows.
.fit_pairwise_grid_domains <- function(data, domains, task, settings, tuning_plan) {
  cache <- tuning_plan$pair_cache
  if (!identical(cache$data, data)) {
    cache$entries <- new.env(parent = emptyenv())
    cache$plans <- withr::with_seed(cache$seed,
      lapply(domains, function(domain) .plan_component(domain, task, settings)))
    cache$data <- data
  }
  get_fit <- function(domain, role, fit) {
    index <- tuning_plan$pair_indices[[paste0("..", domain$pair_id, "_", role)]]
    key <- paste(domain$pair_id, role, index, sep = "/")
    entry <- cache$entries[[key]]
    if (is.null(entry)) {
      before <- tuning_plan$counter$n
      # Each pair/role has its own RNG scope, independent of other grid choices.
      value <- tryCatch(withr::with_seed(cache$seed, fit), error = identity)
      entry <- list(value = value, fits = tuning_plan$counter$n - before)
      cache$entries[[key]] <- entry
    } else cache$reused <- cache$reused + entry$fits
    if (inherits(entry$value, "error")) stop(entry$value)
    entry$value
  }
  fitted <- lapply(seq_along(domains), function(i) {
    component <- get_fit(domains[[i]], "source", {
      sources <- .fit_component_sources(cache$plans[[i]], domains[[i]], task, settings, tuning_plan)
      sources$models <- lapply(sources$models, .memoised_model)
      sources
    })
    .calibrate_component(component, task, settings)
  })
  if (identical(settings$prediction_model, "full")) {
    fitted <- lapply(seq_along(domains), function(i) .use_full_model(fitted[[i]], get_fit(domains[[i]], "full", {
      full <- .fit_component_full(domains[[i]], task, settings, tuning_plan, names(fitted[[i]]$postprocessing))
      full$model <- .memoised_model(full$model)
      full
    })))
  }
  fitted <- lapply(fitted, function(component) {
    for (id in names(component$model_selection))
      component$model_selection[[id]]$chosen <- tuning_plan$pair_grid_id
    component
  })
  list(plans = cache$plans, fitted = fitted)
}

.memoised_model <- function(model) {
  structure(model, class = c("vennabers_memoised", class(model)), memo = new.env(parent = emptyenv()))
}

#' @export
predict.vennabers_memoised <- function(object, new_data, ...) {
  memo <- attr(object, "memo")
  key <- list(new_data, ...)
  if (!identical(memo$key, key)) {
    memo$value <- NextMethod()
    memo$key <- key
  }
  memo$value
}
