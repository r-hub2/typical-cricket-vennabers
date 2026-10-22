# Public only as a parsnip engine entry point; use calibrated_workflow().
#' @export
fit_vennabers <- function(x, y, base_workflow, outcome, parameters, ...) {
  .validate_dots(...)
  if (inherits(parameters, "vennabers_adapter_callbacks")) return(.fit_model_adapter(x, y, parameters))
  raw <- .calibrated_fit_data(x, y, base_workflow, outcome, parameters)
  settings <- parameters$requested
  if (!is.null(settings$seed)) withr::local_seed(settings$seed)
  counter <- NULL
  if (!is.null(settings$tuning)) {
    counter <- new.env(parent = emptyenv())
    counter$n <- 0L
  }
  plan <- .model_tuning_plan(settings, counter)
  joint_pairwise <- identical(settings$model_training, "pairwise") &&
    identical(settings$tuning$strategy, "joint") && !is.null(settings$tuning$model)
  grid_joint <- joint_pairwise && .pairwise_grid_joint(settings$tuning, TRUE, settings$postprocessing)
  if (!is.null(settings$tuning$calibration) || joint_pairwise) {
    chosen <- if (grid_joint) .choose_pairwise_joint(raw, base_workflow, outcome, parameters, counter) else
      if (!is.null(settings$tuning$calibration)) .choose_calibration(raw, base_workflow, outcome, parameters, counter) else
        .assess_pairwise_joint(raw, base_workflow, outcome, parameters, counter)
    raw <- chosen$data
    parameters$requested <- chosen$settings
    plan <- .model_tuning_plan(chosen$settings, counter, chosen$pairs, chosen$pair_grid_id)
    plan$config <- settings$tuning
    plan$calibration <- if (!joint_pairwise || !is.null(settings$tuning$calibration)) chosen$diagnostics else NULL
    plan$joint <- if (joint_pairwise) chosen$diagnostics else NULL
    if (settings$mode != "regression" && any(tabulate(raw[[outcome]], nbins = nlevels(raw[[outcome]])) == 0L)) {
      stop("Final fit F must contain every class after the S/F split.", call. = FALSE)
    }
  }
  .fit_calibrated_components(base_workflow, raw, outcome, parameters, plan)
}

.calibrated_fit_data <- function(x, y, base_workflow, outcome, parameters) {
  .integration_packages()
  settings <- parameters$requested
  # A saved specification can be fitted in a session without its optional packages.
  if (!is.null(settings$tuning)) .tuning_packages()
  .postprocessing_packages(settings$postprocessing)
  regression <- settings$mode == "regression"
  .validate_base_workflow(base_workflow, outcome, settings$classification,
                          allow_tune = !is.null(settings$tuning$model), postprocessing = settings$postprocessing)
  .validate_names(names(x), "input column")
  if (is.data.frame(y) && ncol(y) == 1L) y <- y[[1L]]
  classes <- .validate_truth(y, outcome, settings$classification, settings$mode)
  if (!regression && any(tabulate(y, nbins = length(classes)) == 0L)) {
    stop("Outcome ", sQuote(outcome), " must contain all classes in the development data.", call. = FALSE)
  }
  raw <- as.data.frame(x)
  raw[[outcome]] <- y
  if (!is.null(parameters$schema)) {
    .schema_encoded_data(raw, parameters$schema, truth = TRUE)
    raw <- .bind_schema(raw, parameters$schema)
  }
  raw
}

# Both public fits and proxy candidates enter here with model selection resolved
# to local searches or fixed pair values; this fitter never selects calibration.
.fit_calibrated_components <- function(base_workflow, raw, outcome, parameters, tuning_plan) {
  settings <- parameters$requested
  regression <- settings$mode == "regression"
  multiclass <- !regression && settings$classification != "binary"
  y <- raw[[outcome]]
  classes <- levels(y)
  predictors <- setdiff(names(raw), outcome)
  pre <- workflows::extract_preprocessor(base_workflow)
  if (inherits(pre, "workflow_variables") &&
      !identical(names(tidyselect::eval_select(pre$outcomes, raw)), outcome)) {
    stop("Base workflow outcome must be ", sQuote(outcome), ".", call. = FALSE)
  }
  effective <- settings
  if (multiclass) effective$event_level <- NULL
  if (regression) effective[c("classification", "event_level", "loss", "precision")] <- NULL
  effective$predictors <- predictors
  effective$strata <- if (!regression && is.null(settings$strata)) outcome else FALSE
  if (settings$method == "cross") effective$calibration_prop <- NULL else effective$v <- NULL
  if (.inductive_sized(settings)) effective$calibration_prop <- NULL
  if (!is.null(settings$splits)) {
    effective$calibration_prop <- effective$v <- effective$strata <- NULL
  }
  parameters$effective <- effective
  event <- if (regression || multiclass) NULL else .event_class(classes, settings$event_level)
  # One local RNG scope includes the resolver, split generation and model learning.
  rng_kind <- RNGkind()
  state <- structure(list(parameters = parameters, schema = parameters$schema, outcome = outcome,
    predictors = predictors, columns = names(raw), class_levels = classes,
    event_class = event, loss = effective$loss,
    versions = c(R = as.character(getRversion()),
      vapply(unique(c("vennabers", names(.integration_versions), generics::required_pkgs(base_workflow))),
        function(p) as.character(utils::packageVersion(p)), "")),
    rng_kind = rng_kind,
    tuning = if (!is.null(tuning_plan$config)) list(config = tuning_plan$config,
      calibration = tuning_plan$calibration, models = NULL, seed = settings$seed,
      joint = tuning_plan$joint,
      classes = classes, schema = parameters$schema, method = settings$method,
      classification = settings$classification, model_training = settings$model_training,
      prediction_model = settings$prediction_model) else NULL),
    class = "vennabers_fit")
  if (!is.null(settings$postprocessing)) {
    packages <- names(.postprocessing_versions(settings$postprocessing))
    state$versions <- c(state$versions, vapply(setdiff(packages, names(state$versions)),
      function(p) as.character(utils::packageVersion(p)), ""))
  }
  if (identical(settings$model_training, "pairwise")) {
    return(.fit_pairwise(base_workflow, raw, state, tuning_plan))
  }
  domain <- list(data = raw, workflow = base_workflow, plan = settings$splits,
    schema = parameters$schema,
    outcome = outcome, strata = effective$strata, classes = classes, event = event, y = y,
    # After calibration selection, the local rows are the final-fit F.
    rows_word = if (is.null(tuning_plan$calibration)) "local rows" else "final-fit F rows")
  task <- .task_for(settings)
  result <- .fit_domains(list(domain), task, settings, tuning_plan)
  if (regression) {
    parameters$effective$m <- vapply(result$plans[[1L]], function(f) f$moderation$m, 0)
    parameters$effective$epsilon <- vapply(result$plans[[1L]], function(f) f$moderation$epsilon, 0)
  }
  if (.inductive_sized(settings)) {
    parameters$effective[c("training_size", "calibration_size", "unassigned_size")] <-
      as.list(lengths(result$plans[[1L]]$Inductive1[c("analysis", "assessment", "unassigned")]))
  }
  component <- result$fitted[[1L]]
  state$parameters <- parameters
  state$folds <- component$folds
  state$models <- component$models
  state$sources <- component$sources
  state$calibrator <- component$calibrator
  state$calibration <- component$calibration
  if (!is.null(settings$postprocessing)) {
    state$postprocessing <- component$postprocessing
  }
  if (!is.null(state$tuning)) {
    state$tuning$models <- component$model_selection
    state$tuning$sources <- component$sources
    state$tuning$fit_roles <- component$folds
  }
  if (!is.null(state$tuning)) state$tuning$cost <- list(observed_base_fits = tuning_plan$counter$n)
  state
}

# Every fit is a list of components: one shared component, or one per class
# pair with its recipe-retargeted workflow. Four phase functions, each applied
# in turn, give one path for all task and model-training combinations. Only
# calibration reads the calibration settings, so source and full fits can be
# reused across them. A domain is plain data; a pair domain differs only by its
# `pair_id`, which prefixes its error labels, seeds its split and adds the pair
# class check.
.fit_domains <- function(domains, task, settings, tuning_plan) {
  plans <- lapply(domains, function(domain) .plan_component(domain, task, settings))
  if (!is.null(settings$postprocessing)) {
    plans <- lapply(seq_along(domains), function(i) .plan_postprocessing(plans[[i]], domains[[i]], settings))
  }
  fitted <- lapply(seq_along(domains), function(i) .fit_component_sources(plans[[i]], domains[[i]], task, settings, tuning_plan))
  fitted <- lapply(fitted, .calibrate_component, task = task, settings = settings)
  if (identical(settings$prediction_model, "full")) {
    fitted <- lapply(seq_along(domains), function(i) .use_full_model(fitted[[i]],
      .fit_component_full(domains[[i]], task, settings, tuning_plan, names(fitted[[i]]$postprocessing))))
  }
  list(plans = plans, fitted = fitted)
}

# A fold or model slot ID, prefixed by the pair ID in a pair domain.
.domain_label <- function(domain, slot) paste(c(domain$pair_id, slot), collapse = " ")

.resolve_split <- function(domain, settings, local, cross) {
  plan <- domain$plan
  if (is.function(plan)) return(plan(local))
  if (!is.null(plan)) return(plan)
  outcome <- domain$outcome
  if (cross) {
    if (identical(domain$strata, FALSE)) rsample::vfold_cv(local, v = settings$v) else
      rlang::inject(rsample::vfold_cv(local, v = settings$v, strata = !!outcome, pool = 0))
  } else if (identical(domain$strata, FALSE)) {
    rsample::initial_split(local, prop = 1 - settings$calibration_prop)
  } else {
    rlang::inject(rsample::initial_split(local, prop = 1 - settings$calibration_prop,
                                         strata = !!outcome, pool = 0))
  }
}

.plan_component <- function(domain, task, settings) {
  local <- domain$data
  cross <- settings$method == "cross"
  pair <- !is.null(domain$pair_id)
  # Automatic size plans may leave rows unassigned (U), so they bypass the
  # complete-partition validator that rsample and user plans still pass through.
  sized <- !cross && is.null(domain$plan) && .inductive_sized(settings)
  # Pair split and fold errors name their pair and slot.
  labelled <- function(expr) if (!pair) expr else tryCatch(expr, error = function(e) {
    stop(.domain_label(domain, if (cross) "Cross" else "Inductive1"), ": ", conditionMessage(e), call. = FALSE)
  })
  draw <- function(f) {
    # A fixed local seed makes pair plans independent of preceding model draws.
    if (pair && !is.null(settings$seed)) withr::with_seed(settings$seed, f()) else f()
  }
  split <- if (sized) draw(function() .inductive_size_positions(domain, settings, local)) else labelled({
    if (cross && is.null(domain$plan) && settings$v > nrow(local)) {
      stop("`v` cannot exceed the number of ", domain$rows_word,
           "; the fold count is not reduced.", call. = FALSE)
    }
    resolve <- function() draw(function() .resolve_split(domain, settings, local, cross))
    if (pair) resolve() else tryCatch(resolve(), error = function(e) stop(settings$method,
      " rsample split failed (nonempty training and calibration required): ", conditionMessage(e), call. = FALSE))
  })
  labelled({
    folds <- if (sized) list(Inductive1 = split) else
      if (cross) .cross_splits(split, local) else list(Inductive1 = .inductive_split(split, local))
    for (fold_id in names(folds)) {
      folds[[fold_id]]$training_size <- length(folds[[fold_id]]$analysis)
      folds[[fold_id]]$calibration_size <- length(folds[[fold_id]]$assessment)
    }
    if (pair && !cross) .pair_fold_classes(folds, local[[domain$outcome]])
    task$folds(folds, domain$y, settings)
  })
}

.fit_component_sources <- function(folds, domain, task, settings, tuning_plan) {
  local <- domain$data
  workflow <- domain$workflow
  classes <- domain$classes
  event <- domain$event
  y <- domain$y
  models <- calibration_predictions <- truths <- vector("list", length(folds))
  model_selection <- vector("list", length(folds))
  postprocessing <- if (!is.null(settings$postprocessing)) vector("list", length(folds))
  names(models) <- names(folds)
  names(model_selection) <- names(folds)
  if (!is.null(postprocessing)) names(postprocessing) <- names(folds)
  for (id in names(folds)) {
    analysis <- folds[[id]]$analysis
    assessment <- folds[[id]]$assessment
    if (is.null(settings$postprocessing)) {
      selected <- .choose_model(workflow, local[analysis, , drop = FALSE], domain, settings, tuning_plan, id)
      model_selection[id] <- list(selected$diagnostics)
      models[[id]] <- tryCatch(.fit_base_model(selected$workflow, local[analysis, , drop = FALSE], tuning_plan$counter),
        error = function(e) stop(.domain_label(domain, id), " base workflow training: ", conditionMessage(e), call. = FALSE))
    } else {
      # Each inner B selects its own settings; P and C never enter the search.
      fitted <- .fit_postprocessing(workflow, local, folds[[id]], domain, settings, id, tuning_plan = tuning_plan)
      models[[id]] <- fitted$model
      postprocessing[[id]] <- fitted$diagnostics
      model_selection[id] <- list(fitted$selection)
    }
    k <- match(id, names(folds))
    calibration_predictions[[k]] <- task$pred(models[[id]], local[assessment, , drop = FALSE], .domain_label(domain, id), classes, event)
    truths[[k]] <- task$truth(y, assessment, classes, event)
  }
  calibration <- stats::setNames(lapply(seq_along(folds), function(k) {
    task$record(folds[[k]]$assessment, truths[[k]], calibration_predictions[[k]])
  }), names(folds))
  # `sources` maps every fold slot to the model that predicts it: itself when
  # matched, the shared Full1 refit when full. `calibration` holds the held-out
  # predictions that .calibrate_component() fits on.
  c(domain$fields, list(folds = folds, models = models, model_selection = model_selection,
    sources = stats::setNames(names(folds), names(folds)), calibration = calibration),
    if (!is.null(postprocessing)) list(postprocessing = postprocessing))
}

.calibrate_component <- function(component, task, settings) {
  component$calibrator <- task$calibrate(component$calibration, .va_fold(component$folds, settings), settings)
  component
}

# The full refit on all local rows. `taken` are the source diagnostics IDs the
# full ID must not overwrite; it reads no calibration settings.
.fit_component_full <- function(domain, task, settings, tuning_plan, taken = NULL) {
  local <- domain$data
  # Explicit VA fold names must not be overwritten by the final diagnostics.
  id <- utils::tail(make.unique(c(taken, "Full1")), 1L)
  diagnostics <- NULL
  if (is.null(settings$postprocessing)) {
    selected <- .choose_model(domain$workflow, local, domain, settings, tuning_plan, id)
    model <- tryCatch(.fit_base_model(selected$workflow, local, tuning_plan$counter),
      error = function(e) stop(.domain_label(domain, id), " base workflow training: ", conditionMessage(e), call. = FALSE))
    selection <- selected$diagnostics
  } else {
    # Resolve F only after all held-out sources: its RNG draws cannot alter C.
    plan <- .plan_postprocessing(stats::setNames(list(list(analysis = seq_len(nrow(local)), assessment = integer())), id), domain, settings)
    fitted <- .fit_postprocessing(domain$workflow, local, plan[[id]], domain, settings, id, tuning_plan = tuning_plan)
    model <- fitted$model
    selection <- fitted$selection
    diagnostics <- fitted$diagnostics
  }
  # Check the prediction contract without using full-model values for calibration.
  task$pred(model, local[1L, , drop = FALSE], .domain_label(domain, id), domain$classes, domain$event)
  list(id = id, model = model, selection = selection, diagnostics = diagnostics)
}

# The full model replaces the temporary source models and predicts every fold
# slot; called only after the complete full fit succeeded.
.use_full_model <- function(component, full) {
  if (!is.null(full$diagnostics)) component$postprocessing[[full$id]] <- full$diagnostics
  component$models <- stats::setNames(list(full$model), full$id)
  component$model_selection[full$id] <- list(full$selection)
  component$sources[] <- full$id
  component
}

# The only binary/multiclass/regression switch. Each task owns its fold checks,
# held-out truth and base predictions, calibrator fitting and prediction, and
# the calibration record shape. `folds` runs on every planned component before
# the first recipe/model fit.
.task_for <- function(settings) {
  if (settings$mode == "regression") .regression_task() else
    if (settings$classification != "binary") .multiclass_task() else .binary_task()
}

# Fold IDs of the pooled held-out rows for a cross core; NULL selects the
# inductive core in .va().
.va_fold <- function(folds, settings) {
  if (settings$method == "cross") rep(names(folds), lengths(lapply(folds, `[[`, "assessment")))
}

.va <- function(truth, p, fold = NULL, ..., inductive = venn_abers, cross = cross_venn_abers) {
  if (is.null(fold)) inductive(truth, p, ...) else cross(truth, p, fold, ...)
}

# A cross core scores the sources of all its folds together; an inductive core
# scores its single source.
.va_newdata <- function(calibrator, points, stack) {
  if (is.null(calibrator$folds)) points[[1L]] else stack(points)
}

# Per-fold binary or regression calibrators, named like the sources.
.va_fold_calibrators <- function(calibrator, ids) {
  if (is.null(calibrator$folds)) stats::setNames(list(calibrator), ids) else calibrator$calibrators
}

# One long detail row per new row and fold, in fold order.
.fold_details <- function(calibrators, points, ...) {
  do.call(rbind, lapply(names(calibrators), function(id) {
    n <- length(points[[id]])
    data.frame(row = seq_len(n), fold = rep(id, n), stats::predict(calibrators[[id]], points[[id]], ...))
  }))
}

.binary_task <- function() list(
  folds = function(folds, y, settings) folds,
  pred = function(model, data, id, classes, event) .base_probability(model, data, classes, event, id),
  empty = function(classes) numeric(),
  truth = function(y, assessment, classes, event) y[assessment] == event,
  calibrate = function(records, fold, settings) {
    .va(unlist(lapply(records, `[[`, "truth"), use.names = FALSE),
        unlist(lapply(records, `[[`, "probability"), use.names = FALSE), fold, precision = settings$precision)
  },
  predict = function(calibrator, points, details, classes, event, loss) {
    probability <- stats::predict(calibrator, .va_newdata(calibrator, points, function(p) do.call(cbind, p)),
                                  loss = loss)$probability
    probabilities <- lapply(classes, function(cl) if (cl == event) probability else 1 - probability)
    names(probabilities) <- paste0(".pred_", classes)
    result <- list(prediction = tibble::as_tibble(probabilities, .name_repair = "minimal"))
    if (!details) return(result)
    folds <- .fold_details(.va_fold_calibrators(calibrator, names(points)), points, loss = loss)
    folds$event_class <- rep(event, nrow(folds))
    result$details <- folds[c("row", "fold", "event_class", "lower", "upper", "probability")]
    result
  },
  record = function(rows, truth, prediction) data.frame(row = rows, truth = truth, probability = prediction))

.multiclass_task <- function() list(
  folds = function(folds, y, settings) {
    held <- lapply(folds, `[[`, "assessment")
    fold <- if (settings$method == "cross") rep(names(folds), lengths(held)) else NULL
    .multiclass_eligibility(y[unlist(held, use.names = FALSE)], settings$classification, fold)
    folds
  },
  pred = function(model, data, id, classes, event) .base_probability(model, data, classes, event, id),
  empty = function(classes) matrix(numeric(), 0L, length(classes), dimnames = list(NULL, classes)),
  truth = function(y, assessment, classes, event) y[assessment],
  calibrate = function(records, fold, settings) {
    .va(do.call(c, lapply(records, `[[`, "truth")), do.call(rbind, lapply(records, `[[`, "probability")), fold,
        precision = settings$precision, method = settings$classification)
  },
  predict = function(calibrator, points, details, classes, event, loss) {
    result <- stats::predict(calibrator, .va_newdata(calibrator, points, identity), loss = loss, type = "details")
    prediction <- result$prediction
    colnames(prediction) <- paste0(".pred_", classes)
    prediction <- tibble::as_tibble(prediction, .name_repair = "minimal")
    if (!details) return(list(prediction = prediction))
    method <- if (identical(calibrator$method, "ovr")) "ovr" else "ovo"
    details <- if (method == "ovo") result$pairs else result$classes
    names(details)[names(details) == "non_event"] <- "class_1"
    names(details)[names(details) == "event"] <- if (method == "ovo") "class_2" else "class"
    # Inductive details carry no fold column; their fold is the single source.
    if (is.null(details$fold)) details$fold <- rep(names(points), nrow(details))
    details$method <- rep(method, nrow(details))
    details$event_class <- if (method == "ovo") details$class_2 else details$class
    identity <- if (method == "ovo") c("class_1", "class_2") else "class"
    list(prediction = prediction,
         details = details[c("row", "fold", "method", "event_class", identity, "lower", "upper", "probability")])
  },
  record = function(rows, truth, prediction) list(row = rows, truth = truth, probability = prediction))

.regression_task <- function() list(
  folds = function(folds, y, settings) {
    for (id in names(folds)) {
      k <- length(folds[[id]]$assessment)
      folds[[id]]$moderation <- tryCatch(.regression_moderation(k, settings$m, settings$epsilon),
        error = function(e) stop(id, " calibration size ", k, ": ", conditionMessage(e), call. = FALSE))
    }
    folds
  },
  pred = function(model, data, id, classes, event) .base_numeric(model, data, id),
  empty = function(classes) numeric(),
  truth = function(y, assessment, classes, event) y[assessment],
  calibrate = function(records, fold, settings) {
    .va(unlist(lapply(records, `[[`, "truth"), use.names = FALSE),
        unlist(lapply(records, `[[`, "prediction"), use.names = FALSE), fold,
        m = settings$m, epsilon = settings$epsilon,
        inductive = .venn_abers_regression, cross = .cross_venn_abers_regression)
  },
  predict = function(calibrator, points, details, classes, event, loss) {
    calibrated <- stats::predict(calibrator, .va_newdata(calibrator, points, function(p) do.call(cbind, p)))
    result <- list(prediction = tibble::tibble(.pred = as.double(calibrated$prediction)))
    if (!details) return(result)
    result$details <- .fold_details(.va_fold_calibrators(calibrator, names(points)), points)
    # Only a cross core aggregates several folds into one interval per row.
    if (!is.null(calibrator$folds)) {
      result$aggregate <- tibble::tibble(row = seq_len(nrow(calibrated)), lower = calibrated$lower,
                                         upper = calibrated$upper)
    }
    result
  },
  record = function(rows, truth, prediction) data.frame(row = rows, truth = truth, prediction = prediction))
