# The outer workflow only carries raw variables. Learning belongs to the base
# workflow after the local split, including any recipe preparation.
#' @export
calibrated_workflow <- function(workflow, outcome, ..., predictors = NULL,
                                method = "inductive", classification = "binary",
                                event_level = "first", calibration_prop = 0.2,
                                training_size = NULL, calibration_size = NULL, v = 5L,
                                strata = NULL, splits = NULL, seed = NULL, precision = NULL,
                                loss = "log", m = NULL, epsilon = NULL, schema = NULL,
                                prediction_model = "matched", model_training = "shared",
                                tuning = NULL, postprocessing = NULL) {
  .validate_dots(...)
  .exact_argument_names()
  explicit <- c(predictors = !missing(predictors), method = !missing(method),
                classification = !missing(classification), event_level = !missing(event_level),
                calibration_prop = !missing(calibration_prop), v = !missing(v), strata = !missing(strata),
                splits = !missing(splits), seed = !missing(seed), precision = !missing(precision), loss = !missing(loss),
                m = !missing(m), epsilon = !missing(epsilon), prediction_model = !missing(prediction_model),
                model_training = !missing(model_training), tuning = !missing(tuning),
                postprocessing = !missing(postprocessing))
  if (!is.null(schema)) {
    options <- .schema_options(schema, outcome, classification, event_level, explicit)
    classification <- options$classification
    event_level <- options$event_level
    if (!is.null(predictors) && !setequal(predictors, setdiff(schema$columns, outcome))) {
      stop("Selected predictors must agree with `schema`.", call. = FALSE)
    }
  }
  .validate_names(outcome, "outcome")
  if (length(outcome) != 1L) stop("`outcome` must name one column.", call. = FALSE)
  if (grepl("^(\\.pred($|_)|\\.row$|\\.config$|id([1-9][0-9]*)?$)", outcome)) {
    stop("Outcome ", sQuote(outcome), " is reserved for prediction/resample metadata; rename it explicitly.", call. = FALSE)
  }
  if (!is.null(predictors)) {
    .validate_names(predictors, "predictor")
    if (outcome %in% predictors) stop("Outcome ", sQuote(outcome), " is also a predictor.", call. = FALSE)
  }
  .validate_choice(method, c("inductive", "cross"), "method")
  .validate_choice(classification, c("binary", "ovo", "ovr"), "classification")
  .validate_choice(prediction_model, c("matched", "full"), "prediction_model")
  .validate_choice(model_training, c("shared", "pairwise"), "model_training")
  pairwise <- model_training == "pairwise"
  if (pairwise && classification != "ovo") {
    stop("`model_training = 'pairwise'` requires classification = 'ovo'.", call. = FALSE)
  }
  .integration_packages()
  if (!is.null(tuning)) .tuning_packages()
  postprocessing <- .validate_postprocessing(postprocessing)
  mode <- .validate_base_workflow(workflow, outcome, classification, allow_tune = !is.null(tuning),
                                  postprocessing = postprocessing)
  if (!is.null(postprocessing)) workflow <- .postprocessing_workflow(workflow, postprocessing)
  if (!is.null(schema)) {
    if (!identical(mode, schema$mode)) stop("Base workflow mode must agree with `schema`.", call. = FALSE)
    pre <- workflows::extract_preprocessor(workflow)
    if (inherits(pre, "recipe") && !is.null(pre$template)) .schema_encoded_data(pre$template, schema, truth = TRUE)
  }
  regression <- mode == "regression"
  if ((regression || classification == "binary") && explicit[["model_training"]]) {
    stop("`model_training` requires classification = 'ovo' or 'ovr'.", call. = FALSE)
  }
  .validate_calibration_options(regression, classification, event_level, precision, loss, m, epsilon, explicit)
  sized <- .validate_split_options(method, calibration_prop, training_size, calibration_size, v, strata,
                                   splits, seed, pairwise, explicit)
  tuning <- .validate_tuning(tuning, workflow, mode, model_training,
                             loss, precision, seed, m, epsilon, postprocessing)
  .validate_tuning_compatibility(tuning, splits, postprocessing, sized, training_size, calibration_size, pairwise)
  .register_vennabers()
  requested <- list(mode = mode, predictors = predictors, method = method, classification = classification,
                    event_level = event_level, calibration_prop = calibration_prop, v = v,
                    strata = strata, splits = splits, seed = seed, precision = precision, loss = loss,
                    m = m, epsilon = epsilon, prediction_model = prediction_model,
                    model_training = if (!regression && classification != "binary") model_training else NULL,
                    tuning = tuning, postprocessing = postprocessing)
  # Size fields exist only when requested, so default specifications keep their stored shape.
  if (sized) requested[c("training_size", "calibration_size")] <- list(training_size, calibration_size)
  if (sized) explicit[c("training_size", "calibration_size")] <- c(!is.null(training_size), !is.null(calibration_size))
  spec <- parsnip::new_model_spec("vennabers", args = list(), eng_args = NULL,
                                mode = mode, method = NULL, engine = method)
  spec <- parsnip::set_engine(spec, method, base_workflow = workflow,
                             outcome = outcome, parameters = list(requested = requested, explicit = explicit, schema = schema))
  selection <- if (!is.null(predictors)) predictors else if (!is.null(schema)) setdiff(unname(schema$columns), outcome) else quote(everything())
  carrier <- rlang::inject(workflows::add_variables(
    workflows::workflow(), outcomes = .integration_select_outcome(!!outcome, !!classification, !!mode, !!schema), predictors = !!selection,
    blueprint = hardhat::default_xy_blueprint(intercept = FALSE, allow_novel_levels = TRUE)))
  carrier <- workflows::add_model(carrier, spec)
  # Check schema bindings before workflows::predict() discards input attributes.
  if (!is.null(schema)) class(carrier) <- c("vennabers_workflow", class(carrier))
  carrier
}

# Arguments before ... also need exact names (R otherwise partially matches).
.exact_argument_names <- function() {
  supplied <- names(sys.call(-1L))[-1L]
  bad <- supplied[nzchar(supplied) & !supplied %in% names(formals(sys.function(-1L)))]
  if (length(bad)) stop("Unknown argument: ", sQuote(bad[1L]), ".", call. = FALSE)
}

# Task options that do not fit the base workflow mode or classification. Shared
# by calibrated_workflow() and calibrate_workflow(); `explicit` flags the
# options the caller supplied.
.validate_task_options <- function(regression, classification, explicit) {
  if (regression && any(explicit[c("classification", "event_level", "precision", "loss")])) {
    stop("Regression does not accept explicit classification, event_level, precision or loss options.", call. = FALSE)
  }
  if (!regression && any(explicit[c("m", "epsilon")])) {
    stop("`m` and `epsilon` require a regression workflow.", call. = FALSE)
  }
  if (!regression && classification != "binary" && explicit[["event_level"]]) {
    stop("`event_level` requires classification = 'binary'.", call. = FALSE)
  }
}

# Single-value calibration options; tune() placeholders are checked by .validate_tuning().
.validate_calibration_options <- function(regression, classification, event_level, precision, loss,
                                          m, epsilon, explicit) {
  .validate_task_options(regression, classification, explicit)
  .validate_moderation(m, epsilon, tunable = TRUE)
  .validate_choice(event_level, c("first", "second"), "event_level")
  if (is.null(.tune_id(loss))) .validate_choice(loss, c("log", "brier"), "loss")
  if (is.null(.tune_id(precision))) .validate_precision(precision)
}

# The split plan options; returns whether inductive sizes were requested.
.validate_split_options <- function(method, calibration_prop, training_size, calibration_size, v, strata,
                                    splits, seed, pairwise, explicit) {
  cross <- method == "cross"
  if (!cross && explicit[["v"]]) stop("`v` requires method = 'cross'.", call. = FALSE)
  if (cross && explicit[["calibration_prop"]]) {
    stop("`calibration_prop` requires method = 'inductive'.", call. = FALSE)
  }
  if (!.is_count(v, 2)) {
    stop("`v` must be one finite whole number of at least two folds.", call. = FALSE)
  }
  if (!.is_proportion(calibration_prop)) {
    stop("`calibration_prop` must be one number strictly between 0 and 1.", call. = FALSE)
  }
  sized <- .validate_inductive_sizes(training_size, calibration_size)
  if (sized) {
    if (cross) stop("`training_size` and `calibration_size` require method = 'inductive'.", call. = FALSE)
    if (explicit[["calibration_prop"]]) {
      stop("Supply either `calibration_prop` or `training_size`/`calibration_size`, not both.", call. = FALSE)
    }
    if (!is.null(splits)) {
      stop("Explicit `splits` partition all rows and cannot be combined with `training_size` or `calibration_size`.", call. = FALSE)
    }
  }
  if (!is.null(strata) && !identical(strata, FALSE)) {
    stop("`strata` must be NULL (classification outcome stratification; none for regression) or FALSE.", call. = FALSE)
  }
  if (!is.null(seed) && !.is_count(seed, -.Machine$integer.max)) {
    stop("`seed` must be NULL or one integer in R's seed range.", call. = FALSE)
  }
  if (!is.null(splits)) {
    if (pairwise && !is.function(splits) &&
        (!is.list(splits) || inherits(splits, c("rsplit", "rset", "data.frame")))) {
      stop("Pairwise `splits` must be a function of local pair data or a named list of pair plans.", call. = FALSE)
    }
    if (!pairwise && !inherits(splits, if (cross) "rset" else "rsplit") && !is.function(splits)) {
      stop("`splits` must be NULL, an IID ", if (cross) "Cross rset" else "holdout rsplit",
           ", or a function of local data.", call. = FALSE)
    }
    if (any(explicit[c("calibration_prop", "v", "strata")])) {
      stop("Explicit `splits` cannot be combined with `calibration_prop`, `v` or `strata`.", call. = FALSE)
    }
  }
  sized
}

# Calibration and joint model candidates are calibrated workflows fitted again
# on resample data, so their split plans must be providers or proportions.
.validate_tuning_compatibility <- function(tuning, splits, postprocessing, sized,
                                           training_size, calibration_size, pairwise) {
  resampled <- !is.null(tuning$calibration) || identical(tuning$strategy, "joint") && !is.null(tuning$model)
  if (pairwise && identical(tuning$strategy, "joint") && !is.null(tuning$model) &&
      !is.null(splits) && !is.function(splits)) {
    stop("Pairwise joint model selection requires local pair plan providers, not fixed plans.", call. = FALSE)
  }
  if (!is.null(tuning$calibration) && !is.null(splits) && !is.function(splits)) {
    stop("Calibration grid selection requires a local `splits` provider function or NULL.", call. = FALSE)
  }
  # A static B/P plan cannot bind to resample data.
  if (resampled && !is.null(postprocessing$splits) && !is.function(postprocessing$splits)) {
    stop("Joint or calibration selection with postprocessing requires a `postprocessing$splits` provider function or prop/v.", call. = FALSE)
  }
  # Only proportions carry over when the sizes are resolved again on resample data.
  if (resampled && sized && any(c(training_size, calibration_size) >= 1)) {
    stop("Calibration tuning and joint model tuning fit Venn-Abers candidates on resamples; ",
         "supply `training_size` and `calibration_size` as proportions, not row counts.", call. = FALSE)
  }
}

#' @export
predict.vennabers_workflow <- function(object, new_data, ...) {
  .schema_encoded_data(new_data, workflows::extract_fit_engine(object)$schema)
  NextMethod()
}

.integration_select_outcome <- function(outcome, classification, mode, schema = NULL) {
  .validate_names(tidyselect::peek_vars(), "input column")
  data <- tidyselect::peek_data()
  if (!is.null(schema)) .schema_encoded_data(data, schema, truth = TRUE)
  if (outcome %in% names(data)) .validate_truth(data[[outcome]], outcome, classification, mode)
  outcome
}

.validate_base_workflow <- function(workflow, outcome, classification, allow_tune = FALSE, postprocessing = NULL) {
  if (!inherits(workflow, "workflow")) stop("`workflow` must be an unfitted workflow.", call. = FALSE)
  fitted <- tryCatch(workflows::extract_fit_parsnip(workflow), error = function(e) NULL)
  if (!is.null(fitted)) stop("`workflow` is already fitted.", call. = FALSE)
  post <- tryCatch(workflows::extract_postprocessor(workflow), error = function(e) NULL)
  if (!is.null(post)) {
    if (is.null(postprocessing)) stop("Additional workflow postprocessors require an explicit matching `postprocessing` request.", call. = FALSE)
    .postprocessing_workflow(workflow, postprocessing)
  }
  spec <- workflows::extract_spec_parsnip(workflow)
  if (!spec$mode %in% c("classification", "regression")) {
    stop("The base workflow model mode must be classification or regression.", call. = FALSE)
  }
  if (is.null(spec$engine)) stop("Configure the base workflow model engine before calibration.", call. = FALSE)
  pre <- workflows::extract_preprocessor(workflow)
  if (inherits(pre, "recipe")) {
    if (recipes::fully_trained(pre) || any(vapply(pre$steps, function(x) isTRUE(recipes::is_trained(x)), logical(1)))) {
      stop("The base workflow recipe is already prepped/trained.", call. = FALSE)
    }
    targets <- pre$var_info$variable[pre$var_info$role == "outcome"]
    if (!identical(targets, outcome)) stop("Base recipe outcome must be ", sQuote(outcome), ".", call. = FALSE)
    if (outcome %in% names(pre$template)) .validate_truth(pre$template[[outcome]], outcome, classification, spec$mode)
  } else if (inherits(pre, "formula")) {
    if (length(pre) != 3L || !identical(pre[[2L]], as.name(outcome))) {
      stop("Base formula outcome must be the raw column ", sQuote(outcome), ".", call. = FALSE)
    }
  }
  if (!allow_tune && nrow(generics::tune_args(workflow))) stop("Resolve all tune() placeholders before calibration.", call. = FALSE)
  spec$mode
}
