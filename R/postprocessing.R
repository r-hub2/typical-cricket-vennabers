.validate_postprocessing <- function(post) {
  if (is.null(post)) return(NULL)
  .tuning_names(post, c("method", "strategy", "prop", "v", "splits", "sigmoid_input"), "postprocessing")
  .validate_choice(post$method, c("tailor_logistic", "tailor_multinomial", "tailor_linear", "sigmoid", "isotonic_numeric"), "postprocessing$method")
  .validate_choice(post$strategy, c("holdout", "ensemble"), "postprocessing$strategy")
  if (post$method == "sigmoid") {
    .validate_choice(post$sigmoid_input, c("probability", "glm_link"), "postprocessing$sigmoid_input")
  } else if ("sigmoid_input" %in% names(post)) {
    stop("`postprocessing$sigmoid_input` requires sigmoid.", call. = FALSE)
  }
  if ("splits" %in% names(post)) {
    if (any(c("prop", "v") %in% names(post))) stop("Postprocessing splits cannot be combined with prop/v.", call. = FALSE)
    if (!is.function(post$splits) && !inherits(post$splits, if (post$strategy == "holdout") "rsplit" else "rset")) {
      stop("`postprocessing$splits` must supply a local plan or function(data).", call. = FALSE)
    }
  } else if (post$strategy == "holdout") {
    if ("v" %in% names(post)) stop("Postprocessing holdout does not accept v.", call. = FALSE)
    if (!.is_proportion(post$prop)) {
      stop("`postprocessing$prop` must be one finite numeric P proportion strictly between 0 and 1.", call. = FALSE)
    }
  } else {
    if ("prop" %in% names(post)) stop("Postprocessing ensemble does not accept prop.", call. = FALSE)
    if (!.is_count(post$v, 2)) {
      stop("`postprocessing$v` must be one whole number of at least two.", call. = FALSE)
    }
  }
  .postprocessing_packages(post)
  post
}

.postprocessing_workflow <- function(workflow, postprocessing) {
  if (postprocessing$method == "isotonic_numeric") {
    if (!is.null(tryCatch(workflows::extract_postprocessor(workflow), error = function(e) NULL))) {
      stop("isotonic_numeric requires a base workflow without a postprocessor.", call. = FALSE)
    }
    return(workflow)
  }
  if (postprocessing$method == "sigmoid") {
    .validate_sigmoid_workflow(workflow, postprocessing$sigmoid_input)
    return(workflow)
  }
  .configure_tailor_workflow(workflow, sub("^tailor_", "", postprocessing$method))
}

.plan_postprocessing <- function(folds, domain, settings) {
  post <- settings$postprocessing
  numeric <- is.null(domain$classes)
  k <- if (numeric) 2L else length(domain$classes)
  expected <- if (numeric) "tailor_linear" else if (k == 2L) "tailor_logistic" else "tailor_multinomial"
  # Sigmoid needs a binary component (a shared binary fit or a class pair): one with an event class.
  if (post$method != expected && !(post$method == "isotonic_numeric" && numeric) &&
      !(post$method == "sigmoid" && !is.null(domain$event))) {
    stop(paste(c(domain$pair_id, post$method), collapse = " "), " B/P plan: ", if (numeric) "Regression requires" else paste(k, "model classes require"),
      " explicitly selecting ", if (numeric) "tailor_linear or isotonic_numeric" else expected, ".", call. = FALSE)
  }
  ensemble <- post$strategy == "ensemble"
  for (id in names(folds)) {
    inner_id <- if (ensemble) "ensemble" else "Holdout1"
    folds[[id]]$postprocessing <- tryCatch({
      Q <- folds[[id]]$analysis
      local <- domain$data[Q, , drop = FALSE]
      plan <- post$splits
      if (is.function(plan)) plan <- plan(local)
      if (!"splits" %in% names(post)) {
        outcome <- if (numeric) NULL else domain$outcome
        if (ensemble && post$v > nrow(local)) stop("v cannot exceed local Q rows; folds are not reduced.", call. = FALSE)
        plan <- if (ensemble) rlang::inject(rsample::vfold_cv(local, v = post$v, strata = !!outcome, pool = 0)) else
          rlang::inject(rsample::initial_split(local, prop = 1 - post$prop, strata = !!outcome, pool = 0))
      }
      bindings <- if (inherits(plan, "rset")) plan$splits else list(plan)
      for (split in bindings) if (inherits(split, "rsplit") && !identical(names(split$data), names(local))) {
        stop("Split data binding columns and their order must match local Q.", call. = FALSE)
      }
      positions <- if (ensemble) .cross_splits(plan, local) else list(Holdout1 = .inductive_split(plan, local))
      pairs <- positions
      for (inner_id in names(positions)) {
        pair <- positions[[inner_id]]
        if (!numeric) for (role in c("analysis", "assessment")) {
          rows <- pair[[role]]
          label <- if (role == "analysis") "B" else "P"
          if (any(tabulate(local[[domain$outcome]][rows], nbins = length(domain$classes)) == 0L)) {
            stop(label, " (", length(rows), " rows) must contain all model classes: ", paste(domain$classes, collapse = ", "), ".", call. = FALSE)
          }
        }
        minimum <- if (post$method %in% c("sigmoid", "isotonic_numeric")) 1L else k
        if (length(pair$assessment) <= minimum) {
          stop("P has ", length(pair$assessment), " rows; requires at least ", minimum + 1L,
            if (!post$method %in% c("sigmoid", "isotonic_numeric")) paste0(" (more than design rank ", k, ")"), ".", call. = FALSE)
        }
        pairs[[inner_id]] <- list(Q = Q, B = Q[pair$analysis], P = Q[pair$assessment],
          analysis = pair$analysis, assessment = pair$assessment)
      }
      # A holdout plan stays the bare B/P pair, as stored since v0.9.0.
      if (ensemble) pairs else pairs[[1L]]
    }, error = function(e) stop(.domain_label(domain, id), " / ", inner_id, " ", post$method, " B/P plan: ", conditionMessage(e), call. = FALSE))
  }
  folds
}

.fit_postprocessing <- function(workflow, data, fold, domain, settings, id, tuning_plan = NULL) {
  plans <- if (settings$postprocessing$strategy == "holdout") list(Holdout1 = fold$postprocessing) else
    fold$postprocessing
  pairs <- lapply(names(plans), function(j) {
    .fit_postprocessing_pair(workflow, data, plans[[j]], fold, domain, settings, id, j,
      tuning_plan = tuning_plan)
  })
  names(pairs) <- names(plans)
  if (settings$postprocessing$strategy == "holdout") return(pairs[[1L]])
  diagnostics <- lapply(pairs, `[[`, "diagnostics")
  model <- structure(list(models = lapply(pairs, `[[`, "model"), classes = domain$classes,
    inner_folds = names(pairs), method = settings$postprocessing$method, context = paste(.domain_label(domain, id), settings$postprocessing$method)), class = "vennabers_ensemble")
  selection <- do.call(c, unname(lapply(pairs, `[[`, "selection")))
  identity <- if (!is.null(domain$pair_id)) list(pair_id = domain$pair_id, pair_rows = domain$pair_rows) else list()
  list(model = model, selection = selection, diagnostics = c(identity, list(format_version = 1L,
    method = settings$postprocessing$method, strategy = "ensemble", input = diagnostics[[1L]]$input,
    Q = fold$analysis, C = fold$assessment, VA_fold = id, inner_folds = diagnostics,
    observed_base_fits = length(pairs), observed_postprocessor_fits = length(pairs))))
}

.fit_postprocessing_pair <- function(workflow, data, plan, fold, domain, settings, id, inner_id,
                                     tuning_plan = NULL) {
  identity <- if (!is.null(domain$pair_id)) list(pair_id = domain$pair_id, pair_rows = domain$pair_rows) else list()
  method <- sub("^tailor_", "", settings$postprocessing$method)
  context <- paste0(.domain_label(domain, id), " / ", inner_id, " ", settings$postprocessing$method, " (B=", length(plan$B), ", P=", length(plan$P), ")")
  warnings <- character()
  selection <- NULL
  result <- tryCatch(withCallingHandlers({
    # Local search sees only B; the engine/recipe choice is fixed before P is used.
    # One draw per B gives each search its own stream: data-dependent draws of
    # another B's search (which may contain this P) cannot shift this one.
    search_seed <- if (!is.null(tuning_plan$model)) sample.int(.Machine$integer.max, 1L)
    choose <- function() .choose_model(workflow, data[plan$B, , drop = FALSE], domain, settings, tuning_plan, id)
    selected <- if (is.null(search_seed)) choose() else withr::with_seed(search_seed, choose())
    workflow <- selected$workflow
    if (!is.null(selected$diagnostics)) {
      selection <- stats::setNames(list(c(selected$diagnostics,
        list(fold_id = id, inner_fold = inner_id, B = plan$B, search_seed = search_seed))), inner_id)
    }
    if (!is.null(tuning_plan$counter)) tuning_plan$counter$n <- tuning_plan$counter$n + 1L
    if (method == "isotonic_numeric") {
      base <- generics::fit(workflow, data[plan$B, , drop = FALSE])
      prediction <- .base_numeric(base, data[plan$P, , drop = FALSE], context)
      truth <- data[[domain$outcome]][plan$P]
      iso <- stats::isoreg(prediction, truth)
      sorted_rows <- order(prediction)
      distinct_support <- !duplicated(prediction[sorted_rows])
      knots <- list(x = prediction[sorted_rows][distinct_support], y = iso$yf[distinct_support])
      model <- structure(list(workflow = base, isotonic = knots), class = "vennabers_isotonic_workflow")
      result <- list(model = model)
    } else if (method == "sigmoid") {
      result <- .fit_sigmoid_workflow(workflow, data, plan, domain, settings$postprocessing$sigmoid_input, context)
      model <- result$model$workflow
    } else {
      result <- .fit_tailor_workflow(workflow, data[plan$B, , drop = FALSE],
        data[plan$P, , drop = FALSE], domain, method, context)
      model <- result$model
    }
    engine <- workflows::extract_fit_engine(if (method == "isotonic_numeric") model$workflow else model)
    if (is.list(engine) && (identical(engine$converged, FALSE) ||
        inherits(engine, "nnet") && !isTRUE(engine$convergence == 0))) {
      stop("B engine did not converge.", call. = FALSE)
    }
    result
  }, warning = function(w) {
    message <- paste0(context, ": ", conditionMessage(w))
    warnings <<- c(warnings, message)
    warning(message, call. = FALSE)
    invokeRestart("muffleWarning")
  }), error = function(e) stop(context, ": ", conditionMessage(e), call. = FALSE))
  list(model = result$model, selection = selection, diagnostics = c(plan, identity, list(format_version = 1L,
    method = settings$postprocessing$method, strategy = "holdout",
    input = if (method == "sigmoid") settings$postprocessing$sigmoid_input else if (method %in% c("linear", "isotonic_numeric")) "numeric" else "probability", inner_fold = inner_id,
    VA_fold = id, C = fold$assessment, B_size = length(plan$B), P_size = length(plan$P),
    classes = domain$classes, event_class = domain$event, schema = domain$schema,
    columns = names(data), column_types = lapply(data, function(x) list(class = class(x), typeof = typeof(x), levels = levels(x))),
    seed = settings$seed, rng_kind = RNGkind(), rank = result$rank, coefficients = result$coefficients,
    scale = result$scale, sigmoid_class = if (method == "sigmoid") domain$classes[2L],
    converged = TRUE, warnings = warnings, selected = selection[[1L]]$selected,
    observed_base_fits = 1L, observed_postprocessor_fits = 1L)))
}

#' @export
predict.vennabers_isotonic_workflow <- function(object, newdata, type = "numeric", ...) {
  .validate_choice(type, "numeric", "type")
  if (!nrow(newdata)) return(tibble::tibble(.pred = numeric()))
  prediction <- .base_numeric(object$workflow, newdata, "isotonic_numeric")
  knots <- object$isotonic
  adjusted <- if (length(knots$x) == 1L) rep(knots$y, length(prediction)) else
    stats::approx(knots$x, knots$y, xout = prediction, rule = 2)$y
  tibble::tibble(.pred = adjusted)
}

#' @export
predict.vennabers_ensemble <- function(object, newdata, type = if (is.null(object$classes)) "numeric" else "prob", ...) {
  .validate_dots(...)
  regression <- is.null(object$classes)
  .validate_choice(type, if (regression) "numeric" else "prob", "type")
  if (!identical(names(object$models), object$inner_folds) || length(object$models) < 2L ||
      any(vapply(object$models, is.null, logical(1)))) {
    stop(object$context, " incomplete inner ensemble.", call. = FALSE)
  }
  if (!nrow(newdata) && regression) return(tibble::tibble(.pred = numeric()))
  if (!nrow(newdata)) return(tibble::as_tibble(stats::setNames(
    rep(list(numeric()), length(object$classes)), paste0(".pred_", object$classes))))
  # RDS restores S3 objects without necessarily loading their optional methods.
  if (object$method %in% c("tailor_logistic", "tailor_multinomial", "tailor_linear") && !requireNamespace("tailor", quietly = TRUE)) {
    stop(object$context, " prediction requires optional package 'tailor'.", call. = FALSE)
  }
  if (regression) {
    predictions <- lapply(object$inner_folds, function(j) {
      .base_numeric(object$models[[j]], newdata, paste(object$context, "/", j))
    })
    values <- do.call(cbind, predictions)
    return(tibble::tibble(.pred = vapply(seq_len(nrow(newdata)), function(i) .regression_mean(values[i, ]), 0)))
  }
  probabilities <- lapply(object$inner_folds, function(j) {
    .base_probability(object$models[[j]], newdata, object$classes, NULL, paste(object$context, "/", j))
  })
  averaged <- Reduce(`+`, probabilities) / length(probabilities)
  colnames(averaged) <- paste0(".pred_", object$classes)
  tibble::as_tibble(averaged)
}
