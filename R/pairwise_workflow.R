# Pairwise model training: every class pair is a component with its own local
# rows, recipe-retargeted workflow and binary calibrator; shared phases do the
# fitting and prediction, so this file only retargets recipes, builds pair
# domains and couples pair points into class probabilities.
.fit_pairwise <- function(base_workflow, data, state, tuning_plan) {
  settings <- state$parameters$requested
  pairs <- utils::combn(state$class_levels, 2L, simplify = FALSE)
  ids <- paste0("pair_", seq_along(pairs))
  pair_plans <- settings$splits
  named <- !is.null(pair_plans) && !is.function(pair_plans)
  if (named && (is.null(names(pair_plans)) || anyNA(names(pair_plans)) || anyDuplicated(names(pair_plans)) ||
      !setequal(names(pair_plans), ids))) {
    stop("Pair plans must be named uniquely by exactly the required pair IDs: ",
         paste(ids, collapse = ", "), ".", call. = FALSE)
  }
  # Resolve and validate every pair before any recipe or model learns.
  planned <- lapply(seq_along(pairs), function(i) {
    id <- ids[i]
    pair <- pairs[[i]]
    rows <- which(data[[state$outcome]] %in% pair)
    local <- data[rows, , drop = FALSE]
    local[[state$outcome]] <- factor(local[[state$outcome]], levels = pair,
                                     ordered = is.ordered(data[[state$outcome]]))
    workflow <- base_workflow
    recipe <- workflows::extract_preprocessor(base_workflow)
    if (inherits(recipe, "recipe")) {
      # recipes has no public ptype setter. Change only the outcome domain;
      # retain all steps, roles and blueprint options, without pre-learning them.
      recipe$template <- recipe$template[0L, , drop = FALSE]
      recipe$template[[state$outcome]] <- local[[state$outcome]][0L]
      if (!is.null(recipe$ptype)) recipe$ptype[[state$outcome]] <- local[[state$outcome]][0L]
      workflow <- workflows::update_recipe(workflow, recipe,
        blueprint = base_workflow$pre$actions$recipe$blueprint)
    }
    list(id = id, pair = pair, rows = rows, local = local, workflow = workflow,
      plan = if (is.function(pair_plans)) pair_plans else if (named) pair_plans[[id]] else NULL)
  })
  domains <- lapply(planned, function(entry) {
    pair <- entry$pair
    list(data = entry$local, workflow = entry$workflow, plan = entry$plan, pair_id = entry$id,
      pair_rows = entry$rows, schema = state$schema,
      outcome = state$outcome, strata = settings$strata, classes = pair, event = pair[2L],
      y = entry$local[[state$outcome]], rows_word = "local pair rows",
      fields = list(pair_id = entry$id, class_levels = pair, event_class = pair[2L], rows = entry$rows))
  })
  task <- .binary_task()
  # Finish every source fit first so full-model RNG draws cannot change calibration.
  result <- if (is.null(tuning_plan$pair_cache)) .fit_domains(domains, task, settings, tuning_plan) else
    .fit_pairwise_grid_domains(data, domains, task, settings, tuning_plan)
  components <- stats::setNames(lapply(result$fitted, function(component) component), ids)
  state$pair_models <- components
  # Sizes resolve on each pair's rows, so the resolved counts are named by pair ID.
  if (.inductive_sized(settings)) {
    roles <- c(training_size = "analysis", calibration_size = "assessment", unassigned_size = "unassigned")
    state$parameters$effective[names(roles)] <- lapply(roles, function(role) {
      stats::setNames(vapply(result$plans, function(plan) length(plan$Inductive1[[role]]), 1L), ids)
    })
  }
  if (!is.null(state$tuning)) {
    state$tuning$models <- lapply(components, function(component) {
      selections <- component$model_selection
      for (fold in names(selections)) if (!is.null(selections[[fold]])) {
        if (is.null(settings$postprocessing)) {
          selections[[fold]]$pair_id <- component$pair_id
          selections[[fold]]$fold_id <- fold
        } else for (inner in names(selections[[fold]])) {
          # Inner B searches already carry their fold ID.
          selections[[fold]][[inner]]$pair_id <- component$pair_id
        }
      }
      selections
    })
    state$tuning$sources <- lapply(components, `[[`, "sources")
    state$tuning$fit_roles <- lapply(components, `[[`, "folds")
  }
  if (!is.null(state$tuning)) state$tuning$cost <- list(observed_base_fits = tuning_plan$counter$n)
  state
}

# An inductive pair fold must train and calibrate on both of its classes.
.pair_fold_classes <- function(folds, y) {
  for (part in c("analysis", "assessment")) for (fold_id in names(folds)) {
    if (any(tabulate(y[folds[[fold_id]][[part]]], nbins = 2L) == 0L)) {
      stop("Both pair classes must occur in ", part, " (training/calibration).", call. = FALSE)
    }
  }
}

# Check the whole pair/class mapping before any pair component predicts.
.validate_pair_components <- function(object) {
  pairs <- utils::combn(object$class_levels, 2L, simplify = FALSE)
  ids <- paste0("pair_", seq_along(pairs))
  if (!identical(names(object$pair_models), ids)) {
    stop("Pair/class mapping does not match the fitted global class order.", call. = FALSE)
  }
  cross <- object$parameters$requested$method == "cross"
  for (i in seq_along(pairs)) {
    id <- ids[i]
    pair <- pairs[[i]]
    component <- object$pair_models[[id]]
    fold_ids <- if (cross) component$calibrator$folds else "Inductive1"
    if (!identical(component$pair_id, id) || !identical(component$class_levels, pair) ||
        !identical(component$event_class, pair[2L]) ||
        !identical(names(component$sources), fold_ids) ||
        !identical(names(component$models), unique(unname(component$sources))) ||
        !identical(names(component$folds), fold_ids) ||
        !identical(names(component$calibration), fold_ids) ||
        (cross && (length(fold_ids) < 2L || anyNA(fold_ids) || any(!nzchar(fold_ids)) ||
                   anyDuplicated(fold_ids) || !identical(names(component$calibrator$calibrators), fold_ids)))) {
      stop(id, " model/pair/fold identity does not match the fitted mapping.", call. = FALSE)
    }
  }
  invisible(NULL)
}

# Pair detail rows name their pair and classes; a matched pair's model source
# is the pair itself, a full one its pair-qualified shared refit.
.pair_details <- function(part, details, prediction_model) {
  k <- nrow(details)
  pair <- part$class_levels
  data.frame(row = details$row, pair_id = rep(part$pair_id, k),
    fold = details$fold, method = rep("ovo", k), event_class = details$event_class,
    class_1 = rep(pair[1L], k), class_2 = rep(pair[2L], k),
    details[c("lower", "upper", "probability", "model_training", "prediction_model")],
    model_source = rep(if (identical(prediction_model, "full")) paste(part$pair_id, part$sources[[1L]], sep = "/") else
      part$pair_id, k))
}

# Couple the pair components' event probabilities into class probabilities.
.couple_pairs <- function(object, parts) {
  pairs <- lapply(object$pair_models, `[[`, "class_levels")
  points <- do.call(cbind, lapply(seq_along(parts), function(i) {
    parts[[i]]$prediction[[paste0(".pred_", pairs[[i]][2L])]]
  }))
  prediction <- .couple_ovo(points, unname(pairs), object$class_levels)
  colnames(prediction) <- paste0(".pred_", object$class_levels)
  result <- list(prediction = tibble::as_tibble(prediction, .name_repair = "minimal"))
  details <- lapply(parts, `[[`, "details")
  if (!is.null(details[[1L]])) result$details <- do.call(rbind, unname(details))
  result
}

# Check wrappers as well as every leaf workflow: swapping one inner member
# must not silently mix a different class pair into the global prediction.
.validate_pair_model <- function(model, pair, context) {
  if (inherits(model, "vennabers_ensemble")) {
    if (!identical(names(model$models), model$inner_folds) || length(model$models) < 2L ||
        any(vapply(model$models, is.null, logical(1)))) {
      stop(context, " incomplete inner ensemble.", call. = FALSE)
    }
    for (j in names(model$models)) .validate_pair_model(model$models[[j]], pair, paste(context, "/", j))
    classes <- model$classes
  } else if (inherits(model, c("vennabers_sigmoid", "vennabers_tailor_workflow"))) {
    .validate_pair_model(model$workflow, pair, context)
    classes <- model$classes
  } else {
    classes <- if (inherits(model, "workflow")) levels(workflows::extract_mold(model)$outcomes[[1L]])
  }
  if (!identical(classes, pair)) {
    stop(context, " model/pair identity does not match the fitted mapping.", call. = FALSE)
  }
  invisible(NULL)
}
