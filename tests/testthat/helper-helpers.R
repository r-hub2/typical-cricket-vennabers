register_reference_engine <- function(model, mode, fit, predict, type, engine = "reference", new_model = TRUE) {
  assign(fit, get(fit, envir = parent.frame()), envir = globalenv())
  assign(predict, get(predict, envir = parent.frame()), envir = globalenv())
  if (new_model) parsnip::set_new_model(model)
  parsnip::set_model_mode(model, mode)
  parsnip::set_model_engine(model, mode, engine)
  parsnip::set_fit(model, mode, engine,
    list(interface = "data.frame", protect = c("x", "y"), func = c(fun = fit), defaults = list()))
  parsnip::set_encoding(model, mode, engine,
    list(predictor_indicators = "none", compute_intercept = FALSE, remove_intercept = FALSE, allow_sparse_x = FALSE))
  parsnip::set_pred(model, mode, engine, type,
    list(pre = NULL, post = NULL, func = c(fun = predict),
      args = list(object = quote(object$fit), new_data = quote(new_data))))
}

pairwise_reference <- function(model = "va_pair_reference") {
  audit <- new.env()
  audit$fits <- audit$predictions <- list()
  pair_fit <- function(x, y, fail_size = integer(), bad = FALSE) {
    audit$fits[[length(audit$fits) + 1L]] <- list(x = x, y = y)
    if (nrow(x) %in% fail_size) stop("planned pair failure")
    list(x = x, y = y, levels = levels(y), bad = bad, draw = runif(1))
  }
  pair_prob <- function(object, new_data) {
    stopifnot(nrow(new_data) > 0L)
    audit$predictions[[length(audit$predictions) + 1L]] <- new_data$row_id
    p <- vapply(seq_along(object$levels), function(i) {
      exp(sin(new_data$row_id / (i + 1)) + mean(object$x$x[object$y == object$levels[i]]))
    }, numeric(nrow(new_data)))
    p <- matrix(p, nrow(new_data), length(object$levels))
    p <- as.data.frame(p / rowSums(p))
    names(p) <- paste0(".pred_", object$levels)
    if (object$bad) p[[1L]] <- NA_real_
    p[rev(seq_along(p))]
  }
  register_reference_engine(model, "classification", "pair_fit", "pair_prob", "prob")
  spec <- parsnip::new_model_spec(model, args = list(), eng_args = NULL,
    mode = "classification", method = NULL, engine = "reference")
  make_base <- function(data) workflows::workflow() |> workflows::add_recipe(
    recipes::recipe(truth ~ ., data) |> recipes::step_impute_median(x) |> recipes::step_normalize(x)) |>
    workflows::add_model(spec)
  pair_data <- function(data) setNames(lapply(combn(levels(data$truth), 2, simplify = FALSE), function(pair) {
    local <- data[data$truth %in% pair, , drop = FALSE]
    local$truth <- factor(local$truth, levels = pair)
    local
  }), paste0("pair_", seq_len(choose(nlevels(data$truth), 2))))
  couple <- function(q, classes) {
    pairs <- combn(classes, 2, simplify = FALSE)
    p <- vapply(classes, function(cl) {
      relevant <- which(vapply(pairs, function(pair) cl %in% pair, logical(1)))
      odds <- lapply(relevant, function(j) {
        prob <- if (cl == pairs[[j]][2L]) q[, j] else 1 - q[, j]
        (1 - prob) / prob
      })
      1 / (1 + Reduce(`+`, odds))
    }, numeric(nrow(q)))
    p <- matrix(p, nrow(q), length(classes))
    p / rowSums(p)
  }
  list(audit = audit, pair_fit = pair_fit, pair_prob = pair_prob, spec = spec,
       make_base = make_base, pair_data = pair_data, couple = couple)
}

worker_replay_available <- function() {
  all(vapply(c("tune", "yardstick", "future", "withr"), requireNamespace, logical(1), quietly = TRUE))
}

# Fresh-process script lines for a script holding `s$replay` and `s$snapshot`:
# two multisession workers, each a process other than the script's own, must
# rebuild the saved outer evaluation (or a sequential run of it) and every state.
worker_replay <- function(metrics, packages, untraced) c(
  sprintf('untrace("%s", where = asNamespace("%s"))', untraced, names(untraced)),
  'Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")',
  'r <- s$replay',
  paste0('control <- function(extract) tune::control_resamples(save_pred = TRUE, extract = extract, parallel_over = "resamples", pkgs = ',
    paste(deparse(packages), collapse = ""), ')'),
  paste0('run <- function(extract = identity) withr::with_seed(43, tune::fit_resamples(r$model, r$outer, metrics = ', metrics,
    ', control = control(extract)))'),
  'reference <- if (is.null(r$evaluation)) run() else r$evaluation',
  'future::plan(future::multisession, workers = 2)',
  'parallel <- run(function(x) list(workflow = x, pid = Sys.getpid()))',
  'future::plan(future::sequential)',
  'found <- tune::collect_extracts(parallel)$.extracts; expected <- tune::collect_extracts(reference)$.extracts',
  'pids <- vapply(found, function(x) x$pid, 0L)',
  'stopifnot(nrow(tune::collect_notes(parallel)) == 0L, length(found) == nrow(r$outer), length(expected) == length(found),',
  '  !any(pids == Sys.getpid()), length(unique(pids)) == 2L,',
  '  identical(tune::collect_predictions(parallel), tune::collect_predictions(reference)),',
  '  identical(tune::collect_metrics(parallel), tune::collect_metrics(reference)))',
  'for (i in seq_along(found)) stopifnot(identical(s$snapshot(found[[i]]$workflow, r$new, r$schema), s$snapshot(expected[[i]], r$new, r$schema)))')

# Skips unless every package calibrated_workflow() requires, plus each additional
# optional package named in `...` (for example "tune" or "nnet"), is installed. Minimum versions come from the tables the
# package itself checks (R/registration.R), so a partly available or outdated
# integration skips instead of failing; a complete environment never skips here.
skip_if_no_integration <- function(...) {
  integration <- vennabers:::.integration_versions
  minimums <- c(integration, vennabers:::.tuning_versions, vennabers:::.tailor_versions)
  for (package in union(names(integration), c(...))) {
    testthat::skip_if_not_installed(package, minimum_version = if (package %in% names(minimums)) minimums[[package]])
  }
  invisible(NULL)
}

# Fresh-process script lines that load the installed vennabers under test (the
# copy at `root`) with the caller's library paths, failing instead of using
# any other copy.
installed_vennabers_lines <- function(root = getNamespaceInfo(asNamespace("vennabers"), "path")) c(
  paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
  "library(vennabers)",
  paste0("stopifnot(normalizePath(find.package('vennabers')) == normalizePath(", deparse(root), "))"))
