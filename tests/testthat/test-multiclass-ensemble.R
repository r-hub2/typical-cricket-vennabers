# Exercise the public workflow, stored diagnostics and manual VA seams.
ensemble_multiclass_data <- function(k = 3L) withr::with_seed(123, {
  truth <- factor(rep(c("z", "a", "b", "c")[seq_len(k)], length.out = 240),
    levels = c("z", "a", "b", "c")[seq_len(k)])
  data.frame(x = rnorm(240, sd = 2) + c(-1, 1, -1, 1)[truth],
    z = rnorm(240, sd = 2) + c(-1, -1, 1, 1)[truth],
    w = rnorm(240, sd = 2) + c(1, -1, -1, 1)[truth], truth = truth)
})
ensemble_multiclass_base <- function(d) workflows::workflow() |>
  workflows::add_recipe(recipes::recipe(truth ~ ., d) |>
    recipes::step_normalize(recipes::all_numeric_predictors())) |>
  workflows::add_model(if (nlevels(d$truth) == 2L) parsnip::logistic_reg(engine = "glm") else
    parsnip::multinom_reg(penalty = .1) |> parsnip::set_engine("nnet", trace = FALSE))
ensemble_multiclass_plan <- function(data) {
  rows <- order(ave(seq_len(nrow(data)), data$truth, FUN = seq_along))
  held <- list(rows[1:48], rows[49:108], rows[109:nrow(data)])
  rsample::manual_rset(lapply(held, function(p) rsample::make_splits(
    list(analysis = setdiff(seq_len(nrow(data)), p), assessment = p), data)), c("small", "medium", "large"))
}
# Training rows are calibration scores, where the last bit of a prediction can depend on
# the batch size (OpenBLAS) and move a Venn-Abers step; compare predictions on unseen rows.
ensemble_multiclass_unseen <- function(data) data[c("x", "z", "w")] + .01

testthat::test_that("shared ensembles average vectors before binary VA and multiclass merging", {
  for (pkg in c("workflows", "tailor", "probably", "rsample", "recipes", "nnet")) testthat::skip_if_not_installed(pkg)
  for (k in 2:4) {
    d <- ensemble_multiclass_data(k); base <- ensemble_multiclass_base(d)
    unseen <- ensemble_multiclass_unseen(d)
    classes <- levels(d$truth)
    method <- if (k == 2L) "logistic" else "multinomial"
    probability <- function(model, data) {
      p <- as.matrix(predict(model, data, type = "prob")[paste0(".pred_", classes)])
      colnames(p) <- classes
      p
    }
    for (va in c("inductive", "cross")) for (classification in c("ovo", "ovr")) for (precision in list(NULL, 0L, 2L)) {
      fit <- generics::fit(calibrated_workflow(base, "truth", method = va, classification = classification,
        precision = precision, seed = 24,
        postprocessing = list(method = paste0("tailor_", method), strategy = "ensemble", splits = ensemble_multiclass_plan)), d)
      state <- workflows::extract_fit_engine(fit)
      held <- new <- truth <- ids <- list()
      testthat::expect_identical(state$class_levels, classes)
      testthat::expect_null(state$event_class)
      for (id in names(state$models)) {
        ensemble <- state$models[[id]]; plans <- state$postprocessing[[id]]$inner_folds
        testthat::expect_identical(names(ensemble$models), c("small", "medium", "large"))
        testthat::expect_identical(names(plans), ensemble$inner_folds)
        testthat::expect_identical(state$postprocessing[[id]]$observed_base_fits, 3L)
        testthat::expect_identical(state$postprocessing[[id]]$observed_postprocessor_fits, 3L)
        calibration <- predictions <- list()
        for (j in names(plans)) {
          p <- plans[[j]]
          testthat::expect_identical(p$inner_fold, j)
          testthat::expect_identical(p$rank, k)
          testthat::expect_identical(p$classes, classes)
          testthat::expect_setequal(p$B, setdiff(p$Q, p$P))
          testthat::expect_length(intersect(p$B, p$P), 0L)
          testthat::expect_length(intersect(p$Q, p$C), 0L)
          ref <- generics::fit(workflows::add_tailor(base, tailor::tailor() |>
            tailor::adjust_probability_calibration(method = method, smooth = FALSE)),
            d[p$B, ], data_calibration = d[p$P, ])
          calibration[[j]] <- probability(ref, d[p$C, ])
          predictions[[j]] <- probability(ref, unseen)
          testthat::expect_lte(max(abs(predictions[[j]] - probability(ensemble$models[[j]], unseen))), 1e-10)
        }
        # Complete vectors, equal weights even though P sizes differ.
        held[[id]] <- Reduce(`+`, calibration) / 3
        new[[id]] <- Reduce(`+`, predictions) / 3
        sizes <- vapply(plans, `[[`, 1L, "P_size")
        weighted <- Reduce(`+`, Map(`*`, predictions, sizes / sum(sizes)))
        testthat::expect_gt(max(abs(new[[id]] - weighted)), 1e-6)
        testthat::expect_lte(max(abs(new[[id]] - probability(ensemble, unseen))), 1e-10)
        testthat::expect_lte(max(abs(held[[id]] - state$calibration[[id]]$probability)), 1e-10)
        truth[[id]] <- d$truth[plans[[1]]$C]; ids[[id]] <- rep(id, length(truth[[id]]))
      }
      core <- if (va == "cross") cross_venn_abers(do.call(c, truth), do.call(rbind, held), unlist(ids),
        method = classification, precision = precision) else
        venn_abers(truth[[1]], held[[1]], method = classification, precision = precision)
      for (loss in c("log", "brier")) {
        input <- if (va == "cross") new else new[[1]]
        expected <- predict(core, input, loss = loss)
        actual <- venn_abers_details(fit, unseen, loss = loss)
        testthat::expect_lte(max(abs(as.matrix(actual$prediction) - expected)), 1e-10)
        testthat::expect_identical(actual$prediction, predict(fit, unseen, type = "prob", loss = loss))
        testthat::expect_identical(predict(fit, unseen, type = "class", loss = loss)$.pred_class,
          predict(core, input, type = "class", loss = loss))
        testthat::expect_true(all(vapply(actual$details$inner_fold, identical, logical(1), c("small", "medium", "large"))))
        for (n in 0:1) {
          small <- venn_abers_details(fit, unseen[seq_len(n), ], loss = loss)
          testthat::expect_identical(nrow(small$prediction), n)
          testthat::expect_identical(names(small$details), names(actual$details))
          if (n) testthat::expect_lte(max(abs(as.matrix(small$prediction) - expected[1, , drop = FALSE])), 1e-10)
        }
      }
    }
  }
})

testthat::test_that("shared ensembles count complete fits and predictions and align every pair", {
  for (pkg in c("workflows", "tailor", "probably", "rsample", "recipes", "nnet")) testthat::skip_if_not_installed(pkg)
  original_fit <- getS3method("fit", "workflow", envir = asNamespace("generics"))
  original_predict <- getS3method("predict", "workflow", envir = asNamespace("stats"))
  original_tailor <- getS3method("predict", "tailor", envir = asNamespace("stats"))
  audit <- new.env(); audit$fits <- audit$predictions <- audit$tailor <- 0L; audit$broken <- FALSE
  fit_method <- function(object, ...) {
    if (!inherits(workflows::extract_spec_parsnip(object), "vennabers")) audit$fits <- audit$fits + 1L
    original_fit(object, ...)
  }
  predict_method <- function(object, ...) {
    p <- original_predict(object, ...)
    if (inherits(workflows::extract_spec_parsnip(object), "vennabers")) return(p)
    audit$predictions <- audit$predictions + 1L
    if (audit$predictions %% 3L == 2L) {
      if (audit$broken) p[[1]][1] <- NaN
      p <- p[rev(names(p))]
    }
    p
  }
  tailor_method <- function(object, ...) { audit$tailor <- audit$tailor + 1L; original_tailor(object, ...) }
  registerS3method("fit", "workflow", fit_method, envir = asNamespace("generics"))
  registerS3method("predict", "tailor", tailor_method, envir = asNamespace("stats"))
  on.exit({
    registerS3method("fit", "workflow", original_fit, envir = asNamespace("generics"))
    registerS3method("predict", "workflow", original_predict, envir = asNamespace("stats"))
    registerS3method("predict", "tailor", original_tailor, envir = asNamespace("stats"))
  }, add = TRUE)
  logistic <- probably::cal_estimate_logistic; multinomial <- probably::cal_estimate_multinomial
  audit$post <- 0L
  testthat::local_mocked_bindings(
    cal_estimate_logistic = function(...) { audit$post <- audit$post + 1L; logistic(...) },
    cal_estimate_multinomial = function(...) { audit$post <- audit$post + 1L; multinomial(...) }, .package = "probably")
  for (k in 2:4) for (method in c("inductive", "cross")) for (classification in c("ovo", "ovr")) {
    d <- ensemble_multiclass_data(k)
    spec <- calibrated_workflow(ensemble_multiclass_base(d), "truth", method = method, classification = classification,
      seed = 24, postprocessing = list(method = if (k == 2L) "tailor_logistic" else "tailor_multinomial", strategy = "ensemble", v = 3L))
    audit$fits <- audit$post <- 0L
    fit <- generics::fit(spec, d)
    count <- if (method == "cross") 15L else 3L
    testthat::expect_identical(audit$fits, count)
    testthat::expect_identical(audit$post, count)
    state <- workflows::extract_fit_engine(fit)
    testthat::expect_output(print(state), paste0("Models per nonempty prediction: ", count))
    unseen <- ensemble_multiclass_unseen(d)
    expected <- predict(fit, unseen, type = "prob")
    registerS3method("predict", "workflow", predict_method, envir = asNamespace("stats"))
    for (n in c(0L, 1L, nrow(unseen))) {
      audit$predictions <- audit$tailor <- 0L
      result <- predict(fit, unseen[seq_len(n), ], type = "prob")
      testthat::expect_identical(audit$predictions, if (n) count else 0L)
      testthat::expect_identical(audit$tailor, if (n) count else 0L)
      testthat::expect_equal(result, expected[seq_len(n), ], tolerance = 1e-10)
      testthat::expect_identical(audit$fits, count)
      testthat::expect_identical(audit$post, count)
    }
    ensemble <- state$models[[1]]
    for (remove in c(TRUE, FALSE)) {
      missing <- ensemble
      if (remove) missing$models[[2]] <- NULL else missing$models[2] <- list(NULL)
      testthat::expect_error(predict(missing, d, type = "prob"), "incomplete inner ensemble")
    }
    audit$broken <- TRUE; audit$predictions <- 0L
    testthat::expect_error(predict(fit, d, type = "prob"), "Fold2.*finite probabilities")
    audit$broken <- FALSE
    registerS3method("predict", "workflow", original_predict, envir = asNamespace("stats"))
  }
})

testthat::test_that("every shared inner pair enforces class, size, rank and convergence rules", {
  for (pkg in c("workflows", "tailor", "probably", "rsample", "recipes", "nnet")) testthat::skip_if_not_installed(pkg)
  d <- ensemble_multiclass_data(); base <- ensemble_multiclass_base(d)
  va <- function(data) rsample::make_splits(list(analysis = 1:180, assessment = 181:240), data)
  post <- list(method = "tailor_multinomial", strategy = "ensemble", splits = ensemble_multiclass_plan)
  spec <- calibrated_workflow(base, "truth", classification = "ovo", splits = va, seed = 24, postprocessing = post)
  for (classification in c("ovo", "ovr")) {
    for (k in 2:4) {
      bad <- post; bad$method <- if (k == 2L) "tailor_multinomial" else "tailor_logistic"
      testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", classification = classification,
        splits = va, postprocessing = bad), ensemble_multiclass_data(k)), "require explicitly selecting tailor_")
    }
    for (P in list(1:3, which(d$truth[1:180] != "b"), which(d$truth[1:180] == "b"))) {
      bad <- post
      bad$splits <- function(data) rsample::manual_rset(lapply(list(P, setdiff(1:180, P)), function(p)
        rsample::make_splits(list(analysis = setdiff(1:180, p), assessment = p), data)), c("bad", "rest"))
      testthat::expect_error(generics::fit(calibrated_workflow(base, "truth", classification = classification,
        splits = va, postprocessing = bad), d), "bad tailor_multinomial B/P plan")
    }
  }
  original <- probably::cal_estimate_multinomial
  for (fault in c("none", "missing", "rank", "nonconverged", "nonfinite")) local({
    calls <- 0L
    replacement <- function(...) {
      calls <<- calls + 1L
      if (calls == 2L && fault == "none") return(probably::cal_estimate_none(...))
      result <- original(...)
      if (calls == 2L) {
        if (fault == "missing") result$estimates <- list()
        if (fault == "rank") result$estimates[[1]]$estimate$rank <- 2L
        if (fault == "nonconverged") result$estimates[[1]]$estimate$convergence <- 1L
        if (fault == "nonfinite") result$estimates[[1]]$estimate$wts[] <- NA_real_
      }
      result
    }
    testthat::local_mocked_bindings(cal_estimate_multinomial = replacement, .package = "probably")
    testthat::expect_error(suppressWarnings(generics::fit(spec, d)), "Inductive1 / medium tailor_multinomial.*P")
    testthat::expect_identical(calls, 2L)
  })
  constant <- d; constant$x <- constant$z <- constant$w <- 1
  constant_base <- workflows::workflow() |> workflows::add_formula(truth ~ .) |>
    workflows::add_model(parsnip::multinom_reg(penalty = .1) |> parsnip::set_engine("nnet", trace = FALSE))
  testthat::expect_error(suppressWarnings(generics::fit(calibrated_workflow(constant_base, "truth", classification = "ovo",
    splits = va, seed = 24, postprocessing = post), constant)), "small tailor_multinomial.*full design rank 3")
})

testthat::test_that("shared inner B P and C mutations respect each local learning boundary", {
  for (pkg in c("workflows", "tailor", "probably", "rsample", "recipes", "nnet")) testthat::skip_if_not_installed(pkg)
  d <- ensemble_multiclass_data(); base <- ensemble_multiclass_base(d)
  for (method in c("inductive", "cross")) for (classification in c("ovo", "ovr")) {
    va <- function(data) {
      held <- list(1:60, 61:120, 121:240)
      parts <- lapply(held, function(c) rsample::make_splits(list(analysis = setdiff(1:240, c), assessment = c), data))
      if (method == "cross") rsample::manual_rset(parts, c("one", "two", "three")) else parts[[1]]
    }
    spec <- calibrated_workflow(base, "truth", method = method, classification = classification, splits = va, seed = 24,
      postprocessing = list(method = "tailor_multinomial", strategy = "ensemble", splits = function(data) {
        rsample::manual_rset(lapply(list(seq_len(nrow(data) / 2), seq.int(nrow(data) / 2 + 1, nrow(data))), function(p)
          rsample::make_splits(list(analysis = setdiff(seq_len(nrow(data)), p), assessment = p), data)), c("first", "second"))
      }))
    fit <- function(data) workflows::extract_fit_engine(generics::fit(spec, data))
    original <- fit(d)
    for (j in names(original$models[[1]]$models)) {
      p <- original$postprocessing[[1]]$inner_folds[[j]]
      model <- original$models[[1]]$models[[j]]
      for (role in c("B", "P", "C")) for (variable in c("x", "truth")) {
        changed <- d; rows <- p[[role]]
        if (variable == "x") changed$x[rows] <- changed$x[rows] * 2 + .1 else
          changed$truth[rows] <- levels(d$truth)[as.integer(changed$truth[rows]) %% 3L + 1L]
        after <- fit(changed); other <- after$models[[1]]$models[[j]]
        if (role == "B") {
          testthat::expect_false(identical(coef(workflows::extract_fit_engine(model)), coef(workflows::extract_fit_engine(other))))
        } else {
          testthat::expect_identical(coef(workflows::extract_fit_engine(model)), coef(workflows::extract_fit_engine(other)))
          testthat::expect_identical(workflows::extract_recipe(model)$steps, workflows::extract_recipe(other)$steps)
        }
        if (role == "P") testthat::expect_false(identical(p$coefficients, after$postprocessing[[1]]$inner_folds[[j]]$coefficients))
        if (role == "C") testthat::expect_identical(original$postprocessing[[1]], after$postprocessing[[1]])
      }
    }
  }
})
