testthat::test_that("frozen-schema", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))

packages <- c("parsnip", "workflows", "hardhat", "rsample", "recipes",
              "rlang", "tibble", "withr", "generics", "tidyselect")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) local({
  original <- data.frame(x = sin(1:100), z = cos(1:100),
                         y = factor(rep(c("class", "vaL0001", "vaL0001", "class"), 25),
                                    levels = c("class", "vaL0001")))
  names(original) <- c("vaV0002", "a b", ".pred_class")
  schema <- venn_abers_schema(original[1:40, ], ".pred_class", event_level = "second")
  encoded <- venn_abers_encode(original, schema)
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(vaV0003 ~ ., encoded[1:40, ]) |>
      recipes::step_normalize(recipes::all_numeric_predictors())) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  frozen <- generics::fit(base, encoded[1:40, ])
  fit <- calibrate_workflow(frozen, encoded[41:80, ], schema$encoded_outcome, schema = schema)
  new <- encoded[81:100, 1:2]
  native <- predict(fit, new)
  decoded <- venn_abers_decode(venn_abers_details(fit, new), schema)
  manual <- venn_abers(encoded$vaV0003[41:80] == "vaL0002",
                       predict(frozen, encoded[41:80, 1:2], type = "prob")$.pred_vaL0002)
  assert_conditions(identical(fit$models[[1L]], frozen), identical(fit$schema, schema),
            identical(fit$calibrator, manual), fit$event_class == "vaL0002",
            identical(attr(native, "vennabers_schema"), schema),
            identical(decoded$probability[, "vaL0001"], native$.pred_vaL0002),
            all(decoded$details$event_class == "vaL0001"),
            is.null(decoded$identity), is.null(decoded$aggregate))
  reordered <- encoded[41:80, 3:1]
  reordered$vaV0003 <- factor(reordered$vaV0003, levels = rev(levels(reordered$vaV0003)))
  stopifnot(identical(calibrate_workflow(frozen, reordered, schema$encoded_outcome, schema = schema)$calibrator,
                      fit$calibrator))
})

if (all(vapply(c(packages, "nnet"), requireNamespace, logical(1), quietly = TRUE))) local({
  error <- function(expr, pattern) {
    e <- tryCatch(force(expr), error = identity)
    stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e)))
  }
  close <- function(actual, expected, regression = FALSE) {
    assert_conditions(identical(dim(actual), dim(expected)), length(actual) == length(expected),
              all(is.finite(actual)), all(is.finite(expected)),
              all(abs(actual - expected) <= if (regression) 1e-10 + 1e-10 * abs(expected) else 1e-10))
  }
  outputs <- function(fit, new) {
    native <- if (fit$schema$mode == "regression") list(predict(fit, new), venn_abers_details(fit, new)) else
      unlist(lapply(c("log", "brier"), function(loss) list(
        predict(fit, new, loss = loss), predict(fit, new, type = "class", loss = loss),
        venn_abers_details(fit, new, loss = loss))), recursive = FALSE)
    list(native = native, decoded = lapply(native, venn_abers_decode, schema = fit$schema))
  }
  post_available <- all(vapply(c("tailor", "probably", "mgcv"), requireNamespace, logical(1), quietly = TRUE))
  guards <- list(c("stats", "glm.fit"), c("stats", "lm.fit"), c("nnet", "multinom"),
                 c("workflows", "fit.workflow"), c("recipes", "prep.recipe"))
  if (post_available) guards <- c(guards, list(c("tailor", "fit.tailor"), c("mgcv", "gam")))
  saved <- list()
  for (method in c("binary", "ovo", "ovr", "regression")) {
    for (k in if (method == "regression") 1L else if (method == "binary") 2L else 2:4) {
      for (event in if (method == "binary") c("first", "second") else "first") {
        for (post in if (post_available) c(FALSE, TRUE) else FALSE) local({
          set.seed(905 + k)
          original <- data.frame(x = rnorm(360), z = rnorm(360))
          original$y <- if (k == 1L) 2 + original$x + cos(1:360) else
            factor(rep(c("class", "vaL0001", "a b", ".pred_yes")[seq_len(k)], length.out = 360),
                   levels = c("class", "vaL0001", "a b", ".pred_yes")[seq_len(k)])
          names(original) <- c("vaV0002", "a-b", ".pred_class")
          schema <- if (k == 1L) venn_abers_schema(original[1:120, ], ".pred_class") else
            if (method == "binary") venn_abers_schema(original[1:120, ], ".pred_class", event_level = event) else
              venn_abers_schema(original[1:120, ], ".pred_class", classification = method)
          encoded <- venn_abers_encode(original, schema)
          spec <- if (k == 1L) parsnip::linear_reg(engine = "lm") else if (k == 2L)
            parsnip::logistic_reg(engine = "glm") else
              parsnip::multinom_reg(penalty = .1) |> parsnip::set_engine("nnet", trace = FALSE)
          base <- workflows::workflow() |>
            workflows::add_recipe(recipes::recipe(vaV0003 ~ ., encoded[1:120, ]) |>
              recipes::step_normalize(recipes::all_numeric_predictors())) |>
            workflows::add_model(spec)
          frozen <- generics::fit(base, encoded[1:120, ])
          if (post) {
            base <- workflows::add_tailor(base, if (k == 1L) tailor::tailor() |>
              tailor::adjust_numeric_calibration(method = "linear", smooth = FALSE) else tailor::tailor() |>
                tailor::adjust_probability_calibration(method = if (k == 2L) "logistic" else "multinomial", smooth = (k > 2L)))
            with_post <- generics::fit(base, encoded[1:120, ], data_calibration = encoded[121:220, ])
            type <- if (k == 1L) "numeric" else "prob"
            stopifnot(max(abs(as.matrix(predict(frozen, encoded[321:360, ], type = type)) -
                              as.matrix(predict(with_post, encoded[321:360, ], type = type)))) > 1e-5)
            frozen <- with_post
          }
          for (g in guards) trace(g[2L], where = asNamespace(g[1L]), tracer = quote(stop("Unexpected learning")), print = FALSE)
          on.exit(for (g in guards) untrace(g[2L], where = asNamespace(g[1L])), add = TRUE)
          calibration <- original[221:320, 3:1]
          if (k > 1L) calibration$.pred_class <- factor(calibration$.pred_class, levels = rev(schema$class_levels))
          calibration <- venn_abers_encode(calibration, schema)[3:1]
          new <- venn_abers_encode(original[321:360, 2:1], schema)[2:1]
          final_cal <- predict(frozen, calibration[2:3], type = if (k == 1L) "numeric" else "prob")
          final_new <- predict(frozen, new, type = if (k == 1L) "numeric" else "prob")
          if (k > 1L) {
            final_cal <- as.matrix(final_cal); final_new <- as.matrix(final_new)
            colnames(final_cal) <- colnames(final_new) <- unname(schema$classes)
          }
          options <- if (k == 1L) list(list(m = 1), list(epsilon = .1)) else
            lapply(list(NULL, 0, 2, -1), function(p) list(precision = p, loss = "brier"))
          for (option in options) {
            fit <- do.call(calibrate_workflow, c(list(frozen, calibration, schema$encoded_outcome, schema = schema), option))
            manual <- if (k == 1L) do.call(venn_abers_regression, c(list(calibration$vaV0003, final_cal$.pred), option)) else
              if (method == "binary") venn_abers(calibration$vaV0003 == fit$event_class,
                final_cal[, fit$event_class], precision = option$precision) else
                  venn_abers(calibration$vaV0003, final_cal, method = method, precision = option$precision)
            assert_conditions(identical(fit$calibrator, manual), identical(fit$models[[1L]], frozen), identical(fit$schema, schema),
                      is.null(fit$folds), !any(vapply(fit[setdiff(names(fit), "models")], is.data.frame, logical(1))))
            for (loss in if (k == 1L) "log" else c("log", "brier")) {
              native <- if (k == 1L) venn_abers_details(fit, new) else venn_abers_details(fit, new, loss = loss)
              decoded <- venn_abers_decode(native, schema)
              stopifnot(is.null(decoded$identity), is.null(decoded$aggregate), decoded$outcome == ".pred_class",
                        identical(attr(native$prediction, "vennabers_schema"), schema))
              if (k == 1L) {
                expected <- predict(manual, final_new$.pred)
                close(decoded$prediction, expected$prediction, TRUE)
                close(as.matrix(decoded$details[c("prediction", "lower", "upper")]), as.matrix(expected), TRUE)
                stopifnot(is.null(decoded$probability), is.null(decoded$class))
              } else {
                expected <- if (method == "binary") predict(manual, final_new[, fit$event_class], loss = loss) else
                  predict(manual, final_new, type = "details", loss = loss)
                if (method == "binary") {
                  close(decoded$probability[, schema$event_class], expected$probability)
                  parts <- expected
                  stopifnot(all(decoded$details$event_class == schema$event_class))
                } else {
                  close(unname(decoded$probability), unname(expected$prediction))
                  parts <- if (method == "ovo") expected$pairs else expected$classes
                  stopifnot(identical(decoded$details$event_class, schema$class_levels[match(parts$event, schema$classes)]))
                  if (method == "ovo") stopifnot(identical(decoded$details$class_1, schema$class_levels[match(parts$non_event, schema$classes)]),
                                                 identical(decoded$details$class_2, decoded$details$event_class)) else
                    stopifnot(identical(decoded$details$class, decoded$details$event_class))
                }
                close(as.matrix(decoded$details[c("lower", "upper", "probability")]),
                      as.matrix(parts[c("lower", "upper", "probability")]))
                stopifnot(identical(decoded$details[c("row", "fold")], native$details[c("row", "fold")]),
                          identical(colnames(decoded$probability), schema$class_levels))
              }
              table <- native$prediction
              if (k > 1L) table$.pred_class <- predict(fit, new, type = "class", loss = loss)$.pred_class
              table[[schema$encoded_outcome]] <- encoded$vaV0003[321:360]
              table$.row <- 321:360; table$id <- "assessment"; table$.config <- "frozen"
              all_roles <- venn_abers_decode(table, schema)
              stopifnot(identical(all_roles$truth, original$.pred_class[321:360]),
                        identical(names(all_roles$identity), c(".row", "id", ".config")))
              if (k > 1L) stopifnot(identical(levels(all_roles$class), schema$class_levels))
            }
            for (n in 0:1) {
              result <- outputs(fit, new[seq_len(n), , drop = FALSE])
              de <- result$decoded[[2L]]
              stopifnot(if (k == 1L) length(de$prediction) == n else length(de$class) == n)
            }
          }
          # Persist the final option of each trained procedure for a fresh-process check.
          saved[[length(saved) + 1L]] <<- list(fit = fit, new = new, expected = outputs(fit, new),
            empty = outputs(fit, new[0, ]), single = outputs(fit, new[1, , drop = FALSE]),
            parameters = fit$parameters, schema = schema)
          reordered <- calibration
          if (k > 1L) reordered$vaV0003 <- factor(reordered$vaV0003, levels = rev(levels(reordered$vaV0003)))
          same <- do.call(calibrate_workflow, c(list(frozen, reordered, schema$encoded_outcome, schema = schema), option))
          stopifnot(identical(same$calibrator, fit$calibrator),
                    identical(predict(fit, new), predict(fit, new[2:1])))
          # All detectable schema errors must precede calibrator learning.
          for (name in c("venn_abers", "venn_abers_regression")) trace(name, where = asNamespace("vennabers"),
            tracer = quote(stop("Unexpected calibration learning")), print = FALSE)
          on.exit(for (name in c("venn_abers", "venn_abers_regression")) untrace(name, where = asNamespace("vennabers")), add = TRUE)
          foreign <- schema
          names(foreign$columns)[1L] <- foreign$original_columns[1L] <- "another predictor"
          bound <- calibration; attr(bound, "vennabers_schema") <- foreign
          error(calibrate_workflow(frozen, bound, schema$encoded_outcome, schema = schema), "different schema")
          error(calibrate_workflow(frozen, calibration, schema$encoded_outcome, schema = foreign), "different schema")
          error(calibrate_workflow(frozen, calibration, schema$encoded_outcome), "different schema|exactly two")
          error(calibrate_workflow(frozen, calibration, ".pred_class", schema = schema), "schema")
          error(calibrate_workflow(frozen, calibration[-1L], schema$encoded_outcome, schema = schema), "encoded columns")
          error(calibrate_workflow(frozen, calibration[-2L], schema$encoded_outcome, schema = schema), "encoded columns")
          error(calibrate_workflow(frozen, calibration[0, ], schema$encoded_outcome, schema = schema), "nonempty")
          for (bad in list(bound, calibration[-2L], original[321:360, ], transform(new, unknown = 1))) {
            error(predict(fit, bad), "schema|encoded columns")
            error(venn_abers_details(fit, bad), "schema|encoded columns")
          }
          error(venn_abers_decode(predict(fit, new), foreign), "different schema")
          error(venn_abers_decode(venn_abers_details(fit, new), foreign), "different schema")
          error(calibrate_workflow(frozen, calibration, schema$encoded_outcome, schema = schema,
                                   classification = if (method == "ovo") "ovr" else "ovo"), "schema")
          error(calibrate_workflow(frozen, calibration, schema$encoded_outcome, schema = schema,
                                   event_level = if (event == "first") "second" else "first"), "schema")
          bad <- calibration
          bad$vaV0003[1L] <- NA
          error(calibrate_workflow(frozen, bad, schema$encoded_outcome, schema = schema), "nonmissing|finite")
          if (k > 1L) {
            bad$vaV0003 <- factor(rep("unknown", nrow(bad)))
            error(calibrate_workflow(frozen, bad, schema$encoded_outcome, schema = schema), "class order")
            bad <- original[221:320, ]; bad$.pred_class <- factor(rep("unknown", 100))
            error(venn_abers_encode(bad, schema), "Unknown outcome classes")
          }
        })
      }
    }
  }
  path <- tempfile(fileext = ".rds"); script <- tempfile(fileext = ".R"); log <- tempfile(fileext = ".log")
  on.exit(unlink(c(path, script, log)), add = TRUE)
  saveRDS(saved, path)
  writeLines(c(paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"), "library(vennabers)",
    paste0("outputs <- ", paste(deparse(outputs), collapse = "\n")),
    # The main process and each worker restore the saved fits under the same learning guards.
    paste0("guard <- function() {", paste(vapply(guards, function(g) paste0("trace(", deparse(g[2L]),
      ", where = asNamespace(", deparse(g[1L]), "), tracer = quote(stop('Unexpected learning')), print = FALSE)"), ""),
      collapse = "; "), "; for (name in c('venn_abers', 'venn_abers_regression')) trace(name, where = asNamespace('vennabers'),",
      " tracer = quote(stop('Unexpected calibration learning')), print = FALSE) }"),
    "guard()",
    paste0("saved <- readRDS(", deparse(path), ")"),
    'for (s in saved) stopifnot(identical(s$expected, outputs(s$fit, s$new)), identical(s$empty, outputs(s$fit, s$new[0, ])), identical(s$single, outputs(s$fit, s$new[1, , drop = FALSE])), identical(s$fit$schema, s$schema), identical(s$fit$parameters, s$parameters))',
    'if (requireNamespace("future", quietly = TRUE)) {',
    '  future::plan(future::multisession, workers = 2)',
    '  jobs <- lapply(saved, function(s) future::future({ suppressMessages(guard()); list(pid = Sys.getpid(), expected = outputs(s$fit, s$new), empty = outputs(s$fit, s$new[0, ]), single = outputs(s$fit, s$new[1, , drop = FALSE])) }, globals = list(s = s, outputs = outputs, guard = guard), packages = "vennabers", seed = NULL))',
    '  found <- lapply(jobs, future::value)',
    '  future::plan(future::sequential)',
    '  pids <- vapply(found, function(x) x$pid, 0L)',
    '  stopifnot(!any(pids == Sys.getpid()), length(unique(pids)) == 2L)',
    '  for (i in seq_along(saved)) stopifnot(identical(found[[i]][c("expected", "empty", "single")], saved[[i]][c("expected", "empty", "single")]))',
    '}'), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = log, stderr = log)
  if (status != 0L) stop(paste(readLines(log), collapse = "\n"))
})

if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) local({
  error <- function(expr, pattern) {
    e <- tryCatch(force(expr), error = identity)
    stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e)))
  }
  original <- data.frame(x = sin(1:100), z = cos(1:100),
                         y = factor(rep(c("yes", "no", "no", "yes"), 25), levels = c("yes", "no")))
  names(original) <- c("vaV0002", "a b", ".pred_class")
  schema <- venn_abers_schema(original[1:40, ], ".pred_class")
  encoded <- venn_abers_encode(original, schema)
  legacy <- workflows::workflow() |>
    workflows::add_formula(`.pred_class` ~ `vaV0002` + `a b`) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  legacy <- generics::fit(legacy, original[1:40, ])
  adopted <- calibrate_workflow(legacy, original[41:80, ], ".pred_class")
  stopifnot(identical(adopted$models[[1L]], legacy), is.null(adopted$schema),
            identical(names(predict(adopted, original[81:100, ])), c(".pred_yes", ".pred_no")))
  error(calibrate_workflow(legacy, encoded[41:80, ], schema$encoded_outcome, schema = schema), "outcome.*pred_class")
  error(calibrate_workflow(legacy, original[41:80, ], ".pred_class", schema = schema), "schema")
  error(predict(adopted, encoded[81:100, ]), "different schema")
  detached <- encoded; attr(detached, "vennabers_schema") <- NULL
  error(predict(adopted, detached[81:100, ]), "missing original column.*vaV0002|missing original column.*a b")
  error(venn_abers_details(adopted, detached[81:100, ]), "missing original column")

  # Formula fits can lose bindings. Validate every piece of retained metadata,
  # while documenting that matching detached codes cannot establish provenance.
  base <- workflows::workflow() |> workflows::add_formula(vaV0003 ~ .) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  models <- list()
  changed <- detached
  changed$vaV0003 <- factor(changed$vaV0003, levels = rev(levels(changed$vaV0003)))
  models[[1L]] <- generics::fit(base, changed[1:40, ])
  changed$vaV0003 <- ordered(detached$vaV0003)
  models[[2L]] <- generics::fit(base, changed[1:40, ])
  models[[3L]] <- generics::fit(workflows::update_formula(base, vaV0003 ~ vaV0001), detached[1:40, ])
  changed$vaV0003 <- sin(1:100)
  models[[4L]] <- generics::fit(workflows::update_model(base, parsnip::linear_reg(engine = "lm")), changed[1:40, ])
  correct <- generics::fit(base, detached[1:40, ])
  foreign <- correct; attr(foreign, "vennabers_schema") <- venn_abers_schema(original, ".pred_class", event_level = "second")
  for (name in c("venn_abers", "venn_abers_regression")) trace(name, where = asNamespace("vennabers"),
    tracer = quote(stop("Unexpected calibration learning")), print = FALSE)
  on.exit(for (name in c("venn_abers", "venn_abers_regression")) untrace(name, where = asNamespace("vennabers")), add = TRUE)
  for (model in models) error(calibrate_workflow(model, encoded[41:80, ], schema$encoded_outcome, schema = schema), "must agree with.*schema")
  error(calibrate_workflow(foreign, encoded[41:80, ], schema$encoded_outcome, schema = schema), "different schema")
})

if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) local({
  # A real legacy GLM whose reserved class breaks the public workflow path.
  # A raw GLM response must never become an automatic product fallback.
  d <- data.frame(x = sin(1:100), y = factor(rep(c("class", "other", "other", "class"), 25)))
  base <- workflows::workflow() |> workflows::add_formula(y ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  frozen <- generics::fit(base, d[1:40, ])
  raw <- predict(frozen, d[41:80, "x", drop = FALSE], type = "raw", opts = list(type = "response"))
  stopifnot(length(raw) == 40L, all(is.finite(raw)))
  e <- tryCatch(calibrate_workflow(frozen, d[41:80, ], "y"), error = identity)
  stopifnot(inherits(e, "error"), grepl("Frozen calibration.*class.*x", conditionMessage(e)),
            grepl("probability|unique|duplicat", conditionMessage(e)))
})

if (all(vapply(c(packages, "nnet"), requireNamespace, logical(1), quietly = TRUE))) local({
  # An ordered trained outcome keeps ordered decisions, like the original workflow.
  # Calibration factors are matched by level identity, whatever their order or status.
  set.seed(906)
  for (k in 2:3) for (method in if (k == 2L) c("binary", "ovo", "ovr") else c("ovo", "ovr")) {
    for (encoded in c(FALSE, TRUE)) {
      classes <- c("low", "mid", "high")[seq_len(k)]
      d <- data.frame(x = rnorm(300), z = rnorm(300),
                      truth = factor(rep(classes, length.out = 300), levels = classes, ordered = TRUE))
      schema <- if (encoded) venn_abers_schema(d[1:100, ], "truth", classification = method) else NULL
      data <- if (encoded) venn_abers_encode(d, schema) else d
      outcome <- if (encoded) schema$encoded_outcome else "truth"
      spec <- if (k == 2L) parsnip::logistic_reg(engine = "glm") else
        parsnip::multinom_reg(penalty = .1) |> parsnip::set_engine("nnet", trace = FALSE)
      frozen <- generics::fit(workflows::workflow() |>
        workflows::add_recipe(recipes::recipe(stats::as.formula(paste(outcome, "~ .")), data[1:100, ])) |>
        workflows::add_model(spec), data[1:100, ])
      calibration <- data[101:200, ]
      reordered <- calibration
      reordered[[outcome]] <- factor(reordered[[outcome]], levels = rev(levels(reordered[[outcome]])))
      fit <- calibrate_workflow(frozen, calibration, outcome, classification = method, schema = schema)
      same <- calibrate_workflow(frozen, reordered, outcome, classification = method, schema = schema)
      new <- data[201:300, setdiff(names(data), outcome)]
      native <- predict(frozen, new, type = "class")$.pred_class
      assert_conditions(is.ordered(native), identical(same$calibrator, fit$calibrator))
      for (n in c(0L, 1L, nrow(new))) {
        decision <- predict(fit, new[seq_len(n), ], type = "class")$.pred_class
        assert_conditions(is.ordered(decision), identical(levels(decision), levels(native)), length(decision) == n)
        if (encoded) {
          decoded <- venn_abers_decode(predict(fit, new[seq_len(n), ], type = "class"), schema)$class
          assert_conditions(is.ordered(decoded), identical(levels(decoded), classes))
        }
      }
    }
  }
})

  testthat::expect_true(TRUE)
})
