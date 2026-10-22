testthat::test_that("pairwise postprocessing agrees with independent binary fits before coupling", {
  for (pkg in c("workflows", "recipes", "rsample", "tailor", "probably")) testthat::skip_if_not_installed(pkg)
  reference <- pairwise_reference("va_post_pair_reference")
  d <- data.frame(x = sin(1:240), row_id = 1:240,
    truth = factor(rep(c("z", "a", "b", "four", "z", "a"), 40),
      levels = c("z", "a", "b", "four")))
  base <- workflows::workflow() |> workflows::add_recipe(
    recipes::recipe(truth ~ x, d) |> recipes::step_normalize(x)) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  va <- function(data) rsample::make_splits(list(analysis = seq_len(floor(nrow(data) * .75)),
    assessment = seq.int(floor(nrow(data) * .75) + 1L, nrow(data))), data)
  for (k in 2:4) for (kind in c("tailor_logistic", "probability", "glm_link"))
    for (strategy in c("holdout", "ensemble")) for (method in c("inductive", "cross")) {
    local_data <- droplevels(d[d$truth %in% levels(d$truth)[seq_len(k)], ])
    post <- list(method = if (kind == "tailor_logistic") kind else "sigmoid", strategy = strategy)
    if (kind != "tailor_logistic") post$sigmoid_input <- kind
    if (strategy == "holdout") post$prop <- .3 else post$v <- 2L
    args <- list(base, "truth", classification = "ovo", model_training = "pairwise",
      method = method, seed = 24, precision = 2L, postprocessing = post)
    if (method == "inductive") args$splits <- va else args$v <- 2L
    fit <- generics::fit(do.call(calibrated_workflow, args), local_data)
    state <- workflows::extract_fit_engine(fit)
    new <- d[1:7, c("x", "row_id")]
    q <- list(log = list(), brier = list())
    for (pair_id in names(state$pair_models)) {
      component <- state$pair_models[[pair_id]]
      local <- local_data[component$rows, ]; local$truth <- factor(local$truth, levels = component$class_levels)
      binary_base <- workflows::workflow() |> workflows::add_recipe(
        recipes::recipe(truth ~ x, local) |> recipes::step_normalize(x)) |>
        workflows::add_model(parsnip::logistic_reg(engine = "glm"))
      column <- paste0(".pred_", component$event_class)
      held <- predicted <- truth <- ids <- list()
      for (fold in names(component$folds)) {
        diagnostic <- component$postprocessing[[fold]]
        members <- if (strategy == "ensemble") diagnostic$inner_folds else list(Holdout1 = diagnostic)
        predictions <- calibration <- list()
        for (j in names(members)) {
          p <- members[[j]]
          inner <- post; inner$prop <- inner$v <- NULL; inner$strategy <- "holdout"
          inner$splits <- rsample::make_splits(p[c("analysis", "assessment")], local[p$Q, ])
          outer <- rsample::make_splits(list(analysis = p$Q, assessment = p$C), local)
          ref <- generics::fit(calibrated_workflow(binary_base, "truth", event_level = "second",
            splits = outer, postprocessing = inner), local)
          complete <- workflows::extract_fit_engine(ref)$models[[1]]
          predictions[[j]] <- predict(complete, new, type = "prob")[[column]]
          calibration[[j]] <- predict(complete, local[p$C, ], type = "prob")[[column]]
        }
        predicted[[fold]] <- unname(rowMeans(do.call(cbind, predictions)))
        held[[fold]] <- unname(rowMeans(do.call(cbind, calibration)))
        truth[[fold]] <- local$truth[diagnostic$C] == component$event_class
        ids[[fold]] <- rep(fold, length(diagnostic$C))
        testthat::expect_equal(component$calibration[[fold]]$probability, held[[fold]], tolerance = 1e-10)
        testthat::expect_equal(unname(predict(component$models[[fold]], new, type = "prob")[[column]]), predicted[[fold]], tolerance = 1e-10)
      }
      core <- if (method == "cross") cross_venn_abers(unlist(truth), unlist(held), unlist(ids), precision = 2L) else
        venn_abers(truth[[1]], held[[1]], precision = 2L)
      for (loss in names(q)) q[[loss]][[pair_id]] <- predict(core,
        if (method == "cross") do.call(cbind, predicted) else predicted[[1]], loss = loss)$probability
    }
    for (loss in names(q)) testthat::expect_equal(unname(as.matrix(predict(fit, new, type = "prob", loss = loss))),
      unname(reference$couple(do.call(cbind, q[[loss]]), levels(local_data$truth))), tolerance = 1e-10)
  }
})

testthat::test_that("pair-local diagnostics retain encoded identities and unequal inner plans", {
  for (pkg in c("workflows", "recipes", "rsample", "tailor", "probably")) testthat::skip_if_not_installed(pkg)
  raw <- data.frame(x = sin(1:240), truth = factor(rep(c("class", "a b", "a.b", "class"), 60),
    levels = c("class", "a b", "a.b")))
  names(raw) <- c("input value", ".pred_class")
  schema <- venn_abers_schema(raw, ".pred_class", classification = "ovo")
  d <- venn_abers_encode(raw, schema); outcome <- schema$encoded_outcome
  base <- workflows::workflow() |> workflows::add_recipe(
    recipes::recipe(reformulate(setdiff(names(d), outcome), outcome), d) |>
      recipes::step_normalize(recipes::all_numeric_predictors())) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  va <- function(data) {
    stopifnot(nlevels(data[[outcome]]) == 2L)
    n <- nrow(data)
    held <- list(seq_len(n %/% 3L), seq.int(n %/% 3L + 1L, n))
    rsample::manual_rset(lapply(held, function(h) rsample::make_splits(
      list(analysis = setdiff(seq_len(n), h), assessment = h), data)), c("small", "large"))
  }
  prepared <- list()
  inner <- function(data) {
    stopifnot(nlevels(data[[outcome]]) == 2L)
    plan <- withr::with_seed(37, rlang::inject(rsample::vfold_cv(data,
      v = if (nrow(data) > 80L) 3L else 2L, strata = !!outcome, pool = 0)))
    prepared[[length(prepared) + 1L]] <<- plan
    plan
  }
  for (kind in c("tailor_logistic", "probability", "glm_link")) {
    post <- list(method = if (kind == "tailor_logistic") kind else "sigmoid", strategy = "ensemble", splits = inner)
    if (kind != "tailor_logistic") post$sigmoid_input <- kind
    spec <- calibrated_workflow(base, outcome, schema = schema, method = "cross", model_training = "pairwise",
      splits = va, postprocessing = post)
    prepared <- list()
    fit <- generics::fit(spec, d); state <- workflows::extract_fit_engine(fit)
    new <- d[setdiff(names(d), outcome)]
    details <- venn_abers_details(fit, new)
    decoded <- venn_abers_decode(details, schema)
    testthat::expect_identical(names(details$postprocessing),
      c("pair_1/small", "pair_1/large", "pair_2/small", "pair_2/large", "pair_3/small", "pair_3/large"))
    testthat::expect_identical(unname(lengths(lapply(details$postprocessing, `[[`, "inner_folds"))), c(3L, 2L, 3L, 2L, 2L, 2L))
    for (id in names(details$postprocessing)) {
      p <- details$postprocessing[[id]]
      testthat::expect_identical(p$pair_rows, state$pair_models[[p$pair_id]]$rows)
      testthat::expect_identical(sort(unlist(lapply(p$inner_folds, `[[`, "P"), use.names = FALSE)), sort(p$Q))
      for (j in names(p$inner_folds)) {
        native <- p$inner_folds[[j]]; display <- decoded$postprocessing[[id]]$inner_folds[[j]]
        testthat::expect_identical(native$schema, schema)
        testthat::expect_identical(native$pair_id, p$pair_id)
        testthat::expect_identical(display$classes, schema$class_levels[match(native$classes, schema$classes)])
        testthat::expect_identical(display$column_types[[schema$outcome]]$levels, display$classes)
        testthat::expect_identical(display$columns, names(raw))
      }
    }
    for (n in c(0L, 1L, 7L)) {
      batch <- new[seq_len(n), , drop = FALSE]
      small <- venn_abers_details(fit, batch)
      testthat::expect_identical(names(small$details), names(details$details))
      testthat::expect_identical(small$prediction, predict(fit, batch, type = "prob"))
      testthat::expect_equal(nrow(small$details), n * 6L)
      if (n) {
        testthat::expect_identical(small$details$inner_fold[[1]], names(details$postprocessing[[1]]$inner_folds))
        testthat::expect_identical(small$details$B_size[[1]], vapply(details$postprocessing[[1]]$inner_folds, `[[`, 1L, "B_size"))
      }
    }
    testthat::expect_output(print(state), "Base fits: 14")
    testthat::expect_output(print(state), "Models per nonempty prediction: 14")
    # Data-bound prepared plans replay in reverse pair/fold call order.
    resolver <- unserialize(serialize(local_calibration_splits(prepared, names(d)), NULL))
    for (plan in rev(prepared)) testthat::expect_identical(
      unname(lapply(resolver(plan$splits[[1]]$data)$splits, as.integer, "assessment")),
      lapply(plan$splits, as.integer, "assessment"))
    post$splits <- resolver
    replay <- generics::fit(calibrated_workflow(base, outcome, schema = schema, method = "cross",
      model_training = "pairwise", splits = function(data) va(data)[2:1, ], postprocessing = post), d)
    testthat::expect_equal(predict(fit, new, type = "prob"), predict(replay, new, type = "prob"), tolerance = 1e-10)
  }
})
