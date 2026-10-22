testthat::test_that("local-plans", {
  skip_if_no_integration()
library(vennabers)
local({
  data <- data.frame(x = sin(1:80), truth = factor(rep(c("yes", "no", "no", "yes"), 20)))
  base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  split <- withr::with_seed(42, rsample::initial_split(data, .8, strata = truth, pool = 0))
  # Unselected data must not travel with a serialized prepared resolver.
  wide <- data
  wide$unused <- strrep("x", 1000000L)
  for (prepared in list(rsample::make_splits(list(analysis = 1:60, assessment = 61:80), wide),
                        withr::with_seed(42, rsample::vfold_cv(wide, v = 2)))) {
    compact <- prepared
    if (inherits(compact, "rset")) {
      compact$splits <- lapply(compact$splits, function(s) { s$data <- data; s })
    } else compact$data <- data
    resolver <- local_calibration_splits(list(prepared), names(data))
    reference <- local_calibration_splits(list(compact), names(data))
    serialized <- serialize(resolver, NULL)
    stopifnot(length(serialized) <= length(serialize(reference, NULL)) + 1024L)
    restored <- unserialize(serialized)
    stopifnot(identical(ls(environment(restored), all.names = TRUE), c("columns", "plans", "schemas", "snapshots")),
              identical(parent.env(environment(restored)), environment(local_calibration_splits)))
    actual <- restored(data)
    expected <- reference(data)
    if (!inherits(actual, "rset")) {
      actual <- list(splits = list(actual))
      expected <- list(splits = list(expected))
    }
    for (part in c("analysis", "assessment")) {
      stopifnot(identical(lapply(actual$splits, as.integer, data = part),
                          lapply(expected$splits, as.integer, data = part)))
    }
  }
  automatic <- generics::fit(calibrated_workflow(base, "truth", seed = 42), data)
  null_plan <- generics::fit(calibrated_workflow(base, "truth", splits = NULL, seed = 42), data)
  stopifnot(identical(venn_abers_details(automatic, data), venn_abers_details(null_plan, data)))
  explicit <- generics::fit(calibrated_workflow(base, "truth", splits = split), data)
  stopifnot(identical(venn_abers_details(automatic, data), venn_abers_details(explicit, data)))

  error <- function(expr, pattern) {
    result <- tryCatch(force(expr), error = identity)
    stopifnot(inherits(result, "error"), grepl(pattern, conditionMessage(result)))
  }
  for (arg in list(list(calibration_prop = .2), list(strata = NULL), list(strata = FALSE))) {
    error(do.call(calibrated_workflow, c(list(base, "truth", splits = split), arg)), "cannot be combined")
  }
  error(calibrated_workflow(base, "truth", splits = list(split)), "rsplit")
  for (seed in list(NA_real_, Inf, c(1, 2), .5, "1", .Machine$integer.max + 1)) {
    error(calibrated_workflow(base, "truth", seed = seed), "seed")
  }
  # The selected raw columns are canonical; unused columns/container/row names are not identity.
  foreign <- data
  foreign$x <- foreign$x + 1
  reordered <- data[80:1, ]
  changed_type <- data
  changed_type$x <- as.character(changed_type$x)
  changed_levels <- data
  levels(changed_levels$truth) <- rev(levels(data$truth))
  changed_target <- data
  changed_target$truth[1] <- data$truth[2]
  for (bad in list(foreign, reordered, changed_type, changed_levels, changed_target)) {
    error(generics::fit(calibrated_workflow(base, "truth", splits = split), bad), "data binding")
  }
  same <- tibble::as_tibble(data[c("truth", "x")])
  same$unused <- seq_len(nrow(same))
  fit_same <- generics::fit(calibrated_workflow(base, "truth", predictors = "x", splits = split), same)
  stopifnot(identical(predict(explicit, data, type = "prob"), predict(fit_same, data, type = "prob")))
  bad <- split
  bad$data <- NULL
  error(generics::fit(calibrated_workflow(base, "truth", splits = bad), data), "data binding")
  for (indices in list(c(1, 1), c(0, 1), c(1, 81), c(1, NA), c(1, Inf), c(1, 2.5), numeric())) {
    bad <- split
    bad$in_id <- indices
    bad$out_id <- 41:80
    error(generics::fit(calibrated_workflow(base, "truth", splits = bad), data), "partition")
  }
  for (held in list(40:80, 42:80)) {
    bad <- rsample::make_splits(list(analysis = 1:40, assessment = held), data)
    error(generics::fit(calibrated_workflow(base, "truth", splits = bad), data), "partition")
  }
  bad <- split
  class(bad) <- c("group_initial_split", class(bad))
  error(generics::fit(calibrated_workflow(base, "truth", splits = bad), data), "IID holdout")
  error(generics::fit(calibrated_workflow(base, "truth", splits = function(data) 1:10), data), "IID holdout")

  other <- data
  other$x <- other$x / 2
  other_split <- rsample::make_splits(list(analysis = 1:60, assessment = 61:80), other)
  resolver <- local_calibration_splits(list(split, other_split, split), c("x", "truth"))
  path <- tempfile()
  saveRDS(resolver, path)
  resolver <- readRDS(path)
  unlink(path)
  # Reverse lookup order and serialization do not change a prepared assignment.
  for (d in list(other, data, same)) {
    expected_split <- if (identical(d, other)) other_split else split
    resolved <- generics::fit(calibrated_workflow(base, "truth", predictors = "x", splits = resolver), d)
    expected <- generics::fit(calibrated_workflow(base, "truth", predictors = "x", splits = expected_split), d)
    stopifnot(identical(venn_abers_details(resolved, d), venn_abers_details(expected, d)))
  }
  error(resolver(foreign), "No prepared plan")
  conflict <- rsample::make_splits(list(analysis = 1:60, assessment = 61:80), data)
  error(local_calibration_splits(list(split, conflict), names(data)), "entries 1 and 2")
  # Order of positions is part of the plan, even if their sets agree.
  conflict <- split
  conflict$in_id <- rev(conflict$in_id)
  error(local_calibration_splits(list(split, conflict), names(data)), "Conflicting")
  # Check every entry, including entries unrelated to the first lookup.
  bad <- other_split
  bad$in_id <- c(1L, 1L)
  error(local_calibration_splits(list(split, bad), names(data)), "Prepared entry 2")
  error(local_calibration_splits(list(split), c("x", "missing")), "missing selected")
  error(generics::fit(calibrated_workflow(base, "truth", splits = resolver), data[1:60, ]), "No prepared plan")
  # Equal observations remain separate positions; row names are immaterial.
  duplicates <- data[rep(1:20, 2), ]
  duplicate_split <- rsample::make_splits(list(analysis = 1:20, assessment = 21:40), duplicates)
  rownames(duplicates) <- paste0("row", 1:40)
  duplicate_fit <- generics::fit(calibrated_workflow(base, "truth", splits = duplicate_split), duplicates)
  stopifnot(identical(workflows::extract_fit_engine(duplicate_fit)$folds$Inductive1$assessment, 21:40))

  # Prepared mappings must bind every selected column, including outcome.
  partial <- local_calibration_splits(list(split), "x")
  error(generics::fit(calibrated_workflow(base, "truth", splits = partial), changed_target), "canonical selected")
  partial <- local_calibration_splits(list(split), "truth")
  error(generics::fit(calibrated_workflow(base, "truth", splits = partial), foreign), "canonical selected")
  matrix_data <- data
  matrix_data$x <- cbind(sin(1:80), cos(1:80))
  matrix_split <- rsample::make_splits(list(analysis = 1:60, assessment = 61:80), matrix_data)
  matrix_resolver <- local_calibration_splits(list(matrix_split), c("x", "truth"))
  stopifnot(identical(rsample::analysis(matrix_resolver(matrix_data)), matrix_data[1:60, ]))

  # A stochastic plan and model formula share one RNG scope, including errors.
  stochastic <- workflows::workflow() |> workflows::add_formula(truth ~ jitter(x)) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  random_plan <- function(data) rsample::initial_split(data, .8)
  model <- calibrated_workflow(stochastic, "truth", splits = random_plan, seed = 123)
  set.seed(9)
  before <- .Random.seed
  one <- generics::fit(model, data)
  stopifnot(identical(.Random.seed, before))
  two <- generics::fit(model, data)
  stopifnot(identical(workflows::extract_fit_engine(one)$folds, workflows::extract_fit_engine(two)$folds),
            identical(workflows::extract_fit_engine(workflows::extract_fit_engine(one)$models$Inductive1)$coefficients,
                      workflows::extract_fit_engine(workflows::extract_fit_engine(two)$models$Inductive1)$coefficients))
  generics::fit(calibrated_workflow(base, "truth", splits = random_plan), data)
  stopifnot(!identical(.Random.seed, before))
  broken <- function(data) { runif(1); stop("planned failure") }
  before <- .Random.seed
  error(generics::fit(calibrated_workflow(base, "truth", splits = broken, seed = 123), data), "planned failure")
  stopifnot(identical(.Random.seed, before))
  rm(".Random.seed", envir = .GlobalEnv)
  generics::fit(model, data)
  stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

  testthat::expect_true(TRUE)
})
