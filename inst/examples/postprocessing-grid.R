# Grid selection with learned postprocessing. Every engine/recipe search runs
# inside its own B: each VA source, each inner ensemble pair and each full F
# composition searches separately. P, C, the outer assessment A and a reserved
# option-selection share S never enter a search. Joint candidates are complete
# fixed compositions (no nested search); the first two-stage step scores the raw
# engine/recipe without postprocessor or VA.
source(system.file("examples", "assertions.R", package = "vennabers"))
packages <- c("vennabers", "parsnip", "withr", "workflows", "recipes", "rsample", "tune", "yardstick", "tailor", "probably", "nnet")
cases <- list()
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  library(vennabers)
  checker <- new.env()
  sys.source(system.file("examples", "check-encoded-evaluation.R", package = "vennabers"), checker)
  stratified <- function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)
  plain <- function(d) rsample::vfold_cv(d, v = 2)
  grid <- function(n) list(method = "grid", grid = data.frame(degree = seq_len(n)))
  classes <- withr::with_seed(7, {
    truth <- factor(rep(c("a", "b", "c"), 800))
    data.frame(x = stats::rnorm(2400, sd = 3) + c(-1, 1, -1)[truth],
      z = stats::rnorm(2400, sd = 3) + c(-1, -1, 1)[truth], truth = truth)
  })
  binary <- withr::with_seed(8, {
    x <- stats::rnorm(1200, sd = 1.5)
    data.frame(x = x, z = stats::rnorm(1200),
      truth = factor(ifelse(stats::runif(1200) < stats::plogis(x - .4 * x^2 + .5), "yes", "no")))
  })
  numeric <- withr::with_seed(9, {
    x <- stats::rnorm(1200)
    data.frame(x = x, z = stats::rnorm(1200), truth = 3 + x - .5 * x^2 + stats::rnorm(1200, sd = .5))
  })
  base <- function(d, engine) workflows::workflow() |> workflows::add_recipe(
    recipes::recipe(truth ~ x + z, d) |> recipes::step_poly(x, degree = tune::tune("degree"))) |>
    workflows::add_model(engine)
  glm <- parsnip::logistic_reg(engine = "glm")
  nnet <- parsnip::multinom_reg(penalty = .1) |> parsnip::set_engine("nnet", trace = FALSE, maxit = 1000)
  specs <- list(
    # The v0.9 benchmark case: binary full ensemble, two stages, four candidates per local search.
    binary_sigmoid_ensemble = list(data = binary, spec = calibrated_workflow(base(binary, glm), "truth",
      method = "cross", v = 2, prediction_model = "full", seed = 3,
      tuning = list(strategy = "two_stage", model = c(grid(4), list(resamples = stratified))),
      postprocessing = list(method = "sigmoid", sigmoid_input = "probability", strategy = "ensemble", v = 3))),
    pairwise_tailor_joint = list(data = classes[c("x", "z", "truth")], spec = calibrated_workflow(base(classes, glm), "truth",
      classification = "ovo", model_training = "pairwise", seed = 4,
      tuning = list(strategy = "joint", model = c(grid(2), list(resamples = stratified)), joint_resamples = stratified),
      postprocessing = list(method = "tailor_logistic", strategy = "holdout", prop = .3))),
    ovr_multinomial_options = list(data = classes, spec = calibrated_workflow(base(classes, nnet), "truth",
      classification = "ovr", seed = 5, loss = tune::tune("va_loss"),
      tuning = list(strategy = "two_stage", selection_prop = .5,
        model = c(grid(2), list(resamples = stratified)),
        calibration = list(method = "grid", grid = data.frame(va_loss = c("log", "brier")), resamples = stratified)),
      postprocessing = list(method = "tailor_multinomial", strategy = "holdout", prop = .3))),
    regression_isotonic_joint = list(data = numeric, spec = calibrated_workflow(base(numeric, parsnip::linear_reg()), "truth",
      method = "cross", v = 2, prediction_model = "full", seed = 6,
      tuning = list(strategy = "joint", model = c(grid(2), list(resamples = plain))),
      postprocessing = list(method = "isotonic_numeric", strategy = "holdout", prop = .3))),
    regression_linear_ensemble = list(data = numeric, spec = calibrated_workflow(base(numeric, parsnip::linear_reg()), "truth",
      seed = 7, tuning = list(strategy = "two_stage", model = c(grid(2), list(resamples = plain))),
      postprocessing = list(method = "tailor_linear", strategy = "ensemble", v = 2))))
  # Every local search with its inner B/P plan and, for pairwise, its pair rows.
  searches <- function(state) {
    parts <- if (is.null(state$pair_models)) list(state) else state$pair_models
    unlist(lapply(parts, function(part) {
      local <- if (is.null(part$pair_id)) state$tuning$models else state$tuning$models[[part$pair_id]]
      unlist(lapply(names(part$postprocessing), function(path) {
        post <- part$postprocessing[[path]]
        inner <- if (is.null(post$inner_folds)) list(Holdout1 = post) else post$inner_folds
        lapply(names(inner), function(j) list(pair = part$rows, plan = inner[[j]], search = local[[path]][[j]]))
      }), recursive = FALSE)
    }), recursive = FALSE)
  }
  for (name in names(specs)) {
    d <- specs[[name]]$data
    spec <- specs[[name]]$spec
    regression <- is.numeric(d$truth)
    metrics <- if (regression) yardstick::metric_set(yardstick::rmse) else yardstick::metric_set(yardstick::mn_log_loss)
    # Repeated outer folds: exact (id, id2, .row, .config) coverage, no notes.
    outer <- withr::with_seed(11, rsample::vfold_cv(d, v = 2, repeats = 2))
    result <- tune::fit_resamples(spec, outer, metrics = metrics,
      control = tune::control_resamples(save_pred = TRUE, extract = workflows::extract_fit_engine, allow_par = FALSE))
    states <- checker$check_outer_identities(result, outer, metrics = 1L)$extracts
    for (h in seq_along(states)) {
      state <- states[[h]]
      F_rows <- as.integer(outer$splits[[h]], "analysis"); A_rows <- as.integer(outer$splits[[h]], "assessment")
      # Final rows after an optional option-selection share S.
      final <- if (is.null(state$tuning$calibration)) F_rows else F_rows[state$tuning$calibration$final_rows]
      S_rows <- if (is.null(state$tuning$calibration)) integer() else F_rows[state$tuning$calibration$selection_rows]
      for (local in searches(state)) {
        rows <- if (is.null(local$pair)) final else final[local$pair]
        B <- rows[local$plan$B]
        search <- local$search
        collected <- tune::collect_predictions(search$native)
        assert_conditions(identical(search$B, local$plan$B), identical(local$plan$selected, search$selected),
          !any(B %in% c(A_rows, S_rows, rows[local$plan$P])),
          identical(sort(unique(collected$.row)), seq_along(B)),
          search$chosen %in% search$candidates$id[search$candidates$status == "eligible"])
      }
    }
    n <- nrow(d); A <- seq.int(.8 * n + 1, n)
    split <- rsample::make_splits(list(analysis = seq_len(.8 * n), assessment = A), d)
    last <- tune::last_fit(spec, split, metrics = metrics)
    fit <- tune::extract_workflow(last)
    new <- d[A, setdiff(names(d), "truth")]
    changed <- d
    changed$truth[A] <- if (regression) changed$truth[A] + 100 else rev(changed$truth[A])
    other <- tune::last_fit(spec, rsample::make_splits(list(analysis = seq_len(.8 * n), assessment = A), changed), metrics = metrics)
    stopifnot(identical(venn_abers_details(fit, new), venn_abers_details(tune::extract_workflow(other), new)))
    tuning <- venn_abers_tuning(fit)
    cat(name, ": searches", length(searches(workflows::extract_fit_engine(fit))),
      "; observed base fits", tuning$cost$observed_base_fits, "\n")
    cases[[name]] <- list(fit = fit, new = new)
  }
  stopifnot(length(cases) == 5L)
}
