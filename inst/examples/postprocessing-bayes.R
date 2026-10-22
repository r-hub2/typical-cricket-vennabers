# Stepwise Bayes selection with learned postprocessing. As with Grid, every
# engine/recipe search runs inside its own B (each VA source, inner ensemble
# pair and full F composition, per class pair for pairwise OvO); P, C, the outer
# assessment A and a reserved option-selection share S never enter it. Only
# complete valid candidates feed the next GP. Failed candidates stay in the
# history with their reasons and spend budget; too few valid initial points stop
# the fit with a clear error, never a switch to Grid.
source(system.file("examples", "assertions.R", package = "vennabers"))
packages <- c("vennabers", "parsnip", "rlang", "withr", "workflows", "recipes", "rsample", "tune", "dials",
  "yardstick", "tailor", "probably", "nnet")
cases <- list()
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  library(vennabers)
  checker <- new.env()
  sys.source(system.file("examples", "check-encoded-evaluation.R", package = "vennabers"), checker)
  stratified <- function(d) rsample::vfold_cv(d, v = 2, strata = truth, pool = 0)
  plain <- function(d) rsample::vfold_cv(d, v = 2)
  bayes <- function(upper, resamples, initial = 3L, iter = 3L) list(method = "bayes", initial = initial, iter = iter,
    param_info = dials::parameters(dials::new_quant_param(type = "integer", range = c(1L, upper),
      inclusive = c(TRUE, TRUE), label = c(degree = "Degree"))),
    resamples = resamples, control = tune::control_bayes(seed = 4, save_pred = TRUE, no_improve = 3L))
  classes <- withr::with_seed(7, {
    truth <- factor(rep(c("a", "b", "c"), 800))
    data.frame(x = stats::rnorm(2400, sd = 3) + c(-1, 1, -1)[truth],
      z = stats::rnorm(2400, sd = 3) + c(-1, -1, 1)[truth], truth = truth)
  })
  # z has four distinct values: step_poly(z) fails for degree >= 4 in every B.
  binary <- withr::with_seed(8, {
    x <- stats::rnorm(1200, sd = 1.5)
    z <- sample(1:4, 1200, TRUE)
    data.frame(x = x, z = z, truth = factor(ifelse(stats::runif(1200) < stats::plogis(x - .3 * z + .5), "yes", "no")))
  })
  numeric <- withr::with_seed(9, {
    x <- stats::rnorm(1200)
    data.frame(x = x, z = stats::rnorm(1200), truth = 3 + x - .5 * x^2 + stats::rnorm(1200, sd = .5))
  })
  base <- function(d, engine, term = "x") workflows::workflow() |> workflows::add_recipe(
    recipes::recipe(truth ~ x + z, d) |> recipes::step_poly(!!rlang::sym(term), degree = tune::tune("degree"))) |>
    workflows::add_model(engine)
  glm <- parsnip::logistic_reg(engine = "glm")
  nnet <- parsnip::multinom_reg(penalty = .1) |> parsnip::set_engine("nnet", trace = FALSE, maxit = 1000)
  sigmoid <- list(method = "sigmoid", sigmoid_input = "probability", strategy = "holdout", prop = .3)
  specs <- list(
    # Budget and failure case: degrees 4 to 6 fail in every B, initial or proposed.
    binary_sigmoid_joint = list(data = binary, spec = calibrated_workflow(base(binary, glm, "z"), "truth",
      prediction_model = "full", seed = 3, tuning = list(strategy = "joint", model = bayes(6L, stratified)),
      postprocessing = sigmoid)),
    pairwise_tailor_two_stage = list(data = classes[c("x", "z", "truth")], spec = calibrated_workflow(base(classes, glm), "truth",
      classification = "ovo", model_training = "pairwise", seed = 4,
      tuning = list(strategy = "two_stage", model = bayes(2L, stratified, initial = 2L)),
      postprocessing = list(method = "tailor_logistic", strategy = "ensemble", v = 2))),
    ovr_multinomial_options = list(data = classes, spec = calibrated_workflow(base(classes, nnet), "truth",
      classification = "ovr", seed = 5, loss = tune::tune("va_loss"),
      tuning = list(strategy = "two_stage", selection_prop = .5,
        model = bayes(2L, stratified, initial = 2L),
        calibration = list(method = "bayes", initial = 2L, iter = 1L, resamples = stratified,
          param_info = dials::parameters(dials::new_qual_param(type = "character", values = c("log", "brier"),
            label = c(va_loss = "Loss"))), control = tune::control_bayes(seed = 5, save_pred = TRUE))),
      postprocessing = list(method = "tailor_multinomial", strategy = "holdout", prop = .3))),
    regression_linear_joint = list(data = numeric, spec = calibrated_workflow(base(numeric, parsnip::linear_reg()), "truth",
      method = "cross", v = 2, seed = 6, tuning = list(strategy = "joint", model = bayes(3L, plain)),
      postprocessing = list(method = "tailor_linear", strategy = "holdout", prop = .3))))
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
        eligible <- search$candidates$status == "eligible"
        # The GP only saw complete valid candidates, all on rows of this B.
        assert_conditions(identical(search$B, local$plan$B), identical(local$plan$selected, search$selected),
          !any(B %in% c(A_rows, S_rows, rows[local$plan$P])),
          identical(sort(unique(tune::collect_predictions(search$native)$.row)), seq_along(B)),
          search$chosen %in% search$candidates$id[eligible], !anyDuplicated(search$candidates$id),
          setequal(unique(tune::collect_metrics(search$native)$degree), search$candidates$degree[eligible]),
          search$cost$candidate_evaluations == search$attempts * length(search$resamples))
        if (name == "binary_sigmoid_joint") {
          # Failed initial and continuation candidates keep their IDs and reasons.
          assert_conditions(any(!eligible), identical(!eligible, search$candidates$degree >= 4L),
            grepl("must be less than number of unique points", paste(search$notes$note, collapse = " ")))
        }
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
    local <- searches(workflows::extract_fit_engine(fit))
    cat(name, ": searches", length(local),
      "; attempts", paste(vapply(local, function(x) x$search$attempts, 1L), collapse = "/"),
      "; stop", paste(unique(vapply(local, function(x) x$search$stop_reason, "")), collapse = "/"),
      "; excluded", sum(vapply(local, function(x) x$search$cost$excluded_candidates, 1L)),
      "; observed base fits", venn_abers_tuning(fit)$cost$observed_base_fits, "\n")
    cases[[name]] <- list(fit = fit, new = new)
  }
  # Two initial points (2 and 5) leave one valid candidate: a clear error, no Grid fallback.
  short <- calibrated_workflow(base(binary, glm, "z"), "truth", seed = 3,
    tuning = list(strategy = "joint", model = bayes(6L, stratified, initial = 2L)), postprocessing = sigmoid)
  failure <- tryCatch({
    generics::fit(short, binary)
    ""
  }, error = conditionMessage)
  assert_conditions(grepl("Bayes needs at least 2 complete initial model candidates; found 1", failure),
    grepl("must be less than number of unique points", failure))
  stopifnot(length(cases) == 4L)
}
