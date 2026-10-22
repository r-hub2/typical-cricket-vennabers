# Run with an installed vennabers and its optional tidymodels dependencies.
# Each inner B/P pair learns a complete workflow and postprocessor. Equal
# pair means train VA on C; outer assessment A stays outside all learning.
if (all(vapply(c("vennabers", "tune", "yardstick", "recipes", "tailor", "probably"),
               requireNamespace, logical(1), quietly = TRUE))) local({
  d <- data.frame(x = sin(1:240), row_id = 1:240,
                  truth = factor(rep(c("yes", "no", "no", "yes"), 60)))
  base <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(truth ~ x, d) |> recipes::step_normalize(x)) |>
    workflows::add_model(parsnip::logistic_reg(engine = "glm"))
  # P folds partition each local Q exactly once; each B is the complement.
  # J=3 complete pairs are retained per VA fold. No inner VA is fitted.
  outer <- withr::with_seed(74, rsample::vfold_cv(d, v = 2, repeats = 2, strata = truth, pool = 0))
  metrics <- yardstick::metric_set(yardstick::brier_class, yardstick::mn_log_loss)
  audit <- function(result, plans, event) {
    predictions <- tune::collect_predictions(result, summarize = FALSE)
    scores <- tune::collect_metrics(result, summarize = FALSE)
    stopifnot(nrow(tune::collect_notes(result)) == 0L, all(is.finite(scores$.estimate)),
              nrow(scores) == 2L * nrow(plans))
    ids <- setdiff(names(plans), "splits")
    if (inherits(result, "last_fit")) for (id in ids) scores[[id]] <- plans[[id]][1L]
    expected <- do.call(rbind, lapply(seq_len(nrow(plans)), function(i) {
      rows <- as.integer(plans$splits[[i]], "assessment")
      cbind(as.data.frame(plans[rep(i, length(rows)), ids, drop = FALSE]),
            .row = rows, .config = "pre0_mod0_post0")
    }))
    keys <- c(ids, ".row", ".config")
    stopifnot(!anyDuplicated(predictions[keys]), nrow(predictions) == nrow(expected),
              nrow(merge(predictions[keys], expected, by = keys)) == nrow(expected),
              !anyDuplicated(scores[c(ids, ".metric", ".config")]))
    for (i in seq_len(nrow(plans))) {
      matches <- function(x) Reduce(`&`, lapply(ids, function(id) x[[id]] == plans[[id]][i]))
      p <- predictions[matches(predictions), ]
      m <- scores[matches(scores), ]
      column <- paste0(".pred_", levels(d$truth)[if (event == "first") 1L else 2L])
      values <- c(brier_class = yardstick::brier_class_vec(p$truth, p[[column]], event_level = event),
                  mn_log_loss = yardstick::mn_log_loss_vec(p$truth, p[[column]], event_level = event))
      stopifnot(all(abs(m$.estimate - values[m$.metric]) <= 1e-10))
    }
  }
  for (kind in c("tailor_logistic", "probability", "glm_link")) for (method in c("inductive", "cross")) for (event in c("first", "second")) {
    post <- list(method = if (kind == "tailor_logistic") kind else "sigmoid", strategy = "ensemble", v = 3L)
    if (kind != "tailor_logistic") post$sigmoid_input <- kind
    spec <- vennabers::calibrated_workflow(base, "truth", method = method, event_level = event,
      seed = 42, postprocessing = post)
    fitted <- generics::fit(spec, d)
    details <- vennabers::venn_abers_details(fitted, d[c("x", "row_id")])
    stopifnot(identical(details$prediction, predict(fitted, d, type = "prob")),
              all(details$details$postprocessing_method == post$method),
              all(lengths(details$details$inner_fold) == 3L))
    result <- tune::fit_resamples(spec, outer, metrics = metrics,
      control = tune::control_resamples(save_pred = TRUE, extract = workflows::extract_fit_engine,
        event_level = event, allow_par = FALSE))
    audit(result, outer, event)
    states <- tune::collect_extracts(result)$.extracts
    for (i in seq_along(states)) {
      local <- rsample::analysis(outer$splits[[i]])
      A <- rsample::assessment(outer$splits[[i]])$row_id
      for (ensemble in states[[i]]$postprocessing) for (plan in ensemble$inner_folds) {
        stopifnot(!any(local$row_id[c(plan$B, plan$P, plan$C)] %in% A),
                  !length(intersect(plan$B, plan$P)), !length(intersect(plan$Q, plan$C)))
      }
    }
    split <- rsample::make_splits(list(analysis = 1:180, assessment = 181:240), d)
    last <- tune::last_fit(spec, split, metrics = metrics,
      control = tune::control_last_fit(event_level = event))
    audit(last, rsample::manual_rset(list(split), "train/test split"), event)
    # Changing either A predictors or truth changes evaluation, never learning.
    for (variable in c("x", "truth")) {
      changed <- d
      if (variable == "x") changed$x[181:240] <- 1000 else changed$truth[181:240] <- "yes"
      changed_split <- rsample::make_splits(list(analysis = 1:180, assessment = 181:240), changed)
      other <- tune::last_fit(spec, changed_split, metrics = metrics,
        control = tune::control_last_fit(event_level = event))
      original <- workflows::extract_fit_engine(last$.workflow[[1]])
      after <- workflows::extract_fit_engine(other$.workflow[[1]])
      stopifnot(identical(original$postprocessing, after$postprocessing),
        identical(predict(last$.workflow[[1]], d, type = "prob"), predict(other$.workflow[[1]], d, type = "prob")))
    }
  }
})
