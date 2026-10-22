# Temporary diagnostic for issue #63; never part of a candidate.
testthat::test_that("ATLAS diagnostic: one row against the same row in a batch", {
  for (pkg in c("workflows", "tailor", "probably", "rsample")) testthat::skip_if_not_installed(pkg)
  ulps <- function(a, b) if (a == b) 0 else (a - b) / 2^(floor(log2(abs(b))) - 52)
  lines <- c(paste("BLAS:", extSoftVersion()[["BLAS"]], "| La:", La_library()))
  for (s in c(1e-9, 1, 1e9)) for (shift in c(0, 100)) for (method in c("inductive", "cross")) {
    d <- data.frame(x = sin(1:120), truth = s * (shift + 3 * sin(1:120) + cos(1:120)))
    base <- workflows::workflow() |> workflows::add_formula(truth ~ x) |>
      workflows::add_model(parsnip::linear_reg(engine = "lm"))
    held_rows <- list(1:15, 16:60, 61:120)
    splits <- lapply(held_rows, function(C) rsample::make_splits(list(analysis = setdiff(1:120, C), assessment = C), d))
    plan <- if (method == "cross") rsample::manual_rset(splits, c("a", "b", "c")) else splits[[1L]]
    for (option in list(list(m = 1), list(m = 7), list(epsilon = .875))) {
      fit <- generics::fit(do.call(calibrated_workflow, c(list(workflow = base, outcome = "truth",
        method = method, splits = plan, seed = 42,
        postprocessing = list(method = "tailor_linear", strategy = "holdout", prop = .3)), option)), d)
      state <- workflows::extract_fit_engine(fit)
      fold <- vapply(names(state$models), function(id) {
        m <- state$models[[id]]; cal <- state$calibration[[id]]
        one <- predict(m, d[1, , drop = FALSE])$.pred; all <- predict(m, d)$.pred[1L]
        held <- cal$prediction[cal$row == 1L]
        sprintf("%s[1v120 %+g ulp, 1vC %s, 120vC %s]", id, ulps(one, all),
          if (length(held)) format(ulps(one, held)) else "-", if (length(held)) format(ulps(all, held)) else "-")
      }, "")
      one <- venn_abers_details(fit, d[1L, "x", drop = FALSE])$prediction$.pred
      all <- venn_abers_details(fit, d["x"])$prediction$.pred[1L]
      tol <- 1e-10 * s + 1e-10 * abs(all)
      q <- data.frame(x = .5)
      qd <- abs(venn_abers_details(fit, q)$prediction$.pred -
        venn_abers_details(fit, rbind(d["x"], q))$prediction$.pred[121L])
      lines <- c(lines, sprintf("s=%g shift=%g %s %s=%g | old |1-120|=%.3g tol=%.2g %s | new off-knot |1-121|=%.3g %s | %s",
        s, shift, method, names(option), option[[1L]], abs(one - all), tol, if (abs(one - all) > tol) "OLD-FAIL" else "old-ok",
        qd, if (qd > tol) "NEW-FAIL" else "new-ok", paste(fold, collapse = " ")))
    }
  }
  testthat::fail(paste(c("", lines), collapse = "\n"))
})
