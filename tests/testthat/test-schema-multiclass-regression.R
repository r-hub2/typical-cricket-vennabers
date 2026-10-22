testthat::test_that("schema-multiclass-regression", {
library(vennabers)
source(system.file("examples", "assertions.R", package = "vennabers"))
source("helper-helpers.R")
if (all(vapply(c("parsnip", "workflows", "rsample", "recipes", "withr"),
               requireNamespace, logical(1), quietly = TRUE))) {
  error <- function(expr, pattern = "") {
    e <- tryCatch(force(expr), error = identity)
    stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e)))
  }
  # Public reference engine records the actual training data and distinguishes folds.
  fits <- 0L
  encoded_fit <- function(x, y) {
    fits <<- fits + 1L
    list(x = x, y = y, levels = levels(y))
  }
  encoded_prob <- function(object, new_data) {
    k <- length(object$levels)
    p <- as.matrix(new_data[seq_len(k)])
    p <- p / rowSums(p)
    colnames(p) <- paste0(".pred_", object$levels)
    as.data.frame(p)[rev(seq_len(k))]
  }
  encoded_numeric <- function(object, new_data) mean(object$y) + sd(object$y) * new_data[[1]]
  parsnip::set_new_model("encoded_reference")
  for (mode in c("classification", "regression")) {
    register_reference_engine("encoded_reference", mode, "encoded_fit",
      if (mode == "regression") "encoded_numeric" else "encoded_prob",
      if (mode == "regression") "numeric" else "prob", new_model = FALSE)
  }
  make_base <- function(data, outcome, mode) workflows::workflow() |>
    workflows::add_formula(reformulate(setdiff(names(data), outcome), outcome)) |>
    workflows::add_model(parsnip::new_model_spec("encoded_reference", args = list(), eng_args = NULL,
                                                mode = mode, method = NULL, engine = "reference"))
  # Fixed before comparison: s=1 outcome unit, multiplied when units change.
  close <- function(actual, expected, scale = NULL) {
    assert_conditions(identical(dim(actual), dim(expected)), length(actual) == length(expected),
              all(is.finite(actual)), all(is.finite(expected)),
              all(abs(actual - expected) <= if (is.null(scale)) 1e-10 else 1e-10 * scale + 1e-10 * abs(expected)))
  }
  class_names <- c("class", ".pred_yes", "a b", "a-b", "a.b", "1yes", "ümlaut", "vaL0001")
  for (k in 2:4) for (classification in c("ovo", "ovr")) for (method in c("inductive", "cross")) {
    original <- as.data.frame(matrix(abs(sin(seq_len(120 * k))) + .01, 120, k))
    original$row <- 1:120
    classes <- class_names[seq_len(k) + if (classification == "ovr") 4L else 0L]
    original$y <- factor(rep(classes, length.out = 120), levels = classes)
    names(original)[c(1, k + 1L, k + 2L)] <- c("vaV0001", "row identity", ".pred_class")
    schema <- venn_abers_schema(original, ".pred_class", classification = classification)
    data <- venn_abers_encode(original, schema)
    outcome <- schema$encoded_outcome
    base <- make_base(data, outcome, "classification")
    for (precision in list(NULL, 0, 2, -1)) {
      before <- fits
      fit <- generics::fit(calibrated_workflow(base, outcome, schema = schema, method = method,
                                               seed = 42, precision = precision), data)
      state <- workflows::extract_fit_engine(fit)
      assert_conditions(fits == before + length(state$folds), identical(names(state$models), names(state$folds)),
                identical(state$class_levels, unname(schema$classes)), is.null(state$event_class),
                !state$parameters$explicit[["classification"]], !state$parameters$explicit[["event_level"]])
      new <- data[c(1, 19, 87), setdiff(names(data), outcome)]
      ps <- ys <- fs <- new_p <- list()
      for (id in names(state$folds)) {
        f <- state$folds[[id]]
        engine <- workflows::extract_fit_engine(state$models[[id]])
        stopifnot(identical(as.numeric(engine$x[[k + 1L]]), as.numeric(data[[k + 1L]][f$analysis])),
                  identical(as.vector(engine$y), as.vector(data[[outcome]][f$analysis])))
        p <- predict(state$models[[id]], data[f$assessment, ], type = "prob")
        ps[[id]] <- as.matrix(p[paste0(".pred_", state$class_levels)])
        colnames(ps[[id]]) <- classes
        ys[[id]] <- original[[schema$outcome]][f$assessment]
        fs[[id]] <- rep(id, length(f$assessment))
        p <- predict(state$models[[id]], new, type = "prob")
        new_p[[id]] <- as.matrix(p[paste0(".pred_", state$class_levels)])
        colnames(new_p[[id]]) <- classes
      }
      core <- if (method == "cross") cross_venn_abers(do.call(c, ys), do.call(rbind, ps), unlist(fs),
        precision = precision, method = classification) else
        venn_abers(ys[[1]], ps[[1]], precision = precision, method = classification)
      for (loss in c("log", "brier")) {
        expected <- predict(core, if (method == "cross") new_p else new_p[[1]], loss = loss, type = "details")
        native <- venn_abers_details(fit, new, loss = loss)
        decoded <- venn_abers_decode(native, schema)
        close(decoded$probability, expected$prediction)
        diagnostics <- if (classification == "ovo") expected$pairs else expected$classes
        close(as.matrix(decoded$details[c("lower", "upper", "probability")]),
              as.matrix(diagnostics[c("lower", "upper", "probability")]))
        stopifnot(identical(decoded$details$row, diagnostics$row),
                  identical(decoded$details$event_class, diagnostics$event),
                  identical(decoded$probability, venn_abers_decode(predict(fit, new, type = "prob", loss = loss), schema)$probability))
        if (classification == "ovo") stopifnot(identical(decoded$details$class_1, diagnostics$non_event),
                                                 identical(decoded$details$class_2, diagnostics$event)) else
          stopifnot(identical(decoded$details$class, diagnostics$event))
        if (method == "cross") stopifnot(identical(decoded$details$fold, diagnostics$fold))
        cl <- venn_abers_decode(predict(fit, new, type = "class", loss = loss), schema)$class
        stopifnot(identical(as.character(cl), classes[max.col(expected$prediction, ties.method = "first")]))
        for (n in 0:1) {
          small <- new[seq_len(n), ]
          d <- venn_abers_decode(venn_abers_details(fit, small, loss = loss), schema)
          cl <- venn_abers_decode(predict(fit, small, type = "class", loss = loss), schema)$class
          assert_conditions(identical(dim(d$probability), c(n, k)), is.double(d$probability),
                    identical(levels(cl), classes), length(cl) == n,
                    identical(colnames(d$probability), classes))
        }
      }
      plans <- lapply(state$folds, function(f) rsample::make_splits(f[c("analysis", "assessment")], data))
      plan <- if (method == "cross") rsample::manual_rset(plans, names(state$folds)) else plans[[1]]
      explicit <- generics::fit(calibrated_workflow(base, outcome, schema = schema, method = method,
                                                   splits = plan, precision = precision), data)
      stopifnot(identical(venn_abers_details(explicit, new), venn_abers_details(fit, new)))
    }
    # Early schema/universe failures must precede any model fit.
    before <- fits
    error(calibrated_workflow(base, outcome, schema = schema, event_level = "first"))
    error(calibrated_workflow(base, outcome, schema = schema, classification = "binary"))
    bad <- data; bad[[outcome]] <- factor(bad[[outcome]], levels = rev(state$class_levels))
    model <- calibrated_workflow(base, outcome, schema = schema, method = method, seed = 42)
    error(generics::fit(model, bad), "schema")
    bad <- data; bad[[outcome]][1] <- NA
    error(generics::fit(model, bad), "nonmissing")
    bad <- data[data[[outcome]] != state$class_levels[k], ]
    error(generics::fit(model, bad), "all classes")
    foreign <- venn_abers_schema(original, schema$outcome, classification = if (classification == "ovo") "ovr" else "ovo")
    bad <- data; attr(bad, "vennabers_schema") <- foreign
    error(generics::fit(model, bad), "different schema")
    error(predict(fit, bad, type = "prob"), "different schema")
    error(venn_abers_decode(venn_abers_details(fit, new), foreign), "different schema")
    stopifnot(fits == before)
    if (k == 2L) {
      tied <- data; tied[1:2] <- .5
      fit <- generics::fit(model, tied)
      for (loss in c("log", "brier")) {
        stopifnot(all(predict(fit, tied, type = "prob", loss = loss) == .5),
                  all(venn_abers_decode(predict(fit, tied, type = "class", loss = loss), schema)$class == classes[1]))
      }
    }
  }
}
if (exists("encoded_fit")) {
  outcomes <- c(".pred", ".pred_class", ".row", ".config", "id", "id2", ".pred_vaL0001", "a b")
  for (i in seq_along(outcomes)) for (method in c("inductive", "cross")) {
    original <- data.frame(x = sin(1:60), row = 1:60, y = 3 + sin(1:60) + cos(1:60))
    names(original) <- c("vaV0001", "row identity", outcomes[i])
    for (scale in if (i == 1L) c(1e-9, 1, 1e9) else 1) {
      raw <- original; raw[[3]] <- scale * raw[[3]]
      schema <- venn_abers_schema(raw, outcomes[i])
      data <- venn_abers_encode(raw, schema)
      outcome <- schema$encoded_outcome
      base <- make_base(data, outcome, "regression")
      # Unequal calibration folds expose accidental size weighting/model swaps.
      held <- list(1:7, 8:23, 24:60)
      plans <- lapply(held, function(h) rsample::make_splits(
        list(analysis = setdiff(1:60, h), assessment = h), data))
      plan <- if (method == "cross") rsample::manual_rset(plans, c("a", "b", "c")) else plans[[2]]
      for (option in list(list(m = 2), list(epsilon = .5))) {
        before <- fits
        model <- do.call(calibrated_workflow, c(list(workflow = base, outcome = outcome,
          schema = schema, method = method, splits = plan), option))
        fit <- generics::fit(model, data)
        state <- workflows::extract_fit_engine(fit)
        stopifnot(fits == before + length(state$folds), identical(state$schema, schema),
                  is.null(state$class_levels), is.null(state$event_class))
        ps <- ys <- fs <- new_p <- list()
        new <- data[c(1, 19, 57), 1:2]
        for (id in names(state$folds)) {
          f <- state$folds[[id]]
          engine <- workflows::extract_fit_engine(state$models[[id]])
          stopifnot(identical(as.numeric(engine$x[[2]]), as.numeric(data[[2]][f$analysis])),
                    identical(as.vector(engine$y), data[[outcome]][f$analysis]))
          ps[[id]] <- predict(state$models[[id]], data[f$assessment, ])$.pred
          ys[[id]] <- data[[outcome]][f$assessment]
          fs[[id]] <- rep(id, length(f$assessment))
          new_p[[id]] <- predict(state$models[[id]], new)$.pred
        }
        core <- if (method == "cross") do.call(cross_venn_abers_regression,
          c(list(truth = unlist(ys), prediction = unlist(ps), fold = unlist(fs)), option)) else
          do.call(venn_abers_regression, c(list(truth = ys[[1]], prediction = ps[[1]]), option))
        expected <- if (method == "cross") predict(core, do.call(cbind, new_p), type = "details") else
          predict(core, new_p[[1]])
        native <- venn_abers_details(fit, new)
        decoded <- venn_abers_decode(native, schema)
        points <- if (method == "cross") expected$prediction else expected
        close(decoded$prediction, points$prediction, scale)
        stopifnot(identical(decoded$details, native$details), is.null(decoded$probability), is.null(decoded$class),
                  identical(decoded$prediction, venn_abers_decode(predict(fit, new), schema)$prediction))
        if (method == "cross") {
          stopifnot(identical(decoded$aggregate, native$aggregate))
          close(as.matrix(decoded$aggregate[c("lower", "upper")]), as.matrix(points[c("lower", "upper")]), scale)
          for (id in names(state$folds)) close(as.matrix(decoded$details[decoded$details$fold == id,
            c("prediction", "lower", "upper")]), as.matrix(expected$folds[[id]]), scale)
          close(rowMeans(matrix(decoded$details$prediction, nrow(new))), decoded$prediction, scale)
          if (i == 1L && scale == 1 && !is.null(option$m)) {
            wrong <- do.call(cbind, rev(new_p)); colnames(wrong) <- names(new_p)
            stopifnot(any(abs(predict(core, wrong)$prediction - points$prediction) > 1e-5),
              any(abs(as.vector(matrix(decoded$details$prediction, nrow(new)) %*% lengths(held) / 60) - points$prediction) > 1e-5))
          }
        } else close(as.matrix(decoded$details[c("prediction", "lower", "upper")]), as.matrix(points), scale)
        for (n in 0:1) {
          d <- venn_abers_decode(venn_abers_details(fit, new[seq_len(n), ]), schema)
          assert_conditions(is.double(d$prediction), length(d$prediction) == n,
                    nrow(d$details) == n * length(state$folds), is.integer(d$details$row),
                    is.character(d$details$fold), all(vapply(d$details[3:5], is.double, logical(1))))
          if (method == "cross") stopifnot(nrow(d$aggregate) == n)
        }
      }
      before <- fits
      for (option in list(list(classification = "binary"), list(event_level = "first"), list(loss = "log"),
                          list(precision = NULL))) error(do.call(calibrated_workflow,
        c(list(workflow = base, outcome = outcome, schema = schema), option)))
      bad <- data; bad[[outcome]][1] <- Inf
      error(generics::fit(model, bad), "finite")
      bad <- data; bad[[outcome]] <- factor(bad[[outcome]])
      error(generics::fit(model, bad), "numeric")
      bad <- data; attr(bad, "vennabers_schema") <- NULL; bad[[1]][1] <- 999
      error(generics::fit(model, bad), "data")
      stopifnot(fits == before)
    }
    # Automatic and prepared local plans agree, also for numeric 0/1 truth.
    raw <- original; raw[[3]] <- rep(c(0L, 1L), 30)
    schema <- venn_abers_schema(raw, outcomes[i]); data <- venn_abers_encode(raw, schema)
    base <- make_base(data, schema$encoded_outcome, "regression")
    fit <- generics::fit(calibrated_workflow(base, schema$encoded_outcome, schema = schema,
                                            method = method, seed = 42), data)
    state <- workflows::extract_fit_engine(fit)
    plans <- lapply(state$folds, function(f) rsample::make_splits(f[c("analysis", "assessment")], data))
    plan <- if (method == "cross") rsample::manual_rset(plans, names(state$folds)) else plans[[1]]
    explicit <- generics::fit(calibrated_workflow(base, schema$encoded_outcome, schema = schema,
      method = method, splits = local_calibration_splits(list(plan), names(data))), data)
    stopifnot(identical(venn_abers_details(fit, data), venn_abers_details(explicit, data)))
  }
}
if (exists("encoded_fit")) {
  for (classification in c("ovo", "ovr")) {
    raw <- data.frame(p = rep(.2, 60), q = rep(.3, 60), r = rep(.5, 60),
                      y = factor(rep(c("class", ".pred_yes", "vaL0001"), 20)))
    schema <- venn_abers_schema(raw, "y", classification = classification)
    data <- venn_abers_encode(raw, schema); outcome <- schema$encoded_outcome
    base <- make_base(data, outcome, "classification")
    make_plan <- function(held) rsample::manual_rset(lapply(held, function(h)
      rsample::make_splits(list(analysis = setdiff(1:60, h), assessment = h), data)), paste0("f", seq_along(held)))
    # A late singleton is permitted for OvR; OvO rejects its empty pair before fitting.
    plan <- make_plan(list(2:30, 31:60, 1L))
    model <- calibrated_workflow(base, outcome, schema = schema, method = "cross", splits = plan)
    before <- fits
    if (classification == "ovo") {
      error(generics::fit(model, data), "Empty calibration for pair")
      stopifnot(fits == before)
    } else {
      fit <- generics::fit(model, data)
      stopifnot(identical(workflows::extract_fit_engine(fit)$class_levels, unname(schema$classes)),
                fits == before + 3L)
    }
    absent <- which(data[[outcome]] != schema$classes[3])
    plan <- rsample::make_splits(list(analysis = setdiff(1:60, absent), assessment = absent), data)
    before <- fits
    error(generics::fit(calibrated_workflow(base, outcome, schema = schema, splits = plan), data), "missing")
    stopifnot(fits == before)
  }
  raw <- data.frame(x = sin(1:60), y = 1:60)
  schema <- venn_abers_schema(raw, "y"); data <- venn_abers_encode(raw, schema)
  base <- make_base(data, schema$encoded_outcome, "regression")
  plans <- lapply(list(1:58, 59:60), function(h) rsample::make_splits(
    list(analysis = setdiff(1:60, h), assessment = h), data))
  before <- fits
  error(generics::fit(calibrated_workflow(base, schema$encoded_outcome, schema = schema, method = "cross",
    splits = rsample::manual_rset(plans, c("valid", "invalid"))), data), "invalid calibration size 2")
  stopifnot(fits == before)
}

  testthat::expect_true(TRUE)
})
