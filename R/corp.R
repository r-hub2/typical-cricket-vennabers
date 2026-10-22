#' @export
corp_plot <- function(truth, probability, band = "consistency", level = 0.9,
                      replicates = 1000, seed = NULL, ...) {
  probability <- .validate_corp_input(truth, probability)
  .validate_choice(band, c("consistency", "confidence", "none"), "band")
  if (!.is_proportion(level)) {
    stop("`level` must be one number strictly between 0 and 1.", call. = FALSE)
  }
  if (!.is_count(replicates)) {
    stop("`replicates` must be one finite whole number of at least 1.", call. = FALSE)
  }
  # Each band quantile needs a replicate beyond it. The tolerance keeps exact
  # boundaries such as 20 replicates at level 0.9, where 1 - 0.9 < 0.1.
  if (replicates * (1 - level) / 2 < 1 - 1e-10) {
    stop("`replicates * (1 - level) / 2` must be at least 1; increase `replicates` or lower `level`.",
         call. = FALSE)
  }
  if (!is.null(seed) && !.is_count(seed, -.Machine$integer.max)) {
    stop("`seed` must be NULL or one integer in R's seed range.", call. = FALSE)
  }
  controls <- list(band = band, level = level, replicates = as.integer(replicates),
                   seed = if (is.null(seed)) NA_integer_ else as.integer(seed))
  draw <- function() {
    .corp_panels(truth, probability, "CORP reliability", function(truth, probability, main, class) {
      table <- .corp_table(truth, probability, controls, class)
      .corp_draw(table, utils::modifyList(list(main = main), list(...)))
      table
    })
  }
  # One seed covers all classes, which draw their replicates in turn.
  invisible(if (is.null(seed)) draw() else .with_seed(seed, draw()))
}

#' @export
corp_decomposition_plot <- function(truth, probability, ...) {
  probability <- .validate_corp_input(truth, probability)
  if (is.factor(truth)) {
    classes <- levels(truth)
    rows <- lapply(classes, function(class) .corp_decomposition(truth == class, probability[, class]))
    table <- data.frame(class = factor(classes, levels = classes), do.call(rbind, rows))
    height <- as.matrix(table[c("mcb", "dsc", "unc")])
    # The class legend takes rows of up to five classes in a strip above the
    # tallest bar.
    columns <- min(length(classes), 5L)
    defaults <- list(beside = TRUE, legend.text = classes,
                     col = grDevices::hcl.colors(length(classes), "Dark 3"),
                     args.legend = list(x = "top", ncol = columns, bty = "n"))
    if (max(height) > 0) {
      defaults$ylim <- c(0, (1 + .12 * ceiling(length(classes) / columns)) * max(height))
    }
  } else {
    table <- .corp_decomposition(truth, probability)
    height <- unlist(table[c("mcb", "dsc", "unc")])
    defaults <- list()
  }
  # Groups are the components; multiclass bars within a group are the classes.
  options <- utils::modifyList(c(defaults, list(
    names.arg = c("MCB", "DSC", "UNC"),
    ylab = "Brier score component",
    main = "CORP Brier score decomposition"
  )), list(...))
  options$height <- unname(height)
  do.call(graphics::barplot, options)
  invisible(table)
}

# Truth and probability under the scoring contract; multiclass returns the
# normalized class matrix.
.validate_corp_input <- function(truth, probability) {
  if (is.factor(truth)) return(.validate_multiclass_input(truth, probability))
  .validate_binary_input(truth, probability)
  probability
}

# One table per panel: binary draws one panel, multiclass one per class
# against the rest in its own grid, each table led by `class`. The panel
# receives the class, NULL for binary truth.
.corp_panels <- function(truth, probability, title, panel) {
  if (!is.factor(truth)) return(invisible(panel(truth, probability, title, NULL)))
  classes <- levels(truth)
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old))
  rows <- ceiling(sqrt(length(classes)))
  graphics::par(mfrow = c(rows, ceiling(length(classes) / rows)))
  tables <- lapply(classes, function(class) {
    table <- panel(truth == class, probability[, class], paste0("CORP: ", class), class)
    data.frame(class = factor(rep(class, nrow(table)), levels = classes), table)
  })
  invisible(do.call(rbind, tables))
}

# The CORP estimate: the PAV fit of the events on the probabilities, one row
# per unique probability in increasing order. stats::isoreg orders tied
# probabilities by decreasing event, a run that PAV pools into one block, so
# tied probabilities share one estimate.
.corp_pav <- function(truth, probability) {
  fit <- stats::isoreg(probability, as.numeric(truth))
  sorted <- if (fit$isOrd) probability else probability[fit$ord]
  first <- !duplicated(sorted)
  data.frame(probability = sorted[first], estimate = fit$yf[first],
             count = tabulate(cumsum(first)))
}

.corp_table <- function(truth, probability, controls, class = NULL) {
  table <- .corp_pav(truth, probability)
  n <- nrow(table)
  band <- if (controls$band != "none") .corp_band(truth, probability, table, controls, class)
  if (is.null(band)) {
    band <- list(lower = rep(NA_real_, n), upper = rep(NA_real_, n), valid = rep(NA_integer_, n))
  }
  data.frame(table, band_lower = band$lower, band_upper = band$upper,
             valid_replicates = band$valid, lapply(controls, rep, n))
}

# The pointwise band by resampling as in reliabilitydiag 0.2.1. Each
# replicate draws the pairs with replacement and their events from the
# diagonal (consistency) or the CORP estimate (confidence), refits PAV and
# interpolates it linearly at the unique probabilities, missing outside the
# drawn range. Pairs in the reference order and the same calls to sample()
# and rbinom() give the reference band for equal seeds. NULL when the
# confidence band is omitted.
.corp_band <- function(truth, probability, table, controls, class) {
  truth <- as.numeric(truth)
  if (controls$band == "confidence" && length(unique(truth)) == 1L) {
    warning("Confidence band omitted: `truth` has only one outcome",
            if (!is.null(class)) paste0(" for class \"", class, "\""), ".", call. = FALSE)
    return(NULL)
  }
  x <- probability[order(probability, -truth)]
  center <- if (controls$band == "consistency") x else table$estimate[match(x, table$probability)]
  values <- table$probability
  n <- length(x)
  draws <- vapply(seq_len(controls$replicates), function(replicate) {
    s <- sample(n, replace = TRUE)
    events <- stats::rbinom(n, 1L, center[s])
    fit <- .corp_pav(events, x[s])
    if (nrow(fit) == 1L) {
      ifelse(values == fit$probability, fit$estimate, NA_real_)
    } else {
      stats::approx(fit$probability, fit$estimate, xout = values)$y
    }
  }, numeric(length(values)))
  draws <- matrix(draws, nrow = length(values))
  limits <- apply(draws, 1L, stats::quantile, probs = .5 + c(-.5, .5) * controls$level,
                  na.rm = TRUE, names = FALSE, type = 7)
  band <- list(lower = limits[1L, ], upper = limits[2L, ], valid = as.integer(rowSums(!is.na(draws))))
  if (controls$band == "confidence") {
    band$lower <- .corp_bound_correction(band$lower, table)
    band$upper <- .corp_bound_correction(band$upper, table)
  }
  band
}

# The heuristic of the reference for the confidence band: inner
# probabilities with an estimate of exactly 0 or 1 take limits interpolated
# linearly from their neighbours.
.corp_bound_correction <- function(limit, table) {
  x <- table$probability
  if (length(x) == 1L) return(limit)
  limit[table$estimate %in% c(0, 1) & !(x %in% range(x))] <- NA
  if (sum(!is.na(limit)) < 2L) return(limit)
  stats::approx(x, limit, xout = x)$y
}

.corp_draw <- function(table, graphical) {
  options <- utils::modifyList(list(
    type = "s",
    xlim = c(0, 1),
    ylim = c(0, 1),
    xaxs = "i",
    yaxs = "i",
    xlab = "Predicted probability",
    ylab = "Conditional event probability"
  ), graphical)
  options$x <- table$probability
  options$y <- table$estimate
  band <- table[!is.na(table$band_lower) & !is.na(table$band_upper), ]
  fill <- grDevices::adjustcolor("grey50", alpha.f = .35)
  # The band goes under the step: panel.first draws it once the plot region
  # is set up.
  if (nrow(band) > 0L && is.null(options$panel.first)) {
    options$panel.first <- as.call(list(quote(graphics::polygon),
                                        c(band$probability, rev(band$probability)),
                                        c(band$band_lower, rev(band$band_upper)),
                                        col = fill, border = NA))
  }
  do.call(graphics::plot, options)
  graphics::abline(a = 0, b = 1, lty = 2)
  if (nrow(band) > 0L) {
    level <- paste0(signif(100 * table$level[1L], 10), "%")
    label <- if (table$band[1L] == "consistency") {
      paste(level, "pointwise consistency band")
    } else {
      paste("Approx.", level, "pointwise confidence band")
    }
    graphics::legend("topleft", label, fill = fill, border = NA, bty = "n")
  }
}

# Evaluates `code` after set.seed(seed) and restores the caller's random
# state on exit, also its absence, as withr::with_seed() does. The RNG kind
# stays as it is.
.with_seed <- function(seed, code) {
  env <- globalenv()
  if (exists(".Random.seed", envir = env, inherits = FALSE)) {
    state <- get(".Random.seed", envir = env, inherits = FALSE)
    on.exit(assign(".Random.seed", state, envir = env))
  } else {
    on.exit(if (exists(".Random.seed", envir = env, inherits = FALSE)) rm(".Random.seed", envir = env))
  }
  set.seed(seed)
  code
}

# The CORP decomposition of the binary Brier score into miscalibration (MCB),
# discrimination (DSC) and uncertainty (UNC), against the PAV estimate on the
# same pairs and the constant event rate.
.corp_decomposition <- function(truth, probability) {
  truth <- as.numeric(truth)
  pav <- .corp_pav(truth, probability)
  estimate <- pav$estimate[match(probability, pav$probability)]
  brier <- mean((truth - probability)^2)
  pooled <- mean((truth - estimate)^2)
  unc <- mean((truth - mean(truth))^2)
  data.frame(brier = brier, mcb = brier - pooled, dsc = unc - pooled, unc = unc)
}
