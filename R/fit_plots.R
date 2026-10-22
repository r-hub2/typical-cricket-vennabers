#' @export
plot.venn_abers <- function(x, which = "map", loss = "log", ...) {
  which <- .plot_which(which)
  .validate_choice(loss, c("log", "brier"), "loss")
  if (!is.null(x$classes)) return(.plot_multiclass(x, which, loss, ...))
  .plot_fit_panels(list(.fit_map_table(x, .probability_tables(loss))), "", which, .probability_family(loss), ...)
}

#' @export
plot.cross_venn_abers <- function(x, which = "map", loss = "log", ...) {
  which <- .plot_which(which)
  .validate_choice(loss, c("log", "brier"), "loss")
  if (!is.null(x$classes)) return(.plot_multiclass(x, which, loss, ...))
  .plot_fit_panels(list(.fit_map_table(x, .probability_tables(loss))), "", which, .probability_family(loss), ...)
}

#' @export
plot.venn_abers_regression <- function(x, which = "map", ...) {
  .plot_regression(x, which, ...)
}

#' @export
plot.cross_venn_abers_regression <- function(x, which = "map", ...) {
  .plot_regression(x, which, ...)
}

# Workflow fits show their stored VA calibrators; postprocessing models are
# source models before them. A shared fit delegates to its one calibrator, a
# pairwise fit draws one panel per pair component.
#' @export
plot.vennabers_fit <- function(x, which = "map", loss = x$loss, ...) {
  requested <- x$parameters$requested
  components <- .fit_components(x)
  calibrator_classes <- c("venn_abers", "cross_venn_abers", "venn_abers_regression", "cross_venn_abers_regression")
  if (!length(components) ||
      !all(vapply(components, function(part) inherits(part$calibrator, calibrator_classes), TRUE))) {
    stop("`x` must be a fitted calibrated_workflow() or calibrate_workflow().", call. = FALSE)
  }
  if (identical(requested$mode, "regression")) {
    if (!missing(loss)) stop("Regression does not accept `loss`.", call. = FALSE)
    return(plot(x$calibrator, which, ...))
  }
  if (!identical(requested$model_training, "pairwise")) return(plot(x$calibrator, which, loss, ...))
  which <- .plot_which(which)
  .validate_choice(loss, c("log", "brier"), "loss")
  .validate_pair_components(x)
  # `pair_id` leads the identity columns, followed by the fold and the classes.
  pairs <- lapply(unname(components), `[[`, "class_levels")
  tables <- lapply(unname(components), function(part) {
    .identify_map(.fit_map_table(part$calibrator, .probability_tables(loss)), list(pair_id = part$pair_id),
                  list(class_1 = part$class_levels[1L], class_2 = part$class_levels[2L]))
  })
  .plot_fit_panels(tables, .pair_titles(pairs), which, .probability_family(loss), ...)
}

# The chosen diagrams, always drawn in the order of `choices`: for a fit
# diagnosis map before width.
.plot_which <- function(which, choices = c("map", "width")) {
  if (!is.character(which) || !is.vector(which) || !length(which) ||
      anyNA(which) || anyDuplicated(which) || !all(which %in% choices)) {
    stop("`which` must be one or more distinct values of: ",
         paste(choices, collapse = ", "), ".", call. = FALSE)
  }
  choices[choices %in% which]
}

.stop_stored_map <- function() {
  stop("`x` has no stored calibration map (`knots`, `lower`, `upper`); ",
       "refit it with vennabers 1.0.0.", call. = FALSE)
}

.validate_stored_map <- function(calibrator) {
  if (!all(c("knots", "lower", "upper") %in% names(calibrator))) .stop_stored_map()
}

# The map table of a binary or regression calibrator, Cross fits leading
# with their stored fold IDs and stacking the map of every fold calibrator.
# `map_tables` turns the checked calibrators into one table each.
.fit_map_table <- function(x, map_tables) {
  cross <- inherits(x, c("cross_venn_abers", "cross_venn_abers_regression"))
  if (cross && (is.null(x$folds) || length(x$calibrators) != length(x$folds))) .stop_stored_map()
  calibrators <- if (cross) x$calibrators else list(x)
  for (calibrator in calibrators) .validate_stored_map(calibrator)
  tables <- map_tables(calibrators)
  if (!cross) return(tables[[1L]])
  tables <- lapply(seq_along(x$folds), function(k) {
    data.frame(fold = rep(x$folds[k], nrow(tables[[k]])), tables[[k]])
  })
  do.call(rbind, tables)
}

# One panel per OvR class or OvO pair; identity columns follow the fold and
# precede the values, named as in venn_abers_details(). Each subproblem
# stores its folds in the order they occur among its own rows, so Cross folds
# are put in the stored order of the fit to keep their colors across panels.
.plot_multiclass <- function(x, which, loss, ...) {
  ovr <- identical(x$method, "ovr")
  problems <- if (ovr) as.list(x$classes) else x$pairs
  if (length(x$calibrators) != length(problems)) .stop_stored_map()
  tables <- lapply(seq_along(problems), function(i) {
    table <- .fit_map_table(x$calibrators[[i]], .probability_tables(loss))
    if (!is.null(table$fold)) table <- table[order(match(table$fold, x$folds)), ]
    problem <- problems[[i]]
    .identify_map(table, list(), if (ovr) list(class = problem) else
      list(class_1 = problem[1L], class_2 = problem[2L]))
  })
  .plot_fit_panels(tables, .pair_titles(problems), which, .probability_family(loss), ...)
}

# A map table with identity columns around its fold: `outer` before it,
# `inner` after it, both before the values.
.identify_map <- function(table, outer, inner) {
  n <- nrow(table)
  lead <- names(table) == "fold"
  data.frame(c(lapply(outer, rep, n), table[lead], lapply(inner, rep, n), table[!lead]))
}

# Panel title suffixes, ": <class>" or ": <class_1> vs <class_2>".
.pair_titles <- function(problems) paste0(": ", vapply(problems, paste, "", collapse = " vs "))

# One row per knot plus the range ends, valued by predict() so the knot
# convention holds exactly.
.probability_tables <- function(loss) {
  function(calibrators) lapply(calibrators, function(calibrator) {
    score <- sort(unique(c(0, calibrator$knots, 1)))
    map <- stats::predict(calibrator, score, loss = loss)
    data.frame(score = score, map, width = map$upper - map$lower)
  })
}

# How the maps of a family are drawn: the point column and its value between
# knots, the default axes, and the reference lines of the map.
.probability_family <- function(loss) {
  list(point = "probability", label = "Regularized probability",
       between = function(lower, upper, fold) .regularized_probability(lower, upper, loss),
       axes = list(xlim = c(0, 1), xlab = "Uncalibrated probability"),
       map = list(ylim = c(0, 1), yaxs = "i", ylab = "Calibrated probability"),
       reference = function() graphics::abline(a = 0, b = 1, lty = 3, col = "grey60"))
}

# Regression maps span the knots of all folds with a 5 % margin on each side,
# or a given `xlim`, keeping the knots within. References mark each `a` and
# `b`, the bounds of every VA value.
.plot_regression <- function(x, which, ...) {
  which <- .plot_which(which)
  graphical <- list(...)
  if ("loss" %in% names(graphical)) stop("Regression does not accept `loss`.", call. = FALSE)
  table <- .fit_map_table(x, function(calibrators) {
    ends <- graphical[["xlim"]]
    ends <- if (is.numeric(ends) && length(ends) == 2L && all(is.finite(ends))) sort(ends) else
      .regression_range(unlist(lapply(calibrators, `[[`, "knots")))
    lapply(calibrators, function(calibrator) {
      knots <- calibrator$knots
      score <- sort(unique(c(ends, knots[knots >= ends[1L] & knots <= ends[2L]])))
      map <- stats::predict(calibrator, score)
      data.frame(score = score, map[c("lower", "upper", "prediction")], width = map$upper - map$lower)
    })
  })
  calibrators <- if (inherits(x, "cross_venn_abers_regression")) x$calibrators else list(x)
  references <- unique(unlist(lapply(calibrators, `[`, c("a", "b"))))
  family <- list(
    point = "prediction", label = "Point prediction",
    between = function(lower, upper, fold) {
      calibrator <- calibrators[[if (is.null(fold)) 1L else match(fold, x$folds)]]
      .regression_point(lower, upper, calibrator$a, calibrator$b)
    },
    axes = list(xlim = range(table$score), xlab = "Uncalibrated prediction"),
    map = list(ylim = range(references), ylab = "Calibrated prediction"),
    reference = function() graphics::abline(h = references, lty = 2, col = "grey60"))
  .plot_fit_panels(list(table), "", which, family, ...)
}

# The knot range widened by 5 % of its span on each side, or of the absolute
# value of a single knot (0.05 at zero), within the finite doubles.
.regression_range <- function(knots) {
  ends <- range(knots)
  margin <- (ends[2L] / 2 - ends[1L] / 2) / 10
  if (margin == 0) margin <- if (ends[1L] == 0) 1 / 20 else abs(ends[1L]) / 20
  pmin(pmax(ends + c(-margin, margin), -.Machine$double.xmax), .Machine$double.xmax)
}

# Draws the chosen diagrams of each map table, one panel each, overlaying fold
# maps. Tables fill a grid row by row, and below each row of maps follow the
# widths of the same tables. Returns the stacked tables invisibly whatever
# `which` holds.
.plot_fit_panels <- function(tables, titles, which, family, ...) {
  rows <- ceiling(sqrt(length(tables)))
  columns <- ceiling(length(tables) / rows)
  if (length(tables) * length(which) > 1L) {
    old <- graphics::par(mfrow = c(rows * length(which), columns))
    on.exit(graphics::par(old))
  }
  for (start in seq(1L, length(tables), by = columns)) {
    row <- start:min(start + columns - 1L, length(tables))
    for (diagram in which) {
      for (i in row) .plot_fit_diagram(tables[[i]], titles[i], diagram, family, ...)
      if (diagram != which[length(which)]) for (i in seq_len(columns - length(row))) graphics::plot.new()
    }
  }
  table <- do.call(rbind, tables)
  rownames(table) <- NULL
  invisible(table)
}

.plot_fit_diagram <- function(table, title, diagram, family, ...) {
  maps <- if (is.null(table$fold)) list(table) else split(table, factor(table$fold, unique(table$fold)))
  paths <- lapply(maps, .map_paths, family = family)
  colors <- if (length(maps) > 1L) grDevices::hcl.colors(length(maps), "Dark 3") else "black"
  width <- max(unlist(lapply(paths, function(path) path$width$y)), .Machine$double.eps)
  vertical <- if (diagram == "map") family$map else
    list(ylim = c(0, 1.05 * width), yaxs = "i", ylab = "Upper minus lower")
  options <- utils::modifyList(c(
    list(x = 0:1, y = 0:1, type = "n", xaxs = "i"), family$axes, vertical,
    list(main = paste0(if (diagram == "map") "Calibration map" else "VA width", title))
  ), list(...))
  do.call(graphics::plot, options)
  if (diagram == "map") family$reference()
  for (k in seq_along(paths)) {
    curves <- if (diagram == "map") paths[[k]][c("lower", "upper", family$point)] else paths[[k]]["width"]
    for (curve in names(curves)) {
      graphics::lines(curves[[curve]], col = colors[k], lty = if (curve == family$point) 2 else 1)
    }
  }
  .fit_legend(diagram, family$label, names(maps), colors)
}

# Exact drawing paths of a map table. Between neighbouring rows the lower
# value of the left row and the upper value of the right row hold, so each
# curve is constant there and jumps vertically to its value at a knot.
.map_paths <- function(map, family) {
  last <- nrow(map)
  lower <- map$lower[-last]
  upper <- map$upper[-1L]
  between <- list(lower = lower, upper = upper, width = upper - lower)
  between[[family$point]] <- family$between(lower, upper, map$fold[1L])
  x <- c(rbind(map$score[-last], map$score[-last], map$score[-1L]), map$score[last])
  lapply(stats::setNames(names(between), names(between)), function(curve) {
    list(x = x, y = c(rbind(map[[curve]][-last], between[[curve]], between[[curve]]), map[[curve]][last]))
  })
}

.fit_legend <- function(diagram, point, folds, colors) {
  if (length(folds)) folds <- paste("Fold", folds)
  if (diagram == "map") {
    labels <- c("VA bounds", point)
    graphics::legend("topleft", c(labels, folds), lty = c(1, 2, rep(1, length(folds))),
                     col = c("black", "black", if (length(folds)) colors), bty = "n", cex = .8)
  } else if (length(folds)) {
    graphics::legend("topright", folds, lty = 1, col = colors, bty = "n", cex = .8)
  }
}
