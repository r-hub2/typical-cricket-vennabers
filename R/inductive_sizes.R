.inductive_sized <- function(settings) {
  !is.null(settings$training_size) || !is.null(settings$calibration_size)
}

.validate_inductive_sizes <- function(training_size, calibration_size) {
  sizes <- list(training_size = training_size, calibration_size = calibration_size)
  for (name in names(sizes)) {
    size <- sizes[[name]]
    if (is.null(size)) next
    if (rlang::is_call(size, "tune")) {
      stop("`", name, "` cannot be tuned; supply a fixed proportion or row count.", call. = FALSE)
    }
    if (!.is_proportion(size) && !.is_count(size)) {
      stop("`", name, "` must be NULL, one proportion strictly between 0 and 1, or one whole row count of at least 1.", call. = FALSE)
    }
  }
  if (!is.null(training_size) && !is.null(calibration_size) &&
      training_size < 1 && calibration_size < 1 && training_size + calibration_size > 1) {
    stop("`training_size` and `calibration_size` proportions must not sum to more than 1.", call. = FALSE)
  }
  !is.null(training_size) || !is.null(calibration_size)
}

# sklearn rounding on the local domain: training floors, calibration ceils and
# whole numbers are exact; a missing size is the complement. Totals are exact
# under stratification too (largest remainders per class); U is the rest.
.inductive_size_positions <- function(domain, settings, local) {
  n <- nrow(local)
  pool <- paste0("the ", n, " ", domain$rows_word)
  fail <- function(...) stop(if (!is.null(domain$pair_id)) paste0("Pair [", paste(domain$classes, collapse = ", "), "]: "),
                             ..., call. = FALSE)
  rows <- function(size, round) if (is.null(size)) NULL else if (size < 1) as.integer(round(n * size)) else as.integer(size)
  shown <- function(name, k) paste0("`", name, "` = ", format(settings[[name]]), " (", k, " rows)")
  training <- rows(settings$training_size, floor)
  calibration <- rows(settings$calibration_size, ceiling)
  given <- c(if (!is.null(training)) shown("training_size", training),
             if (!is.null(calibration)) shown("calibration_size", calibration))
  if (sum(training, calibration) > n) {
    fail(paste(given, collapse = " and "), " request ", sum(training, calibration), " rows, but only ", pool, " are available.")
  }
  if (is.null(training)) training <- n - calibration
  if (is.null(calibration)) calibration <- n - training
  for (role in c("training", "calibration")) if (get(role) == 0L) {
    fail(paste(given, collapse = " and "), " leave", if (length(given) == 1L) "s", " 0 ", role,
         " rows of ", pool, "; training and calibration must be nonempty.")
  }
  if (identical(domain$strata, FALSE)) {
    groups <- list(sample.int(n))
    take <- list(c(training, calibration))
  } else {
    y <- local[[domain$outcome]]
    groups <- lapply(levels(y), function(level) {
      members <- which(y == level)
      members[sample.int(length(members))]
    })
    counts <- lengths(groups)
    trained <- .largest_remainders(training, counts, counts)
    take <- Map(c, trained, .largest_remainders(calibration, counts, counts - trained))
  }
  roles <- c("analysis", "assessment", "unassigned")
  parts <- Map(function(members, k) split(members, factor(rep(roles, c(k, length(members) - sum(k))), roles)),
               groups, take)
  stats::setNames(lapply(roles, function(role) {
    sort(as.integer(unlist(lapply(parts, `[[`, role), use.names = FALSE)))
  }), roles)
}

# Proportional integer allocation of `total` over classes; ties keep level order.
.largest_remainders <- function(total, counts, capacity) {
  quota <- total * counts / sum(counts)
  take <- pmin(floor(quota), capacity)
  ranked <- order(floor(quota) - quota, seq_along(quota))
  while (sum(take) < total) for (k in ranked) {
    if (sum(take) < total && take[[k]] < capacity[[k]]) take[[k]] <- take[[k]] + 1
  }
  as.integer(take)
}
