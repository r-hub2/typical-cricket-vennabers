# Selected predictors, in carrier order, then outcome. Row names and container
# classes are deliberately absent; column types/attributes and row order remain.
.local_split_data <- function(data, columns) {
  if (!is.data.frame(data)) stop("Split data binding must be a data frame.", call. = FALSE)
  .validate_names(names(data), "split data column")
  if (!all(columns %in% names(data))) {
    stop("Split data binding is missing selected raw columns: ",
         paste(setdiff(columns, names(data)), collapse = ", "), ".", call. = FALSE)
  }
  lapply(stats::setNames(columns, columns), function(column) data[[column]])
}

.inductive_split <- function(split, data) {
  if (!inherits(split, "rsplit") ||
      !all(class(split) %in% c("rsplit", "initial_split", "mc_split"))) {
    stop("`splits` must supply an Inductive IID holdout from initial_split() or make_splits().", call. = FALSE)
  }
  .split_positions(split, data, "Inductive1")
}

.split_positions <- function(split, data, id) {
  schema <- attr(data, "vennabers_schema", exact = TRUE)
  if (!is.null(schema)) .schema_binding(split$data, schema)
  if (!identical(.local_split_data(split$data, names(data)),
                 .local_split_data(data, names(data)))) {
    stop("Split data binding does not match the local selected raw data (values, types, levels or row order).", call. = FALSE)
  }
  positions <- list(analysis = as.integer(split, data = "analysis"),
                    assessment = as.integer(split, data = "assessment"))
  n <- nrow(data)
  valid <- function(x) is.numeric(x) && is.vector(x) && length(x) > 0L &&
    all(is.finite(x)) && all(x == floor(x)) && all(x >= 1 & x <= n) && !anyDuplicated(x)
  if (!all(vapply(positions, valid, logical(1))) ||
      length(intersect(positions$analysis, positions$assessment)) ||
      !identical(sort(as.integer(unlist(positions, use.names = FALSE))), seq_len(n))) {
    stop(id, " must partition all rows into nonempty, disjoint training and calibration positions with unique valid indices.", call. = FALSE)
  }
  lapply(positions, as.integer)
}

.cross_splits <- function(splits, data) {
  if (!inherits(splits, "rset") ||
      !all(class(splits) %in% c("vfold_cv", "manual_rset", "rset", "tbl_df", "tbl", "data.frame")) ||
      !identical(sort(names(splits)), c("id", "splits")) || nrow(splits) < 2L) {
    stop("`splits` must supply a simple Cross IID rset with at least two folds and only splits/id columns; repeated, group and time plans are unsupported.", call. = FALSE)
  }
  .validate_names(splits$id, "fold ID")
  folds <- lapply(seq_len(nrow(splits)), function(i) {
    split <- splits$splits[[i]]
    if (!inherits(split, "rsplit") || !all(class(split) %in% c("rsplit", "vfold_split"))) {
      stop(splits$id[i], " must be a simple Cross IID split from vfold_cv() or make_splits().", call. = FALSE)
    }
    .split_positions(split, data, splits$id[i])
  })
  held <- unlist(lapply(folds, `[[`, "assessment"), use.names = FALSE)
  if (!identical(sort(held), seq_len(nrow(data)))) {
    stop("Cross assessment folds must partition every local row exactly once.", call. = FALSE)
  }
  stats::setNames(folds, splits$id)
}

# A prepared mapping is checked in full before a resolver can be used for fitting.
#' @export
local_calibration_splits <- function(splits, columns) {
  .integration_packages()
  .validate_names(columns, "selected raw column")
  if (!is.list(splits) || inherits(splits, c("rsplit", "rset")) || !length(splits)) {
    stop("`splits` must be a nonempty list of prepared Inductive rsplit or Cross rset objects.", call. = FALSE)
  }
  snapshots <- plans <- schemas <- vector("list", length(splits))
  for (i in seq_along(splits)) {
    tryCatch({
      cross <- inherits(splits[[i]], "rset")
      schemas[[i]] <- lapply(if (cross) splits[[i]]$splits else list(splits[[i]]),
                            function(s) attr(s$data, "vennabers_schema", exact = TRUE))
      bound <- if (cross) splits[[i]]$splits[[1L]]$data else splits[[i]]$data
      snapshots[[i]] <- .local_split_data(bound, columns)
      data <- bound[columns]
      plans[[i]] <- if (cross) list(method = "cross", folds = .cross_splits(splits[[i]], data)) else
        list(method = "inductive", folds = list(Inductive1 = .inductive_split(splits[[i]], data)))
    }, error = function(e) stop("Prepared entry ", i, ": ", conditionMessage(e), call. = FALSE))
    # ponytail: quadratic setup in plan count; index snapshots if this becomes costly.
    for (j in seq_len(i - 1L)) {
      if (identical(snapshots[[i]], snapshots[[j]]) && !identical(plans[[i]], plans[[j]])) {
        stop("Conflicting prepared plans for identical local data in entries ", j, " and ", i, ".", call. = FALSE)
      }
    }
  }
  # Capture only the selected data and plans, never the construction frame.
  eval(quote(function(data) {
    if (!identical(names(data), columns)) {
      stop("Prepared plan columns must match all canonical selected raw columns, including outcome, in order.", call. = FALSE)
    }
    snapshot <- .local_split_data(data, columns)
    matches <- which(vapply(snapshots, identical, logical(1), y = snapshot))
    if (!length(matches)) stop("No prepared plan matches the local selected raw data.", call. = FALSE)
    for (bindings in schemas[matches]) for (schema in bindings) if (!is.null(schema)) {
      .schema_binding(data, schema)
      data <- .bind_schema(data, schema)
    }
    # All matching plans were proven identical; unique() does not pick a first hit.
    plan <- unique(plans[matches])[[1L]]
    if (plan$method == "cross") {
      rsample::manual_rset(lapply(plan$folds, rsample::make_splits, data = data), names(plan$folds))
    } else rsample::make_splits(plan$folds$Inductive1, data)
  }), envir = list(columns = columns, snapshots = snapshots, plans = plans, schemas = schemas),
  enclos = environment(local_calibration_splits))
}
