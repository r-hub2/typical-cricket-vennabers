#' @export
venn_abers_schema <- function(data, outcome, classification = "binary", event_level = "first") {
  .schema_data_names(data)
  .validate_names(outcome, "outcome")
  if (length(outcome) != 1L || !outcome %in% names(data)) {
    stop("`outcome` must name one column in `data`.", call. = FALSE)
  }
  y <- data[[outcome]]
  regression <- is.numeric(y)
  if (regression) {
    if (!missing(classification) || !missing(event_level)) {
      stop("Regression does not accept explicit classification or event_level options.", call. = FALSE)
    }
    .validate_regression_vector(y, outcome)
    classification <- event_level <- NULL
  } else {
    .validate_choice(classification, c("binary", "ovo", "ovr"), "classification")
    if (classification != "binary" && !missing(event_level)) {
      stop("`event_level` requires classification = 'binary'.", call. = FALSE)
    }
    .validate_choice(event_level, c("first", "second"), "event_level")
    if (!is.factor(y) || nlevels(y) < 2L || classification == "binary" && nlevels(y) != 2L) {
      stop("Schema outcome must be a factor with ",
           if (classification == "binary") "exactly two" else "at least two", " levels.", call. = FALSE)
    }
    .validate_names(levels(y), "class")
    if (classification != "binary") event_level <- NULL
  }
  columns <- stats::setNames(sprintf("vaV%04d", seq_along(data)), names(data))
  classes <- stats::setNames(sprintf("vaL%04d", seq_len(nlevels(y))), levels(y))
  structure(list(version = 1L, outcome = outcome, encoded_outcome = unname(columns[[outcome]]),
                 original_columns = names(data), columns = columns, class_levels = levels(y),
                 classes = classes, ordered = is.ordered(y), mode = if (regression) "regression" else "classification",
                 classification = classification, event_level = event_level,
                 event_class = if (!is.null(event_level)) .event_class(levels(y), event_level)),
            class = "vennabers_schema")
}

.schema_data_names <- function(data) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  .validate_names(names(data), "input column")
}

.validate_schema <- function(schema) {
  fields <- c("version", "outcome", "encoded_outcome", "original_columns", "columns",
              "class_levels", "classes", "ordered", "mode", "classification", "event_level", "event_class")
  invalid <- function() stop("Invalid or incomplete `schema` mapping.", call. = FALSE)
  if (!inherits(schema, "vennabers_schema") || !is.list(schema) ||
      !identical(names(schema), fields) || !identical(schema$version, 1L)) invalid()
  valid_names <- function(x) is.character(x) && is.vector(x) && length(x) > 0L &&
    !anyNA(x) && all(nzchar(x)) && !anyDuplicated(x)
  if (!valid_names(schema$original_columns) || !valid_names(schema$outcome) ||
      length(schema$outcome) != 1L || !schema$outcome %in% schema$original_columns ||
      !identical(schema$columns, stats::setNames(sprintf("vaV%04d", seq_along(schema$original_columns)), schema$original_columns)) ||
      !identical(schema$classes, stats::setNames(sprintf("vaL%04d", seq_along(schema$class_levels)), schema$class_levels)) ||
      !identical(schema$encoded_outcome, unname(schema$columns[[schema$outcome]])) ||
      !is.logical(schema$ordered) || length(schema$ordered) != 1L || is.na(schema$ordered)) invalid()
  if (identical(schema$mode, "regression")) {
    if (!is.null(schema$class_levels) || !is.null(schema$classification) || schema$ordered ||
        !is.null(schema$event_level) || !is.null(schema$event_class)) invalid()
  } else {
    if (!identical(schema$mode, "classification") || !valid_names(schema$class_levels) ||
        length(schema$class_levels) < 2L || !valid_names(schema$classification) ||
        length(schema$classification) != 1L || !schema$classification %in% c("binary", "ovo", "ovr")) invalid()
    if (schema$classification == "binary") {
      if (length(schema$class_levels) != 2L ||
          !identical(schema$event_level, "first") && !identical(schema$event_level, "second")) invalid()
      if (!identical(schema$event_class, .event_class(schema$class_levels, schema$event_level))) invalid()
    } else if (!is.null(schema$event_level) || !is.null(schema$event_class)) invalid()
  }
  invisible(schema)
}

# A schema fixes the encoded outcome, and the classification and event level it
# was built with; explicit options must agree with it. Returns both options.
.schema_options <- function(schema, outcome, classification, event_level, explicit) {
  .validate_schema(schema)
  if (!identical(outcome, schema$encoded_outcome) ||
      explicit[["classification"]] && !identical(classification, schema$classification) ||
      explicit[["event_level"]] && !identical(event_level, schema$event_level)) {
    stop("Outcome, classification and event_level must agree with `schema`.", call. = FALSE)
  }
  if (schema$mode == "classification") classification <- schema$classification
  if (identical(schema$classification, "binary")) event_level <- schema$event_level
  list(classification = classification, event_level = event_level)
}

.schema_binding <- function(x, schema) {
  bound <- attr(x, "vennabers_schema", exact = TRUE)
  if (!is.null(bound) && !identical(bound, schema)) {
    stop("Object is bound to a different schema.", call. = FALSE)
  }
  invisible(NULL)
}

.bind_schema <- function(x, schema) {
  if (!is.null(schema)) attr(x, "vennabers_schema") <- schema
  x
}

#' @export
venn_abers_encode <- function(data, schema) {
  .validate_schema(schema)
  .schema_data_names(data)
  .schema_binding(data, schema)
  unknown <- setdiff(names(data), schema$original_columns)
  required <- setdiff(schema$original_columns, schema$outcome)
  if (length(unknown)) stop("Unknown input columns: ", paste(unknown, collapse = ", "), ".", call. = FALSE)
  if (!all(required %in% names(data))) stop("Input is missing required predictors.", call. = FALSE)
  data <- data[intersect(schema$original_columns, names(data))]
  if (schema$outcome %in% names(data)) {
    y <- data[[schema$outcome]]
    if (schema$mode == "regression") {
      .validate_regression_vector(y, schema$outcome)
    } else {
      if (!is.factor(y)) stop("Schema outcome must be a factor.", call. = FALSE)
      if (any(!levels(y) %in% schema$class_levels)) stop("Unknown outcome classes.", call. = FALSE)
      data[[schema$outcome]] <- factor(unname(schema$classes[as.character(y)]),
                                      levels = unname(schema$classes), ordered = schema$ordered)
    }
  }
  names(data) <- unname(schema$columns[names(data)])
  .bind_schema(data, schema)
}

.schema_encoded_data <- function(data, schema, truth = FALSE) {
  .validate_schema(schema)
  .schema_data_names(data)
  .schema_binding(data, schema)
  required <- unname(schema$columns)
  if (!truth) required <- setdiff(required, schema$encoded_outcome)
  if (!all(required %in% names(data)) || any(!names(data) %in% schema$columns)) {
    stop("Data must contain the schema's encoded columns (truth is required for fit).", call. = FALSE)
  }
  if (schema$mode == "regression" && schema$encoded_outcome %in% names(data)) {
    .validate_regression_vector(data[[schema$encoded_outcome]], schema$outcome)
  } else if (schema$encoded_outcome %in% names(data) &&
      (!is.factor(data[[schema$encoded_outcome]]) ||
       !identical(levels(data[[schema$encoded_outcome]]), unname(schema$classes)) ||
       !identical(is.ordered(data[[schema$encoded_outcome]]), schema$ordered))) {
    stop("Encoded outcome must be a factor with class order and ordered status matching the schema; use venn_abers_encode().", call. = FALSE)
  }
  invisible(NULL)
}

.schema_decode_class <- function(x, schema, role) {
  if ((!is.factor(x) && !is.character(x)) ||
      any(!as.character(x)[!is.na(x)] %in% schema$classes) ||
      is.factor(x) && any(!levels(x) %in% schema$classes)) {
    stop("Unknown encoded classes in ", role, ".", call. = FALSE)
  }
  factor(schema$class_levels[match(as.character(x), schema$classes)],
         levels = schema$class_levels, ordered = schema$ordered)
}

#' @export
venn_abers_decode <- function(x, schema) {
  .validate_schema(schema)
  .schema_binding(x, schema)
  regression <- schema$mode == "regression"
  if (is.list(x) && !is.data.frame(x) &&
      (identical(names(x), c("prediction", "details")) ||
       identical(names(x), c("prediction", "details", "postprocessing")) ||
       regression && (identical(names(x), c("prediction", "details", "aggregate")) ||
                      identical(names(x), c("prediction", "details", "aggregate", "postprocessing"))))) {
    result <- venn_abers_decode(x$prediction, schema)
    details <- x$details
    .schema_data_names(details)
    .schema_binding(details, schema)
    classes <- if (regression) character() else switch(schema$classification,
      binary = "event_class", ovo = c("event_class", "class_1", "class_2"), ovr = c("event_class", "class"))
    required <- c("row", "fold", "lower", "upper", if (regression) "prediction" else "probability", classes)
    if (!all(required %in% names(details))) {
      stop("Incomplete prediction details.", call. = FALSE)
    }
    for (cl in classes) details[[cl]] <- as.character(.schema_decode_class(details[[cl]], schema, cl))
    attr(details, "vennabers_schema") <- NULL
    result$details <- details
    if (!is.null(x$postprocessing)) {
      decode_postprocessing <- function(p) {
        if (!is.null(p$inner_folds)) {
          p$inner_folds <- lapply(p$inner_folds, decode_postprocessing)
          return(p)
        }
        if (!regression) p$classes <- as.character(.schema_decode_class(p$classes, schema, "postprocessing classes"))
        if (!is.null(p$event_class)) p$event_class <- as.character(.schema_decode_class(p$event_class, schema, "postprocessing event"))
        if (!is.null(p$sigmoid_class)) p$sigmoid_class <- as.character(.schema_decode_class(p$sigmoid_class, schema, "sigmoid class"))
        p$columns <- schema$original_columns[match(p$columns, schema$columns)]
        names(p$column_types) <- p$columns
        if (!regression) p$column_types[[schema$outcome]]$levels <- p$classes
        if (!is.null(p$coefficients)) {
          coefficient_names <- if (is.matrix(p$coefficients)) colnames(p$coefficients) else names(p$coefficients)
          coefficient_classes <- match(coefficient_names, paste0(".pred_", schema$classes))
          mapped <- !is.na(coefficient_classes)
          coefficient_names[mapped] <- paste0(".pred_", schema$class_levels[coefficient_classes[mapped]])
          if (is.matrix(p$coefficients)) {
            colnames(p$coefficients) <- coefficient_names
            rownames(p$coefficients) <- as.character(.schema_decode_class(rownames(p$coefficients), schema, "postprocessing coefficients"))
          } else names(p$coefficients) <- coefficient_names
        }
        p
      }
      result$postprocessing <- lapply(x$postprocessing, decode_postprocessing)
    }
    if ("aggregate" %in% names(x)) {
      .schema_data_names(x$aggregate)
      .schema_binding(x$aggregate, schema)
      if (!all(c("row", "lower", "upper") %in% names(x$aggregate))) {
        stop("Incomplete regression aggregate.", call. = FALSE)
      }
      result$aggregate <- x$aggregate
      attr(result$aggregate, "vennabers_schema") <- NULL
    }
    return(result)
  }
  .schema_data_names(x)
  expected <- if (regression) character() else paste0(".pred_", unname(schema$classes))
  probabilities <- setdiff(names(x)[startsWith(names(x), ".pred_")], ".pred_class")
  if (length(probabilities) && !setequal(probabilities, expected)) {
    stop("Incomplete or unknown probability mapping.", call. = FALSE)
  }
  identities <- names(x)[grepl("^(\\.row|\\.config|id([1-9][0-9]*)?)$", names(x))]
  if (any(!names(x) %in% c(expected, if (regression) ".pred" else ".pred_class", schema$encoded_outcome, identities))) {
    stop("Unknown encoded output columns.", call. = FALSE)
  }
  result <- list(outcome = schema$outcome)
  if (regression) {
    if (".pred" %in% names(x)) {
      .validate_regression_vector(x$.pred, "prediction")
      result$prediction <- x$.pred
    }
    if (schema$encoded_outcome %in% names(x)) {
      .validate_regression_vector(x[[schema$encoded_outcome]], "truth")
      result$truth <- x[[schema$encoded_outcome]]
    }
  }
  if (length(probabilities)) {
    if (!all(vapply(x[expected], .valid_probability, logical(1))) || any(abs(rowSums(x[expected]) - 1) > 1e-8)) {
      stop("Invalid class probabilities.", call. = FALSE)
    }
    result$probability <- as.matrix(x[expected])
    storage.mode(result$probability) <- "double"
    colnames(result$probability) <- schema$class_levels
  }
  if (".pred_class" %in% names(x)) result$class <- .schema_decode_class(x$.pred_class, schema, "class")
  if (!regression && schema$encoded_outcome %in% names(x)) result$truth <- .schema_decode_class(x[[schema$encoded_outcome]], schema, "truth")
  if (length(identities)) {
    result$identity <- x[identities]
    attr(result$identity, "vennabers_schema") <- NULL
  }
  result
}
