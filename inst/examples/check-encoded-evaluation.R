# Run against an installed package. Two audits of one fixed outer workflow:
# check_outer_identities() for exact outer coverage, check_encoded_evaluation()
# for decoded names and reconstructed metrics.
source(system.file("examples", "assertions.R", package = "vennabers"))
# Exact outer coverage: no notes, every assessment row once per resample
# (including repeat IDs such as id2) with the single fixed .config, `metrics`
# finite unsummarized estimates per resample, and the required saved extracts
# returned in the row order of `outer`.
check_outer_identities <- function(result, outer, metrics) {
  ids <- setdiff(names(outer), "splits")
  notes <- tune::collect_notes(result)
  native <- as.data.frame(tune::collect_predictions(result, summarize = FALSE))
  scores <- as.data.frame(tune::collect_metrics(result, summarize = FALSE))
  expected <- do.call(rbind, lapply(seq_len(nrow(outer)), function(i) {
    rows <- as.integer(outer$splits[[i]], data = "assessment")
    cbind(as.data.frame(outer[rep(i, length(rows)), ids, drop = FALSE]), .row = rows, .config = "pre0_mod0_post0")
  }))
  key <- function(table, columns) do.call(paste, c(unname(as.list(table[columns])), sep = "\r"))
  have_ids <- all(ids %in% names(native)) && all(ids %in% names(scores))
  assert_conditions(nrow(notes) == 0L, have_ids, nrow(native) == nrow(expected),
    !anyDuplicated(key(native, names(expected))),
    setequal(key(native, names(expected)), key(expected, names(expected))),
    nrow(scores) == metrics * nrow(outer), !anyDuplicated(key(scores, c(ids, ".metric"))),
    setequal(key(scores, ids), key(outer, ids)), all(is.finite(scores$.estimate)),
    all(scores$.config == "pre0_mod0_post0"))
  collected <- as.data.frame(tune::collect_extracts(result))
  position <- match(key(outer, ids), key(collected, ids))
  assert_conditions(nrow(collected) == nrow(outer), !anyNA(position),
    all(collected$.config == "pre0_mod0_post0"), all(lengths(collected$.extracts) > 0L))
  invisible(list(predictions = native, metrics = scores, extracts = collected$.extracts[position]))
}
check_encoded_evaluation <- function(result, outer, schema) {
  notes <- tune::collect_notes(result)
  native <- tune::collect_predictions(result, summarize = FALSE)
  metrics <- tune::collect_metrics(result, summarize = FALSE)
  decoded <- vennabers::venn_abers_decode(native, schema)
  stopifnot(nrow(notes) == 0L, !anyNA(decoded$truth), all(is.finite(metrics$.estimate)))
  ids <- setdiff(names(outer), "splits")
  expected <- do.call(rbind, lapply(seq_len(nrow(outer)), function(i) {
    split <- outer$splits[[i]]
    rows <- as.integer(split, data = "assessment")
    stopifnot(!length(intersect(rows, as.integer(split, data = "analysis"))))
    cbind(as.data.frame(outer[rep(i, length(rows)), ids, drop = FALSE]),
          .row = rows, .config = "pre0_mod0_post0")
  }))
  actual <- decoded$identity
  assert_conditions(!anyNA(actual), !anyDuplicated(actual), nrow(actual) == nrow(expected),
            nrow(merge(actual, expected, by = names(expected))) == nrow(expected),
            !anyDuplicated(metrics[c(ids, ".config", ".metric")]), nrow(metrics) == 2L * nrow(outer),
            all(metrics$.estimator == if (schema$mode == "regression") "standard" else
              if (length(schema$classes) == 2L) "binary" else "multiclass"), all(metrics$.config == "pre0_mod0_post0"))
  for (i in seq_len(nrow(outer))) {
    matches <- function(table) Reduce(`&`, lapply(ids, function(id) table[[id]] == outer[[id]][i]))
    rows <- matches(actual)
    truth <- decoded$truth[rows]
    encoded_truth <- outer$splits[[i]]$data[[schema$encoded_outcome]][actual$.row[rows]]
    if (schema$mode == "regression") {
      stopifnot(identical(truth, encoded_truth))
      prediction <- decoded$prediction[rows]
      reconstructed <- c(rmse = yardstick::rmse_vec(truth, prediction),
                         mae = yardstick::mae_vec(truth, prediction))
    } else {
      stopifnot(identical(as.character(truth), schema$class_levels[match(as.character(encoded_truth), schema$classes)]))
      probability <- decoded$probability[rows, , drop = FALSE]
      event <- if (is.null(schema$event_level)) "first" else schema$event_level
      if (ncol(probability) == 2L) probability <- probability[, if (event == "first") 1L else 2L]
      reconstructed <- c(brier_class = yardstick::brier_class_vec(truth, probability, event_level = event),
                         mn_log_loss = yardstick::mn_log_loss_vec(truth, probability, event_level = event))
    }
    m <- metrics[matches(metrics), ]
    stopifnot(setequal(m$.metric, names(reconstructed)),
              all(abs(m$.estimate - reconstructed[m$.metric]) <= 1e-10))
  }
  invisible(list(decoded = decoded, metrics = metrics, notes = notes))
}
