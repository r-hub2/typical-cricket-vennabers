assert_conditions <- function(...) {
  conditions <- as.list(substitute(list(...)))[-1L]
  environment <- parent.frame()
  for (i in seq_along(conditions)) {
    result <- eval(conditions[[i]], environment)
    if (!(is.logical(result) && !anyNA(result) && all(result))) {
      values <- mget(ls(environment, all.names = TRUE), environment, inherits = FALSE)
      scalars <- vapply(values, function(x) is.null(x) || (is.atomic(x) && length(x) == 1L), logical(1))
      members <- vapply(values, function(x) {
        for (name in names(values)) if (is.list(values[[name]])) {
          index <- which(vapply(values[[name]], identical, logical(1), x))
          if (length(index)) return(paste0(name, "[[", index[1L], "L]]"))
        }
        ""
      }, "")
      context <- setNames(vapply(seq_along(values), function(i) {
        if (nzchar(members[i])) members[i] else paste(deparse(values[[i]]), collapse = "")
      }, ""), names(values))
      context <- context[scalars | nzchar(members)]
      stop(paste0(paste(deparse(conditions[[i]]), collapse = ""), " is not TRUE\nContext: ",
                  paste(paste(names(context), context, sep = "="), collapse = ", ")),
           call. = FALSE)
    }
  }
  invisible()
}
