# Minimum versions match DESCRIPTION Suggests; each optional path checks its own
# packages at the call that needs them and never installs anything.
.integration_versions <- c(parsnip = "1.6.0", workflows = "1.3.0", hardhat = "1.4.3",
                           rsample = "1.3.2", recipes = "1.4.0", rlang = "1.2.0",
                           tibble = "3.3.1", withr = "3.0.3", generics = "0.1.4",
                           tidyselect = "1.2.1")
.tuning_versions <- c(tune = "2.1.0", dials = "1.4.3", yardstick = "1.4.0")
.tailor_versions <- c(tailor = "0.1.0", probably = "1.1.0")

.require_optional <- function(versions, user) {
  for (package in names(versions)) {
    version <- versions[[package]]
    if (!requireNamespace(package, quietly = TRUE) || utils::packageVersion(package) < version) {
      stop(user, " requires optional package ", sQuote(package), " >= ", version,
           "; install or update it explicitly.", call. = FALSE)
    }
  }
  invisible(NULL)
}

.integration_packages <- function() .require_optional(.integration_versions, "Workflow integration")

.tuning_packages <- function() .require_optional(.tuning_versions, "Tuning")

# Sigmoid and numeric isotonic postprocessing need no optional package.
.postprocessing_versions <- function(post) {
  if (is.null(post) || post$method %in% c("sigmoid", "isotonic_numeric")) character() else .tailor_versions
}

.postprocessing_packages <- function(post) .require_optional(.postprocessing_versions(post), "Postprocessing")

.registration <- new.env(parent = emptyenv())

# All package models share the data.frame fit interface, untouched predictors
# and the exported prediction entry points.
.register_model <- function(model, engines, fun) {
  if (model %in% parsnip::get_from_env("models")) {
    if (identical(.registration[[model]], parsnip::get_fit(model))) return(invisible(NULL))
    stop("A foreign parsnip registration already owns model '", model, "'; it will not be overwritten.", call. = FALSE)
  }
  parsnip::set_new_model(model)
  for (mode in c("classification", "regression")) parsnip::set_model_mode(model, mode)
  for (mode in c("classification", "regression")) for (engine in engines) {
    parsnip::set_model_engine(model, mode, engine)
    parsnip::set_dependency(model, engine, pkg = "vennabers", mode = mode)
    parsnip::set_fit(model, mode, engine,
                     list(interface = "data.frame", protect = c("x", "y"),
                          func = c(pkg = "vennabers", fun = fun), defaults = list()))
    parsnip::set_encoding(model, mode, engine,
                          list(predictor_indicators = "none", compute_intercept = FALSE,
                               remove_intercept = FALSE, allow_sparse_x = FALSE))
    for (type in if (mode == "regression") "numeric" else c("prob", "class")) {
      parsnip::set_pred(model, mode, engine, type,
                        list(pre = NULL, post = NULL,
                             func = c(pkg = "vennabers", fun = paste0("predict_vennabers_", type)),
                             args = list(object = quote(object$fit), new_data = quote(new_data))))
    }
  }
  .registration[[model]] <- parsnip::get_fit(model)
  invisible(NULL)
}

.register_vennabers <- function() {
  if (utils::packageVersion("parsnip") < .integration_versions[["parsnip"]]) return(invisible(NULL))
  .register_model("vennabers", c("inductive", "cross"), "fit_vennabers")
  .register_model_adapter()
}

.register_grid_proxy <- function() {
  .register_model("vennabers_grid", "proxy", "fit_vennabers_proxy")
}

# Callback adapters reuse the public fit entry point: their validated callbacks
# arrive as `parameters` (see model_adapter()).
.register_model_adapter <- function() {
  .register_model("vennabers_adapter", "callbacks", "fit_vennabers")
}

#' @export
update.vennabers_grid <- function(object, ..., fresh = FALSE) {
  replacements <- rlang::enquos(...)
  object$eng_args[names(replacements)] <- replacements
  object
}

# tune finalizes every model spec, also one without parameters. Callback
# settings are fixed in the functions, so nothing else may be replaced.
#' @export
update.vennabers_adapter <- function(object, ..., fresh = FALSE) {
  if (...length()) {
    stop("model_adapter() specs have no tunable or replaceable parameters; ",
         "bind fixed settings inside the `fit` and `predict` functions.", call. = FALSE)
  }
  object
}

.onLoad <- function(libname, pkgname) {
  setHook(packageEvent("parsnip", "onLoad"), function(...) .register_vennabers())
  if ("parsnip" %in% loadedNamespaces()) .register_vennabers()
}
