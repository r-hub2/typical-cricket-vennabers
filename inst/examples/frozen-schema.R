# Keep the schema from BEFORE training. Later data are explicitly encoded with it.
source(system.file("examples", "assertions.R", package = "vennabers"))
# Separate earlier model/recipe learning, earlier postprocessing, new VA
# calibration, and independent assessment. Never rebuild a schema from calibration.
packages <- c("workflows", "parsnip", "recipes", "tailor", "probably", "nnet", "mgcv",
              "rsample", "hardhat", "rlang", "tibble", "withr", "generics", "tidyselect")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) local({
  set.seed(615)
  for (method in c("binary", "ovo", "ovr", "regression")) {
    original <- data.frame(x = rnorm(400), z = rnorm(400))
    original$y <- if (method == "regression") 2 + original$x + cos(1:400) else {
      classes <- if (method == "binary") c("class", "vaL0001") else c("class", "a b", "vaL0001")
      factor(rep(classes, length.out = 400), levels = classes)
    }
    names(original) <- c("vaV0002", "a-b", ".pred_class")
    schema <- if (method == "regression") vennabers::venn_abers_schema(original[1:150, ], ".pred_class") else
      if (method == "binary") vennabers::venn_abers_schema(original[1:150, ], ".pred_class", event_level = "second") else
        vennabers::venn_abers_schema(original[1:150, ], ".pred_class", classification = method)
    training <- vennabers::venn_abers_encode(original[1:150, ], schema)
    earlier_postprocessing <- vennabers::venn_abers_encode(original[151:250, ], schema)
    spec <- if (method == "regression") parsnip::linear_reg(engine = "lm") else
      if (method == "binary") parsnip::logistic_reg(engine = "glm") else
        parsnip::multinom_reg(penalty = .1) |> parsnip::set_engine("nnet", trace = FALSE)
    base <- workflows::workflow() |>
      workflows::add_recipe(recipes::recipe(vaV0003 ~ ., training) |>
        recipes::step_normalize(recipes::all_numeric_predictors())) |>
      workflows::add_model(spec) |>
      workflows::add_tailor(if (method == "regression") tailor::tailor() |>
        tailor::adjust_numeric_calibration(method = "linear", smooth = FALSE) else tailor::tailor() |>
          tailor::adjust_probability_calibration(method = if (method == "binary") "logistic" else "multinomial",
                                                smooth = (method != "binary")))
    earlier <- generics::fit(base, training, data_calibration = earlier_postprocessing)
    # Store earlier and schema together if calibration will happen in a later session.
    calibration <- vennabers::venn_abers_encode(original[251:350, ], schema)
    calibrated <- vennabers::calibrate_workflow(earlier, calibration, schema$encoded_outcome, schema = schema)
    new <- vennabers::venn_abers_encode(original[351:400, 1:2], schema)
    native <- predict(calibrated, new)
    decoded <- vennabers::venn_abers_decode(vennabers::venn_abers_details(calibrated, new), schema)
    assert_conditions(identical(calibrated$models[[1L]], earlier), identical(calibrated$schema, schema),
              decoded$outcome == ".pred_class", is.null(decoded$identity), is.null(decoded$aggregate))
    if (method != "regression") {
      decisions <- vennabers::venn_abers_decode(predict(calibrated, new, type = "class"), schema)
      stopifnot(identical(colnames(decoded$probability), schema$class_levels),
                identical(levels(decisions$class), schema$class_levels))
    }
    print(head(native))
  }
  # A workflow trained on original names continues to use those names, with
  # schema = NULL. Retrofitting a schema cannot rename learned state or repair
  # an old full-predict failure. Errors preserve names, stage and cause; there
  # is no automatic raw-response/engine fallback that bypasses postprocessing.
  # Retained bindings and learned metadata are checked, but detached matching
  # codes cannot prove provenance: callers must supply the actual training schema.
})
