# After package installation: Rscript inst/examples/one-hot-classes.R
library(vennabers)

# One-hot matrices represent class decisions the predictor has already made.
# venn_abers_one_hot() never chooses a class from probabilities itself.

# Binary: the logical class decision is TRUE (the event) at probability >= 0.5,
# so an exact 0.5 tie sets the TRUE column.
binary <- venn_abers(c(FALSE, TRUE, FALSE, TRUE), c(.2, .2, .8, .8))
new_scores <- c(.2, .5, .8)
probability <- predict(binary, new_scores)$probability
binary_class <- predict(binary, new_scores, type = "class")
binary_one_hot <- venn_abers_one_hot(binary_class)
print(cbind(probability, binary_one_hot))
stopifnot(probability[[2L]] == .5,
          identical(colnames(binary_one_hot), c("FALSE", "TRUE")),
          identical(binary_one_hot[, "TRUE"], as.double(binary_class)))

# Multiclass: columns follow the stored factor levels, even for classes that
# no row predicts, and an exact tie keeps the first stored level.
classes <- c("cat", "dog", "bird")
truth <- factor(rep(classes, 2), levels = classes)
p <- rbind(c(.7, .2, .1), c(.2, .6, .2), c(.1, .2, .7),
           c(.5, .3, .2), c(.3, .5, .2), c(.2, .3, .5))
colnames(p) <- classes
fit <- venn_abers(truth, p, method = "ovo")
new_p <- rbind(c(.6, .3, .1), c(.1, .4, .5), c(1, 1, 1) / 3)
colnames(new_p) <- classes
multiclass_class <- predict(fit, new_p, type = "class")
multiclass_one_hot <- venn_abers_one_hot(multiclass_class)
print(multiclass_one_hot)
# The third row is an exact tie, so its one sits in the first level "cat".
stopifnot(length(unique(predict(fit, new_p)[3L, ])) == 1L,
          identical(multiclass_one_hot[3L, ], c(cat = 1, dog = 0, bird = 0)))
stopifnot(identical(colnames(multiclass_one_hot), classes),
          all(rowSums(multiclass_one_hot) == 1),
          identical(multiclass_class, factor(classes[max.col(multiclass_one_hot)], levels = classes)))

# Missing decisions give missing rows and keep the row order.
print(venn_abers_one_hot(factor(c("bird", NA, "cat"), levels = classes)))

# Workflow, Grid, parsnip and frozen predictions pass their `.pred_class` table:
#   venn_abers_one_hot(predict(fitted_workflow, new_data, type = "class"))
# Schema-encoded tables keep the encoded class names; decode first for originals:
#   venn_abers_one_hot(venn_abers_decode(prediction, schema)$class)
