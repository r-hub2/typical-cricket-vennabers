\name{vennabers_examples}
\alias{vennabers_examples}
\title{Find and Run Venn-Abers Example Scripts}
\description{
Annotated example collections cover manual calibration with its fit
diagnoses and assessment graphics, training and frozen
workflows, learned postprocessing, and Grid/Bayes selection. Each collection
has a contents list with a separate entry for every example script.
}
\details{
Open the collections with \code{vignette()} or list all guides and examples
with \code{browseVignettes("vennabers")}. Vignettes are available after installing
a built source package. The original scripts are installed under
\code{system.file("examples", package = "vennabers")} even when installing
directly from a checkout without building vignettes.

The articles display the complete scripts without running their fits during
documentation builds. Use \code{file.show()} to read an installed script and
\code{source()} to run it, including its checks. Workflow, postprocessing and
tuning scripts use optional packages named in their code and DESCRIPTION;
install them explicitly before running those examples. Resampling, tuning
and composition matrices can take longer than the introductory guides.
}
\section{manual calibration}{
\code{vignette("examples-manual", package = "vennabers")}.

Scripts: \file{binary-options.R},
\file{cross-venn-abers.R},
\file{inductive-ovo.R},
\file{inductive-ovr.R},
\file{cross-multiclass.R},
\file{one-hot-classes.R},
\file{inductive-regression.R},
\file{cross-regression.R},
\file{calibration-graphics.R}.
}
\section{training and frozen workflows}{
\code{vignette("examples-workflows", package = "vennabers")}.

Scripts: \file{calibrated-workflow.R},
\file{encoded-binary.R},
\file{encoded-multiclass-regression.R},
\file{reference-splits.R},
\file{inductive-sizes.R},
\file{model-adapter.R},
\file{engine-options.R},
\file{full-model-workflow.R},
\file{full-multiclass-workflow.R},
\file{full-regression-workflow.R},
\file{pairwise-inductive.R},
\file{pairwise-cross.R},
\file{frozen-binary.R},
\file{frozen-multiclass-regression.R},
\file{frozen-schema.R}.
}
\section{learned postprocessing}{
\code{vignette("examples-postprocessing", package = "vennabers")}.

Scripts: \file{sigmoid-holdout.R},
\file{tailor-holdout.R},
\file{tailor-multiclass-holdout.R},
\file{tailor-regression-holdout.R},
\file{isotonic-regression-holdout.R},
\file{binary-ensemble.R},
\file{multiclass-ensemble.R},
\file{regression-ensemble.R},
\file{pairwise-postprocessing.R},
\file{full-postprocessing.R}.
}
\section{Grid and Bayes selection}{
\code{vignette("examples-tuning", package = "vennabers")}.

Scripts: \file{binary-grid.R},
\file{shared-multiclass-grid.R},
\file{shared-multiclass-bayes.R},
\file{pairwise-grid.R},
\file{pairwise-bayes.R},
\file{regression-grid.R},
\file{regression-bayes.R},
\file{postprocessing-grid.R},
\file{postprocessing-bayes.R}.
}
\seealso{
\code{\link{vennabers-package}}, \code{\link{venn_abers}},
\code{\link{calibrated_workflow}}, \code{\link{calibrate_workflow}},
\link{workflow_postprocessing}, \code{\link{venn_abers_tuning}}.
}
\examples{
# Locate an installed script without running its fits:
script <- system.file("examples", "binary-options.R", package = "vennabers")
script
# file.show(script)
# source(script)
}
\keyword{documentation}
