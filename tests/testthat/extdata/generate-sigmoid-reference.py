"""Run separately with scikit-learn==1.7.2; no Python needed by R tests."""
import csv
import json
from pathlib import Path

import numpy as np
import scipy
import sklearn
from sklearn.calibration import _SigmoidCalibration

assert sklearn.__version__ == "1.7.2"
root = Path(__file__).resolve().parent
cases = {
    "ordinary": ([.01, .12, .3, .4, .6, .7, .83, .99], [0, 0, 1, 0, 1, 0, 1, 1]),
    "rare": (np.linspace(.01, .99, 40), [0] * 39 + [1]),
    "constant": ([.6] * 8, [0] * 6 + [1] * 2),
    "two": ([.2, .8], [0, 1]),
    "saturated": ([0, 0, .1, .2, .8, .9, 1, 1], [0, 1, 0, 0, 1, 1, 0, 1]),
    # Many L-BFGS-B iterations: the optimizer path, not only its tolerances, matters.
    "concentrated": (.5 + .01 * np.sin(np.arange(1, 31)), [0, 1, 0, 0, 1, 1] * 5),
}
rows = []
for input_type in ["probability", "glm_link"]:
    for name, (x, y) in cases.items():
        x = np.array(x, dtype=float)
        if input_type == "glm_link":
            x = (x - .5) * (1e150 if name == "saturated" else 100)
        model = _SigmoidCalibration().fit(x, y)
        for role, values in [("P", x), ("prediction", np.linspace(min(x), max(x), 11))]:
            for i, (value, prediction) in enumerate(zip(values, model.predict(values))):
                rows.append([input_type, name, role, value, y[i] if role == "P" else "", prediction])
with (root / "sigmoid-reference.csv").open("w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["input", "case", "role", "score", "truth", "probability"])
    w.writerows(rows)
(root / "sigmoid-environment.json").write_text(json.dumps({
    "sklearn": sklearn.__version__, "numpy": np.__version__, "scipy": scipy.__version__,
    "source": "https://github.com/scikit-learn/scikit-learn/blob/1.7.2/sklearn/calibration.py",
    "generator": "tests/testthat/extdata/generate-sigmoid-reference.py", "absolute_tolerance": 1e-7,
}, indent=2) + "\n")
