# Cross Venn-Abers numerical fixtures

`cross-calibration.csv` stores aligned calibration probability/truth by case,
fold identity, and within-fold row. `cross-predictions.csv` stores each matched
fold-model probability and its expected IVAP endpoints by case/fold/new row.
`cross-combined.csv` stores the expected combined event probability by case/new
row. Decimal values retain Python's round-trip float representation: no
precision rounding is applied. The `near_knots` case stores probabilities as
hexadecimal floats so Base R also preserves adjacent values at 0.5 on R 4.1
arm64, whose decimal parser can round the predecessor onto the knot. This
changes only serialization, not the reference inputs, outputs or tolerances.

The independent oracle is Python `venn-abers` **1.5.4**, commit
`0e7b4bc575ce230134576ba1c76688d46bf99ac7`:

- Numerical module SHA-256: `ea7816e8832eef3f0fa129d2bb80dffae4efcb785f37f083fc087136edb720a5`.
- Release sdist SHA-256: `a6bf5383e618ce42d0ecb2335bba0d6b104ecbe979041f10019b57572364bc69`.
- Source: https://github.com/ip200/venn-abers/tree/0e7b4bc575ce230134576ba1c76688d46bf99ac7
- Release: https://pypi.org/project/venn-abers/1.5.4/

Regenerate from the source checkout with:

```sh
python -m venv /tmp/cross-reference-env
/tmp/cross-reference-env/bin/pip install venn-abers==1.5.4
/tmp/cross-reference-env/bin/python tools/generate-cross-fixtures.py
```

For exact environment replication, install the dependency versions recorded
in `cross-environment.json`. That file preserves the actual generating Python,
platform, dependency versions, source identities, fixture hashes, seed, and
explicit coverage results. The generator checks the installed numerical
module hash before using it. Python is optional fixture-generation tooling;
normal installation and R tests require neither Python nor networking.

The generator exercises the pinned `VennAbersCV.predict_proba` path with
classification, `precision=None`, `loss="log"`, and `cv_ensemble=True`. It
injects deterministic calibration arrays and table predictors into that
version's internal state. Thus each new-observation column is matched to its
calibration fold and the actual upstream ensemble aggregation is exercised
without training estimators or creating fresh random splits. This internal
harness is reference tooling, not a proposed public interface.

Five deterministic cases reuse the v0.2 numerical research: unequal folds
with distinct probability scales; crossed aggregate endpoints; singleton
folds; input endpoints, knots, gaps and tails; and adjacent representable
probabilities. They also cover repeated probabilities, exact test ties, and
single-class folds. The unequal case distinguishes equal-weight endpoint
aggregation from averaging fold point predictions, size weighting, and pooled
calibration using the first model's test probabilities. Its alternative
outputs are recorded in the environment asset and checked during generation.

The 100 randomized cases use NumPy `default_rng(20260906)`, 2–10 folds, and 25
new observations per case. Fold sizes cycle through every integer from 1 to
200; continuous/tied probabilities and mixed/single-class truth are explicit
parts of generation. Both generation and ordinary R tests verify coverage
and distinct per-model columns. These ranges are test coverage, not API caps.

Ordinary parity requires maximum absolute error at most `1e-10`, with no
relative component, for fold endpoints and combined probabilities. Tests
check schema, dimensions, ordering, double type, finiteness and range separately.
They use installed public IVAP and cross interfaces without inspecting fitted
object fields.

Many-fold stability is checked separately against exact identical-fold
reductions by invoking the production aggregation helper: `(0.9, 0.91)` at
1,024 folds gives `91/101`; `(0.4, 0.7)` at 4,096 gives `7/13`; `(0.1, 0.2)`
at 2,048 gives `2/11`. The pinned Python product merger yields respectively
1, 0.5, and 0 for these research cases. Those unstable answers are not numerical
oracles and the tolerance is not widened to accommodate them. A one-fold
internal reduction does not change the public minimum of two folds.
