# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
"""Golden tests: tables for fixed inputs must not change silently. After a deliberate change run tests/make_golden.py and review the diff."""
import pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
import numpy as np, pandas as pd, pytest
import make_golden as G


@pytest.mark.parametrize("name", list(G.cases()))
def test_golden(name):
    got = G.cases()[name]
    ref = pd.read_csv(G.OUT / f"{name}.csv")
    assert list(got.columns) == list(ref.columns)
    np.testing.assert_allclose(got.to_numpy(float), ref.to_numpy(float), rtol=1e-8, atol=1e-9)
