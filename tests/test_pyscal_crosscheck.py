# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
"""Compare the native engine with the real pyscal package (skipped when pyscal is not installed)."""
import pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
import pytest

pytest.importorskip("pyscal")
import core as C
from pyscal_bridge import crosscheck_wo


@pytest.mark.parametrize("p", [
    C.WOParams(),
    C.WOParams(swl=0.15, swcr=0.2, sorw=0.3, krwend=0.6, kroend=0.8, krw=C.KrModel("corey", n=2.0), kro=C.KrModel("corey", n=3.5)),
    C.WOParams(krw=C.KrModel("let", L=2.0, E=3.0, T=1.5), kro=C.KrModel("let", L=2.5, E=2.0, T=1.2)),
])
def test_native_matches_pyscal(p):
    diff, msg = crosscheck_wo(p, 0.01)
    assert diff is not None, msg          # the bridge itself failing is a test failure when pyscal is present
    assert diff < 1e-6, msg
