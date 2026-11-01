# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
import pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
import pytest
import core as C, exporters as X, deckcheck as D


def regs():
    p, g, pc = C.WOParams(), C.GOParams(), C.PcParams()
    wo = C.add_pc_wo(C.wateroil(p), p, pc)
    p2 = C.WOParams(swcr=0.1)
    go = C.add_pc_go(C.gasoil(g), g, C.PcParams(pd=0.05, pcmax=2.0))
    return [dict(name="A", wo=wo, go=go), dict(name="B", wo=C.add_pc_wo(C.wateroil(p2), p2, pc), go=go)]


@pytest.mark.parametrize("units", ["METRIC", "FIELD", "SI-kPa"])
@pytest.mark.parametrize("family", ["SWOF/SGOF", "SWFN/SGFN/SOF3", "SWOF/SLGOF"])
def test_eclipse_roundtrip(units, family):
    r = regs()
    checks = D.verify_eclipse(X.export_eclipse(r, units, family), r, units, family)
    bad = [c for c in checks if not c.ok]
    assert not bad, bad


def test_eclipse_detects_tampering():
    r = regs()
    txt = X.export_eclipse(r, "METRIC", "SWOF/SGOF").replace("0.250000", "0.260000", 1)
    assert any(not c.ok for c in D.verify_eclipse(txt, r, "METRIC", "SWOF/SGOF"))


def test_gaswater_roundtrip():
    q = C.GWParams()
    r = [dict(name="gw", gw=C.add_pc_gw(C.gaswater(q), q, C.PcParams(pd=1.0, pcmax=20.0)))]
    for fam in ("SWFN/SGFN/SOF3", "SWOF/SGOF"):
        assert all(c.ok for c in D.verify_eclipse(X.export_eclipse(r, "METRIC", fam), r, "METRIC", fam))


def test_cmg_and_mrst_roundtrip():
    r = regs()
    assert all(c.ok for c in D.verify_cmg(X.export_cmg(r, "SI-kPa"), r, "SI-kPa"))
    assert all(c.ok for c in D.verify_mrst(X.export_mrst(r), r))
