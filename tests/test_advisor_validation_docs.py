# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
import pathlib, re, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
import numpy as np
import pytest
import advisor as A, core as C, validation as V, docs_content as DC, presets as P, deckcheck as D, exporters as X

ROOT = pathlib.Path(__file__).resolve().parent.parent


# ---------------------------------------------------------------- advisor
@pytest.mark.parametrize("key", list(A.PROFILES))
@pytest.mark.parametrize("kind", ["corey", "let"])
def test_profiles_valid_and_ordered_by_displacement(key, kind):
    prof = A.PROFILES[key]
    lo, ba, hi = (A.make_params(prof, c, kind) for c in ("low", "base", "high"))
    for p in (lo, ba, hi):
        p.validate()
    rf = [C.welge(C.wateroil(p, 0.01), 0.5, 3.0, swi=p.swl)["rf_ult"] for p in (lo, ba, hi)]
    assert rf[0] < rf[1] < rf[2], (key, rf)      # low case = pessimistic displacement, high = optimistic


def test_every_profile_has_basis_for_every_parameter():
    for p in A.PROFILES.values():
        assert set(p.ranges) == set(A.PARAMS) == set(p.basis)
        assert all(b in ("Craig 1971", "published", "judgement") for b in p.basis.values())
        assert p.lam[0] < p.lam[1] < p.lam[2] and p.pc_mult[0] > p.pc_mult[1] > p.pc_mult[2]


def test_triple_interpolation_hits_the_three_sets():
    pr = A.PROFILES["mixed_wet_sandstone"]
    lo, ba, hi = (A.make_params(pr, c) for c in ("low", "base", "high"))
    for t, ref in ((0.0, lo), (0.5, ba), (1.0, hi)):
        q = C.interpolate_triple(lo, ba, hi, t)
        assert (q.swl, q.sorw, q.krwend, q.krw.n) == pytest.approx((ref.swl, ref.sorw, ref.krwend, ref.krw.n))
    q = C.interpolate_triple(lo, ba, hi, 0.25)
    assert q.sorw == pytest.approx(0.5 * (lo.sorw + ba.sorw))
    assert C.interpolate_scalar_triple(2.0, 1.0, 0.5, 0.0, log=True) == pytest.approx(2.0)
    assert C.interpolate_scalar_triple(2.0, 1.0, 0.5, 0.75, log=True) == pytest.approx(np.sqrt(0.5))


def test_assess_flags_outside_values_and_let_guidance():
    p = A.make_params(A.PROFILES["water_wet_sandstone"], "base")
    msgs = A.assess(p, 0.5, 3.0, A.PROFILES["oil_wet_sandstone"], 2.0)
    assert any(s == "warn" for s, _ in msgs)
    bad = C.WOParams(krw=C.KrModel("let", L=0.5, E=1.0, T=0.2))
    assert any("LET" in m for s, m in A.assess(bad, 0.5, 3.0))
    ok = A.assess(p, 0.5, 3.0, A.PROFILES["water_wet_sandstone"])
    assert ok[-1][0] == "ok"


def test_presets_agree_with_craig_names():
    expect = {"Water-wet sandstone, waterflood": "water-wet", "Oil-wet carbonate, waterflood": "oil-wet"}
    for name, label in expect.items():
        s = P.PRESETS[name]["state"]
        swl = s["wo_swl"]
        p = C.WOParams(swl, max(s["wo_swcr"], swl), s["wo_sorw"], s["wo_krwend"], s["wo_kroend"], C.KrModel("corey", n=s["wo_krw_n"]), C.KrModel("corey", n=s["wo_kro_n"]))
        w = C.welge(C.wateroil(p, 0.01), s["mu_w"], s["mu_o"], swi=p.swl)
        assert A.craig_votes(p.swl, w["sw_cross"], p.krwend)["label"] == label


# ---------------------------------------------------------------- validation
def test_validation_suite_passes():
    res = V.run_all()
    bad = [(c.id, k.name, k.expected, k.got) for c in res for k in c.checks if not k.ok]
    assert not bad, bad
    assert V.summary(res)["FAIL"] == 0
    kinds = {c.kind for c in res}
    assert {"Published value", "Closed form", "Independent numerical", "Published guidance", "Synthetic recovery"} <= kinds


def test_validation_detects_a_broken_engine(monkeypatch):
    orig = C.KrModel.f
    monkeypatch.setattr(C.KrModel, "f", lambda self, x: orig(self, x) * 1.02)   # 2 % error in every kr shape
    assert V.case_pyscal_doc().status == "FAIL"


def test_fv_solver_conserves_mass_and_approaches_endpoint():
    p = C.WOParams()
    df = C.wateroil(p, 0.01)
    w = C.welge(df, 0.5, 5.0, swi=p.swl)
    r = V.fv_waterflood(df, 0.5, 5.0, p.swl, w["smax"], n_cells=400, checkpoints=(0.2,))
    assert r[0.2] == pytest.approx(0.2, abs=0.01)         # before breakthrough every PV injected displaces one PV of oil


# ---------------------------------------------------------------- documentation
def test_documentation_file_is_in_sync():
    assert (ROOT / "DOCUMENTATION.md").read_text() == DC.render_markdown(), "run: python docs_content.py > DOCUMENTATION.md"


def test_documentation_references_and_symbols_exist():
    import importlib
    for s in DC.SECTIONS:
        assert all(r in DC.REFERENCES for r in s["refs"]), s["id"]
        for c in s["code"]:
            m = re.match(r"^(\w+)\.(\w+)", c)
            if m and not c.split()[0].endswith(".py"):
                assert hasattr(importlib.import_module(m.group(1)), m.group(2)), c
    assert {s["id"] for s in DC.SECTIONS} >= {"kr-normalisation", "pc-models", "displacement", "hysteresis", "lab", "batch", "validation", "limits"}


def test_every_tab_has_a_reading_guide():
    tabs = ["Summary", "Relative permeability", "Capillary pressure", "End-point scaling", "Parameter advisor", "Uncertainty", "Curve-set interpolation",
            "Hysteresis", "Lab interpretation", "Fit SCAL data", "Animations", "Batch mode", "Validation", "Export and report"]
    for t in tabs:
        g = DC.GUIDES[t]
        assert set(g) == {"what", "read", "watch"} and all(len(v) > 40 for v in g.values()), t


def test_documented_formulas_match_the_code_numerically():
    """Spot-check that formulas printed in the documentation are the implemented ones."""
    x = 0.4
    assert C.KrModel("let", L=2.0, E=3.0, T=1.5).f(x) == pytest.approx(x ** 2 / (x ** 2 + 3.0 * (1 - x) ** 1.5))
    p = C.PcParams(kind="van-genuchten", pd=0.3, vg_n=2.5, smooth_cap=False, pcmax=1e6)
    s = 0.6; m = 1 - 1 / 2.5
    assert C.pc_curve(np.array([s]), 0.0, 0.0, p)[0] == pytest.approx(0.3 * (s ** (-1 / m) - 1) ** (1 / 2.5))
    sk = C.PcParams(kind="skjaeveland", cw=0.3, aw=0.5, co=0.15, ao=1.0, pcmax=99, pcmin=99)
    sw, swirr, sor = 0.5, 0.1, 0.2
    ref = 0.3 / ((sw - swirr) / (1 - swirr)) ** 0.5 - 0.15 / ((1 - sw - sor) / (1 - sor)) ** 1.0
    assert C.pc_curve(np.array([sw]), swirr, sor, sk)[0] == pytest.approx(ref)
    wo, go = C.wateroil(C.WOParams()), C.gasoil(C.GOParams())
    sw_, sg_ = 0.3, 0.2
    krow, krw = np.interp(sw_, wo.SW, wo.KROW), np.interp(sw_, wo.SW, wo.KRW)
    krog, krg = np.interp(sg_, go.SG, go.KROG), np.interp(sg_, go.SG, go.KRG)
    e = float(wo.KROW.max())
    assert C.three_phase_kro(wo, go, 0.1, np.array([sw_]), np.array([sg_]), "stone2")[0] == pytest.approx(max(e * ((krow / e + krw) * (krog / e + krg) - (krw + krg)), 0))


# ---------------------------------------------------------------- deck reader
def test_deck_reader_roundtrip_all_families():
    p, g = C.WOParams(), C.GOParams()
    reg = [dict(name="A", wo=C.add_pc_wo(C.wateroil(p), p, C.PcParams()), go=C.add_pc_go(C.gasoil(g), g, C.PcParams(pd=0.05, pcmax=2.0)))]
    for fam, kws in (("SWOF/SGOF", {"SWOF", "SGOF"}), ("SWOF/SLGOF", {"SWOF", "SLGOF"}), ("SWFN/SGFN/SOF3", {"SWFN", "SGFN", "SOF3"})):
        t = D.tables_from_eclipse(X.export_eclipse(reg, "METRIC", fam))
        assert set(t) == kws, (fam, set(t))
    t = D.tables_from_eclipse(X.export_eclipse(reg, "METRIC", "SWOF/SGOF"))
    assert np.allclose(t["SWOF"][0].SW.values, reg[0]["wo"].SW.values, atol=1e-6) and list(t["SWOF"][0].columns) == ["SW", "KRW", "KROW", "PC"]
