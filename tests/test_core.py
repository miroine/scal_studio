# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
import pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
import numpy as np
import pytest
import core as C
import exporters as X


def test_wateroil_endpoints():
    p = C.WOParams(); d = C.wateroil(p)
    assert d.KRW[d.SW <= p.swcr].max() == 0
    assert d.KRW.max() == pytest.approx(p.krwend) and d.KROW.max() == pytest.approx(p.kroend)
    assert d.KROW[d.SW >= 1 - p.sorw].max() == pytest.approx(0, abs=1e-12)
    assert d.KRW.is_monotonic_increasing and d.KROW.is_monotonic_decreasing


def test_let_reduces_to_linear():
    k = C.KrModel("let", L=1, E=1, T=1)
    assert k.f(np.array([0.25]))[0] == pytest.approx(0.25)


def test_gasoil_gaswater():
    g = C.gasoil(C.GOParams()); assert g.KRG.max() == pytest.approx(0.85) and g.KROG.iloc[0] == pytest.approx(0.9)
    w = C.gaswater(C.GWParams()); assert w.KRW.max() == pytest.approx(0.6) and w.KRG.max() == pytest.approx(0.5)


def test_invalid_endpoints_raise():
    with pytest.raises(ValueError):
        C.wateroil(C.WOParams(swl=0.5, swcr=0.5, sorw=0.6))


@pytest.mark.parametrize("m", ["stone1", "stone2", "baker"])
def test_three_phase_limits(m):
    wo, go = C.wateroil(C.WOParams()), C.gasoil(C.GOParams())
    k = C.three_phase_kro(wo, go, 0.1, np.array([0.1, 0.5]), np.array([0.0, 0.5]), m)
    assert k[0] == pytest.approx(0.9, abs=1e-6)      # So = 1 - Swl: full oil end-point
    assert k[1] == 0                                  # So = 0
    sof3 = C.sof3_table(wo, go, 0.1); assert sof3.SO.iloc[0] == 0 and sof3.SO.iloc[-1] == pytest.approx(0.9)


@pytest.mark.parametrize("kind", C.PC_MODELS)
def test_pc_monotone_decreasing(kind):
    sw = np.linspace(0.1, 0.8, 60); pc = C.pc_curve(sw, 0.1, 0.2, C.PcParams(kind=kind))
    assert np.all(np.diff(pc) <= 1e-9)


def test_pc_conversion_and_j_invariance():
    f = C.FLUID_SYSTEMS; ow, hg = f["Oil-water (reservoir)"], f["Mercury-air (lab, MICP)"]
    pc_ow = np.array([1.0, 2.0]); pc_hg = C.pc_convert(pc_ow, ow, hg)
    assert np.all(pc_hg > pc_ow)                                       # mercury needs far higher pressure
    j1 = C.j_function(pc_ow, ow, 100, 0.2); j2 = C.j_function(pc_hg, hg, 100, 0.2)
    assert j1 == pytest.approx(j2)                                     # J-function collapses fluid systems
    assert C.pc_from_j(j1, ow, 100, 0.2) == pytest.approx(pc_ow)
    assert C.pc_convert(1.0, ow, ow, 100, .2, 400, .2)[()] == pytest.approx(0.5)   # 4x k -> half Pc


def test_rescale_endpoints():
    p = C.WOParams(); t = C.wateroil(p); n = C.WOParams(swl=0.15, swcr=0.15, sorw=0.3, krwend=0.5, kroend=0.8)
    r = C.rescale_wo_table(t, p, n)
    assert r.KRW.max() == pytest.approx(0.5) and r.KROW.max() == pytest.approx(0.8)
    r3 = C.rescale_wo_table(t, p, n, "three-point", krwr_new=0.5, krwmax_new=0.8)
    assert r3.KRW.iloc[-1] == pytest.approx(0.8)


def test_rock_family_monotone_in_k():
    lo = C.WOParams(swl=0.25, swcr=0.25, sorw=0.35, krwend=0.2, kroend=0.7); hi = C.WOParams()
    fam = C.rock_type_family(lo, hi, 5, 5, 500, 0.12, 0.28)
    assert [r["params"].swl for r in fam] == sorted([r["params"].swl for r in fam], reverse=True)
    assert fam[0]["pc_scale"] > fam[-1]["pc_scale"]


def test_land_and_scanning_continuity():
    g = C.gasoil(C.GOParams()); h = C.HystSetup(g.SG.values, g.KRG.values, 0.05, 0.75, 0.30)
    assert C.land_trapped(0.75, 0.05, 0.30, 0.75) == pytest.approx(0.30)
    assert C.land_trapped(0.05, 0.05, 0.30, 0.75) == pytest.approx(0.05)
    for scheme in ("killough", "carlson"):
        for s_hy in (0.2, 0.4, 0.6, 0.75):
            c = C.scanning_curve(h, s_hy, scheme)
            assert c.KR.iloc[-1] == pytest.approx(float(h.kr_d(s_hy)), abs=1e-6)   # starts on drainage curve
            assert c.KR.min() >= -1e-12 and c.S.iloc[0] < s_hy
    imb = C.imbibition_gas_table(g, h); assert (imb.KRG <= g.KRG + 1e-9).all()    # hysteresis lowers krg


def test_pc_scanning_continuity():
    sw = np.linspace(0.1, 0.8, 50); d = C.pc_curve(sw, .1, .2, C.PcParams()); i = C.pc_curve(sw, .1, .2, C.PcParams(kind="let-two-sided"))
    s = C.pc_scanning(sw, d, i, 0.3, 0.8)
    assert s[sw >= 0.3][0] == pytest.approx(d[sw >= 0.3][0]) and s[-1] == pytest.approx(i[-1])


def test_welge_and_recovery():
    d = C.wateroil(C.WOParams()); w = C.welge(d, 0.5, 5.0, 0.1)
    assert w["mobility_ratio"] == pytest.approx((0.4 / 0.5) / (0.9 / 5.0))
    rec = C.recovery_curve(w, 3.0)
    assert rec.RF.is_monotonic_increasing and rec.RF.max() <= w["rf_ult"] + 1e-9
    assert np.interp(w["pvi_bt"], rec.PVI, rec.RF) == pytest.approx(w["rf_bt"], abs=5e-3)


def test_fit_recovers_parameters():
    p = C.WOParams(krw=C.KrModel("corey", n=2.0), kro=C.KrModel("corey", n=4.0)); d = C.wateroil(p, 0.03)
    f, rmse = C.fit_relperm(d, "corey", p)
    assert f.krw.n == pytest.approx(2.0, abs=0.05) and f.kro.n == pytest.approx(4.0, abs=0.05) and rmse < 1e-3


def test_exporters_roundtrip():
    wo = C.add_pc_wo(C.wateroil(C.WOParams(), 0.1), C.WOParams(), C.PcParams())
    txt = X.export_eclipse([dict(wo=wo)], "METRIC", "SWOF/SGOF", None)
    rows = [l.split() for l in txt.splitlines() if l.strip() and l.split()[0][0].isdigit()]
    assert len(rows) == len(wo) and float(rows[0][3]) == pytest.approx(5.0)
    assert float(X.export_cmg([dict(wo=wo)], "SI-kPa").split("*SWT")[1].splitlines()[1].split()[3]) == pytest.approx(500.0)


# ---------------------------------------------------------------- Pc upgrade
def test_smooth_cap_hits_pcmax_exactly_and_matches_bc_elsewhere():
    sw = np.linspace(0.1, 0.8, 200)
    sm = C.pc_curve(sw, 0.1, 0.2, C.PcParams(smooth_cap=True))
    hard = C.pc_curve(sw, 0.1, 0.2, C.PcParams(smooth_cap=False))
    assert sm[0] == pytest.approx(5.0) and hard[0] == pytest.approx(5.0)
    assert sm[-1] == pytest.approx(0.3 * ((1 - 0.0) / 1) ** 0 * ((0.8 - 0.1) / 0.9) ** -0.5, rel=0.01)
    assert np.max(np.abs(sm[40:] - hard[40:]) / hard[40:]) < 0.02      # same physics away from the cap


def test_skjaeveland_two_sided_and_clip():
    sw = np.linspace(0.1, 0.8, 100); p = C.PcParams(kind="skjaeveland", pcmax=6.0, pcmin=2.0)
    pc = C.pc_curve(sw, 0.1, 0.2, p)
    assert pc.max() == pytest.approx(6.0) and pc.min() == pytest.approx(-2.0) and (pc > 0).any() and (pc < 0).any()


def test_inverse_pc_and_saturation_height():
    p = C.PcParams(); sw = np.array([0.3, 0.5, 0.7])
    assert C.sw_from_pc(C.pc_curve(sw, 0.1, 0.2, p), 0.1, 0.2, p) == pytest.approx(sw, abs=2e-3)
    h = np.array([0.0, 5.0, 30.0, 100.0]); s = C.saturation_height(h, 250.0, 0.1, 0.2, p)
    assert s[0] == pytest.approx(1.0) and np.all(np.diff(s) <= 1e-12) and s[-1] < 0.2   # Sw falls with height
    assert C.saturation_height([30.0], 250.0, 0.1, 0.2, p, sigma_scale=2.0)[0] > s[2]   # larger sigma*cos(theta) holds more water at the same height -> higher Sw
    assert C.seal_column_height(3.0, 300.0) == pytest.approx(3.0e5 / (300 * 9.81))


def test_micp_conversion_and_fit_recovers_swirr():
    d = C.demo_micp(); ow = C.FLUID_SYSTEMS["Oil-water (reservoir)"]
    sw, pc = C.micp_to_pc(d.P_HG_BAR, d.S_HG, ow, 0.0)
    assert np.all(np.diff(sw) >= 0) and pc.max() < d.P_HG_BAR.max()     # sorted, lower than Hg pressure
    fp, rmse, r2, s0 = C.fit_pc(sw, pc, float(sw.min()), 0.0, "brooks-corey")
    assert abs(s0 - 0.08) < 0.02 and r2 > 0.9                              # true Swirr = 0.08
    _, _, r2_fixed, _ = C.fit_pc(sw, pc, float(sw.min()), 0.0, "brooks-corey", fit_swirr=False)
    assert r2 > r2_fixed + 0.5                                             # fixing Swirr at min(Sw) is far worse


def test_pc_checks_and_gas_oil_ramp():
    sw = np.linspace(0.1, 0.8, 30)
    assert C.pc_checks(sw, C.pc_curve(sw, 0.1, 0.2, C.PcParams()), C.PcParams()) == []
    assert any("monoton" in m for m in C.pc_checks(sw, sw))                # increasing Pc is flagged
    go = C.GOParams(sgcr=0.1); g = C.add_pc_go(C.gasoil(go, 0.02), go, C.PcParams(kind="let", pd=0.05, pcmax=1.5))
    assert g.PC.iloc[0] == 0 and np.all(np.diff(g.PC.values) >= -1e-12)    # starts at 0 and rises without a jump
    first = g.PC[(g.SG > 0) & (g.SG <= 0.02 + 1e-9)].iloc[0]; assert first < 0.5 * g.PC[g.SG >= 0.1].iloc[0]


# ---------------------------------------------------------------- v3: ensembles, interpolation, lab interpretation
def test_ensemble_reproducible_ordered_and_degenerate():
    spec = C.EnsembleSpec(n=80, seed=3)
    a = C.ensemble_wo(C.WOParams(), spec, 0.5, 5.0, 0.02, C.PcParams()); b = C.ensemble_wo(C.WOParams(), spec, 0.5, 5.0, 0.02, C.PcParams())
    assert np.array_equal(a["krw"], b["krw"])                                     # seed -> reproducible
    for key in ("krw", "kro", "pc"):
        p = C.percentile_bands(a[key]); assert np.all(p[10] <= p[50] + 1e-12) and np.all(p[50] <= p[90] + 1e-12)
    z = C.ensemble_wo(C.WOParams(), C.EnsembleSpec(n=4, d_swl=0, d_sorw=0, f_krwend=0, f_kroend=0, d_nw=0, d_no=0, f_pc=0), 0.5, 5.0, 0.02)
    assert np.allclose(z["krw"][0], z["krw"][3]) and z["samples"].rf_ult.nunique() == 1   # no spread -> identical
    pk = C.pick_cases(a); v = a["samples"].rf_ult.values
    assert abs(v[pk[50]] - np.median(v)) < 0.02 and v[pk[10]] < v[pk[50]] < v[pk[90]]


def test_ensemble_let_and_validity():
    base = C.WOParams(krw=C.KrModel("let", L=2, E=1.5, T=2), kro=C.KrModel("let", L=2, E=2, T=1.5))
    e = C.ensemble_wo(base, C.EnsembleSpec(n=40), 0.5, 5.0, 0.02)
    assert len(e["params"]) == 40 and all(p.swcr >= p.swl and p.swcr + p.sorw < 1 for p in e["params"])


def test_interpolation_between_curve_sets():
    a = C.WOParams(); b = C.WOParams(sorw=0.10, krwend=0.6, krw=C.KrModel("corey", n=2.0))
    for t, sorw in ((0.0, 0.20), (0.5, 0.15), (1.0, 0.10)):
        assert C.interpolate_wo(a, b, t).sorw == pytest.approx(sorw)
    with pytest.raises(ValueError):
        C.interpolate_wo(a, C.WOParams(krw=C.KrModel("let"), kro=C.KrModel("let")), 0.5)
    assert C.interpolation_factor(1e-5, 1e-6, 1e-4) == pytest.approx(0.5)
    assert C.interpolation_factor(5000, 1000, 100000, "linear") < 0.05 and C.interpolation_factor(1, 10, 100) == 0.0


def test_jbn_recovers_known_curves():
    true = C.wateroil(C.WOParams(), 0.01); dat = C.synthetic_core_flood(true, 0.5, 5.0, 0.1)
    r = C.jbn(dat, 0.1, 0.5, 5.0, kro_swi=0.9)
    mid = (r.KROW > 0.02) & (r.KRW > 0.01)
    err_w = np.abs(r.KRW[mid] / np.interp(r.SW[mid], true.SW, true.KRW) - 1); err_o = np.abs(r.KROW[mid] / np.interp(r.SW[mid], true.SW, true.KROW) - 1)
    assert mid.sum() >= 8 and np.median(err_w) < 0.15 and np.median(err_o) < 0.15
    with pytest.raises(ValueError):
        C.jbn(dat.iloc[1:], 0.1, 0.5, 5.0)                                          # first row must be PVI = 0


def test_centrifuge_inlet_pressure_and_fit():
    assert C.centrifuge_pc_inlet(4500, 11, 5, 1.0) == pytest.approx(9.44, abs=0.02)
    d = C.demo_centrifuge(); f, rmse, s0 = C.fit_centrifuge(d.RPM, d.SW_AVG, 11, 5, 1.0)
    assert abs(s0 - 0.12) < 0.03 and abs(f.lam - 1.6) / 1.6 < 0.15 and abs(f.pd - 0.5) < 0.1 and rmse < 0.01
    assert np.all(np.diff(C.centrifuge_avg_sw(d.RPM, 11, 5, 1.0, 0.12, 0.0, f)) <= 1e-9)    # Sw_avg falls with speed
