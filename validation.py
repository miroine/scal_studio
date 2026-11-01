# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
"""validation.py - comparison of the engine with published numbers, closed-form solutions and independent calculations.

Each case states WHAT it is compared with (the `kind` field), so a closed-form identity is never presented as a
published measurement. Kinds:
  Published value        a number printed in published documentation (pyscal docs, literature)
  Closed form            an analytic result derived from the cited model
  Independent numerical  a separately coded numerical solution of the same physics
  Published guidance     a published rule or range that the built-in profiles must respect
  Synthetic recovery     known inputs -> synthetic data -> interpretation must recover the inputs
  Self-consistency       internal consistency between two code paths
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

import advisor as A
import core as C

REFS = {
    "pyscal": "https://equinor.github.io/pyscal/introduction.html",
    "petrowiki": "https://petrowiki.spe.org/Relative_permeability_models",
    "lomeland": "https://www.jgmaas.com/SCA/2005/SCA2005-32.pdf",
    "craig": "Craig, F.F. (1971) The Reservoir Engineering Aspects of Waterflooding, SPE Monograph 3",
    "bc64": "Brooks, R.H. and Corey, A.T. (1964) Hydraulic Properties of Porous Media, Colorado State University Hydrology Paper 3",
    "welge": "Welge, H.J. (1952) A simplified method for computing oil recovery by gas or water drive, Trans. AIME 195",
    "land": "Land, C.S. (1968) Calculation of imbibition relative permeability for two- and three-phase flow from rock properties, SPEJ 8(2)",
    "jbn": "Johnson, Bossler and Naumann (1959) Calculation of relative permeability from displacement experiments, Trans. AIME 216",
}


@dataclass
class Check:
    name: str
    expected: float
    got: float
    tol: float          # absolute tolerance
    ok: bool = False

    def __post_init__(self):
        self.ok = bool(abs(self.got - self.expected) <= self.tol)


@dataclass
class Case:
    id: str
    title: str
    kind: str
    reference: str
    checks: list = field(default_factory=list)
    status: str = ""     # PASS | FAIL | SKIP
    note: str = ""

    def finish(self):
        if not self.status:
            self.status = "PASS" if self.checks and all(c.ok for c in self.checks) else "FAIL"
        return self


# ------------------------------------------------------------------ independent numerical solver
def fv_waterflood(df, mu_w, mu_o, swi, smax, n_cells=1000, checkpoints=(0.5, 1.0, 2.0, 3.0), cfl=0.4):
    """Explicit first-order upwind finite-volume solution of the 1-D Buckley-Leverett equation (no capillarity,
    no gravity). Independent of core.welge. Returns {pore volumes injected: oil produced [PV]}."""
    sg = np.linspace(swi, smax, 2000)
    fw = C.fractional_flow(np.interp(sg, df.SW, df.KRW), np.interp(sg, df.SW, df.KROW), mu_w, mu_o)
    fprime = float(np.max(np.gradient(fw, sg)))
    dx = 1.0 / n_cells
    dt = cfl * dx / fprime
    s = np.full(n_cells, swi)
    t, k, out, cps = 0.0, 0, {}, sorted(checkpoints)
    while k < len(cps):
        fl = np.interp(s, sg, fw)
        s = s - dt / dx * (fl - np.concatenate([[1.0], fl[:-1]]))   # inlet: pure water (fw = 1)
        t += dt
        if t >= cps[k]:
            out[cps[k]] = float(np.sum(s - swi) * dx)
            k += 1
    return out


# ------------------------------------------------------------------ cases
def case_pyscal_doc() -> Case:
    """pyscal documentation example: swl=0.1 swcr=0.1 sorw=0.05, Corey nw=2 now=3, krwend=kroend=1, delta_s=0.2.
    The documentation prints 'krw = krow @ sw=0.46856' and the first / last table rows."""
    p = C.WOParams(swl=0.1, swcr=0.1, sorw=0.05, krwend=1.0, kroend=1.0, krw=C.KrModel("corey", n=2.0), kro=C.KrModel("corey", n=3.0))
    fine = C.wateroil(p, 0.001)
    nodes = np.array([0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 1.0])          # pyscal grid for delta_s = 0.2 (steps plus 1 - sorw and 1)
    krw, kro = np.interp(nodes, fine.SW, fine.KRW), np.interp(nodes, fine.SW, fine.KROW)
    d = krw - kro
    i = int(np.where(np.diff(np.sign(d)) != 0)[0][0])
    cross = nodes[i] + (nodes[i + 1] - nodes[i]) * d[i] / (d[i] - d[i + 1])
    c = Case("V1", "pyscal documentation example: Corey water-oil crossover", "Published value", REFS["pyscal"])
    c.checks = [Check("krw = krow at Sw (printed: 0.46856)", 0.46856, float(cross), 2e-5),
                Check("first row krw (printed 0.0000000)", 0.0, float(krw[0]), 1e-9),
                Check("first row krow (printed 1.0000000)", 1.0, float(kro[0]), 1e-9),
                Check("last row krw (printed 1.0000000)", 1.0, float(krw[-1]), 1e-9),
                Check("last row krow (printed 0.0000000)", 0.0, float(kro[-1]), 1e-9)]
    c.note = "Reproduced on the same coarse grid pyscal used (delta_s = 0.2); the exact crossover on a fine grid is 0.4656."
    return c.finish()


def case_bl_closed_form() -> Case:
    """Quadratic kr, equal viscosities, Swc = Sor = 0: Welge tangent gives Sw_f = 1/sqrt(2), PV_bt = 2(sqrt(2) - 1)."""
    p = C.WOParams(swl=0.0, swcr=0.0, sorw=0.0, krwend=1.0, kroend=1.0, krw=C.KrModel("corey", n=2.0), kro=C.KrModel("corey", n=2.0))
    w = C.welge(C.wateroil(p, 0.005), 1.0, 1.0, swi=0.0)
    s2 = 2 ** 0.5
    c = Case("V2", "Buckley-Leverett with quadratic kr: shock saturation and breakthrough", "Closed form", REFS["welge"])
    c.checks = [Check("Shock saturation Swf = 1/sqrt(2)", 1 / s2, w["swf"], 1e-3),
                Check("Fractional flow at shock = 1/(2 - sqrt(2))", 0.5 / (2 - s2), w["fwf"], 1e-3),
                Check("Breakthrough PV = 2 (sqrt(2) - 1)", 2 * (s2 - 1), w["pvi_bt"], 1e-3),
                Check("RF at breakthrough (Swi = 0) = PV_bt", 2 * (s2 - 1), w["rf_bt"], 1e-3)]
    c.note = "Derivation: f = S^2/(S^2 + (1-S)^2); tangent from the origin requires f/S = f'."
    return c.finish()


def case_bl_numerical() -> Case:
    """Oil recovery versus pore volumes injected: Welge solution vs an independent upwind finite-volume code."""
    p = C.WOParams(swl=0.1, swcr=0.1, sorw=0.2, krwend=0.4, kroend=0.9, krw=C.KrModel("corey", n=3.0), kro=C.KrModel("corey", n=2.5))
    df = C.wateroil(p, 0.01)
    w = C.welge(df, 0.5, 5.0, swi=0.1)
    fv = fv_waterflood(df, 0.5, 5.0, 0.1, w["smax"], n_cells=1000)
    c = Case("V3", "Waterflood recovery: Welge vs independent finite-volume simulation (1000 cells)", "Independent numerical", "Self-coded first-order upwind scheme")
    c.checks = [Check(f"Oil produced after {t:g} PV [PV]", C.oil_produced_pv(w, t), v, 0.005) for t, v in fv.items()]
    c.note = "Mobility ratio 5.6; numerical diffusion of the first-order scheme explains differences below 0.001 PV."
    return c.finish()


def case_brooks_corey() -> Case:
    """Brooks and Corey (1964): lambda = 2 reduces to the Corey (1954) exponent 4; Pc = pd * Se^(-1/lambda)."""
    krw2, krn2 = C.brooks_corey_kr(np.array([0.5]), 2.0)
    pc = C.pc_curve(np.array([0.5]), 0.0, 0.0, C.PcParams(kind="brooks-corey", pd=0.3, lam=2.0, smooth_cap=False))
    c = Case("V4", "Brooks-Corey (1964) relations", "Closed form", REFS["bc64"] + "; " + REFS["petrowiki"])
    c.checks = [Check("krw exponent at lambda = 2 is 4", 4.0, C.brooks_corey_krw_exponent(2.0), 1e-12),
                Check("krw(Se = 0.5, lambda = 2) = 0.5^4", 0.5 ** 4, float(krw2[0]), 1e-12),
                Check("kr_nw(Se = 0.5, lambda = 2) = (1-Se)^2 (1-Se^2)", 0.25 * 0.75, float(krn2[0]), 1e-12),
                Check("Pc(Se = 0.5) = pd * 0.5^(-1/2)", 0.3 * 0.5 ** -0.5, float(pc[0]), 1e-12)]
    return c.finish()


def case_land() -> Case:
    """Land (1968): boundary conditions of the trapping relation and its closed form at a mid reversal point."""
    s_crit, s_trap, s_max = 0.05, 0.30, 0.80
    mid = 0.5
    cc = 1 / (s_trap - s_crit) - 1 / (s_max - s_crit)
    expected_mid = s_crit + (mid - s_crit) / (1 + cc * (mid - s_crit))
    c = Case("V5", "Land (1968) trapped saturation", "Closed form", REFS["land"])
    c.checks = [Check("Reversal at maximum saturation gives the maximum trapped saturation", s_trap, C.land_trapped(s_max, s_crit, s_trap, s_max), 1e-12),
                Check("Reversal at critical saturation traps nothing above Sgc", s_crit, C.land_trapped(s_crit, s_crit, s_trap, s_max), 1e-12),
                Check("Mid-point value from the closed form", expected_mid, C.land_trapped(mid, s_crit, s_trap, s_max), 1e-12)]
    return c.finish()


def case_leverett_centrifuge() -> Case:
    """Leverett scaling and the centrifuge inlet pressure checked by hand calculation."""
    ow, hg, ab = C.FLUID_SYSTEMS["Oil-water (reservoir)"], C.FLUID_SYSTEMS["Mercury-air (lab, MICP)"], C.FLUID_SYSTEMS["Air-brine (lab)"]
    j1, j2 = C.j_function(2.0, ow, 150, 0.22), C.j_function(float(C.pc_convert(2.0, ow, hg)), hg, 150, 0.22)
    r1, r2 = 0.11 - 0.05, 0.11                                   # 11 cm outer radius, 5 cm core
    pc_hand = 0.5 * 1000.0 * (2 * np.pi * 3000 / 60) ** 2 * (r2 ** 2 - r1 ** 2) / 1e5
    c = Case("V6", "Leverett J-function invariance and centrifuge inlet pressure", "Closed form", "Leverett (1941); Hassler and Brunner (1945)")
    c.checks = [Check("J-function identical in oil-water and mercury-air", j1, float(j2), 1e-9 * abs(j1)),
                Check("Hg-air to air-brine factor = |sigma cos theta| ratio", abs(72.0 / (480 * np.cos(np.radians(140)))), float(C.pc_convert(1.0, hg, ab)), 1e-12),
                Check("Centrifuge inlet Pc at 3000 rpm, drho 1 g/cc [bar]", pc_hand, float(C.centrifuge_pc_inlet(3000, 11.0, 5.0, 1.0)), 1e-9)]
    return c.finish()


def case_guidance() -> Case:
    """Built-in benchmark profiles must obey the published Craig classification and the LET guidance."""
    c = Case("V7", "Benchmark profiles versus published guidance", "Published guidance", REFS["craig"] + "; " + REFS["lomeland"])
    expect = {"water_wet_sandstone": "water-wet", "tight_sandstone": "water-wet", "oil_wet_sandstone": "oil-wet", "carbonate_mixed_oil_wet": "oil-wet"}
    for key, label in expect.items():
        b = A.make_params(A.PROFILES[key], "base")
        w = C.welge(C.wateroil(b, 0.01), 0.5, 3.0, swi=b.swl)
        got = A.craig_votes(b.swl, w["sw_cross"], b.krwend)["label"]
        c.checks.append(Check(f"{A.PROFILES[key].title}: Craig class is {label}", 1.0, 1.0 if got == label else 0.0, 0.0))
    worst = []
    for n in (1.5, 2, 3, 4, 5):
        m = C.corey_to_let(n)
        worst.append(float(m.L >= 1 or n < 1.5) * float(m.E > 0) * float(m.T >= 0.5))
    c.checks.append(Check("LET equivalents of Corey 1.5-5 satisfy L >= 1, E > 0, T >= 0.5 (Lomeland 2005)", 1.0, min(worst), 0.0))
    x = np.linspace(0, 1, 101)
    err = max(float(np.max(np.abs(C.corey_to_let(n).f(x) - x ** n))) for n in (2.0, 3.0, 4.0))
    c.checks.append(Check("LET equivalent reproduces Corey within 0.03 in kr", 0.0, err, 0.03))
    return c.finish()


def case_jbn() -> Case:
    base = C.wateroil(C.WOParams(), 0.01)
    data = C.synthetic_core_flood(base, 0.5, 5.0, 0.10)
    jb = C.jbn(data, 0.10, 0.5, 5.0, 1.0)
    ref = np.interp(jb.SW, base.SW, base.KROW)
    err = float(np.median(np.abs(jb.KROW - ref) / np.maximum(ref, 1e-3)))
    c = Case("V8", "JBN interpretation recovers kro from noisy synthetic core-flood data", "Synthetic recovery", REFS["jbn"])
    c.checks = [Check("Median relative kro error below 15 %", 0.0, err, 0.15), Check("Points interpreted at least 8", 8.0, float(len(jb)), 1e6)]
    c.checks[1].ok = len(jb) >= 8
    c.note = "Early post-breakthrough points are the most noise sensitive, as is typical for JBN."
    return c.finish()


def case_centrifuge_fit() -> Case:
    d = C.demo_centrifuge()
    fp, rmse, s0 = C.fit_centrifuge(d.RPM.values, d.SW_AVG.values, 11.0, 5.0, 1.0, "brooks-corey", 0.0)
    c = Case("V9", "Centrifuge inversion recovers the Pc curve used to generate the data", "Synthetic recovery", "Hassler and Brunner (1945), exact radial forward model")
    c.checks = [Check("Entry pressure (true 0.5 bar)", 0.5, fp.pd, 0.05), Check("Pore-size index (true 1.6)", 1.6, fp.lam, 0.15),
                Check("Irreducible saturation (true 0.12)", 0.12, s0, 0.02)]
    return c.finish()


def case_pyscal_live() -> Case:
    """Native tables against the real pyscal package (runs only when pyscal is installed)."""
    c = Case("V10", "Native engine versus the installed pyscal package", "Independent numerical", REFS["pyscal"])
    try:
        from pyscal_bridge import PYSCAL_AVAILABLE, crosscheck_wo
    except Exception as exc:  # noqa: BLE001
        c.status, c.note = "SKIP", f"bridge unavailable ({exc})"
        return c
    if not PYSCAL_AVAILABLE:
        c.status, c.note = "SKIP", "pyscal is not installed in this environment; run `pip install pyscal` and re-run."
        return c
    for name, p in (("Corey", C.WOParams()), ("LET", C.WOParams(krw=C.KrModel("let", L=2.0, E=3.0, T=1.5), kro=C.KrModel("let", L=2.5, E=2.0, T=1.2)))):
        d, msg = crosscheck_wo(p, 0.01)
        c.checks.append(Check(f"{name}: max |native - pyscal| over Sw in [Swl, 1 - Sorw]", 0.0, 1.0 if d is None else d, 1e-6))
    return c.finish()


def case_export_roundtrip() -> Case:
    import deckcheck as D
    import exporters as X
    p, g, pc = C.WOParams(), C.GOParams(), C.PcParams()
    reg = [dict(name="A", wo=C.add_pc_wo(C.wateroil(p), p, pc), go=C.add_pc_go(C.gasoil(g), g, C.PcParams(pd=0.05, pcmax=2.0)))]
    bad = 0
    for units in ("METRIC", "FIELD", "SI-kPa"):
        for fam in ("SWOF/SGOF", "SWOF/SLGOF", "SWFN/SGFN/SOF3"):
            bad += sum(not k.ok for k in D.verify_eclipse(X.export_eclipse(reg, units, fam), reg, units, fam))
    bad += sum(not k.ok for k in D.verify_cmg(X.export_cmg(reg, "SI-kPa"), reg, "SI-kPa")) + sum(not k.ok for k in D.verify_mrst(X.export_mrst(reg), reg))
    c = Case("V11", "Exported decks parse back to the same numbers (all families and units)", "Self-consistency", "deckcheck.py independent parsers")
    c.checks = [Check("Failed round-trip checks", 0.0, float(bad), 0.0)]
    return c.finish()


ALL_CASES = [case_pyscal_doc, case_bl_closed_form, case_bl_numerical, case_brooks_corey, case_land, case_leverett_centrifuge,
             case_guidance, case_jbn, case_centrifuge_fit, case_pyscal_live, case_export_roundtrip]


def run_all() -> list:
    return [f() for f in ALL_CASES]


def summary(cases) -> dict:
    return {s: sum(c.status == s for c in cases) for s in ("PASS", "FAIL", "SKIP")}


if __name__ == "__main__":
    rs = run_all()
    for r in rs:
        print(f"{r.status:5s} {r.id:4s} [{r.kind}] {r.title}")
        for k in r.checks:
            print(f"        {'ok ' if k.ok else 'BAD'} {k.name}: expected {k.expected:.6g}, got {k.got:.6g} (tol {k.tol:g})")
    print(summary(rs))
