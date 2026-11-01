# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
"""
core.py - SCAL engine for SCAL Studio.

Pure numpy / pandas / scipy (no Streamlit, no Plotly) so it can be unit tested.
Relative-permeability parametrisations follow pyscal conventions
(Corey and LET, normalised saturations built from swl/swcr/sorw/sgcr/sorg).

Contents
  1. Relative permeability   (water-oil, gas-oil, gas-water, three-phase Stone I/II/Baker)
  2. Capillary pressure       (Brooks-Corey, LET, van Genuchten, LET two-sided, J-function,
                               fluid-system conversion, height above FWL, wettability indices)
  3. Endpoint scaling         (two-point / three-point table rescaling, rock-type families)
  4. Hysteresis               (Land trapping, Killough / Carlson scanning curves, Pc scanning)
  5. Displacement analytics   (fractional flow, Welge shock, recovery vs PVI)
  6. Fitting                  (Corey / LET fit to measured SCAL data)
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
import pandas as pd

BAR_TO_PSI = 14.5037738
BAR_TO_KPA = 100.0
MD_TO_M2 = 9.869233e-16


# --------------------------------------------------------------------------- #
# 1. Relative permeability
# --------------------------------------------------------------------------- #
def _clip01(x):
    return np.clip(np.asarray(x, dtype=float), 0.0, 1.0)


@dataclass
class KrModel:
    """Shape function on a normalised saturation x in [0, 1] (0 -> kr=0, 1 -> kr=1)."""

    kind: str = "corey"  # "corey" | "let"
    n: float = 2.0
    L: float = 2.0
    E: float = 2.0
    T: float = 2.0

    def f(self, x):
        x = _clip01(x)
        if self.kind == "corey":
            return x ** self.n
        num = x ** self.L
        den = num + self.E * (1.0 - x) ** self.T
        return np.divide(num, den, out=np.zeros_like(x), where=den > 0)


@dataclass
class WOParams:
    swl: float = 0.10
    swcr: float = 0.10
    sorw: float = 0.20
    krwend: float = 0.40
    kroend: float = 0.90
    krw: KrModel = field(default_factory=lambda: KrModel("corey", n=3.0))
    kro: KrModel = field(default_factory=lambda: KrModel("corey", n=2.5))

    def validate(self):
        if self.swcr < self.swl - 1e-9:
            raise ValueError("Swcr must be >= Swl")
        if self.swcr + self.sorw >= 1.0 or self.swl + self.sorw >= 1.0:
            raise ValueError("Swcr/Swl + Sorw must be < 1")
        if not (0 < self.krwend <= 1 and 0 < self.kroend <= 1):
            raise ValueError("Kr end-points must be in (0, 1]")


@dataclass
class GOParams:
    swl: float = 0.10
    sgcr: float = 0.05
    sorg: float = 0.15
    krgend: float = 0.85
    kroend: float = 0.90
    krg: KrModel = field(default_factory=lambda: KrModel("corey", n=2.0))
    kro: KrModel = field(default_factory=lambda: KrModel("corey", n=2.5))

    def validate(self):
        if self.sgcr + self.sorg + self.swl >= 1.0:
            raise ValueError("Sgcr + Sorg + Swl must be < 1")


@dataclass
class GWParams:
    swl: float = 0.20
    sgcr: float = 0.05
    sgr: float = 0.25  # residual (trapped) gas, sets where krw reaches its end-point
    krwend: float = 0.60
    krgend: float = 0.50
    krw: KrModel = field(default_factory=lambda: KrModel("corey", n=3.0))
    krg: KrModel = field(default_factory=lambda: KrModel("corey", n=2.0))

    def validate(self):
        if self.swl + self.sgcr + self.sgr >= 1.0:
            raise ValueError("Swl + Sgcr + Sgr must be < 1")


def _grid(lo, hi, h, extra=()):
    n = max(int(round((hi - lo) / h)) + 1, 3)
    g = np.concatenate([np.linspace(lo, hi, n), np.asarray(extra, float)])
    g = g[(g >= lo - 1e-12) & (g <= hi + 1e-12)]
    return np.unique(np.round(g, 10))


def wateroil(p: WOParams, h: float = 0.01) -> pd.DataFrame:
    p.validate()
    sw = _grid(p.swl, 1.0, h, [p.swcr, 1.0 - p.sorw])
    xw = (sw - p.swcr) / (1.0 - p.swcr - p.sorw)
    xo = (1.0 - sw - p.sorw) / (1.0 - p.sorw - p.swl)
    return pd.DataFrame({"SW": sw, "KRW": p.krwend * p.krw.f(xw), "KROW": p.kroend * p.kro.f(xo)})


def gasoil(p: GOParams, h: float = 0.01) -> pd.DataFrame:
    p.validate()
    sg = _grid(0.0, 1.0 - p.swl, h, [p.sgcr, 1.0 - p.swl - p.sorg])
    xg = (sg - p.sgcr) / (1.0 - p.sgcr - p.sorg - p.swl)
    xo = (1.0 - sg - p.sorg - p.swl) / (1.0 - p.sorg - p.swl)
    return pd.DataFrame({"SG": sg, "KRG": p.krgend * p.krg.f(xg), "KROG": p.kroend * p.kro.f(xo)})


def gaswater(p: GWParams, h: float = 0.01) -> pd.DataFrame:
    """Gas-water drainage table as a function of Sw (gas = non-wetting)."""
    p.validate()
    sw = _grid(p.swl, 1.0, h, [1.0 - p.sgr, 1.0 - p.sgcr])
    sg = 1.0 - sw
    xw = (sw - p.swl) / (1.0 - p.swl - p.sgr)
    xg = (sg - p.sgcr) / (1.0 - p.swl - p.sgcr)
    return pd.DataFrame({"SW": sw, "KRW": p.krwend * p.krw.f(xw), "KRG": p.krgend * p.krg.f(xg)})


def three_phase_kro(wo: pd.DataFrame, go: pd.DataFrame, swl: float, sw, sg, model: str = "stone2"):
    """Oil relative permeability in a three-phase system (arrays sw, sg; so = 1-sw-sg >= 0)."""
    sw, sg = np.asarray(sw, float), np.asarray(sg, float)
    krow = np.interp(sw, wo.SW, wo.KROW)
    krw = np.interp(sw, wo.SW, wo.KRW)
    krog = np.interp(sg, go.SG, go.KROG)
    krg = np.interp(sg, go.SG, go.KRG)
    kroend = float(wo.KROW.max())
    if model == "stone2":
        kro = kroend * ((krow / kroend + krw) * (krog / kroend + krg) - (krw + krg))
    elif model == "stone1":
        sorw = float(1.0 - wo.SW[wo.KROW > 1e-9].max()) if (wo.KROW > 1e-9).any() else 0.0
        sorg = float(1.0 - swl - go.SG[go.KROG > 1e-9].max()) if (go.KROG > 1e-9).any() else 0.0
        frac = np.clip(sg / max(1.0 - swl - sorg, 1e-9), 0, 1)
        som = (1 - frac) * sorw + frac * sorg
        so = 1.0 - sw - sg
        d = np.maximum(1.0 - swl - som, 1e-9)
        son, swn, sgn = np.clip((so - som) / d, 0, 1), np.clip((sw - swl) / d, 0, 1), np.clip(sg / d, 0, 1)
        a = np.maximum((1 - swn) * (1 - sgn), 1e-9)
        kro = son * krow * krog / (kroend * a)
    else:  # baker (saturation-weighted linear)
        wgt_w, wgt_g = np.maximum(sw - swl, 0), np.maximum(sg, 0)
        den = wgt_w + wgt_g
        kro = np.where(den > 1e-12, (wgt_w * krow + wgt_g * krog) / np.where(den > 1e-12, den, 1), krow)
    so = 1.0 - sw - sg
    return np.where(so > 1e-12, np.clip(kro, 0.0, kroend), 0.0)


def sof3_table(wo: pd.DataFrame, go: pd.DataFrame, swl: float) -> pd.DataFrame:
    """Eclipse SOF3 table (So, krow, krog) built from the WO and GO tables."""
    so = np.unique(np.round(np.concatenate([1.0 - wo.SW.values, 1.0 - swl - go.SG.values]), 10))
    so = so[(so >= 0) & (so <= 1.0 - swl + 1e-12)]
    krow = np.interp(1.0 - so, wo.SW, wo.KROW)
    krog = np.interp(1.0 - swl - so, go.SG, go.KROG)
    return pd.DataFrame({"SO": so, "KROW": krow, "KROG": krog})


# --------------------------------------------------------------------------- #
# 2. Capillary pressure  (all pressures in bar)
# --------------------------------------------------------------------------- #
PC_MODELS = ["brooks-corey", "let", "van-genuchten", "skjaeveland", "let-two-sided"]


@dataclass
class PcParams:
    kind: str = "brooks-corey"
    pd: float = 0.30  # entry / displacement pressure [bar]
    pcmax: float = 5.0  # pressure at irreducible saturation [bar]
    pcmin: float = 1.0  # |most negative Pc| for the two-sided model [bar]
    lam: float = 2.0  # Brooks-Corey pore-size index
    L: float = 1.5
    E: float = 4.0
    T: float = 1.5
    Ln: float = 1.5
    En: float = 4.0
    Tn: float = 1.5
    vg_n: float = 2.5
    # Skjaeveland et al. (1998): Pc = cw / Sw*^aw - co / So*^ao   (cw, co >= 0)
    cw: float = 0.30
    aw: float = 0.50
    co: float = 0.15
    ao: float = 1.00
    smooth_cap: bool = True  # Brooks-Corey / van Genuchten: reach exactly pcmax at Swirr (no kink)

    def validate(self):
        if self.pcmax <= self.pd and self.kind != "let-two-sided":
            raise ValueError("Pcmax must exceed the entry pressure")


def pc_curve(sw, swirr: float, sorw: float, p: PcParams):
    """Pc(Sw) for the wetting-phase saturation `sw`."""
    sw = np.asarray(sw, float)
    if p.kind == "let-two-sided":
        s = _clip01((sw - swirr) / (1.0 - sorw - swirr))
        om = 1.0 - s
        pos = p.pcmax * om ** p.L / (om ** p.L + p.E * s ** p.T + 1e-300)
        neg = p.pcmin * s ** p.Ln / (s ** p.Ln + p.En * om ** p.Tn + 1e-300)
        return pos - neg
    if p.kind == "skjaeveland":
        sws = _clip01((sw - swirr) / (1.0 - swirr))
        sos = _clip01((1.0 - sw - sorw) / (1.0 - sorw))
        pos = p.cw / np.maximum(sws, 1e-9) ** p.aw
        neg = p.co / np.maximum(sos, 1e-9) ** p.ao
        return np.clip(pos - neg, -p.pcmin, p.pcmax)
    s = _clip01((sw - swirr) / (1.0 - swirr))
    if p.kind == "brooks-corey":
        if p.smooth_cap and p.pcmax > p.pd > 0:
            s_min = (p.pd / p.pcmax) ** p.lam
            s = s_min + (1.0 - s_min) * s
            return p.pd * s ** (-1.0 / p.lam)
        return np.minimum(p.pd * np.maximum(s, 1e-9) ** (-1.0 / p.lam), p.pcmax)
    if p.kind == "let":
        om = 1.0 - s
        return p.pd + (p.pcmax - p.pd) * om ** p.L / (om ** p.L + p.E * s ** p.T + 1e-300)
    if p.kind == "van-genuchten":
        m = 1.0 - 1.0 / p.vg_n
        if p.smooth_cap and p.pd > 0 and p.pcmax > 0:
            s_min = (1.0 + (p.pcmax / p.pd) ** p.vg_n) ** (-m)
            s = s_min + (1.0 - s_min) * s
            return p.pd * (s ** (-1.0 / m) - 1.0) ** (1.0 / p.vg_n)
        return np.minimum(p.pd * (np.maximum(s, 1e-9) ** (-1.0 / m) - 1.0) ** (1.0 / p.vg_n), p.pcmax)
    raise ValueError(f"unknown Pc model {p.kind}")


# Typical interfacial tension (mN/m) and contact angle (deg). Illustrative defaults - edit per study.
FLUID_SYSTEMS = {
    "Oil-water (reservoir)": dict(sigma=30.0, theta=0.0),
    "Gas-oil (reservoir)": dict(sigma=5.0, theta=0.0),
    "Gas-water (reservoir)": dict(sigma=50.0, theta=0.0),
    "CO2-brine (storage)": dict(sigma=30.0, theta=0.0),
    "H2-brine (storage)": dict(sigma=40.0, theta=0.0),
    "Air-brine (lab)": dict(sigma=72.0, theta=0.0),
    "Oil-brine (lab)": dict(sigma=48.0, theta=0.0),
    "Mercury-air (lab, MICP)": dict(sigma=480.0, theta=140.0),
}


def sigma_cos(system: dict) -> float:
    return system["sigma"] * np.cos(np.radians(system["theta"]))


def pc_convert(pc, sys_from: dict, sys_to: dict, k_from=None, phi_from=None, k_to=None, phi_to=None):
    """Leverett scaling of Pc between fluid systems and (optionally) rock properties."""
    # magnitudes: mercury-air injection pressure is positive although cos(140 deg) < 0
    f = abs(sigma_cos(sys_to)) / abs(sigma_cos(sys_from))
    if None not in (k_from, phi_from, k_to, phi_to):
        f *= np.sqrt((phi_to / k_to) / (phi_from / k_from))
    return np.asarray(pc) * f


def j_function(pc_bar, sys: dict, k_md: float, phi: float):
    pa = np.asarray(pc_bar) * 1e5
    return pa * np.sqrt(k_md * MD_TO_M2 / phi) / (abs(sigma_cos(sys)) * 1e-3)


def pc_from_j(j, sys: dict, k_md: float, phi: float):
    return np.asarray(j) * (abs(sigma_cos(sys)) * 1e-3) / np.sqrt(k_md * MD_TO_M2 / phi) / 1e5


def height_above_fwl(pc_bar, drho_kgm3: float):
    """Height above free-water level [m] for density difference drho [kg/m3]."""
    return np.asarray(pc_bar) * 1e5 / (drho_kgm3 * 9.81)


def wettability_indices(sw, pc):
    """Indicative indices from a two-sided Pc curve: Amott-type (spontaneous imbibition share)
    and USBM = log10(A1/A2)."""
    sw, pc = np.asarray(sw, float), np.asarray(pc, float)
    pos, neg = np.clip(pc, 0, None), np.clip(-pc, 0, None)
    a1, a2 = np.trapezoid(pos, sw), np.trapezoid(neg, sw)
    usbm = float(np.log10(a1 / a2)) if a1 > 0 and a2 > 0 else float("nan")
    if (pc < 0).any() and (pc > 0).any():
        i = np.where(np.diff(np.sign(pc)) != 0)[0][0]
        sw0 = sw[i] + (0 - pc[i]) * (sw[i + 1] - sw[i]) / (pc[i + 1] - pc[i])
        amott = float((sw0 - sw[0]) / (sw[-1] - sw[0]))
    else:
        amott = float("nan")
    return {"usbm": usbm, "amott": amott}


def add_pc_wo(df: pd.DataFrame, p: WOParams, pc: PcParams, scale: float = 1.0) -> pd.DataFrame:
    out = df.copy()
    out["PC"] = scale * pc_curve(out.SW.values, p.swl, p.sorw, pc)
    return out


def add_pc_go(df: pd.DataFrame, p: GOParams, pc: PcParams, scale: float = 1.0) -> pd.DataFrame:
    """Pcgo = Pg - Po, increasing with Sg; computed on the liquid saturation."""
    out = df.copy()
    sg = out.SG.values
    pcv = pc_curve(1.0 - sg, p.swl, 0.0, pc)
    if p.sgcr > 0:  # ramp from 0 to the model value at Sgcr: avoids a jump at the first gas saturation
        pcv = pcv * np.minimum(sg / p.sgcr, 1.0)
    pcv[sg <= 0] = 0.0
    out["PC"] = scale * pcv
    return out


def add_pc_gw(df: pd.DataFrame, p: GWParams, pc: PcParams, scale: float = 1.0) -> pd.DataFrame:
    out = df.copy()
    out["PC"] = scale * pc_curve(out.SW.values, p.swl, p.sgr, pc)
    return out


def sw_from_pc(pc_bar, swirr: float, sorw: float, p: PcParams):
    """Inverse of the (monotone decreasing) Pc curve: wetting saturation at a given Pc."""
    sw = np.linspace(swirr, 1.0 - (sorw if p.kind in ("let-two-sided", "skjaeveland") else 0.0), 4000)
    pc = pc_curve(sw, swirr, sorw, p)
    pc = np.minimum.accumulate(pc)  # enforce monotone for np.interp
    return np.interp(-np.asarray(pc_bar, float), -pc, sw)


def saturation_height(heights_m, drho_kgm3: float, swirr: float, sorw: float, p: PcParams, sigma_scale: float = 1.0):
    """Saturation-height function: Sw at heights above the free-water level.
    `sigma_scale` converts the curve's fluid system to the in-situ one (see pc_convert)."""
    pc = np.asarray(heights_m, float) * drho_kgm3 * 9.81 / 1e5
    return sw_from_pc(pc / sigma_scale, swirr, sorw, p)


def seal_column_height(entry_pressure_bar: float, drho_kgm3: float) -> float:
    """Maximum static column [m] a seal with the given entry pressure can hold (capillary seal capacity)."""
    return float(entry_pressure_bar * 1e5 / (drho_kgm3 * 9.81))


def micp_to_pc(p_hg_bar, s_hg, target: dict, sconf: float = 0.0, k_from=None, phi_from=None, k_to=None, phi_to=None):
    """Convert mercury-injection data (Hg pressure, Hg saturation) to Pc vs wetting saturation in `target`.
    `sconf` is the closure / conformance correction (fraction of pore volume filled by surface conformance)."""
    s_hg = np.asarray(s_hg, float)
    s_corr = np.clip((s_hg - sconf) / (1.0 - sconf), 0.0, 1.0)
    sw = 1.0 - s_corr
    pc = pc_convert(p_hg_bar, FLUID_SYSTEMS["Mercury-air (lab, MICP)"], target, k_from, phi_from, k_to, phi_to)
    order = np.argsort(sw)
    return sw[order], np.asarray(pc, float)[order]


def demo_micp(seed: int = 11) -> pd.DataFrame:
    """Synthetic mercury-injection capillary pressure test (Hg pressure [bar], Hg saturation) with noise."""
    rng = np.random.default_rng(seed)
    p_hg = np.geomspace(0.7, 120.0, 26)
    swirr, pd_hg, lam = 0.08, 3.2, 1.7
    s_eff = np.clip((pd_hg / p_hg) ** lam, 0, 1)
    sw = swirr + (1 - swirr) * s_eff
    shg = (1 - sw) + 0.0012 * rng.standard_normal(len(sw))
    shg[p_hg < pd_hg] = 0.012 * rng.random((p_hg < pd_hg).sum())  # surface conformance below entry
    shg = np.maximum.accumulate(np.clip(shg, 0, 0.95))  # injection data are cumulative: monotone
    return pd.DataFrame({"P_HG_BAR": p_hg, "S_HG": shg})


def pc_checks(sw, pc, p: PcParams | None = None) -> list[str]:
    """Plain-language sanity checks on a Pc table. Empty list = no issues."""
    sw, pc, out = np.asarray(sw, float), np.asarray(pc, float), []
    if np.any(np.diff(pc) > 1e-9):
        out.append("Pc is not monotonically decreasing with Sw - simulators may reject or smooth the table.")
    if p is not None and p.kind not in ("let-two-sided", "skjaeveland") and p.pcmax <= p.pd * 1.05:
        out.append("Pc max is close to the entry pressure: the transition zone will be almost flat.")
    if np.nanmax(pc) / max(np.nanmin(np.abs(pc[np.abs(pc) > 0])) if np.any(np.abs(pc) > 0) else 1, 1e-9) > 1e4:
        out.append("Pc spans more than four decades - check units (bar) and the cap value.")
    if p is not None and p.kind == "let-two-sided" and not (pc < 0).any():
        out.append("Two-sided model selected but the curve never goes negative: increase Pc min.")
    return out


_FIT_SPECS = {  # kind: (start, lower, upper) as functions of pmax; swirr is appended separately
    "brooks-corey": lambda pm, p0: ([max(p0, 0.01), 2.0], [1e-4, 0.2], [pm, 12.0]),
    "van-genuchten": lambda pm, p0: ([max(p0, 0.05), 2.5], [1e-4, 1.2], [pm, 12.0]),
    "let": lambda pm, p0: ([max(p0, 0.01), pm, 1.5, 3.0, 1.5], [1e-4, 0.8 * pm, 0.2, 0.05, 0.2], [pm, 30 * pm, 10.0, 40.0, 10.0]),
    "skjaeveland": lambda pm, p0: ([0.06 * pm, 0.6, 0.02 * pm, 1.0], [1e-4, 0.1, 0.0, 0.1], [4 * pm, 4.0, 4 * pm, 4.0]),
    "let-two-sided": lambda pm, p0: ([pm, 0.5, 1.5, 3.0, 1.5, 1.5, 3.0, 1.5],
                                     [0.1, 0.0, 0.2, 0.05, 0.2, 0.2, 0.05, 0.2], [30 * pm, 30 * pm, 10, 40, 10, 10, 40, 10]),
}


def _pc_from_vector(kind: str, th, base: PcParams) -> PcParams:
    q = replace(base, kind=kind)
    if kind == "brooks-corey":
        q.pd, q.lam = th
    elif kind == "van-genuchten":
        q.pd, q.vg_n = th
    elif kind == "let":
        q.pd, q.pcmax, q.L, q.E, q.T = th
    elif kind == "skjaeveland":
        q.cw, q.aw, q.co, q.ao = th
    else:
        q.pcmax, q.pcmin, q.L, q.E, q.T, q.Ln, q.En, q.Tn = th
    return q


def fit_pc(sw, pc, swirr: float, sorw: float, kind: str = "brooks-corey", fit_swirr: bool = True):
    """Least-squares fit of a Pc model to measured (Sw, Pc) points.

    - Swirr is fitted too (default): MICP and centrifuge data only approach it asymptotically, and fixing it
      at the lowest measured Sw biases the pore-size index. Bounds: 0 <= Swirr < min(Sw).
    - Brooks-Corey / van Genuchten / Skjaeveland are fitted as pure functions (no cap), with the cap pushed
      far beyond the data, which is the standard way to fit MICP data.
    - Residuals are scaled by max(|Pc|, 5 % of max|Pc|) so low-pressure points are not drowned by the steep
      end; soft-L1 loss for robustness.
    Returns (PcParams, rmse_bar, r2, swirr_fitted).
    """
    from scipy.optimize import least_squares

    sw, pc = np.asarray(sw, float), np.asarray(pc, float)
    pmax = float(np.max(np.abs(pc)))
    scale = np.maximum(np.abs(pc), 0.05 * pmax)
    far = 30.0 * pmax
    base = PcParams(kind=kind, pcmax=far, pcmin=far, smooth_cap=False)
    x0, lo, hi = _FIT_SPECS[kind](pmax, float(pc[pc > 0].min()) if (pc > 0).any() else 0.1)
    sw_hi = float(sw.min()) - 1e-3
    if fit_swirr:
        x0, lo, hi = x0 + [max(min(swirr, sw_hi) - 0.02, 0.0)], lo + [0.0], hi + [max(sw_hi, 1e-3)]
    else:
        x0, lo, hi = x0, lo, hi
    n = len(x0) - (1 if fit_swirr else 0)

    def unpack(th):
        return _pc_from_vector(kind, th[:n], base), (float(th[-1]) if fit_swirr else swirr)

    def resid(th):
        q, s0 = unpack(th)
        return (pc_curve(sw, s0, sorw, q) - pc) / scale

    r = least_squares(resid, x0, bounds=(lo, hi), loss="soft_l1", max_nfev=600)
    fitted, s_fit = unpack(r.x)
    if kind in ("let", "let-two-sided"):
        pass  # pcmax is a fitted parameter
    model = pc_curve(sw, s_fit, sorw, fitted)
    rmse = float(np.sqrt(np.mean((model - pc) ** 2)))
    ss_res, ss_tot = float(np.sum((model - pc) ** 2)), float(np.sum((pc - pc.mean()) ** 2))
    return fitted, rmse, ((1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan")), float(s_fit)


# --------------------------------------------------------------------------- #
# 3. Endpoint scaling
# --------------------------------------------------------------------------- #
def rescale_wo_table(
    df: pd.DataFrame,
    old: WOParams,
    new: WOParams,
    mode: str = "two-point",
    krwr_new: float | None = None,
    krwmax_new: float = 1.0,
    pc_factor: float = 1.0,
    h: float = 0.01,
) -> pd.DataFrame:
    """Re-scale a water-oil table from `old` to `new` end-points (Eclipse ENDSCALE logic).

    two-point : Swcr/Sorw shift, krw and kro scaled by end-point ratio.
    three-point: additionally pins krw at Sorw (krwr_new) and lets krw rise to krwmax_new at Sw=1.
    """
    sw = _grid(new.swl, 1.0, h, [new.swcr, 1.0 - new.sorw])
    xw = _clip01((sw - new.swcr) / (1.0 - new.swcr - new.sorw))
    xo = _clip01((1.0 - sw - new.sorw) / (1.0 - new.sorw - new.swl))
    sw_old_w = old.swcr + xw * (1.0 - old.swcr - old.sorw)
    sw_old_o = 1.0 - old.sorw - xo * (1.0 - old.sorw - old.swl)
    krw_o = np.interp(sw_old_w, df.SW, df.KRW)
    kro_o = np.interp(sw_old_o, df.SW, df.KROW)
    kro = kro_o * new.kroend / old.kroend
    if mode == "three-point":
        krwr_old = float(np.interp(1.0 - old.sorw, df.SW, df.KRW))
        krwr_new = krwr_new if krwr_new is not None else new.krwend
        krw = krw_o * krwr_new / krwr_old
        above = sw > 1.0 - new.sorw
        frac = (sw[above] - (1.0 - new.sorw)) / max(new.sorw, 1e-9)
        krw[above] = krwr_new + frac * (krwmax_new - krwr_new)
    else:
        krw = krw_o * new.krwend / old.krwend
    out = pd.DataFrame({"SW": sw, "KRW": krw, "KROW": kro})
    if "PC" in df:
        s_old = old.swl + (sw - new.swl) / (1.0 - new.swl) * (1.0 - old.swl)
        out["PC"] = pc_factor * np.interp(s_old, df.SW, df.PC)
    return out


def rock_type_family(low: WOParams, high: WOParams, n_types: int, k_low: float, k_high: float,
                     phi_low: float, phi_high: float, ref_pc: PcParams | None = None):
    """Interpolate end-points and exponents in ln(k) between a poor and a good rock type.

    Returns a list of dicts: {id, k, phi, params, pc_scale}. pc_scale follows Leverett sqrt(phi/k)
    relative to the geometric-mean rock.
    """
    ks = np.geomspace(k_low, k_high, n_types)
    t = (np.log(ks) - np.log(k_low)) / (np.log(k_high) - np.log(k_low)) if n_types > 1 else np.array([0.5])
    phis = phi_low + t * (phi_high - phi_low)
    k_ref, phi_ref = np.sqrt(k_low * k_high), 0.5 * (phi_low + phi_high)
    lerp = lambda a, b, x: a + (b - a) * x  # noqa: E731
    fam = []
    for i, (k, phi, ti) in enumerate(zip(ks, phis, t), start=1):
        mk = lambda a, b: KrModel(a.kind, n=lerp(a.n, b.n, ti), L=lerp(a.L, b.L, ti),  # noqa: E731
                                  E=lerp(a.E, b.E, ti), T=lerp(a.T, b.T, ti))
        swl = lerp(low.swl, high.swl, ti)
        p = WOParams(swl=swl, swcr=max(lerp(low.swcr, high.swcr, ti), swl), sorw=lerp(low.sorw, high.sorw, ti),
                     krwend=lerp(low.krwend, high.krwend, ti), kroend=lerp(low.kroend, high.kroend, ti),
                     krw=mk(low.krw, high.krw), kro=mk(low.kro, high.kro))
        fam.append(dict(id=i, k=float(k), phi=float(phi), params=p,
                        pc_scale=float(np.sqrt((phi / k) / (phi_ref / k_ref)))))
    return fam


# --------------------------------------------------------------------------- #
# 4. Hysteresis (non-wetting phase kr + Pc scanning)
# --------------------------------------------------------------------------- #
def land_trapped(s_hy: float, s_crit: float, s_trap_max: float, s_max: float) -> float:
    """Land (1968) trapped non-wetting saturation for a reversal at s_hy."""
    if s_max <= s_crit + 1e-12 or s_trap_max <= s_crit + 1e-12:
        return s_crit
    c = 1.0 / (s_trap_max - s_crit) - 1.0 / (s_max - s_crit)
    d = max(s_hy - s_crit, 0.0)
    return s_crit + d / (1.0 + c * d)


@dataclass
class HystSetup:
    s: np.ndarray  # non-wetting saturation, ascending (primary drainage)
    kr: np.ndarray  # drainage kr
    s_crit: float  # critical non-wetting saturation (drainage)
    s_max: float  # maximum non-wetting saturation reached
    s_trap_max: float  # maximum trapped non-wetting saturation (boundary imbibition end)

    def kr_d(self, s):
        return np.interp(s, self.s, self.kr)


def scanning_curve(h: HystSetup, s_hy: float, scheme: str = "killough", n: int = 60) -> pd.DataFrame:
    """Imbibition scanning curve starting at (s_hy, kr_d(s_hy)). Columns: S, KR."""
    s_hy = float(np.clip(s_hy, h.s_crit + 1e-6, h.s_max))
    if scheme == "killough":
        s_trap = land_trapped(s_hy, h.s_crit, h.s_trap_max, h.s_max)
        s = np.linspace(s_trap, s_hy, n)
        sstar = h.s_crit + (s - s_trap) * (h.s_max - h.s_crit) / max(s_hy - s_trap, 1e-9)
        kr = h.kr_d(sstar) * h.kr_d(s_hy) / max(float(h.kr_d(h.s_max)), 1e-12)
        return pd.DataFrame({"S": s, "KR": kr})
    # Carlson: boundary imbibition curve translated parallel to pass through the reversal point
    span = (h.s_max - h.s_crit) / max(h.s_max - h.s_trap_max, 1e-9)
    s_b_hy = h.s_trap_max + (s_hy - h.s_crit) / span
    delta = s_b_hy - s_hy
    s_end = h.s_trap_max - delta
    s = np.linspace(s_end, s_hy, n)
    sstar = h.s_crit + (s + delta - h.s_trap_max) * span
    return pd.DataFrame({"S": s, "KR": h.kr_d(sstar)})


def hysteresis_family(h: HystSetup, scheme: str, n_scans: int = 5):
    """Drainage curve, boundary imbibition curve and `n_scans` scanning curves."""
    drain = pd.DataFrame({"S": h.s, "KR": h.kr})
    bound = scanning_curve(h, h.s_max, scheme)
    hys = np.linspace(h.s_crit + 0.2 * (h.s_max - h.s_crit), h.s_max, n_scans + 1)[:-1]
    scans = []
    for s_hy in hys:
        sc = scanning_curve(h, float(s_hy), scheme)
        scans.append(dict(s_hy=float(s_hy), s_trap=float(sc.S.iloc[0]), curve=sc))
    return dict(drainage=drain, boundary=bound, scans=scans)


def pc_scanning(sw, pc_dr, pc_imb, sw_hy: float, sw_end: float, direction: str = "up",
                curvature: float = 0.1):
    """Killough-type Pc scanning curve. `pc_dr`/`pc_imb` are the boundary curves evaluated on `sw`.

    direction "up"  : water imbibition from the drainage curve at sw_hy towards the imbibition curve.
    direction "down": drainage from the imbibition curve at sw_hy towards the drainage curve.
    """
    sw = np.asarray(sw, float)
    eps = curvature
    x = (sw - sw_hy) if direction == "up" else (sw_hy - sw)
    x_end = abs(sw_end - sw_hy)
    f = (1.0 / (np.maximum(x, 0) + eps) - 1.0 / eps) / (1.0 / (x_end + eps) - 1.0 / eps)
    f = np.clip(f, 0.0, 1.0)
    if direction == "up":
        return pc_dr + (pc_imb - pc_dr) * f
    return pc_imb + (pc_dr - pc_imb) * f


def imbibition_gas_table(drain: pd.DataFrame, h: HystSetup, scheme: str = "killough") -> pd.DataFrame:
    """Imbibition SGOF-style table (SG, KRG_IMB) on the drainage Sg grid from the boundary curve."""
    b = scanning_curve(h, h.s_max, scheme, n=200)
    sg = drain.SG.values
    kr = np.where(sg <= b.S.iloc[0], 0.0, np.interp(sg, b.S, b.KR))
    out = drain.copy()
    out["KRG"] = kr
    return out


# --------------------------------------------------------------------------- #
# 5. Displacement analytics
# --------------------------------------------------------------------------- #
def fractional_flow(krw, kro, mu_w: float, mu_o: float):
    krw, kro = np.asarray(krw, float), np.asarray(kro, float)
    den = krw / mu_w + kro / mu_o
    return np.divide(krw / mu_w, den, out=np.zeros_like(den), where=den > 0)


def welge(df: pd.DataFrame, mu_w: float, mu_o: float, swi: float | None = None, n: int = 2000) -> dict:
    """Buckley-Leverett / Welge analytics on a WO table (SW, KRW, KROW)."""
    swi = float(df.SW.iloc[0]) if swi is None else swi
    smax = float(df.SW[df.KROW > 1e-9].max()) if (df.KROW > 1e-9).any() else float(df.SW.iloc[-1])
    sw = np.linspace(swi, smax, n)
    krw, kro = np.interp(sw, df.SW, df.KRW), np.interp(sw, df.SW, df.KROW)
    fw = fractional_flow(krw, kro, mu_w, mu_o)
    dfw = np.gradient(fw, sw)
    sl = np.where(sw > swi + 1e-6, (fw - fw[0]) / np.maximum(sw - swi, 1e-12), 0.0)
    i = int(np.argmax(sl))
    slope = float(sl[i])
    pvi_bt = (1.0 / slope) if slope > 0 else float("inf")
    cross = np.where(np.diff(np.sign(df.KRW.values - df.KROW.values)) != 0)[0]
    sw_cross = float("nan")
    if len(cross):
        j = cross[0]
        a = (df.KRW.iloc[j] - df.KROW.iloc[j]) - (df.KRW.iloc[j + 1] - df.KROW.iloc[j + 1])
        sw_cross = float(df.SW.iloc[j] + (df.SW.iloc[j + 1] - df.SW.iloc[j]) * (df.KRW.iloc[j] - df.KROW.iloc[j]) / a)
    m_end = (float(df.KRW.max()) / mu_w) / (float(df.KROW.max()) / mu_o)
    return dict(sw=sw, fw=fw, dfw=dfw, swi=swi, swf=float(sw[i]), fwf=float(fw[i]), slope=slope,
                pvi_bt=pvi_bt, rf_bt=pvi_bt / (1 - swi), rf_ult=(smax - swi) / (1 - swi),
                mobility_ratio=m_end, sw_cross=sw_cross, smax=smax, idx_f=i)


def oil_produced_pv(w: dict, t: float) -> float:
    """Oil produced [PV] after `t` pore volumes of water injected (Welge, no capillary effects)."""
    if t <= w["pvi_bt"]:
        return float(t)
    i0 = w["idx_f"]
    d_rev = np.maximum.accumulate(w["dfw"][i0:][::-1])
    s_rev = w["sw"][i0:][::-1]
    sw2 = float(np.interp(1.0 / t, d_rev, s_rev))
    fw2 = float(np.interp(sw2, w["sw"], w["fw"]))
    return float(min(sw2 + (1.0 - fw2) * t - w["swi"], w["smax"] - w["swi"]))


def recovery_curve(w: dict, pvi_max: float = 3.0, n: int = 300):
    """Oil recovery factor vs pore volumes injected (Welge, no capillary effects)."""
    t = np.linspace(0, pvi_max, n)
    rf = np.array([oil_produced_pv(w, tk) / (1.0 - w["swi"]) for tk in t])
    return pd.DataFrame({"PVI": t, "RF": rf})


def welge_profile(w: dict, xd, t: float):
    """Water saturation along the core (xd in [0,1]) after t PV injected."""
    xd = np.asarray(xd, float)
    if t <= 0:
        return np.full_like(xd, w["swi"])
    i0 = w["idx_f"]
    d_rev = np.maximum.accumulate(w["dfw"][i0:][::-1])
    s_rev = w["sw"][i0:][::-1]
    v = xd / t
    return np.where(v > w["slope"], w["swi"], np.interp(v, d_rev, s_rev))


# --------------------------------------------------------------------------- #
# 6. Fitting measured SCAL data
# --------------------------------------------------------------------------- #
def demo_scal_data(seed: int = 7) -> pd.DataFrame:
    """Synthetic 'lab' water-oil data with noise, for trying out the fitting workflow."""
    rng = np.random.default_rng(seed)
    p = WOParams(swl=0.14, swcr=0.18, sorw=0.27, krwend=0.32, kroend=0.85,
                 krw=KrModel("let", L=2.2, E=1.6, T=2.5), kro=KrModel("let", L=2.8, E=2.4, T=1.4))
    d = wateroil(p, 0.04)
    d["KRW"] *= 1 + 0.04 * rng.standard_normal(len(d))
    d["KROW"] *= 1 + 0.04 * rng.standard_normal(len(d))
    d["KRW"] = d.KRW.clip(lower=0)
    d["KROW"] = d.KROW.clip(lower=0)
    d["PC"] = pc_curve(d.SW.values, 0.14, 0.27, PcParams("let-two-sided", pcmax=3.0, pcmin=1.2)) * (
        1 + 0.03 * rng.standard_normal(len(d)))
    return d.rename(columns=str.upper)


def detect_endpoints(df: pd.DataFrame, tol: float = 1e-3) -> WOParams:
    sw = df.SW.values
    swl = float(sw.min())
    swcr = float(sw[df.KRW.values <= tol].max()) if (df.KRW.values <= tol).any() else swl
    sor_pts = sw[df.KROW.values > tol]
    sorw = float(1.0 - sor_pts.max()) if len(sor_pts) else 0.2
    krwend = float(df.KRW.max())
    kroend = float(df.KROW.max())
    return WOParams(swl=swl, swcr=max(swcr, swl), sorw=sorw, krwend=min(krwend, 1.0), kroend=min(kroend, 1.0))


def fit_relperm(df: pd.DataFrame, kind: str = "corey", ep: WOParams | None = None):
    """Least-squares fit of Corey exponents or LET parameters; end-points held fixed."""
    from scipy.optimize import least_squares

    ep = ep or detect_endpoints(df)
    sw = df.SW.values

    def build(theta):
        if kind == "corey":
            kw, ko = KrModel("corey", n=theta[0]), KrModel("corey", n=theta[1])
        else:
            kw = KrModel("let", L=theta[0], E=theta[1], T=theta[2])
            ko = KrModel("let", L=theta[3], E=theta[4], T=theta[5])
        return replace(ep, krw=kw, kro=ko)

    def resid(theta):
        t = wateroil(build(theta), 0.01)
        rw = np.interp(sw, t.SW, t.KRW) - df.KRW.values
        ro = np.interp(sw, t.SW, t.KROW) - df.KROW.values
        return np.concatenate([rw, ro])

    if kind == "corey":
        x0, lo, hi = [2.5, 2.5], [0.5, 0.5], [8, 8]
    else:
        x0, lo, hi = [2, 2, 2, 2, 2, 2], [0.3, 0.05, 0.3] * 2, [8, 30, 8] * 2
    r = least_squares(resid, x0, bounds=(lo, hi))
    fitted = build(r.x)
    rmse = float(np.sqrt(np.mean(r.fun ** 2)))
    return fitted, rmse


# --------------------------------------------------------------------------- #
# 7. Uncertainty ensembles (water-oil)
# --------------------------------------------------------------------------- #
@dataclass
class EnsembleSpec:
    n: int = 200
    seed: int = 7
    dist: str = "triangular"  # "triangular" (mode = base case) | "uniform"
    d_swl: float = 0.03  # +/- absolute
    d_sorw: float = 0.05  # +/- absolute
    f_krwend: float = 0.20  # +/- relative
    f_kroend: float = 0.10  # +/- relative
    d_nw: float = 0.5  # Corey exponent +/- absolute (LET: L and T scaled by up to +/- d/3)
    d_no: float = 0.5
    f_pc: float = 0.25  # sigma of ln(Pc multiplier)


def shift_kr(m: KrModel, d: float) -> KrModel:
    """Change curve shape by `d`: Corey exponent += d; LET: L and T scaled by up to +/- d/3 (clipped to 0.5..1.5)."""
    if m.kind == "corey":
        return KrModel("corey", n=max(0.5, m.n + d))
    f = float(np.clip(1.0 + d / 3.0, 0.5, 1.5))
    return KrModel("let", L=m.L * f, E=m.E, T=m.T * f)


def _draw(rng, dist: str, half: float, n: int):
    if half <= 0:
        return np.zeros(n)
    if dist == "uniform":
        return rng.uniform(-half, half, n)
    return rng.triangular(-half, 0.0, half, n)


def ensemble_wo(base: WOParams, spec: EnsembleSpec, mu_w: float, mu_o: float, h: float = 0.01,
                pc: PcParams | None = None) -> dict:
    """Monte-Carlo ensemble of water-oil curves around `base`.

    Returns grid, krw/kro/pc matrices (n x grid), a samples DataFrame (inputs + displacement metrics) and
    the list of WOParams. Percentiles are point-wise: P10 = low, P90 = high.
    """
    rng = np.random.default_rng(spec.seed)
    n = int(spec.n)
    dswl, dsor = _draw(rng, spec.dist, spec.d_swl, n), _draw(rng, spec.dist, spec.d_sorw, n)
    fkw, fko = _draw(rng, spec.dist, spec.f_krwend, n), _draw(rng, spec.dist, spec.f_kroend, n)
    dnw, dno = _draw(rng, spec.dist, spec.d_nw, n), _draw(rng, spec.dist, spec.d_no, n)
    pcs = np.exp(spec.f_pc * rng.standard_normal(n)) if spec.f_pc > 0 else np.ones(n)
    grid = np.linspace(max(base.swl - spec.d_swl, 0.0), 1.0, 161)
    krw, kro = np.zeros((n, grid.size)), np.zeros((n, grid.size))
    pcm = np.zeros((n, grid.size)) if pc is not None else None
    rows, plist = [], []

    shape = shift_kr

    for i in range(n):
        swl = float(np.clip(base.swl + dswl[i], 0.0, 0.6))
        swcr = max(base.swcr + dswl[i], swl)
        sorw = float(np.clip(base.sorw + dsor[i], 0.01, 0.95 - swcr))
        p = WOParams(swl=swl, swcr=swcr, sorw=sorw,
                     krwend=float(np.clip(base.krwend * (1 + fkw[i]), 0.02, 1.0)),
                     kroend=float(np.clip(base.kroend * (1 + fko[i]), 0.1, 1.0)),
                     krw=shape(base.krw, dnw[i]), kro=shape(base.kro, dno[i]))
        t = wateroil(p, h)
        krw[i], kro[i] = np.interp(grid, t.SW, t.KRW), np.interp(grid, t.SW, t.KROW)
        if pc is not None:
            pcm[i] = pcs[i] * pc_curve(grid, swl, sorw, pc)
        w = welge(t, mu_w, mu_o, swi=swl)
        rows.append(dict(swl=swl, swcr=swcr, sorw=sorw, krwend=p.krwend, kroend=p.kroend, pc_mult=float(pcs[i]),
                         rf_ult=w["rf_ult"], pvi_bt=w["pvi_bt"], rf_bt=w["rf_bt"], mobility_ratio=w["mobility_ratio"],
                         sw_cross=w["sw_cross"]))
        plist.append(p)
    return dict(grid=grid, krw=krw, kro=kro, pc=pcm, samples=pd.DataFrame(rows), params=plist, pc_mult=pcs)


def percentile_bands(mat, ps=(10, 50, 90)) -> dict:
    return {p: np.percentile(mat, p, axis=0) for p in ps}


def pick_cases(ens: dict, metric: str = "rf_ult", pcts=(10, 50, 90)) -> dict:
    """Index of the realisation closest to each percentile of `metric` (representative low / mid / high cases)."""
    v = ens["samples"][metric].values
    return {p: int(np.argmin(np.abs(v - np.percentile(v, p)))) for p in pcts}


# --------------------------------------------------------------------------- #
# 8. Interpolation between curve sets (low-salinity, polymer, WAG, capillary desaturation)
# --------------------------------------------------------------------------- #
def interpolation_factor(value: float, lo: float, hi: float, scale: str = "log") -> float:
    """Map a physical variable to t in [0, 1] (0 = set A, 1 = set B), linear or log10 between lo and hi."""
    if scale == "log":
        value, lo, hi = np.log10(max(value, 1e-30)), np.log10(max(lo, 1e-30)), np.log10(max(hi, 1e-30))
    if hi == lo:
        return 0.0
    return float(np.clip((value - lo) / (hi - lo), 0.0, 1.0))


def capillary_number(mu_cp: float, v_m_per_day: float, sigma_mNm: float) -> float:
    """Nc = mu * v / sigma (dimensionless), viscosity [cP], Darcy velocity [m/day], IFT [mN/m]."""
    return (mu_cp * 1e-3) * (v_m_per_day / 86400.0) / (sigma_mNm * 1e-3)


def interpolate_wo(a: WOParams, b: WOParams, t: float) -> WOParams:
    """Linear interpolation of end-points and curve-shape parameters between two curve sets."""
    if a.krw.kind != b.krw.kind or a.kro.kind != b.kro.kind or a.krw.kind != a.kro.kind:
        raise ValueError("Both curve sets must use the same parametrisation (Corey or LET)")
    lerp = lambda x, y: x + (y - x) * t  # noqa: E731
    mk = lambda m1, m2: KrModel(m1.kind, n=lerp(m1.n, m2.n), L=lerp(m1.L, m2.L), E=lerp(m1.E, m2.E), T=lerp(m1.T, m2.T))  # noqa: E731
    swl = lerp(a.swl, b.swl)
    return WOParams(swl=swl, swcr=max(lerp(a.swcr, b.swcr), swl), sorw=lerp(a.sorw, b.sorw), krwend=lerp(a.krwend, b.krwend),
                    kroend=lerp(a.kroend, b.kroend), krw=mk(a.krw, b.krw), kro=mk(a.kro, b.kro))


# --------------------------------------------------------------------------- #
# 9. Lab interpretation: JBN unsteady-state relative permeability, centrifuge Pc
# --------------------------------------------------------------------------- #
def synthetic_core_flood(df: pd.DataFrame, mu_w: float, mu_o: float, swi: float, pvi_max: float = 5.0,
                         n_pts: int = 22, noise: float = 0.004, seed: int = 3) -> pd.DataFrame:
    """Constant-rate waterflood data (PVI, NP [PV], DP [relative]) generated from a known curve set."""
    rng = np.random.default_rng(seed)
    w = welge(df, mu_w, mu_o, swi=swi)
    t = np.concatenate([[0.0], np.geomspace(0.15, pvi_max, n_pts)])
    xd = (np.arange(600) + 0.5) / 600
    lam = lambda s: np.interp(s, df.SW, df.KRW) / mu_w + np.interp(s, df.SW, df.KROW) / mu_o  # noqa: E731
    npv = np.array([oil_produced_pv(w, tk) for tk in t])
    dp = np.array([np.mean(1.0 / lam(welge_profile(w, xd, tk))) for tk in t])
    dp = dp / dp[0]
    npv_n = npv + noise * rng.standard_normal(t.size) * (t > 0)
    dp_n = dp * (1 + 0.01 * rng.standard_normal(t.size) * (t > 0))
    return pd.DataFrame({"PVI": t, "NP": np.maximum.accumulate(npv_n), "DP": dp_n})


def jbn(data: pd.DataFrame, swi: float, mu_w: float, mu_o: float, kro_swi: float = 1.0,
        noise_np: float = 0.004, noise_dp: float = 0.01) -> pd.DataFrame:
    """Johnson-Bossler-Naumann relative permeabilities from a constant-rate unsteady-state waterflood.

    data: PVI (cumulative water injected, PV), NP (oil produced, PV), DP (pressure drop, any consistent units;
    the first row must be PVI = 0 = initial pressure drop). Valid after breakthrough.
    kro_swi scales kro to absolute permeability (effective oil permeability at Swi / absolute).
    Returns SW (outlet saturation), KRW, KROW, FW.
    """
    from scipy.interpolate import UnivariateSpline

    d = data.sort_values("PVI").reset_index(drop=True)
    if d.PVI.iloc[0] > 1e-9:
        raise ValueError("The first row must be PVI = 0 with the initial pressure drop")
    wi, npv, dp = d.PVI.values, d.NP.values, d.DP.values
    m = wi > 0
    if m.sum() < 8:
        raise ValueError("At least 8 post-start data points are needed")
    ir = dp[0] / dp  # relative injectivity
    s_np = UnivariateSpline(wi, npv, k=3, s=len(wi) * noise_np ** 2)
    fo2 = np.clip(s_np.derivative()(wi[m]), 1e-6, 1.0)
    x, y = 1.0 / wi[m], 1.0 / (wi[m] * ir[m])
    o = np.argsort(x)
    s_y = UnivariateSpline(x[o], y[o], k=3, s=float(np.sum((y[o] * noise_dp) ** 2)))
    dydx = np.maximum(s_y.derivative()(x), 1e-9)
    fw2 = 1.0 - fo2
    sw2 = swi + s_np(wi[m]) - wi[m] * fo2
    kro = kro_swi * fo2 / dydx
    krw = kro * (fw2 / fo2) * (mu_w / mu_o)
    out = pd.DataFrame({"SW": sw2, "KRW": krw, "KROW": kro, "FW": fw2, "PVI": wi[m]})
    return out[(out.FW > 0.02) & (out.SW > swi) & (out.KROW > 0)].sort_values("SW").reset_index(drop=True)


def centrifuge_geometry(r2_cm: float, L_cm: float):
    r2 = r2_cm / 100.0
    return r2 - L_cm / 100.0, r2


def centrifuge_pc_inlet(rpm, r2_cm: float, L_cm: float, drho_gcc: float):
    """Pc [bar] at the inlet (inner) face of the core; Pc = 0 at the outlet face."""
    r1, r2 = centrifuge_geometry(r2_cm, L_cm)
    w = 2.0 * np.pi * np.asarray(rpm, float) / 60.0
    return 0.5 * drho_gcc * 1000.0 * w ** 2 * (r2 ** 2 - r1 ** 2) / 1e5


def centrifuge_avg_sw(rpm, r2_cm: float, L_cm: float, drho_gcc: float, swirr: float, sorw: float, p: PcParams,
                      nr: int = 240):
    """Forward model: average Sw after equilibrium at each speed for a given Pc curve (exact radial geometry)."""
    r1, r2 = centrifuge_geometry(r2_cm, L_cm)
    r = np.linspace(r1, r2, nr)
    out = []
    for rp in np.atleast_1d(rpm):
        w = 2.0 * np.pi * rp / 60.0
        pc = 0.5 * drho_gcc * 1000.0 * w ** 2 * (r2 ** 2 - r ** 2) / 1e5
        out.append(float(np.mean(sw_from_pc(pc, swirr, sorw, p))))
    return np.array(out)


def demo_centrifuge(seed: int = 5) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rpm = np.array([300, 500, 800, 1200, 1800, 2500, 3500, 4500], float)
    true = PcParams(kind="brooks-corey", pd=0.5, lam=1.6, pcmax=1e3, smooth_cap=False)
    sw = centrifuge_avg_sw(rpm, 11.0, 5.0, 1.0, 0.12, 0.0, true)
    return pd.DataFrame({"RPM": rpm, "SW_AVG": np.clip(sw + 0.004 * rng.standard_normal(rpm.size), 0, 1)})


def fit_centrifuge(rpm, sw_avg, r2_cm: float, L_cm: float, drho_gcc: float, kind: str = "brooks-corey", sorw: float = 0.0):
    """Invert centrifuge average-saturation data to a Pc curve by fitting a parametric model through the exact
    radial forward model (preferred to the closed-form Hassler-Brunner approximation, which assumes a thin core).
    Returns (PcParams, rmse_in_saturation, swirr_fitted)."""
    from scipy.optimize import least_squares

    rpm, sw_avg = np.asarray(rpm, float), np.asarray(sw_avg, float)
    pmax = float(np.max(centrifuge_pc_inlet(rpm, r2_cm, L_cm, drho_gcc)))
    far = 30.0 * pmax
    base = PcParams(kind=kind, pcmax=far, pcmin=far, smooth_cap=False)
    x0, lo, hi = _FIT_SPECS[kind](pmax, 0.1 * pmax)
    sw_hi = float(sw_avg.min()) - 1e-3
    x0, lo, hi = x0 + [max(sw_hi - 0.03, 0.0)], lo + [0.0], hi + [max(sw_hi, 1e-3)]
    n = len(x0) - 1

    def resid(th):
        q = _pc_from_vector(kind, th[:n], base)
        return centrifuge_avg_sw(rpm, r2_cm, L_cm, drho_gcc, float(th[-1]), sorw, q, nr=120) - sw_avg

    r = least_squares(resid, x0, bounds=(lo, hi), max_nfev=300)
    fitted = _pc_from_vector(kind, r.x[:n], base)
    return fitted, float(np.sqrt(np.mean(r.fun ** 2))), float(r.x[-1])


# --------------------------------------------------------------------------- #
# 10. Parameter helpers: Brooks-Corey/Burdine link, Corey -> LET, low/base/high interpolation
# --------------------------------------------------------------------------- #
def brooks_corey_krw_exponent(lam: float) -> float:
    """Wetting-phase kr exponent implied by the Burdine / Brooks-Corey model: krw = Se^((2 + 3*lam)/lam).
    lam = 2 gives 4 (Corey 1954)."""
    return (2.0 + 3.0 * lam) / lam


def brooks_corey_kr(se, lam: float):
    """Brooks-Corey (1964) kr from the pore-size index lam on effective wetting saturation Se:
    krw = Se^((2+3 lam)/lam);  kr_nw = (1-Se)^2 * (1 - Se^((2+lam)/lam))."""
    se = _clip01(se)
    return se ** ((2.0 + 3.0 * lam) / lam), (1.0 - se) ** 2 * (1.0 - se ** ((2.0 + lam) / lam))


def corey_to_let(n: float, e_fixed: float | None = None) -> KrModel:
    """LET parameters whose curve is closest (least squares on 101 points) to a Corey curve with exponent n.
    Offered as a starting point when you want LET flexibility around a Corey-based estimate."""
    from scipy.optimize import least_squares
    x = np.linspace(0.0, 1.0, 101)
    target = x ** n

    def curve(th):
        L, E, T = (th[0], e_fixed if e_fixed is not None else th[1], th[-1])
        return KrModel("let", L=L, E=E, T=T).f(x)

    x0, lo, hi = ([n, 1.0, 1.0], [0.3, 0.05, 0.3], [10.0, 30.0, 10.0])
    r = least_squares(lambda th: curve(th) - target, x0, bounds=(lo, hi))
    L, E, T = r.x[0], (e_fixed if e_fixed is not None else r.x[1]), r.x[2]
    return KrModel("let", L=float(L), E=float(E), T=float(T))


def interpolate_triple(low: WOParams, base: WOParams, high: WOParams, t: float) -> WOParams:
    """Piecewise-linear interpolation through three curve sets: t = 0 -> low, 0.5 -> base, 1 -> high."""
    t = float(np.clip(t, 0.0, 1.0))
    return interpolate_wo(low, base, 2.0 * t) if t <= 0.5 else interpolate_wo(base, high, 2.0 * t - 1.0)


def interpolate_scalar_triple(lo: float, mid: float, hi: float, t: float, log: bool = False) -> float:
    """Same piecewise rule for a scalar (e.g. a Pc multiplier); log=True interpolates ln(value)."""
    f = np.log if log else (lambda v: v)
    g = np.exp if log else (lambda v: v)
    a, b, c = f(lo), f(mid), f(hi)
    t = float(np.clip(t, 0.0, 1.0))
    return float(g(a + (b - a) * 2 * t) if t <= 0.5 else g(b + (c - b) * (2 * t - 1)))
