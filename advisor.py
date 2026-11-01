# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
"""advisor.py - benchmark ranges and parameter recommendations for relative permeability and Pc.

IMPORTANT: the ranges are generic, textbook-style starting points. Each parameter carries a `basis` tag:
  "Craig 1971"   a published rule of thumb (wettability classification, see documentation)
  "published"    a published range quoted in the reference list (Brooks and Corey 1964, Lomeland et al. 2005)
  "judgement"    engineering judgement consolidated from commonly quoted ranges; NOT a citation
None of them replaces site-specific SCAL data. Low / base / high refer to DISPLACEMENT EFFICIENCY:
the low case is the pessimistic set (higher Sorw, higher water mobility), the high case the optimistic one.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

import core as C

PARAMS = ["swl", "sorw", "krwend", "kroend", "nw", "no"]
LABELS = dict(swl="Swl (connate water)", sorw="Sorw (residual oil)", krwend="krw at Sorw", kroend="kro at Swl",
              nw="Corey nw", no="Corey no")


@dataclass
class Profile:
    key: str
    title: str
    description: str
    ranges: dict          # param -> (low, base, high)
    basis: dict           # param -> tag
    lam: tuple            # Brooks-Corey pore-size index (min, base, max) for Pc
    pc_mult: tuple        # Pc multiplier (low, base, high) used when interpolating
    notes: list = field(default_factory=list)


def _b(**kw):
    return kw


PROFILES = {
    "water_wet_sandstone": Profile(
        "water_wet_sandstone", "Water-wet sandstone",
        "Strongly water-wet clastic reservoir: high bound water, low krw end-point, crossover above Sw = 0.5.",
        dict(swl=(0.28, 0.24, 0.20), sorw=(0.35, 0.28, 0.20), krwend=(0.28, 0.20, 0.12), kroend=(0.80, 0.90, 0.95),
             nw=(3.0, 4.0, 5.0), no=(3.0, 2.5, 2.0)),
        _b(swl="Craig 1971", sorw="judgement", krwend="Craig 1971", kroend="judgement", nw="judgement", no="judgement"),
        (1.8, 2.5, 4.2), (1.5, 1.0, 0.7),
        ["Brooks-Corey lambda: 1.8-3.7 reported for sandpacks, 4.17 for a homogeneous consolidated sandstone (Brooks and Corey 1964, as quoted by PetroWiki).",
         "Craig rules: Swc above 20-25 %, crossover above 50 %, krw at the maximum water saturation below 30 %."]),
    "mixed_wet_sandstone": Profile(
        "mixed_wet_sandstone", "Mixed-wet sandstone",
        "Intermediate wettability: lower bound water and lower Sorw after extended injection, krw end-point in the middle range.",
        dict(swl=(0.18, 0.14, 0.10), sorw=(0.25, 0.18, 0.10), krwend=(0.60, 0.45, 0.30), kroend=(0.80, 0.90, 0.95),
             nw=(2.0, 2.5, 3.5), no=(3.5, 3.0, 2.0)),
        _b(swl="judgement", sorw="judgement", krwend="judgement", kroend="judgement", nw="judgement", no="judgement"),
        (1.5, 2.2, 3.5), (1.4, 1.0, 0.7),
        ["Craig's rules only separate water-wet from oil-wet; intermediate systems give mixed signals by definition.",
         "Sorw depends strongly on throughput in mixed-wet rock: use the throughput the field will actually see."]),
    "oil_wet_sandstone": Profile(
        "oil_wet_sandstone", "Oil-wet sandstone",
        "Oil-wet clastic: low bound water, high krw end-point, crossover below Sw = 0.5, long tail of oil production.",
        dict(swl=(0.12, 0.08, 0.05), sorw=(0.30, 0.22, 0.15), krwend=(0.90, 0.70, 0.55), kroend=(0.80, 0.90, 0.95),
             nw=(1.5, 2.0, 3.0), no=(4.5, 3.5, 2.5)),
        _b(swl="Craig 1971", sorw="judgement", krwend="Craig 1971", kroend="judgement", nw="judgement", no="judgement"),
        (1.2, 2.0, 3.5), (1.5, 1.0, 0.7),
        ["Craig rules: Swc below 10-15 %, crossover below 50 %, krw at the maximum water saturation above 50 % (up to 100 %)."]),
    "carbonate_mixed_oil_wet": Profile(
        "carbonate_mixed_oil_wet", "Carbonate (mixed to oil-wet)",
        "Heterogeneous carbonate: wide Swl and Sorw range, strongly oil-wet tendency, often dual-pore behaviour.",
        dict(swl=(0.15, 0.10, 0.06), sorw=(0.35, 0.25, 0.15), krwend=(0.85, 0.65, 0.45), kroend=(0.80, 0.90, 0.95),
             nw=(1.5, 2.0, 3.0), no=(4.0, 3.0, 2.0)),
        _b(swl="Craig 1971", sorw="judgement", krwend="Craig 1971", kroend="judgement", nw="judgement", no="judgement"),
        (0.8, 1.6, 3.0), (1.8, 1.0, 0.6),
        ["Carbonates are strongly rock-type dependent: use the rock-type family (End-point scaling tab) instead of a single set.",
         "Vuggy or fractured intervals invalidate single-curve SCAL; treat the matrix curves as one input to dual-porosity models."]),
    "tight_sandstone": Profile(
        "tight_sandstone", "Tight sandstone (low permeability)",
        "High bound water, very low water end-point and strongly convex curves; Pc entry pressure is high.",
        dict(swl=(0.45, 0.35, 0.28), sorw=(0.35, 0.28, 0.22), krwend=(0.25, 0.15, 0.08), kroend=(0.60, 0.70, 0.80),
             nw=(4.0, 5.0, 6.0), no=(4.0, 3.0, 2.5)),
        _b(swl="judgement", sorw="judgement", krwend="judgement", kroend="judgement", nw="published", no="judgement"),
        (0.8, 1.3, 2.0), (2.0, 1.0, 0.6),
        ["Wide pore-size distributions give lambda below 2 (Brooks and Corey 1964). Corey exponents of 1 to 6 are the quoted range (PetroWiki).",
         "Capillary end effects and slow imbibition make lab data in tight rock difficult: prefer steady-state or centrifuge methods."]),
}

# gas-oil and gas-water are kept deliberately simple and are marked as judgement only
GO_PROFILE = dict(sgcr=(0.10, 0.05, 0.02), sorg=(0.25, 0.15, 0.08), krgend=(0.60, 0.85, 1.00), kroend=(0.80, 0.90, 0.95),
                  ng=(3.0, 2.0, 1.5), nog=(3.5, 2.5, 2.0))
GW_PROFILE = dict(swl=(0.30, 0.20, 0.12), sgcr=(0.10, 0.05, 0.02), sgr=(0.20, 0.30, 0.40), krwend=(0.30, 0.50, 0.70),
                  krgend=(0.35, 0.55, 0.80), nw=(4.0, 3.0, 2.5), ng=(3.0, 2.0, 1.5))


# ------------------------------------------------------------------ construction
def make_params(profile: Profile, case: str, kind: str = "corey") -> C.WOParams:
    """WOParams for case 'low' | 'base' | 'high' of a profile (Corey, or LET equivalents of the Corey exponents)."""
    i = dict(low=0, base=1, high=2)[case]
    r = {k: v[i] for k, v in profile.ranges.items()}
    if kind == "let":
        krw, kro = C.corey_to_let(r["nw"]), C.corey_to_let(r["no"])
    else:
        krw, kro = C.KrModel("corey", n=r["nw"]), C.KrModel("corey", n=r["no"])
    return C.WOParams(swl=r["swl"], swcr=r["swl"], sorw=r["sorw"], krwend=r["krwend"], kroend=r["kroend"], krw=krw, kro=kro)


def pc_start(profile: Profile, k_md: float, phi: float, j_entry: float = 0.3) -> dict:
    """Starting Pc parameters: Brooks-Corey lambda range and entry pressure from Leverett J (order of magnitude)."""
    sysw = C.FLUID_SYSTEMS["Oil-water (reservoir)"]
    pd_ = float(C.pc_from_j(j_entry, sysw, k_md, phi))
    return dict(lam_min=profile.lam[0], lam=profile.lam[1], lam_max=profile.lam[2], pd_bar=pd_, j_entry=j_entry,
                nw_from_lam=C.brooks_corey_krw_exponent(profile.lam[1]),
                nw_range_from_lam=(C.brooks_corey_krw_exponent(profile.lam[2]), C.brooks_corey_krw_exponent(profile.lam[0])))


def benchmark_table() -> pd.DataFrame:
    rows = []
    for p in PROFILES.values():
        for k in PARAMS:
            lo, ba, hi = p.ranges[k]
            rows.append(dict(profile=p.title, parameter=LABELS[k], low_case=lo, base=ba, high_case=hi, basis=p.basis[k]))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ assessment of a user's case
def craig_votes(swl: float, crossover: float, krw_end: float) -> dict:
    """Craig (1971) rules of thumb. Returns the criteria that point to water-wet / oil-wet and a label."""
    ww, ow = [], []
    if swl > 0.20:
        ww.append("Swc above 20 %")
    if swl < 0.15:
        ow.append("Swc below 15 %")
    if crossover == crossover:  # not NaN
        (ww if crossover > 0.5 else ow).append("crossover above 50 %" if crossover > 0.5 else "crossover below 50 %")
    if krw_end < 0.30:
        ww.append("krw end-point below 30 %")
    if krw_end > 0.50:
        ow.append("krw end-point above 50 %")
    label = "water-wet" if len(ww) >= 2 and len(ow) == 0 else ("oil-wet" if len(ow) >= 2 and len(ww) == 0 else "intermediate / mixed")
    return dict(water_wet=ww, oil_wet=ow, label=label)


def assess(p: C.WOParams, mu_w: float, mu_o: float, profile: Profile | None = None, lam: float | None = None) -> list:
    """Findings for a case: Craig class, mobility ratio, position within a benchmark profile, Pc/kr consistency.
    Returns a list of (severity, message) with severity 'ok' | 'info' | 'warn'."""
    out = []
    t = C.wateroil(p, 0.01)
    w = C.welge(t, mu_w, mu_o, swi=p.swl)
    cv = craig_votes(p.swl, w["sw_cross"], p.krwend)
    out.append(("info", f"Craig rules read this case as {cv['label']} "
                        f"(water-wet signs: {', '.join(cv['water_wet']) or 'none'}; oil-wet signs: {', '.join(cv['oil_wet']) or 'none'})."))
    m = w["mobility_ratio"]
    out.append(("warn" if m > 3 else "info", f"End-point mobility ratio {m:.1f}: " + (
        "strongly unfavourable, expect early breakthrough and a long tail; consider mobility control." if m > 3 else
        "unfavourable but moderate." if m > 1 else "favourable, near piston-like displacement.")))
    if profile is not None:
        cur = dict(swl=p.swl, sorw=p.sorw, krwend=p.krwend, kroend=p.kroend, nw=p.krw.n if p.krw.kind == "corey" else np.nan,
                   no=p.kro.n if p.kro.kind == "corey" else np.nan)
        for k in PARAMS:
            v = cur[k]
            if v != v:
                continue
            lo, hi = min(profile.ranges[k]), max(profile.ranges[k])
            if v < lo - 1e-9 or v > hi + 1e-9:
                out.append(("warn", f"{LABELS[k]} = {v:.2f} lies outside the {profile.title.lower()} range {lo:.2f} to {hi:.2f}: justify it with data."))
        cls = {"water_wet_sandstone": "water-wet", "tight_sandstone": "water-wet", "oil_wet_sandstone": "oil-wet",
               "carbonate_mixed_oil_wet": "oil-wet"}.get(profile.key)
        if cls and cv["label"] != cls and cv["label"] != "intermediate / mixed":
            out.append(("warn", f"Craig rules read the case as {cv['label']} but the benchmark profile is {cls}."))
    if lam is not None and p.krw.kind == "corey":
        n_bc = C.brooks_corey_krw_exponent(lam)
        if abs(p.krw.n - n_bc) > 1.0:
            out.append(("info", f"Pc lambda = {lam:.2f} implies a Brooks-Corey/Burdine krw exponent of about {n_bc:.1f}; "
                                f"the case uses nw = {p.krw.n:.1f}. A mismatch is normal for rock with hysteresis or mixed wettability, but check the pair."))
    if p.krw.kind == "let":
        for name, m_ in (("krw", p.krw), ("kro", p.kro)):
            if m_.L < 1 or m_.E <= 0 or m_.T < 0.5:
                out.append(("warn", f"{name} LET parameters L={m_.L:.2f}, E={m_.E:.2f}, T={m_.T:.2f} violate the guidance L >= 1, E > 0, T >= 0.5 (Lomeland et al. 2005)."))
    if not any(s == "warn" for s, _ in out):
        out.append(("ok", "No parameter is outside the benchmark guidance."))
    return out
