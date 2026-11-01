# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
"""pyscal_bridge.py - optional cross-check of the native engine against the real pyscal package."""
from __future__ import annotations

import numpy as np

import core as C

try:  # pyscal is optional
    import pyscal  # type: ignore

    PYSCAL_AVAILABLE = True
except Exception:  # noqa: BLE001
    pyscal = None
    PYSCAL_AVAILABLE = False


def _cols(df):
    out = df.copy()
    out.columns = [str(c).upper() for c in out.columns]
    return out


def pyscal_wateroil(p: C.WOParams, h: float = 0.01):
    """Water-oil kr table from the real pyscal package (columns SW, KRW, KROW). Signatures follow the pyscal documentation:
    WaterOil(swirr, swl, swcr, sorw, socr, h, ...); add_corey_water(nw, krwend, krwmax); add_LET_water(l, e, t, krwend, krwmax)."""
    if not PYSCAL_AVAILABLE:
        raise RuntimeError("pyscal is not installed (pip install pyscal)")
    wo = pyscal.WaterOil(swl=p.swl, swcr=p.swcr, sorw=p.sorw, h=h)
    if p.krw.kind == "corey":
        wo.add_corey_water(nw=p.krw.n, krwend=p.krwend, krwmax=p.krwend)
        wo.add_corey_oil(now=p.kro.n, kroend=p.kroend)
    else:
        wo.add_LET_water(l=p.krw.L, e=p.krw.E, t=p.krw.T, krwend=p.krwend, krwmax=p.krwend)
        wo.add_LET_oil(l=p.kro.L, e=p.kro.E, t=p.kro.T, kroend=p.kroend)
    return _cols(wo.table)[["SW", "KRW", "KROW"]].reset_index(drop=True)


def pyscal_gasoil(p: C.GOParams, h: float = 0.01):
    """Gas-oil kr table (SG, KRG, KROG): GasOil(swirr, sgcr, h, swl, sorg, ...); add_corey_gas(ng, krgend, krgmax); add_corey_oil(nog, kroend)."""
    if not PYSCAL_AVAILABLE:
        raise RuntimeError("pyscal is not installed (pip install pyscal)")
    go = pyscal.GasOil(sgcr=p.sgcr, swl=p.swl, sorg=p.sorg, h=h)
    if p.krg.kind == "corey":
        go.add_corey_gas(ng=p.krg.n, krgend=p.krgend, krgmax=p.krgend)
        go.add_corey_oil(nog=p.kro.n, kroend=p.kroend)
    else:
        go.add_LET_gas(l=p.krg.L, e=p.krg.E, t=p.krg.T, krgend=p.krgend, krgmax=p.krgend)
        go.add_LET_oil(l=p.kro.L, e=p.kro.E, t=p.kro.T, kroend=p.kroend)
    t = _cols(go.table)
    t = t.rename(columns={c: c for c in t.columns})
    return t[["SG", "KRG", "KROG"]].reset_index(drop=True)


def crosscheck_wo(p: C.WOParams, h: float = 0.01):
    """Return (max abs difference or None, message) between native and pyscal water-oil tables,
    compared on Sw in [Swl, 1 - Sorw] (above 1 - Sorw pyscal applies its own krwmax ramp)."""
    if not PYSCAL_AVAILABLE:
        return None, "pyscal not installed."
    try:
        t = pyscal_wateroil(p, h)
        nat = C.wateroil(p, h)
        m = (t["SW"] >= p.swl - 1e-9) & (t["SW"] <= 1.0 - p.sorw + 1e-9)
        d = max(float(np.max(np.abs(np.interp(t["SW"][m], nat.SW, nat.KRW) - t["KRW"][m]))),
                float(np.max(np.abs(np.interp(t["SW"][m], nat.SW, nat.KROW) - t["KROW"][m]))))
        return d, f"pyscal cross-check: max |native - pyscal| in krw/krow = {d:.2e}"
    except Exception as e:  # noqa: BLE001
        return None, f"pyscal cross-check could not run ({e})."
