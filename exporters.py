# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
"""
exporters.py - write SCAL tables for reservoir simulators.

Supported targets
  eclipse   Eclipse 100/300 keywords. Also read by OPM Flow, tNavigator, INTERSECT (Eclipse-style
            decks) and Petrel (as an include file).
  cmg       CMG IMEX / GEM table keywords (*SWT, *SLT, *KROIL).
  mrst      MATLAB script for MRST (interp1-based fluid functions + raw tables).
  csv / xlsx / json   generic tables.

Each "region" is a dict with optional keys:
  wo  DataFrame [SW, KRW, KROW, PC]
  go  DataFrame [SG, KRG, KROG, PC]
  gw  DataFrame [SW, KRW, KRG, PC]
  name  str
  provenance  optional per-system model parameters and source details
  sampling  optional export sampling details (set by prepare_regions)
  imb_go  optional imbibition gas table [SG, KRG, KROG, PC] (hysteresis)
"""
from __future__ import annotations

import io
import json
from dataclasses import asdict
from datetime import date

import numpy as np
import pandas as pd

from core import BAR_TO_KPA, BAR_TO_PSI, GOParams, GWParams, PcParams, WOParams, sof3_table

PC_FACTOR = {"METRIC": 1.0, "FIELD": BAR_TO_PSI, "SI-kPa": BAR_TO_KPA}


def model_provenance(params: WOParams | GOParams | GWParams, pc: PcParams | None = None, *,
                     resolution: float = 0.01, pc_scale: float = 1.0,
                     source: str = "Explicit parameters", engine: str = "native") -> dict:
    """Capture effective parameters, rather than potentially stale sidebar settings."""
    return dict(engine=engine, pc_engine="native", source=source, resolution_dS=resolution, parameters=asdict(params),
                pc=asdict(pc) if pc is not None else {"kind": "disabled", "value_bar": 0.0},
                pc_scale=pc_scale, parameter_pressure_unit="bar")


def resample_table(df: pd.DataFrame, rows: int, anchors: tuple[float, ...] = ()) -> pd.DataFrame:
    """Linear table interpolation with exact row count and kr plateau boundaries retained."""
    if isinstance(rows, bool) or not isinstance(rows, int) or rows < 2:
        raise ValueError("Export rows must be an integer >= 2")
    sat = df.columns[0]
    d = df.sort_values(sat)
    x = d[sat].to_numpy(float)
    required = {0, len(x) - 1}
    for col in d.columns:
        if col.startswith("KR"):
            changes = np.flatnonzero(np.diff(d[col].to_numpy(float)) != 0)
            if changes.size:
                required.update((int(changes[0]), int(changes[-1]) + 1))
    # Model anchors avoid mistaking roundoff-sized kr tails for shifted residual endpoints.
    knots = np.unique(np.concatenate(([x[0], x[-1]], np.round(anchors, 10)))) if anchors else x[sorted(required)]
    if rows < len(knots):
        raise ValueError(f"Export rows must be >= {len(knots)} to retain all saturation endpoints")
    grid = np.linspace(x[0], x[-1], rows)
    available = list(range(1, rows - 1))
    for knot in knots[1:-1]:
        idx = min(available, key=lambda i: abs(grid[i] - knot))
        grid[idx] = knot
        available.remove(idx)
    grid.sort()
    return pd.DataFrame({col: grid if col == sat else np.interp(grid, x, d[col]) for col in d.columns})


def prepare_regions(regions, rows: int) -> list[dict]:
    """Prepare a separate export copy; plotted/analysed tables remain unchanged."""
    result = []
    for r in regions:
        d = dict(r)
        for key in ("wo", "go", "gw", "imb_go"):
            if key in r:
                p = r.get("provenance", {}).get(key, {}).get("parameters")
                anchors = ()
                if p:
                    if key == "wo":
                        anchors = (p["swl"], p["swcr"], 1.0 - p["sorw"], 1.0)
                    elif key == "go":
                        anchors = (0.0, p["sgcr"], 1.0 - p["swl"] - p["sorg"], 1.0 - p["swl"])
                    elif key == "gw":
                        anchors = (p["swl"], 1.0 - p["sgr"], 1.0 - p["sgcr"], 1.0)
                d[key] = resample_table(r[key], rows, anchors)
        d["sampling"] = dict(rows_per_table=rows, method="linear interpolation of source tables",
                             grid="linspace with nearest interior slots replaced by model endpoints or kr plateau boundaries")
        result.append(d)
    return result


def region_sof3(r) -> pd.DataFrame:
    d = sof3_table(r["wo"], r["go"], float(r["wo"].SW.min()))
    rows = r.get("sampling", {}).get("rows_per_table")
    return resample_table(d, rows) if rows is not None else d


def _comment_lines(values: dict, prefix: str = "") -> list[str]:
    lines = []
    for key, value in values.items():
        name = f"{prefix}.{key}" if prefix else key
        if isinstance(value, dict):
            lines.extend(_comment_lines(value, name))
        else:
            lines.append(f"{name}={json.dumps(value, ensure_ascii=True)}")
    return lines


def _provenance_comments(regions, comment="--") -> str:
    lines = [f"{comment} Reproducibility: kr uses Corey (x^n) or LET (x^L / (x^L + E*(1-x)^T)).",
             f"{comment} Saturation normalization: SCAL Studio core.wateroil/gasoil/gaswater (pyscal conventions).",
             f"{comment} Full parameter configurations follow; only parameters of the selected model are active.",
             f"{comment} Saturations and kr are dimensionless; model Pc parameters below are in bar.",
             f"{comment} Pc values in numeric tables are converted to the declared output units.",
             f"{comment} Numeric table values are rounded to 6 decimal places."]
    for i, r in enumerate(regions, 1):
        name = json.dumps(r.get("name", ""), ensure_ascii=True)
        lines.append(f"{comment} Region {i}: {name}")
        provenance = r.get("provenance")
        if provenance:
            lines.extend(f"{comment} R{i}.{line}" for line in _comment_lines(provenance))
        else:
            lines.append(f"{comment} R{i}: source tables supplied without model parameters; no model inferred.")
        if "sampling" in r:
            lines.extend(f"{comment} R{i}.sampling.{line}" for line in _comment_lines(r["sampling"]))
        for key in ("wo", "go", "gw", "imb_go"):
            if key in r:
                lines.append(f"{comment} R{i}.{key}.rows={len(r[key])}")
    return "\n".join(lines) + "\n"


def _rows(df: pd.DataFrame, cols, pcf: float = 1.0, prec: int = 6) -> str:
    lines = []
    for r in df[cols].itertuples(index=False):
        vals = []
        for c, v in zip(cols, r):
            v = float(v) * (pcf if c == "PC" else 1.0)
            vals.append(f"{v:.{prec}f}")
        lines.append("  " + "  ".join(f"{x:>10}" for x in vals))
    return "\n".join(lines)


def _table(keyword: str, header: str, df, cols, pcf=1.0, trailing_slash=True) -> str:
    return f"{keyword}\n-- {header}\n{_rows(df, cols, pcf)}\n/\n" if trailing_slash else ""


def _banner(title: str, comment: str = "--") -> str:
    return (f"{comment} {'=' * 70}\n{comment} {title}\n{comment} Generated by SCAL Studio on {date.today().isoformat()}\n"
            f"{comment} Made by Merouane Hamdani | Educational use only | MIT License\n"
            f"{comment} Not validated for field development, reserves or investment decisions\n"
            f"{comment} {'=' * 70}\n")


# ------------------------------------------------------------------ Eclipse family
def export_eclipse(regions, units="METRIC", family="SWOF/SGOF", stone="STONE2", hysteresis=None,
                   phases=("oil", "water", "gas")) -> str:
    """family: 'SWOF/SGOF' | 'SWFN/SGFN/SOF3' | 'SWOF/SLGOF'."""
    pcf = PC_FACTOR.get(units, 1.0)
    pcu = {"METRIC": "barsa", "FIELD": "psia", "SI-kPa": "kPa"}.get(units, "barsa")
    out = [_banner(f"Saturation functions - {len(regions)} region(s), Pc in {pcu}")]
    out.append("-- Load in the PROPS section. Table N corresponds to SATNUM = N.\n")
    out.append(f"-- Export settings: family={family}; units={units}; three_phase={stone or 'none'}\n")
    out.append(_provenance_comments(regions))
    if hysteresis:
        out.append(f"-- Hysteresis: {hysteresis['scheme']} scanning curves, Land trapping "
                   f"(max trapped Sg = {hysteresis['s_trap_max']:.3f}).\n"
                   "-- Imbibition tables follow the drainage tables (IMBNUM region N+nreg).\n"
                   "-- Verify EHYSTR items against your simulator manual before running.\n"
                   "EHYSTR\n  0.1 /\n")
    if any("go" in r and "wo" in r for r in regions) and stone:
        out.append(f"{stone}\n")
    has_wo = any("wo" in r for r in regions)
    has_go = any("go" in r for r in regions)
    has_gw = any("gw" in r for r in regions)

    def block(getter, keyword, header, cols, label):
        s = [f"{keyword}\n-- {header}"]
        for i, r in enumerate(regions, start=1):
            df = getter(r)
            if df is None:
                continue
            s.append(f"-- region {i}: {r.get('name', '')}".rstrip())
            s.append(_rows(df, cols, pcf) + "\n/")
        s.append("")
        return "\n".join(s)

    if has_wo and family in ("SWOF/SGOF", "SWOF/SLGOF"):
        out.append(block(lambda r: r.get("wo"), "SWOF", f"Sw  Krw  Krow  Pcow[{pcu}]", ["SW", "KRW", "KROW", "PC"], "wo"))
    if has_go and family == "SWOF/SGOF":
        out.append(block(lambda r: r.get("go"), "SGOF", f"Sg  Krg  Krog  Pcgo[{pcu}]", ["SG", "KRG", "KROG", "PC"], "go"))
    if has_go and family == "SWOF/SLGOF":
        def slgof(r):
            g = r.get("go")
            if g is None:
                return None
            d = g.copy()
            d["SL"] = 1.0 - d.SG
            return d.sort_values("SL")
        out.append(block(slgof, "SLGOF", f"Sl  Krg  Krog  Pcgo[{pcu}]", ["SL", "KRG", "KROG", "PC"], "go"))
    if family == "SWFN/SGFN/SOF3":
        def swfn(r):
            if "wo" in r:
                return r["wo"][["SW", "KRW", "PC"]]
            if "gw" in r:
                return r["gw"][["SW", "KRW", "PC"]]
        out.append(block(swfn, "SWFN", f"Sw  Krw  Pc[{pcu}]", ["SW", "KRW", "PC"], "w"))
        def sgfn(r):
            if "go" in r:
                return r["go"][["SG", "KRG", "PC"]]
            if "gw" in r:
                d = r["gw"].copy()
                d["SG"] = 1.0 - d.SW
                d["PC"] = 0.0
                return d.sort_values("SG")[["SG", "KRG", "PC"]]
        if has_go or has_gw:
            out.append(block(sgfn, "SGFN", f"Sg  Krg  Pcgo[{pcu}]", ["SG", "KRG", "PC"], "g"))
        if has_go and has_wo:
            out.append(block(lambda r: region_sof3(r) if "wo" in r and "go" in r else None,
                             "SOF3", "So  Krow  Krog", ["SO", "KROW", "KROG"], "o"))
    if has_gw and family != "SWFN/SGFN/SOF3":
        out.append("-- Gas-water system: use SWFN/SGFN (Pcgw carried in SWFN)\n")
        out.append(block(lambda r: r["gw"][["SW", "KRW", "PC"]] if "gw" in r else None, "SWFN", f"Sw  Krw  Pcgw[{pcu}]", ["SW", "KRW", "PC"], "w"))
        def sgfn2(r):
            if "gw" not in r:
                return None
            d = r["gw"].copy()
            d["SG"] = 1.0 - d.SW
            d["PC"] = 0.0
            return d.sort_values("SG")[["SG", "KRG", "PC"]]
        out.append(block(sgfn2, "SGFN", "Sg  Krg  Pc", ["SG", "KRG", "PC"], "g"))
    if hysteresis and any("imb_go" in r for r in regions):
        out.append("-- ---- Imbibition (hysteresis) tables: region numbers continue after drainage tables ----\n")
        out.append(block(lambda r: r.get("imb_go"), "SGOF", f"Sg  Krg_imb  Krog  Pcgo[{pcu}]", ["SG", "KRG", "KROG", "PC"], "imb"))
    return "\n".join(out)


# ------------------------------------------------------------------ CMG
def export_cmg(regions, units="SI-kPa", stone="STONE2") -> str:
    pcf = PC_FACTOR.get(units, BAR_TO_KPA)
    pcu = {"METRIC": "bar", "FIELD": "psi", "SI-kPa": "kPa"}.get(units, "kPa")
    out = [_banner("CMG IMEX / GEM rock-fluid data", "**")]
    out.append("** Place inside the *ROCKFLUID section. Table sets are numbered with *RPT.\n")
    out.append(_provenance_comments(regions, "**"))
    for i, r in enumerate(regions, start=1):
        out.append(f"*RPT {i}   ** {r.get('name', '')}".rstrip())
        if "wo" in r:
            out.append(f"*SWT   ** Sw  Krw  Krow  Pcow[{pcu}]")
            out.append(_rows(r["wo"], ["SW", "KRW", "KROW", "PC"], pcf))
        if "go" in r:
            d = r["go"].copy()
            d["SL"] = 1.0 - d.SG
            d = d.sort_values("SL")
            out.append(f"*SLT   ** Sl  Krg  Krog  Pcog[{pcu}]")
            out.append(_rows(d, ["SL", "KRG", "KROG", "PC"], pcf))
        if "gw" in r:
            out.append(f"*SWT   ** Sw  Krw  Krow(=Krg proxy, gas-water)  Pcgw[{pcu}]")
            d = r["gw"].copy()
            out.append(_rows(d, ["SW", "KRW", "KRG", "PC"], pcf))
        if "wo" in r and "go" in r and stone:
            out.append(f"*KROIL *{stone}")
        out.append("")
    out.append("** Hysteresis: define imbibition curves and activate with the *HYSKRG / *HYSKRO family per the\n"
               "** IMEX/GEM manual; the imbibition boundary curve is available in the CSV export.\n")
    return "\n".join(out)


# ------------------------------------------------------------------ MRST
def _mat_vec(name: str, v) -> str:
    return f"{name} = [" + "; ".join(f"{float(x):.8g}" for x in v) + "];"


def export_mrst(regions, units="METRIC") -> str:
    out = [_banner("MRST saturation functions (SI units, Pc in Pa)", "%")]
    out.append("function fluid = scal_fluid(region)\n% Usage: fluid = scal_fluid(1);  then fluid.krW(sw), fluid.krOW(sw), ...\n"
               "switch region")
    for i, r in enumerate(regions, start=1):
        out.append(f"case {i}")
        if "wo" in r:
            w = r["wo"]
            out += [_mat_vec("  sw", w.SW), _mat_vec("  krw", w.KRW), _mat_vec("  krow", w.KROW), _mat_vec("  pcow", w.PC * 1e5)]
            out += ["  fluid.krW  = @(s) interp1(sw, krw,  s, 'linear', 'extrap');",
                    "  fluid.krOW = @(s) interp1(sw, krow, s, 'linear', 'extrap');",
                    "  fluid.pcOW = @(s) interp1(sw, pcow, s, 'linear', 'extrap');"]
        if "go" in r:
            g = r["go"]
            out += [_mat_vec("  sg", g.SG), _mat_vec("  krg", g.KRG), _mat_vec("  krog", g.KROG), _mat_vec("  pcgo", g.PC * 1e5)]
            out += ["  fluid.krG  = @(s) interp1(sg, krg,  s, 'linear', 'extrap');",
                    "  fluid.krOG = @(s) interp1(sg, krog, s, 'linear', 'extrap');",
                    "  fluid.pcGO = @(s) interp1(sg, pcgo, s, 'linear', 'extrap');"]
    out.append("end\nend\n")
    return "\n".join(out)


# ------------------------------------------------------------------ generic
def export_csv(regions) -> bytes:
    frames = []
    for i, r in enumerate(regions, start=1):
        for key in ("wo", "go", "gw", "imb_go"):
            if key in r:
                d = r[key].copy()
                d.insert(0, "SYSTEM", key)
                d.insert(0, "REGION", i)
                frames.append(d)
    return pd.concat(frames, ignore_index=True, sort=False).to_csv(index=False).encode()


def export_xlsx(regions) -> bytes:
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as xw:
        for i, r in enumerate(regions, start=1):
            for key in ("wo", "go", "gw", "imb_go"):
                if key in r:
                    r[key].to_excel(xw, sheet_name=f"R{i}_{key}"[:31], index=False)
    return buf.getvalue()


def export_json(regions, meta: dict | None = None) -> bytes:
    payload = {"meta": meta or {}, "regions": []}
    for i, r in enumerate(regions, start=1):
        item = {"region": i, "name": r.get("name", "")}
        for key in ("provenance", "sampling"):
            if key in r:
                item[key] = r[key]
        for key in ("wo", "go", "gw", "imb_go"):
            if key in r:
                item[key] = r[key].round(8).to_dict(orient="list")
        payload["regions"].append(item)
    return json.dumps(payload, indent=2, default=float).encode()


SIMULATORS = {
    "Eclipse 100 / 300 (.inc)": ("eclipse", "scal_eclipse.inc"),
    "OPM Flow (.inc)": ("eclipse", "scal_opm.inc"),
    "tNavigator (.inc)": ("eclipse", "scal_tnav.inc"),
    "INTERSECT (Eclipse-style .inc)": ("eclipse", "scal_ix.inc"),
    "Petrel include (.inc)": ("eclipse", "scal_petrel.inc"),
    "CMG IMEX / GEM (.dat)": ("cmg", "scal_cmg.inc"),
    "MRST (MATLAB .m)": ("mrst", "scal_fluid.m"),
}
