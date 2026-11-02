# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
"""deckcheck.py - independent round-trip verification of exported decks.

An Eclipse-style keyword reader (and small CMG / MRST readers) parses the exported text back into numbers and
compares it with the tables held in the app, so the export is checked by a different code path than the writer.
Each check returns a Check(name, ok, detail). Pure Python, no Streamlit.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np

from core import sof3_table
from exporters import PC_FACTOR, resample_table

_KW = re.compile(r"^[A-Z][A-Z0-9_]{1,7}$")
ATOL = 6e-6      # exporters write 6 decimals


@dataclass
class Check:
    name: str
    ok: bool
    detail: str = ""


# ------------------------------------------------------------------ Eclipse reader
def parse_eclipse(text: str) -> dict:
    """Return {KEYWORD: [array per slash-terminated table, ...]}; comments (--) are removed, tables may repeat."""
    out: dict = {}
    kw, rows = None, []
    for raw in text.splitlines():
        line = raw.split("--", 1)[0].strip()
        if not line:
            continue
        if _KW.match(line):
            kw, rows = line, []
            out.setdefault(kw, [])
            continue
        if kw is None:
            continue
        closed = "/" in line
        body = line.replace("/", " ").split()
        if body:
            try:
                rows.append([float(x) for x in body])
            except ValueError:
                raise ValueError(f"Cannot parse line under {kw}: '{raw.strip()}'")
        if closed:
            width = max((len(r) for r in rows), default=0)
            out[kw].append(np.array([r + [np.nan] * (width - len(r)) for r in rows]) if rows else np.empty((0, 0)))
            rows = []
    return out


_DECK_COLS = {"SWOF": ["SW", "KRW", "KROW", "PC"], "SGOF": ["SG", "KRG", "KROG", "PC"], "SLGOF": ["SL", "KRG", "KROG", "PC"],
              "SWFN": ["SW", "KRW", "PC"], "SGFN": ["SG", "KRG", "PC"], "SOF3": ["SO", "KROW", "KROG"]}


def tables_from_eclipse(text: str) -> dict:
    """Read saturation-function tables from an Eclipse-style deck: {keyword: [DataFrame per region]} with named columns.
    Keywords: SWOF SGOF SLGOF SWFN SGFN SOF3. Tables with a different number of columns than the keyword defines are skipped."""
    import pandas as pd
    out = {}
    for kw, tabs in parse_eclipse(text).items():
        cols = _DECK_COLS.get(kw)
        if not cols:
            continue
        dfs = [pd.DataFrame(t[:, :len(cols)], columns=cols) for t in tabs if t.ndim == 2 and t.shape[1] >= len(cols) and t.shape[0] > 0]
        if dfs:
            out[kw] = dfs
    return out


def _cmp(name, got, exp, atol=ATOL):
    exp = np.asarray(exp, float)
    if got is None or got.shape != exp.shape:
        return Check(name, False, f"shape {None if got is None else got.shape} vs expected {exp.shape}")
    err = float(np.nanmax(np.abs(got - exp))) if exp.size else 0.0
    tol = atol * max(1.0, float(np.nanmax(np.abs(exp))) if exp.size else 1.0)
    return Check(name, err <= tol, f"max abs diff {err:.2e}")


def _sanity(name, tab, sat_col=0, krcols=(1, 2)):
    sat = tab[:, sat_col]
    ok = bool(np.all(np.diff(sat) > 0)) and all(np.all((tab[:, c] >= -1e-9) & (tab[:, c] <= 1 + 1e-9)) for c in krcols)
    return Check(name + " sane", ok, "saturation strictly increasing, kr within [0,1]" if ok else "saturation not increasing or kr outside [0,1]")


def verify_eclipse(text, regions, units="METRIC", family="SWOF/SGOF", hysteresis=None) -> list:
    pcf = PC_FACTOR.get(units, 1.0)
    P = parse_eclipse(text)
    res = []

    def check(keyword, expected, cols=None):
        got = P.get(keyword, [])
        if not expected:
            if got:
                res.append(Check(keyword, False, f"{len(got)} unexpected table(s)"))
            return
        if len(got) != len(expected):
            res.append(Check(keyword, False, f"{len(got)} tables found, {len(expected)} expected"))
            return
        for i, (g, e) in enumerate(zip(got, expected), 1):
            res.append(_cmp(f"{keyword} table {i}", g, e))
            if g.size and g.shape[1] >= 3 and keyword not in ("SOF3",):
                res.append(_sanity(f"{keyword} table {i}", g, 0, (1, 2) if keyword in ("SWOF", "SGOF", "SLGOF") else (1,)))

    def arr(df, cols):
        a = df[cols].to_numpy(float).copy()
        if "PC" in cols:
            a[:, cols.index("PC")] *= pcf
        return a

    wo = [r for r in regions if "wo" in r]
    go = [r for r in regions if "go" in r]
    gw = [r for r in regions if "gw" in r]
    if family in ("SWOF/SGOF", "SWOF/SLGOF"):
        check("SWOF", [arr(r["wo"], ["SW", "KRW", "KROW", "PC"]) for r in wo])
    if family == "SWOF/SGOF":
        exp = [arr(r["go"], ["SG", "KRG", "KROG", "PC"]) for r in go]
        if hysteresis:
            exp += [arr(r["imb_go"], ["SG", "KRG", "KROG", "PC"]) for r in regions if "imb_go" in r]
        check("SGOF", exp)
    if family == "SWOF/SLGOF":
        exp = []
        for r in go:
            d = r["go"].copy(); d["SL"] = 1.0 - d.SG
            exp.append(arr(d.sort_values("SL"), ["SL", "KRG", "KROG", "PC"]))
        check("SLGOF", exp)
    if family == "SWFN/SGFN/SOF3":
        exp_w = [arr((r["wo"] if "wo" in r else r["gw"])[["SW", "KRW", "PC"]], ["SW", "KRW", "PC"]) for r in regions if "wo" in r or "gw" in r]
        check("SWFN", exp_w)
        exp_g = []
        for r in regions:
            if "go" in r:
                exp_g.append(arr(r["go"], ["SG", "KRG", "PC"]))
            elif "gw" in r:
                d = r["gw"].copy(); d["SG"] = 1.0 - d.SW; d["PC"] = 0.0
                exp_g.append(arr(d.sort_values("SG"), ["SG", "KRG", "PC"]))
        check("SGFN", exp_g)
        both = [r for r in regions if "wo" in r and "go" in r]
        if both:
            exp_o = []
            for r in both:
                d = sof3_table(r["wo"], r["go"], float(r["wo"].SW.min()))
                rows = r.get("sampling", {}).get("rows_per_table")
                if rows is not None:
                    d = resample_table(d, rows)
                exp_o.append(d[["SO", "KROW", "KROG"]].to_numpy(float))
            check("SOF3", exp_o)
    else:
        if gw:
            check("SWFN", [arr(r["gw"], ["SW", "KRW", "PC"]) for r in gw])
            exp = []
            for r in gw:
                d = r["gw"].copy(); d["SG"] = 1.0 - d.SW; d["PC"] = 0.0
                exp.append(arr(d.sort_values("SG"), ["SG", "KRG", "PC"]))
            check("SGFN", exp)
    if hysteresis:
        res.append(Check("EHYSTR present", "EHYSTR" in P and len(P["EHYSTR"]) == 1, "keyword found" if "EHYSTR" in P else "missing"))
    text_up = text.upper()
    res.append(Check("credit banner", "MEROUANE HAMDANI".upper() in text_up and "MIT" in text_up, "author + licence present"))
    return res


# ------------------------------------------------------------------ CMG reader
def parse_cmg(text: str) -> dict:
    """{region: {'SWT': array, 'SLT': array}} from a *RPT based file."""
    out, region, cur = {}, None, None
    for raw in text.splitlines():
        line = raw.split("**", 1)[0].strip()
        if not line:
            continue
        m = re.match(r"^\*RPT\s+(\d+)", line, re.I)
        if m:
            region = int(m.group(1)); out[region] = {}; cur = None
            continue
        m = re.match(r"^\*(SWT|SLT)\b", line, re.I)
        if m and region is not None:
            cur = m.group(1).upper(); out[region][cur] = []
            continue
        if line.startswith("*"):
            cur = None
            continue
        if cur and region is not None:
            out[region][cur].append([float(x) for x in line.split()])
    return {r: {k: np.array(v) for k, v in d.items()} for r, d in out.items()}


def verify_cmg(text, regions, units="SI-kPa") -> list:
    pcf = PC_FACTOR.get(units, 1.0)
    P = parse_cmg(text)
    res = [Check("*RPT count", len(P) == len(regions), f"{len(P)} found, {len(regions)} expected")]
    for i, r in enumerate(regions, 1):
        d = P.get(i, {})
        if "wo" in r:
            e = r["wo"][["SW", "KRW", "KROW", "PC"]].to_numpy(float).copy(); e[:, 3] *= pcf
            res.append(_cmp(f"*SWT region {i}", d.get("SWT"), e))
        if "go" in r:
            g = r["go"].copy(); g["SL"] = 1.0 - g.SG
            e = g.sort_values("SL")[["SL", "KRG", "KROG", "PC"]].to_numpy(float).copy(); e[:, 3] *= pcf
            res.append(_cmp(f"*SLT region {i}", d.get("SLT"), e))
    return res


# ------------------------------------------------------------------ MRST reader
def parse_mrst(text: str) -> dict:
    """{region: {name: array}} from the 'case N' blocks."""
    out, region = {}, None
    for raw in text.splitlines():
        line = raw.split("%", 1)[0].strip()
        m = re.match(r"^case\s+(\d+)", line)
        if m:
            region = int(m.group(1)); out[region] = {}
            continue
        m = re.match(r"^(\w+)\s*=\s*\[([^\]]*)\];", line)
        if m and region is not None:
            out[region][m.group(1)] = np.array([float(x) for x in m.group(2).split(";") if x.strip()])
    return out


def verify_mrst(text, regions) -> list:
    P = parse_mrst(text)
    res = [Check("case count", len(P) == len(regions), f"{len(P)} found, {len(regions)} expected")]
    for i, r in enumerate(regions, 1):
        d = P.get(i, {})
        if "wo" in r:
            for key, col, f in (("sw", "SW", 1), ("krw", "KRW", 1), ("krow", "KROW", 1), ("pcow", "PC", 1e5)):
                e = r["wo"][col].to_numpy(float) * f
                g = d.get(key)
                ok = g is not None and g.shape == e.shape and np.allclose(g, e, rtol=1e-6, atol=1e-9 * max(1, np.abs(e).max()))
                res.append(Check(f"{key} region {i}", bool(ok), "match" if ok else "mismatch"))
    return res


def summarize(checks) -> tuple[int, int]:
    return sum(c.ok for c in checks), len(checks)
