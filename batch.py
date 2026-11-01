# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
"""batch.py - run SCAL Studio without the UI: a YAML file defines the cases, the engine and the exports.

    python scal_batch.py case.yaml --out results/        (see scal_batch.py --help, --template)

The YAML is validated strictly (unknown keys, wrong types and out-of-range values stop the run with the path of the
problem). Exports are parsed back with deckcheck and compared with the tables before the run is reported as successful.
"""
from __future__ import annotations

import io
import json
import zipfile
from dataclasses import dataclass, field, fields
from datetime import date

import numpy as np
import pandas as pd
import yaml

import advisor as A
import core as C
import deckcheck as D
import exporters as X

ENGINES = ("native", "pyscal")
SYSTEM_KEYS = {"water-oil": ("wo",), "gas-oil": ("go",), "gas-water": ("gw",), "three-phase": ("wo", "go")}
FORMATS = ("eclipse", "cmg", "mrst", "csv", "xlsx", "json")
FAMILIES = ("SWOF/SGOF", "SWOF/SLGOF", "SWFN/SGFN/SOF3")
UNITS = ("METRIC", "FIELD", "SI-kPa")

TEMPLATE = """\
# SCAL Studio batch input  (made by Merouane Hamdani, educational use, MIT License)
project: demo_batch
engine: native              # native | pyscal  (pyscal: water-oil and gas-oil kr tables from the pyscal package; Pc from the native engine)
resolution: 0.01            # saturation step dS
units: METRIC               # Pc unit in the files: METRIC (bar) | FIELD (psi) | SI-kPa (kPa)
mu_w: 0.5                   # cP, used for the summary table only
mu_o: 3.0

export:
  formats: [eclipse, cmg, mrst, csv, xlsx, json]
  eclipse_family: SWOF/SGOF # SWOF/SGOF | SWOF/SLGOF | SWFN/SGFN/SOF3
  stone: STONE2             # STONE1 | STONE2 | none  (three-phase cases)
  combine: true             # true: one set of files, case N = SATNUM N;  false: one folder per case
  verify: true              # parse the written decks back and compare with the tables

cases:
  # 1) explicit end-points and curve shapes
  - name: RT1_good
    system: water-oil       # water-oil | gas-oil | gas-water | three-phase
    kr_model: corey         # corey | let
    water_oil:
      swl: 0.15
      swcr: 0.15
      sorw: 0.22
      krwend: 0.35
      kroend: 0.90
      krw: {n: 3.5}         # LET: {L: 2.0, E: 3.0, T: 1.5}
      kro: {n: 2.5}
    pc:
      water_oil: {kind: brooks-corey, pd: 0.25, pcmax: 6.0, lam: 2.2}

  # 2) from a benchmark profile; t = 0 low case, 0.5 base, 1 high (piecewise-linear interpolation)
  - name: RT2_benchmark_t03
    system: water-oil
    kr_model: corey
    advisor: {profile: water_wet_sandstone, t: 0.3}
    pc:
      water_oil: {kind: brooks-corey, pd: 0.35, pcmax: 6.0, lam: 2.5}

  # 3) interpolate between two cases defined above (end-points and shapes; t = 0 -> from, 1 -> to)
  - name: RT3_between
    system: water-oil
    kr_model: corey
    interpolate: {from: RT1_good, to: RT2_benchmark_t03, t: 0.5}
"""


class BatchError(ValueError):
    """Invalid batch input; the message carries the path of the offending key."""


# ------------------------------------------------------------------ validation helpers
def _dict(v, path, allowed=None, required=()):
    if not isinstance(v, dict):
        raise BatchError(f"{path}: expected a mapping")
    if allowed is not None:
        bad = sorted(set(v) - set(allowed))
        if bad:
            raise BatchError(f"{path}: unknown key(s) {bad}; allowed: {sorted(allowed)}")
    miss = [k for k in required if k not in v]
    if miss:
        raise BatchError(f"{path}: missing required key(s) {miss}")
    return v


def _num(v, path, lo=None, hi=None):
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        raise BatchError(f"{path}: expected a number, got {v!r}")
    if (lo is not None and v < lo) or (hi is not None and v > hi):
        raise BatchError(f"{path}: {v} outside [{lo}, {hi}]")
    return float(v)


def _choice(v, path, options):
    if v not in options:
        raise BatchError(f"{path}: {v!r} not one of {list(options)}")
    return v


def _kr(d, path, kind):
    d = _dict(d, path)
    if kind == "corey":
        _dict(d, path, ("n",))
        return C.KrModel("corey", n=_num(d.get("n", 2.0), f"{path}.n", 0.5, 8.0))
    _dict(d, path, ("L", "E", "T"))
    return C.KrModel("let", L=_num(d.get("L", 2.0), f"{path}.L", 0.1, 10.0), E=_num(d.get("E", 2.0), f"{path}.E", 0.05, 50.0),
                     T=_num(d.get("T", 2.0), f"{path}.T", 0.1, 10.0))


_PC_ALLOWED = ("kind", "pd", "pcmax", "pcmin", "lam", "L", "E", "T", "Ln", "En", "Tn", "vg_n", "cw", "aw", "co", "ao", "smooth_cap")


def _pc(d, path):
    d = _dict(d, path, _PC_ALLOWED)
    kw = {"kind": _choice(d.get("kind", "brooks-corey"), f"{path}.kind", C.PC_MODELS)}
    for k in _PC_ALLOWED[1:]:
        if k in d:
            kw[k] = bool(d[k]) if k == "smooth_cap" else _num(d[k], f"{path}.{k}", 0.0, 1000.0)
    p = C.PcParams(**kw)
    try:
        p.validate()
    except ValueError as exc:
        raise BatchError(f"{path}: {exc}") from exc
    return p


@dataclass
class CaseSpec:
    name: str
    system: str
    kr_model: str
    wo: C.WOParams | None = None
    go: C.GOParams | None = None
    gw: C.GWParams | None = None
    pc: dict = field(default_factory=dict)      # 'wo' | 'go' | 'gw' -> PcParams
    pc_scale: float = 1.0
    origin: str = "explicit"


def _wo(d, path, kind):
    allowed = ("swl", "swcr", "sorw", "krwend", "kroend", "krw", "kro")
    d = _dict(d, path, allowed)
    base = C.WOParams()
    swl = _num(d.get("swl", base.swl), f"{path}.swl", 0.0, 0.6)
    p = C.WOParams(swl, max(_num(d.get("swcr", swl), f"{path}.swcr", 0.0, 0.7), swl), _num(d.get("sorw", base.sorw), f"{path}.sorw", 0.0, 0.7),
                   _num(d.get("krwend", base.krwend), f"{path}.krwend", 0.001, 1.0), _num(d.get("kroend", base.kroend), f"{path}.kroend", 0.001, 1.0),
                   _kr(d.get("krw", {}), f"{path}.krw", kind), _kr(d.get("kro", {}), f"{path}.kro", kind))
    return _validated(p, path)


def _go(d, path, kind):
    d = _dict(d, path, ("swl", "sgcr", "sorg", "krgend", "kroend", "krg", "kro"))
    b = C.GOParams()
    p = C.GOParams(_num(d.get("swl", b.swl), f"{path}.swl", 0.0, 0.6), _num(d.get("sgcr", b.sgcr), f"{path}.sgcr", 0.0, 0.6),
                   _num(d.get("sorg", b.sorg), f"{path}.sorg", 0.0, 0.7), _num(d.get("krgend", b.krgend), f"{path}.krgend", 0.001, 1.0),
                   _num(d.get("kroend", b.kroend), f"{path}.kroend", 0.001, 1.0),
                   _kr(d.get("krg", {}), f"{path}.krg", kind), _kr(d.get("kro", {}), f"{path}.kro", kind))
    return _validated(p, path)


def _gw(d, path, kind):
    d = _dict(d, path, ("swl", "sgcr", "sgr", "krwend", "krgend", "krw", "krg"))
    b = C.GWParams()
    p = C.GWParams(_num(d.get("swl", b.swl), f"{path}.swl", 0.0, 0.7), _num(d.get("sgcr", b.sgcr), f"{path}.sgcr", 0.0, 0.6),
                   _num(d.get("sgr", b.sgr), f"{path}.sgr", 0.0, 0.7), _num(d.get("krwend", b.krwend), f"{path}.krwend", 0.001, 1.0),
                   _num(d.get("krgend", b.krgend), f"{path}.krgend", 0.001, 1.0),
                   _kr(d.get("krw", {}), f"{path}.krw", kind), _kr(d.get("krg", {}), f"{path}.krg", kind))
    return _validated(p, path)


def _validated(p, path):
    try:
        p.validate()
    except ValueError as exc:
        raise BatchError(f"{path}: {exc}") from exc
    return p


def load_config(text_or_dict) -> dict:
    """Parse and validate a batch configuration; returns the normalised dict with 'cases' as CaseSpec objects."""
    try:
        raw = yaml.safe_load(text_or_dict) if isinstance(text_or_dict, (str, bytes)) else text_or_dict
    except yaml.YAMLError as exc:
        raise BatchError(f"YAML syntax error: {exc}") from exc
    top = _dict(raw, "input", ("project", "engine", "resolution", "units", "mu_w", "mu_o", "export", "cases"), ("cases",))
    cfg = dict(project=str(top.get("project", "scal_batch")), engine=_choice(top.get("engine", "native"), "engine", ENGINES),
               resolution=_num(top.get("resolution", 0.01), "resolution", 0.002, 0.05), units=_choice(top.get("units", "METRIC"), "units", UNITS),
               mu_w=_num(top.get("mu_w", 0.5), "mu_w", 0.01, 100.0), mu_o=_num(top.get("mu_o", 3.0), "mu_o", 0.01, 5000.0))
    ex = _dict(top.get("export", {}), "export", ("formats", "eclipse_family", "stone", "combine", "verify"))
    fm = ex.get("formats", list(FORMATS))
    if not isinstance(fm, list) or not fm:
        raise BatchError("export.formats: expected a non-empty list")
    cfg["export"] = dict(formats=[_choice(f, "export.formats", FORMATS) for f in fm],
                         family=_choice(ex.get("eclipse_family", "SWOF/SGOF"), "export.eclipse_family", FAMILIES),
                         stone=_choice(ex.get("stone", "STONE2"), "export.stone", ("STONE1", "STONE2", "none", None)),
                         combine=bool(ex.get("combine", True)), verify=bool(ex.get("verify", True)))
    if cfg["export"]["stone"] == "none":
        cfg["export"]["stone"] = None
    cases_raw = top["cases"]
    if not isinstance(cases_raw, list) or not cases_raw:
        raise BatchError("cases: expected a non-empty list")
    done: dict[str, CaseSpec] = {}
    for i, cr in enumerate(cases_raw):
        path = f"cases[{i}]"
        cr = _dict(cr, path, ("name", "system", "kr_model", "water_oil", "gas_oil", "gas_water", "pc", "advisor", "interpolate", "pc_multiplier"), ("name",))
        name = str(cr["name"])
        if name in done:
            raise BatchError(f"{path}.name: duplicate case name '{name}'")
        system = _choice(cr.get("system", "water-oil"), f"{path}.system", SYSTEM_KEYS)
        kind = _choice(cr.get("kr_model", "corey"), f"{path}.kr_model", ("corey", "let"))
        spec = CaseSpec(name, system, kind, pc_scale=_num(cr.get("pc_multiplier", 1.0), f"{path}.pc_multiplier", 0.01, 100.0))
        pcs = _dict(cr.get("pc", {}), f"{path}.pc", ("water_oil", "gas_oil", "gas_water"))
        for key, short in (("water_oil", "wo"), ("gas_oil", "go"), ("gas_water", "gw")):
            if key in pcs:
                spec.pc[short] = _pc(pcs[key], f"{path}.pc.{key}")
        if "advisor" in cr and "interpolate" in cr:
            raise BatchError(f"{path}: use either 'advisor' or 'interpolate', not both")
        if "advisor" in cr:
            ad = _dict(cr["advisor"], f"{path}.advisor", ("profile", "t"), ("profile",))
            prof = _choice(ad["profile"], f"{path}.advisor.profile", list(A.PROFILES))
            t = _num(ad.get("t", 0.5), f"{path}.advisor.t", 0.0, 1.0)
            if "wo" not in SYSTEM_KEYS[system]:
                raise BatchError(f"{path}.advisor: profiles describe water-oil curves; system must be water-oil or three-phase")
            lo, ba, hi = (A.make_params(A.PROFILES[prof], c, kind) for c in ("low", "base", "high"))
            spec.wo = C.interpolate_triple(lo, ba, hi, t)
            spec.origin = f"advisor {prof} t={t:g}"
        elif "interpolate" in cr:
            ip = _dict(cr["interpolate"], f"{path}.interpolate", ("from", "to", "t"), ("from", "to"))
            a, b = done.get(ip["from"]), done.get(ip["to"])
            if a is None or b is None:
                raise BatchError(f"{path}.interpolate: 'from' and 'to' must name cases defined earlier in the list")
            if a.wo is None or b.wo is None:
                raise BatchError(f"{path}.interpolate: both cases need water-oil curves")
            if a.kr_model != b.kr_model:
                raise BatchError(f"{path}.interpolate: both cases must use the same kr_model")
            t = _num(ip.get("t", 0.5), f"{path}.interpolate.t", 0.0, 1.0)
            spec.wo = C.interpolate_wo(a.wo, b.wo, t)
            spec.kr_model = a.kr_model
            for k, v in {**b.pc, **a.pc}.items():
                spec.pc.setdefault(k, v)
            spec.origin = f"interpolate {ip['from']} -> {ip['to']} t={t:g}"
        elif "wo" in SYSTEM_KEYS[system]:
            spec.wo = _wo(cr.get("water_oil", {}), f"{path}.water_oil", kind)
        if "go" in SYSTEM_KEYS[system]:
            spec.go = _go(cr.get("gas_oil", {}), f"{path}.gas_oil", kind)
        if "gw" in SYSTEM_KEYS[system]:
            spec.gw = _gw(cr.get("gas_water", {}), f"{path}.gas_water", kind)
        done[name] = spec
    cfg["cases"] = list(done.values())
    if cfg["engine"] == "pyscal" and any(c.gw is not None for c in cfg["cases"]):
        raise BatchError("engine pyscal supports water-oil and gas-oil only; use engine native for gas-water cases")
    return cfg


# ------------------------------------------------------------------ building tables
def _tables(spec: CaseSpec, engine: str, h: float) -> dict:
    reg = {"name": spec.name}
    if engine == "pyscal":
        import pyscal_bridge as PB
    if spec.wo is not None:
        t = PB.pyscal_wateroil(spec.wo, h) if engine == "pyscal" else C.wateroil(spec.wo, h)
        pc = spec.pc.get("wo")
        reg["wo"] = C.add_pc_wo(t, spec.wo, pc, spec.pc_scale) if pc else t.assign(PC=0.0)
    if spec.go is not None:
        t = PB.pyscal_gasoil(spec.go, h) if engine == "pyscal" else C.gasoil(spec.go, h)
        pc = spec.pc.get("go")
        reg["go"] = C.add_pc_go(t, spec.go, pc, spec.pc_scale) if pc else t.assign(PC=0.0)
    if spec.gw is not None:
        t = C.gaswater(spec.gw, h)
        pc = spec.pc.get("gw")
        reg["gw"] = C.add_pc_gw(t, spec.gw, pc, spec.pc_scale) if pc else t.assign(PC=0.0)
    return reg


@dataclass
class BatchResult:
    ok: bool
    files: dict
    summary: pd.DataFrame
    checks: list
    warnings: list
    manifest: dict

    def zip_bytes(self) -> bytes:
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
            for name, blob in self.files.items():
                z.writestr(name, blob)
        return buf.getvalue()


def _write_set(regions, cfg, prefix) -> tuple[dict, list]:
    ex, units = cfg["export"], cfg["units"]
    files, checks = {}, []
    cmg_units = units if units != "METRIC" else "SI-kPa"
    if "eclipse" in ex["formats"]:
        txt = X.export_eclipse(regions, units, ex["family"], ex["stone"])
        files[f"{prefix}scal_eclipse.inc"] = txt.encode()
        if ex["verify"]:
            checks.append((f"{prefix}scal_eclipse.inc", D.verify_eclipse(txt, regions, units, ex["family"], None)))
    if "cmg" in ex["formats"]:
        txt = X.export_cmg(regions, cmg_units, ex["stone"])
        files[f"{prefix}scal_cmg.inc"] = txt.encode()
        if ex["verify"]:
            checks.append((f"{prefix}scal_cmg.inc", D.verify_cmg(txt, regions, cmg_units)))
    if "mrst" in ex["formats"]:
        txt = X.export_mrst(regions)
        files[f"{prefix}scal_fluid.m"] = txt.encode()
        if ex["verify"]:
            checks.append((f"{prefix}scal_fluid.m", D.verify_mrst(txt, regions)))
    if "csv" in ex["formats"]:
        files[f"{prefix}scal_tables.csv"] = X.export_csv(regions)
    if "xlsx" in ex["formats"]:
        files[f"{prefix}scal_tables.xlsx"] = X.export_xlsx(regions)
    if "json" in ex["formats"]:
        files[f"{prefix}scal_tables.json"] = X.export_json(regions, dict(pc_units=units))
    return files, checks


def run_batch(config, engine_override: str | None = None) -> BatchResult:
    """Run a validated (or raw YAML text) configuration. Raises BatchError for invalid input."""
    cfg = load_config(config) if not (isinstance(config, dict) and "export" in config and isinstance(config.get("cases", [None])[0], CaseSpec)) else config
    if engine_override:
        cfg = dict(cfg, engine=_choice(engine_override, "engine", ENGINES))
    if cfg["engine"] == "pyscal":
        import pyscal_bridge as PB
        if not PB.PYSCAL_AVAILABLE:
            raise BatchError("engine pyscal requested but the pyscal package is not installed (pip install pyscal)")
    regions, rows, warnings = [], [], []
    for spec in cfg["cases"]:
        reg = _tables(spec, cfg["engine"], cfg["resolution"])
        regions.append(reg)
        row = dict(case=spec.name, system=spec.system, kr_model=spec.kr_model, source=spec.origin, engine=cfg["engine"])
        if spec.wo is not None:
            w = C.welge(reg["wo"], cfg["mu_w"], cfg["mu_o"], swi=spec.wo.swl)
            cv = A.craig_votes(spec.wo.swl, w["sw_cross"], spec.wo.krwend)
            row.update(Swl=round(spec.wo.swl, 4), Sorw=round(spec.wo.sorw, 4), krw_end=round(spec.wo.krwend, 4), kro_end=round(spec.wo.kroend, 4),
                       RF_limit=round(w["rf_ult"], 4), breakthrough_PV=round(w["pvi_bt"], 4), mobility_ratio=round(w["mobility_ratio"], 3),
                       crossover_Sw=round(w["sw_cross"], 4), craig_reading=cv["label"])
            if spec.wo.krwend > 0.8 and cv["label"] == "water-wet":
                warnings.append(f"{spec.name}: high krw end-point with a water-wet reading; check the inputs.")
        rows.append(row)
    files, checks = {}, []
    ex = cfg["export"]
    if ex["combine"]:
        f, c = _write_set(regions, cfg, "")
        files.update(f), checks.extend(c)
    else:
        for reg in regions:
            f, c = _write_set([reg], cfg, f"cases/{reg['name']}/")
            files.update(f), checks.extend(c)
    summary = pd.DataFrame(rows)
    files["summary.csv"] = summary.to_csv(index=False).encode()
    check_rows = [dict(file=name, passed=sum(k.ok for k in cs), total=len(cs)) for name, cs in checks]
    ok = all(r["passed"] == r["total"] for r in check_rows)
    manifest = dict(app="SCAL Studio batch", author="Merouane Hamdani (educational use, MIT License)", date=date.today().isoformat(),
                    project=cfg["project"], engine=cfg["engine"], units=cfg["units"], resolution=cfg["resolution"],
                    cases=[r["case"] for r in rows], verification=check_rows, files=sorted(files), ok=ok)
    files["manifest.json"] = json.dumps(manifest, indent=2).encode()
    return BatchResult(ok, files, summary, check_rows, warnings, manifest)
