# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
"""cases.py - save / load a case as JSON with strict validation against the widget bounds.

The schema is derived from the app's default state: every numeric key gets bounds (RULES) that mirror the
widgets in app.py, so loading a file can never push a widget out of range (which would crash Streamlit).
Pure Python: no Streamlit import, fully unit tested.
"""
from __future__ import annotations

import json
import re
from datetime import datetime, timezone

APP_TAG = "SCAL Studio case"
SCHEMA_VERSION = 3

# (regex on key, lo, hi) - first match wins. Keep in sync with the widgets in app.py (a test enforces the defaults).
RULES = [
    (r"^export_rows$", 8, 1000),
    (r"^gw_swl$", 0.0, 0.6), (r"^(wo|go|sc|ip)_swl$", 0.0, 0.5),
    (r"^(wo|sc|ip)_swcr$", 0.0, 0.6), (r"^(wo|sc|ip)_sorw$", 0.0, 0.6),
    (r"^(wo|sc|ip)_krwend$", 0.02, 1.0), (r"^(wo|sc|ip)_kroend$", 0.1, 1.0),
    (r"^go_sgcr$", 0.0, 0.4), (r"^go_sorg$", 0.0, 0.5), (r"^go_krgend$", 0.05, 1.0), (r"^go_kroend$", 0.1, 1.0),
    (r"^gw_sgcr$", 0.0, 0.4), (r"^gw_sgr$", 0.0, 0.6), (r"^gw_(krwend|krgend)$", 0.05, 1.0),
    (r"^sc_krwr$", 0.02, 1.0), (r"^sc_krwmax$", 0.1, 1.0),
    (r"^(wo|go|gw)_(krw|kro|krg)_n$", 0.5, 8.0), (r"^(wo|go|gw)_(krw|kro|krg)_L$", 0.1, 10.0),
    (r"^(wo|go|gw)_(krw|kro|krg)_E$", 0.05, 50.0), (r"^(wo|go|gw)_(krw|kro|krg)_T$", 0.1, 10.0),
    (r"^mu_w$", 0.1, 10.0), (r"^mu_o$", 0.1, 500.0),
    (r"^pc(wo|go|gw)_pd$", 0.0, 50.0), (r"^pc(wo|go|gw)_pcmax$", 0.1, 200.0), (r"^pc(wo|go|gw)_pcmin$", 0.0, 50.0),
    (r"^pc(wo|go|gw)_lam$", 0.2, 10.0), (r"^pc(wo|go|gw)_L$", 0.1, 10.0), (r"^pc(wo|go|gw)_E$", 0.05, 50.0),
    (r"^pc(wo|go|gw)_T$", 0.1, 10.0), (r"^pc(wo|go|gw)_vg$", 1.2, 10.0),
    (r"^pc(wo|go|gw)_(cw|co)$", 0.0, 200.0), (r"^pc(wo|go|gw)_(aw|ao)$", 0.05, 5.0),
    (r"^ens_n$", 20, 1000), (r"^ens_seed$", 0, 99999), (r"^ens_dswl$", 0.0, 0.15), (r"^ens_dsorw$", 0.0, 0.2),
    (r"^ens_(fkrw)$", 0.0, 0.9), (r"^ens_fkro$", 0.0, 0.5), (r"^ens_(dnw|dno)$", 0.0, 2.0), (r"^ens_fpc$", 0.0, 1.0),
    (r"^ip_(dnw|dno)$", -3.0, 3.0), (r"^ip_n$", 2, 9), (r"^ip_(lo|hi|value)$", 1e-12, 1e12), (r"^ip_t$", 0.0, 1.0),
    (r"^adv_t$", 0.0, 1.0), (r"^adv_k$", 0.01, 20000.0), (r"^adv_phi$", 0.02, 0.45), (r"^adv_j$", 0.05, 2.0), (r"^ens_pctl$", 1, 99),
    (r"^ip_mu$", 0.1, 500.0), (r"^ip_v$", 1e-6, 100.0), (r"^ip_sigma$", 0.001, 80.0),
]


def bounds_for(key: str):
    for pat, lo, hi in RULES:
        if re.match(pat, key):
            return lo, hi
    return None


def schema_for(defaults: dict, choices: dict) -> dict:
    """key -> ('num', lo, hi, is_int) | ('choice', options) | ('bool',)."""
    sch = {}
    for k, v in defaults.items():
        if k in choices:
            sch[k] = ("choice", list(choices[k]))
        elif isinstance(v, bool):
            sch[k] = ("bool",)
        elif isinstance(v, (int, float)):
            b = bounds_for(k)
            if b is None:
                raise KeyError(f"no bounds rule for numeric key '{k}' - add one to cases.RULES")
            sch[k] = ("num", b[0], b[1], isinstance(v, int))
        else:
            raise KeyError(f"no choices defined for string key '{k}'")
    return sch


def sanitize(raw: dict, schema: dict):
    """Validate a raw {key: value} dict. Returns (clean, problems). Unknown or out-of-range keys are dropped and reported."""
    clean, problems = {}, []
    for k, v in raw.items():
        spec = schema.get(k)
        if spec is None:
            problems.append(f"{k}: unknown key (ignored)")
            continue
        kind = spec[0]
        if kind == "choice":
            if v in spec[1]:
                clean[k] = v
            else:
                problems.append(f"{k}: '{v}' not one of {spec[1]}")
        elif kind == "bool":
            if isinstance(v, bool):
                clean[k] = v
            else:
                problems.append(f"{k}: expected true/false")
        else:
            try:
                x = float(v)
            except (TypeError, ValueError):
                problems.append(f"{k}: not a number")
                continue
            if not (spec[1] <= x <= spec[2]):
                problems.append(f"{k}: {x} outside [{spec[1]}, {spec[2]}] (ignored)")
                continue
            clean[k] = int(round(x)) if spec[3] else x
    return clean, problems


def state_to_json(state: dict, schema: dict, note: str = "") -> str:
    keep = {k: state[k] for k in schema if k in state}
    payload = {"app": APP_TAG, "schema_version": SCHEMA_VERSION, "author": "Merouane Hamdani (educational use, MIT License)",
               "saved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"), "note": note, "state": keep}
    return json.dumps(payload, indent=2)


def json_to_state(text: str, schema: dict):
    """Parse a case file. Returns (clean_state, problems, meta). Raises ValueError for files that are not case files."""
    try:
        payload = json.loads(text)
    except json.JSONDecodeError as e:
        raise ValueError(f"Not valid JSON: {e}") from e
    if not isinstance(payload, dict) or payload.get("app") != APP_TAG or not isinstance(payload.get("state"), dict):
        raise ValueError("This file is not a SCAL Studio case file")
    clean, problems = sanitize(payload["state"], schema)
    return clean, problems, {k: payload.get(k) for k in ("schema_version", "saved_utc", "note")}
