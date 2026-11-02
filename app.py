# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
"""
SCAL Studio - relative permeability, capillary pressure, hysteresis and simulator export.

Run:  streamlit run app.py
"""
from __future__ import annotations

import io
import json
import sys
import zipfile
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import numpy as np
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

import advisor as A
import batch as B
import cases as K
import charts as ch
import core as C
import deckcheck as D
import exporters as X
import docs_content as DC
import presets as PR
import report as R
import validation as V
from state_defs import CHOICES, DEFAULTS, SYSTEMS
from animation import build_animation_html
from pyscal_bridge import PYSCAL_AVAILABLE, crosscheck_wo

st.set_page_config(page_title="SCAL Studio", page_icon="◼", layout="wide", initial_sidebar_state="expanded")

# ------------------------------------------------------------------ theme
st.markdown(
    """
<style>
 :root{--navy:#051C2C;--blue:#2251FF;--cyan:#00A9F4;--grey:#6B7C93;--soft:#F3F6FB;--line:#D5DEEA}
 html, body, [class*="css"]{font-family:Arial,Helvetica,sans-serif;color:var(--navy)}
 .block-container{padding-top:1.1rem;max-width:1400px}
 .hdr{background:var(--navy);color:#fff;padding:20px 28px;margin:-8px 0 18px 0;border-bottom:4px solid var(--blue)}
 .hdr h1{font-family:Georgia,'Times New Roman',serif;font-weight:400;font-size:30px;margin:0;letter-spacing:.2px;color:#fff}
 .hdr p{margin:4px 0 0;color:#9FB6D6;font-size:13px}
 .act{font-family:Georgia,'Times New Roman',serif;font-size:23px;line-height:1.3;color:var(--navy);margin:6px 0 2px;
      border-top:3px solid var(--navy);padding-top:10px}
 .actsub{color:var(--grey);font-size:13px;margin-bottom:10px}
 .src{color:var(--grey);font-size:11px;margin:2px 0 14px;border-top:1px solid var(--line);padding-top:4px}
 .kpis{display:flex;gap:10px;flex-wrap:wrap;margin:6px 0 14px}
 .kpi{flex:1;min-width:130px;background:var(--soft);border-left:4px solid var(--blue);padding:9px 12px}
 .kpi .v{font-family:Georgia,serif;font-size:24px;color:var(--navy);line-height:1.1}
 .kpi .l{font-size:11px;color:var(--grey);text-transform:uppercase;letter-spacing:.5px;margin-top:2px}
 .note{background:var(--soft);border-left:4px solid var(--cyan);padding:10px 14px;font-size:13px;margin:8px 0}
 section[data-testid="stSidebar"]{background:var(--soft);border-right:1px solid var(--line)}
 section[data-testid="stSidebar"] h3{font-family:Georgia,serif;font-weight:400;color:var(--navy);border-bottom:2px solid var(--navy);padding-bottom:4px}
 .stTabs [data-baseweb="tab-list"]{gap:0;border-bottom:2px solid var(--navy)}
 .stTabs [data-baseweb="tab"]{padding:10px 18px;font-weight:600;color:var(--grey);background:transparent}
 .stTabs [aria-selected="true"]{background:var(--navy);color:#fff!important}
 .stTabs [data-baseweb="tab-highlight"]{background:var(--blue)}
 .stButton>button,.stDownloadButton>button{background:var(--blue);color:#fff;border:0;border-radius:0;font-weight:700;padding:.5rem 1.1rem}
 .stButton>button:hover,.stDownloadButton>button:hover{background:#1a3fd6;color:#fff}
 h2,h3{font-family:Georgia,serif;font-weight:400}
</style>
<div class="hdr"><h1>SCAL Studio</h1><p>Relative permeability &middot; capillary pressure &middot; end-point scaling &middot; hysteresis &middot; simulator export</p><p style="margin-top:8px;color:#fff;font-size:12px;letter-spacing:.3px">Made by <b>Merouane Hamdani</b> &middot; educational purposes only &middot; MIT License</p></div>
""",
    unsafe_allow_html=True,
)


def action_title(text: str, sub: str = ""):
    st.markdown(f'<div class="act">{text}</div>' + (f'<div class="actsub">{sub}</div>' if sub else ""), unsafe_allow_html=True)


def source(text: str = "Source: SCAL Studio; Corey / LET parametrisation as in pyscal; illustrative inputs unless fitted to SCAL data"):
    st.markdown(f'<div class="src">{text}</div>', unsafe_allow_html=True)


def kpis(items):
    html = "".join(f'<div class="kpi"><div class="v">{v}</div><div class="l">{l}</div></div>' for l, v in items)
    st.markdown(f'<div class="kpis">{html}</div>', unsafe_allow_html=True)


def note(text: str):
    st.markdown(f'<div class="note">{text}</div>', unsafe_allow_html=True)


def guide(tab: str):
    """Short reading guide at the top of every section (what it shows, how to read it, what to watch)."""
    g = DC.GUIDES.get(tab)
    if g:
        with st.expander("How to read this section"):
            st.markdown(f"**What it shows.** {g['what']}\n\n**How to read it.** {g['read']}\n\n**Watch out.** {g['watch']}")


# ------------------------------------------------------------------ session defaults (widget keys)
SCHEMA = K.schema_for(DEFAULTS, CHOICES)
for _k, _v in DEFAULTS.items():
    st.session_state.setdefault(_k, _v)
    st.session_state[_k] = st.session_state[_k]   # keep values of widgets that are not rendered in this run

# ------------------------------------------------------------------ cached heavy computations
@st.cache_data(show_spinner=False, max_entries=64)
def cached_fit_pc(swd, pcd, swirr0, sorw, kname, fit_swirr):
    return C.fit_pc(swd, pcd, swirr0, sorw, kname, fit_swirr)


@st.cache_data(show_spinner=False, max_entries=32)
def cached_fit_relperm(data, kind, ep):
    return C.fit_relperm(data, kind, ep)


@st.cache_data(show_spinner="Running ensemble...", max_entries=16)
def cached_ensemble(base, spec, mu_w, mu_o, h, pc):
    return C.ensemble_wo(base, spec, mu_w, mu_o, h, pc)


@st.cache_data(show_spinner=False, max_entries=16)
def cached_flood(df, mu_w, mu_o, swi):
    return C.synthetic_core_flood(df, mu_w, mu_o, swi)


@st.cache_data(show_spinner=False, max_entries=32)
def cached_jbn(data, swi, mu_w, mu_o, kro_swi):
    return C.jbn(data, swi, mu_w, mu_o, kro_swi)


@st.cache_data(show_spinner="Fitting centrifuge data...", max_entries=16)
def cached_fit_centrifuge(rpm, sw, r2, length, drho, kind, sorw):
    return C.fit_centrifuge(rpm, sw, r2, length, drho, kind, sorw)


_PC_UNITS_TO_DECK = {"bar": "METRIC", "psi": "FIELD", "kPa": "SI-kPa", "MPa": "METRIC"}


def show(df: pd.DataFrame) -> pd.DataFrame:
    """Copy of a table for display: Pc column converted to the chosen display unit."""
    if "PC" not in df.columns:
        return df
    out = df.copy()
    out["PC"] = out["PC"] * ch.PCF
    return out.rename(columns={"PC": f"PC [{ch.PCU}]"})


def _scale_pc(p: C.PcParams, f: float) -> C.PcParams:
    return replace(p, pd=p.pd * f, pcmax=p.pcmax * f, pcmin=p.pcmin * f, cw=p.cw * f, co=p.co * f)


def _apply_wo_fit(p, k, set_system=True):
    """Callback: push water-oil curves into the sidebar widgets."""
    ss = st.session_state
    if set_system:
        ss.system = SYSTEMS[0]
    ss.kind = "Corey" if k == "corey" else "LET"
    ss.wo_swl, ss.wo_swcr, ss.wo_sorw = float(round(p.swl, 3)), float(round(p.swcr, 3)), float(round(p.sorw, 3))
    ss.wo_krwend, ss.wo_kroend = float(round(p.krwend, 3)), float(round(p.kroend, 3))
    if k == "corey":
        ss.wo_krw_n, ss.wo_kro_n = float(round(p.krw.n, 2)), float(round(p.kro.n, 2))
    else:
        ss.wo_krw_L, ss.wo_krw_E, ss.wo_krw_T = [float(round(v, 2)) for v in (p.krw.L, p.krw.E, p.krw.T)]
        ss.wo_kro_L, ss.wo_kro_E, ss.wo_kro_T = [float(round(v, 2)) for v in (p.kro.L, p.kro.E, p.kro.T)]


def _apply_advisor(profile_key, case, kind, with_pc):
    """Callback: use a benchmark case (low / base / high) of a profile as the sidebar case."""
    ss = st.session_state
    prof = A.PROFILES[profile_key]
    _apply_wo_fit(A.make_params(prof, case, kind), kind, set_system=False)
    if with_pc:
        ps = A.pc_start(prof, ss.get("adv_k", 150.0), ss.get("adv_phi", 0.22), ss.get("adv_j", 0.3))
        ss["pcwo_kind"] = "brooks-corey"
        ss["pcwo_lam"] = float(min(max(ps["lam"], 0.2), 10.0))
        ss["pcwo_pd"] = float(min(max(ps["pd_bar"], 0.0), 50.0))
        ss["pcwo_pcmax"] = float(min(max(ss.get("pcwo_pcmax", 5.0), 10 * ss["pcwo_pd"], 0.1), 200.0))


class _AdvSkip(Exception):
    """Internal: skip the advisor interpolation block when the sidebar system has no water-oil curves."""


@st.cache_data(show_spinner="Running the validation suite...", max_entries=2)
def cached_validation():
    return V.run_all()


# ------------------------------------------------------------------ case files and presets (callbacks run before widgets)
_KEEP_ON_PRESET = ("ens_", "ip_", "pc_unit", "len_unit")


def _apply_preset():
    name = st.session_state.get("preset_sel")
    if name in PR.PRESETS:
        for k, v in DEFAULTS.items():
            if not k.startswith(_KEEP_ON_PRESET):
                st.session_state[k] = v
        clean, probs = K.sanitize(PR.PRESETS[name]["state"], SCHEMA)
        st.session_state.update(clean)
        st.session_state["_case_msg"] = ("ok", f"Preset applied: {name}" + (f" ({len(probs)} value(s) skipped)" if probs else ""))


def _load_case():
    up = st.session_state.get("case_up")
    if up is None:
        return
    try:
        clean, probs, meta = K.json_to_state(up.getvalue().decode("utf-8"), SCHEMA)
    except (ValueError, UnicodeDecodeError) as exc:
        st.session_state["_case_msg"] = ("err", str(exc))
        return
    st.session_state.update(clean)
    msg = f"Case loaded ({len(clean)} settings)" + (f"; saved {meta['saved_utc']}" if meta.get("saved_utc") else "")
    if probs:
        msg += f". {len(probs)} item(s) skipped: " + "; ".join(probs[:4]) + (" ..." if len(probs) > 4 else "")
    st.session_state["_case_msg"] = ("warn" if probs else "ok", msg)


# ------------------------------------------------------------------ sidebar: case definition
sb = st.sidebar
with sb.expander("Case file, presets and units"):
    st.selectbox("Start from a preset", list(PR.PRESETS), key="preset_sel")
    _pn = st.session_state.get("preset_sel")
    if _pn in PR.PRESETS:
        st.caption(PR.PRESETS[_pn]["description"] + " Illustrative values, not from a specific field.")
    st.button("Apply preset", on_click=_apply_preset, key="preset_btn")
    st.file_uploader("Load a case file (.json)", type=["json"], key="case_up")
    st.button("Load case", on_click=_load_case, key="case_btn")
    _m = st.session_state.pop("_case_msg", None)
    if _m:
        {"ok": st.success, "warn": st.warning, "err": st.error}[_m[0]](_m[1])
    st.download_button("Save case (.json)", K.state_to_json(st.session_state, SCHEMA), "scal_case.json", "application/json", key="case_dl")
    st.caption("Case files hold the sidebar, uncertainty and interpolation settings. Edits to the fluid-configuration table are not saved.")
    st.selectbox("Pressure display unit", CHOICES["pc_unit"], key="pc_unit", help="Inputs stay in bar; plots, tables and downloads follow this unit.")
    st.selectbox("Length display unit", CHOICES["len_unit"], key="len_unit")
ch.set_units(st.session_state["pc_unit"], st.session_state["len_unit"])
sb.markdown("### Case definition")
system = sb.radio("Fluid system", SYSTEMS, key="system")
kind_label = sb.radio("Kr parametrisation", ["Corey", "LET"], horizontal=True, key="kind")
kind = kind_label.lower()
res = sb.select_slider("Table resolution (dS)", [0.02, 0.01, 0.005], key="res")
export_rows = sb.number_input("Export rows per table", 8, 1000, step=1, key="export_rows",
                              help="Exact row count in downloads. Linear interpolation retains kr plateau endpoints; plots and calculations keep the dS resolution.")
need_wo = system in (SYSTEMS[0], SYSTEMS[3])
need_go = system in (SYSTEMS[1], SYSTEMS[3])
need_gw = system == SYSTEMS[2]


def kr_widget(prefix: str, label: str) -> C.KrModel:
    if kind == "corey":
        return C.KrModel("corey", n=st.number_input(f"{label} Corey exponent", 0.5, 8.0, step=0.1, key=f"{prefix}_n"))
    c = st.columns(3)
    return C.KrModel("let", L=c[0].number_input(f"{label} L", 0.1, 10.0, step=0.1, key=f"{prefix}_L"),
                     E=c[1].number_input("E", 0.05, 50.0, step=0.1, key=f"{prefix}_E"),
                     T=c[2].number_input("T", 0.1, 10.0, step=0.1, key=f"{prefix}_T"))


def pc_widget(prefix: str, label: str) -> C.PcParams:
    k = st.selectbox(f"{label} model", C.PC_MODELS, key=f"{prefix}_kind")
    p = C.PcParams(kind=k)
    c = st.columns(2)
    if k in ("brooks-corey", "let", "van-genuchten"):
        p.pd = c[0].number_input("Entry pressure [bar]", 0.0, 50.0, step=0.05, key=f"{prefix}_pd")
    p.pcmax = c[1].number_input("Pc max [bar]", 0.1, 200.0, step=0.5, key=f"{prefix}_pcmax")
    if k == "brooks-corey":
        p.lam = st.number_input("Pore-size index lambda", 0.2, 10.0, step=0.1, key=f"{prefix}_lam")
        p.smooth_cap = st.checkbox("Smooth cap at Swirr", key=f"{prefix}_smooth", help="Curve reaches Pc max exactly at Swirr instead of a flat plateau with a kink.")
    elif k == "van-genuchten":
        p.vg_n = st.number_input("van Genuchten n", 1.2, 10.0, step=0.1, key=f"{prefix}_vg")
        p.smooth_cap = st.checkbox("Smooth cap at Swirr", key=f"{prefix}_smooth")
    elif k == "skjaeveland":
        cc = st.columns(2)
        p.cw = cc[0].number_input("cw [bar]", 0.0, 200.0, step=0.05, key=f"{prefix}_cw")
        p.aw = cc[1].number_input("aw", 0.05, 5.0, step=0.05, key=f"{prefix}_aw")
        cd = st.columns(2)
        p.co = cd[0].number_input("co [bar]", 0.0, 200.0, step=0.05, key=f"{prefix}_co")
        p.ao = cd[1].number_input("ao", 0.05, 5.0, step=0.05, key=f"{prefix}_ao")
        p.pcmin = st.number_input("Pc min (negative clip) [bar]", 0.0, 200.0, step=0.5, key=f"{prefix}_pcmin")
    else:
        cc = st.columns(3)
        p.L = cc[0].number_input("Lp", 0.1, 10.0, step=0.1, key=f"{prefix}_L")
        p.E = cc[1].number_input("Ep", 0.05, 50.0, step=0.1, key=f"{prefix}_E")
        p.T = cc[2].number_input("Tp", 0.1, 10.0, step=0.1, key=f"{prefix}_T")
        if k == "let-two-sided":
            p.pcmin = st.number_input("Pc min (negative side) [bar]", 0.0, 50.0, step=0.1, key=f"{prefix}_pcmin")
    return p


wo_p = go_p = gw_p = None
pc_wo = pc_go = pc_gw = None
try:
    if need_wo:
        with sb.expander("Water-oil end-points & curves", expanded=True):
            swl = st.slider("Swl (connate water)", 0.0, 0.5, 0.01, key="wo_swl") if False else st.slider("Swl (connate water)", 0.0, 0.5, step=0.01, key="wo_swl")
            swcr = st.slider("Swcr (critical water)", 0.0, 0.6, step=0.01, key="wo_swcr")
            sorw = st.slider("Sorw (residual oil)", 0.0, 0.6, step=0.01, key="wo_sorw")
            krwend = st.slider("krw at Sorw", 0.02, 1.0, step=0.01, key="wo_krwend")
            kroend = st.slider("kro at Swl", 0.1, 1.0, step=0.01, key="wo_kroend")
            wo_p = C.WOParams(swl, max(swcr, swl), sorw, krwend, kroend, kr_widget("wo_krw", "krw"), kr_widget("wo_kro", "kro"))
            wo_p.validate()
    if need_go:
        with sb.expander("Gas-oil end-points & curves", expanded=not need_wo):
            g_swl = wo_p.swl if need_wo else st.slider("Swl (gas-oil)", 0.0, 0.5, step=0.01, key="go_swl")
            sgcr = st.slider("Sgcr (critical gas)", 0.0, 0.4, step=0.01, key="go_sgcr")
            sorg = st.slider("Sorg (residual oil to gas)", 0.0, 0.5, step=0.01, key="go_sorg")
            krgend = st.slider("krg at Sorg", 0.05, 1.0, step=0.01, key="go_krgend")
            g_kroend = st.slider("kro at Sgcr=0", 0.1, 1.0, step=0.01, key="go_kroend")
            go_p = C.GOParams(g_swl, sgcr, sorg, krgend, g_kroend, kr_widget("go_krg", "krg"), kr_widget("go_kro", "krog"))
            go_p.validate()
    if need_gw:
        with sb.expander("Gas-water end-points & curves", expanded=True):
            w_swl = st.slider("Swl (irreducible water)", 0.0, 0.6, step=0.01, key="gw_swl")
            w_sgcr = st.slider("Sgcr (critical gas)", 0.0, 0.4, step=0.01, key="gw_sgcr")
            w_sgr = st.slider("Sgr (residual / trapped gas)", 0.0, 0.6, step=0.01, key="gw_sgr")
            w_krwend = st.slider("krw at Sgr", 0.05, 1.0, step=0.01, key="gw_krwend")
            w_krgend = st.slider("krg end-point", 0.05, 1.0, step=0.01, key="gw_krgend")
            gw_p = C.GWParams(w_swl, w_sgcr, w_sgr, w_krwend, w_krgend, kr_widget("gw_krw", "krw"), kr_widget("gw_krg", "krg"))
            gw_p.validate()
except ValueError as e:
    st.error(f"Invalid end-points: {e}")
    st.stop()

with sb.expander("Fluid properties"):
    mu_w = st.number_input("Water viscosity [cP]", 0.1, 10.0, step=0.1, key="mu_w")
    mu_o = st.number_input("Oil viscosity [cP]", 0.1, 500.0, step=0.5, key="mu_o")

pc_on = sb.checkbox("Include capillary pressure in tables", key="pc_on")
if pc_on:
    with sb.expander("Capillary pressure models", expanded=False):
        if need_wo:
            st.caption("Oil-water Pc (Pcow = Po - Pw)")
            pc_wo = pc_widget("pcwo", "Pcow")
        if need_go:
            st.caption("Gas-oil Pc (Pcgo = Pg - Po)")
            pc_go = pc_widget("pcgo", "Pcgo")
        if need_gw:
            st.caption("Gas-water Pc (Pcgw = Pg - Pw)")
            pc_gw = pc_widget("pcgw", "Pcgw")

if PYSCAL_AVAILABLE:
    sb.success("pyscal detected - cross-check available")
else:
    sb.caption("pyscal not installed - native engine active (same Corey / LET definitions). `pip install pyscal` enables the cross-check.")

# ------------------------------------------------------------------ build tables
try:
    wo_df = go_df = gw_df = None
    if need_wo:
        wo_df = C.wateroil(wo_p, res)
        wo_df = C.add_pc_wo(wo_df, wo_p, pc_wo) if pc_on else wo_df.assign(PC=0.0)
    if need_go:
        go_df = C.gasoil(go_p, res)
        go_df = C.add_pc_go(go_df, go_p, pc_go) if pc_on else go_df.assign(PC=0.0)
    if need_gw:
        gw_df = C.gaswater(gw_p, res)
        gw_df = C.add_pc_gw(gw_df, gw_p, pc_gw) if pc_on else gw_df.assign(PC=0.0)
except ValueError as e:
    st.error(str(e))
    st.stop()

wo_w = C.welge(wo_df, mu_w, mu_o, swi=wo_p.swl) if need_wo else None
base_region = {"name": "Base case"}
if wo_df is not None:
    base_region["wo"] = wo_df
if go_df is not None:
    base_region["go"] = go_df
if gw_df is not None:
    base_region["gw"] = gw_df
base_region["provenance"] = {
    key: X.model_provenance(p, pc, resolution=res)
    for key, p, pc in (("wo", wo_p, pc_wo), ("go", go_p, pc_go), ("gw", gw_p, pc_gw))
    if p is not None
}


def wo_region(name, table, params, pc_scale=1.0, source="Derived water-oil case"):
    reg = dict(name=name, wo=table, provenance={
        "wo": X.model_provenance(params, pc_wo, resolution=res, pc_scale=pc_scale, source=source),
    })
    if go_df is not None:
        reg["go"] = go_df
        reg["provenance"]["go"] = base_region["provenance"]["go"]
    return reg

hyst_export = None
region_sources = {"Base case": [base_region]}
ens_result = None

TABS = ["Summary", "Relative permeability", "Capillary pressure", "End-point scaling", "Parameter advisor", "Uncertainty", "Curve-set interpolation",
        "Hysteresis", "Lab interpretation", "Fit SCAL data", "Animations", "Batch mode", "Validation", "Export and report", "Documentation"]
T = dict(zip(TABS, st.tabs(TABS)))

# ================================================================== 0 Summary
with T["Summary"]:
    guide("Summary")
    if need_wo:
        w = wo_w
        verdict = "unfavourable" if w["mobility_ratio"] > 1 else "favourable"
        action_title(
            f"Waterflood can displace at most {w['rf_ult']:.0%} of the initial oil; with a mobility ratio of "
            f"{w['mobility_ratio']:.1f} the front is {verdict} and water breaks through after {w['pvi_bt']:.2f} PV",
            "Displacement efficiency from the SCAL end-points and Buckley-Leverett / Welge analysis (no capillary or gravity effects)")
        kpis([("Connate water Swl", f"{wo_p.swl:.2f}"), ("Residual oil Sorw", f"{wo_p.sorw:.2f}"),
              ("Krw-Kro crossover Sw", f"{w['sw_cross']:.2f}"), ("End-point mobility ratio", f"{w['mobility_ratio']:.1f}"),
              ("Breakthrough [PV]", f"{w['pvi_bt']:.2f}"), ("RF at breakthrough", f"{w['rf_bt']:.0%}"),
              ("Displacement limit", f"{w['rf_ult']:.0%}")])
        c1, c2 = st.columns(2)
        c1.plotly_chart(ch.fig_wo(wo_df), use_container_width=True)
        c2.plotly_chart(ch.fig_recovery(C.recovery_curve(w, 3.0), w), use_container_width=True)
        source()
        cross = w["sw_cross"]
        lean = "water-wet" if cross > 0.55 else ("oil-wet" if cross < 0.45 else "mixed / intermediate-wet")
        extra = ("Consider mobility control (polymer, WAG) - the unfavourable ratio will leave bypassed oil." if w["mobility_ratio"] > 1
                 else "A favourable mobility ratio supports stable, piston-like displacement.")
        note(f"<b>Wettability read-across (Craig rules of thumb, indicative only):</b> crossover at Sw = {cross:.2f} points to a "
             f"{lean} system. {extra}")
    if need_go:
        t = go_df
        action_title(f"Gas injection is limited by a critical gas saturation of {go_p.sgcr:.2f} and {go_p.sorg:.0%} residual oil to gas",
                     "Gas-oil end-points and curve shapes drive gas breakthrough and oil trapping")
        kpis([("Sgcr", f"{go_p.sgcr:.2f}"), ("Sorg", f"{go_p.sorg:.2f}"), ("krg end-point", f"{go_p.krgend:.2f}"),
              ("Max mobile gas sat.", f"{1 - go_p.swl - go_p.sorg:.2f}")])
        st.plotly_chart(ch.fig_go(t), use_container_width=True)
        source()
    if need_gw:
        t = gw_df
        action_title(f"Residual gas trapping of {gw_p.sgr:.0%} bounds the storage security of the gas-water system",
                     "Relevant for gas storage, CO2 and H2 injection: end-points control plume migration and trapping")
        kpis([("Swl", f"{gw_p.swl:.2f}"), ("Sgcr", f"{gw_p.sgcr:.2f}"), ("Residual gas Sgr", f"{gw_p.sgr:.2f}"),
              ("krg end-point", f"{gw_p.krgend:.2f}"), ("krw at Sgr", f"{gw_p.krwend:.2f}")])
        st.plotly_chart(ch.fig_gw(t), use_container_width=True)
        source()
    if PYSCAL_AVAILABLE and need_wo:
        diff, msg = crosscheck_wo(wo_p, res)
        (st.success if diff is not None and diff < 1e-6 else st.info)(msg)

# ================================================================== 1 Relative permeability
with T["Relative permeability"]:
    guide("Relative permeability")
    action_title("Curve shape, not only end-points, decides how quickly mobile phases reach the producer",
                 "Interactive relative permeability curves; use the log scale to inspect the low-kr tail")
    logy = st.toggle("Logarithmic kr axis", value=False)
    if need_wo:
        st.markdown("**Water-oil**")
        st.plotly_chart(ch.fig_wo(wo_df, log=logy), use_container_width=True, key="relperm_wo_chart")
        with st.expander("Table"):
            st.dataframe(show(wo_df), use_container_width=True, height=240)
    if need_go:
        st.markdown("**Gas-oil**")
        st.plotly_chart(ch.fig_go(go_df, log=logy), use_container_width=True, key="relperm_go_chart")
        with st.expander("Table"):
            st.dataframe(show(go_df), use_container_width=True, height=240)
    if need_gw:
        st.markdown("**Gas-water**")
        st.plotly_chart(ch.fig_gw(gw_df, log=logy), use_container_width=True, key="relperm_gw_chart")
        with st.expander("Table"):
            st.dataframe(show(gw_df), use_container_width=True, height=240)
    if need_wo and need_go:
        st.markdown("**Three-phase oil relative permeability**")
        model = st.radio("Three-phase model", ["stone2", "stone1", "baker"], horizontal=True,
                         format_func={"stone2": "Stone II", "stone1": "Stone I", "baker": "Baker (linear)"}.get)
        sw_g, sg_g = np.meshgrid(np.arange(wo_p.swl, 1.0001, 0.025), np.arange(0, 1.0 - wo_p.swl + 1e-9, 0.025))
        sw_f, sg_f = sw_g.ravel(), sg_g.ravel()
        keep = (1 - sw_f - sg_f) >= -1e-9
        sw_f, sg_f = sw_f[keep], sg_f[keep]
        kro3 = C.three_phase_kro(wo_df, go_df, wo_p.swl, sw_f, sg_f, model)
        st.plotly_chart(ch.fig_ternary(sw_f, sg_f, kro3), use_container_width=True)
        source("Source: SCAL Studio; three-phase kro from two-phase WO and GO curves")
    if PYSCAL_AVAILABLE and need_wo:
        st.caption("pyscal cross-check is shown on the Summary tab.")

# ================================================================== 2 Capillary pressure
def _apply_pc_fit(kind, fitted, s_fit):
    """Callback: push a fitted Pc model into the sidebar widgets (runs before the next rerun)."""
    ss = st.session_state
    pre = "pcwo" if need_wo else ("pcgw" if need_gw else "pcgo")
    clip = lambda v, lo, hi: float(min(max(v, lo), hi))  # noqa: E731
    ss[f"{pre}_kind"] = kind
    ss[f"{pre}_pcmax"] = clip(fitted.pcmax, 0.1, 200.0)
    ss[f"{pre}_pcmin"] = clip(fitted.pcmin, 0.0, 200.0)
    ss[f"{pre}_pd"] = clip(fitted.pd, 0.0, 50.0)
    ss[f"{pre}_lam"] = clip(fitted.lam, 0.2, 10.0)
    ss[f"{pre}_vg"] = clip(fitted.vg_n, 1.2, 10.0)
    ss[f"{pre}_L"], ss[f"{pre}_E"], ss[f"{pre}_T"] = clip(fitted.L, 0.1, 10.0), clip(fitted.E, 0.05, 50.0), clip(fitted.T, 0.1, 10.0)
    ss[f"{pre}_cw"], ss[f"{pre}_aw"] = clip(fitted.cw, 0.0, 200.0), clip(fitted.aw, 0.05, 5.0)
    ss[f"{pre}_co"], ss[f"{pre}_ao"] = clip(fitted.co, 0.0, 200.0), clip(fitted.ao, 0.05, 5.0)
    ss[f"{pre}_smooth"] = False  # fitted as a pure function
    if ss.get("pc_set_swl") and pre == "pcwo":
        ss.wo_swl = clip(round(s_fit, 3), 0.0, 0.5)


with T["Capillary pressure"]:
    guide("Capillary pressure")
    action_title("Capillary pressure scales with interfacial tension and contact angle: one rock, several fluid configurations",
                 "Define or fit Pc, translate it between fluid systems, and read it as saturation versus height above the free-water level")
    if not pc_on:
        note("Capillary pressure is switched off in the sidebar. Enable it to generate Pc curves.")
    else:
        views = {}
        if need_wo:
            views["Oil-water (Pcow vs Sw)"] = (wo_df.SW.values, wo_df.PC.values, "Sw", "Oil-water (reservoir)", pc_wo, wo_p.swl, wo_p.sorw)
        if need_go:
            views["Gas-oil (Pcgo vs Sg)"] = (go_df.SG.values, go_df.PC.values, "Sg", "Gas-oil (reservoir)", pc_go, go_p.swl, 0.0)
        if need_gw:
            views["Gas-water (Pcgw vs Sw)"] = (gw_df.SW.values, gw_df.PC.values, "Sw", "Gas-water (reservoir)", pc_gw, gw_p.swl, gw_p.sgr)
        vsel = st.selectbox("Pc curve to analyse", list(views))
        x, pc0, xname, default_sys, view_pc, view_swirr, view_sorw = views[vsel]

        st.markdown("**Fluid configurations** (interfacial tension and contact angle are editable typical values)")
        fs_df = pd.DataFrame([dict(system=k, sigma_mNm=v["sigma"], theta_deg=v["theta"]) for k, v in C.FLUID_SYSTEMS.items()])
        fs_ed = st.data_editor(fs_df, hide_index=True, use_container_width=True, key="fs_editor",
                               column_config={"system": st.column_config.TextColumn(disabled=True),
                                              "sigma_mNm": st.column_config.NumberColumn("sigma [mN/m]", min_value=0.1),
                                              "theta_deg": st.column_config.NumberColumn("theta [deg]", min_value=0.0, max_value=180.0)})
        FS = {r.system: dict(sigma=float(r.sigma_mNm), theta=float(r.theta_deg)) for r in fs_ed.itertuples()}
        c1, c2 = st.columns(2)
        base_sys = c1.selectbox("Curve is defined for", list(FS), index=list(FS).index(default_sys))
        targets = c2.multiselect("Translate to", [s for s in FS if s != base_sys],
                                 default=[s for s in ("CO2-brine (storage)", "Air-brine (lab)") if s != base_sys][:2])
        rock = st.checkbox("Also scale for a different rock (Leverett sqrt(phi/k))", value=False)
        k1 = p1 = k2 = p2 = None
        if rock:
            r1, r2, r3, r4 = st.columns(4)
            k1 = r1.number_input("k reference [mD]", 0.1, 20000.0, 150.0)
            p1 = r2.number_input("phi reference", 0.02, 0.45, 0.22)
            k2 = r3.number_input("k target [mD]", 0.1, 20000.0, 150.0)
            p2 = r4.number_input("phi target", 0.02, 0.45, 0.22)
        curves = {f"{base_sys} (base)": (x, pc0)}
        conv = {}
        for tsys in targets:
            y = C.pc_convert(pc0, FS[base_sys], FS[tsys], k1, p1, k2, p2)
            curves[tsys] = (x, y, "dash")
            conv[tsys] = y
        st.plotly_chart(ch.fig_pc(curves, xt=f"{xname} [-]", hline0=need_wo and pc_wo is not None and pc_wo.kind == "let-two-sided"),
                        use_container_width=True)
        source("Source: SCAL Studio; Pc2 = Pc1 * (sigma cos theta)2 / (sigma cos theta)1 [* sqrt((phi2/k2)/(phi1/k1))]")

        # ---- quality checks
        xs_c, ys_c = (x, pc0) if xname == "Sw" else ((1 - x)[::-1], pc0[::-1])
        issues = C.pc_checks(xs_c, ys_c, view_pc)
        if issues:
            for msg in issues:
                st.warning(msg)
        else:
            st.success("Pc checks passed: monotone, sensible range, consistent with the selected model.")

        t1, t2 = st.columns(2)
        with t1:
            st.markdown("**Saturation-height function**")
            if xname != "Sw":
                st.info("Saturation-height is shown for water-wetting systems (oil-water, gas-water).")
            else:
                insitu = st.selectbox("In-situ fluid system", list(FS), index=list(FS).index(base_sys), key="shf_sys")
                drho = st.number_input("Density difference [kg/m3]", 20.0, 1200.0, 250.0 if "Gas" not in base_sys else 600.0, 10.0, key="shf_drho")
                hmax = st.number_input("Max height above FWL [m]", 10.0, 1000.0, 150.0, 10.0, key="shf_h")
                scale = abs(C.sigma_cos(FS[insitu])) / abs(C.sigma_cos(FS[base_sys]))
                hh = np.linspace(0, hmax, 200)
                swh = C.saturation_height(hh, drho, view_swirr, view_sorw, view_pc, scale)
                st.plotly_chart(ch.fig_shf(swh, hh), use_container_width=True)
                pc_at = lambda sv: max(float(C.pc_curve(sv, view_swirr, view_sorw, view_pc)) * scale, 0.0)  # noqa: E731
                kpis([(f"Entry height [{ch.LENU}]", f"{float(C.height_above_fwl(pc_at(1.0), drho)) * ch.LENF:.1f}"),
                      (f"Height at Sw = 50% [{ch.LENU}]", f"{float(C.height_above_fwl(pc_at(0.5), drho)) * ch.LENF:.1f}")])
        with t2:
            st.markdown("**Leverett J-function collapses the configurations**")
            kj = st.number_input("k for J-function [mD]", 0.1, 20000.0, 150.0, key="kj")
            pj = st.number_input("phi for J-function", 0.02, 0.45, 0.22, key="pj")
            jc = {base_sys: (x, C.j_function(np.clip(pc0, 1e-9, None), FS[base_sys], kj, pj))}
            for tsys, y in conv.items():
                jc[tsys] = (x, C.j_function(np.clip(y, 1e-9, None), FS[tsys], kj if not rock else k2, pj if not rock else p2))
            st.plotly_chart(ch.fig_j(jc), use_container_width=True)
        if need_wo and pc_wo is not None and pc_wo.kind in ("let-two-sided", "skjaeveland"):
            wi = C.wettability_indices(wo_df.SW.values, wo_df.PC.values)
            kpis([("USBM index (indicative)", f"{wi['usbm']:.2f}"), ("Amott-type water index", f"{wi['amott']:.2f}")])
        out = pd.DataFrame({xname: x, f"Pc {base_sys} [{ch.PCU}]": pc0 * ch.PCF, **{f"Pc {k} [{ch.PCU}]": v * ch.PCF for k, v in conv.items()}})
        st.download_button("Download Pc table (CSV)", out.to_csv(index=False).encode(), "pc_configurations.csv", "text/csv")

        # ---- seal capacity
        with st.expander("Seal capacity for gas, CO2 or H2 storage"):
            st.caption("A caprock holds a static column h = Pc_entry / (drho g). Enter the entry (threshold) pressure measured with "
                       "one fluid pair and convert it to the storage fluid pair.")
            q1, q2, q3, q4 = st.columns(4)
            seal_pc = q1.number_input("Seal entry pressure [bar]", 0.1, 500.0, 8.0, 0.5, key="seal_pc")
            seal_from = q2.selectbox("Measured with", list(FS), index=list(FS).index("Mercury-air (lab, MICP)"), key="seal_from")
            seal_to = q3.selectbox("Storage fluid pair", list(FS), index=list(FS).index("CO2-brine (storage)"), key="seal_to")
            seal_drho = q4.number_input("Density difference [kg/m3]", 20.0, 1000.0, 400.0, 10.0, key="seal_drho")
            pc_store = float(C.pc_convert(seal_pc, FS[seal_from], FS[seal_to]))
            hcol = C.seal_column_height(pc_store, seal_drho)
            sf = st.slider("Safety factor", 1.0, 5.0, 2.0, 0.1, key="seal_sf")
            kpis([(f"Entry pressure, storage pair [{ch.PCU}]", f"{pc_store * ch.PCF:.2f}"), (f"Max static column [{ch.LENU}]", f"{hcol * ch.LENF:.0f}"),
                  (f"Column with safety factor {sf:.1f} [{ch.LENU}]", f"{hcol * ch.LENF / sf:.0f}")])

        # ---- Pc fitting
        st.markdown("---")
        st.markdown("**Fit capillary pressure to measured data (MICP, centrifuge, porous plate)**")
        f1, f2, f3 = st.columns(3)
        data_kind = f1.radio("Data type", ["Mercury injection (MICP)", "Pc vs Sw"], key="pcfit_kind_sel")
        p_unit = f2.selectbox("Pressure unit", ["bar", "psi", "MPa"], key="pcfit_unit")
        tgt_sys = f3.selectbox("Convert to fluid system", list(FS), index=list(FS).index(default_sys), key="pcfit_tgt")
        micp_demo = C.demo_micp()
        st.download_button("Download demo MICP data (CSV)", micp_demo.to_csv(index=False).encode(), "demo_micp.csv", "text/csv")
        up_pc = st.file_uploader("Upload CSV (MICP: pressure, Hg saturation; or Sw, Pc)", type=["csv"], key="pc_up")
        raw = None
        if up_pc is not None:
            raw = pd.read_csv(up_pc)
            raw.columns = [c.strip().upper() for c in raw.columns]
            raw = raw.rename(columns={"PRESSURE": "P_HG_BAR", "P": "P_HG_BAR", "SHG": "S_HG", "HG_SAT": "S_HG", "PCOW": "PC"})
        elif st.checkbox("Use demo data", value=True, key="pc_demo"):
            raw = micp_demo
        if raw is not None:
            to_bar = {"bar": 1.0, "psi": 1.0 / C.BAR_TO_PSI, "MPa": 10.0}[p_unit]
            try:
                if data_kind.startswith("Mercury"):
                    if not {"P_HG_BAR", "S_HG"} <= set(raw.columns):
                        raise ValueError("MICP data need columns pressure (P_HG_BAR) and Hg saturation (S_HG)")
                    shg = raw.S_HG.values / (100.0 if raw.S_HG.max() > 1.5 else 1.0)
                    sconf = st.slider("Closure / conformance correction (fraction of pore volume)", 0.0, 0.10, 0.0, 0.005, key="pc_sconf")
                    swd, pcd = C.micp_to_pc(raw.P_HG_BAR.values * to_bar, shg, FS[tgt_sys], sconf)
                    sor_fit = 0.0
                else:
                    if not {"SW", "PC"} <= set(raw.columns):
                        if {"P_HG_BAR", "S_HG"} <= set(raw.columns):  # demo path
                            swd, pcd = C.micp_to_pc(raw.P_HG_BAR.values, raw.S_HG.values, FS[tgt_sys], 0.0)
                        else:
                            raise ValueError("Columns needed: Sw and Pc")
                    else:
                        o = np.argsort(raw.SW.values)
                        swd, pcd = raw.SW.values[o], raw.PC.values[o] * to_bar
                    sor_fit = view_sorw
                m1, m2 = st.columns(2)
                models = m1.multiselect("Models to fit", C.PC_MODELS[:-1] + ["let-two-sided"], default=["brooks-corey", "let", "skjaeveland"], key="pc_models")
                fit_s = m2.checkbox("Fit Swirr as well (recommended for MICP)", value=True, key="pc_fit_swirr")
                rows, fcurves, fits = [], {}, {}
                for kname in models:
                    try:
                        fp, rm, r2, s0 = cached_fit_pc(swd, pcd, float(swd.min()), sor_fit, kname, fit_s)
                    except Exception as exc:  # noqa: BLE001
                        st.warning(f"{kname}: fit failed ({exc})")
                        continue
                    fits[kname] = (fp, s0)
                    grid = np.linspace(max(s0, 0.0) + 1e-3, 1.0 - (sor_fit if kname in ("let-two-sided", "skjaeveland") else 0.0), 200)
                    fcurves[f"{kname}  (R2 {r2:.3f})"] = (grid, C.pc_curve(grid, s0, sor_fit, fp))
                    par = {"brooks-corey": f"pd={fp.pd:.3f}, lambda={fp.lam:.2f}", "van-genuchten": f"pd={fp.pd:.3f}, n={fp.vg_n:.2f}",
                           "let": f"pd={fp.pd:.3f}, pcmax={fp.pcmax:.2f}, L/E/T={fp.L:.1f}/{fp.E:.1f}/{fp.T:.1f}",
                           "skjaeveland": f"cw={fp.cw:.3f}, aw={fp.aw:.2f}, co={fp.co:.3f}, ao={fp.ao:.2f}",
                           "let-two-sided": f"pcmax={fp.pcmax:.2f}, pcmin={fp.pcmin:.2f}"}[kname]
                    rows.append(dict(model=kname, R2=round(r2, 3), RMSE_bar=round(rm, 3), Swirr=round(s0, 3), parameters=par))
                if rows:
                    best = min(rows, key=lambda r: r["RMSE_bar"])
                    logfit = st.toggle("Logarithmic Pc axis", value=True, key="pc_logfit")
                    st.plotly_chart(ch.fig_pc_fit(swd, pcd, fcurves, logfit), use_container_width=True)
                    st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
                    note(f"<b>Best fit by RMSE:</b> {best['model']} (R2 {best['R2']}, Swirr {best['Swirr']}). "
                         "Fitting Swirr removes the bias caused by treating the lowest measured Sw as irreducible.")
                    st.checkbox("Also set Swl to the fitted Swirr (water-oil)", value=False, key="pc_set_swl")
                    pick = st.selectbox("Apply this model to the case", [r["model"] for r in rows], index=[r["model"] for r in rows].index(best["model"]))
                    st.button("Use fitted Pc in the sidebar", on_click=_apply_pc_fit, args=(pick, fits[pick][0], fits[pick][1]))
            except ValueError as exc:
                st.error(str(exc))

# ================================================================== 3 End-point scaling
with T["End-point scaling"]:
    guide("End-point scaling")
    action_title("End-point scaling turns one SCAL curve set into a family that follows rock quality",
                 "Re-scale the base curves to new end-points, or generate permeability-ranked rock types for SATNUM regions")
    if not need_wo:
        note("End-point scaling and rock-type families are available for systems that include water-oil. Select Water-oil or Three-phase in the sidebar.")
    else:
        st.markdown("**A. Re-scale the base curve to new end-points**")
        a, b, c = st.columns(3)
        n_swl = a.slider("New Swl", 0.0, 0.5, step=0.01, key="sc_swl")
        n_swcr = a.slider("New Swcr", 0.0, 0.6, step=0.01, key="sc_swcr")
        n_sorw = b.slider("New Sorw", 0.0, 0.6, step=0.01, key="sc_sorw")
        n_krwend = b.slider("New krw at Sorw", 0.02, 1.0, step=0.01, key="sc_krwend")
        n_kroend = c.slider("New kro at Swl", 0.1, 1.0, step=0.01, key="sc_kroend")
        mode = c.radio("Scaling", ["two-point", "three-point"], horizontal=True)
        krwr = krwmax = None
        if mode == "three-point":
            krwr = c.slider("krw at Sorw (3-pt)", 0.02, 1.0, step=0.01, key="sc_krwr")
            krwmax = c.slider("krw at Sw = 1 (3-pt)", 0.1, 1.0, step=0.01, key="sc_krwmax")
        pcf = st.number_input("Pc multiplier (e.g. sqrt(phi/k) change)", 0.1, 10.0, 1.0, step=0.05)
        try:
            new_p = C.WOParams(n_swl, max(n_swcr, n_swl), n_sorw, n_krwend, n_kroend, wo_p.krw, wo_p.kro)
            new_p.validate()
            scaled = C.rescale_wo_table(wo_df, wo_p, new_p, mode, krwr_new=krwr, krwmax_new=krwmax or 1.0, pc_factor=pcf, h=res)
            ws = C.welge(scaled, mu_w, mu_o, swi=new_p.swl)
            action_title(f"Re-scaled curves move the displacement limit from {wo_w['rf_ult']:.0%} to {ws['rf_ult']:.0%} "
                         f"and breakthrough from {wo_w['pvi_bt']:.2f} to {ws['pvi_bt']:.2f} PV")
            st.plotly_chart(ch.fig_wo(scaled, ref=wo_df, title_ref="base"), use_container_width=True)
            if pc_on:
                st.plotly_chart(ch.fig_pc({"Base Pc": (wo_df.SW, wo_df.PC), "Scaled Pc": (scaled.SW, scaled.PC, "dash")}),
                                use_container_width=True)
            st.session_state["scaled_wo"] = scaled
            st.download_button("Download scaled table (CSV)", scaled.to_csv(index=False).encode(), "scaled_wo.csv", "text/csv")
        except ValueError as e:
            st.error(f"Invalid scaled end-points: {e}")

        st.markdown("---")
        st.markdown("**B. Rock-type family (SATNUM regions) from permeability range**")
        st.caption("The current base case is the good-quality rock. Define the poor-quality rock; end-points and exponents are interpolated in ln(k).")
        r1, r2, r3 = st.columns(3)
        n_rt = r1.slider("Number of rock types", 2, 8, 4)
        k_low = r1.number_input("k poor rock [mD]", 0.1, 5000.0, 5.0)
        k_high = r1.number_input("k good rock [mD]", 1.0, 20000.0, 500.0)
        phi_low = r2.number_input("phi poor rock", 0.02, 0.4, 0.12)
        phi_high = r2.number_input("phi good rock", 0.05, 0.45, 0.28)
        swl_poor = r2.slider("Swl poor rock", 0.0, 0.6, min(wo_p.swl + 0.12, 0.6), 0.01)
        sorw_poor = r3.slider("Sorw poor rock", 0.0, 0.6, min(wo_p.sorw + 0.08, 0.6), 0.01)
        krw_poor = r3.slider("krw end-point poor rock", 0.02, 1.0, max(wo_p.krwend * 0.5, 0.05), 0.01)
        kro_poor = r3.slider("kro end-point poor rock", 0.1, 1.0, max(wo_p.kroend * 0.8, 0.1), 0.01)
        try:
            poor = replace(wo_p, swl=swl_poor, swcr=max(swl_poor, wo_p.swcr), sorw=sorw_poor, krwend=krw_poor, kroend=kro_poor)
            poor.validate()
            fam = C.rock_type_family(poor, wo_p, n_rt, k_low, k_high, phi_low, phi_high)
            fam_rows = []
            for r in fam:
                t = C.wateroil(r["params"], res)
                t = C.add_pc_wo(t, r["params"], pc_wo, scale=r["pc_scale"]) if pc_on else t.assign(PC=0.0)
                r["table"] = t
                p = r["params"]
                fam_rows.append(dict(RT=r["id"], k_mD=round(r["k"], 1), phi=round(r["phi"], 3), Swl=round(p.swl, 3), Sorw=round(p.sorw, 3),
                                     krw_end=round(p.krwend, 3), kro_end=round(p.kroend, 3), Pc_scale=round(r["pc_scale"], 2)))
            st.plotly_chart(ch.fig_family(fam), use_container_width=True)
            st.dataframe(pd.DataFrame(fam_rows), hide_index=True, use_container_width=True)
            fam_regs = []
            for r in fam:
                reg = wo_region(f"RT{r['id']} k={r['k']:.0f}mD", r["table"], r["params"], r["pc_scale"],
                                source=f"Rock-type family; log(k) parameter interpolation; k_mD={r['k']!r}; phi={r['phi']!r}")
                fam_regs.append(reg)
            region_sources["Rock-type family (SATNUM regions)"] = fam_regs
            st.caption("The family is available as a region source on the Export tab.")
        except ValueError as e:
            st.error(f"Invalid rock-type end-points: {e}")

# ================================================================== 4 Hysteresis
with T["Hysteresis"]:
    guide("Hysteresis")
    action_title("Hysteresis traps the non-wetting phase: the scanning curve you follow decides how much gas or oil stays behind",
                 "Land trapping with Killough or Carlson scanning curves; Killough-type interpolation for capillary pressure")
    options = {}
    if need_go:
        options["Gas-oil: gas (non-wetting)"] = ("go", go_df.SG.values, go_df.KRG.values, go_p.sgcr, 1 - go_p.swl - go_p.sorg)
    if need_gw:
        o = np.argsort(1 - gw_df.SW.values)
        options["Gas-water: gas (non-wetting)"] = ("gw", (1 - gw_df.SW.values)[o], gw_df.KRG.values[o], gw_p.sgcr, 1 - gw_p.swl)
    if need_wo:
        o = np.argsort(1 - wo_df.SW.values)
        options["Water-oil: oil (non-wetting)"] = ("wo", (1 - wo_df.SW.values)[o], wo_df.KROW.values[o], wo_p.sorw, 1 - wo_p.swl)
    sel = st.selectbox("Non-wetting phase", list(options))
    tag, s_arr, kr_arr, s_crit, s_max = options[sel]
    h1, h2, h3 = st.columns(3)
    scheme = h1.radio("Scanning-curve scheme", ["killough", "carlson"], format_func=str.capitalize, horizontal=True)
    n_scans = h2.slider("Number of scanning curves", 2, 8, 5)
    lo_t, hi_t = s_crit + 0.02, s_max - 0.02
    trap_default = float(np.clip(s_crit + 0.45 * (s_max - s_crit), lo_t, hi_t))
    s_trap_max = h3.slider("Maximum trapped saturation", float(lo_t), float(hi_t), trap_default, 0.01, key=f"hy_trap_{tag}")
    hs = C.HystSetup(s_arr, kr_arr, float(s_crit), float(s_max), float(s_trap_max))
    fam_h = C.hysteresis_family(hs, scheme, n_scans)
    st.plotly_chart(ch.fig_hyst(fam_h), use_container_width=True)
    source("Source: SCAL Studio; Land (1968) trapping; Killough (1976) / Carlson (1981)-type scanning curves")
    cA, cB = st.columns(2)
    with cA:
        sh = np.linspace(s_crit + 0.01, s_max, 40)
        trap = [C.land_trapped(v, s_crit, s_trap_max, s_max) for v in sh]
        f = ch.fig_pc({"Land trapped saturation": (sh, trap), "No-hysteresis limit (1:1)": (sh, sh, "dot")},
                      xt="Reversal saturation [-]", yt="Trapped saturation [-]", hline0=False)
        st.plotly_chart(f, use_container_width=True)
    with cB:
        st.dataframe(pd.DataFrame([dict(reversal=round(s["s_hy"], 3), trapped=round(s["s_trap"], 3),
                                        mobile_share=round(1 - s["s_trap"] / s["s_hy"], 2)) for s in fam_h["scans"]]),
                     hide_index=True, use_container_width=True)
    if tag == "go":
        imb = C.imbibition_gas_table(go_df, hs, scheme)
        imb["PC"] = go_df.PC.values
        hyst_export = dict(scheme=scheme, s_trap_max=float(s_trap_max), imb_go=imb)
        st.caption("The imbibition gas table is added to the Eclipse export when hysteresis export is enabled on the Export tab.")

    st.markdown("---")
    st.markdown("**Capillary-pressure hysteresis (Killough-type scanning)**")
    if need_wo and pc_wo is not None:
        q1, q2, q3 = st.columns(3)
        imb_pcmax = q1.number_input("Imbibition Pc max [bar]", 0.1, 50.0, max(pc_wo.pcmax * 0.5, 0.1), 0.1)
        imb_pcmin = q1.number_input("Imbibition Pc min (negative) [bar]", 0.0, 20.0, 1.0, 0.1)
        eps = q2.slider("Curvature parameter", 0.02, 0.5, 0.1, 0.01)
        sw_hy = q2.slider("Reversal Sw", float(wo_p.swl + 0.05), float(1 - wo_p.sorw - 0.05), float(wo_p.swl + 0.3), 0.01)
        direction = q3.radio("Scanning direction", ["up", "down"], format_func={"up": "Imbibition (Sw up)", "down": "Drainage (Sw down)"}.get)
        swv = np.linspace(wo_p.swl, 1 - wo_p.sorw, 120)
        pdr = C.pc_curve(swv, wo_p.swl, wo_p.sorw, pc_wo) if pc_wo.kind != "let-two-sided" else C.pc_curve(swv, wo_p.swl, wo_p.sorw, C.PcParams(pd=0.3, pcmax=pc_wo.pcmax))
        pim = C.pc_curve(swv, wo_p.swl, wo_p.sorw, C.PcParams(kind="let-two-sided", pcmax=imb_pcmax, pcmin=imb_pcmin))
        if direction == "up":
            seg = swv >= sw_hy
            sc = C.pc_scanning(swv[seg], pdr[seg], pim[seg], sw_hy, swv[-1], "up", eps)
        else:
            seg = swv <= sw_hy
            sc = C.pc_scanning(swv[seg], pdr[seg], pim[seg], sw_hy, swv[0], "down", eps)
        figp = ch.fig_pc({"Primary drainage": (swv, pdr), "Boundary imbibition": (swv, pim, "dash"), "Scanning curve": (swv[seg], sc)})
        st.plotly_chart(figp, use_container_width=True)
    else:
        note("Pc hysteresis is shown for systems that include water-oil with capillary pressure switched on.")


# ================================================================== Parameter advisor
with T["Parameter advisor"]:
    guide("Parameter advisor")
    action_title("Start from a benchmark range for your rock type, then interpolate between the low, base and high case to export one specific curve",
                 "Generic ranges (Craig 1971 rules, Lomeland et al. 2005 LET guidance, judgement); replace them with your own SCAL data whenever you have it")
    pk = st.selectbox("Rock and wettability profile", list(A.PROFILES), format_func=lambda k: A.PROFILES[k].title, key="adv_profile")
    prof = A.PROFILES[pk]
    st.caption(prof.description)
    kind_a = st.radio("Curve shape", ["corey", "let"], horizontal=True, format_func=str.upper, key="adv_kind_radio")

    cur = {}
    if need_wo and wo_p is not None:
        cur = dict(swl=wo_p.swl, sorw=wo_p.sorw, krwend=wo_p.krwend, kroend=wo_p.kroend,
                   nw=wo_p.krw.n if wo_p.krw.kind == "corey" else float("nan"), no=wo_p.kro.n if wo_p.kro.kind == "corey" else float("nan"))
    g1, g2 = st.columns(2)
    g1.plotly_chart(ch.fig_benchmark(prof.ranges, A.LABELS, ["swl", "sorw", "krwend", "kroend"], cur, "End-points"), use_container_width=True)
    g2.plotly_chart(ch.fig_benchmark(prof.ranges, A.LABELS, ["nw", "no"], cur, "Corey exponents"), use_container_width=True)
    bt = A.benchmark_table()
    st.dataframe(bt[bt["profile"] == prof.title].drop(columns="profile"), hide_index=True, use_container_width=True)
    st.caption("Basis: 'Craig 1971' = published rule of thumb; 'published' = literature range; 'judgement' = typical practice, not a measured dataset.")
    for line in prof.notes:
        st.markdown(f"- {line}")

    st.markdown("**Capillary-pressure starting values**")
    q1, q2, q3 = st.columns(3)
    adv_k = q1.number_input("Permeability [mD]", 0.01, 20000.0, step=10.0, key="adv_k")
    adv_phi = q2.number_input("Porosity [-]", 0.02, 0.45, step=0.01, key="adv_phi")
    adv_j = q3.number_input("J at the entry point [-]", 0.05, 2.0, step=0.05, key="adv_j")
    ps = A.pc_start(prof, adv_k, adv_phi, adv_j)
    kpis([("Brooks-Corey lambda", f"{ps['lam_min']:.1f} / {ps['lam']:.1f} / {ps['lam_max']:.1f}"),
          ("Entry pressure", f"{ps['pd_bar'] * ch.PCF:.3g} {ch.PCU}"),
          ("Implied nw (Burdine)", f"{ps['nw_from_lam']:.1f}"),
          ("nw range from lambda", f"{ps['nw_range_from_lam'][0]:.1f} to {ps['nw_range_from_lam'][1]:.1f}")])

    if need_wo and wo_p is not None:
        st.markdown("**Assessment of your current case**")
        for sev, msg in A.assess(wo_p, mu_w, mu_o, prof, ps["lam"]):
            {"ok": st.success, "info": st.info, "warn": st.warning}[sev](msg)

    a1, a2 = st.columns(2)
    a_pc = a1.checkbox("Also set Pc to Brooks-Corey from the starting values", value=False, key="adv_with_pc")
    a1.button("Use the base case as the sidebar case", on_click=_apply_advisor, args=(pk, "base", kind_a, a_pc), key="adv_apply_base")

    st.markdown("**Low / base / high cases** (edit the table; low = pessimistic displacement, high = optimistic)")
    adv_ok = need_wo and not (pc_on and pc_wo is None)
    if not adv_ok:
        note("The low / base / high interpolation works on water-oil curves. Select Water-oil or Three-phase in the sidebar.")
    tbl0 = pd.DataFrame({k: [prof.ranges[k][0], prof.ranges[k][1], prof.ranges[k][2]] for k in A.PARAMS}, index=["low", "base", "high"])
    tbl = st.data_editor(tbl0, key=f"adv_tbl_{pk}", use_container_width=True) if adv_ok else None
    try:
        if not adv_ok:
            raise _AdvSkip()
        trip = []
        for nm in ("low", "base", "high"):
            r = tbl.loc[nm]
            krw_a = C.corey_to_let(float(r["nw"])) if kind_a == "let" else C.KrModel("corey", n=float(r["nw"]))
            kro_a = C.corey_to_let(float(r["no"])) if kind_a == "let" else C.KrModel("corey", n=float(r["no"]))
            pp = C.WOParams(float(r["swl"]), float(r["swl"]), float(r["sorw"]), float(r["krwend"]), float(r["kroend"]), krw_a, kro_a)
            pp.validate()
            trip.append(pp)
        t_a = st.slider("Interpolation (0 = low, 0.5 = base, 1 = high)", 0.0, 1.0, step=0.01, key="adv_t")
        p_t = C.interpolate_triple(*trip, t_a)
        pcm = C.interpolate_scalar_triple(*prof.pc_mult, t_a, log=True)

        def _adv_tab(p, scale=1.0):
            t = C.wateroil(p, res)
            return C.add_pc_wo(t, p, pc_wo, scale=scale) if pc_on else t.assign(PC=0.0)

        tabs_a = {"Low": _adv_tab(trip[0], prof.pc_mult[0]), "Base": _adv_tab(trip[1], prof.pc_mult[1]), "High": _adv_tab(trip[2], prof.pc_mult[2])}
        t_int = _adv_tab(p_t, pcm)
        tabs_a[f"Interpolated t={t_a:.2f}"] = t_int
        st.plotly_chart(ch.fig_interp(tabs_a), use_container_width=True)
        rows_a = []
        for nm, tb_, pp_ in (("Low", tabs_a["Low"], trip[0]), ("Base", tabs_a["Base"], trip[1]), ("High", tabs_a["High"], trip[2]), (f"Interpolated t={t_a:.2f}", t_int, p_t)):
            w_ = C.welge(tb_, mu_w, mu_o, swi=pp_.swl)
            rows_a.append({"Case": nm, "Swl": pp_.swl, "Sorw": pp_.sorw, "krw end": pp_.krwend, "kro end": pp_.kroend,
                           "Crossover Sw": w_["sw_cross"], "Mobility ratio": w_["mobility_ratio"],
                           "Breakthrough [PV]": w_["pvi_bt"], "Ultimate RF": w_["rf_ult"]})
        st.dataframe(pd.DataFrame(rows_a).round(3), hide_index=True, use_container_width=True)
        st.caption(f"Pc multiplier at this t: {pcm:.2f} (log interpolation of {prof.pc_mult[0]} / {prof.pc_mult[1]} / {prof.pc_mult[2]}). Applied only when Pc is included in the sidebar.")

        reg_i = wo_region(f"Advisor interpolated t={t_a:.2f}", t_int, p_t, pcm,
                          source=f"Advisor profile={pk}; t={t_a!r}; effective parameters from edited low/base/high table")
        regs3 = [wo_region(nm, tabs_a[nm], pp, scale, source=f"Advisor profile={pk}; edited {nm} case")
                 for nm, pp, scale in zip(("Low", "Base", "High"), trip, prof.pc_mult)]
        region_sources["Advisor: interpolated curve (1 region)"] = [reg_i]
        region_sources["Advisor: low / base / high (3 regions)"] = regs3
        d1, d2, d3 = st.columns(3)
        advisor_export = X.prepare_regions([reg_i], export_rows)
        d1.download_button("Download the interpolated curve (CSV)", show(advisor_export[0]["wo"]).to_csv(index=False).encode(), "advisor_curve.csv", "text/csv", key="adv_csv")
        advisor_keywords = "SWOF + SGOF" if "go" in reg_i else "SWOF"
        d2.download_button(f"Download the interpolated curve (Eclipse {advisor_keywords})",
                           X.export_eclipse(advisor_export, _PC_UNITS_TO_DECK[st.session_state["pc_unit"]]),
                           "advisor_swof.inc", "text/plain", key="adv_swof")
        d3.button("Use the interpolated curve as the sidebar case", on_click=_apply_wo_fit, args=(p_t, kind_a, False), key="adv_apply_t")
        st.caption("The interpolated curve and the low / base / high set are also listed as region sources on the Export tab.")
        if go_df is None:
            st.caption("To export both SWOF and SGOF, select Three-phase (oil-water-gas) in the sidebar. Advisor changes water-oil only; gas-oil curves are carried unchanged.")
    except _AdvSkip:
        pass
    except ValueError as exc:
        st.error(f"The edited table is not a valid curve set: {exc}")

# ================================================================== Uncertainty
with T["Uncertainty"]:
    guide("Uncertainty")
    action_title("P10 / P50 / P90 curve bands show how much SCAL uncertainty matters for recovery",
                 "Monte-Carlo sampling of end-points, curve shape and Pc around the base case; low / mid / high cases export as separate regions for history matching")
    if not need_wo:
        note("Ensembles are available for water-oil curves. Select Water-oil or Three-phase in the sidebar.")
    else:
        u1, u2, u3 = st.columns(3)
        e_n = u1.number_input("Realisations", 20, 1000, step=10, key="ens_n")
        e_seed = u1.number_input("Random seed", 0, 99999, step=1, key="ens_seed")
        e_dist = u1.radio("Distribution (mode = base case)", CHOICES["ens_dist"], horizontal=True, key="ens_dist")
        e_dswl = u2.slider("Swl, +/- absolute", 0.0, 0.15, step=0.005, key="ens_dswl")
        e_dsorw = u2.slider("Sorw, +/- absolute", 0.0, 0.2, step=0.005, key="ens_dsorw")
        e_fkrw = u2.slider("krw end-point, +/- relative", 0.0, 0.9, step=0.01, key="ens_fkrw")
        e_fkro = u3.slider("kro end-point, +/- relative", 0.0, 0.5, step=0.01, key="ens_fkro")
        e_dnw = u3.slider("krw shape (exponent), +/-", 0.0, 2.0, step=0.1, key="ens_dnw")
        e_dno = u3.slider("kro shape (exponent), +/-", 0.0, 2.0, step=0.1, key="ens_dno")
        e_fpc = u3.slider("Pc multiplier, sigma of ln", 0.0, 1.0, step=0.05, key="ens_fpc")
        spec = C.EnsembleSpec(n=int(e_n), seed=int(e_seed), dist=e_dist, d_swl=e_dswl, d_sorw=e_dsorw, f_krwend=e_fkrw,
                              f_kroend=e_fkro, d_nw=e_dnw, d_no=e_dno, f_pc=e_fpc)
        try:
            ens = cached_ensemble(wo_p, spec, mu_w, mu_o, max(res, 0.01), pc_wo if pc_on else None)
        except ValueError as exc:
            st.error(f"Ensemble failed: {exc}")
            st.stop()
        ens_result = ens
        S = ens["samples"]
        q = {m: np.percentile(S[m].values, [10, 50, 90]) for m in ("rf_ult", "pvi_bt", "rf_bt", "mobility_ratio")}
        action_title(f"The displacement limit spans {q['rf_ult'][0]:.0%} to {q['rf_ult'][2]:.0%} (P10 to P90) and breakthrough "
                     f"{q['pvi_bt'][0]:.2f} to {q['pvi_bt'][2]:.2f} PV across {len(S)} realisations")
        kpis([("Displacement limit P10 / P50 / P90", " / ".join(f"{v:.0%}" for v in q["rf_ult"])),
              ("Breakthrough [PV] P10 / P50 / P90", " / ".join(f"{v:.2f}" for v in q["pvi_bt"])),
              ("RF at breakthrough P10 / P50 / P90", " / ".join(f"{v:.0%}" for v in q["rf_bt"])),
              ("Mobility ratio P10 / P50 / P90", " / ".join(f"{v:.1f}" for v in q["mobility_ratio"]))])
        g1, g2 = st.columns(2)
        g1.plotly_chart(ch.fig_ens_kr(ens, wo_df), use_container_width=True)
        if pc_on and ens["pc"] is not None:
            g2.plotly_chart(ch.fig_ens_pc(ens), use_container_width=True)
        else:
            g2.plotly_chart(ch.fig_hist(S.rf_ult.values, "Displacement limit", "RF [-]"), use_container_width=True)
        h1, h2 = st.columns(2)
        h1.plotly_chart(ch.fig_hist(S.rf_ult.values, "Displacement limit (1 - Sorw, normalised)", "RF [-]",
                                    {f"P{p}": v for p, v in zip((10, 50, 90), q["rf_ult"])}), use_container_width=True)
        h2.plotly_chart(ch.fig_hist(S.pvi_bt.values, "Breakthrough", "PV injected",
                                    {f"P{p}": v for p, v in zip((10, 50, 90), q["pvi_bt"])}), use_container_width=True)
        source("Source: SCAL Studio Monte-Carlo ensemble; P10 = 10th percentile (low), P90 = 90th percentile (high); percentiles are taken point-wise")
        metric = st.selectbox("Pick low / mid / high cases by", CHOICES["ens_metric"], key="ens_metric",
                              format_func={"rf_ult": "Displacement limit", "pvi_bt": "Breakthrough PV", "rf_bt": "RF at breakthrough",
                                           "mobility_ratio": "Mobility ratio"}.get)
        def _ens_region(idx, label):
            pp = ens["params"][idx]
            tb = C.wateroil(pp, res)
            tb = C.add_pc_wo(tb, pp, pc_wo, scale=float(ens["pc_mult"][idx])) if pc_on else tb.assign(PC=0.0)
            return wo_region(f"{label} {metric}={S.loc[idx, metric]:.3f}", tb, pp, float(ens["pc_mult"][idx]),
                             source=f"Uncertainty seed={e_seed}; realisation={idx}; metric={metric}; distribution={e_dist}")

        picks = C.pick_cases(ens, metric)
        regs, rows = [], []
        for pctl, idx in picks.items():
            pp = ens["params"][idx]
            regs.append(_ens_region(idx, f"P{pctl}"))
            rows.append(dict(case=f"P{pctl}", realisation=idx, Swl=round(pp.swl, 3), Sorw=round(pp.sorw, 3), krw_end=round(pp.krwend, 3),
                             kro_end=round(pp.kroend, 3), Pc_multiplier=round(float(ens["pc_mult"][idx]), 2), RF_limit=round(S.loc[idx, "rf_ult"], 3),
                             breakthrough_PV=round(S.loc[idx, "pvi_bt"], 3)))
        region_sources["Uncertainty cases (P10 / P50 / P90 regions)"] = regs
        pc_custom = st.slider("Custom percentile: export one specific realisation", 1, 99, key="ens_pctl")
        idx_c = C.pick_cases(ens, metric, (pc_custom,))[pc_custom]
        reg_c = _ens_region(idx_c, f"P{pc_custom}")
        region_sources["Uncertainty: custom percentile (1 region)"] = [reg_c]
        st.caption(f"Realisation {idx_c} is closest to P{pc_custom} of the chosen metric: Swl {ens['params'][idx_c].swl:.3f}, Sorw {ens['params'][idx_c].sorw:.3f}, "
                   f"krw end {ens['params'][idx_c].krwend:.3f}. It is listed as a region source on the Export tab.")
        st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
        st.caption("The three representative realisations are available as SATNUM regions on the Export tab (history-matching sets). "
                   "Ensembles cover water-oil curves; gas-oil tables are carried unchanged.")
        st.download_button("Download ensemble samples (CSV)", S.to_csv(index=False).encode(), "scal_ensemble_samples.csv", "text/csv", key="ens_dl")

# ================================================================== Curve-set interpolation
with T["Curve-set interpolation"]:
    guide("Curve-set interpolation")
    action_title("Interpolating between two curve sets represents low-salinity, polymer, WAG and capillary-desaturation effects",
                 "Set A is the sidebar base case; define set B and drive the blend with a manual factor, a physical variable or the capillary number")
    if not need_wo:
        note("Curve-set interpolation is available for water-oil curves. Select Water-oil or Three-phase in the sidebar.")
    else:
        b1, b2, b3 = st.columns(3)
        i_swl = b1.slider("Swl (set B)", 0.0, 0.5, step=0.01, key="ip_swl")
        i_swcr = b1.slider("Swcr (set B)", 0.0, 0.6, step=0.01, key="ip_swcr")
        i_sorw = b2.slider("Sorw (set B)", 0.0, 0.6, step=0.01, key="ip_sorw")
        i_krwend = b2.slider("krw at Sorw (set B)", 0.02, 1.0, step=0.01, key="ip_krwend")
        i_kroend = b3.slider("kro at Swl (set B)", 0.1, 1.0, step=0.01, key="ip_kroend")
        i_dnw = b3.slider("krw shape change, B vs A", -3.0, 3.0, step=0.1, key="ip_dnw")
        i_dno = b3.slider("kro shape change, B vs A", -3.0, 3.0, step=0.1, key="ip_dno")
        mode_i = st.radio("What drives the blend", CHOICES["ip_var"], horizontal=True, key="ip_var")
        t_val = 0.5
        if mode_i == "Manual":
            t_val = st.slider("Interpolation factor t (0 = set A, 1 = set B)", 0.0, 1.0, step=0.01, key="ip_t")
        else:
            m1, m2, m3 = st.columns(3)
            if mode_i == "Capillary number":
                mu_i = m1.number_input("Displacing viscosity [cP]", 0.1, 500.0, step=0.1, key="ip_mu")
                v_i = m1.number_input("Darcy velocity [m/day]", 1e-6, 100.0, step=0.1, format="%.4f", key="ip_v")
                sg_i = m1.number_input("Interfacial tension [mN/m]", 0.001, 80.0, step=1.0, key="ip_sigma")
                value_i = C.capillary_number(mu_i, v_i, sg_i)
                m1.caption(f"Nc = {value_i:.2e}")
            else:
                value_i = m1.number_input("Current value (e.g. salinity, polymer, WAG ratio)", 1e-12, 1e12, step=None, format="%.3e", key="ip_value")
            lo_i = m2.number_input("Value at set A", 1e-12, 1e12, step=None, format="%.3e", key="ip_lo")
            hi_i = m2.number_input("Value at set B", 1e-12, 1e12, step=None, format="%.3e", key="ip_hi")
            sc_i = m3.radio("Interpolate on", CHOICES["ip_scale"], horizontal=True, key="ip_scale")
            t_val = C.interpolation_factor(value_i, lo_i, hi_i, sc_i)
            m3.caption(f"t = {t_val:.2f}")
        try:
            pA = wo_p
            pB = C.WOParams(i_swl, max(i_swcr, i_swl), i_sorw, i_krwend, i_kroend, C.shift_kr(wo_p.krw, i_dnw), C.shift_kr(wo_p.kro, i_dno))
            pB.validate()

            def _tab(p, scale=1.0):
                t = C.wateroil(p, res)
                return C.add_pc_wo(t, p, pc_wo, scale=scale) if pc_on else t.assign(PC=0.0)

            tB = _tab(pB)
            tI = _tab(C.interpolate_wo(pA, pB, t_val))
            wA, wB, wI = wo_w, C.welge(tB, mu_w, mu_o, swi=pB.swl), C.welge(tI, mu_w, mu_o, swi=C.interpolate_wo(pA, pB, t_val).swl)
            action_title(f"At t = {t_val:.2f} the displacement limit is {wI['rf_ult']:.0%} (A {wA['rf_ult']:.0%}, B {wB['rf_ult']:.0%}) "
                         f"and breakthrough {wI['pvi_bt']:.2f} PV (A {wA['pvi_bt']:.2f}, B {wB['pvi_bt']:.2f})")
            st.plotly_chart(ch.fig_interp({"A (base)": wo_df, "B": tB, f"Interpolated t={t_val:.2f}": tI}), use_container_width=True)
            n_ip = st.slider("Regions to export along A to B", 2, 9, key="ip_n")
            regs_i = []
            for kk, tk in enumerate(np.linspace(0.0, 1.0, int(n_ip)), start=1):
                params = C.interpolate_wo(pA, pB, float(tk))
                reg = wo_region(f"IP{kk} t={tk:.2f}", _tab(params), params,
                                source=f"Curve-set parameter interpolation A to B; t={float(tk)!r}")
                regs_i.append(reg)
            region_sources["Interpolated curve sets (A to B regions)"] = regs_i
            reg_now = wo_region(f"Interpolated t={t_val:.2f}", tI, C.interpolate_wo(pA, pB, t_val),
                                source=f"Curve-set parameter interpolation A to B; t={t_val!r}; driver={mode_i}")
            region_sources["Interpolated curve at the current t (1 region)"] = [reg_now]
            blended_export = X.prepare_regions([reg_now], export_rows)
            st.download_button("Download the blended curve (CSV)", show(blended_export[0]["wo"]).to_csv(index=False).encode(), "interpolated_curve.csv", "text/csv", key="ip_csv")
            blended_keywords = "SWOF + SGOF" if "go" in reg_now else "SWOF"
            st.download_button(f"Download the blended curve (Eclipse {blended_keywords})", X.export_eclipse(blended_export, _PC_UNITS_TO_DECK[st.session_state["pc_unit"]]),
                               "interpolated_swof.inc", "text/plain", key="ip_swof")
            note("The blended curves are exported as SATNUM regions (Export tab). Switch regions in your simulator with its own mechanism "
                 "(user-defined region changes or built-in saturation-function interpolation, whose keywords differ per simulator). "
                 "Both sets must use the same Corey or LET parametrisation.")
        except ValueError as exc:
            st.error(f"Invalid set B: {exc}")

# ================================================================== Lab interpretation
with T["Lab interpretation"]:
    guide("Lab interpretation")
    action_title("Raw core-flood and centrifuge data can be turned into relative permeability and capillary pressure",
                 "JBN interpretation of an unsteady-state waterflood, and Pc from centrifuge average-saturation data through an exact radial forward model")
    lt1, lt2 = st.tabs(["Core flood (JBN)", "Centrifuge Pc"])
    with lt1:
        base_tbl = wo_df if need_wo else C.wateroil(C.WOParams(), 0.01)
        swi_demo = float(wo_p.swl) if need_wo else 0.10
        demo_flood = cached_flood(base_tbl, mu_w, mu_o, swi_demo)
        st.download_button("Download demo core-flood data (CSV)", demo_flood.to_csv(index=False).encode(), "demo_coreflood.csv", "text/csv", key="flood_dl")
        up_f = st.file_uploader("Core-flood CSV (columns PVI, NP, DP; first row PVI = 0)", type=["csv"], key="jbn_up")
        flood, is_demo = demo_flood, True
        if up_f is not None:
            raw_f = pd.read_csv(up_f)
            raw_f.columns = [c.strip().upper() for c in raw_f.columns]
            if {"PVI", "NP", "DP"} <= set(raw_f.columns):
                flood, is_demo = raw_f[["PVI", "NP", "DP"]].astype(float), False
            else:
                st.error("Columns needed: PVI (water injected, PV), NP (oil produced, PV), DP (pressure drop).")
        j1, j2 = st.columns(2)
        swi_in = j1.number_input("Initial water saturation Swi", 0.0, 0.6, swi_demo, 0.01, key="jbn_swi")
        kro_swi = j2.number_input("kro at Swi (effective / absolute)", 0.1, 1.0, 1.0, 0.05, key="jbn_kroswi")
        st.caption("Constant-rate, 1-D, incompressible, negligible capillary end effect. The demo is generated from the sidebar base case "
                   "with added noise, so the interpretation can be compared with the known answer." if is_demo else
                   "Uploaded data. JBN is valid after breakthrough; early points are the most sensitive to noise.")
        try:
            jb = cached_jbn(flood, swi_in, mu_w, mu_o, kro_swi)
            st.plotly_chart(ch.fig_flood(flood), use_container_width=True)
            if is_demo and len(jb):
                err = np.median(np.abs(jb.KROW - np.interp(jb.SW, base_tbl.SW, base_tbl.KROW)) / np.maximum(np.interp(jb.SW, base_tbl.SW, base_tbl.KROW), 1e-3))
                kpis([("Points interpreted", f"{len(jb)}"), ("Saturation range", f"{jb.SW.min():.2f} - {jb.SW.max():.2f}"),
                      ("Median kro error vs known truth", f"{err:.0%}")])
            else:
                kpis([("Points interpreted", f"{len(jb)}"), ("Saturation range", f"{jb.SW.min():.2f} - {jb.SW.max():.2f}" if len(jb) else "n/a")])
            fk = st.radio("Fit smooth curves to the JBN points (end-points of the sidebar case are kept)", ["corey", "let"], horizontal=True,
                          format_func=str.upper, key="jbn_fit_kind")
            ep_f = wo_p if need_wo else C.WOParams()
            fp_f, rm_f = cached_fit_relperm(jb[["SW", "KRW", "KROW"]], fk, ep_f)
            st.plotly_chart(ch.fig_jbn(jb, ref=C.wateroil(fp_f, 0.01)), use_container_width=True)
            note(f"<b>Fit ({fk.upper()}):</b> RMSE {rm_f:.4f} in kr over the interpreted range. JBN gives kr only between the breakthrough "
                 "saturation and the outlet saturation reached in the test, so extrapolating to the end-points relies on the fitted shape.")
            st.button("Use these fitted curves as the base case", on_click=_apply_wo_fit, args=(fp_f, fk), key="jbn_apply")
            st.download_button("Download JBN points (CSV)", jb.to_csv(index=False).encode(), "jbn_relperm.csv", "text/csv", key="jbn_dl")
        except ValueError as exc:
            st.error(f"JBN interpretation failed: {exc}")
    with lt2:
        FSD = C.FLUID_SYSTEMS
        cc1, cc2, cc3 = st.columns(3)
        r2_cm = cc1.number_input("Outer radius r2 [cm]", 3.0, 40.0, 11.0, 0.5, key="cf_r2")
        len_cm = cc2.number_input("Core length [cm]", 1.0, 15.0, 5.0, 0.5, key="cf_len")
        drho_g = cc3.number_input("Density difference [g/cc]", 0.05, 1.5, 1.0, 0.05, key="cf_drho")
        demo_cf = C.demo_centrifuge()
        st.download_button("Download demo centrifuge data (CSV)", demo_cf.to_csv(index=False).encode(), "demo_centrifuge.csv", "text/csv", key="cf_dl")
        up_c = st.file_uploader("Centrifuge CSV (columns RPM, SW_AVG)", type=["csv"], key="cf_up")
        cf = demo_cf
        if up_c is not None:
            raw_c = pd.read_csv(up_c)
            raw_c.columns = [c.strip().upper() for c in raw_c.columns]
            raw_c = raw_c.rename(columns={"SW": "SW_AVG", "SAVG": "SW_AVG"})
            if {"RPM", "SW_AVG"} <= set(raw_c.columns):
                cf = raw_c[["RPM", "SW_AVG"]].astype(float).sort_values("RPM").reset_index(drop=True)
            else:
                st.error("Columns needed: RPM and SW_AVG (average water saturation).")
        k1c, k2c, k3c = st.columns(3)
        cf_model = k1c.selectbox("Pc model to fit", ["brooks-corey", "van-genuchten", "let", "skjaeveland"], key="cf_model")
        lab_sys = k2c.selectbox("Measured with", list(FSD), index=list(FSD).index("Air-brine (lab)"), key="cf_lab")
        tgt_c = k3c.selectbox("Convert to", list(FSD), index=list(FSD).index("Oil-water (reservoir)" if need_wo else "Gas-water (reservoir)"), key="cf_tgt")
        try:
            fp_c, rm_c, s0_c = cached_fit_centrifuge(cf.RPM.values, cf.SW_AVG.values, r2_cm, len_cm, drho_g, cf_model, 0.0)
            rpm_g = np.linspace(float(cf.RPM.min()), float(cf.RPM.max()), 80)
            sw_g = C.centrifuge_avg_sw(rpm_g, r2_cm, len_cm, drho_g, s0_c, 0.0, fp_c)
            pin = C.centrifuge_pc_inlet(cf.RPM.values, r2_cm, len_cm, drho_g)
            kpis([("Fitted Swirr", f"{s0_c:.3f}"), ("RMSE in Sw", f"{rm_c:.4f}"), (f"Max inlet Pc [{ch.PCU}]", f"{float(pin.max()) * ch.PCF:.1f}")])
            st.plotly_chart(ch.fig_centrifuge(cf.RPM.values, cf.SW_AVG.values, {"Fitted model (radial forward model)": (sw_g, rpm_g)}),
                            use_container_width=True)
            grid_c = np.linspace(s0_c + 1e-3, 1.0, 240)
            pc_c = C.pc_curve(grid_c, s0_c, 0.0, fp_c)
            keep_c = pc_c <= 1.3 * float(pin.max())
            st.plotly_chart(ch.fig_pc({f"Fitted Pc, {lab_sys}": (grid_c[keep_c], pc_c[keep_c])}, hline0=False), use_container_width=True)
            st.dataframe(pd.DataFrame({"RPM": cf.RPM, f"Inlet Pc [{ch.PCU}]": pin * ch.PCF, "SW_AVG measured": cf.SW_AVG,
                                       "SW_AVG fitted": C.centrifuge_avg_sw(cf.RPM.values, r2_cm, len_cm, drho_g, s0_c, 0.0, fp_c)}).round(4),
                         hide_index=True, use_container_width=True)
            fconv = float(C.pc_convert(1.0, FSD[lab_sys], FSD[tgt_c]))
            note(f"<b>Conversion lab to reservoir fluids:</b> Pc is multiplied by {fconv:.3f} when applied to the case. "
                 "Centrifuge data only constrain Pc up to the inlet pressure of the fastest speed; the curve beyond that is extrapolated by the model.")
            st.button("Use this Pc in the sidebar (converted)", on_click=_apply_pc_fit, args=(cf_model, _scale_pc(fp_c, fconv), s0_c), key="cf_apply")
        except (ValueError, RuntimeError) as exc:
            st.error(f"Centrifuge fit failed: {exc}")

# ================================================================== 5 Animations
with T["Animations"]:
    guide("Animations")
    action_title("Animating the displacement shows why mobility ratio and curve shape control sweep",
                 "Interactive Buckley-Leverett waterflood and a hysteresis loop walker - drag the sliders, press play")
    if wo_df is not None:
        wo_anim = wo_df
    else:
        wo_anim = C.wateroil(C.WOParams(), 0.01)
        note("No water-oil system selected: the waterflood animation uses a default water-oil curve set.")
    hyst_payload = dict(label=sel, s_crit=float(hs.s_crit), s_max=float(hs.s_max),
                        drain=dict(s=fam_h["drainage"].S.round(5).tolist(), kr=fam_h["drainage"].KR.round(6).tolist()),
                        bound=dict(s=fam_h["boundary"].S.round(5).tolist(), kr=fam_h["boundary"].KR.round(6).tolist()),
                        scans=[dict(s_hy=s["s_hy"], s_trap=s["s_trap"], s=s["curve"].S.round(5).tolist(), kr=s["curve"].KR.round(6).tolist())
                               for s in fam_h["scans"]])
    html = build_animation_html(dict(sw=wo_anim.SW.round(5).tolist(), krw=wo_anim.KRW.round(6).tolist(), kro=wo_anim.KROW.round(6).tolist()),
                                mu_o / mu_w, hyst_payload)
    components.html(html, height=1130, scrolling=True)
    source("Source: SCAL Studio; Welge construction, no capillary or gravity effects")

# ================================================================== 6 Fit SCAL data
with T["Fit SCAL data"]:
    guide("Fit SCAL data")
    action_title("Fitting measured SCAL data turns lab points into smooth, simulator-ready curves",
                 "Upload relative permeability data (columns Sw, krw, krow; optional Pc), detect end-points and fit Corey or LET")
    up = st.file_uploader("SCAL data (CSV)", type=["csv"])
    demo = C.demo_scal_data()
    st.download_button("Download demo data (CSV)", demo.to_csv(index=False).encode(), "demo_scal.csv", "text/csv")
    data = None
    with st.expander("Read an existing Eclipse-style deck (SWOF, SGOF, SLGOF, SWFN, SGFN, SOF3)"):
        up_deck = st.file_uploader("Deck or include file", type=["inc", "include", "txt", "data", "dat"], key="deck_up")
        deck_tabs = {}
        use_deck = False
        deck_data = None
        if up_deck is not None:
            try:
                deck_tabs = D.tables_from_eclipse(up_deck.getvalue().decode("utf-8", errors="replace"))
            except ValueError as exc:
                st.error(f"Could not read the deck: {exc}")
        if up_deck is not None and not deck_tabs:
            st.warning("No saturation-function tables found. Supported keywords: " + ", ".join(D._DECK_COLS))
        if deck_tabs:
            dk = st.selectbox("Keyword", list(deck_tabs), key="deck_kw")
            di = st.selectbox("Table (region)", list(range(1, len(deck_tabs[dk]) + 1)), key="deck_i")
            dtab = deck_tabs[dk][di - 1]
            st.dataframe(show(dtab), hide_index=True, use_container_width=True, height=220)
            if dk == "SWOF":
                st.plotly_chart(ch.fig_wo(dtab), use_container_width=True)
            elif dk in ("SGOF", "SLGOF"):
                st.plotly_chart(ch.fig_go(dtab.rename(columns={"SL": "SG"}) if dk == "SLGOF" else dtab), use_container_width=True)
            st.download_button("Download this table (CSV)", dtab.to_csv(index=False).encode(), f"{dk}_{di}.csv", "text/csv", key="deck_csv")
            if dk == "SWOF":
                use_deck = st.checkbox("Use this SWOF table as the measured data for fitting below", key="deck_use")
                if use_deck:
                    deck_data = dtab.rename(columns={c: c for c in dtab.columns})[[c for c in ("SW", "KRW", "KROW", "PC") if c in dtab.columns]].copy()
    if use_deck and deck_data is not None:
        data = deck_data.sort_values("SW").reset_index(drop=True)
    elif up is not None:
        raw = pd.read_csv(up)
        raw.columns = [c.strip().upper() for c in raw.columns]
        raw = raw.rename(columns={"KRO": "KROW", "SATURATION": "SW", "KR_W": "KRW", "KR_O": "KROW", "PCOW": "PC"})
        if {"SW", "KRW", "KROW"} <= set(raw.columns):
            data = raw.sort_values("SW").reset_index(drop=True)
        else:
            st.error("Columns needed: Sw, krw, krow (kro is accepted as an alias).")
    elif st.checkbox("Use demo data", value=True):
        data = demo
    if data is not None:
        ep = C.detect_endpoints(data)
        kpis([("Swl", f"{ep.swl:.3f}"), ("Swcr", f"{ep.swcr:.3f}"), ("Sorw", f"{ep.sorw:.3f}"),
              ("krw end-point", f"{ep.krwend:.2f}"), ("kro end-point", f"{ep.kroend:.2f}")])
        fit_kind = st.radio("Fit model", ["corey", "let"], horizontal=True, format_func=str.upper)
        try:
            fitted_p, rmse = cached_fit_relperm(data, fit_kind, ep)
            fitted = C.wateroil(fitted_p, 0.01)
            st.plotly_chart(ch.fig_fit(data, fitted), use_container_width=True)
            m = fitted_p
            desc = (f"krw n = {m.krw.n:.2f}, kro n = {m.kro.n:.2f}" if fit_kind == "corey" else
                    f"krw LET = ({m.krw.L:.2f}, {m.krw.E:.2f}, {m.krw.T:.2f}); kro LET = ({m.kro.L:.2f}, {m.kro.E:.2f}, {m.kro.T:.2f})")
            note(f"<b>Fit result</b> ({fit_kind.upper()}): {desc}. RMSE = {rmse:.4f} in kr.")
            if "PC" in data.columns:
                st.plotly_chart(ch.fig_pc({"Measured Pc": (data.SW, data.PC)}), use_container_width=True)

            st.button("Use fitted curves as the base case", on_click=_apply_wo_fit, args=(fitted_p, fit_kind))
        except Exception as e:  # noqa: BLE001
            st.error(f"Fit failed: {e}")

# ================================================================== Batch mode
with T["Batch mode"]:
    guide("Batch mode")
    action_title("Describe the cases in one YAML file, pick the engine and download the finished curves without touching the other tabs",
                 "The same file runs from the command line: python scal_batch.py config.yaml --engine native --out results")
    src_mode = st.radio("YAML input", ["Edit the template", "Upload a file"], horizontal=True, key="batch_src")
    engines = ["native"] + (["pyscal"] if PYSCAL_AVAILABLE else [])
    eng = st.selectbox("Engine", engines, key="batch_engine")
    if not PYSCAL_AVAILABLE:
        st.caption("The pyscal engine is offered when the pyscal package is installed (pip install pyscal). It supports water-oil and gas-oil only.")
    st.download_button("Download the template (YAML)", B.TEMPLATE.encode(), "scal_batch_template.yaml", "text/yaml", key="batch_tpl")
    if src_mode == "Upload a file":
        up_y = st.file_uploader("YAML file", type=["yaml", "yml", "txt"], key="batch_up")
        batch_text = up_y.getvalue().decode("utf-8", errors="replace") if up_y is not None else ""
        if up_y is not None:
            with st.expander("File content"):
                st.code(batch_text, language="yaml")
    else:
        batch_text = st.text_area("YAML", B.TEMPLATE, height=420, key="batch_text")
    if st.button("Run batch", type="primary", key="batch_run"):
        if not batch_text.strip():
            st.session_state["_batch_res"] = ("error", "No YAML provided.")
        else:
            try:
                st.session_state["_batch_res"] = ("ok", B.run_batch(batch_text, engine_override=eng))
            except B.BatchError as exc:
                st.session_state["_batch_res"] = ("error", str(exc))
    out_b = st.session_state.get("_batch_res")
    if out_b:
        if out_b[0] == "error":
            st.error(out_b[1])
        else:
            br = out_b[1]
            (st.success if br.ok else st.warning)(f"{len(br.summary)} case(s) written, {len(br.files)} file(s); verification {'passed' if br.ok else 'reported issues'}.")
            st.dataframe(br.summary, hide_index=True, use_container_width=True)
            if br.checks:
                st.markdown("**Deck verification**")
                st.dataframe(pd.DataFrame(br.checks), hide_index=True, use_container_width=True)
            for wmsg in br.warnings:
                st.warning(wmsg)
            st.download_button("Download all results (ZIP)", br.zip_bytes(), "scal_batch_results.zip", "application/zip", key="batch_zip")
            with st.expander("Files in the archive"):
                st.write(sorted(br.files))

# ================================================================== Validation
with T["Validation"]:
    guide("Validation")
    action_title("The engine is checked against published numbers, closed-form solutions and an independent numerical solver",
                 "Kinds: published value, closed form, independent numerical, published guidance, self-consistency")
    if st.button("Run the validation suite", type="primary", key="val_run"):
        st.session_state["_val_res"] = cached_validation()
    vr = st.session_state.get("_val_res")
    if vr:
        sm = V.summary(vr)
        kpis([("Passed", str(sm["PASS"])), ("Failed", str(sm["FAIL"])), ("Skipped", str(sm["SKIP"]))])
        for c in vr:
            with st.expander(f"{c.status}  {c.id}  [{c.kind}]  {c.title}"):
                st.markdown(f"**Reference.** {c.reference}")
                if c.checks:
                    st.dataframe(pd.DataFrame([dict(check=k.name, expected=k.expected, got=k.got, tolerance=k.tol, ok=k.ok) for k in c.checks]),
                                 hide_index=True, use_container_width=True)
                if c.note:
                    st.caption(c.note)
    else:
        note("Press the button to run the suite (a few seconds). Skipped cases need an optional package such as pyscal.")

# ================================================================== 7 Export
with T["Export and report"]:
    guide("Export and report")
    action_title("One case, every simulator: export consistent saturation tables for the whole toolchain",
                 "Eclipse-style decks are read by OPM Flow, tNavigator and INTERSECT; CMG and MRST have native writers")
    src_name = st.selectbox("Region source (SATNUM tables to write)", list(region_sources))
    regions = region_sources[src_name]
    e1, e2, e3 = st.columns(3)
    _pc_idx = {"bar": 0, "psi": 1, "kPa": 2, "MPa": 0}[st.session_state["pc_unit"]]
    units = e1.selectbox("Pc units in files", ["METRIC", "FIELD", "SI-kPa"], index=_pc_idx,
                         format_func={"METRIC": "bar (METRIC)", "FIELD": "psi (FIELD)", "SI-kPa": "kPa"}.get)
    fam_choice = e2.selectbox("Eclipse table family", ["SWOF/SGOF", "SWOF/SLGOF", "SWFN/SGFN/SOF3"])
    stone = e3.selectbox("Three-phase model keyword", ["STONE2", "STONE1", "none"])
    stone = None if stone == "none" else stone
    use_hyst = False
    if hyst_export is not None:
        use_hyst = st.checkbox("Include hysteresis (imbibition gas tables + EHYSTR)", value=False)
        if use_hyst:
            regions = [dict(r, imb_go=hyst_export["imb_go"], provenance={
                **r.get("provenance", {}),
                "imb_go": dict(drainage=base_region["provenance"]["go"], scheme=hyst_export["scheme"],
                               s_crit=hs.s_crit, s_max=hs.s_max, s_trap_max=hs.s_trap_max,
                               boundary_curve_points=200, pc="unchanged from base gas-oil drainage"),
            }) if "go" in r else r for r in regions]
    regions = X.prepare_regions(regions, export_rows)
    st.caption(f"{len(regions)} region(s), {export_rows} rows per table; source: {src_name}.")
    if need_wo and not need_go:
        st.info("Water-oil produces SWOF only. Select Three-phase (oil-water-gas) in the sidebar for both SWOF and SGOF.")

    files: dict[str, bytes] = {}
    eclipse_txt = X.export_eclipse(regions, units, fam_choice, stone, hyst_export if use_hyst else None)
    cmg_txt = X.export_cmg(regions, units if units != "METRIC" else "SI-kPa", stone)
    mrst_txt = X.export_mrst(regions)
    chosen = st.multiselect("Simulators", list(X.SIMULATORS), default=list(X.SIMULATORS)[:1] + ["CMG IMEX / GEM (.dat)", "MRST (MATLAB .m)"])
    for name in chosen:
        kind_x, fname = X.SIMULATORS[name]
        files[fname] = {"eclipse": eclipse_txt, "cmg": cmg_txt, "mrst": mrst_txt}[kind_x].encode()
    files["scal_tables.csv"] = X.export_csv(regions)
    files["scal_tables.xlsx"] = X.export_xlsx(regions)
    files["scal_tables.json"] = X.export_json(regions, dict(system=system, kr_model=kind, pc_units=units))

    cols = st.columns(4)
    for i, (fname, blob) in enumerate(files.items()):
        cols[i % 4].download_button(f"{fname}", blob, fname, key=f"dl_{fname}")
    zbuf = io.BytesIO()
    with zipfile.ZipFile(zbuf, "w", zipfile.ZIP_DEFLATED) as z:
        for fname, blob in files.items():
            z.writestr(fname, blob)
    st.download_button("Download everything (ZIP)", zbuf.getvalue(), "scal_export.zip", "application/zip")
    # ---- independent round-trip verification of what was written
    with st.expander("Verify the exported decks (parsed back and compared with the tables)", expanded=True):
        chk_units_cmg = units if units != "METRIC" else "SI-kPa"
        decks = {"Eclipse-style keywords": D.verify_eclipse(eclipse_txt, regions, units, fam_choice, hyst_export if use_hyst else None),
                 "CMG *SWT / *SLT": D.verify_cmg(cmg_txt, regions, chk_units_cmg),
                 "MRST arrays": D.verify_mrst(mrst_txt, regions)}
        for dname, chks in decks.items():
            n_ok, n_tot = D.summarize(chks)
            (st.success if n_ok == n_tot else st.error)(f"{dname}: {n_ok} of {n_tot} checks passed")
            if n_ok != n_tot:
                st.dataframe(pd.DataFrame([dict(check=c.name, status="PASS" if c.ok else "FAIL", detail=c.detail) for c in chks if not c.ok]),
                             hide_index=True, use_container_width=True)
        st.caption("Each deck is read back by an independent parser (keyword tables, *RPT blocks, MATLAB arrays) and compared "
                   "number by number with the in-app tables, including the Pc unit factor.")
    prev = st.selectbox("Preview", ["Eclipse / OPM / tNavigator / INTERSECT", "CMG IMEX / GEM", "MRST"])
    st.code({"Eclipse / OPM / tNavigator / INTERSECT": eclipse_txt, "CMG IMEX / GEM": cmg_txt, "MRST": mrst_txt}[prev], language="text")
    note("Check keyword placement and units against your simulator version before running; hysteresis keywords (EHYSTR items, "
         "CMG *HYSKR*) differ between versions and are provided as a starting template.")

    # ---- one-click PDF summary
    st.markdown("---")
    st.markdown("**One-click PDF summary (A4, two pages)**")
    if not need_wo:
        note("The PDF summary is built around the water-oil displacement. Select Water-oil or Three-phase in the sidebar.")
    else:
        case_name = st.text_input("Case name for the report", "SCAL Studio case", key="case_name")

        def _build_pdf():
            rows = [("Fluid system", system), ("Kr parametrisation", kind_label), ("Swl / Swcr / Sorw", f"{wo_p.swl:.3f} / {wo_p.swcr:.3f} / {wo_p.sorw:.3f}"),
                    ("krw at Sorw / kro at Swl", f"{wo_p.krwend:.2f} / {wo_p.kroend:.2f}"),
                    ("Water / oil viscosity [cP]", f"{mu_w:g} / {mu_o:g}"),
                    ("Capillary pressure", (f"{pc_wo.kind}, Pc max {pc_wo.pcmax * ch.PCF:.1f} {ch.PCU}" if pc_on and pc_wo is not None else "not included")),
                    ("Table resolution dS", f"{res}")]
            rows.append(("Krw / kro shape", f"n = {wo_p.krw.n:.2f} / {wo_p.kro.n:.2f}" if wo_p.krw.kind == "corey" else
                         f"LET krw ({wo_p.krw.L:.1f}, {wo_p.krw.E:.1f}, {wo_p.krw.T:.1f}); kro ({wo_p.kro.L:.1f}, {wo_p.kro.E:.1f}, {wo_p.kro.T:.1f})"))
            if ens_result is not None:
                rows.append(("Uncertainty ensemble", f"{len(ens_result['samples'])} realisations, P10-P90 band shown"))
            rows.append(("Regions in export", f"{len(regions)} ({src_name})"))
            notes = ["Curves are generated from parametric Corey / LET models; calibrate against measured SCAL data."]
            if ens_result is not None:
                notes.append("The shaded kr bands are the point-wise P10 to P90 of the Monte-Carlo ensemble.")
            return R.build_report(case_name, wo_df, mu_w, mu_o, rows, pc_unit=ch.PCU, ensemble=ens_result, notes=notes)

        if st.button("Prepare PDF report", key="pdf_btn"):
            st.session_state["_pdf"] = _build_pdf()
        if st.session_state.get("_pdf"):
            st.download_button("Download PDF report", st.session_state["_pdf"], "scal_summary.pdf", "application/pdf", key="pdf_dl")
            st.caption("Prepared from the case as it was when you pressed the button; prepare it again after changing inputs.")

# ================================================================== About / licence / footer
with st.sidebar.expander("About and licence"):
    st.markdown(
        "**SCAL Studio** was made by **Merouane Hamdani** for **educational purposes only**, under the **MIT License**.\n\n"
        "Results are not validated for field development, reserves booking or investment decisions. "
        "Always check tables and keywords against your simulator manual and your laboratory data.")
    _lic = Path(__file__).parent / "LICENSE"
    if _lic.exists():
        st.code(_lic.read_text(), language="text")

st.markdown(
    '<div class="src" style="margin-top:28px">SCAL Studio &middot; made by Merouane Hamdani &middot; for educational purposes only &middot; '
    'MIT License &middot; Not for use in field development, reserves or investment decisions without independent validation</div>',
    unsafe_allow_html=True)

# ================================================================== Documentation
with T["Documentation"]:
    action_title("Every formula in the engine, with its assumptions, code location and reference",
                 "Search a keyword, then open a section; the same text is available as DOCUMENTATION.md")
    q = st.text_input("Search the documentation", key="doc_query").strip().lower()
    shown = 0
    for sec in DC.SECTIONS:
        blob = (sec["title"] + sec["text"] + " ".join(c for c, _ in sec["formulas"])).lower()
        if q and q not in blob:
            continue
        shown += 1
        with st.expander(sec["title"], expanded=bool(q)):
            st.markdown(sec["text"])
            for cap, tex in sec["formulas"]:
                st.caption(cap)
                st.latex(tex)
            if sec.get("code"):
                st.caption("Implemented in: " + ", ".join(sec["code"]))
            if sec.get("refs"):
                st.caption("References: " + "; ".join(sec["refs"]))
    if not shown:
        st.info("No section matches this search.")
    with st.expander("References"):
        for k, v in DC.REFERENCES.items():
            st.markdown(f"- **{k}**: {v}")
    st.download_button("Download DOCUMENTATION.md", DC.render_markdown().encode(), "DOCUMENTATION.md", "text/markdown", key="doc_md")
    st.caption("Made by Merouane Hamdani for educational purposes only, MIT License.")
