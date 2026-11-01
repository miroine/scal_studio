# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
"""Headless smoke test: runs app.py against stub streamlit/plotly modules for several configurations.
(Real Streamlit/Plotly are not needed for this; it validates the data/logic path of every tab.)"""
import sys, types, runpy, itertools, pathlib, json
ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

class AttrDict(dict):
    __getattr__ = dict.get
    def __setattr__(self, k, v): self[k] = v

STATE = AttrDict()
CAPTURED = {}
WIDGET_BOUNDS = {}
FORCE_BUTTONS = set()
SELECT_INDEX = {}
ALLOW_ERROR = False
ERRORS = []

class Ctx:
    def __enter__(self): return self
    def __exit__(self, *a): return False

class FakeSt(Ctx):
    session_state = STATE
    def __getattr__(self, name):          # any unknown call is a no-op (markdown, caption, error, ...)
        return lambda *a, **k: None
    @property
    def sidebar(self): return FakeSt()
    def stop(self): raise SystemExit("st.stop")
    def error(self, m, *a, **k):
        if ALLOW_ERROR: ERRORS.append(str(m)); return
        raise AssertionError("st.error: " + str(m))
    def columns(self, n, **k):
        n = n if isinstance(n, int) else len(n); return [FakeSt() for _ in range(n)]
    def tabs(self, labels): return [FakeSt() for _ in labels]
    def expander(self, *a, **k): return FakeSt()
    def _s(self, key, default):
        return STATE[key] if key is not None and key in STATE else default
    def _choice(self, label, options, key, default):
        options = list(options)
        v = self._s(key, default)
        if key is not None and key in STATE and v not in options:
            raise AssertionError(f"state {key}={v!r} not among options of '{label}'")
        return v
    def radio(self, label, options, key=None, **k): return self._choice(label, options, key, list(options)[0])
    def select_slider(self, label, options, key=None, **k): return self._choice(label, options, key, list(options)[0])
    def selectbox(self, label, options, index=0, key=None, **k):
        options = list(options)
        for pre, i in SELECT_INDEX.items():
            if label.startswith(pre): index = min(i, len(options) - 1)
        return self._choice(label, options, key, options[index])
    def cache_data(self, *a, **k):
        if a and callable(a[0]): return a[0]
        return lambda f: f
    def text_input(self, label, value="", key=None, **k): return self._s(key, value)
    def multiselect(self, label, options, default=None, **k): return default or []
    def checkbox(self, label, value=False, key=None, **k): return self._s(key, value)
    toggle = checkbox
    def _bounded(self, label, key, lo, hi, v):
        if key is not None:
            WIDGET_BOUNDS[key] = (lo, hi)
            if key in STATE and isinstance(STATE[key], (int, float)) and not (lo - 1e-12 <= STATE[key] <= hi + 1e-12):
                raise AssertionError(f"state {key}={STATE[key]} outside widget bounds [{lo}, {hi}] of '{label}'")
        return self._s(key, v)
    def slider(self, label, *a, key=None, **k):
        return self._bounded(label, key, a[0], a[1], k.get("value", a[2] if len(a) > 2 else a[0]))
    def number_input(self, label, *a, key=None, **k):
        v = k.get("value", a[2] if len(a) > 2 and isinstance(a[2], (int, float)) else a[0])
        return self._bounded(label, key, a[0], a[1], v)
    def data_editor(self, df, **k): return df
    def text_area(self, label, value="", key=None, **k): return self._s(key, value)
    def file_uploader(self, *a, **k): return None
    def button(self, *a, key=None, on_click=None, **k): return key in FORCE_BUTTONS
    def download_button(self, label, data, *a, **k):
        CAPTURED[label] = data; return False
    def code(self, body, **k): CAPTURED["code"] = body
    def set_page_config(self, **k): pass
    class column_config:
        TextColumn = NumberColumn = staticmethod(lambda *a, **k: None)

st_mod = types.ModuleType("streamlit"); fake = FakeSt()
for n in dir(FakeSt):
    if not n.startswith("__"): pass
st_mod.__getattr__ = lambda name: getattr(fake, name)
st_mod.session_state = STATE
st_mod.sidebar = fake.sidebar
st_mod.column_config = FakeSt.column_config
comp = types.ModuleType("streamlit.components"); comp1 = types.ModuleType("streamlit.components.v1")
comp1.html = lambda html, **k: CAPTURED.__setitem__("anim_html", html)
sys.modules.update({"streamlit": st_mod, "streamlit.components": comp, "streamlit.components.v1": comp1})
st_mod.components = comp; comp.v1 = comp1

class Any:
    def __init__(self, *a, **k): pass
    def __getattr__(self, n): return lambda *a, **k: self
go = types.ModuleType("plotly.graph_objects"); go.Figure = Any; go.Scatter = Any; go.Scatterternary = Any; go.Histogram = Any
pl = types.ModuleType("plotly"); pl.graph_objects = go
sys.modules.update({"plotly": pl, "plotly.graph_objects": go})

import importlib
import cases as K
def run(**state):
    STATE.clear(); CAPTURED.clear()
    STATE.update(state)
    ns = runpy.run_path(str(ROOT / "app.py"), run_name="__main__")
    return ns

SYS = ["Water-oil", "Gas-oil", "Gas-water", "Three-phase (oil-water-gas)"]
n = 0
for system, kind, pc_on in itertools.product(SYS, ["Corey", "LET"], [True, False]):
    run(system=system, kind=kind, pc_on=pc_on, pcwo_kind="let-two-sided" if (pc_on and kind == "LET") else "brooks-corey")
    n += 1
    assert "scal_tables.csv" in CAPTURED or any("csv" in k for k in CAPTURED), CAPTURED.keys()
    assert "anim_html" in CAPTURED and "const D =" in CAPTURED["anim_html"]
    print(f"ok  {system:32s} {kind:5s} pc={pc_on}")
for pk in ["let", "van-genuchten", "skjaeveland", "let-two-sided"]:
    run(system="Three-phase (oil-water-gas)", pc_on=True, pcwo_kind=pk, pcgo_kind=pk)
    print("ok  pc model", pk)
# resolution + the three Eclipse table families
for res in (0.02, 0.005):
    run(system="Three-phase (oil-water-gas)", res=res); print("ok  res", res)
run(system="Water-oil", pcfit_kind_sel="Pc vs Sw"); print("ok  Pc-vs-Sw fit path")
run(system="Gas-water", pcgw_kind="skjaeveland"); print("ok  gas-water Skjaeveland")
# new features: uncertainty, interpolation (all three drivers), units, region sources, PDF, lab tabs
FORCE_BUTTONS.add("pdf_btn")
for kw in (dict(system="Water-oil"), dict(system="Three-phase (oil-water-gas)", kind="LET", pc_on=True),
           dict(system="Water-oil", ens_dist="uniform", ens_metric="pvi_bt", ens_n=60),
           dict(system="Water-oil", ip_var="Variable (salinity, polymer, ...)"), dict(system="Water-oil", ip_var="Capillary number", ip_scale="linear"),
           dict(system="Water-oil", pc_unit="psi", len_unit="ft"), dict(system="Water-oil", pc_unit="MPa", pc_on=False),
           dict(system="Gas-water", pc_unit="kPa"), dict(system="Gas-oil", len_unit="ft")):
    run(**kw)
    if kw.get("system") in ("Water-oil", "Three-phase (oil-water-gas)"):
        assert "Download PDF report" in CAPTURED and CAPTURED["Download PDF report"][:5] == b"%PDF-", "PDF not produced"
        assert "Download ensemble samples (CSV)" in CAPTURED
    print("ok  new-feature scenario", kw)
FORCE_BUTTONS.clear()
# v4: advisor (every profile, both shapes, end and middle t), custom percentile, batch mode, validation, documentation
import advisor as ADV, batch as BAT
for pk in ADV.PROFILES:
    for kind_r, t in (("corey", 0.0), ("let", 0.25), ("corey", 0.5), ("let", 1.0)):
        ns = run(system="Water-oil", adv_profile=pk, adv_kind_radio=kind_r, adv_t=t, pc_on=True)
        rs = ns["region_sources"]
        assert "Advisor: interpolated curve (1 region)" in rs and len(rs["Advisor: low / base / high (3 regions)"]) == 3
        assert "Download the interpolated curve (Eclipse SWOF)" in CAPTURED and "SWOF" in CAPTURED["Download the interpolated curve (Eclipse SWOF)"]
    print("ok  advisor profile", pk)
ns = run(system="Gas-oil"); assert "Advisor: interpolated curve (1 region)" not in ns["region_sources"]; print("ok  advisor skipped for gas-oil")
for pc_ in (1, 50, 99):
    ns = run(system="Water-oil", ens_pctl=pc_)
    assert "Uncertainty: custom percentile (1 region)" in ns["region_sources"]
print("ok  custom percentile")
FORCE_BUTTONS.update({"batch_run", "val_run", "adv_apply_base", "adv_apply_t"})
for eng_ in ("native",):
    ns = run(system="Water-oil", batch_engine=eng_, batch_text=BAT.TEMPLATE)
    res_ = STATE["_batch_res"]; assert res_[0] == "ok" and res_[1].ok, res_
    assert "Download all results (ZIP)" in CAPTURED
vals = STATE["_val_res"]; assert vals and all(c.status in ("PASS", "SKIP") for c in vals), [c.id for c in vals if c.status == "FAIL"]
ALLOW_ERROR = True; ERRORS.clear()
ns = run(system="Water-oil", batch_text="cases: [unclosed"); assert STATE["_batch_res"][0] == "error" and any("YAML syntax" in e for e in ERRORS), ERRORS
ALLOW_ERROR = False
assert "Download DOCUMENTATION.md" in CAPTURED
FORCE_BUTTONS.clear()
print("ok  batch mode, validation, documentation")
# advisor callbacks through the real function, then the app must run on the resulting state (bounds enforced)
for pk in ADV.PROFILES:
    for kind_r in ("corey", "let"):
        ns = run(system="Water-oil")
        ns["_apply_advisor"](pk, "high", kind_r, True)
        run(**dict(STATE))
print("ok  advisor apply callbacks for all profiles")
# presets and case files through the real callbacks; the app must then run on the resulting state (widget bounds enforced by FakeSt)
import presets as PRS
class Upload:
    def __init__(self, text): self.text = text
    def getvalue(self): return self.text.encode()
for pname in PRS.PRESETS:
    ns = run()
    STATE["preset_sel"] = pname; ns["_apply_preset"]()
    assert STATE["_case_msg"][0] == "ok", STATE["_case_msg"]
    st_after = {k: v for k, v in STATE.items() if not k.startswith("_")}
    ns = run(**st_after); print("ok  preset", pname)
    saved = K.state_to_json(STATE, ns["SCHEMA"])
    STATE.clear(); ns = run()
    STATE["case_up"] = Upload(saved); ns["_load_case"]()
    assert STATE["_case_msg"][0] == "ok" and STATE["system"] == st_after["system"], STATE["_case_msg"]
ns = run(); STATE["case_up"] = Upload('{"app":"SCAL Studio case","state":{"wo_swl":9,"kind":"Nope","zzz":1}}'); ns["_load_case"]()
assert STATE["_case_msg"][0] == "warn" and STATE["wo_swl"] == 0.10, STATE["_case_msg"]
ns = run(); STATE["case_up"] = Upload("garbage"); ns["_load_case"](); assert STATE["_case_msg"][0] == "err"
print("ok  case-file callbacks (valid, invalid values, garbage)")
# every region source x every Eclipse family must round-trip through the deck verifier (st.error would raise)
for system in ("Water-oil", "Three-phase (oil-water-gas)"):
    for src_i in range(4):
        for fam_i in range(3):
            SELECT_INDEX.clear(); SELECT_INDEX.update({"Region source": src_i, "Eclipse table family": fam_i})
            run(system=system, pc_on=True)
print("ok  all region sources x Eclipse families verified")
SELECT_INDEX.clear()
# widget bounds must contain the case-file rules (otherwise a loaded value could crash a widget)
import cases as K
bad = []
for key, (lo, hi) in WIDGET_BOUNDS.items():
    rb = K.bounds_for(key)
    if rb and not (isinstance(lo, (int, float)) and isinstance(hi, (int, float))): continue
    if rb and (rb[0] < lo - 1e-12 or rb[1] > hi + 1e-12): bad.append((key, rb, (lo, hi)))
assert not bad, bad
print("rule/widget bounds consistent for", sum(1 for k in WIDGET_BOUNDS if K.bounds_for(k)), "keys")
print("code preview head:", CAPTURED["code"][:200].replace("\n", " | "))
print(f"{n} configuration runs passed")
