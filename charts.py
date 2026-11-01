# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
"""charts.py - Plotly figure builders with a McKinsey-style template (white, navy, blue accents)."""
from __future__ import annotations

import numpy as np
import plotly.graph_objects as go

NAVY, BLUE, CYAN, GREY, LIGHT, DEEP = "#051C2C", "#2251FF", "#00A9F4", "#6B7C93", "#D5DEEA", "#1F40E6"
PALETTE = [BLUE, NAVY, CYAN, "#7C8DA6", "#0B5CAD", "#8FB4FF", "#3B4A5C", "#00D1FF"]


# display units (set once per run by the app; all computation stays in bar / metres)
PCU, PCF = "bar", 1.0
LENU, LENF = "m", 1.0
_PC_FACTORS = {"bar": 1.0, "psi": 14.5037738, "kPa": 100.0, "MPa": 0.1}


def set_units(pc_unit: str = "bar", len_unit: str = "m"):
    global PCU, PCF, LENU, LENF
    PCU, PCF = pc_unit, _PC_FACTORS.get(pc_unit, 1.0)
    LENU, LENF = len_unit, (3.280839895 if len_unit == "ft" else 1.0)


def style(fig, xt="", yt="", height=430, legend=True):
    fig.update_layout(
        template="plotly_white", height=height, margin=dict(l=60, r=24, t=24, b=56),
        font=dict(family="Arial, Helvetica, sans-serif", size=13, color=NAVY),
        plot_bgcolor="white", paper_bgcolor="white", hovermode="x unified",
        legend=dict(orientation="h", y=-0.2, x=0, font=dict(size=12)) if legend else None, showlegend=legend,
        colorway=PALETTE,
    )
    fig.update_xaxes(title=xt, showline=True, linecolor=NAVY, linewidth=1.4, gridcolor="#EEF2F8", zeroline=False, ticks="outside")
    fig.update_yaxes(title=yt, showline=True, linecolor=NAVY, linewidth=1.4, gridcolor="#EEF2F8", zeroline=False, ticks="outside")
    return fig


def _l(name, x, y, color, dash=None, width=3):
    return go.Scatter(x=x, y=y, name=name, mode="lines", line=dict(color=color, width=width, dash=dash))


def fig_wo(df, ref=None, log=False, title_ref="Reference"):
    f = go.Figure()
    if ref is not None:
        f.add_trace(_l(f"krw ({title_ref})", ref.SW, ref.KRW, "#8FB4FF", "dot", 2))
        f.add_trace(_l(f"krow ({title_ref})", ref.SW, ref.KROW, "#B8C7DB", "dot", 2))
    f.add_trace(_l("krw  water", df.SW, df.KRW, BLUE))
    f.add_trace(_l("krow  oil", df.SW, df.KROW, NAVY))
    style(f, "Water saturation, Sw [-]", "Relative permeability [-]")
    if log:
        f.update_yaxes(type="log", range=[-4, 0.1])
    return f


def fig_go(df, log=False):
    f = go.Figure()
    f.add_trace(_l("krg  gas", df.SG, df.KRG, CYAN))
    f.add_trace(_l("krog  oil", df.SG, df.KROG, NAVY))
    style(f, "Gas saturation, Sg [-]", "Relative permeability [-]")
    if log:
        f.update_yaxes(type="log", range=[-4, 0.1])
    return f


def fig_gw(df, log=False):
    f = go.Figure()
    f.add_trace(_l("krw  water", df.SW, df.KRW, BLUE))
    f.add_trace(_l("krg  gas", df.SW, df.KRG, CYAN))
    style(f, "Water saturation, Sw [-]", "Relative permeability [-]")
    if log:
        f.update_yaxes(type="log", range=[-4, 0.1])
    return f


def fig_ternary(sw, sg, kro):
    f = go.Figure(go.Scatterternary(
        a=sg, b=1 - sw - sg, c=sw, mode="markers",
        marker=dict(size=9, color=kro, colorscale=[[0, "#EAF0FF"], [0.5, CYAN], [1, NAVY]], showscale=True,
                    colorbar=dict(title="kro", thickness=12), line=dict(width=0)),
        hovertemplate="Sg %{a:.2f}<br>So %{b:.2f}<br>Sw %{c:.2f}<br>kro %{marker.color:.3f}<extra></extra>"))
    f.update_layout(
        ternary=dict(sum=1, aaxis=dict(title="Gas", min=0, linecolor=NAVY), baxis=dict(title="Oil", min=0, linecolor=NAVY),
                     caxis=dict(title="Water", min=0, linecolor=NAVY), bgcolor="white"),
        height=470, margin=dict(l=30, r=30, t=20, b=20), paper_bgcolor="white",
        font=dict(family="Arial", color=NAVY))
    return f


def fig_fw(w):
    f = go.Figure()
    f.add_trace(_l("Fractional flow fw", w["sw"], w["fw"], NAVY))
    f.add_trace(_l("Welge tangent", [w["swi"], w["swi"] + 1 / w["slope"]], [w["fw"][0], w["fw"][0] + 1], CYAN, "dash", 2))
    f.add_trace(go.Scatter(x=[w["swf"]], y=[w["fwf"]], mode="markers", name="Shock front",
                           marker=dict(size=11, color=BLUE, line=dict(color="white", width=2))))
    style(f, "Water saturation, Sw [-]", "Fractional flow [-]")
    f.update_yaxes(range=[0, 1.02])
    return f


def fig_recovery(rec, w):
    f = go.Figure()
    f.add_trace(_l("Oil recovery factor", rec.PVI, rec.RF, BLUE))
    f.add_hline(y=w["rf_ult"], line=dict(color=GREY, dash="dot"), annotation_text=f"Displacement limit {w['rf_ult']:.0%}",
                annotation_position="bottom right")
    f.add_vline(x=w["pvi_bt"], line=dict(color=CYAN, dash="dash"), annotation_text=f"Breakthrough {w['pvi_bt']:.2f} PV",
                annotation_position="top right")
    style(f, "Pore volumes injected [PV]", "Recovery factor [-]", legend=False)
    f.update_yaxes(range=[0, min(1, w["rf_ult"] * 1.15)], tickformat=".0%")
    return f


def fig_pc(curves: dict, xt="Water saturation, Sw [-]", yt=None, secondary=None, hline0=True):
    """curves: name -> (x, y[, dash]); with the default yt the y values are bar and are shown in the display unit."""
    f = go.Figure()
    k = PCF if yt is None else 1.0
    yt = f"Capillary pressure [{PCU}]" if yt is None else yt
    for i, (name, c) in enumerate(curves.items()):
        f.add_trace(_l(name, c[0], np.asarray(c[1], float) * k, PALETTE[i % len(PALETTE)], c[2] if len(c) > 2 else None, 3))
    if hline0:
        f.add_hline(y=0, line=dict(color=GREY, width=1))
    style(f, xt, yt)
    return f


def fig_pc_height(sw, h, name="Pc as height above FWL"):
    f = go.Figure(_l(name, sw, np.asarray(h, float) * LENF, BLUE))
    style(f, "Water saturation, Sw [-]", f"Height above free-water level [{LENU}]", legend=False)
    return f


def fig_j(curves: dict):
    f = go.Figure()
    for i, (n, (x, y)) in enumerate(curves.items()):
        f.add_trace(_l(n, x, y, PALETTE[i % len(PALETTE)]))
    style(f, "Water saturation, Sw [-]", "Leverett J-function [-]")
    f.update_yaxes(type="log")
    return f


def fig_family(fam, key="KROW"):
    f = go.Figure()
    for i, r in enumerate(fam):
        t = r["table"]
        shade = i / max(len(fam) - 1, 1)
        col = f"rgb({int(5 + 29 * shade)},{int(28 + 53 * shade)},{int(44 + 211 * shade)})"
        f.add_trace(go.Scatter(x=t.SW, y=t.KRW, mode="lines", line=dict(color=col, width=2.5),
                               name=f"RT{r['id']}  k={r['k']:.0f} mD", legendgroup=str(r["id"])))
        f.add_trace(go.Scatter(x=t.SW, y=t.KROW, mode="lines", line=dict(color=col, width=2.5, dash="dash"),
                               showlegend=False, legendgroup=str(r["id"])))
    style(f, "Water saturation, Sw [-]", "Relative permeability [-]   (solid krw, dashed krow)")
    return f


def fig_hyst(fa, y="Relative permeability [-]", xt="Non-wetting saturation, Snw [-]", sel=None):
    f = go.Figure()
    f.add_trace(_l("Primary drainage", fa["drainage"].S, fa["drainage"].KR, BLUE))
    f.add_trace(_l("Boundary imbibition", fa["boundary"].S, fa["boundary"].KR, GREY, "dash", 2.5))
    for i, s in enumerate(fa["scans"]):
        f.add_trace(go.Scatter(x=s["curve"].S, y=s["curve"].KR, mode="lines", name=f"Scanning from {s['s_hy']:.2f}",
                               line=dict(color=CYAN, width=1.8), opacity=0.9, showlegend=(i == 0), legendgroup="scan"))
    style(f, xt, y)
    return f


def fig_fit(data, fitted, cols=("KRW", "KROW")):
    f = go.Figure()
    f.add_trace(go.Scatter(x=data.SW, y=data.KRW, mode="markers", name="krw data", marker=dict(color=BLUE, size=8)))
    f.add_trace(go.Scatter(x=data.SW, y=data.KROW, mode="markers", name="krow data", marker=dict(color=NAVY, size=8, symbol="square")))
    f.add_trace(_l("krw fit", fitted.SW, fitted.KRW, BLUE))
    f.add_trace(_l("krow fit", fitted.SW, fitted.KROW, NAVY))
    style(f, "Water saturation, Sw [-]", "Relative permeability [-]")
    return f


def fig_shf(sw, h):
    f = go.Figure(_l("Saturation-height function", sw, np.asarray(h, float) * LENF, BLUE))
    style(f, "Water saturation, Sw [-]", f"Height above free-water level [{LENU}]", legend=False)
    f.update_xaxes(range=[0, 1.0])
    return f


def fig_pc_fit(sw, pc, curves: dict, log=False):
    f = go.Figure()
    f.add_trace(go.Scatter(x=sw, y=np.asarray(pc, float) * PCF, mode="markers", name="Measured", marker=dict(color=NAVY, size=8, symbol="square")))
    for i, (name, (x, y)) in enumerate(curves.items()):
        f.add_trace(_l(name, x, np.asarray(y, float) * PCF, PALETTE[(i * 2) % len(PALETTE)] if i else BLUE, None, 2.8))
    style(f, "Water saturation, Sw [-]", f"Capillary pressure [{PCU}]")
    if log:
        f.update_yaxes(type="log")
    return f


# ------------------------------------------------------------------ v3: uncertainty, interpolation, lab interpretation
def _band(f, x, mat, color, name, k=1.0, dash=None):
    lo, mid, hi = (np.percentile(mat, p, axis=0) * k for p in (10, 50, 90))
    rgb = tuple(int(color[i:i + 2], 16) for i in (1, 3, 5))
    f.add_trace(go.Scatter(x=x, y=hi, mode="lines", line=dict(width=0), showlegend=False, hoverinfo="skip"))
    f.add_trace(go.Scatter(x=x, y=lo, mode="lines", line=dict(width=0), fill="tonexty", fillcolor=f"rgba({rgb[0]},{rgb[1]},{rgb[2]},0.18)",
                           name=f"{name} P10-P90", hoverinfo="skip"))
    f.add_trace(_l(f"{name} P50", x, mid, color, dash, 3))


def fig_ens_kr(ens, base=None):
    f = go.Figure()
    _band(f, ens["grid"], ens["krw"], BLUE, "krw")
    _band(f, ens["grid"], ens["kro"], NAVY, "kro")
    if base is not None:
        f.add_trace(_l("Base case krw", base.SW, base.KRW, BLUE, "dot", 1.5))
        f.add_trace(_l("Base case kro", base.SW, base.KROW, NAVY, "dot", 1.5))
    style(f, "Water saturation, Sw [-]", "Relative permeability [-]")
    return f


def fig_ens_pc(ens):
    f = go.Figure()
    _band(f, ens["grid"], np.clip(ens["pc"], 0, None), BLUE, "Pc", PCF)
    style(f, "Water saturation, Sw [-]", f"Capillary pressure [{PCU}]")
    return f


def fig_hist(values, title, xt, marks=None):
    f = go.Figure(go.Histogram(x=values, marker=dict(color=BLUE), nbinsx=30, opacity=0.9))
    for name, v in (marks or {}).items():
        f.add_vline(x=v, line=dict(color=NAVY, width=1.5, dash="dash"), annotation_text=name, annotation_position="top")
    style(f, xt, "Realisations", legend=False, height=330)
    f.update_layout(title=dict(text=title, font=dict(size=13)))
    return f


def fig_interp(sets: dict, t_curve=None):
    """sets: name -> WO table. Dotted ends A/B, solid interpolated."""
    f = go.Figure()
    for i, (name, d) in enumerate(sets.items()):
        dash = None if name.startswith("Interpolated") else "dot"
        col = BLUE if name.startswith("Interpolated") else PALETTE[(i * 2 + 2) % len(PALETTE)]
        f.add_trace(_l(f"krw {name}", d.SW, d.KRW, col, dash, 3 if dash is None else 2))
        f.add_trace(_l(f"kro {name}", d.SW, d.KROW, NAVY if dash is None else "#7C8DA6", dash, 3 if dash is None else 2))
    style(f, "Water saturation, Sw [-]", "Relative permeability [-]")
    return f


def fig_jbn(jb, ref=None):
    f = go.Figure()
    f.add_trace(go.Scatter(x=jb.SW, y=jb.KRW, mode="markers", name="krw (JBN)", marker=dict(color=BLUE, size=8)))
    f.add_trace(go.Scatter(x=jb.SW, y=jb.KROW, mode="markers", name="kro (JBN)", marker=dict(color=NAVY, size=8, symbol="square")))
    if ref is not None:
        f.add_trace(_l("krw reference", ref.SW, ref.KRW, BLUE, "dot", 1.5)); f.add_trace(_l("kro reference", ref.SW, ref.KROW, NAVY, "dot", 1.5))
    style(f, "Outlet water saturation, Sw [-]", "Relative permeability [-]")
    return f


def fig_flood(data):
    f = go.Figure()
    f.add_trace(go.Scatter(x=data.PVI, y=data.NP, mode="markers+lines", name="Oil produced [PV]", line=dict(color=BLUE), marker=dict(size=6)))
    f.add_trace(go.Scatter(x=data.PVI, y=data.DP / data.DP.iloc[0], mode="markers+lines", name="Pressure drop / initial", yaxis="y2",
                           line=dict(color=NAVY, dash="dot"), marker=dict(size=6)))
    style(f, "Water injected [PV]", "Oil produced [PV]")
    f.update_layout(yaxis2=dict(title="dP / dP0", overlaying="y", side="right", showgrid=False))
    return f


def fig_centrifuge(rpm, sw, fit_curves: dict):
    f = go.Figure()
    f.add_trace(go.Scatter(x=sw, y=np.asarray(rpm, float), mode="markers", name="Measured average Sw", marker=dict(color=NAVY, size=9, symbol="square")))
    for i, (name, (x, y)) in enumerate(fit_curves.items()):
        f.add_trace(_l(name, x, y, BLUE, None, 2.8))
    style(f, "Average water saturation, Sw [-]", "Centrifuge speed [rpm]")
    return f


def fig_benchmark(ranges: dict, labels: dict, params: list, current: dict, title: str = ""):
    """Horizontal range chart: benchmark low-to-high range, base value and the user's current case."""
    f = go.Figure()
    names = [labels[k] for k in params]
    for i, k in enumerate(params):
        lo, ba, hi = ranges[k]
        a, b = min(lo, hi), max(lo, hi)
        f.add_trace(go.Scatter(x=[a, b], y=[labels[k]] * 2, mode="lines", line=dict(color="#B8C7DB", width=14), name="Benchmark range", showlegend=(i == 0)))
        f.add_trace(go.Scatter(x=[ba], y=[labels[k]], mode="markers", marker=dict(color=NAVY, size=11, symbol="diamond"), name="Base", showlegend=(i == 0)))
        v = current.get(k)
        if v is not None and v == v:
            f.add_trace(go.Scatter(x=[v], y=[labels[k]], mode="markers", marker=dict(color=BLUE, size=12, symbol="circle"), name="Your case", showlegend=(i == 0)))
    style(f, "Value [-]", "", height=300)
    f.update_yaxes(categoryorder="array", categoryarray=names[::-1])
    f.update_layout(title=dict(text=title, font=dict(size=13)))
    return f
