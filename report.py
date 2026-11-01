# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
"""report.py - one-click, printable A4 PDF summary (McKinsey-style: action title, KPIs, charts, assumptions).
Uses matplotlib only (no Streamlit / Plotly dependency) so it also runs in CI."""
from __future__ import annotations

import io
from datetime import date

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

import core as C

NAVY, BLUE, CYAN, GREY, LIGHT, RED = "#051C2C", "#2251FF", "#00A9F4", "#6B7C8F", "#E8EEF5", "#C0392B"
CREDIT = "Made by Merouane Hamdani  |  For educational purposes only  |  MIT License"
DISCLAIMER = "Illustrative results from user inputs. Not validated for field development, reserves or investment decisions."
PC_UNITS = {"bar": 1.0, "psi": C.BAR_TO_PSI, "kPa": C.BAR_TO_KPA, "MPa": 0.1}


def _page(title, subtitle, page_no, n_pages):
    fig = plt.figure(figsize=(8.27, 11.69))
    fig.patch.set_facecolor("white")
    band = fig.add_axes([0, 0.925, 1, 0.075]); band.set_facecolor(NAVY); band.set_xticks([]); band.set_yticks([])
    for s in band.spines.values():
        s.set_visible(False)
    band.text(0.05, 0.62, "SCAL Studio", color="white", fontsize=9, fontweight="bold", transform=band.transAxes)
    band.text(0.05, 0.22, "Relative permeability and capillary pressure summary", color=CYAN, fontsize=7.5, transform=band.transAxes)
    band.text(0.95, 0.4, date.today().isoformat(), color="white", fontsize=8, ha="right", transform=band.transAxes)
    fig.text(0.06, 0.885, title, fontsize=12.5, fontweight="bold", color=NAVY, va="top", wrap=True)
    fig.text(0.06, 0.832, subtitle, fontsize=8.5, color=GREY, va="top")
    fig.add_artist(plt.Line2D([0.06, 0.94], [0.03, 0.03], color=BLUE, lw=0.8))
    fig.text(0.06, 0.018, CREDIT, fontsize=6.5, color=GREY)
    fig.text(0.94, 0.018, f"Page {page_no} of {n_pages}", fontsize=6.5, color=GREY, ha="right")
    fig.text(0.06, 0.006, DISCLAIMER, fontsize=5.8, color=GREY)
    return fig


def _style(ax, xlabel, ylabel, title):
    ax.set_title(title, loc="left", fontsize=8.5, fontweight="bold", color=NAVY)
    ax.set_xlabel(xlabel, fontsize=7, color=GREY); ax.set_ylabel(ylabel, fontsize=7, color=GREY)
    ax.tick_params(labelsize=7, colors=GREY)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#B8C4D2")
    ax.grid(True, color=LIGHT, lw=0.7)


def _kpi(fig, x, y, label, value, unit=""):
    ax = fig.add_axes([x, y, 0.205, 0.075]); ax.set_xticks([]); ax.set_yticks([]); ax.set_facecolor(LIGHT)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.add_patch(plt.Rectangle((0, 0), 0.025, 1, color=BLUE, transform=ax.transAxes))
    ax.text(0.09, 0.62, value, fontsize=13, fontweight="bold", color=NAVY, transform=ax.transAxes, va="center")
    ax.text(0.09, 0.2, f"{label}{(' (' + unit + ')') if unit else ''}", fontsize=6.5, color=GREY, transform=ax.transAxes, va="center")


def _table(fig, rect, rows, header, col_w=None):
    ax = fig.add_axes(rect); ax.axis("off")
    if not rows:
        return ax
    t = ax.table(cellText=rows, colLabels=header, loc="upper left", cellLoc="left", colWidths=col_w)
    t.auto_set_font_size(False); t.set_fontsize(7); t.scale(1, 1.25)
    for (r, c), cell in t.get_celld().items():
        cell.set_edgecolor("white")
        if r == 0:
            cell.set_facecolor(NAVY); cell.get_text().set_color("white"); cell.get_text().set_fontweight("bold")
        else:
            cell.set_facecolor(LIGHT if r % 2 else "white"); cell.get_text().set_color(NAVY)
    return ax


def build_report(case_name, wo, mu_w, mu_o, rows_assumptions, pc_unit="bar", ensemble=None, go=None, gw=None,
                 wo_endpoints=None, notes=None) -> bytes:
    """Return PDF bytes. `wo` = water-oil table (SW, KRW, KROW, PC). `ensemble` optional dict from core.ensemble_wo."""
    pcf = PC_UNITS.get(pc_unit, 1.0)
    w = C.welge(wo, mu_w, mu_o)
    rec = C.recovery_curve(w, 3.0)
    rf2 = float(np.interp(2.0, rec.PVI, rec.RF))
    m = w["mobility_ratio"]
    verdict = "favourable" if m < 1 else ("near-unit" if m < 2 else "unfavourable")
    n_pages = 2
    buf = io.BytesIO()
    with PdfPages(buf) as pdf:
        # -------- page 1
        title = (f"Water breakthrough after {w['pvi_bt']:.2f} PV; displacement is {verdict} "
                 f"(end-point mobility ratio {m:.1f})")
        fig = _page(title, f"Case: {case_name}   |   Welge analytics without capillary effects", 1, n_pages)
        _kpi(fig, 0.06, 0.74, "Breakthrough", f"{w['pvi_bt']:.2f}", "PV injected")
        _kpi(fig, 0.285, 0.74, "RF at breakthrough", f"{100 * w['rf_bt']:.0f} %")
        _kpi(fig, 0.51, 0.74, "RF at 2 PV", f"{100 * rf2:.0f} %")
        _kpi(fig, 0.735, 0.74, "Ultimate RF (1-Sorw)", f"{100 * w['rf_ult']:.0f} %")
        ax = fig.add_axes([0.09, 0.40, 0.38, 0.28])
        ax.plot(wo.SW, wo.KRW, color=BLUE, lw=1.8, label="krw"); ax.plot(wo.SW, wo.KROW, color=NAVY, lw=1.8, label="kro")
        if ensemble is not None:
            g = ensemble["grid"]
            for key, col in (("krw", BLUE), ("kro", NAVY)):
                lo, hi = np.percentile(ensemble[key], [10, 90], axis=0)
                ax.fill_between(g, lo, hi, color=col, alpha=0.15, lw=0)
        _style(ax, "Sw", "kr", "Relative permeability" + (" (P10-P90 band)" if ensemble is not None else ""))
        ax.legend(fontsize=7, frameon=False)
        ax = fig.add_axes([0.57, 0.40, 0.37, 0.28])
        ax.plot(w["sw"], w["fw"], color=BLUE, lw=1.8)
        ax.plot([w["swi"], w["swf"]], [0, w["fwf"]], color=CYAN, lw=1, ls="--")
        ax.plot(w["swf"], w["fwf"], "o", color=NAVY, ms=4)
        _style(ax, "Sw", "fw", "Fractional flow and Welge tangent")
        ax = fig.add_axes([0.09, 0.08, 0.38, 0.24])
        ax.plot(rec.PVI, 100 * rec.RF, color=BLUE, lw=1.8); ax.axvline(w["pvi_bt"], color=GREY, ls=":", lw=1)
        _style(ax, "Pore volumes injected", "Oil recovery (% of OIIP)", "Recovery factor vs PVI")
        ax = fig.add_axes([0.57, 0.08, 0.37, 0.24])
        if "PC" in wo and np.nanmax(np.abs(wo.PC)) > 0:
            ax.plot(wo.SW, wo.PC * pcf, color=NAVY, lw=1.8); ax.axhline(0, color=GREY, lw=0.6)
            _style(ax, "Sw", f"Pc [{pc_unit}]", "Capillary pressure")
        else:
            ax.axis("off"); ax.text(0, 0.5, "Capillary pressure not enabled", fontsize=8, color=GREY)
        pdf.savefig(fig); plt.close(fig)
        # -------- page 2
        fig = _page("Assumptions, end-points and how to read these results",
                    f"Case: {case_name}", 2, n_pages)
        header = ["Parameter", "Value"]
        rows = [[a, b] for a, b in rows_assumptions]
        _table(fig, [0.06, 0.45, 0.42, 0.33], rows[:18], header, [0.62, 0.38])
        ep = wo_endpoints or {}
        right = [["Swl", f"{wo.SW.iloc[0]:.3f}"], ["Krw max", f"{wo.KRW.max():.3f}"], ["Kro max", f"{wo.KROW.max():.3f}"],
                 ["Crossover Sw", "n/a" if np.isnan(w["sw_cross"]) else f"{w['sw_cross']:.3f}"],
                 ["Mu_w / Mu_o [cP]", f"{mu_w:g} / {mu_o:g}"], ["Table points", f"{len(wo)}"]] + [[k, str(v)] for k, v in ep.items()]
        _table(fig, [0.54, 0.45, 0.40, 0.33], right, ["Derived", "Value"], [0.6, 0.4])
        if "PC" in wo and np.nanmax(np.abs(wo.PC)) > 0:
            msgs = C.pc_checks(wo.SW.values, wo.PC.values) if hasattr(C, "pc_checks") else []
            txt = "Pc checks: " + ("all passed" if not msgs else "; ".join(msgs))
        else:
            txt = "Pc checks: capillary pressure not enabled"
        fig.text(0.06, 0.40, txt, fontsize=7.5, color=NAVY, wrap=True, va="top")
        lines = (notes or []) + [
            "Reading guide: breakthrough and recovery use Welge/Buckley-Leverett theory (1-D, incompressible, no gravity or capillary effects).",
            "Use these curves as a starting point; calibrate against core data and history-match before field use.",
        ]
        y = 0.36
        fig.text(0.06, y, "Notes", fontsize=9, fontweight="bold", color=NAVY, va="top"); y -= 0.025
        for ln in lines:
            fig.text(0.06, y, "• " + ln, fontsize=7.3, color=NAVY, va="top", wrap=True); y -= 0.03
        pdf.savefig(fig); plt.close(fig)
        info = pdf.infodict(); info["Title"] = f"SCAL Studio summary - {case_name}"; info["Author"] = "Merouane Hamdani (SCAL Studio, MIT, educational)"
    return buf.getvalue()
