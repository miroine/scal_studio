# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
"""presets.py - illustrative starting cases. Generic textbook-style values, NOT from a specific field:
replace them with your own SCAL data before drawing conclusions."""
from __future__ import annotations

_WO = "Water-oil"
_GO = "Gas-oil"
_GW = "Gas-water"
_3P = "Three-phase (oil-water-gas)"

PRESETS = {
    "Water-wet sandstone, waterflood": dict(
        description="Strongly water-wet clastic: high Swcr, low krw end-point, crossover above Sw = 0.5.",
        state=dict(system=_WO, kind="Corey", mu_w=0.5, mu_o=3.0, wo_swl=0.22, wo_swcr=0.22, wo_sorw=0.28, wo_krwend=0.25, wo_kroend=0.85,
                   wo_krw_n=4.0, wo_kro_n=2.5, pc_on=True, pcwo_kind="brooks-corey", pcwo_pd=0.35, pcwo_pcmax=6.0, pcwo_lam=2.0)),
    "Mixed-wet sandstone, waterflood": dict(
        description="Mixed-wet: low Swcr, higher krw end-point, lower Sorw, Pc crossing zero (two-sided).",
        state=dict(system=_WO, kind="Corey", mu_w=0.5, mu_o=4.0, wo_swl=0.12, wo_swcr=0.08, wo_sorw=0.18, wo_krwend=0.55, wo_kroend=0.85,
                   wo_krw_n=2.5, wo_kro_n=3.0, pc_on=True, pcwo_kind="let-two-sided", pcwo_pcmax=3.0, pcwo_pcmin=1.5,
                   pcwo_L=1.5, pcwo_E=3.0, pcwo_T=1.5)),
    "Oil-wet carbonate, waterflood": dict(
        description="Oil-wet carbonate: low Swcr, high krw end-point, crossover below Sw = 0.5, mostly negative Pc.",
        state=dict(system=_WO, kind="Corey", mu_w=0.6, mu_o=2.5, wo_swl=0.08, wo_swcr=0.05, wo_sorw=0.20, wo_krwend=0.75, wo_kroend=0.80,
                   wo_krw_n=1.8, wo_kro_n=3.5, pc_on=True, pcwo_kind="skjaeveland", pcwo_pcmax=4.0, pcwo_pcmin=5.0,
                   pcwo_cw=0.08, pcwo_aw=0.6, pcwo_co=0.5, pcwo_ao=0.9)),
    "Heavy-oil waterflood (unfavourable mobility)": dict(
        description="Viscous oil: end-point mobility ratio well above 1, early breakthrough, long tail.",
        state=dict(system=_WO, kind="Corey", mu_w=0.7, mu_o=80.0, wo_swl=0.20, wo_swcr=0.20, wo_sorw=0.25, wo_krwend=0.30, wo_kroend=0.80,
                   wo_krw_n=3.0, wo_kro_n=2.0, pc_on=True, pcwo_kind="brooks-corey", pcwo_pd=0.15, pcwo_pcmax=3.0, pcwo_lam=1.5)),
    "Gas cap / gas injection (gas-oil)": dict(
        description="Gas-oil displacement: small critical gas saturation, moderate residual oil to gas.",
        state=dict(system=_GO, kind="Corey", go_swl=0.12, go_sgcr=0.05, go_sorg=0.15, go_krgend=0.80, go_kroend=0.90,
                   go_krg_n=2.0, go_kro_n=2.5, pc_on=True, pcgo_kind="let", pcgo_pd=0.02, pcgo_pcmax=1.5)),
    "Tight gas sandstone (gas-water)": dict(
        description="High irreducible water, low krw end-point, strong residual gas trapping on imbibition.",
        state=dict(system=_GW, kind="Corey", gw_swl=0.45, gw_sgcr=0.05, gw_sgr=0.20, gw_krwend=0.30, gw_krgend=0.60,
                   gw_krw_n=4.0, gw_krg_n=2.0, pc_on=True, pcgw_kind="brooks-corey", pcgw_pd=3.0, pcgw_pcmax=40.0, pcgw_lam=1.2)),
    "CO2 storage, saline aquifer (gas-water)": dict(
        description="Drainage curves for CO2 injection into brine; use the Hysteresis tab for residual (trapped) CO2.",
        state=dict(system=_GW, kind="Corey", gw_swl=0.25, gw_sgcr=0.05, gw_sgr=0.30, gw_krwend=0.55, gw_krgend=0.45,
                   gw_krw_n=3.5, gw_krg_n=2.0, pc_on=True, pcgw_kind="brooks-corey", pcgw_pd=0.12, pcgw_pcmax=6.0, pcgw_lam=2.5)),
    "WAG sandstone (three-phase)": dict(
        description="Oil-water plus gas-oil sets for water-alternating-gas studies; check hysteresis for the gas phase.",
        state=dict(system=_3P, kind="Corey", mu_w=0.5, mu_o=2.0, wo_swl=0.15, wo_swcr=0.15, wo_sorw=0.22, wo_krwend=0.35, wo_kroend=0.90,
                   wo_krw_n=3.0, wo_kro_n=2.5, go_sgcr=0.05, go_sorg=0.15, go_krgend=0.85, go_kro_n=2.5, go_krg_n=2.0,
                   pc_on=True, pcwo_kind="brooks-corey", pcwo_pd=0.30, pcwo_pcmax=5.0, pcwo_lam=2.0)),
}
