# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
"""state_defs.py - default widget state and the allowed choices (pure Python, shared by the app, case files and tests)."""
from __future__ import annotations

import advisor as A
import core as C

SYSTEMS = ["Water-oil", "Gas-oil", "Gas-water", "Three-phase (oil-water-gas)"]
DEFAULTS = dict(
    system=SYSTEMS[0], kind="Corey", res=0.01, export_rows=51,
    mu_w=0.5, mu_o=5.0,
    wo_swl=0.10, wo_swcr=0.10, wo_sorw=0.20, wo_krwend=0.40, wo_kroend=0.90,
    wo_krw_n=3.0, wo_kro_n=2.5, wo_krw_L=2.0, wo_krw_E=2.0, wo_krw_T=2.0, wo_kro_L=2.0, wo_kro_E=2.0, wo_kro_T=2.0,
    go_swl=0.10, go_sgcr=0.05, go_sorg=0.15, go_krgend=0.85, go_kroend=0.90,
    go_krg_n=2.0, go_kro_n=2.5, go_krg_L=2.0, go_krg_E=2.0, go_krg_T=2.0, go_kro_L=2.0, go_kro_E=2.0, go_kro_T=2.0,
    gw_swl=0.20, gw_sgcr=0.05, gw_sgr=0.25, gw_krwend=0.60, gw_krgend=0.50,
    gw_krw_n=3.0, gw_krg_n=2.0, gw_krw_L=2.0, gw_krw_E=2.0, gw_krw_T=2.0, gw_krg_L=2.0, gw_krg_E=2.0, gw_krg_T=2.0,
    pc_on=True,
    pcwo_kind="brooks-corey", pcwo_pd=0.30, pcwo_pcmax=5.0, pcwo_pcmin=1.0, pcwo_lam=2.0, pcwo_L=1.5, pcwo_E=4.0, pcwo_T=1.5, pcwo_vg=2.5, pcwo_cw=0.3, pcwo_aw=0.5, pcwo_co=0.15, pcwo_ao=1.0, pcwo_smooth=True,
    pcgo_kind="let", pcgo_pd=0.02, pcgo_pcmax=1.5, pcgo_pcmin=0.5, pcgo_lam=2.0, pcgo_L=1.5, pcgo_E=4.0, pcgo_T=1.5, pcgo_vg=2.5, pcgo_cw=0.05, pcgo_aw=0.5, pcgo_co=0.02, pcgo_ao=1.0, pcgo_smooth=True,
    pcgw_kind="brooks-corey", pcgw_pd=0.50, pcgw_pcmax=8.0, pcgw_pcmin=1.0, pcgw_lam=2.0, pcgw_L=1.5, pcgw_E=4.0, pcgw_T=1.5, pcgw_vg=2.5, pcgw_cw=0.5, pcgw_aw=0.5, pcgw_co=0.25, pcgw_ao=1.0, pcgw_smooth=True,
    sc_swl=0.15, sc_swcr=0.15, sc_sorw=0.28, sc_krwend=0.55, sc_kroend=0.80, sc_krwr=0.55, sc_krwmax=0.85,
)

DEFAULTS.update(
    pc_unit="bar", len_unit="m",
    ens_n=200, ens_seed=7, ens_dist="triangular", ens_dswl=0.03, ens_dsorw=0.05, ens_fkrw=0.20, ens_fkro=0.10,
    ens_dnw=0.5, ens_dno=0.5, ens_fpc=0.25, ens_metric="rf_ult",
    ip_var="Manual", ip_t=0.5, ip_scale="log", ip_lo=1e-7, ip_hi=1e-3, ip_value=1e-5, ip_mu=1.0, ip_v=1.0, ip_sigma=30.0, ip_n=5,
    ip_swl=0.10, ip_swcr=0.05, ip_sorw=0.12, ip_krwend=0.55, ip_kroend=0.90, ip_dnw=-1.0, ip_dno=0.0,
    adv_profile="water_wet_sandstone", adv_t=0.5, adv_k=150.0, adv_phi=0.22, adv_j=0.3, ens_pctl=50,
)

CHOICES = dict(
    system=SYSTEMS, kind=["Corey", "LET"], res=[0.02, 0.01, 0.005],
    pcwo_kind=C.PC_MODELS, pcgo_kind=C.PC_MODELS, pcgw_kind=C.PC_MODELS,
    pc_unit=["bar", "psi", "kPa", "MPa"], len_unit=["m", "ft"],
    ens_dist=["triangular", "uniform"], ens_metric=["rf_ult", "pvi_bt", "rf_bt", "mobility_ratio"],
    adv_profile=list(A.PROFILES), ip_var=["Manual", "Variable (salinity, polymer, ...)", "Capillary number"], ip_scale=["log", "linear"],
)
