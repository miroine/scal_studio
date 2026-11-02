# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
"""docs_content.py - single source for the in-app Documentation tab, the per-tab reading guides and DOCUMENTATION.md.

Every formula below is the one implemented in core.py (function names are listed under 'Implemented in').
Regenerate the Markdown with:  python docs_content.py > DOCUMENTATION.md   (a test keeps the file in sync).
"""
from __future__ import annotations

REFERENCES = {
    "Corey 1954": "Corey, A.T. (1954) The interrelation between gas and oil relative permeabilities, Producers Monthly 19(1), 38-41.",
    "Brooks and Corey 1964": "Brooks, R.H. and Corey, A.T. (1964) Hydraulic Properties of Porous Media, Colorado State University, Hydrology Paper 3.",
    "Burdine 1953": "Burdine, N.T. (1953) Relative permeability calculations from pore size distribution data, Trans. AIME 198, 71-78.",
    "van Genuchten 1980": "van Genuchten, M.T. (1980) A closed-form equation for predicting the hydraulic conductivity of unsaturated soils, Soil Sci. Soc. Am. J. 44, 892-898.",
    "Lomeland et al. 2005": "Lomeland, F., Ebeltoft, E. and Thomas, W.H. (2005) A new versatile relative permeability correlation, SCA2005-32. https://www.jgmaas.com/SCA/2005/SCA2005-32.pdf",
    "Skjaeveland et al. 1998": "Skjaeveland, S.M., Siqveland, L.M., Kjosavik, A., Hammervold, W.L. and Virnovsky, G.A. (1998) Capillary pressure correlation for mixed-wet reservoirs, SPE 60900.",
    "Stone 1970, 1973": "Stone, H.L. (1970) Probability model for estimating three-phase relative permeability, JPT 22(2); Stone, H.L. (1973) Estimation of three-phase relative permeability and residual oil data, J. Can. Pet. Tech. 12(4).",
    "Baker 1988": "Baker, L.E. (1988) Three-phase relative permeability correlations, SPE 17369.",
    "Leverett 1941": "Leverett, M.C. (1941) Capillary behavior in porous solids, Trans. AIME 142, 152-169.",
    "Hassler and Brunner 1945": "Hassler, G.L. and Brunner, E. (1945) Measurement of capillary pressures in small core samples, Trans. AIME 160, 114-123.",
    "Land 1968": "Land, C.S. (1968) Calculation of imbibition relative permeability for two- and three-phase flow from rock properties, SPE J. 8(2), 149-156.",
    "Killough 1976": "Killough, J.E. (1976) Reservoir simulation with history-dependent saturation functions, SPE J. 16(1), 37-48.",
    "Carlson 1981": "Carlson, F.M. (1981) Simulation of relative permeability hysteresis to the nonwetting phase, SPE 10157.",
    "Buckley and Leverett 1942": "Buckley, S.E. and Leverett, M.C. (1942) Mechanism of fluid displacement in sands, Trans. AIME 146, 107-116.",
    "Welge 1952": "Welge, H.J. (1952) A simplified method for computing oil recovery by gas or water drive, Trans. AIME 195, 91-98.",
    "JBN 1959": "Johnson, E.F., Bossler, D.P. and Naumann, V.O. (1959) Calculation of relative permeability from displacement experiments, Trans. AIME 216, 370-372.",
    "Craig 1971": "Craig, F.F. (1971) The Reservoir Engineering Aspects of Waterflooding, SPE Monograph Vol. 3.",
    "Amott 1959": "Amott, E. (1959) Observations relating to the wettability of porous rock, Trans. AIME 216, 156-162.",
    "Donaldson et al. 1969": "Donaldson, E.C., Thomas, R.D. and Lorenz, P.B. (1969) Wettability determination and its effect on recovery efficiency, SPE J. 9(1), 13-20 (USBM index).",
    "pyscal documentation": "Equinor pyscal documentation. https://equinor.github.io/pyscal/",
    "PetroWiki": "PetroWiki, Relative permeability models. https://petrowiki.spe.org/Relative_permeability_models",
}

# Each section: id, title, text (markdown), formulas [(caption, latex)], code (implemented in), refs (keys of REFERENCES)
SECTIONS = [
    dict(id="conventions", title="1. Conventions and units", refs=["pyscal documentation"], code=["core.py (whole module)"],
         text="""All saturations are fractions of pore volume. Pressures inside the engine are in **bar**, lengths in **metres**, permeability in **mD**,
viscosity in **cP**, interfacial tension in **mN/m**. The display unit switch changes only what is shown, written to files or downloaded.

Phase and curve naming follows the Eclipse/pyscal convention: `Swl` connate (irreducible) water, `Swcr` critical water saturation,
`Sorw` residual oil to water, `Sorg` residual oil to gas, `Sgcr` critical gas, `Sgr` residual (trapped) gas, `krwend`, `kroend`, `krgend` end-point relative permeabilities.""",
         formulas=[]),
    dict(id="kr-normalisation", title="2. Relative permeability: normalised saturations and end-points", refs=["pyscal documentation", "Corey 1954", "Lomeland et al. 2005"],
         code=["core.wateroil", "core.gasoil", "core.gaswater", "core.KrModel.f"],
         text="""Each curve is an end-point multiplied by a shape function of a **normalised saturation** x in [0, 1] (clipped). The normalisations follow pyscal.
Above 1 - Sorw the water curve is held at its end-point (pyscal can ramp linearly to a separate `krwmax`; this engine does not, so the cross-check compares Sw in [Swl, 1 - Sorw]).""",
         formulas=[("Water-oil", r"x_w=\frac{S_w-S_{wcr}}{1-S_{wcr}-S_{orw}},\quad x_o=\frac{1-S_w-S_{orw}}{1-S_{orw}-S_{wl}},\quad k_{rw}=k_{rw}^{end}\,f_w(x_w),\quad k_{ro}=k_{ro}^{end}\,f_o(x_o)"),
                   ("Gas-oil", r"x_g=\frac{S_g-S_{gcr}}{1-S_{gcr}-S_{org}-S_{wl}},\quad x_{o}=\frac{1-S_g-S_{org}-S_{wl}}{1-S_{org}-S_{wl}}"),
                   ("Gas-water (drainage)", r"x_w=\frac{S_w-S_{wl}}{1-S_{wl}-S_{gr}},\quad x_g=\frac{S_g-S_{gcr}}{1-S_{wl}-S_{gcr}}"),
                   ("Corey shape", r"f(x)=x^{n}"),
                   ("LET shape (Lomeland et al. 2005)", r"f(x)=\frac{x^{L}}{x^{L}+E\,(1-x)^{T}},\qquad L\ge 1,\; E>0,\; T\ge 0.5\ \text{(published guidance)}")]),
    dict(id="three-phase", title="3. Three-phase oil relative permeability", refs=["Stone 1970, 1973", "Baker 1988"], code=["core.three_phase_kro", "core.sof3_table"],
         text="""Oil kr in the presence of water and gas is built from the two-phase tables (krow(Sw), krog(Sg)). Values are set to zero when So = 0 and clipped to [0, kro_end].""",
         formulas=[("Stone II (as implemented)", r"k_{ro}=k_{ro}^{end}\left[\left(\frac{k_{row}}{k_{ro}^{end}}+k_{rw}\right)\left(\frac{k_{rog}}{k_{ro}^{end}}+k_{rg}\right)-(k_{rw}+k_{rg})\right]"),
                   ("Stone I", r"k_{ro}=\frac{S_o^{*}\,k_{row}\,k_{rog}}{k_{ro}^{end}(1-S_w^{*})(1-S_g^{*})},\quad S_o^{*}=\frac{S_o-S_{om}}{1-S_{wl}-S_{om}},\; S_w^{*}=\frac{S_w-S_{wl}}{1-S_{wl}-S_{om}},\; S_g^{*}=\frac{S_g}{1-S_{wl}-S_{om}},\; S_{om}=(1-\alpha)S_{orw}+\alpha S_{org},\; \alpha=\frac{S_g}{1-S_{wl}-S_{org}}"),
                   ("Baker (saturation-weighted)", r"k_{ro}=\frac{(S_w-S_{wl})\,k_{row}+S_g\,k_{rog}}{(S_w-S_{wl})+S_g}")]),
    dict(id="pc-models", title="4. Capillary pressure models", refs=["Brooks and Corey 1964", "van Genuchten 1980", "Skjaeveland et al. 1998", "Lomeland et al. 2005"],
         code=["core.pc_curve", "core.PcParams", "core.add_pc_wo", "core.add_pc_go", "core.add_pc_gw"],
         text="""All models return Pc in bar as a function of the wetting-phase saturation. With S = (Sw - Swirr) / (1 - Swirr) (the two-sided LET model uses 1 - Sorw - Swirr in the denominator).
**Smooth cap:** Brooks-Corey and van Genuchten reach Pc max exactly at Swirr by remapping S to s_min + (1 - s_min) S, instead of a flat plateau with a kink.
**Gas-oil Pc** is evaluated on the liquid saturation 1 - Sg and multiplied by min(Sg / Sgcr, 1), so it is zero at Sg = 0 and has no jump at the first gas saturation.""",
         formulas=[("Brooks-Corey", r"P_c=p_d\,S^{-1/\lambda},\qquad s_{min}=\left(\frac{p_d}{p_{c,max}}\right)^{\lambda}\ \text{(smooth cap)}"),
                   ("LET", r"P_c=p_d+(p_{c,max}-p_d)\,\frac{(1-S)^{L}}{(1-S)^{L}+E\,S^{T}}"),
                   ("van Genuchten", r"P_c=p_d\left[S^{-1/m}-1\right]^{1/n},\quad m=1-\tfrac{1}{n},\qquad s_{min}=\left(1+(p_{c,max}/p_d)^{n}\right)^{-m}"),
                   ("Skjaeveland (two-sided, mixed-wet)", r"P_c=\frac{c_w}{S_w^{*\,a_w}}-\frac{c_o}{S_o^{*\,a_o}},\quad S_w^{*}=\frac{S_w-S_{wirr}}{1-S_{wirr}},\; S_o^{*}=\frac{1-S_w-S_{or}}{1-S_{or}},\quad -p_{c,min}\le P_c\le p_{c,max}"),
                   ("LET two-sided", r"P_c=p_{c,max}\frac{(1-s)^{L}}{(1-s)^{L}+E s^{T}}-p_{c,min}\frac{s^{L_n}}{s^{L_n}+E_n(1-s)^{T_n}},\quad s=\frac{S_w-S_{wirr}}{1-S_{or}-S_{wirr}}")]),
    dict(id="pc-scaling", title="5. Fluid-system conversion, J-function and height above FWL", refs=["Leverett 1941"],
         code=["core.pc_convert", "core.j_function", "core.pc_from_j", "core.height_above_fwl", "core.saturation_height", "core.seal_column_height", "core.micp_to_pc"],
         text="""Pc is translated between fluid systems with the ratio of |sigma cos(theta)| (magnitudes, because cos 140 deg is negative for mercury-air injection) and, optionally,
between rocks with the Leverett factor. The J-function collapses the curves of one rock type measured with different fluids. Saturation-height inverts the Pc curve at Pc = drho g h.
Seal capacity is the static column a given entry pressure can hold. MICP data (Hg saturation S_Hg) are converted with a closure/conformance correction S_c.""",
         formulas=[("Conversion", r"P_{c,2}=P_{c,1}\,\frac{|\sigma_2\cos\theta_2|}{|\sigma_1\cos\theta_1|}\sqrt{\frac{\phi_2/k_2}{\phi_1/k_1}}"),
                   ("J-function (SI)", r"J=\frac{P_c}{\sigma\cos\theta}\sqrt{\frac{k}{\phi}}"),
                   ("Height above free-water level / seal column", r"h=\frac{P_c}{\Delta\rho\, g},\qquad h_{seal}=\frac{p_{entry}}{\Delta\rho\, g}"),
                   ("MICP to Pc", r"S_w=1-\frac{S_{Hg}-S_c}{1-S_c}")]),
    dict(id="pc-fit", title="6. Capillary pressure fitting and wettability indices", refs=["Amott 1959", "Donaldson et al. 1969"], code=["core.fit_pc", "core.wettability_indices"],
         text="""The fit minimises soft-L1 residuals scaled by max(|Pc|, 5 % of max |Pc|), so low-pressure points are not drowned by the steep end. **Swirr is fitted** (0 <= Swirr < min Sw): fixing it at the lowest measured Sw biases
the pore-size index (on the demo data R2 fell from 0.98 to a negative value). Brooks-Corey, van Genuchten and Skjaeveland are fitted as pure functions with the cap far beyond the data.
The indices are indicative: USBM from the areas under the positive and negative parts of a two-sided Pc curve, and an Amott-type index from the zero crossing.""",
         formulas=[("Residual", r"r_i=\frac{P_c(S_{w,i};\theta)-P_{c,i}}{\max(|P_{c,i}|,\,0.05\max|P_c|)}"),
                   ("USBM index", r"I_{USBM}=\log_{10}\frac{A_1}{A_2}"),
                   ("Amott-type water index", r"I_{A}=\frac{S_{w0}-S_{w,min}}{S_{w,max}-S_{w,min}},\quad P_c(S_{w0})=0")]),
    dict(id="scaling", title="7. End-point scaling and rock-type families", refs=["pyscal documentation"], code=["core.rescale_wo_table", "core.rock_type_family"],
         text="""**Two-point scaling** maps the new saturation range back onto the old table and scales kr by the end-point ratio. **Three-point scaling** additionally pins krw at Sorw and lets krw rise linearly to krwmax at Sw = 1.
**Rock-type families** interpolate end-points and shape parameters in ln(k) between a poor and a good rock and scale Pc with the Leverett factor relative to the geometric-mean rock.""",
         formulas=[("Scaled water kr (two-point)", r"k_{rw}^{new}(S_w)=k_{rw}^{old}\!\left(S_{wcr}^{o}+x_w^{new}(1-S_{wcr}^{o}-S_{orw}^{o})\right)\frac{k_{rw}^{end,new}}{k_{rw}^{end,old}}"),
                   ("Rock-type interpolation", r"t_i=\frac{\ln k_i-\ln k_{low}}{\ln k_{high}-\ln k_{low}},\qquad p_i=p_{low}+t_i\,(p_{high}-p_{low})"),
                   ("Pc scale of a rock type", r"m_i=\sqrt{\frac{\phi_i/k_i}{\phi_{ref}/k_{ref}}},\quad k_{ref}=\sqrt{k_{low}k_{high}},\ \phi_{ref}=\tfrac12(\phi_{low}+\phi_{high})")]),
    dict(id="hysteresis", title="8. Hysteresis", refs=["Land 1968", "Killough 1976", "Carlson 1981"], code=["core.land_trapped", "core.scanning_curve", "core.hysteresis_family", "core.pc_scanning", "core.imbibition_gas_table"],
         text="""The non-wetting phase is trapped on saturation reversal. **Land** fixes the trapped saturation; **Killough** maps the imbibition scanning curve onto the drainage curve; **Carlson** translates the boundary
imbibition curve parallel to itself through the reversal point. Capillary-pressure scanning curves blend the drainage and imbibition boundary curves with a Killough-type curvature function.""",
         formulas=[("Land trapping", r"S_{nt}=S_{nc}+\frac{d}{1+C\,d},\quad d=S_{n,hy}-S_{nc},\quad C=\frac{1}{S_{nt}^{max}-S_{nc}}-\frac{1}{S_{n}^{max}-S_{nc}}"),
                   ("Killough scanning kr", r"k_{r}^{imb}(S)=k_r^{dr}(S^{*})\,\frac{k_r^{dr}(S_{n,hy})}{k_r^{dr}(S_n^{max})},\quad S^{*}=S_{nc}+(S-S_{nt})\frac{S_n^{max}-S_{nc}}{S_{n,hy}-S_{nt}}"),
                   ("Pc scanning weight", r"f=\frac{1/(x+\varepsilon)-1/\varepsilon}{1/(x_{end}+\varepsilon)-1/\varepsilon},\qquad P_c=P_c^{dr}+(P_c^{imb}-P_c^{dr})\,f\ \ (\text{imbibition})")]),
    dict(id="displacement", title="9. Displacement analytics", refs=["Buckley and Leverett 1942", "Welge 1952", "Craig 1971"],
         code=["core.fractional_flow", "core.welge", "core.oil_produced_pv", "core.recovery_curve", "core.welge_profile"],
         text="""One-dimensional, incompressible, no gravity or capillary effects. The shock saturation comes from the tangent to the fractional-flow curve drawn from the initial water saturation (Welge construction);
breakthrough is the inverse of its slope. After breakthrough the outlet saturation S_w2 satisfies f'(S_w2) = 1/t, and the produced oil follows from the average saturation.
The end-point mobility ratio compares the water and oil end-point mobilities. **Craig rules of thumb** (Swc above 20-25 % and crossover above 50 % and krw end-point below 30 % indicate water-wet; Swc below 10-15 %, crossover below 50 % and krw end-point above 50 % indicate oil-wet) are used for the indicative wettability reading.""",
         formulas=[("Fractional flow", r"f_w=\frac{k_{rw}/\mu_w}{k_{rw}/\mu_w+k_{ro}/\mu_o}"),
                   ("Welge tangent and breakthrough", r"\left.\frac{df_w}{dS_w}\right|_{S_{wf}}=\frac{f_w(S_{wf})-f_w(S_{wi})}{S_{wf}-S_{wi}}=\frac{1}{PV_{bt}}"),
                   ("Oil produced after t pore volumes", r"N_p(t)=S_{w2}+(1-f_{w2})\,t-S_{wi},\qquad f_w'(S_{w2})=\frac1t"),
                   ("Recovery factors", r"RF_{bt}=\frac{PV_{bt}}{1-S_{wi}},\qquad RF_{ult}=\frac{S_{w,max}-S_{wi}}{1-S_{wi}}"),
                   ("End-point mobility ratio", r"M=\frac{k_{rw}^{end}/\mu_w}{k_{ro}^{end}/\mu_o}")]),
    dict(id="uncertainty", title="10. Uncertainty ensembles", refs=["Lomeland et al. 2005"], code=["core.EnsembleSpec", "core.ensemble_wo", "core.percentile_bands", "core.pick_cases"],
         text="""Monte-Carlo sampling around the base case (triangular with the base value as mode, or uniform): Swl and Swcr shift together, Sorw is clipped, end-points are scaled by (1 + relative change),
Corey exponents shift by an absolute change (LET: L and T scale by clip(1 + d/3, 0.5, 1.5)) and Pc is multiplied by a log-normal factor. P10 / P50 / P90 are **point-wise percentiles** (P10 = low, P90 = high).
Representative cases are the realisations closest to the chosen percentile of a displacement metric, and can be exported as regions.""",
         formulas=[("Sampled end-points", r"S_{wl}^{(i)}=S_{wl}+\delta_i,\quad k_{rw}^{end,(i)}=k_{rw}^{end}(1+\varepsilon_i),\quad n^{(i)}=n+\eta_i"),
                   ("Pc multiplier", r"m^{(i)}=\exp(\sigma Z_i),\quad Z_i\sim N(0,1)"),
                   ("Representative case", r"i_p=\arg\min_i\,|y_i-\text{percentile}_p(y)|")]),
    dict(id="interpolation", title="11. Interpolation between curve sets", refs=["Craig 1971"], code=["core.interpolate_wo", "core.interpolate_triple", "core.interpolation_factor", "core.capillary_number", "core.shift_kr"],
         text="""Two curve sets A and B are blended parameter by parameter (end-points and shape parameters, with Swcr kept at or above Swl), which keeps every intermediate curve a valid Corey/LET set.
The factor t comes from a manual slider, a physical variable (linear or log scale between its value at A and at B) or the capillary number. **Low / base / high** sets are blended piecewise-linearly
(t = 0.5 is the base case); a Pc multiplier can be blended in ln-space. Both sets must use the same parametrisation.""",
         formulas=[("Blend", r"p(t)=(1-t)\,p_A+t\,p_B"),
                   ("Factor from a variable", r"t=\mathrm{clip}\!\left(\frac{g(v)-g(v_A)}{g(v_B)-g(v_A)},0,1\right),\quad g=\log_{10}\ \text{or identity}"),
                   ("Capillary number", r"N_c=\frac{\mu\,v}{\sigma}"),
                   ("Low / base / high", r"p(t)=\begin{cases}p_{low}+2t\,(p_{base}-p_{low}) & t\le 0.5\\ p_{base}+(2t-1)(p_{high}-p_{base}) & t>0.5\end{cases}")]),
    dict(id="lab", title="12. Lab interpretation: JBN and centrifuge", refs=["JBN 1959", "Hassler and Brunner 1945"], code=["core.jbn", "core.synthetic_core_flood", "core.centrifuge_pc_inlet", "core.centrifuge_avg_sw", "core.fit_centrifuge"],
         text="""**JBN** interprets a constant-rate unsteady-state waterflood after breakthrough from the oil produced N_p(W_i) (spline derivative) and the relative injectivity I_r = dP_0 / dP, both versus the injected pore volumes W_i.
The outlet fractional flow of oil is f_o2 = dN_p/dW_i. **Centrifuge Pc** is obtained by fitting a parametric Pc model through the exact radial forward model (Pc zero at the outlet face, maximum at the inlet), which is preferable to the closed-form
Hassler-Brunner approximation for thick cores. Both are valid only within the saturation range the experiment actually covered.""",
         formulas=[("JBN outlet saturation", r"S_{w2}=S_{wi}+N_p-W_i f_{o2},\qquad f_{o2}=\frac{dN_p}{dW_i}"),
                   ("JBN kro and krw", r"k_{ro}=\frac{f_{o2}}{d\!\left(1/(W_iI_r)\right)/d(1/W_i)}\,k_{ro}(S_{wi}),\qquad k_{rw}=k_{ro}\frac{f_{w2}}{f_{o2}}\frac{\mu_w}{\mu_o}"),
                   ("Centrifuge capillary pressure", r"P_c(r)=\tfrac12\,\Delta\rho\,\omega^{2}\,(r_2^{2}-r^{2}),\quad \omega=\frac{2\pi\,\mathrm{rpm}}{60}"),
                   ("Average saturation", r"\bar S_w(\omega)=\frac{1}{r_2-r_1}\int_{r_1}^{r_2}S_w\!\left(P_c(r)\right)dr")]),
    dict(id="benchmarks", title="13. Benchmarks and the parameter advisor", refs=["Craig 1971", "Brooks and Corey 1964", "Burdine 1953", "Lomeland et al. 2005", "PetroWiki"], code=["advisor.PROFILES", "advisor.assess", "advisor.pc_start", "core.brooks_corey_krw_exponent", "core.corey_to_let"],
         text="""Benchmark profiles are **starting points, not data**. Each parameter carries a basis tag: *Craig 1971* (a published rule), *published* (a published range) or *judgement* (consolidated engineering judgement, not a citation).
Low / base / high refer to displacement efficiency (the low case is pessimistic). The advisor compares a case with the profile ranges, reads it with the Craig rules, reports the end-point mobility ratio and checks LET guidance (L >= 1, E > 0, T >= 0.5).
Published anchors used: Brooks-Corey lambda of 1.8-3.7 for sandpacks, 4.17 for a homogeneous consolidated sandstone and 7.3 for uniform glass beads; wide pore-size distributions give lambda below 2; modified Brooks-Corey exponents n between 1 and 6 (PetroWiki).
The **Burdine/Brooks-Corey link** relates the Pc pore-size index to the wetting-phase kr exponent, so a Pc fit can inform nw. The entry-pressure starting value uses a Leverett J at entry of about 0.3 (order of magnitude; calibrate with MICP).""",
         formulas=[("Burdine / Brooks-Corey wetting-phase kr", r"k_{rw}=S_e^{(2+3\lambda)/\lambda},\quad k_{r,nw}=(1-S_e)^2\left(1-S_e^{(2+\lambda)/\lambda}\right)\quad(\lambda=2\Rightarrow n_w=4)"),
                   ("Entry pressure from J", r"p_d=J_{entry}\,\sigma\cos\theta\,\sqrt{\phi/k}")]),
    dict(id="export", title="14. Simulator export and verification", refs=["pyscal documentation"], code=["exporters.py", "deckcheck.py"],
         text="""Eclipse-style keywords (SWOF, SGOF, SLGOF, SWFN, SGFN, SOF3, STONE1/2, EHYSTR template), CMG IMEX/GEM (*SWT, *SLT, *RPT, *KROIL), MRST arrays, CSV, XLSX and JSON.
Pc in the files uses the chosen unit (bar, psi, kPa). Every deck is parsed back by an independent reader and compared number by number with the tables (Verify panel and batch mode).
Select **Three-phase (oil-water-gas)** to write both SWOF and SGOF, including in the advisor and blended-curve direct INC downloads. Gas-oil curves are carried unchanged; Water-oil mode writes SWOF only.
**Export rows per table** sets an exact count (8-1000, default 51) for Export and advisor/blended-curve downloads, including SOF3 and imbibition tables. Linear interpolation retains the saturation range and critical/residual kr plateau boundaries within that count. Plots and calculations keep their source dS resolution; more export rows do not add model accuracy. The count is saved in case JSON; batch mode keeps its YAML resolution.
INC comments record the engine, source, full kr/Pc configuration, Pc multiplier, source dS, output settings, sampling method and row counts for each region. Derived cases use their effective parameters, not sidebar base values. Only the selected model's parameters are active (Corey n; LET L/E/T). Parameter pressures stay in bar; table pressures use the output unit. Eclipse/CMG include comments; JSON includes provenance and sampling metadata. Missing provenance is explicitly identified, not inferred.
Hysteresis keywords differ between simulator versions: the export is a starting template.""",
         formulas=[("Pc unit factors", r"1\ \text{bar}=14.5038\ \text{psi}=100\ \text{kPa}")]),
    dict(id="batch", title="15. Batch mode (YAML)", refs=[], code=["batch.py", "scal_batch.py"],
         text="""`python scal_batch.py input.yaml --out results/` validates the YAML strictly, builds every case with the chosen engine, writes the selected formats, verifies the decks and writes a summary and a manifest.
Top-level keys: `project, engine (native|pyscal), resolution, units, mu_w, mu_o, export, cases`. Case keys: `name, system, kr_model, water_oil, gas_oil, gas_water, pc, pc_multiplier, advisor {profile, t}, interpolate {from, to, t}`.
`--template` prints a complete example, `--check` validates only, `--list-profiles` lists the benchmark profiles. The pyscal engine supports water-oil and gas-oil kr tables (Pc always from the native engine) and needs `pip install pyscal`.
Exit codes: 0 success, 1 verification failed, 2 invalid input.""",
         formulas=[]),
    dict(id="validation", title="16. Validation against published and independent results", refs=["pyscal documentation", "Welge 1952", "Brooks and Corey 1964", "Land 1968", "Craig 1971", "Lomeland et al. 2005"], code=["validation.py"],
         text="""Each case states what it is compared with: a **published value** (e.g. the crossover printed in the pyscal documentation example), a **closed form**, an **independent numerical** solution (a separately coded upwind finite-volume
Buckley-Leverett solver), **published guidance** (Craig rules, LET ranges), a **synthetic recovery** test (known inputs, noisy synthetic data, interpretation must recover the inputs) or an internal **self-consistency** check.
A closed-form identity is never presented as a measurement. The comparison with the installed pyscal package runs only when pyscal is available and is skipped otherwise.""",
         formulas=[("Closed-form Buckley-Leverett check (quadratic kr, M = 1, Swc = Sor = 0)", r"S_{wf}=\tfrac{1}{\sqrt2},\quad f_{wf}=\frac{1}{2-\sqrt2},\quad PV_{bt}=2(\sqrt2-1)")]),
    dict(id="limits", title="17. Limitations", refs=[], code=[],
         text="""Not validated for field development, reserves or investment decisions. Ensembles cover water-oil curves. Displacement analytics ignore gravity, capillary end effects and compressibility.
JBN is sensitive to noise in early post-breakthrough points; centrifuge data constrain Pc only up to the highest inlet pressure. Benchmark ranges are generic and must be replaced by site SCAL data.
Capillary-pressure hysteresis and simulator hysteresis keywords are templates. The Streamlit interface has been exercised with a stand-in for Streamlit and Plotly in the development environment; run the included AppTest suite in your own environment.""",
         formulas=[]),
]

GUIDES = {
    "Summary": dict(what="Headline results of the current case: displacement limit, breakthrough, mobility ratio and a Craig-type wettability reading.",
                    read="Read the action title first, then the KPI row. The curves show krw and kro; the recovery curve flattens at the displacement limit (1 - Sorw normalised).",
                    watch="Mobility ratio above 1 means the front is unstable in the analytic model; the Craig reading is indicative only and cannot identify intermediate wettability."),
    "Relative permeability": dict(what="Interactive kr curves for water-oil, gas-oil, gas-water and the three-phase oil kr on a ternary diagram.",
                                  read="Use the log axis to see the low-kr tail that controls late-time production. On the ternary plot colours show kro for each (Sw, Sg, So) combination.",
                                  watch="Stone I, Stone II and Baker can differ strongly in the middle of the triangle; pick the one used by your simulator."),
    "Capillary pressure": dict(what="Pc models, fluid-system conversion, J-function, saturation-height, seal capacity and fitting of MICP or Pc-versus-Sw data.",
                               read="The top chart shows the base curve and its translations; J-function curves should collapse onto one line. Quality checks flag non-monotone or implausible curves.",
                               watch="Fix the fluid table before converting; mercury-air to reservoir conversion assumes the same pore structure. Fit Swirr together with the model."),
    "End-point scaling": dict(what="Re-scale the base curve to new end-points, or generate a permeability-ranked rock-type family.",
                              read="Compare the dotted base curves with the scaled ones. The rock-type table lists each type's end-points and Pc scale; the family can be exported as SATNUM regions.",
                              watch="Interpolating in ln(k) assumes a smooth trend between the poor and good rock; check intermediate types against core data."),
    "Parameter advisor": dict(what="Benchmark ranges by rock and wettability, an assessment of your current case, and a low-base-high interpolation that you can export as a single curve.",
                              read="The range chart shows low-to-high ranges with the base value and your case. Findings are ranked: warnings mean a value is outside the benchmark or conflicts with published guidance.",
                              watch="Ranges are generic starting points (basis tags tell you which are published rules and which are judgement). Replace them with your SCAL data."),
    "Uncertainty": dict(what="Monte-Carlo curve ensembles with P10/P50/P90 bands for kr and Pc and the spread in recovery and breakthrough.",
                        read="Bands are point-wise percentiles (P10 low, P90 high). The histograms show the metric distributions; pick low/mid/high cases by any metric or by a custom percentile.",
                        watch="Point-wise percentile curves are not themselves a realisation; export representative realisations for history matching."),
    "Curve-set interpolation": dict(what="Blend two curve sets with a manual factor, a physical variable or the capillary number, and export the blended curve or a set of regions.",
                                    read="Dotted curves are sets A and B, the solid curve is the blend at the current factor. The title reports displacement limit and breakthrough for A, B and the blend.",
                                    watch="Both sets must use the same Corey or LET parametrisation. Check how your simulator switches regions or interpolates saturation functions."),
    "Hysteresis": dict(what="Land trapping with Killough or Carlson scanning curves for the non-wetting phase, and Pc scanning curves.",
                       read="The family chart shows drainage, boundary imbibition and scanning curves; the table lists trapped saturation for each reversal point.",
                       watch="Simulator hysteresis keywords differ by version; the export is a template to adapt."),
    "Lab interpretation": dict(what="JBN interpretation of an unsteady-state core flood and fitting of centrifuge data to a Pc curve.",
                               read="JBN points should follow the reference curve once the demo is used; apply the fitted curves to the sidebar case. For centrifuge data compare measured and fitted average saturations.",
                               watch="JBN is valid after breakthrough and noise-sensitive early on. Centrifuge Pc is only constrained up to the highest inlet pressure."),
    "Fit SCAL data": dict(what="Upload measured kr (and optional Pc) data, detect end-points and fit Corey or LET; read existing Eclipse decks.",
                          read="End-points are detected first and held fixed during the fit; the RMSE in kr is reported. A deck reader shows tables from an uploaded Eclipse-style file.",
                          watch="End-point detection assumes the data span the full saturation range."),
    "Animations": dict(what="JavaScript animation of a Buckley-Leverett waterflood and a hysteresis loop walker.",
                       read="Move the viscosity ratio and press play; the front speed follows the Welge solution. The walker steps through drainage and imbibition reversals.",
                       watch="Animation uses the analytic solution without capillarity or gravity."),
    "Batch mode": dict(what="Run many cases from a YAML file without the interface; choose the engine and download all exports.",
                       read="Start from the template, edit, press Run. The summary table lists key metrics per case; verification shows how many keyword checks passed per file.",
                       watch="Validation errors name the key path. The pyscal engine requires the pyscal package and covers water-oil and gas-oil kr tables."),
    "Validation": dict(what="Comparison of the engine with published numbers, closed-form results and independent calculations.",
                       read="Each row names the kind of reference. Expand a case to see expected values, computed values and tolerances.",
                       watch="Closed-form checks verify the implementation of a model, not the model against nature. Skipped cases need packages that are not installed."),
    "Export and report": dict(what="Choose which curves to write (base, rock types, uncertainty cases, interpolated sets), export to simulators and build a PDF summary.",
                              read="The Verify panel parses each written deck back and compares it with the tables. The preview shows the first part of each file.",
                              watch="Check keyword placement and units against your simulator version before running."),
}


def render_markdown() -> str:
    """DOCUMENTATION.md text (LaTeX in $$ blocks)."""
    out = ["# SCAL Studio documentation", "",
           "*Made by Merouane Hamdani for educational purposes only, MIT License. Generated from `docs_content.py`; do not edit by hand.*", "",
           "## Contents", ""]
    out += [f"- [{s['title']}](#{s['id']})" for s in SECTIONS] + ["- [References](#references)", ""]
    for s in SECTIONS:
        out += [f'<a id="{s["id"]}"></a>', f"## {s['title']}", "", s["text"].strip(), ""]
        for cap, tex in s["formulas"]:
            out += [f"**{cap}**", "", f"$$ {tex} $$", ""]
        if s["code"]:
            out += ["*Implemented in:* " + ", ".join(f"`{c}`" for c in s["code"]), ""]
        if s["refs"]:
            out += ["*References:* " + "; ".join(s["refs"]), ""]
    out += ['<a id="references"></a>', "## References", ""]
    out += [f"- **{k}**: {v}" for k, v in REFERENCES.items()]
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    print(render_markdown(), end="")
