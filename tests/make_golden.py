# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
"""Regenerate the golden tables:  python tests/make_golden.py   (only when a deliberate numerical change is intended)."""
import pathlib, sys
ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import core as C

OUT = ROOT / "tests" / "golden"


def cases():
    wo_c = C.WOParams()
    wo_l = C.WOParams(krw=C.KrModel("let", L=2.0, E=3.0, T=1.5), kro=C.KrModel("let", L=2.5, E=2.0, T=1.2))
    pc_bc, pc_ts = C.PcParams(kind="brooks-corey", pd=0.3, pcmax=5.0, lam=2.0), C.PcParams(kind="let-two-sided", pcmax=3.0, pcmin=1.5)
    g, w = C.GOParams(), C.GWParams()
    return {
        "wo_corey_pc_bc": C.add_pc_wo(C.wateroil(wo_c, 0.02), wo_c, pc_bc),
        "wo_let_pc_twosided": C.add_pc_wo(C.wateroil(wo_l, 0.02), wo_l, pc_ts),
        "go_corey": C.gasoil(g, 0.02),
        "gw_corey": C.gaswater(w, 0.02),
    }


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    for name, df in cases().items():
        df.to_csv(OUT / f"{name}.csv", index=False, float_format="%.10g")
        print("wrote", name, df.shape)
