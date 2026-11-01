# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
import pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
import core as C, report as R


def test_report_pdf():
    p, pc = C.WOParams(), C.PcParams()
    wo = C.add_pc_wo(C.wateroil(p), p, pc)
    pdf = R.build_report("test", wo, 0.5, 3.0, [("Swl", "0.10"), ("Sorw", "0.20")], pc_unit="psi")
    assert pdf[:5] == b"%PDF-" and len(pdf) > 10_000
    wo0 = C.wateroil(p)   # no Pc column must also work
    assert R.build_report("nopc", wo0.assign(PC=0.0), 0.5, 3.0, [])[:5] == b"%PDF-"
