# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pytest

import core as C
import deckcheck as D
import exporters as X


def region():
    p = C.WOParams(swl=0.13, swcr=0.21, sorw=0.23)
    g = C.GOParams(swl=p.swl, sgcr=0.07, sorg=0.19)
    pc = C.PcParams(kind="let", L=1.7, E=3.1, T=2.2)
    return dict(
        name="Derived case",
        wo=C.add_pc_wo(C.wateroil(p), p, pc, 1.4),
        go=C.add_pc_go(C.gasoil(g), g, pc),
        provenance={
            "wo": X.model_provenance(p, pc, resolution=0.01, pc_scale=1.4, source="Advisor t=0.37"),
            "go": X.model_provenance(g, pc, resolution=0.01),
        },
    )


@pytest.mark.parametrize("count", [8, 23, 100])
@pytest.mark.parametrize("family", ["SWOF/SGOF", "SWOF/SLGOF", "SWFN/SGFN/SOF3"])
@pytest.mark.parametrize("units", ["METRIC", "FIELD", "SI-kPa"])
def test_exact_rows_and_roundtrip(count, family, units):
    original = region()
    r = X.prepare_regions([original], count)
    txt = X.export_eclipse(r, units=units, family=family)
    parsed = D.tables_from_eclipse(txt)
    assert all(len(df) == count for tables in parsed.values() for df in tables)
    assert all(c.ok for c in D.verify_eclipse(txt, r, units=units, family=family))
    for key, sat, endpoints in (
        ("wo", "SW", [0.13, 0.21, 0.77, 1.0]),
        ("go", "SG", [0.0, 0.07, 0.68, 0.87]),
    ):
        table = r[0][key]
        assert len(table) == count
        assert np.all(np.diff(table[sat]) > 0)
        for endpoint in endpoints:
            assert np.any(np.isclose(table[sat], endpoint, atol=1e-10))
        for col in table.columns[1:]:
            np.testing.assert_allclose(table[col], np.interp(table[sat], original[key][sat], original[key][col]))
    assert "sampling" not in original
    assert len(original["wo"]) != count
    assert all(c.ok for c in D.verify_cmg(X.export_cmg(r), r))
    assert all(c.ok for c in D.verify_mrst(X.export_mrst(r), r))


def test_reproducibility_comments_and_json():
    r = X.prepare_regions([region()], 23)
    txt = X.export_eclipse(r, units="FIELD")
    assert "SGOF\n" in txt and "SWOF\n" in txt
    assert "family=SWOF/SGOF" in txt and "units=FIELD" in txt
    assert 'krw.kind="corey"' in txt and "krw.n=3.0" in txt
    assert 'pc.kind="let"' in txt and "pc.L=1.7" in txt
    assert "pc_scale=1.4" in txt and "swcr=0.21" in txt
    assert "Advisor t=0.37" in txt and "rows_per_table=23" in txt
    payload = json.loads(X.export_json(r))
    assert payload["regions"][0]["provenance"] == r[0]["provenance"]


def test_gaswater_and_imbibition_sampling():
    p = C.GWParams()
    gw = C.add_pc_gw(C.gaswater(p), p, C.PcParams())
    r = X.prepare_regions([dict(gw=gw)], 17)
    assert all(len(t) == 17 for ts in D.tables_from_eclipse(X.export_eclipse(r)).values() for t in ts)
    assert all(c.ok for c in D.verify_eclipse(X.export_eclipse(r), r))
    original = region()
    hs = C.HystSetup(original["go"].SG.values, original["go"].KRG.values, 0.07, 0.68, 0.3)
    original["imb_go"] = C.imbibition_gas_table(original["go"], hs, "killough").assign(PC=0.0)
    r = X.prepare_regions([original], 17)
    hyst = dict(scheme="killough", s_trap_max=0.3)
    txt = X.export_eclipse(r, hysteresis=hyst)
    assert [len(t) for t in D.parse_eclipse(txt)["SGOF"]] == [17, 17]
    assert all(c.ok for c in D.verify_eclipse(txt, r, hysteresis=hyst))


@pytest.mark.parametrize("sorw", [0.17, 0.21, 0.27, 0.31])
def test_model_endpoints_survive_floating_point_tails(sorw):
    p = C.WOParams(swcr=0.18, sorw=sorw)
    original = dict(wo=C.wateroil(p).assign(PC=0.0), provenance={"wo": X.model_provenance(p)})
    r = X.prepare_regions([original], 23)[0]["wo"]
    assert np.any(np.isclose(r.SW, 1 - sorw, atol=1e-10, rtol=0))
    assert np.any(np.isclose(r.SW, p.swcr, atol=1e-10, rtol=0))


@pytest.mark.parametrize("count", [0, 1, 2.5, True])
def test_invalid_row_counts(count):
    with pytest.raises(ValueError, match="rows"):
        X.prepare_regions([region()], count)
