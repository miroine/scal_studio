# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
import pathlib, sys, io, zipfile, json, subprocess
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
import pytest
import batch as B, core as C, advisor as A, deckcheck as D

ROOT = pathlib.Path(__file__).resolve().parent.parent


def test_template_runs_and_verifies():
    r = B.run_batch(B.TEMPLATE)
    assert r.ok and len(r.summary) == 3 and set(r.files) >= {"scal_eclipse.inc", "scal_cmg.inc", "scal_fluid.m", "scal_tables.csv", "summary.csv", "manifest.json"}
    assert all(c["passed"] == c["total"] for c in r.checks)
    with zipfile.ZipFile(io.BytesIO(r.zip_bytes())) as z:
        assert "manifest.json" in z.namelist()
    assert json.loads(r.files["manifest.json"])["ok"] is True


def test_advisor_t_matches_direct_interpolation():
    y = "cases:\n  - {name: a, advisor: {profile: oil_wet_sandstone, t: 0.25}}\n"
    cfg = B.load_config(y)
    lo, ba, hi = (A.make_params(A.PROFILES["oil_wet_sandstone"], c) for c in ("low", "base", "high"))
    ref = C.interpolate_triple(lo, ba, hi, 0.25)
    assert cfg["cases"][0].wo.sorw == pytest.approx(ref.sorw) and cfg["cases"][0].wo.krw.n == pytest.approx(ref.krw.n)


def test_batch_inc_records_model_and_parameters():
    r = B.run_batch("""
cases:
  - name: reproducible
    system: three-phase
    kr_model: let
    gas_oil: {sgcr: 0.07, krg: {L: 1.7, E: 3.2, T: 2.4}}
    pc_multiplier: 1.4
    pc: {gas_oil: {kind: let, pcmax: 2.3}}
""")
    txt = r.files["scal_eclipse.inc"].decode()
    assert r.ok
    assert "SGOF\n" in txt and "SWOF\n" in txt
    assert 'go.parameters.krg.kind="let"' in txt
    assert "go.parameters.krg.E=3.2" in txt
    assert "go.pc.pcmax=2.3" in txt
    assert "go.pc_scale=1.4" in txt


def test_interpolate_between_cases_and_endpoints():
    y = """
cases:
  - {name: lo, water_oil: {swl: 0.1, sorw: 0.3, krwend: 0.6}}
  - {name: hi, water_oil: {swl: 0.1, sorw: 0.1, krwend: 0.2}}
  - {name: mid, interpolate: {from: lo, to: hi, t: 0.5}}
"""
    c = B.load_config(y)["cases"]
    assert c[2].wo.sorw == pytest.approx(0.2) and c[2].wo.krwend == pytest.approx(0.4)


@pytest.mark.parametrize("units,factor", [("METRIC", 1.0), ("FIELD", C.BAR_TO_PSI), ("SI-kPa", C.BAR_TO_KPA)])
def test_units_applied_in_files(units, factor):
    y = f"units: {units}\nexport: {{formats: [eclipse]}}\ncases:\n  - name: a\n    pc: {{water_oil: {{kind: brooks-corey, pd: 0.3, pcmax: 5.0}}}}\n"
    r = B.run_batch(y)
    assert r.ok
    tab = D.parse_eclipse(r.files["scal_eclipse.inc"].decode())["SWOF"][0]
    assert tab[0, 3] == pytest.approx(5.0 * factor, rel=1e-5)


def test_three_phase_gas_oil_gas_water_and_split_output():
    y = """
export: {combine: false, formats: [eclipse, csv], eclipse_family: SWFN/SGFN/SOF3}
cases:
  - {name: tp, system: three-phase, kr_model: let}
  - {name: gw, system: gas-water, gas_water: {sgr: 0.3}}
  - {name: go, system: gas-oil}
"""
    r = B.run_batch(y)
    assert r.ok and "cases/tp/scal_eclipse.inc" in r.files and "cases/gw/scal_tables.csv" in r.files


@pytest.mark.parametrize("bad,frag", [
    ("cases: []", "non-empty"), ("cases:\n  - {name: a, water_oil: {swl: 0.1, bogus: 1}}", "unknown key"),
    ("cases:\n  - {name: a, water_oil: {swl: 2.0}}", "outside"), ("cases:\n  - {name: a, system: foo}", "not one of"),
    ("engine: magic\ncases:\n  - {name: a}", "engine"), ("cases:\n  - {name: a}\n  - {name: a}", "duplicate"),
    ("cases:\n  - {name: a, water_oil: {swl: 0.5, swcr: 0.5, sorw: 0.6}}", "must be"),
    ("cases:\n  - {name: a, advisor: {profile: nope}}", "not one of"),
    ("cases:\n  - {name: a, interpolate: {from: x, to: y}}", "defined earlier"),
    ("cases:\n  - {name: a, water_oil: {krw: {L: 2}}, kr_model: corey}", "unknown key"),
    ("cases:\n  - {name: a, pc: {water_oil: {kind: brooks-corey, pd: 9, pcmax: 1}}}", "Pcmax"),
    ("cases: [unclosed", "YAML syntax"),
])
def test_invalid_inputs_are_rejected_with_path(bad, frag):
    with pytest.raises(B.BatchError) as e:
        B.load_config(bad)
    assert frag in str(e.value)


def test_pyscal_engine_errors_cleanly_when_missing():
    import pyscal_bridge as PB
    if PB.PYSCAL_AVAILABLE:
        pytest.skip("pyscal installed")
    with pytest.raises(B.BatchError, match="pyscal"):
        B.run_batch("engine: pyscal\ncases:\n  - {name: a}")
    with pytest.raises(B.BatchError, match="gas-water"):
        B.load_config("engine: pyscal\ncases:\n  - {name: a, system: gas-water}")


def test_cli(tmp_path):
    y = tmp_path / "in.yaml"; y.write_text(B.TEMPLATE)
    r = subprocess.run([sys.executable, str(ROOT / "scal_batch.py"), str(y), "--out", str(tmp_path / "o")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr + r.stdout
    assert (tmp_path / "o" / "scal_eclipse.inc").exists() and "verify scal_eclipse.inc: 7/7" in r.stdout
    y.write_text("cases: []")
    r = subprocess.run([sys.executable, str(ROOT / "scal_batch.py"), str(y), "--check"], capture_output=True, text=True)
    assert r.returncode == 2 and "ERROR" in r.stderr
