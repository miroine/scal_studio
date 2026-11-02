# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
"""Real Streamlit run with streamlit.testing.v1.AppTest (needs streamlit + plotly installed; skipped otherwise)."""
import pathlib
import pytest

pytest.importorskip("streamlit")
pytest.importorskip("plotly")
from streamlit.testing.v1 import AppTest

APP = str(pathlib.Path(__file__).resolve().parent.parent / "app.py")


def fresh():
    return AppTest.from_file(APP, default_timeout=90).run()


def test_default_run_has_no_exception():
    at = fresh()
    assert not at.exception, at.exception
    assert len(at.tabs) >= 11


@pytest.mark.parametrize("system", ["Water-oil", "Gas-oil", "Gas-water", "Three-phase (oil-water-gas)"])
def test_each_system_runs(system):
    at = fresh()
    at.radio(key="system").set_value(system).run()
    assert not at.exception, at.exception


def test_units_switch_and_pdf():
    at = fresh()
    at.selectbox(key="pc_unit").set_value("psi").run()
    assert not at.exception, at.exception
    at.button(key="pdf_btn").click().run()
    assert not at.exception, at.exception
    assert at.session_state["_pdf"][:5] == b"%PDF-"


def test_preset_button():
    at = fresh()
    at.selectbox(key="preset_sel").set_value("Oil-wet carbonate, waterflood").run()
    at.button(key="preset_btn").click().run()
    assert not at.exception, at.exception
    assert at.session_state["wo_krwend"] == pytest.approx(0.75)


@pytest.mark.parametrize("kind", ["Corey", "LET"])
def test_direct_downloads_include_gas_and_exact_rows(monkeypatch, kind):
    import exporters as X
    import deckcheck as D

    written = []
    export = X.export_eclipse

    def capture(regions, *args, **kwargs):
        txt = export(regions, *args, **kwargs)
        written.append((regions, txt))
        return txt

    monkeypatch.setattr(X, "export_eclipse", capture)
    at = fresh()
    at.radio(key="system").set_value("Three-phase (oil-water-gas)")
    at.radio(key="kind").set_value(kind)
    at.number_input(key="export_rows").set_value(23)
    written.clear()
    at.run()
    assert not at.exception, at.exception
    direct = [(r, txt) for r, txt in written if r[0]["name"].startswith(("Advisor interpolated", "Interpolated t="))]
    assert len(direct) == 2
    for regions, txt in direct:
        tabs = D.tables_from_eclipse(txt)
        assert set(tabs) == {"SWOF", "SGOF"}
        assert all(len(table) == 23 for ts in tabs.values() for table in ts)
        assert all(c.ok for c in D.verify_eclipse(txt, regions))
        assert regions[0]["provenance"]["go"]["parameters"]["krg"]["kind"] == kind.lower()
        assert "pc_scale" in txt and "resolution_dS" in txt
    assert any("SWOF + SGOF" in button.label for button in at.get("download_button"))
    preview = next(block.value for block in at.code if "Saturation functions -" in block.value)
    assert "\nSWOF\n" in preview and "\nSGOF\n" in preview

    written.clear()
    at.checkbox(key="pc_on").uncheck().run()
    assert not at.exception, at.exception
    assert all(r[0]["provenance"]["wo"]["pc"]["kind"] == "disabled" for r, _ in written)


def test_derived_exports_have_effective_parameters(monkeypatch):
    import exporters as X

    written = []
    export = X.export_eclipse

    def capture(regions, *args, **kwargs):
        written.append(regions)
        return export(regions, *args, **kwargs)

    monkeypatch.setattr(X, "export_eclipse", capture)
    at = fresh()
    source = next(s for s in at.selectbox if s.label.startswith("Region source"))
    for option in source.options:
        written.clear()
        source.set_value(option).run()
        assert not at.exception, at.exception
        regions = written[-1]
        for r in regions:
            p = r["provenance"]["wo"]["parameters"]
            assert r["wo"]["SW"].min() == pytest.approx(p["swl"])
            assert r["wo"]["KRW"].max() == pytest.approx(p["krwend"])
            assert r["wo"]["KROW"].max() == pytest.approx(p["kroend"])
            assert len(r["wo"]) == at.session_state["export_rows"]
