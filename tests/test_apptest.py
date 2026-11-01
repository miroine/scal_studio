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
