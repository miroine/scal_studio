# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
import pathlib, subprocess, sys


def test_headless_smoke():
    """Runs app.py in a subprocess against stub streamlit/plotly across systems, models, units, presets and export families."""
    r = subprocess.run([sys.executable, str(pathlib.Path(__file__).parent / "smoke_app.py")], capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout[-1500:] + r.stderr[-3000:]
    assert "configuration runs passed" in r.stdout
