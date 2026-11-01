# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
import pathlib, sys, json
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
import pytest
import cases as K, presets as P

D = dict(system="Water-oil", kind="Corey", wo_swl=0.1, wo_krw_n=3.0, pc_on=True, pcwo_pd=0.3, ens_n=200)
CH = dict(system=["Water-oil", "Gas-oil"], kind=["Corey", "LET"])


def test_roundtrip_and_clamp():
    sch = K.schema_for(D, CH)
    clean, probs, meta = K.json_to_state(K.state_to_json(D, sch, "hello"), sch)
    assert clean == D and not probs and meta["note"] == "hello"


def test_bad_values_dropped():
    sch = K.schema_for(D, CH)
    clean, probs = K.sanitize(dict(wo_swl=9.0, kind="Bogus", pc_on="yes", nope=1, ens_n=250.0, wo_krw_n="x"), sch)
    assert clean == dict(ens_n=250) and len(probs) == 5 and isinstance(clean["ens_n"], int)


def test_non_case_file_rejected():
    sch = K.schema_for(D, CH)
    for txt in ("not json", json.dumps({"app": "other", "state": {}})):
        with pytest.raises(ValueError):
            K.json_to_state(txt, sch)


def test_missing_rule_is_loud():
    with pytest.raises(KeyError):
        K.schema_for(dict(mystery_value=1.0), {})


def test_presets_use_valid_values():
    """Every preset key must be a numeric/choice key with an in-bounds value (guards the widgets against bad presets)."""
    for name, pr in P.PRESETS.items():
        for k, v in pr["state"].items():
            if isinstance(v, (bool, str)):
                continue
            b = K.bounds_for(k)
            assert b is not None, f"{name}: no bounds for {k}"
            assert b[0] <= v <= b[1], f"{name}: {k}={v} outside {b}"


def test_app_defaults_and_presets_match_schema():
    import state_defs as S
    sch = K.schema_for(S.DEFAULTS, S.CHOICES)           # raises if any key lacks a bounds rule or choices
    clean, probs = K.sanitize(S.DEFAULTS, sch)
    assert not probs and len(clean) == len(S.DEFAULTS)
    for name, pr in P.PRESETS.items():
        clean, probs = K.sanitize(pr["state"], sch)
        assert not probs and len(clean) == len(pr["state"]), (name, probs)
