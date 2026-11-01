#!/usr/bin/env python
# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
"""Command-line batch runner:  python scal_batch.py input.yaml --out results/ [--engine native|pyscal]"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import batch as B  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="SCAL Studio batch mode: YAML in, simulator-ready saturation tables out.")
    ap.add_argument("config", nargs="?", help="YAML input file")
    ap.add_argument("--out", default="scal_out", help="output directory (default scal_out)")
    ap.add_argument("--engine", choices=B.ENGINES, help="override the engine set in the YAML")
    ap.add_argument("--template", action="store_true", help="print an example YAML and exit")
    ap.add_argument("--check", action="store_true", help="validate the YAML only; write nothing")
    ap.add_argument("--list-profiles", action="store_true", help="list benchmark profiles usable in 'advisor:' and exit")
    a = ap.parse_args(argv)
    if a.template:
        print(B.TEMPLATE)
        return 0
    if a.list_profiles:
        import advisor
        for k, p in advisor.PROFILES.items():
            print(f"{k:28s} {p.title}")
        return 0
    if not a.config:
        ap.error("config file required (or --template)")
    try:
        text = Path(a.config).read_text()
        cfg = B.load_config(text)
        if a.check:
            print(f"OK: {len(cfg['cases'])} case(s), engine {a.engine or cfg['engine']}")
            return 0
        res = B.run_batch(text, a.engine)
    except (B.BatchError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    out = Path(a.out)
    for name, blob in res.files.items():
        dest = out / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(blob)
    print(res.summary.to_string(index=False))
    for w in res.warnings:
        print("warning:", w)
    for r in res.checks:
        print(f"verify {r['file']}: {r['passed']}/{r['total']}")
    print(f"wrote {len(res.files)} file(s) to {out}/" + ("" if res.ok else "  -- VERIFICATION FAILED"))
    return 0 if res.ok else 1


if __name__ == "__main__":
    sys.exit(main())
