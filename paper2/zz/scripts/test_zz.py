#!/usr/bin/env python3
"""Re-analyse the stored points and logs and compare with the stored results (exact match):
results_zz.json from analyse_zz.py, summary.md and summary.json from make_summary.py,
results_windows.json from analyse_windows.py, results_forward.json from forward_curve.py and
results_wide.json from wide_curve.py (with the points in points_wide/ next to POINTS).

Usage: test_zz.py POINTS POINTS_ROUTE_B LOGS RESULTS_JSON SUMMARY_MD"""
import json
import pathlib
import subprocess
import sys
import tempfile

here = pathlib.Path(__file__).parent
points, routeb, logs = sys.argv[1:4]
stored, summary = pathlib.Path(sys.argv[4]), pathlib.Path(sys.argv[5])
ok = True
with tempfile.TemporaryDirectory() as d:
    out = pathlib.Path(d) / "r.json"
    subprocess.run([sys.executable, str(here / "analyse_zz.py"), points, str(out)], check=True,
                   stdout=subprocess.DEVNULL)
    same = json.loads(out.read_text()) == json.loads(stored.read_text())
    print("OK: the stored points reproduce results_zz.json" if same else "FAILED: results_zz.json differs")
    ok &= same
    md = pathlib.Path(d) / "summary.md"
    subprocess.run([sys.executable, str(here / "make_summary.py"), str(out), points, routeb, logs, str(md)],
                   check=True, stdout=subprocess.DEVNULL)
    win, fwd = pathlib.Path(d) / "results_windows.json", pathlib.Path(d) / "results_forward.json"
    wide = pathlib.Path(d) / "results_wide.json"
    subprocess.run([sys.executable, str(here / "wide_curve.py"), points, str(pathlib.Path(points).parent / "points_wide"),
                    str(wide)], check=True, stdout=subprocess.DEVNULL)
    subprocess.run([sys.executable, str(here / "analyse_windows.py"), points, str(win)], check=True,
                   stdout=subprocess.DEVNULL)
    subprocess.run([sys.executable, str(here / "forward_curve.py"), points, str(out), str(fwd)], check=True,
                   stdout=subprocess.DEVNULL)
    for new, old in ((md, summary), (md.with_suffix(".json"), summary.with_suffix(".json")),
                     (win, stored.with_name("results_windows.json")), (fwd, stored.with_name("results_forward.json")),
                     (wide, stored.with_name("results_wide.json"))):
        same = new.read_text() == old.read_text()
        print(f"OK: {old.name} reproduced" if same else f"FAILED: {old.name} differs")
        ok &= same
sys.exit(0 if ok else 1)
