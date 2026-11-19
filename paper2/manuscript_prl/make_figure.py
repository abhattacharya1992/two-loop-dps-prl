#!/usr/bin/env python3
# The figures of the two documents.
#
# Supplement: the figures with results, drawn by the figure script of the full manuscript
# (../_build/make_figure.py) from the result files listed in its SOURCES, with the output
# directed to figures/ here.
#
# Letter, Fig. 1: the two-loop DPS configuration for gluon and quark splittings, drawn with
# TikZ-Feynman. The TikZ source is figures/dps_diagram_tikz.tex, compiled inside the wrapper
# figures/dps_diagram.tex.
#
#   python3 make_figure.py
import importlib.util
import pathlib
import subprocess
import sys
import argparse

HERE = pathlib.Path(__file__).resolve().parent
FULL = HERE.parent/'_build'

sys.path.insert(0, str(FULL))        # make_figure.py imports make_tables from its own folder
spec = importlib.util.spec_from_file_location('make_figure_full', FULL/'make_figure.py')
mf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mf)
mf.OUT = HERE/'figures'
mf.OUT.mkdir(exist_ok=True)


def diagram_figure():
    """Build figures/dps_diagram.pdf from the retained local TikZ source."""
    picture = HERE/'figures'/'dps_diagram_tikz.tex'
    if not picture.exists():
        sys.exit(f'Required local figure source is missing: {picture}')
    run = subprocess.run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error', 'dps_diagram.tex'],
                         cwd=HERE/'figures', capture_output=True, text=True)
    if run.returncode != 0:
        sys.exit(run.stdout[-2000:])
    for by_product in ('dps_diagram.aux', 'dps_diagram.log'):
        (HERE/'figures'/by_product).unlink(missing_ok=True)
    print('wrote figures/dps_diagram.pdf')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reproduce the PRL v2 figures')
    parser.add_argument('--diagnostics', action='store_true',
                        help='also render forward-approach and ZZ-coefficient diagnostic plots')
    args = parser.parse_args()
    if args.diagnostics:
        mf.forward_figure()
        mf.zz_figure()
    mf.zz_prediction_figure()
    mf.three_process_comparison_figure()
    mf.windows_figure()
    diagram_figure()
    print('Publication figures regenerated'
          + (' with optional diagnostic plots' if args.diagnostics else ''))
