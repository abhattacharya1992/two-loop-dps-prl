# DPS splitting diagrams

The [DPS splitting diagram](dps_sketch.pdf) shows the two configurations considered in the paper:

- **Gluon splitting:** The two incoming gluons split into gluon pairs, which contribute to $gg\to HH$ and $gg\to AA$ through the effective Higgs and pseudoscalar vertices.
- **Quark splitting:** The incoming gluons split into quark-antiquark pairs, which contribute to $gg\to ZZ$ through massless quark loops.

The diagram is provided as a reference for the common DPS configuration. It is **not included in the PRL Letter or Supplemental Material**, and it is excluded from the arXiv submission archive.

The original drawing is in [`dps_sketch.tex`](dps_sketch.tex). To reproduce the PDF, run in this directory:

```bash
pdflatex -interaction=nonstopmode -halt-on-error dps_sketch.tex
```

The source uses TikZ and the `standalone` document class.
