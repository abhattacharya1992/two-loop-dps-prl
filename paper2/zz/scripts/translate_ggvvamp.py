#!/usr/bin/env python3
"""Translate the generated C++ expression files of ggvvamp 1.0 (von Manteuffel, Tancredi,
arXiv:1503.08835) into Python code objects that run with mpmath at any precision.

Usage: translate_ggvvamp.py GGVVAMP_SOURCE_DIR OUTPUT_DIR

The generated code has three parts:
  p0002.cpp        letters, function arguments and the polylogarithms (via GiNaC)
  p0003..p2463.cpp straight-line arithmetic on the array r[] (master integrals, coefficients)
  p2464.cpp        assignment of the amplitude coefficients A[j][loop][re/im]
p0002 is translated into one Python function that calls mpmath for the special functions;
p0003..p2463 into one function per file; p2464 is read into a table of r[] indices.
Nothing from the source directory is executed or imported.
"""
import json
import marshal
import pathlib
import re
import sys

src, out = map(pathlib.Path, sys.argv[1:3])
out.mkdir(parents=True, exist_ok=True)

STMT = re.compile(r"^r\[(\d+)\] = (.*);$")
RAT = re.compile(r"n\((-?\d+),(\d+)\)")
INT = re.compile(r"n\((-?\d+)\)")
SAFE = re.compile(r"^[-+*/() 0-9rQ\[\]_]*$")


def arith(expr, rationals):
    """C++ arithmetic on r[] with n(a), n(a,b) -> Python with ints and Q[k] rationals."""
    def rat(m):
        key = f"{m.group(1)}/{m.group(2)}"
        if key not in rationals:
            rationals[key] = len(rationals)
        return f"Q[{rationals[key]}]"
    e = RAT.sub(rat, expr)
    e = INT.sub(lambda m: f"({m.group(1)})", e)
    if not SAFE.match(e):
        raise ValueError(f"unexpected syntax: {expr[:200]}")
    return e


rationals = {}

# --- p0002: letters, arguments, functions -------------------------------------------------
LI22 = re.compile(r"^n\(Li\(lst\(2,2\),lst\(tonumeric\((.*)\),tonumeric\((.*)\)\)\)\)$")
LIN = re.compile(r"^n\(Li\(([234]),tonumeric\((.*)\)\)\)$")
LOG = re.compile(r"^n\(log\(tonumeric\((.*)\)\)\)$")
ZETA = re.compile(r"^n\(zeta\(([234])\)\)$")
lines = ["def p0002(r, Q, F):"]
nfun = 0
for line in (src / "p0002.cpp").read_text().splitlines():
    m = STMT.match(line.strip())
    if not m:
        continue
    i, e = m.groups()
    if (f := LI22.match(e)):
        rhs = f"F.li22({arith(f.group(1), rationals)}, {arith(f.group(2), rationals)})"
        nfun += 1
    elif (f := LIN.match(e)):
        rhs = f"F.lin({f.group(1)}, {arith(f.group(2), rationals)})"
        nfun += 1
    elif (f := LOG.match(e)):
        rhs = f"F.log({arith(f.group(1), rationals)})"
        nfun += 1
    elif (f := ZETA.match(e)):
        rhs = f"F.zeta({f.group(1)})"
    else:
        rhs = arith(e, rationals)
    lines.append(f"    r[{i}] = {rhs}")
p0002_src = "\n".join(lines) + "\n"
(out / "p0002.py").write_text(p0002_src)

# --- p0003..p2463: arithmetic --------------------------------------------------------------
codes = []
nstmt = 0
for k in range(3, 2464):
    body = []
    for line in (src / f"p{k:04d}.cpp").read_text().splitlines():
        m = STMT.match(line.strip())
        if m:
            body.append(f"    r[{m.group(1)}] = {arith(m.group(2), rationals)}")
            nstmt += 1
        elif line.startswith("r["):
            raise ValueError(f"unparsed statement in p{k:04d}: {line[:120]}")
    text = f"def p{k:04d}(r, Q):\n" + ("\n".join(body) if body else "    pass") + "\n"
    codes.append(compile(text, f"p{k:04d}", "exec"))
(out / "algebra.marshal").write_bytes(marshal.dumps(codes))

# --- p2464: which r[] entries are the coefficients ----------------------------------------
table = {}
for line in (src / "p2464.cpp").read_text().splitlines():
    m = re.match(r"^A\[(\d+)\]\[(\d)\]\[(\d)\] = (.*);$", line.strip())
    if not m:
        continue
    j, loop, reim, e = int(m.group(1)), int(m.group(2)), int(m.group(3)), m.group(4)
    pi = e.startswith("n(GiNaC::Pi)*")
    e = e.removeprefix("n(GiNaC::Pi)*").strip("()")
    terms = {}
    for t in e.split("+"):
        f = re.fullmatch(r"(?:(CF|CA)\*)?r\[(\d+)\]", t.strip())
        if not f:
            raise ValueError(f"unexpected coefficient form: {line}")
        terms[f.group(1) or "1"] = int(f.group(2))
    table[f"{j},{loop},{reim}"] = {"pi": pi, "terms": terms}
(out / "coefficients.json").write_text(json.dumps(table, indent=1))
(out / "rationals.json").write_text(json.dumps(rationals))
print(f"p0002: {nfun} function calls; algebra: {nstmt} statements in {len(codes)} files; "
      f"{len(rationals)} distinct rationals; {len(table)} coefficient entries")
