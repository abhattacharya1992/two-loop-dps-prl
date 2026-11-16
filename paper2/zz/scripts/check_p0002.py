#!/usr/bin/env python3
"""Check of the translation of p0002.cpp, the polylogarithm part of ggvvamp, into Python.

Usage: check_p0002.py SRC BUILD_GGV
       SRC        the original code (build/ggvvamp-1.0)
       BUILD_GGV  the translation (build/ggv)

p0002.cpp is a list of statements r[i] = <expression>; with rational numbers n(k) and the GiNaC
calls Li(lst(2,2),lst(tonumeric(a),tonumeric(b))) for Li22(a, b), Li(n,tonumeric(x)) for Li_n(x),
log(tonumeric(x)) and zeta(n). This script transcribes every statement on its own, with regular
expressions independent of translate_ggvvamp.py, into the form used by the generated p0002.py
(F.li22, F.lin, F.log, F.zeta), and compares the two statement by statement, spaces removed.
The other checks do not cover p0002 completely: the compiled C++ code (check_cxx.py) takes its
values from the translation, and the FORM file (route B) checks only the values that enter the
one-loop and C_F terms. The functions themselves are mpmath's, and Li22 is checked by
check_li22_values.py.
"""
import collections
import pathlib
import re
import sys

LI22 = re.compile(r"^n\(Li\(lst\(2,2\),lst\(tonumeric\((.*?)\),tonumeric\((.*?)\)\)\)\)$")
LIN = re.compile(r"^n\(Li\((\d+),tonumeric\((.*)\)\)\)$")
LOG = re.compile(r"^n\(log\(tonumeric\((.*)\)\)\)$")
ZETA = re.compile(r"^n\(zeta\((\d+)\)\)$")
RATIONAL = re.compile(r"n\((-?\d+)\)")


def transcribe(rhs):
    """C++ right-hand side -> (kind, Python right-hand side)"""
    for kind, pat, fmt in (("Li22", LI22, "F.li22({0}, {1})"), ("Li_n", LIN, "F.lin({0}, {1})"),
                           ("log", LOG, "F.log({0})"), ("zeta", ZETA, "F.zeta({0})")):
        m = pat.match(rhs)
        if m:
            return kind, fmt.format(*(RATIONAL.sub(r"(\1)", g) for g in m.groups()))
    if any(name in rhs for name in ("Li(", "log(", "zeta(")):
        raise ValueError(f"Unrecognised special-function statement: {rhs}")
    return "algebraic", RATIONAL.sub(r"(\1)", rhs)


def main():
    src, build = map(pathlib.Path, sys.argv[1:3])
    cpp, kinds = {}, collections.Counter()
    for line in (src / "p0002.cpp").read_text().splitlines():
        m = re.match(r"^r\[(\d+)\] = (.*);$", line.strip())
        if m:
            kind, py = transcribe(m.group(2))
            index = int(m.group(1))
            if index in cpp:
                raise ValueError(f"Duplicate C++ assignment r[{index}]")
            cpp[index] = py
            kinds[kind] += 1
        elif re.match(r"\s*r\s*\[", line):
            raise ValueError(f"Unparsed C++ assignment: {line}")
    gen = {}
    for line in (build / "p0002.py").read_text().splitlines():
        m = re.match(r"^\s+r\[(\d+)\] = (.*)$", line)
        if m:
            index = int(m.group(1))
            if index in gen:
                raise ValueError(f"Duplicate Python assignment r[{index}]")
            gen[index] = m.group(2)
        elif re.match(r"\s*r\s*\[", line):
            raise ValueError(f"Unparsed Python assignment: {line}")
    expected = 1613  # p0002.cpp in the checked ggvvamp 1.0 source
    if len(cpp) != expected or len(gen) != expected:
        raise ValueError(f"Expected {expected} statements; found {len(cpp)} C++ and {len(gen)} Python")
    norm = lambda s: s.replace(" ", "")
    bad = sorted(k for k in set(cpp) | set(gen) if norm(cpp.get(k, "?")) != norm(gen.get(k, "!")))
    for k in bad[:10]:
        print(f"r[{k}]: C++ -> {cpp.get(k)}\n        Python   {gen.get(k)}")
    print(f"p0002.cpp: {len(cpp)} statements ({', '.join(f'{v} {k}' for k, v in sorted(kinds.items()))}); "
          f"p0002.py: {len(gen)} statements; statements that differ: {len(bad)}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
