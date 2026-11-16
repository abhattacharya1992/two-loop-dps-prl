#!/usr/bin/env python3
"""Cross-check of the translated ggvvamp code (route A) against the FORM file of the
coefficients and the master integrals in the Wolfram Language file (route B).

Usage: crosscheck_form.py BUILD INPUTS s t ma2 mb2 DPS [coefficients...]

Route B: A_j from gg_Aj_A_qT-1.0.inc, with the master integrals M(fam,idx,w) evaluated from
vvpmasters_opt-1.0.m. The two routes share only the kinematic map, mpmath and the Li22 routine
of ggvv.py (checked separately by check_li22_values.py).
"""
import json
import argparse
import pathlib
import re
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import mpmath as mp
from formeval import Evaluator, letters
from ggvv import GGVV, Functions, li22
from numeric_checks import require_close
from zz_master_families import FAMILY

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('coefficients',nargs='*')
parser.add_argument('--allow-unresolved-ca',action='store_true',
                    help='report the documented CA discrepancy without certifying agreement')
parser.add_argument('--rtol',default='1e-20')
parser.add_argument('--output',type=pathlib.Path,help='write the compared coefficients and master values')
if len(sys.argv)<8:raise SystemExit(__doc__)
# Preserve signed/scientific-notation kinematics as positional strings.
build,inputs=map(pathlib.Path,sys.argv[1:3])
s_,t_,ma2_,mb2_=sys.argv[3:7]
dps=int(sys.argv[7])
args=parser.parse_args(sys.argv[8:])
wanted=args.coefficients or ['1l']
mp.mp.dps = dps
s, t, ma2, mb2 = (mp.mpf(v) for v in (s_, t_, ma2_, mb2_))

# ---------------------------------------------------------------- route A
g = GGVV(build)
A_route_a = {}
for label, (tt, a2, b2) in {"": (t, ma2, mb2), "x12": (ma2 + mb2 - s - t, ma2, mb2),
                            "x34": (ma2 + mb2 - s - t, mb2, ma2), "x12x34": (t, mb2, ma2)}.items():
    A_route_a[label] = g.compute(s, tt, a2, b2)
CROSS = {3: ("x12", 2, 1), 6: ("x12", 5, 1), 7: ("x12", 4, 1), 10: ("x34", 8, -1), 11: ("x34", 9, -1),
         12: ("x12", 9, 1), 13: ("x12", 8, 1), 14: ("x12x34", 9, -1), 15: ("x12x34", 8, -1),
         19: ("x12", 18, 1), 20: ("x12", 17, 1)}


def route_a(j, part):
    if j in CROSS:
        lab, k, sign = CROSS[j]
        return sign * A_route_a[lab][k][part]
    return A_route_a[""][j][part]


# ---------------------------------------------------------------- route B
u, x, y, z = g.variables(s, t, ma2, mb2)
L = letters(x, y, z)
root = mp.sqrt((s - ma2 - mb2) ** 2 - 4 * ma2 * mb2)

text = (inputs / "gg_Aj_A_qT-1.0.inc").read_text()
heads = [(m.start(), m.group(1)) for m in re.finditer(r"^g \[(A\d+_A_\dl_qT)\] =", text, re.M)]
exprs = {}
for k, (pos, name) in enumerate(heads):
    end = heads[k + 1][0] if k + 1 < len(heads) else len(text)
    body = text[pos:end].split("=", 1)[1]
    exprs[name] = body[: body.rindex(";")]

opt = (inputs / "vvpmasters_opt-1.0.m").read_text()
rules = {}
for m in re.finditer(r'M\["(\w+)", (\d+), (\d+)\] -> (.*?)\n,\n', opt + "\n,\n", re.S):
    rules[(m.group(1), int(m.group(2)), int(m.group(3)))] = m.group(4)


class MasterFunctions:
    def __init__(self):
        self.cache = {}

    def __call__(self, name, args):
        key = (name,mp.mp.prec,*args)
        if key in self.cache:
            return self.cache[key]
        if name == "Log":
            v = mp.log(args[0])
        elif name == "PolyLog":
            v = mp.polylog(int(mp.re(args[0])), args[1])
        elif name == "Li":
            if not (args[0] == 2 and args[1] == 2):
                raise ValueError("unexpected Li")
            v = li22(mp.re(args[2]), mp.re(args[3]))
        elif name == "Zeta":
            v = mp.zeta(int(mp.re(args[0])))
        else:
            raise ValueError(name)
        if not mp.isfinite(v):raise ValueError('Nonfinite FORM master function')
        self.cache[key] = v
        return v


sym_m = dict(L, I=mp.mpc(0, 1), Pi=+mp.pi)
mf = MasterFunctions()


def master(fam, idx, w):
    key = (FAMILY[fam], idx, w)
    if key not in rules:
        raise KeyError(f"no rule for {key}")
    return Evaluator(sym_m, mf).evaluate(rules[key])


class FormFunctions:
    def __call__(self, name, args):
        if name == "Den":
            return 1 / mp.mpmathify(args[0])
        if name == "Num":
            return args[0]
        if name == "Log":
            return mp.log(args[0])
        if name == "PolyLog":
            return mp.polylog(int(mp.re(args[0])), args[1])
        raise ValueError(name)


def route_b(name, CA, CF):
    body = exprs[name]
    syms = dict(L, t=t, u=u, ma2=ma2, mb2=mb2, root=root, CA=mp.mpf(CA), CF=mp.mpf(CF),
                Pi=+mp.pi, i_=mp.mpc(0, 1))
    for m in set(re.findall(r"M\(\s*(\w+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)", body)):
        syms[f"M__{m[0]}__{m[1]}__{m[2]}"] = master_values[m]
    body = re.sub(r"M\(\s*(\w+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)", r"M__\1__\2__\3", body)
    return Evaluator(syms, FormFunctions()).evaluate(body)


master_values = {}
report = []
for spec in wanted:
    if spec == "1l":
        names = [f"A{j}_A_1l_qT" for j in range(1, 21)]
    else:
        names = [f"A{spec}_A_2l_qT"]
    for name in names:
        j = int(name[1:name.index("_")])
        if name.endswith("2l_qT"):
            t0 = time.time()
            for m in set(re.findall(r"M\(\s*(\w+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)", exprs[name])):
                if m not in master_values:
                    master_values[m] = master(m[0], int(m[1]), int(m[2]))
            tm = time.time() - t0
            for part, (CA, CF) in (("2l_CA", (1, 0)), ("2l_CF", (0, 1))):
                b = route_b(name, CA, CF)
                a = route_a(j, part)
                rel = abs(a-b)/abs(a) if a else (mp.inf if b else mp.mpf(0))
                if part=='2l_CA' and args.allow_unresolved_ca:
                    if not all(mp.isfinite(v) for v in (a,b,rel)):
                        raise ValueError('Nonfinite unresolved CA comparison')
                    print(f'UNRESOLVED CA: {name}; this comparison is diagnostic only',flush=True)
                else:
                    require_close(b,a,rtol=args.rtol,label=f'{name} {part}')
                report.append((name, part, a, b, rel))
                print(f"{name} {part}: route A {mp.nstr(a, 15)}  route B {mp.nstr(b, 15)}  rel. diff {mp.nstr(rel, 3)}"
                      f"  (masters {tm:.0f} s)", flush=True)
            b0 = route_b(name, 0, 0)
            require_close(b0,0,rtol=0,atol=mp.mpf(10)**(-dps+5),label='colour-free part')
            print(f"{name} colour-free part (CA = CF = 0): {mp.nstr(b0, 5)}", flush=True)
        else:
            b = route_b(name, 0, 0)
            a = route_a(j, "1l")
            rel = abs(a-b)/abs(a) if a else (mp.inf if b else mp.mpf(0))
            require_close(b,a,rtol=args.rtol,label=name)
            report.append((name, "1l", a, b, rel))
            print(f"{name}: route A {mp.nstr(a, 15)}  route B {mp.nstr(b, 15)}  rel. diff {mp.nstr(rel, 3)}", flush=True)
print("largest relative difference:", mp.nstr(max(r[4] for r in report), 3))
if args.output:
    pair=lambda z:[mp.nstr(mp.re(z),dps),mp.nstr(mp.im(z),dps)]
    data={'decimal_inputs':[s_,t_,ma2_,mb2_],'dps':dps,
          'enforced_comparisons_passed':True,
          'unresolved_ca_allowed':args.allow_unresolved_ca,
          'comparisons':[{'coefficient':name,'part':part,'route_a':pair(a),'route_b':pair(b),
                          'relative_difference':mp.nstr(rel,dps)} for name,part,a,b,rel in report],
          'master_values':{','.join(map(str,k)):pair(v) for k,v in master_values.items()}}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(data,indent=2)+'\n')
