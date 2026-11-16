import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import mpmath as mp
from ggvv import GGVV
build = sys.argv[1]
ns = {}
exec(compile(pathlib.Path(build, "p0002.py").read_text(), "p0002", "exec"), ns)
g = GGVV.__new__(GGVV)
mp.mp.dps = 30
class Probe:
    def __init__(s): s.pairs = []
    def li22(s, a, b): s.pairs.append((a, b)); return mp.mpf(0)
    def lin(s, n, x): return mp.mpf(0)
    def log(s, x): return mp.mpf(0)
    def zeta(s, n): return mp.mpf(0)
import json
rat = json.loads(pathlib.Path(build, "rationals.json").read_text())
Q = [None]*len(rat)
for k, v in rat.items():
    a, b = k.split("/"); Q[v] = mp.mpf(int(a))/int(b)
for (s, t, ma2, mb2) in [("15625", "-2384.89727896833", "6467.215561", "625"), ("10", "-2.1", "1", "1"), ("10", "-1.0001", "1", "1")]:
    s, t, ma2, mb2 = map(mp.mpf, (s, t, ma2, mb2))
    u, x, y, z = GGVV.variables(g, s, t, ma2, mb2)
    r = [None]*2000
    r[1], r[2], r[3], r[4], r[5], r[6], r[7] = ma2, t/ma2, u/ma2, mb2/ma2, x, y, z
    P = Probe(); ns["p0002"](r, Q, P)
    given = sum(1 for a, b in P.pairs if abs(a) < 1 and abs(a*b) < 1)
    swapped = sum(1 for a, b in P.pairs if abs(b) < 1 and abs(a*b) < 1)
    print(f"s={s} t={t} ma2={ma2} mb2={mb2}: {len(P.pairs)} Li22 calls; convergent as given: {given}; convergent if swapped: {swapped}; max|ab| = {mp.nstr(max(abs(a*b) for a,b in P.pairs), 6)}, max|a| = {mp.nstr(max(abs(a) for a,b in P.pairs), 6)}")
