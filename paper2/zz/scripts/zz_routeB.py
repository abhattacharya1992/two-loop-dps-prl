#!/usr/bin/env python3
"""Route B for forward gg -> ZZ points: the coefficients A_j from the FORM file
gg_Aj_A_qT-1.0.inc with the master integrals of vvpmasters_opt-1.0.m, for comparison with
the translated ggvvamp code (route A, zz_points.py).

Usage: zz_routeB.py INPUTS OUTDIR BETA K DPS PARTS [--recompute]
       PARTS: comma-separated subset of 1l,2l_CA,2l_CF

Writes the 20 coefficients and the 36 helicity amplitudes of the requested parts.
Both the point and the split-input cache require matching provenance. Use a fresh
OUTDIR to preserve an older run, or --recompute to rebuild the split cache and point.
"""
import argparse
import json
import pathlib
import re
import sys
import tempfile
import time

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import mpmath as mp
from formeval import Evaluator, letters
from ggvv import GGVV, li22
from zz_points import kinematics, helicity_amplitudes
from zz_cache import cache_provenance, form_inputs, form_provenance, validate_point
from zz_master_families import FAMILY

MREF = re.compile(r"M\(\s*(\w+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)")


def split_inputs(inputs, cache, *, recompute=False):
    """Split actual input bytes once, validating provenance and every cached output."""
    inputs, cache = pathlib.Path(inputs), pathlib.Path(cache)
    metadata = form_provenance(inputs, "zz-form-split", {"format": 1})
    manifest = cache / "manifest.json"
    if cache.exists() and any(cache.iterdir()) and not recompute:
        if not manifest.is_file():
            raise cache_provenance.CacheMismatch(
                f"{cache}: untagged legacy split cache; use a fresh output directory "
                "or explicitly pass --recompute.")
        record = cache_provenance.load_result(manifest, metadata)
        outputs = record.get("outputs")
        if not isinstance(outputs, dict) or "rules.json" not in outputs:
            raise cache_provenance.CacheMismatch(f"{cache}: incomplete split manifest; pass --recompute.")
        actual_names = {p.name for p in cache.iterdir()} - {manifest.name}
        if actual_names != set(outputs):
            raise cache_provenance.CacheMismatch(f"{cache}: changed split files; pass --recompute.")
        actual = cache_provenance.file_hashes({name: cache / name for name in outputs})
        if actual != outputs:
            raise cache_provenance.CacheMismatch(f"{cache}: changed split contents; pass --recompute.")
        return metadata
    text = (inputs / "gg_Aj_A_qT-1.0.inc").read_text()
    heads = [(m.start(), m.group(1)) for m in re.finditer(r"^g \[(A\d+_A_\dl_qT)\] =", text, re.M)]
    if not heads or len({name for _, name in heads}) != len(heads):
        raise ValueError("FORM input must contain distinct amplitude coefficient definitions")
    files = {}
    for k, (pos, name) in enumerate(heads):
        end = heads[k + 1][0] if k + 1 < len(heads) else len(text)
        body = text[pos:end].split("=", 1)[1]
        body = body[: body.rindex(";")]
        body = MREF.sub(r"M__\1__\2__\3", body)
        files[f"{name}.txt"] = body
    opt = (inputs / "vvpmasters_opt-1.0.m").read_text()
    rules = {}
    for m in re.finditer(r'M\["(\w+)", (\d+), (\d+)\] -> (.*?)\n,\n', opt + "\n,\n", re.S):
        rules[f"{m.group(1)},{m.group(2)},{m.group(3)}"] = m.group(4)
    files["rules.json"] = json.dumps(rules)
    if metadata["inputs"] != cache_provenance.file_hashes(form_inputs(inputs)):
        raise cache_provenance.CacheMismatch("FORM inputs changed while the split cache was being built")
    cache.parent.mkdir(parents=True, exist_ok=True)
    # Commit a complete new directory. A parse/write failure leaves the old cache intact.
    with tempfile.TemporaryDirectory(prefix=cache.name + ".", dir=cache.parent) as temporary:
        staged = pathlib.Path(temporary) / "new"
        staged.mkdir()
        for name, contents in files.items():
            (staged / name).write_text(contents)
        record = {"provenance": metadata,
                  "outputs": cache_provenance.file_hashes({name: staged / name for name in files})}
        cache_provenance.atomic_write(staged / "manifest.json", record)
        previous = pathlib.Path(temporary) / "previous"
        if cache.exists():
            cache.rename(previous)
        try:
            staged.rename(cache)
        except BaseException:
            if previous.exists():
                previous.rename(cache)
            raise
    return metadata


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
            v = li22(mp.re(args[2]), mp.re(args[3]))
        elif name == "Zeta":
            v = mp.zeta(int(mp.re(args[0])))
        else:
            raise ValueError(name)
        if not mp.isfinite(v):raise ValueError('Nonfinite FORM master function')
        self.cache[key] = v
        return v


def form_functions(name, args):
    if name == "Den":
        return 1 / mp.mpmathify(args[0])
    if name == "Num":
        return args[0]
    if name == "Log":
        return mp.log(args[0])
    if name == "PolyLog":
        return mp.polylog(int(mp.re(args[0])), args[1])
    raise ValueError(name)


def evaluate(inputs, outdir, beta, k, dps, parts, *, recompute=False):
    inputs, outdir = pathlib.Path(inputs), pathlib.Path(outdir)
    if not parts or len(parts) != len(set(parts)) or set(parts) - {"1l", "2l_CA", "2l_CF"}:
        raise ValueError("parts must be a nonempty, distinct subset of 1l,2l_CA,2l_CF")
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"zzB_b{beta.replace('/', 'o')}_d{k}_p{dps}.json"
    metadata = form_provenance(inputs, "zz-form-amplitude",
        {"beta": beta, "delta_exp": k, "dps": dps, "parts": parts,
         "ma2": 1, "mb2": 1, "scheme": "qT", "mu2": "s"})
    validator = lambda record: validate_point(record, route_b=True)
    if cache_provenance.load_result(path, metadata, recompute=recompute, validator=validator) is not None:
        return path, "cached"
    cache = outdir / "_split"
    split_metadata = split_inputs(inputs, cache, recompute=recompute)
    if split_metadata["inputs"] != metadata["inputs"]:
        raise cache_provenance.CacheMismatch("FORM inputs changed before point evaluation")
    rules = json.loads((cache / "rules.json").read_text())
    mp.mp.dps = dps
    kin = kinematics(beta, k)
    t0 = time.time()
    g = GGVV.__new__(GGVV)
    u, x, y, z = GGVV.variables(g, kin["s"], kin["t"], mp.mpf(1), mp.mpf(1))
    L = letters(x, y, z)
    root = mp.sqrt((kin["s"] - 2) ** 2 - 4)
    mf = MasterFunctions()
    msym = dict(L, I=mp.mpc(0, 1), Pi=+mp.pi)
    masters = {}
    A = {p: {} for p in parts}
    for j in range(1, 21):
        for p in parts:
            loop = "1l" if p == "1l" else "2l"
            body = (cache / f"A{j}_A_{loop}_qT.txt").read_text()
            syms = dict(L, t=kin["t"], u=u, ma2=mp.mpf(1), mb2=mp.mpf(1), root=root,
                        CA=mp.mpf(1 if p == "2l_CA" else 0), CF=mp.mpf(1 if p == "2l_CF" else 0),
                        Pi=+mp.pi, i_=mp.mpc(0, 1))
            for name in set(re.findall(r"M__\w+?__\d+__\d+", body)):
                if name not in masters:
                    fam, idx, w = name[3:].rsplit("__", 2)
                    masters[name] = Evaluator(msym, mf).evaluate(rules[f"{FAMILY[fam]},{idx},{w}"])
                syms[name] = masters[name]
            A[p][j] = Evaluator(syms, form_functions).evaluate(body)
    H = helicity_amplitudes({p: A[p] for p in parts} | {q: {j: mp.mpc(0) for j in range(1, 21)}
                                                        for q in ("1l", "2l_CF", "2l_CA") if q not in parts}, kin)
    c = lambda zz: [mp.nstr(mp.re(zz), dps), mp.nstr(mp.im(zz), dps)]
    rec = {"beta": beta, "delta_exp": k, "dps": dps, "route": "B", "parts": parts,
           "pT2": mp.nstr(kin["pT2"], dps), "G": c(kin["G"]),
           "M": {h: {p: c(d[p]) for p in parts} for h, d in H.items()},
           "A": {p: {str(j): c(v) for j, v in A[p].items()} for p in parts},
           "sec": round(time.time() - t0, 1), "provenance": metadata}
    cache_provenance.validate_result(rec, validator, path)
    cache_provenance.atomic_write(path, rec)
    return path, rec["sec"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs")
    parser.add_argument("outdir")
    parser.add_argument("beta")
    parser.add_argument("k", type=int)
    parser.add_argument("dps", type=int)
    parser.add_argument("parts")
    parser.add_argument("--recompute", action="store_true",
                        help="rebuild split inputs and replace the evaluated output")
    args = parser.parse_args()
    path, sec = evaluate(args.inputs, args.outdir, args.beta, args.k, args.dps,
                         args.parts.split(","), recompute=args.recompute)
    print(path, sec)


if __name__ == "__main__":
    main()
