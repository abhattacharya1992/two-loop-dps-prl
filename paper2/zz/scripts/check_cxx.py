#!/usr/bin/env python3
"""Check of the Python translation of ggvvamp (route A) against the original C++ code.

Usage: check_cxx.py BUILD_GGV SRC WORKDIR POINT...
       BUILD_GGV  the translation (build/ggv), SRC  the original code (build/ggvvamp-1.0),
       WORKDIR    where the C++ files are copied and compiled (build/cxx),
       POINT      'bench' (the benchmark point of evalcmdlinegg.cpp) or BETA:K (a point of the
                  forward line, delta = 10^-K, equal masses m = 1, as in zz_points.py)

The generated files p0003.cpp ... p2464.cpp, ggvvamp.xh and ggvvamp.xi are copied unchanged,
except that p2464.cpp takes C_F and C_A from the harness instead of 4/3 and 3, so that the two
two-loop parts come out separately. They are compiled with g++ in the quad precision of the
original code (libquadmath), with scripts/cxx/ggvvamp.h in place of the original header, which
needs GiNaC and CLN. The polylogarithms (p0002.cpp, through GiNaC) are not compiled: their
values come from the Python translation at 50 digits. The check therefore compares the
translation of the algebraic part and of the assignment, which is where the C_F and C_A terms
differ, with the compiled original, for the nine computed coefficients A_j, at one loop and
separately for the C_F and C_A terms at two loops.

Copied sources are refreshed by content. Objects and the executable are reused only with a
matching build manifest covering sources, headers, compiler, flags and output hashes.
Stale or untagged build products are rebuilt; a failed build preserves the previous products.
"""
import concurrent.futures
import argparse
import hashlib
import json
import os
import pathlib
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import mpmath as mp
from ggvv import GGVV, Functions, RSIZE
from zz_points import kinematics
from zz_cache import cache_provenance
from numeric_checks import require_close

HERE = pathlib.Path(__file__).resolve().parent
CXX = os.environ.get("CXX", "g++")


def write_if_changed(path, contents):
    """Refresh a copied/generated source without trusting its timestamp."""
    if path.is_file() and cache_provenance.sha256_file(path) == hashlib.sha256(contents).hexdigest():
        return
    handle, temporary = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    try:
        with os.fdopen(handle, "wb") as stream:
            stream.write(contents)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def prepare(src, work):
    src, work = pathlib.Path(src).resolve(), pathlib.Path(work).resolve()
    s = work / "src"
    s.mkdir(parents=True, exist_ok=True)
    files = [f"p{k:04d}.cpp" for k in range(3, 2465)]
    for name in files[:-1] + ["ggvvamp.xh", "ggvvamp.xi"]:
        write_if_changed(s / name, (src / name).read_bytes())
    text = (src / "p2464.cpp").read_text()
    patched = text.replace("T CF = n(4,3);", "T CF = harness_CF;").replace("T CA = n(3);", "T CA = harness_CA;")
    assert patched.count("harness_") == 2, "p2464.cpp: Casimir lines not found"
    write_if_changed(s / "p2464.cpp", patched.encode())
    write_if_changed(s / "main.cpp", (HERE / "cxx" / "main.cpp").read_bytes())
    calls = "".join(f"  g.p{k:04d}();\n" for k in range(3, 2464))
    write_if_changed(s / "calls.inc", calls.encode())
    return s, files + ["main.cpp"]


def compiler_identity():
    command = shlex.split(CXX)
    if not command or not shutil.which(command[0]):
        raise FileNotFoundError(f"C++ compiler not found: {CXX}")
    command[0] = str(pathlib.Path(shutil.which(command[0])).resolve())
    inputs = {"compiler/driver": pathlib.Path(command[0])}
    # CXX may contain a launcher such as 'ccache g++'. Hash executable arguments too.
    for index, argument in enumerate(command[1:], 1):
        executable = shutil.which(argument)
        if executable:
            inputs[f"compiler/command-{index}"] = pathlib.Path(executable).resolve()

    def query(option):
        return subprocess.run([*command, option], check=True, capture_output=True, text=True).stdout.strip()

    details = {"command": command, "version": query("--version"), "target": query("-dumpmachine")}
    for name in ("cc1plus", "collect2", "ld"):
        reported = query("-print-prog-name=" + name)
        resolved = pathlib.Path(reported) if pathlib.Path(reported).is_file() else None
        if resolved is None and shutil.which(reported):
            resolved = pathlib.Path(shutil.which(reported))
        details[name] = reported
        if resolved is not None:
            inputs[f"compiler/{name}"] = resolved.resolve()
    library = query("-print-file-name=libquadmath.so")
    if not pathlib.Path(library).is_file():
        raise FileNotFoundError(f"Compiler did not locate libquadmath.so: {library}")
    inputs["compiler/libquadmath.so"] = pathlib.Path(library).resolve()
    return command, details, inputs


def read_dependencies(path, cwd):
    text = path.read_text().replace("\\\n", " ")
    if ":" not in text:
        raise ValueError(f"Missing compiler dependency target in {path}")
    # Preserve include symlinks: retargeting one must change the bytes hashed next time.
    return {os.path.abspath(cwd / name) for name in shlex.split(text.split(":", 1)[1])}


def compile_all(s, files, *, flags=None):
    s = pathlib.Path(s).resolve()
    flags = list(flags) if flags is not None else ["-O1", "-std=gnu++11", "-w", "-I.", f"-I{HERE / 'cxx'}"]
    link_flags = ["-lquadmath"]
    compiler, details, compiler_inputs = compiler_identity()
    inputs = {f"source/{name}": s / name for name in files}
    inputs.update({f"header/{name}": s / name for name in ("ggvvamp.xh", "ggvvamp.xi", "calls.inc")})
    inputs.update({f"harness/{p.relative_to(HERE / 'cxx')}": p for p in (HERE / "cxx").rglob("*.h")})
    inputs.update(compiler_inputs)
    environment = {key: os.environ.get(key) for key in (
        "CPATH", "CPLUS_INCLUDE_PATH", "C_INCLUDE_PATH", "LIBRARY_PATH", "GCC_EXEC_PREFIX",
        "COMPILER_PATH", "LD_LIBRARY_PATH", "LD_PRELOAD", "SOURCE_DATE_EPOCH", "LANG", "LC_ALL")}

    def metadata(dependencies):
        all_inputs = dict(inputs)
        all_inputs.update({f"dependency/{name}": pathlib.Path(name) for name in dependencies})
        return cache_provenance.fingerprint("zz-cxx-build", inputs=all_inputs,
            sources={"check_cxx.py": pathlib.Path(__file__), "dps/provenance.py": cache_provenance.__file__},
            parameters={"compiler": details, "compile_flags": flags, "dependency_flags": ["-MD", "-MF"],
                        "link_flags": link_flags, "environment": environment})

    manifest = s / "build_manifest.json"
    outputs = [name[:-4] + ".o" for name in files] + ["ggvv_cxx"]
    base_metadata = metadata([])
    try:
        old = json.loads(manifest.read_text())
        dependencies = old["dependencies"]
        if not isinstance(dependencies, list) or not all(isinstance(p, str) for p in dependencies):
            raise ValueError("invalid dependency list")
        expected = metadata(dependencies)
        record = cache_provenance.load_result(manifest, expected)
        if set(record["outputs"]) != set(outputs):
            raise ValueError("incomplete build manifest")
        actual = cache_provenance.file_hashes({name: s / name for name in outputs})
        if actual != record["outputs"]:
            raise ValueError("changed or incomplete compiled output")
        return s / "ggvv_cxx"
    except (OSError, ValueError, KeyError, TypeError, cache_provenance.CacheMismatch):
        # Build products are disposable. Never bless or reuse legacy objects on this path.
        pass

    with tempfile.TemporaryDirectory(prefix="cxx-build.", dir=s.parent) as temporary:
        staged = pathlib.Path(temporary)

        def one(name):
            obj = staged / (name[:-4] + ".o")
            depfile = staged / (name[:-4] + ".d")
            subprocess.run([*compiler, *flags, "-MD", "-MF", str(depfile),
                            "-c", name, "-o", str(obj)], cwd=s, check=True)
            if not obj.is_file() or not obj.stat().st_size:
                raise RuntimeError(f"Compiler produced no object for {name}")
            return read_dependencies(depfile, s)

        with concurrent.futures.ThreadPoolExecutor(os.cpu_count() or 1) as ex:
            dependencies = sorted(set().union(*ex.map(one, files)))
        exe = staged / "ggvv_cxx"
        subprocess.run([*compiler, "-o", str(exe),
                        *(str(staged / (name[:-4] + ".o")) for name in files), *link_flags], cwd=s, check=True)
        if not exe.is_file() or not exe.stat().st_size:
            raise RuntimeError("Linker produced no executable")
        if metadata([]) != base_metadata:
            raise cache_provenance.CacheMismatch("C++ build inputs changed during compilation; retry the build")
        record = {"provenance": metadata(dependencies), "dependencies": dependencies,
                  "outputs": cache_provenance.file_hashes({name: staged / name for name in outputs})}
        for name in outputs:
            os.replace(staged / name, s / name)
        cache_provenance.atomic_write(manifest, record)
    return s / "ggvv_cxx"


def point(spec):
    if spec == "bench":
        return "benchmark point", (mp.mpf("15625"), mp.mpf("-2384.89727896833"), mp.mpf("6467.215561"), mp.mpf("625"))
    beta, k = spec.split(":")
    kin = kinematics(beta, int(k))
    return f"beta = {beta}, delta = 10^-{k}", (kin["s"], kin["t"], mp.mpf(1), mp.mpf(1))


def compare_output(py, output, *, rtol, atol=0):
    """Require all nine coefficients and three colour/loop parts exactly once."""
    parts={'1l','2l_CF','2l_CA'}
    if len(py)!=9 or any(set(v)!=parts for v in py.values()):
        raise ValueError('Incomplete Python coefficient set')
    expected={(j,p) for j in py for p in parts};seen=set();worst={}
    for line in output.splitlines():
        if not line.strip():continue
        j,part,re_,im_=line.split();key=(int(j),part)
        if key not in expected or key in seen:
            raise ValueError(f'Unexpected or duplicate compiled coefficient {key}')
        seen.add(key)
        value=mp.mpc(mp.mpf(re_),mp.mpf(im_))
        deviation=require_close(value,py[key[0]][part],rtol=rtol,atol=atol,label=f'A{j} {part}')
        worst[part]=max(worst.get(part,0),float(deviation))
    if seen!=expected:raise ValueError(f'Missing compiled coefficients: {sorted(expected-seen)}')
    return worst


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build',type=pathlib.Path);parser.add_argument('src',type=pathlib.Path)
    parser.add_argument('work',type=pathlib.Path);parser.add_argument('points',nargs='+')
    parser.add_argument('--rtol',help='override relative tolerance (default: 1e-25 at bench, 1e-14 forward)')
    parser.add_argument('--atol',default='0',help='absolute tolerance, including zero references')
    args=parser.parse_args()
    build,src,work=args.build,args.src,args.work
    s, files = prepare(src, work)
    exe = compile_all(s, files)
    g = GGVV(build)
    worst_all = {}
    for spec in args.points:
        mp.mp.dps = 50
        label, (S, T, ma2, mb2) = point(spec)
        u, x, y, z = g.variables(S, T, ma2, mb2)
        Q = [mp.mpf(int(k.split("/")[0])) / int(k.split("/")[1]) for _, k in g.rationals]
        r = [None] * (RSIZE + 1)
        r[1], r[2], r[3], r[4], r[5], r[6], r[7] = ma2, T / ma2, u / ma2, mb2 / ma2, x, y, z
        F = Functions()
        g.p0002(r, Q, F)
        rfile = s / f"r_{spec.replace('/', 'o').replace(':', '_')}.txt"
        rfile.write_text("".join(f"{i} {mp.nstr(v, 45)}\n" for i, v in enumerate(r) if v is not None and i > 0))
        out = subprocess.run([str(exe), str(rfile)], check=True, capture_output=True, text=True).stdout
        py = g.compute(S, T, ma2, mb2, functions=F)       # the polylogarithms come from the cache of F
        tolerance=args.rtol or ('1e-25' if spec=='bench' else '1e-14')
        worst=compare_output(py,out,rtol=tolerance,atol=args.atol)
        print(f"{label}: largest relative difference between the Python translation and the compiled "
              f"C++ code: " + ", ".join(f"{p} {mp.nstr(v, 3)}" for p, v in worst.items()), flush=True)
        for p, v in worst.items():
            worst_all[p] = max(worst_all.get(p, 0), v)
    print("largest relative difference over all points: " + ", ".join(f"{p} {mp.nstr(v, 3)}" for p, v in worst_all.items()))


if __name__ == "__main__":
    main()
