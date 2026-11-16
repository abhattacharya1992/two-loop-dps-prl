"""Provenance shared by the ZZ point drivers and their translated-input loader.

Legacy numerical records remain analysis inputs. Only records with matching
actual input bytes, evaluator sources, parameters and runtime can resume a run.
"""
import pathlib
import sys
from itertools import product

HERE = pathlib.Path(__file__).resolve().parent
for parent in HERE.parents:
    for package in (parent / "common/dps", parent / "dps"):
        if (package / "dps/provenance.py").is_file():
            sys.path.insert(0, str(package))
            break
    else:
        continue
    break
from dps import provenance as cache_provenance

TRANSLATED = ("p0002.py", "algebra.marshal", "rationals.json", "coefficients.json")
FORM_INPUTS = ("gg_Aj_A_qT-1.0.inc", "vvpmasters_opt-1.0.m")


def translated_inputs(build):
    return {name: pathlib.Path(build) / name for name in TRANSLATED}


def form_inputs(inputs):
    return {name: pathlib.Path(inputs) / name for name in FORM_INPUTS}


def source_files(*names):
    files = {name: HERE / name for name in ("zz_cache.py", *names)}
    files["dps/provenance.py"] = pathlib.Path(cache_provenance.__file__)
    return files


def translated_provenance(g, stage, parameters, *, wide=False):
    names = ["ggvv.py", "translate_ggvvamp.py", "zz_points.py"]
    if wide:
        names.append("zz_wide.py")
    metadata = cache_provenance.fingerprint(
        stage, inputs=translated_inputs(g.build), sources=source_files(*names),
        parameters=parameters, libraries=("mpmath",))
    if metadata["inputs"] != g.input_hashes:
        raise cache_provenance.CacheMismatch(
            "Translated ZZ inputs changed after GGVV was loaded; construct a new GGVV "
            "instance before evaluating or resuming points.")
    return metadata


def form_provenance(inputs, stage, parameters):
    return cache_provenance.fingerprint(
        stage, inputs=form_inputs(inputs),
        sources=source_files("zz_routeB.py", "formeval.py", "ggvv.py", "zz_points.py",
                             "zz_master_families.py"),
        parameters=parameters, libraries=("mpmath",))


def validate_point(record, *, wide=False, route_b=False):
    """Require the complete finite payload, not merely a matching metadata tag."""
    required = {"beta", "dps", "pT2", "G", "M", "sec", "provenance"}
    required.update({"j", "sin2theta", "delta", "x", "s", "t"} if wide else {"delta_exp", "A"})
    if route_b:
        required.update({"route", "parts"})
    else:
        required.update({"li22_max_rate", "li22_bad_order", "x", "s", "t"})
    cache_provenance.require_keys(record, required)
    parameters = record["provenance"]["parameters"]
    for name in ("beta", "dps", "j" if wide else "delta_exp"):
        if record[name] != parameters[name]:
            raise ValueError(f"{name} differs from the recorded parameters")
    parts = parameters["parts"]
    if route_b and (record["route"] != "B" or record["parts"] != parts):
        raise ValueError("route or parts differs from the recorded parameters")
    for name in ("pT2", "sec"):
        cache_provenance.require_finite(record[name])
    for name in ("x", "s", "t", "delta", "sin2theta", "li22_max_rate"):
        if name in required:
            cache_provenance.require_finite(record[name])
    if not route_b and (type(record["li22_bad_order"]) is not int or record["li22_bad_order"] < 0):
        raise ValueError("li22_bad_order must be a nonnegative integer")
    cache_provenance.require_pair(record["G"])
    helicities = ["".join(h) for h in product("+-", "+-", "+0-", "+0-")]
    cache_provenance.require_keys(record["M"], helicities, exact=True)
    for amplitudes in record["M"].values():
        cache_provenance.require_keys(amplitudes, parts, exact=True)
        for value in amplitudes.values():
            cache_provenance.require_pair(value)
    if not wide:
        cache_provenance.require_keys(record["A"], parts, exact=True)
        for coefficients in record["A"].values():
            cache_provenance.require_keys(coefficients, [str(j) for j in range(1, 21)], exact=True)
            for value in coefficients.values():
                cache_provenance.require_pair(value)
