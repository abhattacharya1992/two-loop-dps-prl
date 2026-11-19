#!/usr/bin/env python3
"""Verify or recover archived sector-decomposition performance caches.

These files accelerate costly original region integrations. The included
Laurent coefficients and all normal publication tests do not require them.
Recover only files from the verified project Git history.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "checks" / "region_cache_manifest.json"


def git(*args):
    result = subprocess.run(["git", "-C", str(ROOT), *args],
                            capture_output=True, check=True)
    return result.stdout


def inspect_record(entry, revision):
    path = entry["path"]
    blob = git("rev-parse", f"{revision}:{path}").decode().strip()
    if blob != entry["git_sha"]:
        raise ValueError(f"{path}: Git object differs from the archived manifest")
    raw = git("cat-file", "blob", blob)
    if len(raw) != entry["size"] or hashlib.sha256(raw).hexdigest() != entry["sha256"]:
        raise ValueError(f"{path}: archived data size or SHA-256 mismatch")
    return raw


def verify(manifest):
    revision = manifest["baseline"]
    for entry in manifest["entries"]:
        inspect_record(entry, revision)
    return len(manifest["entries"])


def restore(manifest, *, force=False):
    recovered = []
    for entry in manifest["entries"]:
        path = ROOT / entry["path"]
        if path.exists() and not force:
            existing = path.read_bytes()
            if hashlib.sha256(existing).hexdigest() != entry["sha256"]:
                raise ValueError(f"{path}: existing file differs; use --force to replace")
            continue
        raw = inspect_record(entry, manifest["baseline"])
        path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=path.parent, prefix=".region-cache-",
                                         delete=False) as stream:
            temp = Path(stream.name)
            try:
                stream.write(raw)
                stream.flush()
            except BaseException:
                temp.unlink(missing_ok=True)
                raise
        temp.replace(path)
        recovered.append(entry["path"])
    return recovered


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true",
                        help="verify each archived Git blob against the recorded SHA-256")
    parser.add_argument("--restore", action="store_true",
                        help="restore optional performance caches from the frozen Git commit")
    parser.add_argument("--force", action="store_true",
                        help="replace a conflicting local cache (with --restore)")
    args = parser.parse_args()
    if args.force and not args.restore:
        parser.error("--force requires --restore")
    data = json.loads(MANIFEST.read_text())
    if not (args.verify or args.restore):
        args.verify = True
    if args.verify:
        print(f"Verified {verify(data)} archived sector-decomposition caches")
    if args.restore:
        restored = restore(data, force=args.force)
        print(f"Restored {len(restored)} optional caches (others already match)")


if __name__ == "__main__":
    main()
