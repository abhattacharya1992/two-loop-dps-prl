#!/usr/bin/env python3
"""Validate retained blind-input provenance against Git blob hashes."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def blob_hash(data):
    return hashlib.sha1(b'blob ' + str(len(data)).encode('ascii') + b'\0' + data).hexdigest()

def main():
    manifest = json.loads((ROOT / 'checks/blind_points_manifest.json').read_text())
    assert manifest['hash_algorithm'] == 'git-blob-sha1'
    assert len(manifest['files']) == 46
    hh = zz = 0
    for path, sha in manifest['files'].items():
        file = ROOT / path
        assert file.is_file(), f'Absent: {path}'
        assert blob_hash(file.read_bytes()) == sha, f'Unexpected changes: {path}'
        hh += '/points_hh/' in path
        zz += '/points_zz/' in path
    assert (hh, zz) == (32, 9)
    print('PASS: 32 HH and 9 ZZ blind points, and 5 source/results files, match upstream Git blobs')
    print('SCOPE: no new calculation from original external amplitudes')

if __name__ == '__main__':
    main()
