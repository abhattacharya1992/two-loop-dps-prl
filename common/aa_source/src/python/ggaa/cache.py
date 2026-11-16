"""Content provenance and atomic writes for resumable numerical evaluations."""
import hashlib
import json
import os
from pathlib import Path
import platform
import tempfile

import mpmath as mp


def fingerprint(paths, **parameters):
    hashes = {}
    for name, path in paths.items():
        with Path(path).open('rb') as stream:
            digest = hashlib.sha256()
            for block in iter(lambda: stream.read(1024 * 1024), b''):
                digest.update(block)
            hashes[name] = digest.hexdigest()
    runtime={'python':platform.python_version(),'implementation':platform.python_implementation(),
             'machine':platform.machine(),'mpmath':mp.__version__,'backend':mp.libmp.BACKEND}
    if mp.libmp.BACKEND=='gmpy':
        import gmpy2
        runtime['gmpy2']=gmpy2.version()
    return {'schema':1,'files':hashes,'parameters':parameters,'runtime':runtime}


def load(path, provenance):
    if path is None or not Path(path).exists():
        return {}
    try:
        record = json.loads(Path(path).read_text())
    except (ValueError, OSError):
        return {}
    if (not isinstance(record, dict) or record.get('provenance') != provenance or
            not isinstance(record.get('values'), dict)):
        return {}
    return record['values']


def write(path, provenance, values):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=path.name + '.', dir=path.parent)
    try:
        with os.fdopen(handle, 'w') as stream:
            json.dump({'provenance': provenance, 'values': values}, stream, indent=1)
            stream.write('\n')
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
