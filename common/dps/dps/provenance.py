"""Content-based provenance for resumable scientific calculations.

Stored datasets may be analysed without this metadata. Resuming a calculation,
however, requires the same actual inputs, source files, parameters and runtime.
An untagged or incompatible result is preserved until recomputation is requested.
"""
import hashlib
from decimal import Decimal, InvalidOperation
import importlib.metadata
import importlib
import json
import os
import pathlib
import pickle
import platform
import tempfile


class CacheMismatch(RuntimeError):
    """A result cannot safely be reused by the requested calculation."""


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def file_hashes(paths):
    """Hash actual bytes; keys are stable labels, independent of checkout paths."""
    return {str(name): sha256_file(path) for name, path in sorted(paths.items())}


def fingerprint(stage, *, inputs=None, sources=None, parameters=None, libraries=()):
    runtime = {'python': platform.python_version(),
               'implementation': platform.python_implementation(),
               'machine': platform.machine()}
    runtime.update({name: importlib.metadata.version(name) for name in libraries})
    if 'mpmath' in libraries:
        runtime['mpmath_backend'] = importlib.import_module('mpmath.libmp').BACKEND
        if runtime['mpmath_backend'] == 'gmpy':
            runtime['gmpy2'] = importlib.metadata.version('gmpy2')
    metadata = {'schema': 1, 'stage': stage, 'inputs': file_hashes(inputs or {}),
                'sources': file_hashes(sources or {}), 'parameters': parameters or {},
                'runtime': runtime}
    # Round-trip also normalises tuples and integer dictionary keys as JSON does.
    metadata = json.loads(json.dumps(metadata, sort_keys=True, allow_nan=False))
    metadata['fingerprint'] = hashlib.sha256(
        json.dumps(metadata, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    return metadata


def require_match(record, expected, path):
    stored = record.get('provenance') if isinstance(record, dict) else None
    if stored != expected:
        reason = 'untagged legacy result' if not stored else 'changed inputs, code, parameters or runtime'
        raise CacheMismatch(f'{path}: {reason}; refusing to resume. Preserve the old result and '
                            'use a fresh output directory, or explicitly pass --recompute.')


def require_finite(value):
    try:
        valid = not isinstance(value, bool) and Decimal(str(value)).is_finite()
    except (InvalidOperation, ValueError):
        valid = False
    if not valid:
        raise ValueError(f'not a finite real number: {value!r}')


def require_pair(value):
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise ValueError('complex value must be a real/imaginary pair')
    for component in value:
        require_finite(component)


def require_keys(value, required, *, exact=False):
    if not isinstance(value, dict) or (set(value) != set(required) if exact else not set(required) <= set(value)):
        raise ValueError(f'missing or unexpected fields; expected {sorted(required)}')


def validate_result(record, validator, path):
    try:
        validator(record)
    except (KeyError, TypeError, ValueError, ArithmeticError, AssertionError) as exc:
        raise CacheMismatch(f'{path}: malformed or incomplete scientific result ({exc}); '
                            'preserve it and use a fresh directory, or pass --recompute.') from exc


def load_result(path, expected, *, recompute=False, binary=False, validator=None):
    path = pathlib.Path(path)
    if recompute or not path.exists():
        return None
    try:
        with path.open('rb' if binary else 'r') as stream:
            record = pickle.load(stream) if binary else json.load(stream)
    except (ValueError, EOFError, pickle.UnpicklingError) as exc:
        raise CacheMismatch(f'{path}: unreadable cache; preserve it or pass --recompute.') from exc
    require_match(record, expected, path)
    if validator is not None:
        validate_result(record, validator, path)
    return record


def atomic_write(path, record, *, binary=False):
    """A failed calculation/write must not leave a truncated checkpoint."""
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=path.name + '.', dir=path.parent)
    try:
        with os.fdopen(handle, 'wb' if binary else 'w') as stream:
            if binary:
                pickle.dump(record, stream, protocol=4)
            else:
                json.dump(record, stream, indent=1, allow_nan=False)
                stream.write('\n')
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
