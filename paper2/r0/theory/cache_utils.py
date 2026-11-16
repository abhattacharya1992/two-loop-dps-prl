"""Provenance and checkpoint handling shared by the region calculations."""
import argparse
import os
from pathlib import Path
import pickle
import sys
import warnings

HERE = Path(__file__).resolve().parent
for parent in HERE.parents:
    for package in (parent / 'common/dps', parent / 'dps'):
        if (package / 'dps/provenance.py').is_file():
            sys.path.insert(0, str(package))
            break
    else:
        continue
    break
from dps import provenance

CacheMismatch = provenance.CacheMismatch
atomic_write = provenance.atomic_write

INTEGRAND_SOURCES = {
    'gqq': ('run_gqq.py', 'driver.py', 'gqq_diagrams.py', 'lc_engine.py'),
    'gg': ('run_gg.py', 'gg_driver.py', 'gg_diagrams.py', 'vec_engine.py', 'lc_engine.py'),
    'n4': ('run_n4.py', 'gg_driver.py', 'gg_diagrams.py', 'n4_diagrams.py', 'vec_engine.py', 'lc_engine.py'),
    'fdh': ('run_n4.py', 'gg_driver.py', 'gg_diagrams.py', 'n4_diagrams.py', 'vec_engine.py', 'lc_engine.py'),
}
ASSEMBLY_SOURCES = ('assemble_lib.py', 'assemble_gg.py', 'assemble_gqq.py', 'assemble_n4.py',
                    'secdec.py', 'factors.py', 'lc_engine.py',
                    'gg_driver.py', 'vec_engine.py', 'gg_diagrams.py', 'n4_diagrams.py', 'gqq_diagrams.py')


def sources(names):
    return {**{name: HERE / name for name in names}, 'cache_utils.py': Path(__file__),
            'dps/provenance.py': Path(provenance.__file__)}


def arguments(description, *, fdh=False, assembly=False):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('x', nargs='?', default='3/10')
    parser.add_argument('names', nargs='?', help='comma-separated diagram names')
    parser.add_argument('--recompute', action='store_true',
                        help='replace the existing output; discard all old diagram results')
    if fdh:
        parser.add_argument('--fdh', action='store_true')
    if assembly:
        parser.add_argument('--allow-legacy-input', action='store_true',
                            help='use supplied untagged integrands, recorded as legacy/unverified')
    args = parser.parse_args()
    args.only = args.names.split(',') if args.names else None
    return args


def diagram_names(kind):
    if kind == 'gqq':
        from gqq_diagrams import DIAGRAMS
    elif kind == 'n4':
        from n4_diagrams import DIAGRAMS
    else:
        from gg_diagrams import DIAGRAMS
        if kind == 'fdh':
            DIAGRAMS = {name: diag for name, diag in DIAGRAMS.items() if name.startswith('V')}
    return sorted(DIAGRAMS)


def integrand_provenance(kind, x):
    return provenance.fingerprint(f'{kind}-integrands', sources=sources(INTEGRAND_SOURCES[kind]),
        parameters={'x': str(x), 'scheme': 'FDH' if kind == 'fdh' else 'CDR',
                    'expected_diagrams': diagram_names(kind), 'arithmetic': 'symbolic'},
        libraries=('sympy', 'numpy', 'mpmath'))


def new_integrands(path, kind, x, *, recompute=False):
    metadata = integrand_provenance(kind, x)
    record = provenance.load_result(path, metadata, recompute=recompute, binary=True,
                                    validator=lambda r: validate_integrands(r, kind))
    return record if record is not None else {'x': x, 'results': {}, 'provenance': metadata,
        'expected_diagrams': diagram_names(kind), 'complete': False, 'failed_diagrams': []}


def checkpoint(path, record, key):
    record['complete'] = (set(record[key]) == set(record['expected_diagrams'])
                          and not record.get('failed_diagrams'))
    if key == 'results':
        validate_integrands(record, record['provenance']['stage'].split('-')[0])
    else:
        validate_laurent(record)
    atomic_write(path, record, binary=(key == 'results'))


def selected_names(only, names):
    unknown = set(only or ()) - set(names)
    if unknown:
        raise ValueError(f'unknown diagram names: {sorted(unknown)}')
    return list(only) if only is not None else list(names)


def load_integrands(path, kind, x, *, allow_legacy=False):
    with open(path, 'rb') as stream:
        data = pickle.load(stream)
    expected = diagram_names(kind)
    provenance.validate_result(data, lambda r: validate_integrands(r, kind, legacy='provenance' not in r), path)
    if str(data.get('x')) != str(x) or set(data.get('results', {})) != set(expected):
        raise CacheMismatch(f'{path}: wrong x or incomplete diagram set; regenerate all integrands first.')
    if 'provenance' not in data:
        if not allow_legacy:
            raise CacheMismatch(f'{path}: untagged legacy integrands. Regenerate them, or explicitly '
                                'use --allow-legacy-input for a calculation from supplied data.')
        warnings.warn(f'{path}: using legacy integrands without verified producer provenance', stacklevel=2)
        data['_input_provenance_status'] = 'legacy-unverified'
    else:
        provenance.require_match(data, integrand_provenance(kind, x), path)
        if not data.get('complete') or data.get('failed_diagrams'):
            raise CacheMismatch(f'{path}: incomplete or failed integrand calculation; resume it first.')
        data['_input_provenance_status'] = 'verified'
    return data


def quadrature():
    nqs = os.environ.get('NQ')
    neps = int(os.environ.get('NEPS', '32'))
    counts = [64, 48, 28] if not nqs else list(map(int, nqs.split(',')))
    if len(counts) != 3 or min(counts) < 1 or neps < 1:
        raise ValueError('NQ must contain three positive orders; NEPS must be positive')
    suffix = ('' if not nqs else f'_nq{nqs.replace(",", "-")}') + ('' if neps == 32 else f'_neps{neps}')
    return dict(zip((1, 2, 3), counts)), neps, suffix


def laurent_provenance(data, input_path, script, nq, neps):
    import mpmath as mp
    return provenance.fingerprint('region-laurent', inputs={Path(input_path).name: input_path},
        sources=sources(ASSEMBLY_SOURCES + (script,)),
        parameters={'x': str(data['x']), 'nq': nq, 'neps': neps, 'mpmath_dps': mp.mp.dps,
                    'eps_radius': 0.2, 'eta_radius_fraction': 0.05, 'n_eta': 8,
                    'eta_orders': [-1, 0], 'eps_orders': [-2, -1, 0],
                    'input_provenance_status': data['_input_provenance_status'],
                    'expected_diagrams': sorted(data['results'])},
        libraries=('sympy', 'numpy', 'scipy', 'mpmath'))


def new_laurent(path, data, metadata, *, recompute=False):
    out = provenance.load_result(path, metadata, recompute=recompute, validator=validate_laurent)
    return out if out is not None else {'x': str(data['x']), 'diagrams': {},
        'provenance': metadata, 'expected_diagrams': sorted(data['results']), 'complete': False}


def validate_integrands(record, kind, *, legacy=False):
    import sympy as sp
    expected = diagram_names(kind)
    if not isinstance(record['results'], dict) or not set(record['results']) <= set(expected):
        raise ValueError('missing results mapping or unknown diagrams')
    if not legacy:
        if str(record['x']) != record['provenance']['parameters']['x']:
            raise ValueError('x differs from provenance')
        if sorted(record['expected_diagrams']) != expected:
            raise ValueError('wrong expected diagram set')
        if (not isinstance(record['complete'], bool) or not isinstance(record['failed_diagrams'], list)
                or not set(record['failed_diagrams']) <= set(expected)
                or set(record['failed_diagrams']) & set(record['results'])):
            raise ValueError('invalid completion or failure status')
        if bool(record['complete']) != (set(record['results']) == set(expected) and not record['failed_diagrams']):
            raise ValueError('inconsistent completion status')
    structures = ('ek', 'ke') if kind == 'gqq' else ('S1', 'S2', 'S3', 'T4')
    if kind != 'gqq':
        provenance.require_keys(record['tree'], structures, exact=True)
        if any(sp.sympify(value).is_finite is not True for value in record['tree'].values()):
            raise ValueError('nonfinite tree normalisation')
    for pieces in record['results'].values():
        if not isinstance(pieces, (list, tuple)):
            raise ValueError('integrand pieces must be a sequence')
        for interval, parts in pieces:
            if len(interval) != 2 or not all(sp.sympify(v).is_finite for v in interval):
                raise ValueError('invalid integration interval')
            for parameters, values in parts:
                if not isinstance(parameters, (list, tuple)):
                    raise ValueError('invalid integration parameters')
                provenance.require_keys(values, structures, exact=True)
                if any(sp.sympify(v).has(sp.nan, sp.oo, -sp.oo, sp.zoo) for v in values.values()):
                    raise ValueError('nonfinite integrand')


def validate_laurent(record):
    import sympy as sp
    expected = record['provenance']['parameters']['expected_diagrams']
    if sorted(record['expected_diagrams']) != sorted(expected) or not isinstance(record['diagrams'], dict):
        raise ValueError('invalid diagram registry')
    if not set(record['diagrams']) <= set(expected):
        raise ValueError('unknown cached diagram')
    if not isinstance(record['complete'], bool) or record['complete'] != (set(record['diagrams']) == set(expected)):
        raise ValueError('inconsistent completion status')
    if record['x'] != record['provenance']['parameters']['x']:
        raise ValueError('x differs from provenance')
    is_gqq = any(name.startswith('gqq_integrands_') for name in record['provenance']['inputs'])
    is_n4 = any(name.startswith('n4_integrands_') for name in record['provenance']['inputs'])
    colour_keys = ('CF', 'CA') if is_gqq else ('ns',) if is_n4 else ('CA', 'nf')
    structures = ('ek', 'ke') if is_gqq else ('S1', 'S2', 'S3', 'T4')
    allowed_orders = {f'{i},{j}' for i in (-1, 0) for j in (-2, -1, 0)}
    for diagram in record['diagrams'].values():
        provenance.require_keys(diagram, ('pref' if is_gqq else 'colour', 'laurent'))
        provenance.require_keys(diagram['laurent'], structures, exact=True)
        weights = diagram['pref' if is_gqq else 'colour']
        provenance.require_keys(weights, colour_keys, exact=True)
        if is_gqq and any(sp.sympify(value).is_finite is not True for value in weights.values()):
            raise ValueError('nonfinite colour coefficients')
        if not is_gqq:
            for value in diagram['colour'].values():
                provenance.require_pair(value)
        for coefficients in diagram['laurent'].values():
            if coefficients:
                provenance.require_keys(coefficients, allowed_orders, exact=True)
            elif not isinstance(coefficients, dict):
                raise ValueError('missing Laurent coefficient mapping')
            for value in coefficients.values():
                provenance.require_pair(value)
    if not is_gqq:
        provenance.require_keys(record['tree'], structures, exact=True)
        for value in record['tree'].values():
            provenance.require_pair(value)


def sector_cache(path):
    """Only source-compatible sector decompositions are reusable.

    This is a performance cache, not a numerical result. Legacy or incompatible
    entries are ignored and rebuilt, with a warning. Symbolic expressions and
    integration parameters are additionally part of each entry's key.
    """
    metadata = provenance.fingerprint('sector-decomposition', sources=sources(ASSEMBLY_SOURCES),
                                      libraries=('sympy',))
    try:
        record = provenance.load_result(path, metadata, binary=True)
    except CacheMismatch as exc:
        warnings.warn(f'{exc} Rebuilding the sector-decomposition performance cache.', stacklevel=2)
        record = None
    return metadata, {} if record is None else record['entries']
