"""Read complete Laurent data for downstream physics calculations.

Legacy stored records remain usable when their diagram inventory and numerical
payload are complete. This reader does not assign producer provenance to them.
Partial checkpoints are for resuming assembly, not for forming a physical kernel.
"""
from fractions import Fraction
import json
import math

from cache_utils import diagram_names


STRUCTURES = {
    'gqq': ('ek', 'ke'),
    'gg': ('S1', 'S2', 'S3', 'T4'),
    'n4': ('S1', 'S2', 'S3', 'T4'),
    'fdh': ('S1', 'S2', 'S3', 'T4'),
}
COLOURS = {'gqq': ('CF', 'CA'), 'gg': ('CA', 'nf'), 'fdh': ('CA', 'nf'), 'n4': ('ns',)}
ORDERS = {f'{i},{j}' for i in (-1, 0) for j in (-2, -1, 0)}


def _keys(value, expected, label):
    if not isinstance(value, dict) or set(value) != set(expected):
        raise ValueError(f'{label} must contain exactly {sorted(expected)}')


def _pair(value, label):
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise ValueError(f'{label} must be a real/imaginary pair')
    if any(isinstance(v, bool) or not math.isfinite(float(v)) for v in value):
        raise ValueError(f'{label} must contain finite numbers')


def load_laurent(path, kind, *, expected_x=None):
    """Return validated data without modifying the file or inventing provenance."""
    with open(path) as stream:
        record = json.load(stream)
    try:
        structures = STRUCTURES[kind]
        expected = diagram_names(kind)
        _keys(record['diagrams'], expected, 'diagram inventory')
        if ('complete' in record or 'provenance' in record) and record.get('complete') is not True:
            raise ValueError('calculation is marked incomplete or has no completion status')
        if record.get('failed_diagrams'):
            raise ValueError('calculation contains failed diagrams')
        if 'expected_diagrams' in record and sorted(record['expected_diagrams']) != expected:
            raise ValueError('recorded expected diagram inventory is inconsistent')
        x = Fraction(str(record['x']))
        if expected_x is not None and x != Fraction(str(expected_x)):
            raise ValueError(f'recorded x={record["x"]} does not match requested x={expected_x}')
        if kind != 'gqq':
            _keys(record['tree'], structures, 'tree structures')
            for structure, value in record['tree'].items():
                _pair(value, f'tree {structure}')
            if complex(*record['tree']['S1']) == 0:
                raise ValueError('tree normalisation S1 must be nonzero')
        for name, diagram in record['diagrams'].items():
            colour_field = 'pref' if kind == 'gqq' else 'colour'
            _keys(diagram[colour_field], COLOURS[kind], f'{name} colour coefficients')
            for colour, value in diagram[colour_field].items():
                if kind == 'gqq':
                    if not math.isfinite(float(Fraction(str(value)))):
                        raise ValueError(f'{name} {colour} must be finite')
                else:
                    _pair(value, f'{name} {colour}')
            _keys(diagram['laurent'], structures, f'{name} helicity structures')
            for structure, coefficients in diagram['laurent'].items():
                # Assembly represents an identically zero structure by {}.
                if coefficients:
                    _keys(coefficients, ORDERS, f'{name} {structure} Laurent orders')
                elif not isinstance(coefficients, dict):
                    raise ValueError(f'{name} {structure} has no Laurent coefficient mapping')
                for order, value in coefficients.items():
                    _pair(value, f'{name} {structure} coefficient {order}')
    except (KeyError, TypeError, ValueError, OverflowError, ZeroDivisionError) as exc:
        raise ValueError(f'{path}: invalid or incomplete {kind} Laurent data: {exc}') from exc
    return record
