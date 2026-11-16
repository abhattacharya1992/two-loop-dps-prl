"""Finite-value and tolerance gates for region quadrature diagnostics."""
import numpy as np
import sympy as sp
import scipy.special as ss
import mpmath as mp


def close(actual, expected, *, atol=1e-9, rtol=1e-4, label='quadrature'):
    if not (np.all(np.isfinite(actual)) and np.all(np.isfinite(expected))):
        raise ValueError(f'{label}: nonfinite comparison')
    if np.any(np.abs(actual-expected) > atol + rtol*np.maximum(np.abs(actual),np.abs(expected))):
        raise ValueError(f'{label}: {actual!r} vs {expected!r} exceeds atol={atol:g}, rtol={rtol:g}')


def finite_array(value):
    value=np.asarray(value,dtype=complex)
    if not np.all(np.isfinite(value)):
        raise ValueError('Nonfinite integrand at an interior quadrature node')
    return value


def interior_values(expr, variables, coordinates):
    """Factor before numerical evaluation to avoid cancelling endpoint factors.

    If floating-point evaluation still loses a small denominator, recompute
    those nodes at higher precision. Never replace nonfinite values by zero.
    """
    expression=sp.factor(expr,deep=True)
    function=sp.lambdify(variables,expression,[{'gamma':ss.gamma},'numpy'])
    grids=np.broadcast_arrays(*coordinates)
    with np.errstate(all='ignore'):
        values=np.asarray(function(*coordinates),dtype=complex)*np.ones_like(grids[0],dtype=complex)
    missing=~np.isfinite(values)
    if np.any(missing):
        accurate=sp.lambdify(variables,expression,'mpmath')
        with mp.workdps(70):
            for index in zip(*np.nonzero(missing)):
                arguments=[mp.mpc(complex(grid[index])) for grid in grids]
                values[index]=complex(accurate(*arguments))
    return finite_array(values)
