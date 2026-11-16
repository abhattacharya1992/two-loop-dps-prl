"""Numerical comparison gates; unavailable or nonfinite results cannot pass."""
import mpmath as mp


def require_close(actual, reference, *, rtol, atol=0, label='comparison'):
    actual,reference=mp.mpmathify(actual),mp.mpmathify(reference)
    rtol,atol=mp.mpf(rtol),mp.mpf(atol)
    if not all(mp.isfinite(v) for v in (actual,reference,rtol,atol)) or min(rtol,atol)<0:
        raise ValueError(f'{label}: nonfinite value or invalid tolerance')
    error=abs(actual-reference)
    if error>atol+rtol*abs(reference):
        raise ValueError(f'{label}: difference {mp.nstr(error,6)} exceeds rtol={rtol}, atol={atol}')
    return error/abs(reference) if reference else error
