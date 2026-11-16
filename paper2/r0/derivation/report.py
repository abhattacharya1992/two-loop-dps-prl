"""Bookkeeping for the checks. Each check prints one line; a script ends with status 1 if any check fails."""
import sys

_failed = []


def check(label, deviation, tolerance):
    """pass if |deviation| <= tolerance"""
    ok = abs(deviation) <= tolerance
    print(f"  [{'ok' if ok else 'FAIL'}] {label}: {abs(deviation):.1e} (tolerance {tolerance:.0e})")
    if not ok: _failed.append(label)


def require(label, condition):
    print(f"  [{'ok' if condition else 'FAIL'}] {label}")
    if not condition: _failed.append(label)


def finish():
    if _failed:
        print(f'{len(_failed)} check(s) failed: ' + '; '.join(_failed))
        sys.exit(1)
    print('All checks passed.')
