"""Evaluate the FORM (gg_Aj_A_qT) and Wolfram Language (vvpmasters_opt) expressions of
VVamp with mpmath, as an independent route to the coefficients A_j^[A].

The evaluator is an operator-precedence (shunting-yard) pass over the text, as in the
Paper 1 package (ggaa.twoloop): only a stack of numbers is kept. Function calls are written
f(...) in FORM and f[...] in the Wolfram Language; both are accepted.
"""
import re

import mpmath as mp

_TOKEN = re.compile(r"\s*(?:(\d+)|([A-Za-z][A-Za-z0-9_]*)([\(\[])?|(.))", re.S)
_PREC = {"+": 1, "-": 1, "*": 2, "/": 2, "neg": 3, "^": 4}


def _apply(op, v):
    if op == "neg":
        v.append(-v.pop())
        return
    b, a = v.pop(), v.pop()
    if op == "+":
        v.append(a + b)
    elif op == "-":
        v.append(a - b)
    elif op == "*":
        v.append(a * b)
    elif op == "/":
        v.append(mp.mpf(a) / b if isinstance(a, int) and isinstance(b, int) else a / b)
    else:
        if isinstance(b, int) or (mp.im(b) == 0 and mp.re(b) == int(mp.re(b))):
            n = int(mp.re(b)) if not isinstance(b, int) else b
            v.append((mp.mpf(a) if isinstance(a, int) else a) ** n)
        else:
            v.append(a ** b)


class Evaluator:
    def __init__(self, symbols, functions):
        self.symbols = symbols
        self.functions = functions

    def evaluate(self, text):
        values, ops, arity = [], [], []
        expect = True
        for m in _TOKEN.finditer(text):
            number, name, call, other = m.groups()
            if number is not None:
                values.append(int(number))
                expect = False
            elif name is not None and call:
                ops.append(("call", name, "]" if call == "[" else ")"))
                arity.append(1)
                expect = True
            elif name is not None:
                values.append(self.symbols[name])
                expect = False
            elif other is None:
                break
            elif other in "([":
                ops.append(other)
                expect = True
            elif other in ")]":
                # a closing bracket ends either a group or a call
                while not (ops[-1] in ("(", "[") or isinstance(ops[-1], tuple)):
                    _apply(ops.pop(), values)
                top = ops.pop()
                if isinstance(top, tuple):
                    n = arity.pop()
                    args = values[-n:]
                    del values[-n:]
                    values.append(self.functions(top[1], args))
                expect = False
            elif other == ",":
                while not isinstance(ops[-1], tuple):
                    _apply(ops.pop(), values)
                arity[-1] += 1
                expect = True
            elif other in "+-*/^":
                if expect:
                    if other == "-":
                        ops.append("neg")
                    continue
                p = _PREC[other]
                while ops and isinstance(ops[-1], str) and ops[-1] in _PREC:
                    q = _PREC[ops[-1]]
                    if q > p or (q == p and other != "^"):
                        _apply(ops.pop(), values)
                    else:
                        break
                ops.append(other)
                expect = True
            elif other.isspace() or other == ";":
                continue
            else:
                raise ValueError(f"unexpected character {other!r}")
        while ops:
            _apply(ops.pop(), values)
        if len(values) != 1:
            raise ValueError("malformed expression")
        return values[0]


def letters(x, y, z):
    """The 19 letters of eqs. (5.11), (5.13) of arXiv:1503.04812."""
    return {
        "l1": x, "l2": 1 + x, "l3": y, "l4": 1 - y, "l5": z, "l6": 1 - z, "l7": -y + z,
        "l8": 1 + y - z, "l9": 1 + x * y, "l10": 1 + x * z, "l11": x * y + z,
        "l12": 1 + y + x * y - z, "l13": 1 + x + x * y - x * z,
        "l14": 1 + y + 2 * x * y - z + x**2 * y * z,
        "l15": 2 * x * y + x**2 * y + x**2 * y**2 + z - x**2 * y * z,
        "l16": 1 + x + y + x * y + x * y**2 - z - x * z - x * y * z,
        "l17": 1 + y + x * y + y**2 + x * y**2 - z - y * z - x * y * z,
        "l18": -x * y + z + x * z + x * y * z, "l19": -y + z + y * z + x * y * z}
