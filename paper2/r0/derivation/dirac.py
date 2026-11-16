"""Gamma matrices in even dimension d, metric (+,-,...,-), with the light cone vectors n = (1,0,...,0,1) and
nbar = (1,0,...,0,-1): vslash = (v+/2) nslash + (v-/2) nbarslash - v_perp.Gamma, v+ = v0 + v(d-1), v- = v0 - v(d-1),
and Euclidean transverse vectors of d-2 components. P = nslash nbarslash/4 projects onto the good components."""
import numpy as np

_S = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]


class Space:
    def __init__(self, d):
        m = d//2
        def kron(ops):
            out = np.eye(1, dtype=complex)
            for o in ops: out = np.kron(out, o)
            return out
        E = []                                   # Hermitian matrices with square 1
        for j in range(m):
            for a in (0, 1):
                E.append(kron([_S[2]]*j + [_S[a]] + [np.eye(2)]*(m - j - 1)))
        self.d = d; self.n = d - 2
        self.g0 = E[0]; self.gs = [1j*e for e in E[1:]]          # gamma^0 and the d-1 spatial gamma matrices
        self.gt = self.gs[:-1]                                    # transverse
        self.nsl = self.g0 - self.gs[-1]; self.nbsl = self.g0 + self.gs[-1]
        self.P = self.nsl @ self.nbsl/4; self.Pb = self.nbsl @ self.nsl/4
        self.one = np.eye(2**m)

    def perp(self, v): return -sum(vi*g for vi, g in zip(v, self.gt))
    def slash(self, plus, minus, v): return plus/2*self.nsl + minus/2*self.nbsl + self.perp(v)
    def contract(self, X): return self.g0 @ X @ self.g0 - sum(g @ X @ g for g in self.gs)     # gamma^a X gamma_a
    def basis(self, r):
        """the 2(d-2) points +-r e_i: their mean reproduces the angular average of a polynomial of degree three or less"""
        return [sg*r*np.eye(self.n)[i] for i in range(self.n) for sg in (1, -1)]


def charge_conjugation(sp):
    """C with C gamma^mu C^-1 = -(gamma^mu)^T, four dimensions"""
    assert sp.d == 4
    gam = [sp.g0] + sp.gs
    for trial in ([1, 3], [0, 2]):               # product of the gamma matrices that are symmetric (or antisymmetric)
        C = np.eye(4, dtype=complex)
        for i in trial: C = C @ gam[i]
        if all(abs(C @ g @ np.linalg.inv(C) + g.T).max() < 1e-12 for g in gam): return C
    raise ValueError('no charge conjugation matrix found')


dot = lambda a, b: float(np.dot(a, b))
mx = lambda M: float(abs(M).max())
