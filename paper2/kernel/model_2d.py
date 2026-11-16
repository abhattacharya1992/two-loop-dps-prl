#!/usr/bin/env python3
"""Stage B2, step 1d: the transverse integral of the cross-pair model,
    K_cross(lam) = (1/4) int d^2 l_perp F(l_perp),   F(l_perp) = int dl^+ dl^- (...) from model_fast.py,
with k = lam khat, q = lam qhat, eps0 = eps0_rel lam^2.

Singular points of F in the transverse plane (probed numerically, 16 Oct 2026):
  P0 = 0 (soft gluon): F = c0 + c1 ln|l_perp| (integrable; large values, |F| ~ 1e12 in |l_perp| < 0.05 lam);
  P1 = q - k (t = l_perp + k - q = 0, the D4 line without transverse momentum): F ~ (k.t/|t|^2) G + regular,
       the leading term is odd in t and integrates to zero over the angle;
  elsewhere smooth (after the correction of the side of the z6' pole; with the old side, F grew like
  1/distance from the half line l_perp = s (q - k), s > 1, a fake pinch).  At large |l_perp|, F ~ 1/|l_perp|^5.

Integration: polar coordinates about P0 everywhere except a 'polar rectangle' R1 = {r in [d(1-rho), d(1+rho)],
|phi - phi0| <= alpha} around P1 (d = |q - k|), which is done in polar coordinates about P1 with the ray
length to the boundary of R1.  Radii are log-spaced (geometric segments), angles split at the directions
of the corners of R1 (about P1) and at phi0 +- alpha (about P0); Gauss-Legendre on every segment.
The disc |l_perp| < r_in about P0 and the disc |t| < r1_in about P1 are dropped (their contributions
scale like r^2 ln r and r^2).  Points are evaluated in parallel (multiprocessing).

    python3 model_2d.py            # convergence tests at lam = 0.03 (several minutes)
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from multiprocessing import Pool
import model_fast as MF
PI = np.pi

def gl_segments(edges, n):
    """Gauss-Legendre nodes and weights on consecutive segments [edges[i], edges[i+1]]"""
    gx, gw = leggauss(n)
    a = np.asarray(edges[:-1]); b = np.asarray(edges[1:])
    mid = 0.5*(a+b); half = 0.5*(b-a)
    return (mid[:, None] + half[:, None]*gx[None, :]).ravel(), (half[:, None]*gw[None, :]).ravel()

def log_edges(r0, r1, per_decade):
    """geometric breakpoints from r0 to r1 with about per_decade segments per decade (at least one segment)"""
    n = max(1, int(np.ceil(per_decade*np.log10(r1/r0))))
    return np.geomspace(r0, r1, n+1)

def ray_length(px, py, ux, uy, r1, r2, phi0, alpha):
    """distance from P1 = (px, py) along the unit vector u to the boundary of the polar rectangle
    {r in [r1, r2], |phi - phi0| <= alpha} (about the origin), P1 inside."""
    hits = []
    pu = px*ux + py*uy; d2 = px*px + py*py
    def ang_ok(X, Y):
        dphi = (np.arctan2(Y, X) - phi0 + PI) % (2*PI) - PI
        return abs(dphi) <= alpha + 1e-12
    # outer arc r = r2
    disc = pu*pu + r2*r2 - d2
    if disc >= 0:
        s = -pu + np.sqrt(disc)
        if s > 0 and ang_ok(px+s*ux, py+s*uy): hits.append(s)
    # inner arc r = r1
    disc = pu*pu + r1*r1 - d2
    if disc >= 0:
        for s in (-pu - np.sqrt(disc), -pu + np.sqrt(disc)):
            if s > 0 and ang_ok(px+s*ux, py+s*uy): hits.append(s)
    # radial sides phi = phi0 +- alpha
    for sgn in (1, -1):
        vx, vy = np.cos(phi0+sgn*alpha), np.sin(phi0+sgn*alpha)
        cr = ux*vy - uy*vx
        if abs(cr) > 1e-14:
            s = -(px*vy - py*vx)/cr
            if s > 0:
                t = (px+s*ux)*vx + (py+s*uy)*vy
                if r1 - 1e-12 <= t <= r2 + 1e-12: hits.append(s)
    assert hits, 'ray does not leave the polar rectangle'
    return min(hits)

def nodes_2d(lam, khat, qhat, rho=0.3, alpha=0.3, r_in_rel=1e-4, r1_in_rel=1e-3, rmax_rel=300.0,
             per_decade=2, n_r=8, n_phi=12, n_phi_big=40):
    """all (x, y, weight) nodes of the 2D quadrature; weight includes the Jacobian r"""
    kx, ky = lam*khat[0], lam*khat[1]; qx, qy = lam*qhat[0], lam*qhat[1]
    px, py = qx - kx, qy - ky; d = np.hypot(px, py); phi0 = np.arctan2(py, px)
    r1, r2 = d*(1-rho), d*(1+rho)
    X, Y, W, TAG = [], [], [], []
    # ---- about P0: angular segments [phi0-alpha, phi0+alpha] (n_phi) and the rest (3 segments, n_phi_big/3 each)
    phi_edges_small = [phi0-alpha, phi0+alpha]
    phi_edges_big = [phi0+alpha, phi0+PI/2, phi0+3*PI/2, phi0+2*PI-alpha]
    ph_s, wph_s = gl_segments(phi_edges_small, n_phi)
    ph_b, wph_b = gl_segments(phi_edges_big, max(4, n_phi_big//3))
    def add_polar0(redges, which, tag):
        rr, wr = gl_segments(redges, n_r)
        for ph, wph in which:
            R, P = np.meshgrid(rr, ph, indexing='ij'); WR, WP = np.meshgrid(wr, wph, indexing='ij')
            X.append((R*np.cos(P)).ravel()); Y.append((R*np.sin(P)).ravel()); W.append((WR*WP*R).ravel())
            TAG.append(np.full(R.size, tag))
    both = [(ph_s, wph_s), (ph_b, wph_b)]
    add_polar0(log_edges(r_in_rel*lam, r1, per_decade), both, 0)          # inner region (contains the soft point)
    add_polar0(log_edges(r1, r2, per_decade), [(ph_b, wph_b)], 1)        # shell, outside the window of R1
    add_polar0(log_edges(r2, rmax_rel*lam, per_decade), both, 2)         # outer region
    # ---- R1 about P1
    corners = [np.arctan2(r*np.sin(phi0+s*alpha) - py, r*np.cos(phi0+s*alpha) - px) for r in (r1, r2) for s in (1, -1)]
    corners = sorted(c % (2*PI) for c in corners)
    sect_edges = corners + [corners[0] + 2*PI]
    gx, gw = leggauss(n_phi)
    for i in range(4):
        a, b = sect_edges[i], sect_edges[i+1]
        for g, w in zip(gx, gw):
            ph = 0.5*(a+b) + 0.5*(b-a)*g; wph = 0.5*(b-a)*w
            ux, uy = np.cos(ph), np.sin(ph)
            rmax = ray_length(px, py, ux, uy, r1, r2, phi0, alpha)
            rr, wr = gl_segments(log_edges(r1_in_rel*lam, rmax, per_decade), n_r)
            X.append(px + rr*ux); Y.append(py + rr*uy); W.append(wr*wph*rr); TAG.append(np.full(rr.size, 3))
    return np.concatenate(X), np.concatenate(Y), np.concatenate(W), np.concatenate(TAG)

_ctx = {}
def _init(kx, ky, qx, qy, x, eps0):
    _ctx.update(kx=kx, ky=ky, qx=qx, qy=qy, x=x, eps0=eps0)

def _F(pt):
    lx, ly = pt
    c = _ctx
    L = max(1000.0, 100.0*(lx*lx + ly*ly))
    return MF.lplus_integral(lx, ly, c['kx'], c['ky'], c['qx'], c['qy'], c['x'], c['eps0'], L=L)

def kernel_cross(lam, khat=(0.6, 0.3), qhat=(1.0, 0.0), x=0.75, eps0_rel=1e-4, nproc=2, return_parts=False, **grid):
    kx, ky = lam*khat[0], lam*khat[1]; qx, qy = lam*qhat[0], lam*qhat[1]
    X, Y, W, TAG = nodes_2d(lam, khat, qhat, **grid)
    with Pool(nproc, initializer=_init, initargs=(kx, ky, qx, qy, x, eps0_rel*lam*lam)) as pool:
        vals = np.array(pool.map(_F, list(zip(X, Y)), chunksize=64))
    tot = 0.25*np.sum(W*vals)
    if return_parts:
        parts = {t: 0.25*np.sum((W*vals)[TAG == t]) for t in range(4)}
        return tot, parts, len(X)
    return tot

if __name__ == '__main__':
    import time, sys
    lam = 0.03
    print(f'lam = {lam}, khat = (0.6, 0.3), qhat = (1, 0), x = 0.75, eps0 = 1e-4 lam^2')
    tests = [dict(),                                     # default grid
             dict(n_r=12, n_phi=18, n_phi_big=60),       # more nodes
             dict(per_decade=3),                         # finer radial segments
             dict(rho=0.2, alpha=0.2),                   # different partition
             dict(r_in_rel=1e-5, r1_in_rel=1e-4, rmax_rel=1000.0)]   # cut-offs
    for g in tests:
        t0 = time.time()
        v, parts, npts = kernel_cross(lam, return_parts=True, **g)
        print(f'{str(g):60s} K_cross = {v:.6e}  [{npts} points, {time.time()-t0:.0f} s]')
        print('    parts: inner', f'{parts[0]:.4e}', 'shell', f'{parts[1]:.4e}', 'outer', f'{parts[2]:.4e}', 'R1', f'{parts[3]:.4e}')
        sys.stdout.flush()
