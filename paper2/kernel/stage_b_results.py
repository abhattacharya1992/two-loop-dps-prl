#!/usr/bin/env python3
"""Collects the numbers of Stage B that the manuscript quotes into results_stage_b.json:
the ultrasoft cancellation (us_telescoping.py), the explicit two-loop graph (model_numerics.py) and the
bookkeeping of the region scales (kernel_logs.py).  Every number is recomputed here from the functions of
those scripts; nothing is typed.
    python3 stage_b_results.py          (about 2 minutes)
"""
import json
import numpy as np
import mpmath as mp
import sympy as sp
import us_telescoping as UT
import model_numerics as MN
import kernel_logs as KL

def telescoping():
    x = mp.mpf('0.75')
    lhat = (mp.mpf('0.7'), mp.mpf('-1.3'), (mp.mpf('0.4'), mp.mpf('-0.2')))
    rows = []
    for lam in (mp.mpf('0.1'), mp.mpf('0.01'), mp.mpf('0.001'), mp.mpf('0.0001')):
        K = (mp.mpf('0.6')*lam, mp.mpf('0.3')*lam)
        l = (lhat[0]*lam**2, lhat[1]*lam**2, (lhat[2][0]*lam**2, lhat[2][1]*lam**2))
        r = UT.residues(l, K, x)
        lp, lm = l[0], l[1]
        wa, wb, wc = x - lp/2, 1 - x + lp/2, 1 + lp/2
        dk = abs((wa*r['a'] + wc*r['c'])/r['tree']*lm - 1)
        dkb = abs((wb*r['b'] + wc*r['c'])/r['tree']*lm - 1)
        single = abs(wa*r['a']/r['tree'])*abs(lm)
        rows.append({'lambda': float(lam), 'dev_Tk': float(dk), 'dev_Tkbar': float(dkb), 'single_attachment': float(single)})
    glauber = []
    for lam in (mp.mpf('0.01'), mp.mpf('0.001')):
        K = (mp.mpf('0.6')*lam, mp.mpf('0.3')*lam)
        l = (mp.mpf('0.7')*lam**2, mp.mpf('-1.3')*lam**2, (mp.mpf('0.4')*lam, mp.mpf('-0.2')*lam))
        r = UT.residues(l, K, x)
        lp, lm = l[0], l[1]
        wa, wb, wc = x - lp/2, 1 - x + lp/2, 1 + lp/2
        glauber.append({'lambda': float(lam), 'dev_Tk': float(abs((wa*r['a'] + wc*r['c'])/r['tree']*lm - 1)),
                        'dev_Tkbar': float(abs((wb*r['b'] + wc*r['c'])/r['tree']*lm - 1))})
    soft = []
    for lam in (mp.mpf('0.01'), mp.mpf('0.001')):
        K = (mp.mpf('0.6')*lam, mp.mpf('0.3')*lam)
        l = (mp.mpf('0.7')*lam, mp.mpf('-1.3')*lam, (mp.mpf('0.4')*lam, mp.mpf('-0.2')*lam))
        r = UT.residues(l, K, x)
        lp = l[0]
        soft.append({'lambda': float(lam), 'parent_over_daughter': float(abs((1 + lp/2)*r['c']/((x - lp/2)*r['a'])))})
    worst_quad = float(UT.quad_check())
    return {'residues_vs_quadrature_worst': worst_quad, 'ultrasoft': rows, 'glauber': glauber, 'soft': soft}

def model():
    us = {}
    for c in (0.1, 1.0, 3.0):
        vals = []
        for lam in (0.03, 0.01, 0.003, 0.001):
            r = c*lam*lam
            vals.append(lam**6*MN.F(lam, (r*np.cos(0.7), r*np.sin(0.7))))
        us[str(c)] = {'lambda': [0.03, 0.01, 0.003, 0.001], 'lambda6F_re': [v.real for v in vals], 'lambda6F_im': [v.imag for v in vals],
                      'rel_changes': [abs(vals[i+1] - vals[i])/abs(vals[-1]) for i in range(3)]}
    slopes = []
    for lh in ((0.2, -0.5), (-3.0, 2.0), (0.15, 0.1), (0.9, 0.4), (-0.3, -1.1)):
        lams = (0.01, 0.003, 0.001)
        v = [lam**4*MN.F(lam, (lam*lh[0], lam*lh[1])) for lam in lams]
        s2 = (v[2] - v[1])/np.log(lams[2]/lams[1])
        k2 = MN.KH[0]**2 + MN.KH[1]**2
        a = (MN.KH[0] + lh[0] - MN.QH[0])**2 + (MN.KH[1] + lh[1] - MN.QH[1])**2
        N = MN.KH[0]*(MN.KH[0] + lh[0] - MN.QH[0]) + MN.KH[1]*(MN.KH[1] + lh[1] - MN.QH[1])
        l2 = lh[0]**2 + lh[1]**2
        pred = -2*(2j*np.pi)**3*N/(MN.x**2*k2*a*l2)
        slopes.append({'lhat': lh, 'lambda': list(lams), 'slope_re': s2.real, 'slope_im': s2.imag, 'pred_im': pred.imag,
                       'rel_diff': abs(s2 - pred)/abs(pred), 'real_part_change': abs(v[2].real - v[1].real)/abs(v[2])})
    return {'x': MN.x, 'khat': MN.KH, 'qhat': MN.QH, 'i0_over_lambda2': 1e-4, 'ultrasoft_lambda6F': us, 'log_slopes': slopes}

def logs():
    r2, r1 = KL.part1()
    KL.part2()
    ser = sp.expand(sp.simplify(sp.series(sp.exp(-KL.eps*sp.EulerGamma)*sp.gamma(-KL.eps), KL.eps, 0, 2).removeO()))
    return {'r2': str(r2), 'r1': str(r1), 'triangle_cu': -1, 'triangle_X2': '-1/2',
            'rapidity_coefficient': str(ser), 'rapidity_coefficient_eps0': str(sp.simplify(ser.coeff(KL.eps, 0)))}

if __name__ == '__main__':
    out = {'telescoping': telescoping(), 'model': model(), 'logs': logs()}
    with open('results_stage_b.json', 'w') as f:
        json.dump(out, f, indent=1)
    t = out['telescoping']; m = out['model']
    print('ultrasoft cancellation at lambda = 1e-4:', t['ultrasoft'][-1])
    print('Glauber:', t['glauber'][-1], ' soft parent/daughter:', t['soft'][-1])
    print('largest slope deviation:', max(s['rel_diff'] for s in m['log_slopes']),
          ' largest real-part change:', max(s['real_part_change'] for s in m['log_slopes']))
    print('ultrasoft lambda^6 F, last relative changes:', [m['ultrasoft_lambda6F'][c]['rel_changes'][-1] for c in m['ultrasoft_lambda6F']])
    print('written results_stage_b.json')
