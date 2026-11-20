#!/usr/bin/env python3
# Figures 2-5 of the paper, drawn from the result files listed in SOURCES of make_tables.py,
# so that the figures and the tables always read the same files. Figure 1, the sketch of the
# DPS configuration, is TikZ source (figures/dps_sketch.tex), compiled by `make figure`.
#
#   Fig. 2  figures/forward_approach.pdf  exact one-loop gg -> AA and gg -> ZZ amplitudes from 90
#                                         degrees to the forward limit, in units of their DPS
#                                         terms (Sec. 2.3)
#   Fig. 3  figures/hh_coefficients.pdf   two-loop gg -> HH coefficients and the predictions (Sec. 3.3)
#   Fig. 4  figures/zz_coefficients.pdf   two-loop gg -> ZZ coefficients and the predictions (Sec. 6.2)
#   Fig. 5  figures/fit_windows.pdf       deviation from the prediction against the fit window and
#                                         the fit model (App. B)
#
# Every series has its own marker shape as well as its own colour, so the figures also read in
# black and white. The categorical colours are used in a fixed order (blue, orange, aqua); in
# Fig. 5 each fit model keeps its colour in all panels, and the fits without p_T^2 terms are grey.
# In Fig. 2 the real parts are solid and the imaginary parts dotted.
#
#   python3 make_figure.py
import logging
import pathlib
from fractions import Fraction
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from make_tables import load, preferred_beta, zz_predictions, lin_coeffs, ZZ_LOG, ZZ_K

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE/'figures'
mp.mp.dps = 40
# fontTools warns that the time stamps in the header of the cmr10 font file are implausibly
# early whenever it embeds the font; the warning is harmless
logging.getLogger('fontTools').setLevel(logging.ERROR)

INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
C1, C2, C3, C5 = '#2a78d6', '#eb6834', '#1baf7a', '#e87ba4'      # categorical slots 1-3 and 5
MARKS = ('o', 's', '^')
plt.rcParams.update({'font.size': 9, 'font.family': 'serif', 'font.serif': ['cmr10'],
                     'mathtext.fontset': 'cm', 'axes.formatter.use_mathtext': True,
                     'axes.unicode_minus': False, 'pdf.fonttype': 42,   # Computer Modern, embedded as TrueType
                     'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.linewidth': 0.6,
                     'legend.frameon': False})


def frame(a, title=None):
    a.grid(True, color=GRID, lw=0.5)
    a.set_axisbelow(True)
    if title:
        a.set_title(title, fontsize=9, loc='left', color=INK)


def save(fig, name):
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT/f'{name}.pdf', metadata={'CreationDate': None})    # no date, so reruns give the same file
    plt.close(fig)
    print(f'wrote figures/{name}.pdf')


def betas_of(d):
    return sorted(d, key=Fraction)


def cplx(v):
    return mp.mpc(mp.mpf(v[0]), mp.mpf(v[1]))


# ---------------------------------------------------------------- Fig. 2: from 90 degrees to the forward limit
def forward_figure():
    """The exact one-loop amplitudes along the lines of fixed beta, from 90 degrees (on the right)
    towards the forward limit, in units of their DPS terms (Sec. 2.3), against
    p_T^2/m^2. Top: the derivative of the amplitude with Delta J_z = 0 with respect to G at fixed
    beta, in units of C_1, which tends to 1 if the amplitude contains C_1 G plus terms that remain
    finite (the difference quotient between neighbouring points, at the geometric mean of their
    p_T^2). Bottom: the amplitude with |Delta J_z| = 2 in units of K, which tends to 1. Real parts
    solid (filled symbols), imaginary parts dotted (open symbols). gg -> AA from the remainders of
    arXiv:1909.08993 as evaluated by P1, on a dense grid (lines); gg -> ZZ at the evaluated points
    (symbols)."""
    sets = ((load('aa_curve'), 'AA', (r'$\mathcal{M}_1^{(1)}$', r'$\mathcal{M}_2^{(1)}$')),
            (load('zz_wide'), 'ZZ', (r'$\mathcal{M}_{{+}{+}{+}{+}}^{(1)}$', r'$\mathcal{M}_{{+}{-}{+}{+}}^{(1)}$')))
    assert betas_of(sets[0][0]) == betas_of(sets[1][0])
    fig, ax = plt.subplots(2, 2, figsize=(7.0, 5.5), sharex=True, sharey='row')
    dot = (0, (1.0, 1.6))
    colours = (C1, C2, C3)
    for col, (data, proc, names) in enumerate(sets):
        for i, b in enumerate(betas_of(data)):
            rows = data[b]['rows']                     # from large to small p_T^2
            pT2 = [mp.mpf(r['pT2']) for r in rows]
            G = [cplx(r['G']) for r in rows]
            m1 = [cplx(r['M1_over_C1_minus_G'])+g for r, g in zip(rows, G)]          # M_1/C_1
            d1 = [(m1[k+1]-m1[k])/(G[k+1]-G[k]) for k in range(len(rows)-1)]
            pm = [mp.sqrt(pT2[k]*pT2[k+1]) for k in range(len(rows)-1)]
            r2 = [cplx(r['M2_over_K']) for r in rows]
            for a, xs, ys in ((ax[0, col], pm, d1), (ax[1, col], pT2, r2)):
                x = [float(v) for v in xs]
                re, im = [float(mp.re(v)) for v in ys], [float(mp.im(v)) for v in ys]
                if proc == 'AA':
                    a.plot(x, re, color=colours[i], lw=1.5, zorder=3-0.1*i)
                    a.plot(x, im, color=colours[i], lw=1.3, ls=dot, zorder=3-0.1*i)
                else:
                    a.plot(x, re, color=colours[i], lw=0.9, marker=MARKS[i], ms=5, mfc=colours[i], mec='white',
                           mew=0.6, zorder=3)
                    a.plot(x, im, color=colours[i], lw=0.9, ls=dot, marker=MARKS[i], ms=4.6, mfc='white',
                           mec=colours[i], mew=1.0, zorder=3)
        for row, name in enumerate(names):
            what = (r'$\Delta J_z=0$', r'$|\Delta J_z|=2$')[row]
            frame(ax[row, col], f'({"abcd"[2*row+col]}) $gg\\to {proc}$, {what}: {name}')
    for a in ax.flat:
        a.set_xscale('log')
        for y in (0, 1):
            a.axhline(y, color=INK2, lw=0.8, ls=(0, (5, 3)), zorder=1)
    ax[0, 0].set_xlim(1e-8, 1.2)
    ax[0, 0].set_xticks([10.0**e for e in range(-8, 1, 2)])
    for a in ax.flat:
        a.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    for a in ax[1]:
        a.set_xlabel(r'$p_T^2/m^2$')
    ax[0, 0].set_ylabel(r'$(\partial\mathcal{M}^{(1)}/\partial\mathcal{G})\,/\,C_1$')
    ax[1, 0].set_ylabel(r'$\mathcal{M}^{(1)}/K$')
    handles = [Line2D([], [], color=c, lw=1.5, marker=mk, ms=5, mec='white', mew=0.6, label=f'$\\beta={b}$')
               for c, mk, b in zip(colours, MARKS, betas_of(sets[0][0]))]
    handles += [Line2D([], [], color=INK, lw=1.5, marker='o', ms=4.6, mfc=INK, mec='white', mew=0.6,
                       label='real part'),
                Line2D([], [], color=INK, lw=1.3, ls=dot, marker='o', ms=4.4, mfc='white', mec=INK, mew=1.0,
                       label='imaginary part'),
                Line2D([], [], color=INK2, lw=0.8, ls=(0, (5, 3)), label='DPS value')]
    fig.legend(handles=handles, loc='upper center', ncol=6, fontsize=8, handlelength=2.2,
               columnspacing=1.3, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, 'forward_approach')


# ---------------------------------------------------------------- Fig. 3: gg -> HH coefficients
def hh_figure():
    """Top: the N_c = 3, n_f = 5 coefficients at mu_R = m_H, divided by their colour and beta0
    factors, on the one-loop weight curves F_S(x) and W(x), x = (1+beta)/2. Bottom: relative
    deviation of every colour term with a nonzero prediction (task B)."""
    res = load('colour')
    BETAS = betas_of(res)
    b = np.array([float(Fraction(x)) for x in BETAS])
    N, nf = 3, 5
    weight = {'N2': N**2, 'N_Nf': N*nf, 'Nf_over_N': nf/N, 'Nf2': nf**2}
    beta0 = (11*N-2*nf)/3

    def total(beta, coeff):
        rows = [r for r in res[beta]['rows'] if r['coefficient'] == coeff]
        return sum(weight[r['colour']]*float(r['re']) for r in rows)

    def dev(beta, coeff, colour):
        return float([r for r in res[beta]['rows'] if r['coefficient'] == coeff and r['colour'] == colour][0]['deviation'])

    norm = {'spin0 G^3': 8/3*np.pi*N**2, 'spin0 G^2': 4/3*np.pi*N*(11*N-2*nf),
            'spin2 G^2': -4*np.pi*N**2, 'spin2 G^1': (4*N-4*beta0)*np.pi*N}
    w = np.linspace(0.5, 0.97, 200)
    FS = 1/(w*(1-w))+2*w*(1-w)
    Kt = (w**4+(1-w)**4)/(w*(1-w))
    bb = 2*w-1

    fig, ax = plt.subplots(2, 2, figsize=(7.0, 5.2), gridspec_kw={'height_ratios': [1.25, 1]})
    for a in ax.flat:
        frame(a)
        a.set_xlim(0.15, 0.75)

    a = ax[0, 0]
    a.plot(bb, FS, color=INK2, lw=1.2, ls='--', label=r'$F_S((1+\beta)/2)$')
    a.plot(b, [total(x, 'spin0 G^3')/norm['spin0 G^3'] for x in BETAS], 'o', ms=10, mfc='none', mec=C1,
           mew=1.6, label=r'$c_{13}/[\frac{4}{3} C_A\cdot2\pi C_A]$')
    a.plot(b, [total(x, 'spin0 G^2')/norm['spin0 G^2'] for x in BETAS], 's', ms=4.5, color=C2,
           mec=C2, label=r'$c_{12}/[2\beta_0\cdot2\pi C_A]$')
    a.set_ylim(4, 10)
    a.set_ylabel('normalised coefficient')
    a.set_title(r'(a) equal helicities, $\mathcal{M}_1$', fontsize=9, loc='left', color=INK)
    a.legend(loc='upper left', fontsize=8)

    a = ax[0, 1]
    a.plot(bb, Kt, color=INK2, lw=1.2, ls='--', label=r'$W((1+\beta)/2)$')
    a.plot(b, [total(x, 'spin2 G^2')/norm['spin2 G^2'] for x in BETAS], 'o', ms=10, mfc='none', mec=C1,
           mew=1.6, label=r'$c_{22}/[-4C_A\cdot\pi C_A]$')
    a.plot(b, [total(x, 'spin2 G^1')/norm['spin2 G^1'] for x in BETAS], 's', ms=4.5, color=C2,
           mec=C2, label=r'$c_{21}/[(4C_A-4\beta_0)\pi C_A]$')
    a.set_ylim(0, 6)
    a.set_title(r'(b) opposite helicities, $\mathcal{M}_2$', fontsize=9, loc='left', color=INK)
    a.legend(loc='upper left', fontsize=8)

    for a, coeffs, title in ((ax[1, 0], (('spin0 G^3', 'N2', 'o', C1, r'$c_{13}$, $N^2$'),
                                         ('spin0 G^2', 'N2', 's', C2, r'$c_{12}$, $N^2$'),
                                         ('spin0 G^2', 'N_Nf', '^', C3, r'$c_{12}$, $Nn_f$')),
                              r'(c) relative deviation, $\mathcal{M}_1$'),
                             (ax[1, 1], (('spin2 G^2', 'N2', 'o', C1, r'$c_{22}$, $N^2$'),
                                         ('spin2 G^1', 'N2', 's', C2, r'$c_{21}$, $N^2$'),
                                         ('spin2 G^1', 'N_Nf', '^', C3, r'$c_{21}$, $Nn_f$')),
                              r'(d) relative deviation, $\mathcal{M}_2$')):
        for i, (ck, col, mk, c, lab) in enumerate(coeffs):
            a.plot(b+(i-1)*0.012, [dev(x, ck, col) for x in BETAS], mk, ms=6, color=c, mec='white', mew=0.8, label=lab)
        a.set_yscale('log')
        a.set_ylim(1e-15, 1e-5)
        a.set_xlabel(r'$\beta$')
        a.set_title(title, fontsize=9, loc='left', color=INK)
        a.legend(loc='upper left', fontsize=8, ncol=1)
    ax[1, 0].set_ylabel(r'$|\mathrm{extracted}/\mathrm{predicted}-1|$')

    fig.tight_layout()
    save(fig, 'hh_coefficients')


# A compact comparison for the PRL Letter: predictions and independent VVamp fits.
def zz_prediction_figure():
    """Predicted leading two-loop ZZ coefficients and fits at three velocities.

    The plotted c_j/C_1 equals 2*f_j/C_1 from the input's alpha_s/(2*pi)
    convention, since the common amplitude expansion uses a_s=alpha_s/(4*pi).
    """
    data=load('zz'); pred=zz_predictions()
    betas=betas_of(data['betas'])
    grid=np.linspace(0.20,0.72,180)
    xg=(1+grid)/2
    Lx=np.log(xg*(1-xg))
    bf=np.array([float(Fraction(b)) for b in betas])
    fig,ax=plt.subplots(figsize=(3.45,2.52))
    for power,marker in ((3,'o'),(2,'s')):
        for part,colour in (('CF',C1),('CA',C2)):
            a0,a1=(float(q) for q in lin_coeffs(pred[('log',power,part)]))
            curve=a0+a1*Lx
            vals=[]
            for beta in betas:
                rec=data['betas'][beta]
                hel=ZZ_LOG[0]
                ratio=2*cplx(rec['hel'][hel]['2l_'+part]['f'][power])/cplx(rec['hel'][hel]['c1'])
                vals.append(float(mp.re(ratio)))
            style='-' if a1 else '--'
            ax.plot(grid,curve,color=colour,ls=style,lw=1.2,zorder=2)
            ax.plot(bf,vals,marker,ls='none',color=colour,ms=4.5,
                    mfc='white' if power==3 else colour,mew=1.1,zorder=3,
                    label=rf'$C_{{{part[-1]}}},\;\mathcal{{G}}^{power}$')
    frame(ax)
    ax.set_xlim(.18,.74)
    ax.set_xlabel(r'$\beta$')
    ax.set_ylabel(r'$c_j/C_1$')   # c_j = 2 f_j, converting published alpha_s/(2pi) to a_s
    ax.legend(fontsize=7,loc='lower left',ncol=2,handletextpad=.3,columnspacing=.85)
    fig.tight_layout(pad=.6)
    save(fig,'zz_prediction')



# The Letter compares the two splitting channels and all three tested hard processes.
def three_process_comparison_figure():
    """Leading and next-to-leading isotropic two-loop DPS logarithms.

    The reported HH/AA amplitudes use a_s=alpha_s/(4*pi); the published ZZ
    amplitudes instead use alpha_s/(2*pi), giving a factor two in the
    two-loop/one-loop ratios. The plotted points are *unconstrained* amplitude
    fits. Prediction curves are obtained solely from Delta_C and B; no
    plotted fit value is used in their construction.
    """
    hh, one, aa, zz = (load(k) for k in ('colour', 'one_loop', 'dps_aa', 'zz'))
    aa = aa['betas']; zz = zz['betas']
    betas = ['1/4', '1/2', '2/3']
    assert all(set(betas).issubset(d) for d in (hh, one, aa, zz))
    Nc, nf = 3, 5
    CA, CF = mp.mpf(Nc), mp.mpf(Nc*Nc-1)/(2*Nc)
    beta0 = (11*CA - 2*nf)/3
    colour_weights = {'N2': Nc*Nc, 'N_Nf': Nc*nf,
                      'Nf_over_N': mp.mpf(nf)/Nc, 'Nf2': nf*nf}
    out = {p: {j: [] for j in (3, 2)} for p in ('HH', 'AA', 'ZZ')}
    theory = {p: {j: [] for j in (3, 2)} for p in out}
    for b in betas:
        # The HH one-loop slope in results_one_loop.json is per Nc.
        ch = CA*mp.mpf(one[b]['one_loop']['DPS_minus_2piFS'])
        # The AA analytic one-loop slope in results_aa.json is also per Nc.
        ca = CA*mp.mpf(aa[b]['one_loop']['C1/N']['analytic'])
        cz = cplx(zz[b]['hel']['++++']['c1'])
        for j in (3, 2):
            h = sum(colour_weights[r['colour']]*cplx((r['re'],r['im']))
                    for r in hh[b]['rows'] if r['coefficient']==f'spin0 G^{j}')/ch
            a = sum(colour_weights[r['colour']]*cplx(r['value'])
                    for r in aa[b]['leading_logs']
                    if r['component']==1 and r['power']==j)/ca
            z = 2*sum(c*cplx(zz[b]['hel']['++++']['2l_'+label]['f'][j])
                      for label,c in (('CF',CF),('CA',CA)))/cz
            Lx = mp.log((1-(mp.mpf(Fraction(b).numerator)/Fraction(b).denominator)**2)/4)
            for proc, value, expected in (
                ('HH',h, -mp.mpf(4)/3*CA if j==3 else -2*beta0),
                ('AA',a, -mp.mpf(4)/3*CA if j==3 else -2*beta0),
                ('ZZ',z, -mp.mpf(4)/3*(2*CF-CA) if j==3 else -6*CF+4*(CA-CF)*Lx),
            ):
                if abs(value-expected)>mp.mpf('1e-8'):
                    raise AssertionError(f'{proc} beta={b}, power={j}: {value} != {expected}')
                out[proc][j].append(float(mp.re(value)))
                theory[proc][j].append(float(expected))
    bf=np.array([float(Fraction(b)) for b in betas])
    bg=np.linspace(0.18,0.73,201)
    gluon_c3=float(-mp.mpf(4)*CA/3)
    quark_c3=float(-mp.mpf(4)*(2*CF-CA)/3)
    gluon_c2=float(-2*beta0)
    quark_c2=(-6*float(CF)+4*float(CA-CF)*np.log((1-bg*bg)/4))
    fig, axes=plt.subplots(1,2,figsize=(6.75,2.6))
    for ax in axes:
        frame(ax)
        ax.set_xlim(0.18,0.73)
        ax.set_xticks([0.25,0.5,2/3]); ax.set_xticklabels([r'$1/4$',r'$1/2$',r'$2/3$'])
        ax.set_xlabel(r'$\beta$')
    axes[0].axhline(gluon_c3,lw=1.3,color=INK,ls='-',zorder=2)
    axes[0].axhline(quark_c3,lw=1.3,color=INK2,ls='--',zorder=2)
    axes[1].axhline(gluon_c2,lw=1.3,color=INK,ls='-',zorder=2)
    axes[1].plot(bg,quark_c2,lw=1.3,color=INK2,ls='--',zorder=2)
    process=(('HH',C1,'s',-0.012),('AA',C2,'D',0.012),('ZZ',C3,'o',0.0))
    for p,col,mark,dx in process:
        for j,ax in ((3,axes[0]),(2,axes[1])):
            ax.plot(bf+dx,out[p][j],linestyle='none',marker=mark,ms=5.1,
                    color=col,mfc=('white' if p=='AA' else col),mec=col,
                    mew=1.3,zorder=5,label=rf'$gg\to {p}$' if j==3 else None)
    axes[0].set_ylabel(r'$c_3/C_1$')
    axes[1].set_ylabel(r'$c_2/C_1$')
    axes[0].set_ylim(-4.65,1.1)
    axes[1].set_ylim(-22.4,-14.0)
    axes[0].set_title(r'(a) Cubic logarithm',fontsize=9,color=INK,loc='left')
    axes[1].set_title(r'(b) Quadratic logarithm',fontsize=9,color=INK,loc='left')
    processes=[Line2D([],[],marker=mk,color=c,mfc=('white' if p=='AA' else c),
                      lw=0,ms=5.5,label=rf'$gg\to {p}$') for p,c,mk,_ in process]
    prediction_lines=[Line2D([],[],color=INK,lw=1.3,label=r'$g\to gg$ prediction'),
                      Line2D([],[],color=INK2,lw=1.3,ls='--',label=r'$g\to q\bar q$ prediction')]
    fig.legend(handles=processes+prediction_lines,loc='upper center',ncol=5,
               bbox_to_anchor=(0.5,1.0),fontsize=8,handletextpad=.4,
               columnspacing=1.05,handlelength=1.6)
    fig.subplots_adjust(left=.10,right=.985,bottom=.23,top=.76,wspace=.34)
    save(fig, 'three_process_comparison')


# ---------------------------------------------------------------- Fig. 4: gg -> ZZ coefficients
def spread_labels(ys, gap):
    """positions for labels at heights ys (any order), at least `gap` apart, as close as possible
    to ys; returned in the original order"""
    order = sorted(range(len(ys)), key=lambda k: ys[k])
    pos = [ys[k] for k in order]
    for _ in range(50):
        moved = False
        for n in range(1, len(pos)):
            if pos[n]-pos[n-1] < gap:
                d = (gap-(pos[n]-pos[n-1]))/2
                pos[n-1] -= d; pos[n] += d
                moved = True
        if not moved:
            break
    out = [0.0]*len(ys)
    for k, y in zip(order, pos):
        out[k] = y
    return out


def zz_figure():
    """Top: the coefficients of eq. (6.3), 2 f_j/C_1 for ++++ and 2 f_j/K for +-++ (the factor 2
    converts the expansion in alpha_s/2pi of the source to a_s), against beta, with the
    predictions of eq. (6.3), and dotted the predictions for the coefficients that depend on L_x
    if the scale were the gluon value Q_*^2 = m^4/s, which has no terms in L_x. Bottom: absolute
    deviation from the prediction, the largest over the four amplitudes of each group (as in
    Table 7)."""
    zz = load('zz')
    pred = zz_predictions()
    betas = betas_of(zz['betas'])
    bf = np.array([float(Fraction(b)) for b in betas])
    B0, B1 = 0.15, 0.75
    grid = np.linspace(B0, B1, 200)
    xg = (1+grid)/2
    Lg = np.log(xg*(1-xg))
    groups = ((ZZ_LOG, 'log', 'c1', (3, 2), r'(a) $\Delta J_z=0$: $\mathcal{M}_{{+}{+}{+}{+}}$, units of $C_1$',
               r'(c) $\Delta J_z=0$: largest deviation'),
              (ZZ_K, 'K', 'K', (2, 1), r'(b) $|\Delta J_z|=2$: $\mathcal{M}_{{+}{-}{+}{+}}$, units of $K$',
               r'(d) $|\Delta J_z|=2$: largest deviation'))
    gname = {3: r'\mathcal{G}^3', 2: r'\mathcal{G}^2', 1: r'\mathcal{G}'}
    pname = {'CF': 'C_F', 'CA': 'C_A'}
    colour = {'CF': C1, 'CA': C2}
    fig, ax = plt.subplots(2, 2, figsize=(7.0, 5.8), gridspec_kw={'height_ratios': [1.4, 1]})
    for col, (hels, kind, key, powers, title, dtitle) in enumerate(groups):
        top, bot = ax[0, col], ax[1, col]
        labels, alts, handles, ys = [], {}, [], []
        for k, j in enumerate(powers):
            for i, part in enumerate(('CF', 'CA')):
                a0, a1 = (float(q) for q in lin_coeffs(pred[(kind, j, part)]))
                vals, devs = [], []
                for b in betas:
                    r = zz['betas'][b]
                    q = Fraction(b)
                    x = (1+mp.mpf(q.numerator)/q.denominator)/2
                    p = a0+a1*mp.log(x*(1-x))
                    ratios = [2*cplx(r['hel'][h]['2l_'+part]['f'][j])/cplx(r['hel'][h][key]) for h in hels]
                    vals.append(float(mp.re(ratios[0])))
                    devs.append(float(max(abs(v-p) for v in ratios)))
                first = k == 0                       # highest power: open circles, as in Fig. 3
                mk = dict(marker='o', ms=8, mfc='none', mec=colour[part], mew=1.4) if first else \
                    dict(marker='s', ms=4.5, mfc=colour[part], mec=colour[part])
                curve = a0+a1*Lg
                top.plot(grid, curve, ls='-' if a1 else (0, (5, 2.5)), color=colour[part], lw=1.1)
                top.plot(bf, vals, ls='none', **mk)
                ys += list(curve)+vals
                labels.append((curve[-1], f'${gname[j]}$, ${pname[part]}$'))
                if a1:                               # the same coefficient for the gluon scale (no L_x)
                    alts.setdefault(a0, []).append((j, part))
                    ys.append(a0)
                bot.plot(bf+(2*k+i-1.5)*0.010, devs, ls='none',
                         **(dict(marker='o', ms=6, mfc='white', mec=colour[part], mew=1.3) if first else
                            dict(marker='s', ms=5, mfc=colour[part], mec='white', mew=0.6)))
                handles.append(Line2D([], [], ls='none', label=f'${gname[j]}$, ${pname[part]}$',
                                      **(dict(marker='o', ms=6, mfc='white', mec=colour[part], mew=1.3) if first else
                                         dict(marker='s', ms=5, mfc=colour[part], mec='white', mew=0.6))))
        y0, y1 = min(ys), max(ys)
        pad = 0.07*(y1-y0)
        top.set_ylim(y0-pad, y1+pad)
        span = y1-y0+2*pad
        for a0, parts in alts.items():               # coinciding lines alternate their colours
            for n, (j, part) in enumerate(parts):
                top.plot(grid, np.full_like(grid, a0), color=colour[part], lw=1.3,
                         ls=(1.8*n, (1.0, 1.8*len(parts)-1.0)) if len(parts) > 1 else (0, (1.0, 1.6)))
            powers_here = {j for j, _ in parts}
            assert len(powers_here) == 1
            who = ' and '.join(f'${pname[part]}$' for _, part in parts)
            text = f'${gname[powers_here.pop()]}$, {who}, for $Q_*^2=m^4/s$'
            # label at the right end of the dotted line, on the side with more room over the width
            # of the label (beta > 0.4)
            near = [c for c in top.get_lines() if len(c.get_xdata()) == len(grid) and c.get_ydata()[0] != a0]
            right = grid > 0.4
            room = lambda s: min(min(abs(c.get_ydata()[right]-(a0+s*f*span)).min() for f in (0.02, 0.05, 0.08))
                                 for c in near)
            s = 1 if room(1) >= room(-1) else -1
            top.text(B1-0.01, a0+s*0.022*span, text, fontsize=7, color=INK,
                     va='bottom' if s > 0 else 'top', ha='right')
        for (y, text), yl in zip(labels, spread_labels([y for y, _ in labels], 0.065*span)):
            top.text(B1+0.012, yl, text, fontsize=8, color=INK, va='center', ha='left', clip_on=False)
        frame(top, title)
        frame(bot, dtitle)
        for a in (top, bot):
            a.set_xlim(B0, B1)
        bot.set_yscale('log')
        bot.set_xlabel(r'$\beta$')
        bot.legend(handles=handles, loc='upper left', fontsize=7.5, ncol=2, handletextpad=0.3, columnspacing=1.0)
    ax[0, 0].set_ylabel('coefficient')
    ax[1, 0].set_ylabel(r'$|\mathrm{extracted}-\mathrm{predicted}|$')
    lo = min(min(l.get_ydata()) for a in ax[1] for l in a.get_lines())
    hi = max(max(l.get_ydata()) for a in ax[1] for l in a.get_lines())
    for a in ax[1]:
        a.set_ylim(10**np.floor(np.log10(lo)-0.3), 10**np.ceil(np.log10(hi)+1.3))
    fig.tight_layout(w_pad=4.5)
    save(fig, 'zz_coefficients')


# ---------------------------------------------------------------- Fig. 5: fit windows
def windows_figure():
    """Relative deviation of a coefficient of the next-to-highest power of G from its prediction at
    beta = 1/2, for fits through windows of consecutive deltas, against the smallest delta of the
    window, for every fit model: no p_T^2 terms, or p_T^2 times a polynomial of degree P in G.
    The last point of the model of the paper (P = 2 for gg -> HH, P = 3 for two-loop gg -> ZZ) is
    the fit of Tables 2 and 7, and the point before it the fit with the window moved one step up
    in delta, from which the numerical uncertainty is estimated. The models keep their colour in all panels."""
    hh, zz = load('windows'), load('zz_windows')
    b = preferred_beta(list(hh))
    models = {'plain': ('no $p_T^2$ terms', 'v', INK2), 'pT2': ('$P=0$', 'D', C1), 'pT2 G': ('$P=1$', '^', C2),
              'pT2 G^2': ('$P=2$', 's', C3), 'pT2 G^3': ('$P=3$', 'o', C5)}
    panels = ((hh[b]['spin0 G^2']['rows'], 'pT2 G^2', r'(a) $gg\to HH$: $c_{12}$, $N^2$ term'),
              (hh[b]['spin2 G^1']['rows'], 'pT2 G^2', r'(b) $gg\to HH$: $c_{21}$, $N^2$ term'),
              ([r for r in zz[b]['++++ 2l_CF']['rows'] if r['power'] == 2], 'pT2 G^3',
               r'(c) $gg\to ZZ$: $\mathcal{G}^2$ in $\mathcal{M}_{{+}{+}{+}{+}}$, $C_F$ term'),
              ([r for r in zz[b]['+-++ 2l_CF']['rows'] if r['power'] == 1], 'pT2 G^3',
               r'(d) $gg\to ZZ$: $\mathcal{G}$ in $\mathcal{M}_{{+}{-}{+}{+}}$, $C_F$ term'))
    fig, ax = plt.subplots(2, 2, figsize=(7.0, 5.2), sharex=True, sharey=True)
    ks = set()
    for a, (rows, paper, title) in zip(ax.flat, panels):
        for m, (lab, mk, col) in models.items():
            sel = [r for r in rows if r['model'] == m]
            if not sel:
                continue
            k = [r['delta_exp_min'] for r in sel]      # the window ends at delta = 10^-k
            ks |= set(k)
            main = m == paper
            a.semilogy(k, [float(mp.mpf(r['rel_dev'])) for r in sel], marker=mk, color=col,
                       lw=2.6 if main else 1.2, ms=7 if main else 4.8, mec='white', mew=0.7,
                       zorder=3 if main else 2)
        frame(a, title)
    kk = list(range(min(ks), max(ks)+1))
    ax[1, 0].set_xticks(kk)
    ax[1, 0].set_xticklabels([f'$10^{{-{k}}}$' for k in kk])
    ax[1, 0].set_xlim(kk[0]-0.4, kk[-1]+0.4)
    lo = min(l.get_ydata().min() for a in ax.flat for l in a.get_lines())
    hi = max(l.get_ydata().max() for a in ax.flat for l in a.get_lines())
    e0, e1 = int(np.floor(np.log10(lo)-0.3)), int(np.ceil(np.log10(hi)+0.3))
    ax[0, 0].set_ylim(10.0**e0, 10.0**e1)
    ax[0, 0].set_yticks([10.0**e for e in range(e0 + (e0 % 2), e1+1, 2)])
    ax[0, 0].minorticks_off()
    for a in ax.flat:
        a.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    for a in ax[1]:
        a.set_xlabel(r'smallest $\delta$ in the fit window')
    for a in ax[:, 0]:
        a.set_ylabel('relative deviation')
    handles = [Line2D([], [], marker=mk, color=col, lw=1.2, ms=4.8, mec='white', mew=0.7, label=lab)
               for lab, mk, col in models.values()]
    handles.append(Line2D([], [], color=INK, lw=2.6, label='used in the text'))
    fig.legend(handles=handles, loc='upper center', ncol=6, fontsize=8, handlelength=2.2,
               columnspacing=1.3, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, 'fit_windows')


if __name__ == '__main__':
    forward_figure()
    hh_figure()
    zz_figure()
    windows_figure()
