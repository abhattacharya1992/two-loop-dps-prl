# Unequal-virtuality test: fixed predictions and conventions

This is supporting research for the `unequal-virtuality` branch. It does not
change the publication Letter or Supplement and does not establish a complete
factorisation theorem. The test concerns the highest two logarithms in
isotropic helicity amplitudes with a nonzero one-loop DPS logarithm.

## 1. External kinematics

Set `a = ma2 = p3^2 > 0`, `b = mb2 = p4^2 > 0`, and
`s > (sqrt(a)+sqrt(b))^2`. All three inputs have mass-squared units.
The massless incoming momenta are p1 and p2, with

    s=(p1+p2)^2, t=(p1-p3)^2, u=(p2-p3)^2,
    s+t+u=a+b.

Define

    D = sqrt((s-a-b)^2-4*a*b),
    t_plus  = (a+b-s+D)/2,
    t_minus = (a+b-s-D)/2,
    t_plus*t_minus = a*b.

The forward endpoint is `t=t_plus`, `u=t_minus`. The implementation evaluates
`t_plus=a*b/t_minus` to avoid subtraction of nearly equal numbers.
Choose a dimensionful distance `d=sqrt(a*b)*10^(-k)` and approach the endpoint by

    t=t_plus-d, u=t_minus+d,
    pT2=(t*u-a*b)/s = d*(D-d)/s.

The last product is used numerically; direct evaluation of `t*u-a*b` is avoided.
The scan requires `0<d<D/2`, so it stays in the forward half of the physical region.
The exactly singular endpoint is never evaluated.

In the centre-of-mass frame:

    Ein=sqrt(s)/2, E3=(s+a-b)/(2*sqrt(s)), E4=(s+b-a)/(2*sqrt(s)),
    |p3|=|p4|=D/(2*sqrt(s)),
    cos(theta)=1-2*d/D, sin(theta)=sqrt(pT2)/|p3|.

The outgoing polarization vectors have different energies and masses. They
satisfy p.epsilon=0 and epsilon.epsilon*=-1; the incoming polarizations also
satisfy epsilon1.p2=epsilon2.p1=0, as required by the published tensor basis.
The code preserves the equal-mass implementation's phases and outgoing
complex conjugation. Timelike virtual-boson polarizations are projections of
the current tensor, not a new approximation to a fully decayed cross section.

At the forward pinch, write `p3=x1*p1+x2*p2`. Then

    x1=(s+a-b+D)/(2*s), x2=(s+a-b-D)/(2*s),
    a=s*x1*x2, b=s*(1-x1)*(1-x2), D=s*(x1-x2).

For equal masses, `x2=1-x1`; unequal masses remove this constraint. The
initial benchmark a=1, b=2, s=15/2 has D=7/2, t_plus=-1/2,
t_minus=-4, x1=2/3 and x2=1/5. The fractions are away from their endpoints.

## 2. Logarithmic convention

Use the reference scale `Q0^2=a*b/s` and the complex logarithm

    G=log(pT/Q0)+i*pi/2 = 0.5*log(s*pT2/(a*b))+i*pi/2,
    L=log(sqrt(a*b)/s).

Both arguments are dimensionless. At a=b=1 this is exactly the original
`G=0.5*log(1-t*u)` on the principal physical branch, since t*u>1.
Under s,a,b -> c*s,c*a,c*b, pT2 -> c*pT2 and G and L are unchanged.

For quark daughters, the Letter gives

    DC=2*CF-CA, B=6*CF,
    DC*log(Qstar^2)=CF*log(a)+CF*log(b)-CA*log(s).

Only ratios of scales enter differences of these logarithms. Direct algebra gives

    Gq-G = 0.5*log(Q0^2/Qstar^2) = -(CA-CF)*L/DC.

For SU(Nc), Qstar^2 = s^(Nc^2)/(a*b)^((Nc^2-1)/2). This expression
has mass-squared dimension, despite its large exponent.

## 3. Conditional two-loop predictions

The prediction uses the SAME boundary condition r1=r2=0 as the existing
Letter. Its validity for unequal virtualities is part of this numerical test;
it is not inferred merely from the renormalisation group equation.

In the Letter's alpha_s/(4*pi) convention the top two logarithms are

    M2 = C1*(-4*DC*Gq^3/3 - B*Gq^2) + O(Gq).

The VVamp evaluator has alpha_s/(2*pi) per loop. If c1 is its one-loop
coefficient of G, the extra relative factor of 1/2 gives

    f3_CF/c1 = -4/3,       f2_CF/c1 = -2*L-3,
    f3_CA/c1 =  2/3,       f2_CA/c1 =  2*L.

Here M2_CF=f3_CF*G^3+f2_CF*G^2+... and similarly for CA, with the colour
factors NOT yet multiplied into the stored coefficients. The change from
Gq to G contributes to f2 and must not be omitted. These are ratios to the
same helicity's one-loop logarithmic coefficient, not ratios to its full
one-loop amplitude. The extraction fits the one-loop coefficient independently;
no guessed unequal-mass closed form for c1 is imposed.

At the benchmark, L=-1.66832943026229210187015638796144849205984..., so

| coefficient/c1 | CF | CA |
|---|---:|---:|
| G^3 | -1.3333333333333333... | +0.6666666666666667... |
| G^2 | +0.33665886052458420374... | -3.33665886052458420374... |

These targets and the fit plan were committed in `predictions_stage1.json`
before extracting any unequal-virtuality amplitude. The prediction module
has no dependency on the two-loop numerical results.

## 4. Published amplitude and crossings

The source is the massless-quark class-A qT finite remainder of arXiv:1503.08835,
at mu^2=s, with electroweak factors and the overall colour tensor stripped.
This is not a massive-top prediction or a complete physical ZZ/WW cross section.

The source computes nine of the twenty tensor coefficients. Its
`ggvvampprec.cpp::employ_crossing_relations` requires four configurations:

    A      : (s,t,a,b)
    A12    : (s,u,a,b)
    A34    : (s,u,b,a)
    A12x34 : (s,t,b,a)

The missing coefficients are

    A3=A12_2, A6=A12_5, A7=A12_4,
    A10=-A34_8, A11=-A34_9,
    A12=A12_9, A13=A12_8,
    A14=-A12x34_9, A15=-A12x34_8,
    A19=A12_18, A20=A12_17.

The equal-mass shortcuts A34=A12 and A12x34=A are NOT used at unequal
virtualities. Tests include a deliberate shortcut error. The shared special-
function cache reduces repeated work without identifying different arguments.
The existing arbitrary-precision translation is reused unchanged; extending
its kinematics is not an independent implementation of the source amplitude.

## 5. Extraction and validation

At each fixed (s,a,b), fit each complex amplitude WITHOUT coefficient constraints:

    sum(f_j*G^j,j=0..J) + rho*sum(h_l*G^l,l=0..P),
    rho=pT2/sqrt(a*b).

J=1 for one loop, J=3 for two loops. One loop uses P=2. Two-loop fits
compare P=-1 (no power term), 0, 1, 2, 3. Central interpolation uses the
deepest J+P+2 points, with a one-point-shifted window; one extra point gives
an overdetermined complex least-squares check. Column rescaling is purely
numerical and does not impose theory. A deliberately wrong synthetic cubic
coefficient is recovered as wrong rather than forced toward its prediction.

The primary P=3 two-loop interpolation uses eight points; nine allow its
shifted-window comparison. Working amplitude precision is 30+5*k decimal
digits. A repeat at k=10 and 95 digits checks the 80-digit point. Higher fit
precision does not restore information lost in the input amplitudes.

Window shifts are numerical stability diagnostics, not rigorous uncertainty
bounds. They cannot exclude a common truncation bias. Model dependence,
precision comparisons and residuals must be assessed together. No fit
claims a probability of agreement or a statistically calibrated confidence interval.

## 6. Scope and references

This extension tests a prediction; it does not assume a universal two-variable
finite matching r0, derive new two-loop integrals, prove operator factorisation,
or change the Letter's scope. Additional virtuality ratios and energies are
needed before treating the single benchmark as a general validation.

Primary amplitude reference: A. von Manteuffel and L. Tancredi,
*The two-loop helicity amplitudes for gg -> V1 V2 -> 4 leptons*,
arXiv:1503.08835, https://arxiv.org/abs/1503.08835.
One-loop DPS reference: J. R. Gaunt and W. J. Stirling,
*Double Parton Scattering Singularity in One-Loop Integrals*,
arXiv:1103.1888, https://arxiv.org/abs/1103.1888.
Internal theoretical input: `../manuscript_prl/letter.tex`, Eqs. labelled
`eq:gammaK`, `eq:kernel`, `eq:rg`, at base main 3d59dd953735bfbe2e3f8ca48dfd12b9e7b86195.
