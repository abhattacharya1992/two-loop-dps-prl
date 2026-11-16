# Analytic derivation of the finite matching functions r0

This folder derives the one-loop collinear function $J_1$ of the DPS kernel in closed form, diagram by diagram, for $g\to q\bar q$ and $g\to gg$. Together with the soft function and the daughter self-energies it gives the finite matching functions of Eqs. (7) and (8) of the Letter,

```math
\begin{aligned}
r_0^{g\to gg}&=C_A\Big[\tfrac{134}{9}-5\zeta_2-4\ln x\ln(1-x)\Big]-\tfrac{40}{9}T_Fn_f+\frac{2(C_A-2T_Fn_f)}{3F_\nu(x)},\\
r_0^{g\to q\bar q}&=C_F\Big[12-6\zeta_2+2\ln^2\frac{x}{1-x}\Big]+C_A\Big[2+\zeta_2-2\ln^2x-2\ln^2(1-x)\Big].
\end{aligned}
```

The calculation follows the section "One-loop correction to the DPS kernel" of the Supplement and uses its conventions. The gluon colour factors are taken from the existing graph calculation in
`paper2/r0/theory/gg_diagrams.py` and checked against its stored results.
The tensor numerators, residues and collinear integrals are evaluated here.
The numerical evaluation in `paper2/r0/theory` does not enter the derivation. Its stored Laurent coefficients (`gqq_laurent_x*.json`, `gg_laurent_x*.json`) serve as a comparison for every diagram.

## How to run

```
make -C paper2/r0/derivation test
```

runs the six scripts in turn and stops at the first failure. A single script takes the repository root as its argument:

```
python3 paper2/r0/derivation/check_gg.py .
```

The scripts need Python 3 with numpy and scipy 1.10 or later, and nothing else. They were run with Python 3.13, numpy 2.5 and scipy 1.18. Every check prints one line marked `[ok]` or `[FAIL]`, and a script exits with status 1 if any check fails. A complete run has 163 checks and takes about three minutes; two and a half of these are the quadratures in two variables in `check_qq_Vggq_W1g_W2g.py`. The folder `outputs` holds the output of each script, and `make -C paper2/r0/derivation outputs` regenerates it.

## Conventions

The bare one-loop correction, normalised to the lowest-order kernel, is

```math
\widetilde{\mathcal K}^{(1)}_{\rm bare}/\widetilde{\mathcal K}_0=J_1+J_2+S+\sum_{\rm daughters}\Sigma+G ,
```

with $p_1^+=\mu=\nu_{\rm r}=1$, $b=b_0$ and $a_s=\alpha_s/(4\pi)$. We work in Feynman gauge in $d=4-2\epsilon$ dimensions with the rapidity regulator $\eta$ of the Supplement, and take $\eta\to0$ before $\epsilon\to0$. The function $r_0$ is the coefficient of $\eta^0\epsilon^0$ in $J_1+J_2+S+\sum\Sigma$ for the helicity structures present at tree level. The Glauber contribution $G$ is imaginary and is not treated here.

The parent gluon has momentum $p_1=(1,0,\boldsymbol 0)$ and polarisation $\varepsilon$. The daughters carry $K=(x,K^-,\boldsymbol k)$ and $p_1-K$, and the loop momentum is $\ell=(z,\ell^-,\boldsymbol\ell)$. Every diagram is evaluated in the same order:

1. Numerator. Dirac algebra in $d$ dimensions for $g\to q\bar q$, contraction of the Lorentz tensors for $g\to gg$.
2. Residues in $K^-$ and $\ell^-$. This leaves a one-loop light-front wave function with at most two transverse denominators.
3. Transverse integral in $2-2\epsilon$ dimensions with one Feynman parameter.
4. Parameter integrals, expanded in $\eta$ and then in $\epsilon$.
5. Transform to impact parameter space at $b=b_0$. A term linear in $\boldsymbol k$ acquires the factor $F_\epsilon=e^{-\epsilon\gamma_E}\Gamma(1-2\epsilon)/[\Gamma(1+\epsilon)\Gamma(1-\epsilon)]$.

The second region of the loop plus momentum, and the diagrams with the two operators interchanged, follow from charge conjugation for $g\to q\bar q$ and from the exchange of the two daughters for $g\to gg$. Both replace $x$ by $1-x$ and interchange the two quark helicity structures, or $S_2$ and $S_3$.

Below, $L=\ln[x(1-x)]$.

## Results for g → q q̄

Each diagram equals its colour factor times $A\,T+\tfrac12B\,(\varepsilon\cdot\boldsymbol k)$, where $T=\tfrac x2\,\not\varepsilon\,\not k-\tfrac{1-x}2\,\not k\,\not\varepsilon$ is the tree splitting. The coefficients of the two quark helicity structures $\not\varepsilon\,\not k$ and $\not k\,\not\varepsilon$, normalised to the tree splitting, are $A-B/(2x)$ and $A+B/(2(1-x))$.

Wilson line couplings to the quark and to the antiquark:

```math
\begin{aligned}
\text{W1q}\ \ (C_F):\quad &A=G(\epsilon)\,(1-x)^{\epsilon},\qquad B=0,\\
\text{W2qb}\ \ (C_F):\quad &A=G(\epsilon)\,x^{\epsilon},\qquad B=0,\\
&G(\epsilon)=-2F_\epsilon\,\frac{\Gamma(\epsilon)\Gamma(-\epsilon)\Gamma(2-\epsilon)}{\Gamma(2-2\epsilon)}=\frac2{\epsilon^2}+\frac2\epsilon+4+\zeta_2+O(\epsilon),\\
\text{W2q}\ \ (C_F-\tfrac12C_A):\quad &A=\frac2{\eta\epsilon}-\frac2{\epsilon^2}-\frac{2L}\epsilon+2\,\mathrm{Li}_2(1-x)-3\zeta_2,\qquad B=4x\ln x,\\
\text{W1qb}\ \ (C_F-\tfrac12C_A):\quad &A=\frac2{\eta\epsilon}-\frac2{\epsilon^2}-\frac{2L}\epsilon+2\,\mathrm{Li}_2(x)-3\zeta_2,\qquad B=-4(1-x)\ln(1-x).
\end{aligned}
```

In W1q and W2qb the gluon from the Wilson line of an operator attaches to the daughter created by the same operator. In W2q and W1qb it attaches to the daughter created by the other operator, and only these two diagrams have a rapidity divergence.

Abelian and non-Abelian vertex corrections, and Wilson line couplings to the parent gluon:

```math
\begin{aligned}
\text{Vqq}\ \ (C_F-\tfrac12C_A):\quad &A=\frac1\epsilon-L,\qquad B=-4x\ln x+4(1-x)\ln(1-x),\\
\text{Vggq}\ \ (\tfrac12C_A):\quad &A=\frac{2-2L}\epsilon+4-\ln^2x-\ln^2(1-x),\\
&B=\frac{-2x\ln x+2(1-x)\ln(1-x)}\epsilon-x\ln^2x+(1-x)\ln^2(1-x),\\
\text{W1g}\ \ (\tfrac12C_A):\quad &A=-\frac2{\epsilon^2}-\frac{2\ln(1-x)+x}\epsilon-3\zeta_2-\ln(1-x)-\ln^2(1-x)+2\,\mathrm{Li}_2(x)-2x,\\
&B=\frac{2x(1-x)-2(1-x)\ln(1-x)}\epsilon+4x(1-x)-(1-x)\ln^2(1-x),\\
\text{W2g}\ \ (\tfrac12C_A):\quad &A_{\rm W2g}(x)=A_{\rm W1g}(1-x),\qquad B_{\rm W2g}(x)=-B_{\rm W1g}(1-x).
\end{aligned}
```

In the sum, $B$ cancels separately in the $C_F$ and $C_A$ parts, so that $J_1$ is proportional to the tree splitting:

```math
J_1^{g\to q\bar q}=C_F\Big[\frac4{\eta\epsilon}+\frac{5-2L}\epsilon+8-2\zeta_2+L+\ln^2\frac{x}{1-x}\Big]
+C_A\Big[-\frac2{\eta\epsilon}+1-\ln^2x-\ln^2(1-x)\Big].
```

## Results for g → gg

The helicity structures are

```math
S_1=(\varepsilon\cdot\boldsymbol k)(e_1\cdot e_2),\quad
S_2=(\varepsilon\cdot e_1)(e_2\cdot\boldsymbol k),\quad
S_3=(\varepsilon\cdot e_2)(e_1\cdot\boldsymbol k),\quad
T_4=\frac{(\varepsilon\cdot\boldsymbol k)(e_1\cdot\boldsymbol k)(e_2\cdot\boldsymbol k)}{\boldsymbol k^2},
```

and the tree splitting is proportional to $S_1-S_2/(1-x)-S_3/x$. We write $c_1,c_2,c_3$ for the coefficients of $S_1,S_2,S_3$ normalised to the tree splitting, and $c_4$ for the coefficient of $T_4$ in units of the tree coefficient of $S_1$. Where $c_4$ is not listed it vanishes.

Wilson line couplings to the daughter created by the same operator, to the parent gluon, and to the other daughter:

```math
\begin{aligned}
\text{W1k}\ \ (C_A):\quad &c_1=c_2=c_3=\frac2{\epsilon^2}+\frac{1+2\ln(1-x)}\epsilon+2+\zeta_2+\ln(1-x)+\ln^2(1-x),\\
\text{W1g}\ \ (\tfrac12C_A):\quad &c_1=-\frac2{\epsilon^2}-\frac{3\ln(1-x)}\epsilon-3\zeta_2+2\,\mathrm{Li}_2(x)-\tfrac12\ln(1-x)-\tfrac32\ln^2(1-x),\\
&c_2=-\frac2{\epsilon^2}-\frac{2\ln(1-x)+x}\epsilon-3\zeta_2+2\,\mathrm{Li}_2(x)-\ln(1-x)-\ln^2(1-x)-2x,\\
&c_3=-\frac2{\epsilon^2}-\frac{\tfrac52\ln(1-x)+\tfrac12x}\epsilon-3\zeta_2+2\,\mathrm{Li}_2(x)+\tfrac34(1-x)\ln(1-x)-\tfrac54\ln^2(1-x)-\tfrac54x,\\
\text{W1kb}\ \ (-\tfrac12C_A):\quad &c_1=-\frac2{\eta\epsilon}+\frac2{\epsilon^2}+\frac{2L}\epsilon+3\zeta_2-2\,\mathrm{Li}_2(x)+\ln(1-x),\\
&c_2=-\frac2{\eta\epsilon}+\frac2{\epsilon^2}+\frac{2L+\ln(1-x)+\tfrac12x}\epsilon+3\zeta_2-2\,\mathrm{Li}_2(x)+\tfrac12\ln(1-x)+\tfrac12\ln^2(1-x)+x,\\
&c_3=-\frac2{\eta\epsilon}+\frac2{\epsilon^2}+\frac{2L}\epsilon+3\zeta_2-2\,\mathrm{Li}_2(x)+\big(\tfrac74-\tfrac34x\big)\ln(1-x)-\tfrac14x.
\end{aligned}
```

The coefficient of $1/\eta$ in W1kb is $2e^{-\epsilon\gamma_E}\Gamma(-\epsilon)$ times the tree splitting to all orders in $\epsilon$, because the numerator at the eikonal endpoint is proportional to the tree splitting.

Four-gluon vertex bubbles on the first daughter, and the gluon, ghost and quark triangles:

```math
\begin{aligned}
\text{Vbuba1, 2, 3}\ \ (-C_A,\,-\tfrac12C_A,\,\tfrac12C_A):\quad &c_1=c_2=\kappa\Big[\frac1\epsilon+2+\ln(1-x)\Big],\quad \kappa=\tfrac34,\ \tfrac38,\ -\tfrac38,\qquad c_3=0,\\
\text{Vtri}\ \ (\tfrac12C_A):\quad &c_1=\frac{\tfrac54-L}\epsilon+\tfrac73-\tfrac32L-\tfrac12\big[\ln^2x+\ln^2(1-x)\big],\\
&c_2=\frac{x-\tfrac12-\tfrac32\ln x-\ln(1-x)}\epsilon+\tfrac54x-\tfrac5{12}-\tfrac{17}4\ln x-\tfrac34\ln^2x-\tfrac32\ln(1-x)-\tfrac12\ln^2(1-x),\\
&c_3(x)=c_2(1-x),\qquad c_4=\tfrac32,\\
\text{Vgh1}+\text{Vgh2}\ \ (C_A):\quad &c_1=-\tfrac1{24}\Big(\frac1\epsilon+L\Big)-\tfrac19,\quad c_2=c_1-\tfrac1{24}(1-x),\quad c_3=c_1-\tfrac1{24}x,\quad c_4=-\tfrac1{12},\\
\text{Vq1}+\text{Vq2}\ \ (T_Fn_f):\quad &c_1=\tfrac43\Big(\frac1\epsilon+L\Big)+\tfrac{20}9,\quad c_2=c_1-\tfrac23(1-x),\quad c_3=c_1-\tfrac23x,\quad c_4=-\tfrac43.
\end{aligned}
```

The three entries Vbuba1, 2, 3 are the three colour and Lorentz parts of the four-gluon vertex. The entries W2g, W2kb, W2k and Vbubb1, 2, 3 follow from W1g, W1k, W1kb and Vbuba1, 2, 3 by exchanging the two daughters. The six entries in which the gluon from the Wilson line enters a four-gluon vertex vanish. The colour sum of the 23 entries is

```math
\begin{aligned}
c_s&=C_A\Big[\frac2{\eta\epsilon}+\frac{\tfrac13-2L}\epsilon+\tfrac59-2\zeta_2-\tfrac53L-2\ln x\ln(1-x)+\tfrac13\rho_s\Big]
+T_Fn_f\Big[\frac4{3\epsilon}+\tfrac{20}9+\tfrac43L-\tfrac23\rho_s\Big],\\
c_4&=\tfrac23\big(C_A-2T_Fn_f\big),\qquad \rho_1=0,\quad\rho_2=1-x,\quad\rho_3=x .
\end{aligned}
```

## From the collinear function to r0

At $\nu_{\rm r}=1$ and $b=b_0$ the soft function of the Supplement is $S=\Delta_C\,[-4/(\eta\epsilon)+2/\epsilon^2-\zeta_2]+O(\epsilon)$. The rapidity divergences cancel: $J_1$ and $J_2$ each contribute $2\Delta_C/(\eta\epsilon)$, with $\Delta_C=2C_F-C_A$ for $g\to q\bar q$ and $\Delta_C=C_A$ for $g\to gg$. The self-energies of the two daughters of one splitting are

```math
\Sigma_q+\Sigma_{\bar q}=C_F\Big[-\frac2\epsilon-2-L\Big],\qquad
\Sigma_g+\Sigma_g=C_A\Big[\frac{10}{3\epsilon}+\tfrac{62}9+\tfrac53L\Big]+T_Fn_f\Big[-\frac8{3\epsilon}-\tfrac{40}9-\tfrac43L\Big].
```

The function $J_2$ is $J_1$ at $1-x$ with the helicity structures interchanged, and the sums above are symmetric, so the finite parts of $J_1$ and $J_2$ agree. For $g\to q\bar q$ this gives

```math
r_0^{g\to q\bar q}=2\,J_1^{g\to q\bar q}\Big|_{\epsilon^0}-\zeta_2\,(2C_F-C_A)+2C_F\,(-2-L),
```

which is Eq. (8).

For $g\to gg$ the terms $\rho_s$ and $c_4$ depend on the helicities. Two steps turn them into the last term of Eq. (7).

Transform of $T_4$. The tensor $T_4-(S_1+S_2+S_3)/4$ is symmetric and traceless in its three polarisation indices, so it is the third harmonic in the direction of $\boldsymbol k$, while the tree splitting is the first. Both have the radial dependence $1/|\boldsymbol k|$, for which the transform of the harmonic $m$ is $i^m/b$ times a factor common to all $m$. Relative to the tree splitting the third harmonic changes sign, so that at $\epsilon=0$

```math
T_4\ \to\ \tfrac12(S_1+S_2+S_3)-T_4 ,\qquad \boldsymbol k\to\boldsymbol b .
```

Contraction of the two splittings. Each hard vertex joins one daughter of each splitting, and the daughter that meets the one with fraction $x$ carries $1-x$. The two daughters at a vertex are contracted with $\delta_{ij}$ for $gg\to HH$ and with the antisymmetric symbol of the two transverse dimensions for $gg\to AA$. The sum over the daughter polarisations has the form $\alpha\,(\varepsilon\cdot\varepsilon')+\beta\,(\varepsilon\cdot\boldsymbol b)(\varepsilon'\cdot\boldsymbol b)/b^2$, with harmonics $\alpha+\beta/2$ for $\nu=0$ and $\beta$ for $\nu=2$. The tree overlaps are

```math
gg\to HH:\quad \frac{F_0}{2x(1-x)},\ \ \frac{F_2}{x(1-x)};\qquad\qquad
gg\to AA:\quad -\frac{F_A}{2x(1-x)},\ \ \frac{F_2}{x(1-x)}
```

for $\nu=0$ and $\nu=2$. Dividing the one-loop overlap by the tree overlap and adding the soft function and the four self-energies gives Eq. (7) for both harmonics. For $gg\to AA$, $F_A$ replaces $F_0$ and the last term changes sign for $\nu=2$, as stated below Eq. (7) of the Letter.

The single pole of $J_1+J_2+\sum\Sigma$ that does not depend on $x$ is $B/\epsilon$, with $B=6C_F$ for $g\to q\bar q$ and $B=2\beta_0$ for $g\to gg$, the coefficients of the Letter.

## What is checked

| Script | Content | Time |
|---|---|---|
| `check_qq_W1q_W2qb.py` | W1q, W2qb | 1 s |
| `check_qq_W1qb_W2q.py` | W2q, W1qb | 1 s |
| `check_qq_Vqq.py` | Vqq | 15 s |
| `check_qq_Vggq_W1g_W2g.py` | Vggq, W1g, W2g | 160 s |
| `check_qq_sum.py` | sum of the eight diagrams, Eq. (8) | 1 s |
| `check_gg.py` | the 23 entries of the gluon channel, colour sum, contraction of the two splittings, Eq. (7) | 5 s |

There are two kinds of check.

Equalities that hold for all $x$. In `check_gg.py` and `check_qq_sum.py` the momentum fraction is a symbol, and every coefficient is a rational function of $x$ multiplying $1$, $\ln x$, $\ln(1-x)$, their squares and product, $\mathrm{Li}_2(x)$ or $\zeta_2$. The gluon channel is computed in this way from the numerators to the Laurent coefficients. The closed forms above, the colour sum, the single poles and Eq. (7) are tested as equalities of these coefficients. For Eq. (7) this includes the term with $F_\nu$ for $\nu=0$ and $2$, with the contraction of the two splittings done in exact arithmetic, and the corresponding statement for $gg\to AA$. The sum of the eight quark diagrams and Eq. (8) are tested in the same way from the closed forms of the individual diagrams.

Numerical checks with a stated tolerance. These cover the Dirac algebra of the quark diagrams in $d=4,6,8,10$ with random momenta (the result is linear in $d$), the residues in $K^-$ and $\ell^-$ against quadrature in two variables, the parameter integrals against quadrature at finite $\epsilon$ and $\eta$, and the comparison with the stored Laurent coefficients at $x=0.1,0.3,0.5,0.7,0.9$. For $g\to q\bar q$, the individual parameter-integral Laurent expansions were obtained analytically and checked numerically at finite $\epsilon$.
The sum of the closed diagram expressions is tested as an exact equality in $x$.

The closed forms agree with the stored coefficients to $10^{-9}$ or better for the quark diagrams, except Vqq, whose stored values are accurate to $10^{-7}$. In the gluon channel the largest deviation is $1.4\times10^{-8}$, except for the gluon triangle: $8.7\times10^{-7}$ with the default quadrature orders, and $1.1\times10^{-8}$ at $x=1/2$ with the orders 96, 72 and 40.

## Files

| File | Content |
|---|---|
| `algebra.py` | polynomials, rational functions of $x$, Laurent series in $\epsilon$, Lorentz tensors in light-cone components; standard library only |
| `gluon_diagrams.py` | tree splitting and the 23 one-loop entries of `paper2/r0/theory/gg_diagrams.py`, and their numerators after the residues |
| `integrals.py` | transverse integral, reduction to $S_1,S_2,S_3,T_4$, parameter integrals |
| `gluon_channel.py` | each entry of the gluon channel as a Laurent series in $\eta$ and $\epsilon$ |
| `dirac.py` | gamma matrices in even dimension, for the quark diagrams |
| `supplement_formulas.py` | $F_\epsilon$, the soft function and the daughter self-energy of the Supplement |
| `report.py` | pass or fail line of each check, exit status |
| `check_*.py` | the six scripts of the table above |
| `outputs/` | output of the six scripts |

## Not covered

- The Glauber contribution $G$ and the imaginary part of the kernel.
- A separate calculation of $J_2$. It is obtained from $J_1$ by charge conjugation or by the exchange of the two splittings.
- A derivation of the daughter self-energies. Their expression is taken from the Supplement and expanded. For the soft function the rapidity integral and the expansion are checked.
