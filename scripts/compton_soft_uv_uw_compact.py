"""Collect the independently integrated UW sector using exact shuffle identities."""
from compton_artifacts import read_artifact, write_artifact
import json
from pathlib import Path
import sympy as S

folder = Path(__file__).resolve().parents[1]
X, Y = S.symbols('X Y')
source = json.loads(read_artifact(folder / 'uw_inverse_basis.json'))
assert source['exact_span']
c = {tuple(row['word']): S.sympify(row['coefficient'])
     for row in source['coefficients_X']}
K = c[('-4', 'r', '-1')]
for word, multiple in [(('0', 'r', '-1'), 2), (('r', '-1', '-1'), 1),
                       (('r', '-1', '0'), 2), (('r', '0', '-1'), 1)]:
    assert S.cancel(c[word] - multiple*K) == 0
assert S.cancel(2*c[('-2', '-2')] - c[('-2', '0')]) == 0
coeff = {
    'a_0': c[()], 'a_1': c[('-1',)], 'a_L': c[('0',)],
    'a_{11}': c[('-1', '-1')]/2, 'a_{10}': c[('-1', '0')],
    'a_D': c[('-1', '0')] - c[('0', '-1')],
    'a_2': c[('-2',)], 'a_{D2}': c[('-2', '0')] - c[('0', '-2')],
    'a_{20}': c[('-2', '0')], 'a_R': c[('r', '-1')], 'K': K,
    'a_{011}': c[('0', '-1', '-1')], 'a_{RR}': c[('r', 'r', '-1')],
}
coeff = {k: S.factor(v) for k, v in coeff.items()}
reconstructed = {
    (): coeff['a_0'], ('-1',): coeff['a_1'], ('0',): coeff['a_L'],
    ('-1', '-1'): 2*coeff['a_{11}'],
    ('-1', '0'): coeff['a_{10}'],
    ('0', '-1'): coeff['a_{10}']-coeff['a_D'],
    ('-2',): coeff['a_2'],
    ('0', '-2'): coeff['a_{20}']-coeff['a_{D2}'],
    ('-2', '0'): coeff['a_{20}'],
    ('-2', '-2'): coeff['a_{20}']/2,
    ('r', '-1'): coeff['a_R'],
    ('-4', 'r', '-1'): K, ('0', 'r', '-1'): 2*K,
    ('r', '-1', '-1'): K, ('r', '-1', '0'): 2*K,
    ('r', '0', '-1'): K,
    ('0', '-1', '-1'): coeff['a_{011}'],
    ('r', 'r', '-1'): coeff['a_{RR}'],
}
assert set(reconstructed) == set(c)
assert all(S.cancel(reconstructed[w]-c[w]) == 0 for w in c)
out = [r'\subsection{The integrated UW contribution}',
 r'Put $H=X>0$ and use $Y$ and $q$ defined at the beginning of this chapter. The completed $UW$ integral can be expressed in the same functions as the virtual contribution. Define the real logarithms',
 r'\begin{equation} L_0=\ln X,\qquad L_1=\ln(1+X),\qquad L_2=\ln(1+X/2).\end{equation}',
 r'In the following formula every $G$ has endpoint $X$. Its rational letters have kernels $\dd z/(z-a)$, and the letter $r$ has kernel $\dd z/\sqrt{z(z+4)}$, with positive square root for $z>0$. Define',
 r'\begin{equation} \mathcal U(X)=G(-4,r,-1;X)+G(r,-1,-1;X)-G(r,0,-1;X).\end{equation}',
 r'The unnormalized finite contribution is',
 r'\begin{equation}\begin{aligned}',
 r' F_{UW}^{\eta}(X)={}&a_0+a_1L_1+a_LL_0+a_{11}L_1^2+a_{10}L_0L_1\\',
 r' &+a_D\operatorname{Li}_2(-X)+a_2L_2+a_{D2}\operatorname{Li}_2(-X/2)\\',
 r' &+a_{20}(L_0L_2+L_2^2/4)+(a_R+2KL_0)G(r,-1;X)\\',
 r' &+K\mathcal U(X)+a_{011}G(0,-1,-1;X)+a_{RR}G(r,r,-1;X).',
 r'\end{aligned}\label{eq:comp-uw-integrated-compact}\end{equation}',
 r'Its contribution to $C_1$ is $F_{UW}^{\eta}/[4X(1+X)]$, including the second Bose-related sector. The coefficients are as follows.']
for name, value in coeff.items():
    num, den = S.fraction(value)
    terms = S.Add.make_args(S.expand(num))
    terms = sorted(terms, key=lambda term: S.degree(term, X), reverse=True)
    if len(terms) > 5:
        pieces = []
        for i in range(0, len(terms), 4):
            piece = S.latex(sum(terms[i:i+4]))
            if i and not piece.startswith('-'):
                piece = '+' + piece
            pieces.append('&' + piece)
        rhs = r'\frac{\begin{aligned}' + r'\\'.join(pieces) + r'\end{aligned}}{' + S.latex(den) + '}'
    else:
        rhs = S.latex(value)
    out += [r'\begin{equation}', name + '=' + rhs + '.', r'\end{equation}']
out += [
 r'\paragraph{Integration and endpoint constants.} Starting from the eight-function kernel above, rational integration by parts gives',
 r'\begin{equation} \Phi_{UW}(t)=\frac{\dd}{\dd t}\sum_i Q_i(t)f_i(t)+\sum_i T_i(t)f_i(t).\end{equation}',
 r'Here $Q_i$ are rational functions, and $T_i$ have only simple poles. The identity follows by substituting the derivatives of the eight functions and matching their rational coefficients. At each endpoint the negative powers in the Laurent expansion of the sum $\sum_iQ_if_i$ cancel. Its constant term must include products of the poles of $Q_i$ with the Taylor coefficients of $f_i$. Thus the integrated boundary is the difference of the two constant terms, not the difference of the constant coefficients of the rational functions alone.',
 r'''The simple-pole terms integrate to GPLs with endpoint $t=1$.
Set $X=(1-q)^2/q$ and $\delta=1-q$, with $0<q<1$. The moving letters are
\begin{equation}
 a=-\frac{q}{(1-q)^2},\qquad c=-\frac{2q}{(1-q)^2},\qquad
 r_-=-\frac{q}{1-q},\qquad r_+=\frac1{1-q}.
\end{equation}
They lie outside the interval $[0,1]$. Together with the fixed letters
$0,1$, they specify the endpoint functions whose constants are needed.

For a word $w$, define $g_w(q)$ as the constant coefficient in the
$t\to1^-$ expansion of $G(w;t)$ in powers of $\ln(1-t)$.
The endpoint is regularized at fixed $q$ before taking $q\to1$.
The remaining terms vanish as $(1-t)$ times powers of its logarithm,
so this constant-term operation respects shuffle products.
In particular, $g_{\varnothing}=1$, $g_0=g_1=0$, and all nonempty
pure-zero or pure-one words have value zero. The nonzero fixed-letter
values through weight three are
\begin{align}
 g_{0,1}&=-\zeta(2),&g_{1,0}&=\zeta(2),\\
 g_{0,0,1}&=-\zeta(3),&g_{0,1,0}&=2\zeta(3),\\
 g_{1,0,0}&=-\zeta(3),&g_{0,1,1}&=\zeta(3),\\
 g_{1,0,1}&=-2\zeta(3),&g_{1,1,0}&=\zeta(3).
\end{align}
Three convergent integrals fix this table:
\begin{align}
 G(0,1;1)&=\int_0^1\frac{\ln(1-t)}t\,\dd t=-\zeta(2),\\
 G(0,0,1;1)&=-\operatorname{Li}_3(1)=-\zeta(3),\\
 G(0,1,1;1)&=\frac12\int_0^1\frac{\ln^2(1-t)}t\,\dd t
 =\frac12\int_0^1\frac{\ln^2u}{1-u}\,\dd u=\zeta(3).
\end{align}
The first two follow from the convergent polylogarithm series.
For the third, expanding $1/(1-u)=\sum_{n=0}^{\infty}u^n$ and using
$\int_0^1u^n\ln^2u\,\dd u=2/(n+1)^3$ gives the stated value.
The remaining entries follow from the regularized shuffle relations
\begin{align}
 0&=g_{0,1}+g_{1,0},\\
 0&=2g_{0,0,1}+g_{0,1,0},\\
 0&=g_{1,0,0}+g_{0,1,0}+g_{0,0,1},\\
 0&=g_{1,0,1}+2g_{0,1,1},\\
 0&=g_{0,1,1}+g_{1,0,1}+g_{1,1,0}.
\end{align}
These are the products $g_0g_1$, $g_0g_{0,1}$, $g_1g_{0,0}$,
$g_1g_{0,1}$ and $g_0g_{1,1}$, respectively.

The words containing a moving letter have zero boundary value as
$q\to1$. To establish this before interchanging limits, first use
shuffle products to remove initial letters $1$ and terminal letters $0$.
The factors $g_1$ and $g_0$ vanish. Every surviving word therefore
starts with a letter other than $1$ and ends with a letter other than
$0$, and still contains each moving letter of the original word.
Its ordered integral is absolutely convergent at both endpoints.
For $0<\delta<1/4$, each moving letter $\alpha$ satisfies
\begin{equation}
 \left|\frac1{t-\alpha(q)}\right|\leq\frac2{|\alpha(q)|},
 \qquad 0\leq t\leq1.
\end{equation}
Replacing the moving kernels by these bounds leaves only the fixed
kernels $1/t$, $1/(1-t)$ and the constant kernel in the ordered integral.
The endpoint conditions just established make this integral finite.
Consequently, for each surviving word $w$,
\begin{equation}
 |G(w;1)|\leq C_w
 \prod_{j\,\mathrm{moving}}\frac2{|\alpha_j(q)|},
 \qquad C_w<\infty,
\end{equation}
where $C_w$ is independent of $q$. The factors $1/a$ and $1/c$
are $O(\delta^2)$, while $1/r_\pm=O(\delta)$. The bound therefore
tends to zero, and so does the finite shuffle combination defining
$g_w$. This proves the moving-letter boundary value through every
weight required here. It applies to each special function before
multiplication by possibly singular rational energy coefficients.

Differentiating the endpoint GPLs by the total-differential recursion
of Section~\ref{sec:comp-virtual-energy-basis} now gives logarithmic
differentials in a fixed alphabet with the explicitly fixed constants
above. Add the complete Hermite boundary difference, including its
Laurent--Taylor cross terms, before returning to the $X$ kernels.
The result contains eighteen words with rational coefficients in $X$
and $Y$. Exact substitution of the kernel maps reproduces every
coefficient of the fixed-alphabet expression.''',
 r'\paragraph{Collection of the functions.} The shuffle identity gives',
 r'\begin{equation}G(0,r,-1;X)=L_0G(r,-1;X)-G(r,0,-1;X)-G(r,-1,0;X).\end{equation}',
 r'The five root words occur with coefficients proportional to $1,2,1,2,1$, in the order $(-4,r,-1)$, $(0,r,-1)$, $(r,-1,-1)$, $(r,-1,0)$, $(r,0,-1)$. Applying this identity reduces their sum to $K[\mathcal U+2L_0G(r,-1;X)]$. For the rational letters, use',
 r'\begin{align}G(-1,0;X)&=L_0L_1+\operatorname{Li}_2(-X),\\',
 r'G(-2,0;X)&=L_0L_2+\operatorname{Li}_2(-X/2),\\',
 r'G(-j,-j;X)&=\tfrac12\ln^2(1+X/j),\qquad j=1,2.\end{align}',
 r'These are shuffle identities and the defining dilogarithm integral, so no continuation across a cut is needed. All kernels and logarithms in the final expression are real for $X>0$.',

]
write_artifact(folder / 'tex/uw_compact_formula.tex' , '\n'.join(out) + '\n')
write_artifact(folder / 'uw_compact_coefficients.json' , json.dumps({
    'source': 'own UW fixed-letter integral and exact shuffle identities',
    'coefficients': {k: str(v) for k, v in coeff.items()},
    'normalization': 'one finite UW sector before 1/(8*s*H)',
    'shuffle_reconstruction': 'exact zero for all eighteen word coefficients',
}, indent=2) + '\n')
