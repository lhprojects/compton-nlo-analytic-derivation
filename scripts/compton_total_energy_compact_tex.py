"""Render the product form without expanding its elementary log polynomial."""
from compton_artifacts import read_artifact, write_artifact
import json
import re
from pathlib import Path
import sympy as S
from compton_total_regular import coefficients, certify_real_pair
folder=Path(__file__).resolve().parents[1]
X,Y,L,M=S.symbols('X Y L M')
data=json.loads(read_artifact(folder/'total_energy_compact.json'))
terms=dict(coefficients())
certify_real_pair()
poly=S.Poly(S.expand(terms.pop(())),L,M)
elementary=[(powers,min([S.factor(c),S.apart(c,X)],key=lambda a:(S.count_ops(a),len(str(a))))) for powers,c in poly.terms()]
assert S.cancel(terms[('0','-1','0')]-terms[('0','1','0')])==0
terms.pop(('0','1','0'))
out=[r'\subsection{A product form of the total}',
 r'\label{sec:comp-total-energy-product}',
 r'Products of elementary logarithms can make the collected result shorter without introducing additional polylogarithms. Define $L=\ln(1+X)$ and $M=\ln X$. The two GPLs with the same coefficient combine into',
 r'\begin{equation}',
 r'\mathcal P(X)=G(0,-1,0;X)+G(0,1,0;X)=\tfrac12\left[\operatorname{Re}\operatorname{Li}_3(X^2)-M\operatorname{Re}\operatorname{Li}_2(X^2)\right].',
 r'\end{equation}',
 r'For the derivation, first take $0<X<1$. The shuffle relation $G(0)G(0,a)=2G(0,0,a)+G(0,a,0)$, with $G(0,a;X)=-\operatorname{Li}_2(X/a)$ and $G(0,0,a;X)=-\operatorname{Li}_3(X/a)$, gives',
 r'\begin{equation}',
 r' G(0,a,0;X)=2\operatorname{Li}_3(X/a)-M\operatorname{Li}_2(X/a),\qquad a=\pm1.',
 r'\end{equation}',
 r'Summing the two signs and using $\operatorname{Li}_n(z)+\operatorname{Li}_n(-z)=2^{1-n}\operatorname{Li}_n(z^2)$ gives the stated result. Both shuffle multiplicities and this duplication factor are needed for its factor $1/2$. For $X>1$, continue all functions on the same edge of the principal cut. On either edge, the imaginary parts of $\operatorname{Li}_2(X)$ and $\operatorname{Li}_3(X)$ are respectively $\pm\pi\ln X$ and $\pm\pi\ln^2X/2$, with the same sign. These follow by integrating their derivatives from $X=1$, and they cancel in $2\operatorname{Li}_3(X)-M\operatorname{Li}_2(X)$. The analogous cancellation holds for the paired expression at argument $X^2$. Thus both edges give the same real function and the continuation adds no branch constant. The following equivalent form uses only real polylogarithms with arguments in $(0,1]$:',
 r'\begin{equation}',
 r'\mathcal P(X)=\begin{cases}\tfrac12[\operatorname{Li}_3(X^2)-M\operatorname{Li}_2(X^2)],&0<X\leq1,\\ \tfrac12[\operatorname{Li}_3(X^{-2})+M\operatorname{Li}_2(X^{-2})]+\tfrac13M^3+\tfrac16\pi^2M,&X\geq1.\end{cases}',
 r'\end{equation}',
 r'To derive the second line, put $t=\ln X>0$, $D(t)=\operatorname{Li}_2(e^{-2t})$ and $T(t)=\operatorname{Li}_3(e^{-2t})$. Their derivatives are $D^{\prime}(t)=2\ln(1-e^{-2t})$ and $T^{\prime}(t)=-2D(t)$. The identity $\ln(e^{2t}-1)=2t+\ln(1-e^{-2t})$ then gives',
 r'\begin{align}',
 r'\operatorname{Re}\operatorname{Li}_2(e^{2t})&=\frac{\pi^2}{3}-2t^2-D(t),\\',
 r'\operatorname{Re}\operatorname{Li}_3(e^{2t})&=T(t)-\frac43t^3+\frac{2\pi^2}{3}t.',
 r'\end{align}',
 r'Indeed, the derivatives of the two sides are respectively $-2\ln(e^{2t}-1)$ and $2\operatorname{Re}\operatorname{Li}_2(e^{2t})$. Write $\zeta(n)=\sum_{k=1}^{\infty}k^{-n}$ for the Riemann zeta function at integers $n>1$. Their values at $t=0$ agree because $\operatorname{Li}_2(1)=\pi^2/6$ and $\operatorname{Li}_3(1)=\zeta(3)$. These derivative and boundary conditions fix both identities and give $\mathcal P(1)=\zeta(3)/2$.',
 r'The sole apparent rational pole at $X=1$ can also be removed before evaluation. Define the continuous elementary function',
 r'\begin{equation}',
 r'R_{\log}(X)=\begin{cases}\ln X/(X-1),&X\ne1,\\1,&X=1,\end{cases}\qquad R_{\log}(1+h)=1-\frac h2+\frac{h^2}{3}+O(h^3).',
 r'\end{equation}',
 r'The series follows by dividing the Taylor series of $\ln(1+h)$ by $h$, which converges for $|h|<1$. Extracting the common power $X^{-3}$ gives the exact total coefficient as',
 r'\begin{equation}',r'C_1(X)=\frac{1}{X^3}\left[R_{\log}(X)+\mathcal E(X,L,M)+\sum_{j=1}^{13}b_j(X,Y,L,M)\mathcal F_j(X)\right].',r'\label{eq:comp-total-energy-product}\end{equation}',
 r'The elementary part is',r'\begin{equation}',
 r'\mathcal E='+ '+'.join(f'e_{{{i}}}'+('L' if p[0]==1 else f'L^{{{p[0]}}}' if p[0] else '')+('M' if p[1]==1 else f'M^{{{p[1]}}}' if p[1] else '') for i,(p,c) in enumerate(elementary))+'.',r'\end{equation}']
def rhs(value):
    # Preserve exact partial fractions when they are shorter than a single fraction.
    poly=S.Poly(S.expand(value),L,M,Y)
    result=0
    for (l,m,y),v in poly.terms():
        best=min([S.factor(v),S.apart(v,X)],key=lambda a:(S.count_ops(a),len(str(a))))
        result+=best*L**l*M**m*Y**y
    return S.latex(result)
for i,(powers,c) in enumerate(elementary):
    out += [r'\begin{equation}',f'e_{{{i}}}='+rhs(c)+'.',r'\end{equation}']
out += [r'The thirteen remaining functions and their coefficients are']
for i,(word,c) in enumerate(terms.items(),1):
    fun=r'\mathcal P(X)' if word==('0','-1','0') else 'G('+','.join(word)+';X)'
    if c.has(L) or c.has(M):
        # Cancel each log coefficient separately before recombining. A common
        # denominator obscures the exact cancellation in the L coefficient.
        pieces=S.Poly(S.expand(c),L,M).terms()
        constant=S.factor(c.subs({L:0,M:0}))
        logs=[(powers,S.factor(value)) for powers,value in pieces if powers!=(0,0)]
        assert S.cancel(constant+sum(value*L**powers[0]*M**powers[1] for powers,value in logs)-c)==0
        logpart=sum(value*L**powers[0]*M**powers[1] for powers,value in logs)
        out += [r'\begin{equation}',f'\mathcal F_{{{i}}}={fun}.',r'\end{equation}',
                r'\begin{equation}',f'b_{{{i}}}='+rhs(constant)+S.latex(logpart)+'.',r'\end{equation}']
    elif word==('-2','-1'):
        out += [r'\begin{equation}',f'\mathcal F_{{{i}}}={fun}.',r'\end{equation}',
                r'\begin{equation}',f'b_{{{i}}}='+rhs(c)+'.',r'\end{equation}']
    else:
        out += [r'\begin{equation}',f'\mathcal F_{{{i}}}={fun},\qquad b_{{{i}}}='+rhs(c)+'.',r'\end{equation}']
out += [r'This form follows from the complete coefficient table using only the identities',
 r'\begin{align}',
 r'G(-1,0)&=LM-G(0,-1),\\',
 r'G(-1,0,-1)&=L\,G(0,-1)-2G(0,-1,-1),\\',
 r'G(-1,-1,0)&=\tfrac12ML^2-L\,G(0,-1)+G(0,-1,-1),\\',
 r'G(0,r,-1)&=M\,G(r,-1)-G(r,0,-1)-G(r,-1,0).',
 r'\end{align}',
 r'The common endpoint is $X$. In addition, repeated letters $-1$ give powers of $L$ divided by factorials. Expanding every product back by the shuffle identity reproduces all twenty-three original coefficients exactly. The two words combined into $\mathcal P$ have the same coefficient $c_{15}=c_{16}$, so this last grouping replaces their sum without another factor of two.',
 r'The extraction of $R_{\log}$ is fixed by the coefficient $c_2$ of $M$ in the collected table. Direct polynomial division gives',
 r'\begin{equation}',
 r' X^3c_2(X,Y)=\frac{1}{X-1}+e_5(X),\qquad X^3c_2(X,Y)M=R_{\log}(X)+e_5(X)M.',
 r'\end{equation}',
 r'The second equality at $X=1$ means its continuous limit. All other coefficients are multiplied by $X^3$ before the common factor $X^{-3}$ is placed outside the bracket. Thus this extraction changes neither the normalization of $C_1$ nor the coefficient of its removable term. Every displayed rational coefficient is regular at $X=1$, as follows by exact polynomial cancellation. The function $G(1,0;X)=\zeta(2)-\operatorname{Li}_2(1-X)$ is also explicitly real and regular for every $X>0$. Its derivative is $\ln X/(X-1)$ and its limit as $X\to0^+$ is zero, which fixes the identity without numerical continuation.'
]
# Share the complete result between the opening overview and its derivation.
cut=out.index(r'This form follows from the complete coefficient table using only the identities')
body='\n'.join(out[2:cut]).replace(r'\label{eq:comp-total-energy-product}',r'\label{\comptonproductlabel}')
write_artifact(folder/'tex/total_energy_product_body.tex' , body+'\n')
wrapper=out[:2]+[r'\begingroup',r'\def\comptonproductlabel{eq:comp-total-energy-product}',
 r'\input{total_energy_product_body}',r'\endgroup']+out[cut:]
write_artifact(folder/'tex/total_energy_compact_formula.tex' , '\n'.join(wrapper)+'\n')

# README and both document occurrences use this same coefficient rendering.
# Keep complete formulas here rather than a manually copied summary.
begin='## The answer'
end='## Simplification of our result'
readme_path=folder/'README.md'
readme=readme_path.read_text()
if readme.count(begin)!=1 or readme.count(end)!=1:
    raise ValueError('README must contain exactly one answer heading and one following simplification heading')
prefix,marked=readme.split(begin,1)
_,suffix=marked.split(end,1)

def displayed_source(path, label):
    source=path.read_text()
    for equation in re.findall(r'\\begin\{equation\}(.*?)\\end\{equation\}',source,re.S):
        marker=r'\label{'+label+'}'
        if marker in equation:
            return equation.replace(marker,'').strip()
    raise ValueError(f'Missing source equation {label}')

total_source=folder/'tex/total_simplification.tex'
total_equation=displayed_source(total_source,'eq:comp-simplified-total')
born_equation=displayed_source(total_source,'eq:comp-simplified-born')
basis_source=(folder/'tex/virtual_energy_basis.tex').read_text()
kernel_equations=re.search(r'\\begin\{align\}(.*?)\\end\{align\}',basis_source,re.S).group(1).strip()
kernel_equations=kernel_equations.replace(r'\dd',r'\mathrm{d}')
markdown=[
    'For a photon of energy $\\omega$ in the electron rest frame, define '
    '$X=2\\omega/m=s_{\\mathrm{phys}}/m^2-1>0$ and $Y=\\sqrt{X/(X+4)}$. '
    'Here $m$ is the electron pole mass and $\\alpha_{\\mathrm T}$ is the Thomson coupling. '
    'The inclusive total through next-to-leading order is',
    '$$\n'+total_equation+'\n$$',
    'With $r_\\omega=X/2$, the Born coefficient is',
    '$$\n'+born_equation+'\n$$',
    'For an ordered list of letters $w$, with $\\varnothing$ denoting the empty list, '
    'the iterated integrals use the following kernels:',
    '$$\n'+r'\begin{aligned}'+'\n'+kernel_equations+'\n'+r'\end{aligned}'+'\n$$',
    'The symbol $r$ labels the positive square-root kernel. Pure-zero words have '
    '$G(0,\\ldots,0;X)=\\ln^nX/n!$ for $n$ zeros. The combinations containing '
    'the letter $1$ have removable singularities, with their continuous real values fixed in the derivation.',
    'The ordinary polylogarithm is $\\operatorname{Li}_n(z)=\\sum_{k=1}^{\\infty}z^k/k^n$ '
    'for $|z|<1$, continued on its principal branch elsewhere.',
    'Define $L=\\ln(1+X)$, $M=\\ln X$ and $\\zeta(3)=\\sum_{k=1}^{\\infty}k^{-3}$.',
    'The paired function and its real representations are given below; the latter use polylogarithm arguments in $(0,1]$:',
]
starts=(r'\mathcal P(X)',r'R_{\log}(X)',r'C_1(X)',r'\mathcal E=',
        'e_{',r'\mathcal F_{','b_{')
for equation in re.findall(r'\\begin\{equation\}(.*?)\\end\{equation\}',body,re.S):
    equation=equation.strip()
    if not equation.startswith(starts):
        continue
    equation=equation.replace(r'\label{\comptonproductlabel}','').strip()
    # Keep README piecewise definitions in separate displays: some Markdown
    # renderers misparse cases row separators before passing them to math.
    if equation.startswith(r'\mathcal P(X)=\begin{cases}'):
        branches=equation.split(r'\begin{cases}',1)[1].split(r'\end{cases}',1)[0]
        for branch in branches.split(r'\\'):
            value,domain=branch.strip().split('&',1)
            markdown.append('$$\n'+r'\mathcal P(X)='+value.rstrip(',')
                            +r',\qquad '+domain.rstrip(',.')+'.\n$$')
        continue
    if equation.startswith(r'R_{\log}(X)=\begin{cases}'):
        cases,tail=equation.split(r'\end{cases}',1)
        branches=cases.split(r'\begin{cases}',1)[1].split(r'\\')
        first,domain=branches[0].split('&',1)
        value,at=branches[1].split('&',1)
        markdown.append('$$\n'+r'R_{\log}(X)='+first.rstrip(',')
                        +r',\qquad '+domain.rstrip(',')+'.\n$$')
        markdown.append('$$\n'+r'R_{\log}(1)='+value.rstrip(',')+','+tail+'\n$$')
        continue
    if equation.startswith(r'\mathcal E='):
        markdown.append('The elementary logarithmic polynomial and its coefficients are:')
    elif equation.startswith(r'\mathcal F_{1}='):
        markdown.append('The thirteen grouped functions and their coefficients are:')
    markdown.append('$$\n'+equation+'\n$$')
markdown.append('The shuffle identities and the analytic derivative and boundary arguments '
                'in [the derivation](tex/total_energy_compact_formula.tex) establish this '
                'form, including its real continuation through $X=1$. Expanding the products '
                'recovers the original twenty-three coefficients exactly.')
# Fenced math prevents Markdown from interpreting TeX underscores or row breaks.
# Use basic supported commands for polylogarithms and real parts.
rendered='\n\n'.join(markdown)
rendered=rendered.replace(r'\operatorname{Li}',r'\mathrm{Li}')
rendered=rendered.replace(r'\operatorname{Re}',r'\Re')
rendered=re.sub(r'\$\$\n(.*?)\n\$\$',lambda m: '```math\n'+m.group(1)+'\n```',rendered,flags=re.S)
readme_path.write_text(prefix+begin+'\n\n'+rendered+'\n\n'+end+suffix)
