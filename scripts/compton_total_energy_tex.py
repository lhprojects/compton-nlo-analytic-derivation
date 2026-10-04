"""Write the full independently collected energy formula with all coefficients."""
from compton_artifacts import read_artifact, write_artifact
import json
from pathlib import Path
import sympy as S

folder=Path(__file__).resolve().parents[1]
data=json.loads(read_artifact(folder/'total_energy_collected.json'))
X,Y=S.symbols('X Y')
out=[r'\section{The collected total in one function basis}',
 r'\label{sec:comp-total-energy-collected}',
 r'All contributions can now be combined in the same energy variable. This section gives the complete coefficient of $\alpha_{\mathrm T}^3/m^2$, including the virtual correction and the integrated two-photon final state. It is obtained from the independently evaluated sectors above. Define $X=s_{\mathrm{phys}}/m^2-1>0$ and $Y=\sqrt{X/(X+4)}$. Use the GPL kernels and real branches in Section~\ref{sec:comp-virtual-energy-basis}, including the removable combinations containing the letter $1$. The complete result is',
 r'\begin{equation}',
 r'\sigma_{\mathrm{inc}}^{\mathrm{NLO}}=\frac{\pi\alpha_{\mathrm T}^2}{m^2}B(X)+\frac{\alpha_{\mathrm T}^3}{m^2}C_1(X)+\frac1{m^2}O(\alpha_{\mathrm T}^4),',
 r'\end{equation}',
 r'where $B$ is given in \eqref{eq:comp-simplified-born} and',
 r'\begin{equation}\boxed{C_1(X)=\sum_{j=0}^{22}c_j(X,Y)F_j(X).}\label{eq:comp-total-energy-collected}\end{equation}',
 r'The functions and their full coefficients follow. No unevaluated phase-space or Feynman-parameter integral occurs in this expression.']
for i,row in enumerate(data['terms']):
    word=row['word'];func='1' if not word else 'G('+','.join(word)+';X)'
    value=S.sympify(row['coefficient']);num,den=S.fraction(value)
    terms=sorted(S.Add.make_args(S.expand(num)),key=lambda t:S.degree(t,X),reverse=True)
    if len(terms)>5:
        pieces=[]
        for n in range(0,len(terms),4):
            p=S.latex(sum(terms[n:n+4]))
            if n and not p.startswith('-'):p='+'+p
            pieces.append('&'+p)
        rhs=r'\frac{\begin{aligned}'+r'\\'.join(pieces)+r'\end{aligned}}{'+S.latex(den)+'}'
    else:rhs=S.latex(value)
    out += [r'\begin{equation}',f' F_{{{i}}}={func},\qquad c_{{{i}}}='+rhs+'.',r'\end{equation}']
out += [
 r'\paragraph{Collection of the real terms.} Put $L=\ln(1+X)$, $M=\ln X$, and $J=\ln(1+X/2)$. The classical functions in the four-sector real expression satisfy',
 r'\begin{align}',
 r'\operatorname{Li}_2\!\left(\frac1{1+X}\right)&=\zeta(2)+G(-1,0;X)-G(-1,-1;X),\\',
 r'\operatorname{Li}_2\!\left(\frac2{2+X}\right)+\ln2\,J&=\zeta(2)+G(-2,0;X)-G(-2,-2;X),\\',
 r'\zeta(3)-\operatorname{Li}_3\!\left(\frac1{1+X}\right)&=\zeta(2)G(-1;X)+G(-1,-1,0;X)-G(-1,-1,-1;X).',
 r'\end{align}',
 r'Each identity follows by differentiating with respect to $X$. The constants are fixed as $X\to0^+$, where the dilogarithms and trilogarithm take their values at one. Every nonconstant GPL on these right-hand sides vanishes in this limit. In particular, the logarithmic lower endpoint in $G(-1,0;X)$ and $G(-2,0;X)$ is integrable because $\int_0^X|\ln t|\,\mathrm{d}t=O(X[1+|\ln X|])$. The remaining products of $L$, $M$, and $J$ are expanded with the shuffle identity, counting each interleaving even when repeated letters give the same word. All classical arguments above lie in $(0,1)$ for $X>0$, so this conversion introduces no continuation across a cut.',
 r'\paragraph{Assembly.} Substitute these identities, the completed $UW$ expression, and \eqref{eq:comp-ck-energy-compact} into \eqref{eq:comp-simplified-master}. The finite virtual coefficient is its constant term in $\ell_\lambda=\ln(\lambda_{\mathrm{phys}}/m)$ after the previously proved infrared cancellation. Taking this coefficient does not set the photon mass to $m$. The finite conversion term remains $-A(1+X)[\ln X+\tfrac12\ln(1+X)]$. Expanding this product gives the last coefficients needed in the sum above.',
 r'There are five contributions to this collection: $\widehat C_{V+\mathrm{UV}}$, $F_{UW}^{\eta}/[4X(1+X)]$, $\mathcal R_4(1+X)/[8X(1+X)]$, $-\mathcal C_K(1+X)/(1+X)$, and the finite conversion term. The factor $1/4$ for $UW$ follows from $\mathcal R_{\mathrm{rem}}=\mathcal R_4+2F_{UW}^{\eta}$, so its Bose-related second copy is already included. Expand each contribution in the same words and add their rational coefficients. This is an algebraic reconstruction and does not require a claim that the functions are linearly independent. The infrared coefficient used in the last contribution is \eqref{eq:comp-ir-short}, with $L^2=2G(-1,-1;X)$ and $\mathcal G(X)=G(r,-1;X)$.',
 r'The coefficients of five words vanish identically in the complete sum. In particular the coefficient of $G(-1,r,-1;X)$ from the finite soft term cancels its coefficient from the infrared conversion term. All coefficients proportional to $\pi^2$ outside the GPLs also cancel. The remaining coefficients are the rational functions listed above.',
 r'\paragraph{The apparent point $X=1$.} The factor $1/(X-1)$ multiplies $G(0;X)=\ln X$, so it has a removable limit. The product involving this coefficient tends to $+1$. The two GPLs containing the rational letter $+1$ are also real and finite:',
 r'\begin{align}',
 r'G(1,0;X)&=\ln X\ln|1-X|+\operatorname{Re}\operatorname{Li}_2(X),\\',
 r'G(0,1,0;X)&=2\operatorname{Re}\operatorname{Li}_3(X)-\ln X\operatorname{Re}\operatorname{Li}_2(X).',
 r'\end{align}',
 r'Their derivatives are $\ln X/(X-1)$ and $G(1,0;X)/X$. At $X=1$ their values are $\zeta(2)$ and $2\zeta(3)$. Thus the apparent inner pole does not require a contour choice for the physical answer.',

]
write_artifact(folder/'tex/total_energy_formula.tex' , '\n'.join(out)+'\n')
