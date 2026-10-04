"""Write the complete compact matched real formula with every coefficient."""
from compton_artifacts import read_artifact, write_artifact
import json
from pathlib import Path
import sympy as S
folder=Path(__file__).resolve().parents[1]
s,Ls,Lh,Lp=S.symbols('s Ls Lh Lp');Li=S.Function('L2');Li3=S.Function('Li3')
d=json.loads(read_artifact(folder/'real_four_matched_short.json'));atoms=list(map(S.sympify,d['atoms']))
e=sum(S.sympify(t['coefficient'])*S.prod(a**p for a,p in zip(atoms,t['powers'])) for t in d['terms'])
L,J,M,D,V,T,P,c,Z=S.symbols('L J M D V T P c Z')
e=e.subs({Li(1/s):D,Li(2/(s+1)):V-c*J,Li3(1/s):Z-T,S.zeta(3):Z,S.pi**2:P,S.log(2):c,Ls:L,Lh:M,Lp:J+c},simultaneous=True)
e=S.expand(e);poly=S.Poly(e,L,J,M,D,V,T,P,c,Z)
terms={k:S.factor(S.cancel(v)) for k,v in poly.terms() if S.cancel(v)!=0}
assert all(k[-1]==k[-2]==0 for k in terms)
def coefficient(mon):
    target=S.Poly(mon,L,J,M,D,V,T,P,c,Z).monoms()[0]
    return terms.pop(target,S.S.Zero)
coeff={**{f'A_{i}':coefficient(L**i) for i in range(4)},
       'B_0':coefficient(P),'B_1':coefficient(P*L),
       'C_0':coefficient(D),'C_1':coefficient(D*L),
       'K':coefficient(T),'N':coefficient(V),'E':coefficient(J*J),'F':coefficient(J),
       'H_0':coefficient(M),'H_1':coefficient(M*L),'H_2':coefficient(M*L*L),'H_3':coefficient(M*J)}
assert not terms,terms
out=[r'\subsection{The matched real contribution}',
 r'The four real sectors and their finite matching shifts can now be combined into a formula with two dilogarithms and one trilogarithm. Define',
 r'\begin{equation}',r' \mathcal R_4(s)=F_{\mathrm{rat}}^z+F_W^z+F_{\log}^z+F_D^z',
 r'                    +\Delta_{\mathrm{rat}}+\Delta_{\log}+\Delta_D.',r'\end{equation}',
 r'Here $s>1$, $H=s-1$ and $r=\sqrt{s}$. Use the real logarithms',
 r'\begin{equation}',r' L=\ln s,\qquad M=\ln(s-1),\qquad J=\ln\frac{s+1}{2},',r'\end{equation}',
 r'and the following three special-function combinations:',r'\begin{align}',
 r' D(s)&=\operatorname{Li}_2(1/s),\\',
 r' V(s)&=\operatorname{Li}_2\!\left(\frac2{s+1}\right)+\ln2\,J,\\',
 r' T(s)&=\zeta(3)-\operatorname{Li}_3(1/s).',r'\end{align}',
 r'Both dilogarithm arguments and the trilogarithm argument lie strictly between zero and one. The complete matched expression is',r'\begin{equation}',r'\begin{aligned}',
 r' \mathcal R_4(s)={}&A_0+A_1L+A_2L^2+A_3L^3',r'                  +\pi^2(B_0+B_1L)\\',
 r' &+(C_0+C_1L)D+KT+NV+EJ^2+FJ\\',r' &+M(H_0+H_1L+H_2L^2+H_3J).',r'\end{aligned}',
 r'\label{eq:comp-real-four-matched-compact}',r'\end{equation}',
 r'Every coefficient is a rational function of $s$. Explicitly,']
def latex(e):return S.latex(e).replace(r'\left(s - 1\right)','H')
for name,value in coeff.items():
    n,den=S.fraction(value);nterms=S.Add.make_args(S.expand(n));nterms=sorted(nterms,key=lambda t:S.degree(t,s),reverse=True)
    if len(nterms)>5:
        pieces=[]
        for i in range(0,len(nterms),4):
            piece=latex(sum(nterms[i:i+4]))
            if i and not piece.startswith('-'):piece='+'+piece
            pieces.append('&'+piece)
        rhs=r'\frac{\begin{aligned}'+r'\\'.join(pieces)+r'\end{aligned}}{'+latex(den)+'}'
    else:rhs=latex(value)
    out.extend([r'\begin{equation}',name+'='+rhs+'.',r'\end{equation}'])
out.extend([
 r'This is the unnormalized real finite contribution. Its relation to the remainder in \eqref{eq:comp-simplified-master} is',
 r'\begin{equation}',r' \mathcal R_{\mathrm{rem}}(X)=\mathcal R_4(1+X)+2F_{UW}^{\eta}(X).',r'\end{equation}',
 r'Thus \eqref{eq:comp-real-four-matched-compact} replaces the four separate endpoint evaluations and their three matching shifts. The UV term is already included with the virtual contribution in that master formula. The remaining $UW$ and finite soft terms are separate.',
 r'\paragraph{Reduction of the function basis.} After the endpoint identities above, real reflection and inversion relations put the dilogarithm arguments in $(0,1)$. Two five-term identities reduce their number further:',
 r'\begin{align}',r' \operatorname{Li}_2\!\left(\frac2{r+1}\right)+\operatorname{Li}_2\!\left(\frac{r+1}{2r}\right)',
 r' &=\operatorname{Li}_2(1/r)+\operatorname{Li}_2\!\left(\frac1{r+1}\right)+\frac{\pi^2}{12}-\frac12\ln^22+\ln2\ln\frac{r+1}{r}.',r'\end{align}',
 r'The same identity with $r$ replaced by $s$ reduces the half-argument function in the matching shift. For the second identity, define $\mathcal L(z)=\operatorname{Li}_2(z)+\tfrac12\ln z\ln(1-z)$ for $0<z<1$. Then',
 r'\begin{equation}',r' \mathcal L\!\left(\frac{r+1}{r^2+1}\right)-\mathcal L\!\left(\frac{r-1}{r(r+1)}\right)',
 r' =\mathcal L(1/r)+\mathcal L\!\left(\frac1{r^2+1}\right)-\mathcal L\!\left(\frac1{r+1}\right).',r'\end{equation}',
 r'Both identities follow by differentiation using $\dd\operatorname{Li}_2(z)/\dd z=-\ln(1-z)/z$, with the integration constant fixed as $r\to\infty$. The remaining square-root arguments occur in a common-coefficient difference. Landen transformation and duplication give',
 r'\begin{equation}',r' \operatorname{Li}_2\!\left(\frac1{r+1}\right)-\operatorname{Li}_2(1/r)',
 r' =-\frac12\operatorname{Li}_2(1/r^2)-\frac12\ln^2\!\left(\frac{r+1}{r}\right).',r'\end{equation}',
 r'After this substitution, all coefficients of $\ln(r+1)$ cancel exactly when $\ln(r-1)=\ln(s-1)-\ln(r+1)$ is used. The result therefore has no square-root argument. Adding the matching shift also cancels the entire coefficient of $\operatorname{Li}_2(1/(s+1))$.',
 r'The matching shifts are reduced with the explicit identities',
 r'\begin{align}',
 r' \operatorname{Li}_2(1-s)&=D-LM+\tfrac12L^2-\frac{\pi^2}{6},\\',
 r' \operatorname{Li}_2\!\left(-\frac{(s-1)^2}{4s}\right)',
 r' &=4\operatorname{Li}_2\!\left(\frac2{s+1}\right)-2D-2\operatorname{Li}_2\!\left(\frac1{s+1}\right)\\',
 r' &\quad+J^2+2J(L+\ln2)-\tfrac32L^2+2M(L-2J)-\frac{\pi^2}{6}-\ln^22.',
 r'\end{align}',
 r'For the second line, put $v=(s-1)/(s+1)$, apply Landen transformation to the negative argument, and then use $\operatorname{Li}_2(v^2)=2\operatorname{Li}_2(v)+2\operatorname{Li}_2(-v)$. Reflection and the preceding half-argument identity give the displayed expression. All intermediate arguments used in these identities are real and their logarithms have the signs fixed by $s>1$.',
 r'Equivalently, denote the two right-hand sides by $Q_1(s)$ and $Q_2(s)$. Differentiation and the positive-factor logarithm identities for $s>1$ give',
 r'\begin{equation}',
 r' Q_1\prime(s)=-\frac{\ln s}{s-1},\qquad',
 r' Q_2\prime(s)=-\frac{2(s+1)}{s(s-1)}\ln\frac{s+1}{2\sqrt s}.',
 r'\end{equation}',
 r'These are precisely the derivatives of the two left-hand sides. At $s\to1^+$ both sides vanish. Here the coefficient of $M=\ln(s-1)$ vanishes at least linearly, and $\operatorname{Li}_2(1/2)=\pi^2/12-\ln^22/2$ fixes the remaining constant. Equality of derivatives and this boundary value establish the identities on the connected domain $s>1$.',
])
write_artifact(folder/'tex/real_four_matched_formula.tex' , '\n'.join(out)+'\n')
write_artifact(folder/'real_four_matched_coefficients.json' , json.dumps({k:str(v) for k,v in coeff.items()},indent=2)+'\n')
print('Wrote compact formula with',len(coeff),'explicit rational coefficients')
