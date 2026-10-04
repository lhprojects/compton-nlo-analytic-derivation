"""Checks the published inclusive expression against independent evaluations.

The amplitude tests are separate in check_compton_amplitudes.py. This file
checks integration accuracy, elementary reductions and Denner--Dittmaier
Table 3 (hep-ph/9805443), not an independent numerical loop integration.
"""
from compton_artifacts import write_artifact
import json
import math
from pathlib import Path
import mpmath as mp
import compton_nlo as c


def main():
    g=c.iterated_integrals(2e4)
    tight=c.iterated_integrals(2e4,rtol=3e-14,atol=2e-15,start=-40)
    worst=0.
    for x in [.07,.5,1.,10.,2e4]:
        with mp.workdps(40):
            X=mp.mpf(x);z=4/(mp.sqrt(X+4)+mp.sqrt(X))**2
            checks={(-2,-1):mp.polylog(2,-X-1)+mp.log(X+1)*mp.log(X+2)+mp.pi**2/12,
                (-1,-2,-1):2*mp.polylog(3,-X-1)-mp.polylog(2,-X-1)*mp.log(X+1)+mp.pi**2/12*mp.log(X+1)+3*mp.zeta(3)/2,
                ('r',-1):mp.polylog(2,-z**3)/3-mp.polylog(2,-z)+mp.log(z)**2/2-mp.pi**2/18,
                ('r','r',-1):mp.polylog(3,-z)-mp.polylog(3,-z**3)/9-mp.log(z)**3/6+mp.pi**2/18*mp.log(z)+2*mp.zeta(3)/3}
            for w,v in checks.items():worst=max(worst,abs(g(x)[w]-float(v))/max(1.,abs(float(v))))
    convergence=max(abs(c.correction(w,g)-c.correction(w,tight))/abs(c.correction(w,tight)) for w in [4.4,5.,5.47])
    olda,oldm=c.ALPHA,c.MASS;c.ALPHA=1/137.0359895;c.MASS=.00051099906
    table=[]
    for E,reference in [(4.,.0023),(6.,.0045),(8.,.0071)]:
        omega=1.165e-9*(E+math.sqrt(E*E-c.MASS*c.MASS))/c.MASS
        value=100*c.correction(omega,tight)/c.born(omega)
        assert abs(value-reference)<.00005
        table.append(dict(electron_GeV=E,laser_eV=1.165,delta_percent=value,published_percent=reference))
    c.ALPHA,c.MASS=olda,oldm
    high=c.iterated_integrals(1/c.MASS**2)
    omega=(1-c.MASS**2)/(2*c.MASS)
    ecm=100*c.correction(omega,high)/c.born(omega)
    assert abs(ecm-8.2)<.05 and worst<2e-11 and convergence<1e-7
    report=dict(iterated_integral_scaled_error=worst,primex_convergence_relative=convergence,
                denner_table3=table,delta_percent_at_ecm_1_GeV=ecm)
    print(json.dumps(report,indent=2))
    write_artifact(Path(__file__).resolve().parents[1]/'inclusive_checks.json' , json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
