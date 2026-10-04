"""Three-way Compton total comparison, with 50-digit independent analytic values.

Reference values only enter the difference column. No missing analytic or
numerical result is replaced by a reference value. Errors are arithmetic and
quadrature estimates, not experimental errors or confidence intervals.
"""
from compton_artifacts import read_artifact, write_artifact
from pathlib import Path
import json
import mpmath as mp
mp.mp.dps=160
BASE=Path(__file__).resolve().parents[1]
ENERGIES=['1.0002','1.02','1.6','10','100','2001','200001']

def number(value,digits=3):
    s=mp.nstr(mp.mpf(value),digits,strip_zeros=False)
    if 'e' in s:
        mantissa,exponent=s.split('e')
        return mantissa+r'\times10^{'+str(int(exponent))+'}'
    return s


def main():
    analytic={}
    for path in (BASE/'reports').glob('independent_closed_s*.json'):
        d=json.loads(read_artifact(path))
        if d.get('complete'):
            analytic.setdefault(mp.mpf(d['s']),[]).append((d,path.name))
    virtual={mp.mpf(d['s']):d for d in json.loads(read_artifact(BASE/'compton_soft_uv_numeric_virtual.json'))}
    real={}
    for name in ['compton_real_analytic_numerical_audit.json','compton_real_analytic_numerical_wide_audit.json']:
        data=json.loads(read_artifact(BASE/name))
        real.update({mp.mpf(k):v for k,v in data['energies'].items()})
    reference={mp.mpf(d['s']):d for d in json.loads(read_artifact(BASE/'compton_box_analytic_reference_check.json'))['points']}
    fresh=json.loads(read_artifact(BASE/'total_compact_table_values.json'))
    assert fresh['complete'], 'fresh compact evaluation incomplete'
    compact={mp.mpf(row['s']):row for row in fresh['points']}
    lines=[];checks=[]
    for energy in ENERGIES:
        ss=mp.mpf(energy);runs=sorted(analytic[ss],key=lambda x:x[0]['dps'])
        assert len(runs)>=2,('second precision run missing',energy)
        d,path=runs[-1];own=mp.mpf(d['components']['total'])
        printed=mp.nstr(own,50,strip_zeros=False)
        compatible=[old for old,_ in runs[:-1]
                    if mp.nstr(mp.mpf(old['components']['total']),50,strip_zeros=False)==printed]
        assert compatible,('50 digits not stable',energy)
        current=compact[ss]
        assert current['analytic_50']==printed, ('compact and original differ',energy)
        own=mp.mpf(current['value'])
        vs=virtual[ss];rr=real[ss];ref=reference[ss]
        assert rr['all_requested_tolerances_met'],('real quadrature not converged',energy)
        num=mp.mpf(vs['values'][-1])+mp.mpf(str(rr['value']))
        en=mp.mpf(vs['estimated_absolute_error'])+mp.mpf(str(rr['estimated_error']))
        er=mp.mpf(ref['error_estimate']);lit=mp.mpf(ref['C_ref'])
        omega=(ss-1)/2
        lines.extend([
            '$'+number(omega,3)+'$ & Our analytic & $'+number(own,50)+r'$ \\',
            ' & Numerical difference & $'+number(num-own)+r'\;\pm\;'+number(en,2)+r'$ \\',
            ' & Reference difference & $'+number(lit-own)+r'\;\pm\;'+number(er,2)+r'$ \\'])
        if energy!=ENERGIES[-1]:lines.append(r'\addlinespace[3pt]')
        checks.append({'s':energy,'omega_over_m':mp.nstr(omega,12),'analytic_50':printed,
          'analytic_report':'total_compact_table_values.json',
          'original_analytic_report':path,'precision_pair':current['precision_pair'],
          'numerical_difference':mp.nstr(num-own,15),'numerical_error':mp.nstr(en,15),
          'reference_difference':mp.nstr(lit-own,15),'reference_error':mp.nstr(er,15)})
    tex=r'''\begin{table}[!htbp]
\centering
\comptonprecisiontable
\begin{tabular}{rll}
\toprule
$\omega/m$ & Calculation & Coefficient or difference\\
\midrule
'''+ '\n'.join(lines)+r'''
\bottomrule
\end{tabular}
\caption{Three independent evaluations of the NLO Compton correction
$\sigma_1=(\alpha_{\mathrm T}^3/m^2)C_1$, from low to high incident
energy in the electron rest frame. Our analytic coefficient $C_1$ is
evaluated from the collected formula and printed to 50 significant digits. Each difference subtracts this
coefficient from the indicated calculation. The $\pm$ entries estimate
numerical error, not experimental uncertainty. The reference is
Equation (6) of Lee, Schwartz and Zhang~\cite{qed-processes-lsz}.}
\label{tab:comp-independent-total-50}
\end{table}
'''
    write_artifact(BASE/'tex/independent_total_table.tex' , tex)
    write_artifact(BASE/'independent_total_table_check.json' , json.dumps(checks,indent=2)+'\n')
    for row in checks:print(row,flush=True)

if __name__=='__main__':main()
