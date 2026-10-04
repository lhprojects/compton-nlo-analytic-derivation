"""Generate the experimental comparison from our independent analytic total.

Only GPL evaluation is used for the NLO theory. No literature total is imported.
Every run evaluates every theory point afresh. Saved reports are outputs only.
"""
from compton_artifacts import read_artifact, write_artifact
from pathlib import Path
import argparse, json, math
from concurrent.futures import ProcessPoolExecutor, as_completed
import mpmath as mp
from compton_total_independent_closed import independent_total
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'tex/figures'
REPORT=ROOT/'comparison_checks.json'
SAMPLE_REPORT=ROOT/'experimental_curve_analytic.json'
RAW=ROOT/'data/primex_carbon_raw.dat'
ALPHA=mp.mpf(1)/mp.mpf('137.035999084')
MASS=mp.mpf('0.00051099895')
GEV2_TO_MB=mp.mpf('0.3893793721')

def born(w):
    r=mp.mpf(w)/MASS
    return 2*mp.pi*ALPHA**2/MASS**2*((1+r)/r**2*(2*(1+r)/(1+2*r)-mp.log1p(2*r)/r)+mp.log1p(2*r)/(2*r)-(1+3*r)/(1+2*r)**2)

def write(name,header,rows,comment):
    write_artifact(OUT/name , '# '+comment+'\n'+header+'\n'+''.join(' '.join(f'{float(v):.7g}' for v in row)+'\n' for row in rows))

def evaluate_sample(key,digits):
    mp.mp.dps=digits+40
    w=mp.mpf(key);ss=1+2*w/mp.mpf("0.00051099895")
    result=independent_total(ss,dps=digits)
    return key,{"dps":digits,"s":mp.nstr(ss,digits+15),"components":{k:mp.nstr(v,digits+10) for k,v in result.items()}}

def main():
    p=argparse.ArgumentParser();p.add_argument('--digits',type=int,default=15);p.add_argument('--samples',type=int,default=9);p.add_argument('--workers',type=int,default=4);a=p.parse_args()
    mp.mp.dps=a.digits+40
    global ALPHA,MASS,GEV2_TO_MB
    ALPHA=1/mp.mpf('137.035999084');MASS=mp.mpf('0.00051099895');GEV2_TO_MB=mp.mpf('0.3893793721')
    samples={}
    energies=[mp.mpf('4.35')+mp.mpf('1.2')*i/(a.samples-1) for i in range(a.samples)]
    if mp.mpf(5) not in energies:energies.append(mp.mpf(5))
    keys=[mp.nstr(w,16) for w in energies]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futures={pool.submit(evaluate_sample,key,a.digits):key for key in keys}
        for future in as_completed(futures):
            key,value=future.result();samples[key]=value
            print('Completed independent sample',key,'GeV',flush=True)
    write_artifact(SAMPLE_REPORT , json.dumps(samples,indent=2)+'\n')
    theory=[]
    for w in sorted(energies):
        key=mp.nstr(w,16)
        correction=ALPHA**3/MASS**2*mp.mpf(samples[key]['components']['total'])*GEV2_TO_MB
        b=born(w)*GEV2_TO_MB
        theory.append((w,b,correction,b+correction,100*correction/b))
    write('compton_nlo_total.dat','energy born correction nlo percent',theory,'Independent analytic amplitude integration; GeV and mb/electron; alpha^-1=137.035999084; m=0.00051099895 GeV. Virtual plus double Compton, no charged-pair final state.')
    data=[]
    for line in read_artifact(RAW).splitlines():
        if line and line[0] not in '#p':data.append(list(map(float,line.split())))
    for period in (1,2):
        rows=[]
        for per,w,sig,stat,syst in data:
            if per==period:
                b=float(born(w)*GEV2_TO_MB);err=math.hypot(stat,syst)
                rows.append((w,sig,stat,syst,err,sig/b,err/b))
        write(f'compton_nlo_primex_{period}.dat','energy sigma stat syst error ratio_born error_born',rows,'PrimEx Table 2, PLB 797 (2019) 134884; unchanged Chapter 9 raw data; errors=hypot(stat,syst).')
    five=next(row for row in theory if row[0]==5)
    report={'source':'independent_total; no literature formula','points':len(data),'analytic_samples':len(theory),'minimum_working_digits':a.digits,'curve_significant_digits':7,'five_gev':dict(zip(('energy','born','correction','nlo','percent'),map(float,five))),'percent_range':[float(min(r[4] for r in theory)),float(max(r[4] for r in theory))]}
    write_artifact(REPORT , json.dumps(report,indent=2)+'\n')
    write_artifact(OUT/'compton_nlo_values.tex' , '% Generated from independent analytic total\n'+f'\\newcommand{{\\comptonBornFive}}{{{float(five[1]):.6f}}}\n'+f'\\newcommand{{\\comptonNloFive}}{{{float(five[3]):.6f}}}\n'+f'\\newcommand{{\\comptonCorrectionFive}}{{{float(five[4]):.3f}}}\n')
    print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':main()
