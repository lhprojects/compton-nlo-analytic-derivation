"""Reproduce the converged extreme-energy real-emission quadrature audit.

Requires Cuba (QFT_CUBA_PREFIX contains cuba.h and libcuba.a). The recorded
refinement error rule is fixed before comparison with analytic expressions.
"""
from compton_artifacts import write_artifact
import concurrent.futures,json,os,subprocess,time
from pathlib import Path
from compton_real_analytic_numerical_audit import digest,ensure_engine,SOURCE,ENGINE as DOUBLE,CUBA,CACHE
BASE=Path(__file__).resolve().parent
OUT=BASE.parent/'compton_real_analytic_numerical_wide_audit.json'
SRC=BASE/'compton_real_analytic_numerical_longdouble.cpp';EXTENDED=str(CACHE/'compton_real_longdouble')

def one(s,tol,engine,maxeval=5000000,cores=1):
 cmd=[engine,s,str(tol),str(maxeval)];start=time.time()
 p=subprocess.run(cmd,capture_output=True,text=True,check=True,env={**os.environ,'CUBACORES':str(cores)})
 r=json.loads(p.stdout.strip().splitlines()[-1]);r.update(command=cmd,requested_relative_tolerance=tol,requested_absolute_tolerance=tol*.01,elapsed_seconds=time.time()-start,cubacores=cores)
 print(json.dumps({'s':s,'level':r}),flush=True);return r

def run(s):
 levels=[]
 if s=='1.0002':tols=(1e-7,1e-8,1e-9,1e-10,1e-11)
 elif s=='1.02':tols=(1e-5,1e-6,1e-7)
 elif s=='2001':tols=(1e-4,3e-5,1e-5)
 else:tols=(1e-3,3e-4,1e-4,3e-5)
 for tol in tols:
  high=s=='200001';r=one(s,tol,DOUBLE if high else EXTENDED,15000000 if high and tol==3e-5 else 5000000,2 if high and tol==3e-5 else 1)
  if levels:r['change_from_previous']=r['real_subtracted']-levels[-1]['real_subtracted']
  levels.append(r)
 a,b=levels[-2:];out={'levels':levels,'value':b['real_subtracted'],
  'estimated_error':2*max(a['error'],b['error'],abs(a['real_subtracted']-b['real_subtracted'])),
  'all_requested_tolerances_met':all(a['fail']==0 for a in levels),
  'arithmetic':'long double amplitudes, double Cuba weights' if s!='200001' else 'double with independent long double repeat of penultimate level'}
 if s=='200001':
  out['longdouble_precision_check']=one(s,1e-4,EXTENDED)
  out['double_longdouble_difference']=out['longdouble_precision_check']['real_subtracted']-levels[-2]['real_subtracted']
 return out

if __name__=='__main__':
 commands=[ensure_engine(SOURCE,DOUBLE),ensure_engine(SRC,EXTENDED)]
 result={'source':str(SRC),'source_sha256':digest(SRC),'double_source_sha256':digest(SOURCE),
  'engine_sha256':digest(EXTENDED),'double_engine_sha256':digest(DOUBLE),'cuba_library_sha256':digest(CUBA),
  'compile_commands':commands,'energies':{},'error_rule':'2 max(last two Cuhre estimates, last refinement change)'}
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
  futures={pool.submit(run,s):s for s in ['1.0002','1.02','2001','200001']}
  for f in concurrent.futures.as_completed(futures):
   result['energies'][futures[f]]=f.result();write_artifact(OUT , json.dumps(result,indent=2)+'\n')
