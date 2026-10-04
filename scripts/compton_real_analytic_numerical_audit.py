"""Reproducible quadrature audit of the original six-diagram real amplitude.

This runner records requested Cuhre tolerances, evaluations, source hashes,
and refinement changes. It never imports the analytic result or a reference.
The quadrature's reported errors are estimates, not rigorous bounds.
"""
from compton_artifacts import read_artifact_bytes, write_artifact
import concurrent.futures,hashlib,json,os,subprocess,time
from pathlib import Path
BASE=Path(__file__).resolve().parent
SOURCE=BASE/'compton_real_integral.cpp'
OUTPUT=BASE.parent/'compton_real_analytic_numerical_audit.json'
CACHE=Path(os.environ.get('XDG_CACHE_HOME',Path.home()/'.cache'))/'qftbook/compton-integration'
ENGINE=str(CACHE/'compton_real_integral')
CUBA_PREFIX=Path(os.environ.get('QFT_CUBA_PREFIX',CACHE/'Cuba-4.2.2'))
CUBA=str(CUBA_PREFIX/'libcuba.a')

def digest(path):return hashlib.sha256(read_artifact_bytes(Path(path))).hexdigest()
def ensure_engine(source,engine):
 Path(engine).parent.mkdir(parents=True,exist_ok=True)
 if not Path(CUBA).is_file():raise SystemExit('Set QFT_CUBA_PREFIX to a Cuba installation containing cuba.h and libcuba.a')
 command=['g++','-O3','-std=c++17','-I'+str(CUBA_PREFIX),str(source),CUBA,'-lm','-o',engine]
 subprocess.run(command,check=True)
 return command

def run_energy(s):
 levels=[]
 for tol in (1e-5,1e-6,1e-7):
  command=[ENGINE,str(s),str(tol),'15000000'];start=time.time()
  proc=subprocess.run(command,text=True,capture_output=True,check=True,
                      env={**os.environ,'CUBACORES':'1'})
  result=json.loads(proc.stdout.strip().splitlines()[-1]);result.update(
   requested_relative_tolerance=tol,requested_absolute_tolerance=tol*.01,
   maxeval=15000000,elapsed_seconds=time.time()-start,command=command)
  if levels:result['change_from_previous']=result['real_subtracted']-levels[-1]['real_subtracted']
  levels.append(result)
  print(json.dumps({'energy':s,'level':result}),flush=True)
 fine,prior=levels[-1],levels[-2]
 estimate=2*max(fine['error'],prior['error'],abs(fine['real_subtracted']-prior['real_subtracted']))
 return {'levels':levels,'value':fine['real_subtracted'],
  'estimated_error':estimate,'error_rule':'2 max(last two reported Cuhre errors, last refinement change)',
  'all_requested_tolerances_met':all(r['fail']==0 for r in levels)}

if __name__=='__main__':
 compile_command=ensure_engine(SOURCE,ENGINE)
 result={'source':str(SOURCE),'source_sha256':digest(SOURCE),
  'engine':ENGINE,'engine_sha256':digest(ENGINE),'cuba_library_sha256':digest(CUBA),
  'compile_command':compile_command,
  'meaning':'CM soft-subtracted real coefficient R_sub = F_eta - A log(2 sqrt(s) (s-1)); sigma_R=alpha^3 R_sub/m^2',
  'cubacores':1,'energies':{}}
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
  futures={pool.submit(run_energy,s):s for s in (1.6,10,100)}
  for f in concurrent.futures.as_completed(futures):
   s=futures[f];result['energies'][str(s)]=f.result();write_artifact(OUTPUT , json.dumps(result,indent=2)+'\n')
 print('Saved '+str(OUTPUT),flush=True)
