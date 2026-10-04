"""Check that generators reproduce TeX and README formulas without writing files."""
from compton_artifacts import read_artifact
from pathlib import Path
import runpy,sys,difflib
repo=Path(__file__).resolve().parents[1]
failures=[]
sys.path.insert(0,str(repo/'scripts'))
names=['compton_total_energy_tex.py','compton_total_energy_compact_tex.py','compton_total_simplification_tex.py','compton_real_equivalence_tex.py','compton_real_four_matched_tex.py','compton_soft_uv_uw_compact.py','compton_real_analytic_tex.py','compton_real_analytic_gpl_tex.py','compton_real_analytic_logs_tex.py','generate_compton_closed_table.py']
def capture(self,text,*args,**kwargs):
 if self.suffix=='.tex' or self==repo/'README.md':
  old=read_artifact(self) if self.exists() else ''
  diff=''.join(difflib.unified_diff(old.splitlines(True),text.splitlines(True),fromfile=str(self),tofile='generated'))
  if diff:
   failures.append(str(self))
   print(diff)
 return len(text)
Path.write_text=capture
for name in names:
 print('GENERATOR',name,flush=True)
 try:runpy.run_path(str(repo/'scripts'/name),run_name='__main__')
 except Exception as e:
  failures.append(name)
  print('ERROR',type(e).__name__,str(e))

if failures:
 raise SystemExit(f"Generator mismatches or failures: {failures}")
print(f"All {len(names)} generators reproduce their TeX sources and README formulas.")
