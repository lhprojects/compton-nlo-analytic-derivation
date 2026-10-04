"""Independent analytic Compton NLO total assembled from amplitude integrals.

Production uses rational/log/Li2/Li3/GPL formulas and convergent series only.
No literature total and no numerical phase-space quadrature is imported.
The physical cross section is sigma_Born + alpha_T**3/m**2 * coefficient.
"""
from compton_artifacts import write_artifact
import argparse,json
from pathlib import Path
import mpmath as mp
from compton_real_analytic_sector_eval import real_sector_zcut
from compton_box_analytic_bose_eval import bose_zcut
from compton_real_endpoint_matching import endpoint_matching_closed,combine_real_finite
from compton_soft_uv_check import real_uv_closed
from compton_soft_uv_uw_finite import real_uw_closed
from compton_soft_uv_interface import soft_angular_closed
from compton_soft_analytic import ir_angular_closed
from compton_virtual_angular_analytic import coefficient_closed
from compton_box_analytic_closed import box_angular_closed


def independent_total(s,lam='1e-5',dps=40,progress=None):
    """Return components and C1(s-1), with s=s_phys/m**2 and lam=lambda_phys/m."""
    with mp.workdps(dps+40):
        s=mp.mpf(s);lam=mp.mpf(lam);h=s-1
        if h<=0 or lam<=0:raise ValueError('s>1 and lambda>0 are required')
        parts={}
        def record(k,v):
            v=mp.mpc(v)
            if abs(v.imag)>mp.mpf(10)**(-dps)*max(1,abs(v.real)):
                raise ArithmeticError(f'{k}: uncancelled imaginary part {v.imag}')
            parts[k]=v.real
            if progress:progress(k,v.real,parts)
            return v.real
        vt=record('virtual_vertex_self',coefficient_closed(s,lam))
        box=record('virtual_box',box_angular_closed(s,lam,dps=dps))
        ck=record('soft_K',soft_angular_closed(s,dps=dps))
        A=record('soft_pole',ir_angular_closed(s)/s)
        uv=real_uv_closed(h,dps=dps);uw=real_uw_closed(h,dps=dps)
        for name,data in [('UV',uv),('UW',uw)]:
            for key,value in data.items():record(name+'_'+key,value)
        finite={}
        for name,code in [('rational','rational'),('logs','logs'),('LW','W')]:
            value=real_sector_zcut(h,code,dps=dps)
            finite[name]=record(name+'_zfinite',value['finite_zcut'])
            record(name+'_pole',value['pole'])
        value=bose_zcut(h,dps=dps)
        finite['bose']=record('bose_zfinite',value['finite_zcut'])
        record('bose_pole',value['pole'])
        matching=endpoint_matching_closed(h,dps=dps)
        for name,data in matching.items():record(name+'_matching',data['finite_shift'])
        realpole=(uv['pole']+2*uw['pole']+sum(parts[n+'_pole'] for n in ['rational','logs','LW','bose']))/(8*s*h)
        residual=record('ir_pole_residual',realpole-A)
        if abs(residual)>mp.mpf(10)**(-dps+5)*max(1,abs(A)):
            raise ArithmeticError('Real and virtual infrared coefficients disagree')
        fr=record('real_finite',combine_real_finite(finite,uv['finite'],uw['finite'],h,dps=dps))
        scale=record('regulator_conversion',A*mp.log(1/(mp.sqrt(s)*h*lam)))
        record('total',vt+box+fr+scale-ck/s)
        return parts


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--s',default='1.6');p.add_argument('--lambda',dest='lam',default='1e-5')
    p.add_argument('--dps',type=int,default=35)
    p.add_argument('--report',help='Optional JSON output; otherwise print results only')
    args=p.parse_args();path=Path(args.report) if args.report else None
    report={'s':args.s,'lambda':args.lam,'dps':args.dps,'complete':False,'components':{}}
    def progress(k,v,parts):
        report['components']={n:mp.nstr(x,args.dps+10) for n,x in parts.items()}
        print(k,mp.nstr(v,args.dps),flush=True)
    independent_total(args.s,args.lam,args.dps,progress)
    report['complete']=True
    if path is not None:
        write_artifact(path , json.dumps(report,indent=2)+'\n')
