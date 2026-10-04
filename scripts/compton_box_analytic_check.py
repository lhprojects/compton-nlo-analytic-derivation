"""Cross-check a closed box against independently derived parameter integrals.

Quadrature here is a verification only. The evaluated answer comes from
compton_box_analytic.AnalyticBox, which performs no quadrature.
"""
import argparse,json
import numpy as np
from compton_box_analytic import AnalyticBox
from compton_physical_virtual import PhysicalLoop,cm_event
from check_compton_amplitudes import dot,slash,pol

def run(order):
    report=[]
    for s,c in [(1.6,.2),(10.,-.4),(100.,.7)]:
        p,k,kp,pp=cm_event(s,c)
        for Q in (p+k,p-kp):
            closed=AnalyticBox(p,pp,Q,1e-5)
            residual=[]
            for e in pol(k):
                for ep in pol(kp):
                    A,B=slash(e),slash(ep)
                    result=closed.matrix(A,B)
                    check=PhysicalLoop(order,1e-5).spanning(p,pp,Q,A,B)
                    residual.append(np.linalg.norm(result-check)/max(1,np.linalg.norm(result)))
            report.append(dict(s=s,c=c,Q_squared=float(dot(Q,Q)),quadrature_order=order,
                               max_relative_matrix_difference=float(max(residual)),
                               rank_three_uv_residue=float(np.linalg.norm(closed.tensor((0,1,2,3),3)[0]))))
    return report
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--order',type=int,default=96)
    args=parser.parse_args()
    rows=run(args.order)
    print(json.dumps(rows,indent=2))
    for row in rows:
        assert row['max_relative_matrix_difference'] < 1e-9, row
        assert row['rank_three_uv_residue'] < 1e-11, row
