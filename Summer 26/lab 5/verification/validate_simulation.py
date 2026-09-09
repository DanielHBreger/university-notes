"""Independent short-time ODE check and long-time statistical step check."""
from pathlib import Path
import sys,json
import numpy as np
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'chua'))
from simulate import make_g,fitted_offset,integrate,rhs,R0

results=[]
for offset in [0,fitted_offset()]:
    g=make_g(i0=offset)
    for r in [316,660,773]:
        row=dict(rpot_ohm=r,removed_offset_A=offset)
        t,Y=integrate([r],1,.0005,.5e-6,0,g)
        reference=solve_ivp(lambda t,y:rhs(y[:,None],np.array([R0+r]),0,g)[:,0],
                            [0,t[-1]],[1,0,0],method='DOP853',rtol=1e-10,atol=1e-12,t_eval=t,
                            max_step=.25e-6)
        row['max_v1_error_vs_DOP853_V']=float(np.max(abs(Y[:,0,0]-reference.y[0])))
        row['reference_success']=reference.success
        results.append(row)
    t,Y=integrate([316,660,773],1,.05,.25e-6,0,g,t_skip=.02)
    _,coarse=integrate([316,660,773],1,.05,.5e-6,0,g,t_skip=.02)
    for k,row in enumerate(results[-3:]):
        fine=Y[:,0,k];base=coarse[:,0,k]
        row['fine_v1_mean_std_V']=[float(fine.mean()),float(fine.std())]
        row['coarse_v1_mean_std_V']=[float(base.mean()),float(base.std())]
        row['fine_v1_quantiles_V']=np.quantile(fine,[.01,.5,.99]).tolist()
        row['coarse_v1_quantiles_V']=np.quantile(base,[.01,.5,.99]).tolist()
    print('Completed calibration',offset,flush=True)
(ROOT/'verification/simulation_validation.json').write_text(json.dumps(results,indent=2))
print(json.dumps(results,indent=2),flush=True)
