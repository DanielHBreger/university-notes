"""Revalidate every retained raw recording and its published derived peaks."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scope_data import read_scope, R0
from rpot import fit_divider
from lorenz_map import maxima_from_arrays

HERE=Path(__file__).resolve().parent

def main():
    rows=[]
    for sw in ['forward','back']:
        tab=pd.read_csv(HERE/f'{sw}_rpot.csv').set_index('filename')
        peaks=pd.read_csv(HERE/f'{sw}_bifurcation_points.csv')
        byfile={fn:g.max_v.to_numpy() for fn,g in peaks.groupby('filename',sort=False)}
        for path in sorted((HERE/sw).glob('*.csv'),key=lambda p:int(p.stem[5:])):
            t,v=read_scope(path,('CH1','CH2','CH3'))
            row=dict(sweep=sw,filename=path.name,n_samples=len(t),duration_s=t[-1]-t[0])
            try:
                r,res=fit_divider(v)
                row.update(rpot_ohm=r,residual_pct=res*100,status='ok' if 0<=r<=1000 and res<=.05 else 'rejected')
                row['divider_difference_ohm']=r-float(tab.loc[path.name,'rpot_ohm'])
                if row['status']=='ok':
                    m,info=maxima_from_arrays(t,v[:,0])
                    old=byfile.get(path.name,np.empty(0))
                    row.update(maxima=len(m),peak_count_matches=len(m)==len(old),
                               max_peak_difference_V=float(max(abs(m-old))) if len(m)==len(old) and len(m) else None,
                               clip_pct=info['clip_pct'])
            except ValueError as e:row.update(status='rejected',reason=str(e))
            rows.append(row)
            if len(rows)%30==0:print(f'checked {len(rows)}',flush=True)
    df=pd.DataFrame(rows);df.to_csv(HERE/'measurement_audit.csv',index=False)
    summary=dict(raw_records=len(rows),admitted=int(df.status.eq('ok').sum()),
        counts=df.groupby(['sweep','status']).size().to_dict(),
        max_divider_difference_ohm=float(df.divider_difference_ohm.abs().max()),
        peak_count_mismatches=int((df.peak_count_matches==False).sum()),
        max_peak_difference_V=float(df.max_peak_difference_V.max()))
    summary['counts']={f'{a}/{b}':int(n) for (a,b),n in summary['counts'].items()}
    (HERE/'measurement_audit.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(summary)

if __name__=='__main__':main()
