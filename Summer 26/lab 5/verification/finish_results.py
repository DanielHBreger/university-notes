"""Refresh map uncertainties and resistance labels using audited peak caches.

Previously computed direct estimates are retained; no direct curve is
invented or reused for a different recording.
"""
from pathlib import Path
import sys,json,csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'chua'))
from lyapunov import lyapunov,CSV_FIELDS,num,rpot_from_name
from simulate import load_measured_bifurcation,plot_bifurcation
AUDIT=json.loads((ROOT/'verification/data_audit.json').read_text())
INDEX={r['path']:r for r in AUDIT}


def refresh(folder,stem):
    path=folder/(stem+'.csv'); old=pd.read_csv(path,keep_default_na=False)
    rows=[]
    for _,record in old.iterrows():
        raw=(folder/('' if stem=='records_lyapunov' else stem.removesuffix('_lyapunov'))/record['filename']).relative_to(ROOT).as_posix()
        r=INDEX[raw]
        data=np.load(ROOT/'verification/cache'/f"{r['peak_cache']}.npz")
        res=lyapunov(data['M'],data['t_peaks'],quantum=r['quantum_V'])
        for csv_name,key in [('lambda_rosenstein_per_s','lam_ros'),('lambda_rosenstein_err_per_s','lam_ros_err'),('rosenstein_fit_r2','ros_fit_r2')]:
            if record.get(csv_name,'')!='': res[key]=float(record[csv_name])
        res['ros_status']=record.get('rosenstein_status','not requested')
        rows.append(dict(res,sweep=record['sweep'],filename=record['filename'],rpot=r['rpot_ohm'],
                         rpot_source='divider fit (R0=990 ohm)',
                         nominal_rpot=rpot_from_name(record['filename']) if stem=='records_lyapunov' else np.nan,
                         mean_T_us=res['mean_T']*1e6,median_T_us=res['median_T']*1e6))
    with path.open('w',newline='') as fh:
        w=csv.writer(fh);w.writerow([c for c,_,_ in CSV_FIELDS])
        w.writerows([[num(r[key],nd) for _,key,nd in CSV_FIELDS] for r in rows])
    fig,ax=plt.subplots(figsize=(10,5.5),layout='constrained')
    for key,err,label,marker,color in [('lam','lam_err','Return map','o','#0072B2'),
                                       ('lam_ros','lam_ros_err','Direct time series','s','#D55E00')]:
        sel=[r for r in rows if np.isfinite(r[key])]
        if sel:
            ax.errorbar([r['rpot'] for r in sel],[r[key] for r in sel],yerr=[r[err] for r in sel],
                        fmt=marker,ms=4,capsize=2,elinewidth=.6,color=color,
                        mfc='none' if marker=='s' else color,label=f'{label} ({len(sel)} records)')
    ax.axhline(0,color='0.5',ls='--',lw=.8)
    ax.set(xlabel='$R_{pot}$ (Ω)',ylabel='Lyapunov estimate (s$^{-1}$)',
           title=stem.replace('_lyapunov','')+' · Lyapunov estimates\nError bars show precision; model and calibration bias are additional')
    if ax.get_legend_handles_labels()[0]: ax.legend(fontsize=9)
    else: ax.text(.5,.5,'No resolved map estimates; see the CSV for reasons',ha='center',transform=ax.transAxes)
    ax.grid(alpha=.15);fig.savefig(folder/(stem+'.png'),dpi=220);plt.close(fig)
    return rows


def redraw_simulated_overlay():
    folder=ROOT/'chua';measured=load_measured_bifurcation(folder)
    for tag in ['', '_corrected']:
        data=pd.read_csv(folder/f'simulated_bifurcation{tag}_points.csv')
        rs=np.sort(data.rpot_ohm.unique())
        maps=[]
        for start in [1,-1]:
            groups=data[data.initial_v1_V==start].groupby('rpot_ohm')
            maps.append([groups.get_group(r).max_v.to_numpy() if r in groups.groups else np.empty(0) for r in rs])
        note='Offset removed (assumption)' if tag else 'Measured current offset retained'
        plot_bifurcation(rs,*maps,measured,folder/f'simulated_bifurcation{tag}.png',note)


if __name__=='__main__':
    for folder in [ROOT/'chua',ROOT/'chua/set 2']:
        for sweep in ['sweep forward','sweep back']:
            refresh(folder,sweep+'_lyapunov')
    refresh(ROOT/'chua','records_lyapunov')
    redraw_simulated_overlay()
    print('Final numerical summaries and overlays refreshed.')
