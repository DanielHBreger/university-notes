"""Reproducible partial uncertainty audit. python uncertainty/calculate.py"""
import sys,json,csv
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'uncertainty'
sys.path.insert(0,str(ROOT/'chua'))
from scope_data import read_scope
from uncertainty import (voltage_quantum,scope_vertical_limit,timebase_limit_fraction,
    supply_voltage_limit,current_covariance,divider_block_bootstrap,delta_covariance)
from cascade_periods import load_sidecar,lag_distances,period_of,bifurcation_points
C=json.loads((OUT/'instrument_inputs.json').read_text())


def channel_stats(path,t,d,channels):
    key=str(path.relative_to(ROOT)); out={}
    for j,ch in enumerate(channels):
        q,quality=voltage_quantum(d[:,j]); settings=dict(C['scope_channels_default'][ch])
        settings.update(C['record_channel_overrides'].get(key,{}).get(ch,{}))
        limit=scope_vertical_limit(settings['volts_per_div_at_tip'],settings['offset_V_at_tip'],settings['probe_attenuation'])
        out[ch]=dict(quantum_V=q,quantum_grid_consistency=quality,u_quantization_V=q/np.sqrt(12),
                     dc_reading_limit_V=limit,dc_limit_status='unknown settings' if limit is None else 'conditional typical DC envelope; probe accuracy additional',**settings)
    return out


def diode_audit():
    t,d=read_scope(ROOT/'trace1.csv'); dt=float(np.median(np.diff(t)))
    stats=channel_stats(ROOT/'trace1.csv',t,d,['CH1','CH2'])
    vv,vi=d.T; current=(vv-vi)/C['shunt_ohm']; cfg=json.loads((ROOT/'diode_fit.json').read_text())
    nblocks=32; blocks=np.array_split(np.arange(len(t)),nblocks)
    rng=np.random.default_rng(20260915); weights=rng.multinomial(nblocks,np.full(nblocks,1/nblocks),size=1000)
    keep=(vv>=cfg['source_fit_range_V'][0])&(vv<=cfg['source_fit_range_V'][1]); segments=[]
    for k,s in enumerate(cfg['segments']):
        mask=keep&(vi>=s['v_lo'])&((vi<s['v_hi']) if k<4 else (vi<=s['v_hi']))
        sums=[]
        for block in blocks:
            idx=block[mask[block]];x=vi[idx];y=current[idx]
            sums.append([len(x),x.sum(),y.sum(),(x*x).sum(),(x*y).sum()])
        samples=weights@np.array(sums);n,sx,sy,sxx,sxy=samples.T
        den=n*sxx-sx*sx;valid=(n>10)&(den>1e-15)
        a=(n[valid]*sxy[valid]-sx[valid]*sy[valid])/den[valid];b=(sy[valid]-a*sx[valid])/n[valid]
        segments.append(dict(label=s['label'],slope_S=s['slope_S'],intercept_A=s['intercept_A'],
                             u_slope_block_S=float(a.std(ddof=1)),u_intercept_block_A=float(b.std(ddof=1)),
                             cov_slope_intercept=np.cov(a,b).tolist(),valid_draws=int(valid.sum()),
                             slope_ci95_block_S=np.quantile(a,[.025,.975]).tolist(),
                             v_lo=s['v_lo'],v_hi=s['v_hi'],total_uncertainty=None))
    cov=current_covariance(stats['CH1']['u_quantization_V'],stats['CH2']['u_quantization_V'],C['shunt_ohm'])
    out=dict(channels=stats,sample_dt_s=dt,cov_VI_I_quantization=cov.tolist(),u_current_quantization_A=float(np.sqrt(cov[1,1])),
             segments=segments,block_s=(t[-1]-t[0])/nblocks,
             note='Paired time-block bootstrap, fixed nominal segment boundaries; calibration, EIV bias, breakpoint selection and shunt uncertainty excluded.')
    (OUT/'diode_uncertainty.json').write_text(json.dumps(out,indent=2)+'\n')


def sweep_audit(sweep):
    table=pd.read_csv(ROOT/'chua'/f'{sweep}_rpot.csv'); table=table[table.status=='ok']
    rows=[]
    for j,row in enumerate(table.itertuples()):
        path=ROOT/'chua'/sweep/row.filename; t,d=read_scope(path,('CH1','CH2','CH3'))
        dt=float(np.median(np.diff(t))); stats=channel_stats(path,t,d,['CH1','CH2','CH3'])
        rec=dict(sweep=sweep,filename=row.filename,rpot_ohm=row.rpot_ohm,dt_s=dt,n_samples=len(t),
                 u_r_total_ohm=None,uncertainty_status='partial: calibration/probe/R0 terms unresolved')
        for ch,st in stats.items():
            for key in ['quantum_V','quantum_grid_consistency','u_quantization_V','dc_reading_limit_V']:
                rec[f'{ch}_{key}']=st[key]
        try:
            results=[divider_block_bootstrap(d,dt,block_s=b,draws=400,seed=20260915+j,r0=C['r0_ohm']) for b in [.005,.010]]
            rec['u_r_block_5ms_ohm']=results[0]['u_r_ohm']; rec['u_r_block_10ms_ohm']=results[1]['u_r_ohm']
            selected=max(results,key=lambda x:x['u_r_ohm'])
            rec['u_r_block_ohm']=selected['u_r_ohm'];rec['block_ci95_low_ohm'],rec['block_ci95_high_ohm']=selected['ci95_ohm']
            rec['block_duration_selected_s']=selected['block_s']
            if C['r0_standard_uncertainty_ohm'] is not None and C['relative_gain_ratio_standard_uncertainty'] is not None:
                rec['u_r_total_ohm']=float(np.sqrt(selected['u_r_ohm']**2+(row.rpot_ohm/C['r0_ohm']*C['r0_standard_uncertainty_ohm'])**2+(row.rpot_ohm*C['relative_gain_ratio_standard_uncertainty'])**2))
                rec['uncertainty_status']='fit + supplied calibration terms; other model effects not included'
        except (ValueError,np.linalg.LinAlgError) as e:
            rec['uncertainty_status']='unidentifiable bootstrap: '+str(e)
        rows.append(rec)
        if j%40==0: print(sweep,j+1,'/',len(table),flush=True)
    pd.DataFrame(rows).to_csv(OUT/f'{sweep}_uncertainty.csv',index=False)
    return pd.DataFrame(rows)


def cascade_audit():
    records=load_sidecar('forward',str(ROOT/'chua'));rows=[]
    for name,(rp,m) in records.items():
        if 744<=rp<=790:
            p,s,f,dr=period_of(lag_distances(m))
            if dr<2.5: rows.append((rp,p,s))
    points=bifurcation_points(rows); r=np.array([points[k][0] for k in [1,2,3]])
    half=np.array([points[k][1] for k in [1,2,3]])
    delta,sens=delta_covariance(r,np.diag(half**2))
    # Actual observed brackets, uniform location within each interval. These
    # give a separate interval-only model, not a re-fit or formal confidence CI.
    rng=np.random.default_rng(20260915);brackets=np.array([points[k][2] for k in [1,2,3]])
    sample=rng.uniform(brackets[:,0],brackets[:,1],size=(200000,3))
    ds=(sample[:,0]-sample[:,1])/(sample[:,1]-sample[:,2])
    out=dict(nominal_transition_R_ohm=r.tolist(),half_brackets_ohm=half.tolist(),brackets_ohm=brackets.tolist(),
             nominal_delta=delta,delta_half_bracket_sensitivity=sens,
             note='Half-brackets are resolution descriptors, NOT standard errors. Correct covariance includes R2 in both gaps.',
             bracket_uniform_delta_median=float(np.median(ds)),bracket_uniform_delta_95_interval=np.quantile(ds,[.025,.975]).tolist(),
             bracket_uniform_delta_std=float(ds.std(ddof=1)),
             interval_note='Conditional on independent uniform locations in the measured brackets. Does not include R calibration, drift, or classification uncertainty.')
    (OUT/'cascade_uncertainty.json').write_text(json.dumps(out,indent=2)+'\n')


def main():
    diode_audit()
    tables={s:sweep_audit(s) for s in ['forward','back']}
    cascade_audit()
    for sweep,table in tables.items():
        ly=pd.read_csv(ROOT/'chua'/f'{sweep}_lyapunov.csv'); ly=ly.merge(table[['filename','u_r_block_ohm','dt_s']],on='filename',how='left')
        frac=timebase_limit_fraction(C['scope_years_since_calibration'])
        ly['u_lambda_timebase_per_s']=np.nan if frac is None else ly.lambda_rosenstein_per_s.abs()*frac/np.sqrt(3)
        ly['lambda_timebase_25ppm_component_per_s']=ly.lambda_rosenstein_per_s.abs()*25e-6/np.sqrt(3)
        ly['lambda_error_status']='Original fit uncertainty only; clock aging and estimator systematic uncertainty unresolved'
        ly.to_csv(OUT/f'{sweep}_lyapunov_uncertainty.csv',index=False)
    supply=dict(setting_V=C['supply_actual_setting_V'],reported_current_setting_A=C['supply_reported_current_setting_A'],
                programming_or_readback_limit_V=supply_voltage_limit(C['supply_actual_setting_V']),
                standard_u_rectangular_V=supply_voltage_limit(C['supply_actual_setting_V'])/np.sqrt(3),
                line_regulation_limit_V=.0001*C['supply_actual_setting_V']+.003,
                load_regulation_limit_V=.0001*C['supply_actual_setting_V']+.003,ripple_limit_Vrms=.001,
                current_setting_limit_A=.003*C['supply_reported_current_setting_A']+.010,
                note='CH1/CH2 at 9 V each, 0.03 A reported setting. Conditional on CV operation, warm-up, temperature. Tracking mode and actual load current not established. Regulation terms are not independent calibration errors and are not RSS-added. Ripple is not added to every scope sample.')
    (OUT/'supply_uncertainty.json').write_text(json.dumps(supply,indent=2)+'\n')
    print(json.dumps(supply,indent=2),flush=True)

if __name__=='__main__': main()
