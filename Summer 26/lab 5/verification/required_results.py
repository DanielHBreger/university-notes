"""Produce the M2 resistance comparison, M3 gallery, and M4/M5 examples."""
from pathlib import Path
import sys,json,csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'chua'),str(ROOT)]
from lorenz_map import draw, frame
from lyapunov import lyapunov, diagnostic
from rosenstein import load_channels, rosenstein
from simulate import make_g, fitted_offset, equilibria, stability, integrate, R0, C1, C2, L
OUT=ROOT/'verification'; AUDIT=json.loads((OUT/'data_audit.json').read_text())


def peaks(r):
    with np.load(OUT/'cache'/f"{r['peak_cache']}.npz") as data:
        return data['M'].copy(),data['t_peaks'].copy()


def required_examples():
    candidates=[]
    for r in AUDIT:
        if str(Path(r['path']).parent).replace('\\','/')!='chua/sweep forward': continue
        if r['resistance_status']!='ok' or r.get('clip_pct',100)>2 or r.get('n_maxima',0)<100: continue
        M,tp=peaks(r)
        res=lyapunov(M,tp,quantum=r['quantum_V'])
        if res['status']=='ok' and res['n_branches']>=2:
            candidates.append((r,res))
    chosen=[]
    for target in [450,550,650]:
        pool=[p for p in candidates if p[0]['path'] not in [q[0]['path'] for q in chosen]]
        if pool: chosen.append(min(pool,key=lambda p:abs(p[0]['rpot_ohm']-target)))
    if len(chosen)<3: raise RuntimeError('fewer than three resolved double-scroll maps')
    fig,axes=plt.subplots(1,3,figsize=(14,5),layout='constrained')
    allM=np.concatenate([peaks(r)[0] for r,_ in chosen]); pad=.05*np.ptp(allM)
    lim=(allM.min()-pad,allM.max()+pad)
    summaries=[]
    for ax,(r,res) in zip(axes,chosen):
        M,tp=peaks(r)
        t,v1,v2=load_channels(ROOT/r['path'])
        ros=rosenstein(v1,v2,r['dt_s'],r['period_samples'])
        alternatives=[]
        for lo,hi in [(0.5,2),(0.5,2.5),(1,3)]:
            time_in_periods=ros['t']/(r['dt_s']*r['period_samples'])
            take=(time_in_periods>=lo)&(time_in_periods<=hi)&np.isfinite(ros['curve'])
            alternatives.append(float(np.polyfit(ros['t'][take],ros['curve'][take],1)[0]))
        title=f"{Path(r['path']).name} · Rpot = {r['rpot_ohm']:.1f} Ω\n{r['duration_s']*1e3:.0f} ms · {len(M)} maxima · resolution {r['quantum_V']*1e3:.0f} mV"
        draw(ax,M,title,lim,size=5,alpha=.5); ax.set_xlabel('$M_n$ (V)');ax.set_ylabel('$M_{n+1}$ (V)')
        diagnostic(OUT/(Path(r['path']).stem+'_M5.png'),M,tp,res,
                   f"M5 · {Path(r['path']).name}, Rpot = {r['rpot_ohm']:.1f} Ω",.05,ros,r['dt_s']*r['period_samples'])
        summary=dict(path=r['path'],rpot_ohm=r['rpot_ohm'],duration_s=r['duration_s'],n_maxima=len(M),
                     residual_pct=r['residual_pct'],quantum_V=r['quantum_V'],
                     mean_T_us=res['mean_T']*1e6,median_T_us=res['median_T']*1e6,
                     mean_median_difference_pct=res['T_diff_pct'],lambda_map_per_s=res['lam'],
                     lambda_map_err_per_s=res['lam_err'],map_r2=res['r2'],map_coverage=1-res['cluster_pct']/100,
                     lambda_rosenstein_per_s=ros['lam'],lambda_rosenstein_block_sem=ros['lam_err'],
                     rosenstein_fit_r2=ros['fit_r2'],rosenstein_fit_sensitivity_per_s=alternatives)
        summaries.append(summary)
    fig.suptitle('M4 · Three long double-scroll records · same axes for comparison')
    fig.savefig(OUT/'M4_return_maps.png',dpi=220);plt.close(fig)
    (OUT/'M4_M5_results.json').write_text(json.dumps(summaries,indent=2))
    with (OUT/'M4_M5_results.csv').open('w',newline='') as fh:
        w=csv.DictWriter(fh,list(summaries[0]));w.writeheader();w.writerows(summaries)
    print('M4/M5 examples:',[(x['path'],x['lambda_map_per_s'],x['lambda_rosenstein_per_s']) for x in summaries],flush=True)


def gallery():
    from scipy.signal import savgol_filter
    names=['limitcycle1 - 773.0 ohm.csv','limitcycle2 - 764.6 ohm.csv','chaos - 716.9 ohm.csv',
           'doublescroll start - 666.0 ohm.csv','doublescroll 2 - 656.7ohm.csv','doublescroll end - 316.1ohm.csv']
    labels=['Recorded limit cycle','Recorded period 2','Recorded single scroll',
            'Recorded double-scroll start','Recorded double scroll','Recorded double-scroll end']
    fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
    named=[]
    for ax,name,label in zip(axes.flat,names,labels):
        r=next(r for r in AUDIT if r['path']=='chua/'+name)
        t,v1,v2=load_channels(ROOT/r['path'])
        stride=max(1,len(t)//100000)
        window=max(5,(r['period_samples']//20)|1)
        ax.scatter(v1[::stride],v2[::stride],s=3,alpha=.06,lw=0,color='0.5',rasterized=True)
        smooth1,smooth2=savgol_filter(v1,window,3),savgol_filter(v2,window,3)
        ax.scatter(smooth1[::stride],smooth2[::stride],s=.8,alpha=.12,lw=0,color='#0072B2',rasterized=True)
        ax.set(xlabel='$V_1$ (V)',ylabel='$V_2$ (V)',title=f"{label}\nRpot = {r['rpot_ohm']:.1f} Ω; {r['duration_s']*1000:.0f} ms; SG {window} samples")
        ax.grid(alpha=.15); ax.ticklabel_format(useOffset=False)
        named.append({k:r[k] for k in ['path','rpot_ohm','residual_pct','duration_s','quantum_V','clip_pct']})
    fig.suptitle('M3 / N2 · Recorded phase portraits · grey: raw codes; blue: light smoothing\n'
                 'Acquisition labels are retained; the final record appears nearly periodic')
    fig.savefig(OUT/'M3_phase_portraits.png',dpi=220);plt.close(fig)
    (OUT/'named_record_resistances.json').write_text(json.dumps(named,indent=2))


def resistance_check():
    fit=json.loads((ROOT/'diode_fit.json').read_text())
    scenarios=[]
    for label,i0 in [('Measured offset retained',0),('Offset removed (assumption)',fitted_offset())]:
        g=make_g(i0=i0)
        for rpot in [315.6318450984,642.874204691,654.765717346,760.123984374]:
            roots=equilibria(g,R0+rpot,0)
            for v in roots:
                eig,G=stability(g,v,R0+rpot,0)
                real=eig[np.abs(eig.imag)<1e-5]
                pair=eig[eig.imag>1e-5]
                sigma=float(pair[0].real) if len(pair) else None
                gamma=float(real[0].real) if len(real)==1 else None
                scenarios.append(dict(calibration=label,rpot_ohm=rpot,total_ohm=R0+rpot,
                                      v1_equilibrium_V=float(v),slope_S=float(G),
                                      eigenvalues=[str(x) for x in eig],
                                      sigma_per_s=sigma,gamma_per_s=gamma,
                                      shilnikov_magnitude_ratio=abs(sigma/gamma) if sigma is not None and gamma else None))
    hopf=[]
    for branch in [fit['segments'][1],fit['segments'][3]]:
        G=branch['slope_S']
        total=-G*(1+C1/C2)/(G*G+C1*C1/(C2*L))
        # Characteristic eigenfrequency at the slope-only Hopf condition.
        g=lambda v:G*np.asarray(v)
        eig,_=stability(g,0,total,0)
        hopf.append(dict(branch=branch['label'],total_ohm=total,pot_ohm=total-R0,
                         frequency_hz=float(np.max(np.abs(eig.imag))/(2*np.pi))))
    low=fit['ranges'][0]['total_lo_ohm']
    result=dict(coupling_resistors_ohm=[990,'Rpot'],shunt_ohm=216,shunt_in_chaos_path=False,
                three_equilibria_slope_total_lower_ohm=low,three_equilibria_slope_pot_lower_ohm=low-R0,
                recorded_end_pot_divider_ohm=315.6318450984,recorded_end_total_ohm=R0+315.6318450984,
                total_shortfall_ohm=low-(R0+315.6318450984),
                offset_retained_A=fitted_offset(),equilibria=scenarios,hopf_slope_only=hopf)
    (OUT/'resistance_check.json').write_text(json.dumps(result,indent=2))
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
    for ax,shift,title in zip(axes,[0,R0],['Potentiometer resistance','Total coupling resistance (990 Ω + Rpot)']):
        for i,r in enumerate(fit['ranges']):
            lo=r['pot_lo_ohm']+shift;hi=r['pot_hi_ohm']+shift
            ax.plot([lo,hi],[i,i],lw=8,alpha=.5,color=f'C{i}',label=f"Slope prediction: {r['side']} shoulder")
        ax.plot(np.array([315.6318450984,654.765717346])+shift,[-.6,-.6],'o-',color='black',lw=1.2,
                label='Records labelled double-scroll end/start')
        ax.axvspan(shift,shift+1000,color='0.6',alpha=.08,label='Physical potentiometer range')
        ax.set(yticks=[],xlabel='Resistance (Ω)',title=title,ylim=(-1,1.6));ax.grid(axis='x',alpha=.2)
    axes[0].legend(fontsize=8,loc='upper left')
    fig.suptitle('M2 / M3 · Same comparison on both resistance scales · shunt excluded in CHAOS mode')
    fig.savefig(OUT/'M2_M3_resistance_comparison.png',dpi=220);plt.close(fig)
    print('Resistance decomposition checked:',result['total_shortfall_ohm'],'ohm remaining slope-bound shortfall',flush=True)


def convergence():
    rows=[]
    for i0 in [0,fitted_offset()]:
        g=make_g(i0=i0)
        ts,ys=integrate([316,660,773],1,.003,.5e-6,0,g)
        tf,yf=integrate([316,660,773],1,.003,.25e-6,0,g)
        # Compare before exponential separation can invalidate pointwise tests.
        err=np.max(np.abs(ys-yf[::2]),axis=0)
        for k,r in enumerate([316,660,773]):
            rows.append(dict(offset_removed_A=i0,rpot_ohm=r,
                             max_v1_difference_V=float(err[0,k]),max_v2_difference_V=float(err[1,k]),
                             max_iL_difference_A=float(err[2,k])))
    (OUT/'convergence.json').write_text(json.dumps(rows,indent=2))
    print('Time-step comparison:',rows,flush=True)


if __name__=='__main__':
    resistance_check(); gallery();required_examples();convergence()
