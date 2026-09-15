"""Rebuild presentation figures from project data. Run from any directory."""
import sys
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'presentation_figures'
sys.path.insert(0, str(ROOT / 'chua'))
from scope_data import read_scope
from gallery import nearest
BLUE, ORANGE, PURPLE, INK = '#147da3', '#d87530', '#8257a6', '#243444'
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':17,
    'axes.labelsize':18,'axes.titlesize':20,'xtick.labelsize':15,'ytick.labelsize':15,
    'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'#85919b',
    'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK,
    'legend.frameon':False,'legend.fontsize':15,'svg.fonttype':'none',
    'savefig.facecolor':'white'})
manifest = []

def canvas(title, subtitle=''):
    fig = plt.figure(figsize=(13.333,7.5))
    fig.text(.075,.93,title,fontsize=27,weight='bold')
    if subtitle: fig.text(.075,.875,subtitle,fontsize=15,color='#586673')
    return fig

def footer(fig, text):
    fig.text(.075,.035,text,fontsize=11,color='#586673')

def style(ax): ax.grid(alpha=.15); ax.set_axisbelow(True)

def save(fig, name, use, source):
    fig.savefig(OUT / (name+'.png'),dpi=180)
    fig.savefig(OUT / (name+'.svg'))
    manifest.append(dict(file=name,use=use,source=source))
    plt.close(fig)
    print('Wrote',name,flush=True)

def record(r,sweep='forward',window=.02):
    name,rp=nearest(str(ROOT/'chua'/sweep),r)
    t,d=read_scope(ROOT/'chua'/sweep/name)
    mask=t-t[0]<=window
    return t[mask]-t[0],d[mask],rp,f'{sweep}/{name}'

# Title image, a measured portrait with minimal annotation.
t,d,r,src=record(550,window=.05)
fig=plt.figure(figsize=(10,6))
ax=fig.add_axes([.05,.1,.9,.82]); ax.plot(d[:,0],d[:,1],lw=.5,color=BLUE,alpha=.65)
ax.axis('off'); fig.text(.5,.035,f'Measured double scroll   •   Rpot = {r:.1f} Ω',ha='center',fontsize=16)
save(fig,'01_title_double_scroll','Title slide image',src+'; first 50 ms')

# Static I–V data, retaining original offset so the fit and measurements match.
cfg=json.loads((ROOT/'diode_fit.json').read_text())
t,vi=read_scope(ROOT/'trace1.csv')
v=vi[:,1]; current=(vi[:,0]-v)/cfg['shunt_ohm']*1000
fig=canvas('The nonlinear element has a negative-slope region','Measured current–voltage characteristic and the five-segment fit')
ax=fig.add_axes([.11,.17,.82,.62]); ax.scatter(v[::4],current[::4],s=3,color=BLUE,alpha=.22,rasterized=True,label='Measurements')
for j,s in enumerate(cfg['segments']):
    x=np.linspace(s['v_lo'],s['v_hi'],100)
    ax.plot(x,1000*(s['slope_S']*x+s['intercept_A']),color=ORANGE,lw=3,label='Segment fits' if j==0 else None)
ax.axhline(0,color='#77838b',lw=.7); ax.set(xlabel='Voltage across nonlinear element (V)',ylabel='Current (mA)',xlim=(-9,8),ylim=(-7,8)); style(ax)
ax.annotate('Negative differential resistance',xy=(1,-1.2),xytext=(-4,5.8),fontsize=17,arrowprops=dict(arrowstyle='->',color=INK))
ax.legend(loc='lower left',markerscale=3)
footer(fig,'Current = (CH1 − CH2) / 216 Ω. Raw offset retained here; simulations remove the fitted current offset.')
save(fig,'02_nonlinear_element','Nonlinear element characterization','trace1.csv; diode_fit.json')

# Four main small-attractor regimes, with explicit shared scales.
fig=canvas('Measured regimes as resistance decreases','V1 and V2 are the voltages across C1 and C2; each panel shows 20 ms')
for j,(lab,target) in enumerate([('Period 1',806),('Period 2',760),('Single scroll',700),('Double scroll',550)]):
    t,d,r,src=record(target)
    ax=fig.add_axes([.105+j*.222,.22,.195,.52]); ax.plot(d[:,0],d[:,1],lw=.55,color=BLUE,alpha=.75)
    ax.set(xlim=(-4.8,2.8),ylim=(-1.05,1.05),xlabel='V1 (V)',title=f'{lab}\n{r:.1f} Ω')
    ax.set_xticks([-4,-2,0,2]);
    if j==0: ax.set_ylabel('V2 (V)')
    else: ax.tick_params(labelleft=False)
    style(ax)
footer(fig,'Forward sweep. All four panels use identical voltage scales. R labels refer to Rpot.')
save(fig,'03_measured_regimes','Phase-portrait overview','forward traces nearest 806, 760, 700, 550 Ω')

fwd=pd.read_csv(ROOT/'chua/forward_bifurcation_points.csv')
back=pd.read_csv(ROOT/'chua/back_bifurcation_points.csv')
# Near-DC early records have unreliable fitted R; exclude as documented in RESULTS.md.
fwd=fwd[~fwd.filename.isin([f'trace{i}.csv' for i in range(14,19)])]
fig=canvas('Period doubling precedes the chaotic regimes','Each vertical slice contains the local maxima of V1 from one recording')
a=fig.add_axes([.09,.19,.48,.61]); b=fig.add_axes([.68,.19,.28,.61])
for ax,xlim,ylim in [(a,(900,0),(-3.5,8)),(b,(805,735),(-1.1,.35))]:
    df=fwd[(fwd.rpot_ohm>=min(xlim))&(fwd.rpot_ohm<=max(xlim))]
    ax.scatter(df.rpot_ohm,df.max_v,s=.35,color=BLUE,alpha=.35,lw=0,rasterized=True)
    ax.set(xlim=xlim,ylim=ylim,xlabel='Rpot (Ω)',ylabel='Local maxima of V1 (V)'); style(ax)
a.set_title('Full downward sweep'); b.set_title('Period-doubling detail')
for r,txt,y in [(774.6,'1 → 2',.23),(755.6,'2 → 4',.08),(750.5,'4 → 8',-.07)]:
    b.axvline(r,color=ORANGE,lw=1,alpha=.6); b.text(r+1,y,txt,ha='right',fontsize=14,color=ORANGE)
footer(fig,'Resistance decreases to the right. Transition estimates: 774.6 ± 1.5, 755.6 ± 0.8, 750.5 ± 0.4 Ω.')
save(fig,'04_bifurcation','Bifurcation and period doubling','forward_bifurcation_points.csv; cascade_periods_forward.txt')

fig=canvas('The observed attractor depends on sweep direction','Different orbits persist over an overlapping range of resistance')
a=fig.add_axes([.09,.20,.46,.60])
for df,col,lab in [(fwd,BLUE,'Decreasing R'),(back,PURPLE,'Increasing R')]:
    a.scatter(df.rpot_ohm,df.max_v,s=.5,color=col,alpha=.35,lw=0,rasterized=True,label=lab)
a.set(xlim=(900,0),ylim=(-3.5,8),xlabel='Rpot (Ω)',ylabel='Local maxima of V1 (V)'); style(a)
leg=a.legend(loc='upper left',markerscale=7)
for r in [328,675]: a.axvline(r,color='#6d7781',ls='--',lw=1)
for j,(sweep,col,lab) in enumerate([('forward',BLUE,'Decreasing R'),('back',PURPLE,'Increasing R')]):
    t,d,r,src=record(550,sweep)
    ax=fig.add_axes([.66,.55-j*.36,.30,.24]); ax.plot(d[:,0],d[:,1],color=col,lw=.65)
    ax.set(xlim=(-8,8),ylim=(-7,7),xlabel='V1 (V)',ylabel='V2 (V)',title=f'{lab}: {r:.1f} Ω'); style(ax)
footer(fig,'Example portraits use nearby measured resistances, not exactly equal R. Dashed lines: 328 Ω and 675 Ω.')
save(fig,'05_hysteresis','History dependence at nearly the same R','forward/back bifurcation CSVs; forward trace270 and back trace21')

ly=pd.read_csv(ROOT/'chua/forward_lyapunov.csv')
ly=ly[(ly.rosenstein_status=='ok')&~ly.filename.isin([f'trace{i}.csv' for i in range(14,19)])]
fig=canvas('Positive Lyapunov estimates in the chaotic regime','Direct Rosenstein estimates from the measured time series')
ax=fig.add_axes([.11,.19,.83,.61])
ax.axhline(0,color=INK,lw=1.5)
ax.errorbar(ly.rpot_ohm,ly.lambda_rosenstein_per_s/1000,yerr=ly.lambda_rosenstein_err_per_s/1000,fmt='o',ms=3.5,color=BLUE,alpha=.8,elinewidth=.5)
ax.set(xlim=(900,0),ylim=(-.3,2.8),xlabel='Rpot (Ω)',ylabel='Largest Lyapunov estimate (10³ s⁻¹)'); style(ax)
p=ly[ly.filename=='trace270.csv'].iloc[0]
ax.annotate(f'Rpot = {p.rpot_ohm:.1f} Ω\nλ ≈ {p.lambda_rosenstein_per_s/1000:.2f} × 10³ s⁻¹',xy=(p.rpot_ohm,p.lambda_rosenstein_per_s/1000),xytext=(825,2.35),fontsize=18,arrowprops=dict(arrowstyle='->',color=INK))
footer(fig,'Error bars show reported fit uncertainty; systematic estimator bias is additional. Positive values indicate local divergence.')
save(fig,'06_lyapunov','Quantitative evidence for chaos','forward_lyapunov.csv; direct Rosenstein estimates only')

identified=json.loads((ROOT/'chua/identified.json').read_text())
ind=pd.DataFrame(identified['inductor']['per_record'])
fig=canvas('The inductor changes with oscillation amplitude','Effective inductance and loss inferred from periodic measured orbits')
for j,(key,label) in enumerate([('L_eff_mH','Effective inductance (mH)'),('r_eff_ohm','Effective resistance (Ω)')]):
    ax=fig.add_axes([.10+j*.46,.19,.37,.61]); ax.scatter(ind.I_amp_mA,ind[key],s=30,color=BLUE,alpha=.75)
    ax.set(xlabel='Inductor current amplitude (mA)',ylabel=label); style(ax)
footer(fig,'These are effective values per orbit. The bench model uses a current-dependent Rayleigh law, plus amplifier dynamics.')
save(fig,'07_identified_inductor','Physical reason to improve the model','identified.json: inductor.per_record')

fig=canvas('The identified model improves the transition resistances','Selected thresholds from the project’s continuation-sweep analysis')
ax=fig.add_axes([.35,.21,.58,.56])
labels=['First period doubling','Double-scroll onset\n(decreasing R)','Double-scroll end\n(decreasing R)','Large-cycle upper limit\n(increasing R)']
vals=np.array([[774.6,772,776],[668,579,697],[328,519,344],[675,769,683]])
for i,row in enumerate(vals): ax.plot([min(row),max(row)],[i,i],color='#c8cfd5',lw=2,zorder=0)
for j,(col,mark,lab) in enumerate([(BLUE,'o','Measurements'),('#919aa3','s','Ideal model'),(ORANGE,'D','Identified model')]):
    ax.scatter(vals[:,j],np.arange(4),s=100,c=col,marker=mark,label=lab,zorder=3)
ax.set_yticks(range(4),labels); ax.set(xlim=(290,820),ylim=(3.5,-.65),xlabel='Transition Rpot (Ω)'); style(ax)
ax.legend(loc='upper center',bbox_to_anchor=(.37,1.14),ncol=3)
footer(fig,'Ideal: C1 = 11.5 nF. Bench: identified components with inductor and amplifier dynamics. Values from RESULTS.md.')
save(fig,'08_model_transitions','Quantitative model comparison','RESULTS.md: Models against the bench; rounded continuation thresholds')

if '--animation' not in sys.argv and (OUT/'manifest.json').exists():
    manifest.extend(item for item in json.loads((OUT/'manifest.json').read_text()) if item['file'].startswith(('09_', '10_')))
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

# Reuse the existing exact-resistance trajectory pipeline, simplifying its renderer.
if '--animation' in sys.argv:
    import animate_resistance as sweep
    def render(measured,modeled,rows,full_limits,model='bench'):
        fig=plt.figure(figsize=(12.8,7.2),dpi=90)
        fig.text(.06,.93,'Measurements and identified model at equal R',fontsize=26,weight='bold')
        label=fig.text(.06,.86,'',fontsize=20)
        axes=[fig.add_axes([.09,.20,.37,.55]),fig.add_axes([.58,.20,.37,.55])]
        lines=[]
        for ax,title,col in zip(axes,['Measurements','Identified bench model'],[BLUE,ORANGE]):
            ax.set(title=title,xlabel='V1 (V)',ylabel='V2 (V)'); style(ax)
            line,=ax.plot([],[],lw=.6,color=col,alpha=.8); lines.append(line)
        footer(fig,'20 ms per portrait. Downward continuation sweep; 100 ms settling per simulated resistance.')
        note=fig.text(.5,.085,'',ha='center',fontsize=14,color='#586673')
        def update(k):
            r=rows[k][1]; wide=r<360
            label.set_text(f'Rpot = {r:.2f} Ω     Rtotal = {r+992:.2f} Ω')
            for ax,line,collection in zip(axes,lines,[measured,modeled]):
                xy=collection[k]; line.set_data(xy[:,0],xy[:,1])
                ax.set_xlim((-8,8) if wide else (-5,3)); ax.set_ylim((-7,7) if wide else (-1.15,1.15))
            note.set_text('Full voltage range' if wide else 'Small-attractor view · both panels have identical scales')
        anim=FuncAnimation(fig,update,frames=len(rows),blit=False)
        anim.save(OUT/'09_equal_r_bench.gif',writer=PillowWriter(fps=8))
        for target,name in [(550,'09_equal_r_bench_still'),(700,'10_model_onset_difference')]:
            update(min(range(len(rows)),key=lambda i:abs(rows[i][1]-target)))
            save(fig,name,'Paused model comparison',f'Identified bench continuation; near {target} Ω')
        print('Wrote comparison animation',flush=True)
    sweep.render_comparison=render
    sys.argv=['animate_resistance.py','--model','bench','--comparison-only']
    sweep.main()
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
