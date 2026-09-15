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
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

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
UNC = {s: pd.read_csv(ROOT/'uncertainty'/f'{s}_uncertainty.csv').set_index('filename') for s in ['forward','back']}
DU = json.loads((ROOT/'uncertainty/diode_uncertainty.json').read_text())

def uncertainty_row(src):
    sweep, name = src.split('/')
    return UNC[sweep].loc[name]

def error_caption(ax, src):
    u=uncertainty_row(src)
    ax.text(.5,-.24,f"R fit: ±{u.u_r_block_ohm:.2f} Ω",
            transform=ax.transAxes,ha='center',fontsize=10,color='#586673')

def canvas(title, subtitle=''):
    fig = plt.figure(figsize=(40/3,7.5))
    fig.text(.075,.93,title,fontsize=27,weight='bold')
    if subtitle: fig.text(.075,.875,subtitle,fontsize=15,color='#586673')
    return fig

def footer(fig, text):
    fig.text(.075,.035,text,fontsize=10,color='#586673')

def style(ax): ax.grid(alpha=.15); ax.set_axisbelow(True)

def save(fig, name, use, source):
    fig.savefig(OUT / (name+'.png'),dpi=144)
    fig.savefig(OUT / (name+'.svg'))
    svg=OUT/(name+'.svg')
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
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
fig=canvas('A simple circuit produces a double-scroll attractor')
ax=fig.add_axes([.10,.12,.80,.68]); ax.plot(d[:,0],d[:,1],lw=.5,color=BLUE,alpha=.65)
ax.axis('off'); fig.text(.5,.035,f'Measured double scroll   •   Rpot = {r:.1f} Ω',ha='center',fontsize=16)
save(fig,'01_title_double_scroll','Title slide image',src+'; first 50 ms')

# Static I–V data, retaining original offset so the fit and measurements match.
cfg=json.loads((ROOT/'diode_fit.json').read_text())
t,vi=read_scope(ROOT/'trace1.csv')
v=vi[:,1]
current=(vi[:,0]-v)/cfg['shunt_ohm']*1000
fig=canvas('The measured I–V curve defines the nonlinear model','Measured current–voltage characteristic and the five-segment fit')
ax=fig.add_axes([.11,.17,.82,.62])
ax.scatter(v[::4],current[::4],s=3,color=BLUE,alpha=.22,rasterized=True,label='Measurements')
for j,s in enumerate(cfg['segments']):
    x=np.linspace(s['v_lo'],s['v_hi'],100)
    ax.plot(x,1000*(s['slope_S']*x+s['intercept_A']),color=ORANGE,lw=3,label='Segment fits' if j==0 else None)
ax.axhline(0,color='#77838b',lw=.7)
ax.set(xlabel='Voltage across nonlinear element (V)',ylabel='Current (mA)',xlim=(-9,8),ylim=(-4.5,4)); style(ax)
ax.legend(markerscale=3)

for j,seg in enumerate(DU['segments']):
    xx=np.linspace(seg['v_lo'],seg['v_hi'],100); X=np.column_stack([xx,np.ones(100)])
    up=np.sqrt(np.einsum('ni,ij,nj->n',X,np.array(seg['cov_slope_intercept']),X))*1000
    yy=(seg['slope_S']*xx+seg['intercept_A'])*1000
    ax.fill_between(xx,yy-up,yy+up,color=ORANGE,alpha=.20)
cov=np.array(DU['cov_VI_I_quantization'])*np.outer([1,1000],[1,1000]); val,vec=np.linalg.eigh(cov)
angle=np.degrees(np.arctan2(vec[1,-1],vec[0,-1]))
footer(fig,'Shaded bands show the uncertainty of the segment fits (one standard deviation).')
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
    style(ax); error_caption(ax,src)
footer(fig,'All panels use the same voltage scales. Resistance error labels show one standard deviation of the divider fit.')
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
a.add_patch(Rectangle((735,-1.1),70,1.45,fill=False,edgecolor=ORANGE,
                      linewidth=2.3,linestyle='--',zorder=6))
a.text(715,-1.7,'Zoom at right',color=ORANGE,fontsize=14)
for spine in b.spines.values():
    spine.set_visible(True); spine.set_color(ORANGE); spine.set_linestyle('--'); spine.set_linewidth(1.5)
for r,txt,y in [(774.6,'1 → 2',.23),(755.6,'2 → 4',.08),(750.5,'4 → 8',-.07)]:
    b.axvline(r,color=ORANGE,lw=1,alpha=.6); b.text(r+1,y,txt,ha='right',fontsize=14,color=ORANGE)
footer(fig,'R decreases to the right. Transition intervals: 774.55–777.62 Ω, 754.84–756.35 Ω and 750.29–751.11 Ω.')
save(fig,'04_bifurcation','Bifurcation and period doubling','forward_bifurcation_points.csv; cascade_periods_forward.txt')

fig=canvas('The observed attractor depends on sweep direction','Different orbits persist over an overlapping range of resistance')
a=fig.add_axes([.09,.20,.46,.60])
for df,col,lab in [(fwd,BLUE,'Decreasing R'),(back,PURPLE,'Increasing R')]:
    a.scatter(df.rpot_ohm,df.max_v,s=2.8 if lab=='Increasing R' else .7,
              marker='^' if lab=='Increasing R' else 'o',color=col,
              alpha=.8 if lab=='Increasing R' else .4,lw=0,rasterized=True,label=lab)
a.set(xlim=(900,0),ylim=(-3.5,8),xlabel='Rpot (Ω)',ylabel='Local maxima of V1 (V)'); style(a)
a.legend(handles=[Line2D([],[],color=c,marker=m,linestyle='none',markersize=10,label=l)
                  for c,m,l in [(BLUE,'o','Decreasing R'),(PURPLE,'^','Increasing R')]],loc='upper left')
for r in [328,675]: a.axvline(r,color='#6d7781',ls='--',lw=1)
for j,(sweep,col,lab) in enumerate([('forward',BLUE,'Decreasing R'),('back',PURPLE,'Increasing R')]):
    t,d,r,src=record(550,sweep)
    ax=fig.add_axes([.66,.55-j*.36,.30,.24]); ax.plot(d[:,0],d[:,1],color=col,lw=.65)
    ax.set(xlim=(-8,8),ylim=(-7,7),xlabel='V1 (V)',ylabel='V2 (V)',title=f'{lab}: {r:.1f} Ω'); style(ax)
footer(fig,'Resistance fit uncertainty: ±0.26 Ω (down) and ±0.05 Ω (up). Dashed lines mark 328 Ω and 675 Ω.')
save(fig,'05_hysteresis','History dependence at nearly the same R','forward/back bifurcation CSVs; forward trace270 and back trace21')

ly=pd.read_csv(ROOT/'chua/forward_lyapunov.csv')
ly=ly[(ly.rosenstein_status=='ok')&~ly.filename.isin([f'trace{i}.csv' for i in range(14,19)])]
fig=canvas('Positive Lyapunov estimates in the chaotic regime','Nearby trajectories separate in the chaotic range')
ax=fig.add_axes([.10,.20,.47,.58]); detail=fig.add_axes([.69,.20,.27,.58])
for panel, data, lims in [(ax,ly,((900,0),(-.3,2.8))),
                          (detail,ly[ly.rpot_ohm.between(540,560)],((562,539),(.90,1.65)))]:
    panel.axhline(0,color=INK,lw=1.5)
    panel.errorbar(data.rpot_ohm,data.lambda_rosenstein_per_s/1000,
                   xerr=data.filename.map(UNC['forward'].u_r_block_ohm),
                   yerr=data.lambda_rosenstein_err_per_s/1000,fmt='o',
                   ms=4 if panel is ax else 5,color=BLUE,ecolor=INK,
                   elinewidth=1.1 if panel is ax else 1.7,capsize=2 if panel is ax else 4,
                   capthick=1.3,alpha=.9)
    panel.set(xlim=lims[0],ylim=lims[1],xlabel='Rpot (Ω)'); style(panel)
ax.set_ylabel('Largest Lyapunov estimate (10³ s⁻¹)')
ax.set_title('Full resistance sweep'); detail.set_title('Error bars: enlarged view')
ax.add_patch(Rectangle((540,.9),20,.75,fill=False,edgecolor=ORANGE,linewidth=2,linestyle='--'))
detail.set_xticks([560,550,540])
footer(fig,'Error bars: one standard deviation of the fits. The enlarged view shows 540–560 Ω with unchanged error values.')
save(fig,'06_lyapunov','Quantitative evidence for chaos','forward_lyapunov.csv; direct Rosenstein estimates only')

identified=json.loads((ROOT/'chua/identified.json').read_text())
ind=pd.DataFrame(identified['inductor']['per_record'])
fig=canvas('The inductor changes with oscillation amplitude','Inductance and loss estimated from the measured oscillations')
for j,(key,label) in enumerate([('L_eff_mH','Effective inductance (mH)'),('r_eff_ohm','Effective resistance (Ω)')]):
    ax=fig.add_axes([.10+j*.46,.19,.37,.61]); ax.scatter(ind.I_amp_mA,ind[key],s=30,color=BLUE,alpha=.75)
    ax.set(xlabel='Inductor current amplitude (mA)',ylabel=label); style(ax)
footer(fig,'Each point is an estimate from one recording. The improved model includes the inductor’s changing behavior.')
save(fig,'07_identified_inductor','Physical reason to improve the model','identified.json: inductor.per_record')

fig=canvas('The identified model improves the transition resistances','Measured transition points compared with the two circuit models')
ax=fig.add_axes([.35,.21,.58,.56])
labels=['First period doubling','Double-scroll onset\n(decreasing R)','Double-scroll end\n(decreasing R)','Large-cycle upper limit\n(increasing R)']
vals=np.array([[774.6,772,776],[668,579,697],[328,519,344],[675,769,683]])
for i,row in enumerate(vals): ax.plot([min(row),max(row)],[i,i],color='#c8cfd5',lw=2,zorder=0)
for j,(col,mark,lab) in enumerate([(BLUE,'o','Measurements'),('#919aa3','s','Ideal model'),(ORANGE,'D','Identified model')]):
    ax.scatter(vals[:,j],np.arange(4),s=100,c=col,marker=mark,label=lab,zorder=3)

for i, target in enumerate(vals[:,0]):
    tab=UNC['back' if i==3 else 'forward']; item=tab.iloc[np.argmin(np.abs(tab.rpot_ohm-target))]
    ax.errorbar(target,i,xerr=item.u_r_block_ohm,fmt='none',color=BLUE,capsize=4,lw=1.5,zorder=4)
ax.set_yticks(range(4),labels); ax.set(xlim=(290,820),ylim=(3.5,-.65),xlabel='Transition Rpot (Ω)'); style(ax)
ax.legend(loc='upper center',bbox_to_anchor=(.37,1.14),ncol=3)
footer(fig,'Measured error bars show the resistance fit uncertainty (one standard deviation). Model points show the predicted values.')
save(fig,'08_model_transitions','Quantitative model comparison','RESULTS.md: Models against the bench; rounded continuation thresholds')

# A second view makes the model comparison readable without following connected dots.
fig=canvas('Measured components improve the predicted transitions',
           'Distance from the measured resistance: shorter bars mean better agreement')
ax=fig.add_axes([.35,.20,.59,.59])
deviations=np.abs(vals[:,1:]-vals[:,[0]])
for j,(color,label,offset) in enumerate([('#919aa3','Ideal model',-.17),
                                        (ORANGE,'Model with measured components',.17)]):
    yy=np.arange(4)+offset
    ax.barh(yy,deviations[:,j],height=.28,color=color,label=label)
    for y,value in zip(yy,deviations[:,j]):
        ax.text(value+3,y,f'{value:g} Ω',va='center',fontsize=16,color=INK)
ax.set_yticks(range(4),labels)
ax.set(xlim=(0,225),ylim=(3.55,-.6),xlabel='Distance from measured transition (Ω)')
ax.grid(axis='x',alpha=.15); ax.set_axisbelow(True)
ax.legend(loc='upper center',bbox_to_anchor=(.40,1.13),ncol=2,fontsize=14)
footer(fig,'Uses the same rounded transition values as the original comparison. Bars show model differences, not uncertainty.')
save(fig,'12_model_transitions_deviation','Simplified model comparison; alternative to 08',
     'Absolute differences of the same rounded thresholds plotted in 08_model_transitions')

if '--animation' not in sys.argv and (OUT/'manifest.json').exists():
    manifest.extend(item for item in json.loads((OUT/'manifest.json').read_text()) if item['file'].startswith(('09_', '10_')))
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

# Validated period-1,2,4,8 records, with peak sequences to make close branches visible.
from cascade_periods import load_sidecar, period_of, lag_distances
period_records=load_sidecar('forward', str(ROOT/'chua'))
fig=canvas('Period doubling creates orbits with 1, 2, 4 and 8 cycles',
           'Measured phase portraits above; 16 successive voltage peaks below')
chosen=[(1,'trace58.csv'),(2,'trace82.csv'),(4,'trace90.csv'),(8,'trace92.csv')]
for j,(period,name) in enumerate(chosen):
    rp,peaks=period_records[name]
    found,*_=period_of(lag_distances(peaks))
    assert found==period, (name,found,period)
    tt,dd=read_scope(ROOT/'chua'/'forward'/name)
    start=tt[0]+.25*(tt[-1]-tt[0]); use=(tt>=start)&(tt<=start+.02)
    ax=fig.add_axes([.105+j*.222,.45,.195,.32])
    ax.plot(dd[use,0],dd[use,1],lw=.5,color=BLUE,alpha=.7)
    ax.set(xlim=(-4.5,.5),ylim=(-1,1),xlabel='V1 (V)',title=f'Period {period}\n{rp:.2f} Ω')
    ax.set_xticks([-4,-2,0]);ax.set_yticks([-1,0,1]);style(ax)
    if j==0: ax.set_ylabel('V2 (V)')
    else: ax.tick_params(labelleft=False)
    bx=fig.add_axes([.105+j*.222,.15,.195,.17])
    values=peaks[len(peaks)//4:len(peaks)//4+16]
    bx.plot(np.arange(1,17),values,'o-',ms=4,lw=.8,color=BLUE)
    bx.set(xlabel='Peak number',ylim=(-1,.15),xlim=(.5,16.5))
    bx.set_xticks([1,8,16]);style(bx)
    if j==0: bx.set_ylabel('V1 peak (V)',fontsize=14)
    else: bx.tick_params(labelleft=False)
footer(fig,'Peak heights come from the smoothed traces used in the bifurcation analysis. All columns share the same scales.')
save(fig,'11_period_doubling_orbits','Measured 1,2,4,8 cycle comparison','forward traces 58,82,90,92; period labels verified by lag distances')
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

# Reuse the existing exact-resistance trajectory pipeline, simplifying its renderer.
if '--animation' in sys.argv:
    import animate_resistance as sweep
    def render(measured,modeled,rows,full_limits,model='bench'):
        fig=plt.figure(figsize=(40/3,7.5),dpi=96)
        fig.text(.06,.93,'The circuit model follows the measured attractors',fontsize=26,weight='bold')
        label=fig.text(.06,.86,'',fontsize=20)
        axes=[fig.add_axes([.09,.20,.37,.55]),fig.add_axes([.58,.20,.37,.55])]
        lines=[]
        for ax,title,col in zip(axes,['Measurements','Model with measured components'],[BLUE,ORANGE]):
            ax.set(title=title,xlabel='V1 (V)',ylabel='V2 (V)'); style(ax)
            line,=ax.plot([],[],lw=.6,color=col,alpha=.8); lines.append(line)
        footer(fig,'20 ms per portrait. The model settles for 100 ms at each resistance. R error labels show the fit uncertainty.')
        note=fig.text(.5,.085,'',ha='center',fontsize=14,color='#586673')
        def update(k):
            r=rows[k][1]; wide=r<360
            u=UNC['forward'].loc[rows[k][0]]
            label.set_text(f'Rpot = {r:.2f} Ω     Fit uncertainty: ±{u.u_r_block_ohm:.2f} Ω')
            for ax,line,collection in zip(axes,lines,[measured,modeled]):
                xy=collection[k]; line.set_data(xy[:,0],xy[:,1])
                ax.set_xlim((-8,8) if wide else (-5,3)); ax.set_ylim((-7,7) if wide else (-1.15,1.15))
            view='Full voltage range' if wide else 'Small-attractor view'
            note.set_text(f'{view}. Both panels use the same voltage scales.')
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
