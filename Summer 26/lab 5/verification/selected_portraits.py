"""Show the phase portraits and complete time records behind M4/M5."""
from pathlib import Path
import sys,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'chua'))
from scope_data import read_scope
audit={r['path']:r for r in json.loads((ROOT/'verification/data_audit.json').read_text())}
examples=json.loads((ROOT/'verification/M4_M5_results.json').read_text())
fig,axes=plt.subplots(2,3,figsize=(14,7),layout='constrained')
for col,r in enumerate(examples):
    t,ch=read_scope(ROOT/r['path']); info=audit[r['path']]
    window=max(5,(info['period_samples']//20)|1)
    smooth=savgol_filter(ch,window,3,axis=0)
    step=max(1,len(t)//100000)
    ax=axes[0,col]
    ax.scatter(ch[::step,0],ch[::step,1],s=2,lw=0,color='0.5',alpha=.04,rasterized=True)
    ax.scatter(smooth[::step,0],smooth[::step,1],s=.8,lw=0,color='#0072B2',alpha=.1,rasterized=True)
    ax.set(xlabel='$V_1$ (V)',ylabel='$V_2$ (V)',
           title=f"{Path(r['path']).name} · Rpot = {r['rpot_ohm']:.1f} Ω\nDivider residual {r['residual_pct']:.2f}%; SG {window} samples")
    ax=axes[1,col]
    # Min/max envelopes retain short excursions over the full half-second.
    chunks=len(t)//5000
    n=chunks*5000
    lo=ch[:n,0].reshape(5000,chunks).min(axis=1)
    hi=ch[:n,0].reshape(5000,chunks).max(axis=1)
    tx=(t[:n].reshape(5000,chunks).mean(axis=1)-t[0])*1000
    ax.fill_between(tx,lo,hi,color='#0072B2',alpha=.75,lw=0)
    ax.set(xlabel='Time (ms)',ylabel='$V_1$ (V)',title='Full 500 ms record · min/max envelope')
for ax in axes.flat: ax.grid(alpha=.12)
fig.suptitle('M3 · Three long records used for M4 and M5\nPhase portraits: grey raw codes, blue light smoothing')
fig.savefig(ROOT/'verification/M3_selected_records.png',dpi=220)
plt.close(fig)
