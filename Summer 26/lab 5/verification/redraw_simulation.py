"""Redraw the saved numerical trajectories with explicit calibration labels."""
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'chua'))
from simulate import plot_nr,plot_portraits,plot_timeseries,find_records,make_g

for tag in ['', '_corrected']:
    folder=ROOT/'chua'
    with np.load(folder/f'simulated_runs{tag}.npz') as data:
        t=data['t']; y=data['states']; r=data['rpot']; offset=float(data['removed_offset_A'])
        runs=[(r[k],t,y[:,:,k],y[:,:,k+1]) for k in range(0,len(r),2)]
        note='Offset removed (assumption)' if offset else 'Measured current offset retained'
        plot_nr(make_g(i0=offset),r[::2],folder/f'simulated_nr{tag}.png',i0=offset,rL=float(data['rL']))
        plot_portraits(runs,find_records(folder),folder/f'simulated_portraits{tag}.png',200000,note)
        plot_timeseries(runs,folder/f'simulated_timeseries{tag}.png',model_note=note)
    print('Redrawn', tag or 'offset retained',flush=True)
