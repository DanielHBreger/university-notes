"""Rebuild sweep figures using the peaks just audited from every raw record.

This harness calls the same public command-line entry points. Its temporary
in-memory peak provider avoids parsing the 14.6 GB dataset repeatedly. It
never changes the analysis modules' normal file-reading behavior on disk.
"""
from pathlib import Path
import sys,json,copy,contextlib
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'chua'))
import lorenz_map
AUDIT=json.loads((ROOT/'verification/data_audit.json').read_text())
INDEX={str((ROOT/r['path']).resolve()):r for r in AUDIT if 'peak_cache' in r}
original_maxima=lorenz_map.maxima


def cached_maxima(path,ch='CH1',prominence=.02,period=None):
    r=INDEX.get(str(Path(path).resolve()))
    if r and ch=='CH1' and prominence==.02 and period is None:
        data=np.load(ROOT/'verification/cache'/f"{r['peak_cache']}.npz")
        info=copy.deepcopy(r['peak_info']); info['t_peaks']=data['t_peaks'].copy()
        return data['M'].copy(),info
    return original_maxima(path,ch,prominence,period)


lorenz_map.maxima=cached_maxima
import bifurcation,lyapunov


def run(module,args,name):
    print(f'Building {name}',flush=True)
    with (ROOT/'verification'/f'{name}.log').open('w',encoding='utf-8',buffering=1) as fh:
        with contextlib.redirect_stdout(fh):
            sys.argv=[module.__file__,*[str(a) for a in args]]
            module.main()


def main():
    folders=[ROOT/'chua'/s for s in ['sweep forward','sweep back','set 2/sweep forward','set 2/sweep back']]
    for i,folder in enumerate(folders):
        run(lorenz_map,[folder],f'lorenz_{i+1}')
        run(bifurcation,[folder],f'bifurcation_{i+1}')
    for i,pair in enumerate([folders[:2],folders[2:]]):
        run(bifurcation,[*pair,'-o',pair[0].parent/'hysteresis_bifurcation.png'],f'hysteresis_{i+1}')
    for i,folder in enumerate(folders):
        options=['--rosenstein'] if i<2 else []
        run(lyapunov,[folder,*options],f'lyapunov_{i+1}')
    named=[ROOT/'chua'/name for name in [
        'limitcycle1 - 773.0 ohm.csv','limitcycle2 - 764.6 ohm.csv',
        'chaos - 716.9 ohm.csv','doublescroll start - 666.0 ohm.csv',
        'doublescroll 2 - 656.7ohm.csv','doublescroll end - 316.1ohm.csv']]
    run(lyapunov,[*named,'--rosenstein','--each','-o',ROOT/'chua/records_lyapunov.png'],'lyapunov_records')
    print('All existing sweep figure families rebuilt.',flush=True)


if __name__=='__main__': main()
