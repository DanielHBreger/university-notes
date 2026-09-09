"""Audit every original scope record; cache peaks for reproducible figure rebuilds."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import sys, json, hashlib, csv
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'chua'))
from scope_data import read_scope, validate_time
from rpot import fit_divider
from lorenz_map import maxima_from_arrays
from sweeplib import is_scope_csv, natural_key
CACHE=ROOT/'verification/cache'; CACHE.mkdir(exist_ok=True)


def audit(path):
    relative=path.relative_to(ROOT).as_posix()
    result=dict(path=relative,status='ok',resistance_status='not measured')
    try:
        digest=hashlib.file_digest(path.open('rb'),'sha256').hexdigest()
        head=pd.read_csv(path,nrows=0)
        channels=tuple(c.split('(')[0] for c in head.columns[1:])
        t,y=read_scope(path,channels,min_samples=100)
        result.update(sha256=digest,n_samples=len(t),duration_s=float(t[-1]-t[0]),
                      dt_s=validate_time(t),n_channels=y.shape[1],
                      ch1_min=float(y[:,0].min()),ch1_max=float(y[:,0].max()),
                      ch2_min=float(y[:,1].min()),ch2_max=float(y[:,1].max()))
        if y.shape[1]>=3:
            try:
                r,res=fit_divider(y[:,:3])
                result.update(rpot_ohm=float(r),residual_pct=float(100*res),
                              resistance_status='ok' if 0<=r<=1000 and res<=.05 else 'CHECK')
            except ValueError as exc:
                result['resistance_status']=str(exc)
        if relative!='trace1.csv':
            M,info=maxima_from_arrays(t,y[:,0])
            key=hashlib.sha256(relative.encode()).hexdigest()[:24]
            np.savez_compressed(CACHE/f'{key}.npz',M=M,t_peaks=info.pop('t_peaks'))
            result.update(n_maxima=len(M),period_samples=info['period_samples'],
                          f0_hz=info['f0_hz'],quantum_V=info['quantum'],clip_pct=info['clip_pct'],
                          peak_status=info['status'],peak_cache=key,peak_info=info)
        return result
    except Exception as exc:
        result['status']=f'{type(exc).__name__}: {exc}'
        return result


def main():
    paths=[p for p in ROOT.rglob('*.csv') if 'verification' not in p.parts
           and '.venv' not in p.parts and is_scope_csv(p)]
    paths.sort(key=lambda p:natural_key(str(p)))
    results=[]
    with ThreadPoolExecutor(max_workers=3) as pool:
        for i,r in enumerate(pool.map(audit,paths),1):
            results.append(r)
            if i%20==0 or i==len(paths):
                print(f'{i}/{len(paths)} records audited',flush=True)
    (ROOT/'verification/data_audit.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    fields=['path','status','n_samples','duration_s','dt_s','n_channels','rpot_ohm','residual_pct',
            'resistance_status','n_maxima','period_samples','f0_hz','quantum_V','clip_pct','peak_status','sha256']
    with (ROOT/'verification/data_audit.csv').open('w',newline='') as fh:
        w=csv.DictWriter(fh,fields,extrasaction='ignore');w.writeheader();w.writerows(results)
    for folder in sorted({p.parent for p in paths if 'sweep ' in p.parent.name or p.parent.name=='comparisons'}):
        rows=[r for r in results if (ROOT/r['path']).parent==folder]
        out=folder.parent/(folder.name+'_rpot.csv')
        with out.open('w',newline='') as fh:
            w=csv.writer(fh);w.writerow(['filename','rpot_ohm','residual_pct','status'])
            for r in rows:
                w.writerow([Path(r['path']).name,round(r['rpot_ohm'],2) if 'rpot_ohm' in r else '',
                            round(r['residual_pct'],3) if 'residual_pct' in r else '',r['resistance_status']])
    errors=[r for r in results if r['status']!='ok']
    print(f'Finished: {len(results)} scope records; {len(errors)} structural errors.',flush=True)
    for r in errors: print(r['path'],r['status'])


if __name__=='__main__': main()
