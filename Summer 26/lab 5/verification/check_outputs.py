"""Check exported data against independently saved peaks and original hashes."""
from pathlib import Path
import hashlib,json
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
audit=json.loads((ROOT/'verification/data_audit.json').read_text())
index={r['path']:r for r in audit}
checks=[]
for folder in [ROOT/'chua'/p for p in ['sweep forward','sweep back','set 2/sweep forward','set 2/sweep back']]:
    stem=folder.parent/folder.name
    maxima=pd.read_csv(str(stem)+'_bifurcation_points.csv')
    summary=pd.read_csv(str(stem)+'_lorenz/lorenz_summary.csv',keep_default_na=False)
    pairs=pd.read_csv(str(stem)+'_lorenz/lorenz_pairs.csv')
    accepted=set(summary.loc[summary.status=='ok','filename'])
    assert set(maxima.filename)==accepted,folder
    assert set(pairs.filename)==accepted,folder
    for name,group in maxima.groupby('filename'):
        r=index[(folder/name).relative_to(ROOT).as_posix()]
        cached=np.load(ROOT/'verification/cache'/f"{r['peak_cache']}.npz")
        expected=cached['M']
        np.testing.assert_allclose(group.max_v,expected,atol=5.01e-7,rtol=0)
        np.testing.assert_allclose(group.rpot_ohm,r['rpot_ohm'],atol=.00501,rtol=0)
        actual=pairs[pairs.filename==name]
        np.testing.assert_allclose(actual.m_n,expected[:-1],atol=5.01e-7,rtol=0)
        np.testing.assert_allclose(actual.m_next,expected[1:],atol=5.01e-7,rtol=0)
    checks.append(dict(folder=folder.relative_to(ROOT).as_posix(),admitted_records=len(accepted),
                       maxima=len(maxima),return_pairs=len(pairs),status='passed'))


def unchanged(r):
    with (ROOT/r['path']).open('rb') as fh:
        return hashlib.file_digest(fh,'sha256').hexdigest()==r['sha256']


with ThreadPoolExecutor(max_workers=3) as pool:
    matches=list(pool.map(unchanged,audit))
assert all(matches),'A raw measurement changed after the audit'
result=dict(raw_records_checked=len(matches),all_raw_hashes_unchanged=all(matches),sweep_checks=checks)
(ROOT/'verification/output_checks.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
