"""
Run rpot.py over every CSV in a folder and write the results to one CSV.

The folder can be given on the command line; with no argument a directory
picker opens beside this script (falls back to a prompt if Tk is unavailable).

Output columns: filename, rpot_ohm, residual_pct, status

Usage:
    python batch_rpot.py                              # pick the folder
    python batch_rpot.py "sweep forward"
    python batch_rpot.py "sweep forward" --g21 0.9760 -o forward.csv
    python batch_rpot.py "sweep forward" --recursive
"""
import argparse
import csv
import os

from rpot import rpot
from sweeplib import (CLEAN_DIVIDER_MAX_PCT, list_csvs, resolve_folders,
                      sibling)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('folder', nargs='?', help='folder of scope CSVs')
    p.add_argument('-o', '--out', help='output CSV (default: <folder>_rpot.csv)')
    p.add_argument('--recursive', action='store_true', help='include subfolders')
    p.add_argument('--r0', type=float, default=992.0, help='fixed resistor, ohm')
    p.add_argument('--g21', type=float, default=1.0,
                   help='gain ratio of the v_C2 channel to the v_C1 channel')
    p.add_argument('--c1', default='CH1', help='column with v_C1')
    p.add_argument('--c2', default='CH2', help='column with v_C2')
    p.add_argument('--mid', default='CH3', help='column with the midpoint')
    p.add_argument('--max-files', type=int, help='stop after this many records')
    a = p.parse_args()

    folder = resolve_folders([a.folder] if a.folder else [], multi=False)[0]
    files = list_csvs(folder, a.recursive, full=True)
    if a.max_files:
        files = files[:a.max_files]
    if not files:
        raise SystemExit(f'no CSV files in {folder}')

    out = a.out or sibling(folder, '_rpot.csv')
    # Never let the output land in the input set of a later re-run.
    files = [f for f in files if os.path.abspath(f) != os.path.abspath(out)]

    rows = []
    for i, f in enumerate(files, 1):
        name = os.path.relpath(f, folder)
        print(f'[{i}/{len(files)}] {name}', end=' ', flush=True)
        try:
            r, res = rpot(f, a.r0, a.g21, a.c1, a.c2, a.mid)
        except Exception as e:
            print(f'-> ERROR: {e}')
            rows.append({'filename': name, 'rpot_ohm': '', 'residual_pct': '',
                         'status': f'error: {e}'})
            continue
        pct = 100 * res
        status = ('ok' if pct <= CLEAN_DIVIDER_MAX_PCT
                  else 'CHECK: not a clean divider')
        print(f'-> {r:7.1f} ohm  ({pct:.2f} %)')
        rows.append({'filename': name, 'rpot_ohm': round(r, 2),
                     'residual_pct': round(pct, 3), 'status': status})

    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, ['filename', 'rpot_ohm', 'residual_pct', 'status'])
        w.writeheader()
        w.writerows(rows)

    bad = sum(1 for r in rows if r['status'] != 'ok')
    print(f'\nwrote {len(rows)} rows to {out}'
          + (f'   ({bad} flagged)' if bad else ''))


if __name__ == '__main__':
    main()
