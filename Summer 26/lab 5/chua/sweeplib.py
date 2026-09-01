"""
Shared plumbing for the Chua sweep scripts.

Everything here is about finding and naming files, not about the physics:
where a sweep folder is, which records are in it, and where the artefacts
derived from it live. Kept to the standard library so the measurement stage
(batch_rpot.py) can use it without pulling in numpy, scipy or matplotlib.

The chain writes its artefacts BESIDE the sweep folder, named after it:

    sweep forward/            the records
    sweep forward_rpot.csv    batch_rpot.py
    sweep forward_lorenz/     lorenz_map.py
    sweep forward_bifurcation.png, _lyapunov.png

`sibling` is the one place that convention is spelled out, so the producer
and the two scripts that read the sidecar back cannot drift apart.
"""
import csv
import os
import re
import sys

# rpot.py calls a divider fit above this residual "not a clean divider".
# batch_rpot.py writes that verdict out; the plotting scripts re-derive it
# from the numeric residual, so the threshold belongs in one place.
CLEAN_DIVIDER_MAX_PCT = 5.0


def natural_key(name):
    """Sort trace2 before trace10."""
    return [int(t) if t.isdigit() else t.lower()
            for t in re.split(r'(\d+)', name)]


def script_dir():
    """Where the running script lives - the sweep folders sit beside it."""
    d = os.path.dirname(os.path.abspath(sys.argv[0]))
    return d or os.getcwd()


def sibling(folder, suffix):
    """Path of an artefact derived from `folder`, named after it."""
    return os.path.join(os.path.dirname(folder),
                        os.path.basename(folder) + suffix)


def list_csvs(folder, recursive=False, full=False):
    """
    The records in a sweep folder, in natural order.

    Names by default, full paths with `full`. Artefacts this chain writes are
    named "<folder>_*.csv" and live outside the folder, so nothing here has to
    exclude them - except a --recursive walk, which can reach them.
    """
    if recursive:
        found = [os.path.join(dirpath, n)
                 for dirpath, _, names in os.walk(folder)
                 for n in names if n.lower().endswith('.csv')]
        found.sort(key=natural_key)
        return found if full else [os.path.relpath(f, folder) for f in found]
    names = sorted((n for n in os.listdir(folder)
                    if n.lower().endswith('.csv')), key=natural_key)
    return [os.path.join(folder, n) for n in names] if full else names


def pick_folder(title='Folder of scope CSVs', initialdir=None):
    """Open a directory chooser; fall back to a prompt without a display."""
    initialdir = initialdir or script_dir()
    try:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        folder = filedialog.askdirectory(title=title, initialdir=initialdir)
        root.destroy()
    except Exception:
        folder = input(f'{title} (under {initialdir}): ').strip()
        # Bare names are read as sitting beside the script, like the dialog.
        if folder and not os.path.isabs(folder) and not os.path.isdir(folder):
            beside = os.path.join(initialdir, folder)
            if os.path.isdir(beside):
                folder = beside
    return folder


def pick_folders():
    """
    Pick any number of sweep folders, one dialog at a time.

    The first opens beside the script; each later one opens beside the folder
    just chosen, since sweeps of the same run live together. Cancel to finish.
    """
    picked, start = [], script_dir()
    while True:
        title = ('Another sweep folder (Cancel when done)' if picked
                 else 'Sweep folder')
        f = pick_folder(title, start)
        if not f:
            break
        f = os.path.abspath(os.path.expanduser(f))
        if not os.path.isdir(f):
            print(f'not a folder, skipped: {f}')
            continue
        if f in picked:
            print(f'already added: {os.path.basename(f)}')
            continue
        picked.append(f)
        start = os.path.dirname(f)
        print(f'added {os.path.basename(f)}')
    return picked


def resolve_folders(named, multi=True):
    """Named folders if given, otherwise ask; exits if none survive."""
    if named:
        folders = [os.path.abspath(os.path.expanduser(f)) for f in named if f]
    elif multi:
        folders = pick_folders()
    else:
        one = pick_folder()
        folders = [os.path.abspath(os.path.expanduser(one))] if one else []
    if not folders:
        sys.exit('no folder selected')
    for f in folders:
        if not os.path.isdir(f):
            sys.exit(f'not a folder: {f}')
    return folders


def load_rpot(folder, required=False):
    """
    filename -> (Rpot in ohm, fit residual in %) from batch_rpot.py's output.

    Empty when there is no sidecar, unless `required` - the scripts that put
    Rpot on an axis cannot proceed without it, the ones that only label with
    it can.
    """
    path = sibling(folder, '_rpot.csv')
    if not os.path.exists(path):
        if required:
            sys.exit(f'need {path}\nrun:  python batch_rpot.py "{folder}"')
        return {}
    out = {}
    with open(path, newline='') as fh:
        for row in csv.DictReader(fh):
            try:
                out[row['filename']] = (float(row['rpot_ohm']),
                                        float(row['residual_pct']))
            except (KeyError, ValueError):
                continue
    return out


def sweep_colors(n):
    """One distinct colour per sweep folder, for any number of them."""
    from matplotlib import pyplot as plt  # kept local: this module is stdlib
    cmap = plt.get_cmap('tab10' if n <= 10 else 'tab20')
    return [cmap(i % cmap.N) for i in range(n)]
