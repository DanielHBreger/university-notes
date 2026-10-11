"""Journal-style figures for the two-column lab report.

Figures are drawn at the size they are printed, so 9 pt text in Python is 9 pt on
the page (the report body is 10 pt). Include them in LaTeX with width=\\columnwidth (default) or, for a
figure* spanning both columns, width=\\textwidth. If the template changes, print
\\the\\columnwidth and \\the\\textwidth in the document and update the two widths.

Text is typeset by LaTeX when it is installed, with LATEX_PREAMBLE chosen to match
the report's body font; otherwise the STIX fonts bundled with matplotlib are used.
"""
import shutil
from pathlib import Path

import matplotlib as mpl

PT = 1 / 72.27                 # TeX points per inch
COLUMN_WIDTH = 251.0 * PT      # \columnwidth: A4, 1.5 cm margins, twocolumn, 10 pt columnsep
TEXT_WIDTH = 512.1 * PT        # \textwidth, for figure* spanning both columns
DPI = 600                      # raster previews; the PDF is vector
# The report sets TeX Gyre Pagella for text and math; this is the pdflatex equivalent
LATEX_PREAMBLE = r'\usepackage[T1]{fontenc}\usepackage{tgpagella}\usepackage{newpxmath}'

# Okabe-Ito, distinguishable for colour-blind readers and in greyscale by marker
COLORS = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#000000']


def apply():
    usetex = shutil.which('latex') is not None
    mpl.rcParams.update({
        'text.usetex': usetex,
        'text.latex.preamble': LATEX_PREAMBLE,
        'font.family': 'serif',
        'font.serif': ['TeX Gyre Pagella', 'Palatino Linotype', 'STIXGeneral', 'DejaVu Serif'],
        'mathtext.fontset': 'stix',
        'font.size': 9, 'axes.labelsize': 9, 'axes.titlesize': 9,
        'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8,
        # frame and ticks: inward, on all four sides, with minor ticks
        'axes.linewidth': 0.5, 'axes.grid': False,
        'xtick.direction': 'in', 'ytick.direction': 'in',
        'xtick.top': True, 'ytick.right': True,
        'xtick.minor.visible': True, 'ytick.minor.visible': True,
        'xtick.major.size': 3.5, 'ytick.major.size': 3.5,
        'xtick.minor.size': 1.8, 'ytick.minor.size': 1.8,
        'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
        'xtick.minor.width': 0.4, 'ytick.minor.width': 0.4,
        'xtick.major.pad': 3, 'ytick.major.pad': 3,
        'axes.labelpad': 2,
        'axes.prop_cycle': mpl.cycler(color=COLORS),
        'lines.linewidth': 0.9, 'lines.markersize': 2.5, 'lines.markeredgewidth': 0.5,
        'legend.frameon': False, 'legend.handlelength': 1.4, 'legend.handletextpad': 0.4,
        'legend.borderaxespad': 0.3, 'legend.labelspacing': 0.25,
        'savefig.dpi': DPI, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.01,
        'pdf.fonttype': 42, 'ps.fonttype': 42,
    })


def figsize(full_width=False, aspect=0.75):
    """(width, height) in inches; aspect is height/width."""
    w = TEXT_WIDTH if full_width else COLUMN_WIDTH
    return (w, w * aspect)


def save(fig, path):
    """Save as vector PDF for LaTeX and as a 600 dpi PNG for previews."""
    path = Path(path).with_suffix('')
    fig.savefig(path.with_suffix('.pdf'))
    fig.savefig(path.with_suffix('.png'))
