from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
from dual_agn.settings import FIGURE_ROOT

mpl.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'cm', 'font.size': 14,
    'axes.labelsize': 16, 'axes.titlesize': 16, 'axes.linewidth': 1.6,
    'xtick.direction': 'in', 'ytick.direction': 'in'})

def save_figure(figure, stem, group):
    folder = FIGURE_ROOT / group
    folder.mkdir(parents=True, exist_ok=True)
    for extension in ('png', 'pdf'):
        figure.savefig(folder / f'{stem}.{extension}', dpi=180, bbox_inches='tight')
    print('Figure:', folder / stem, flush=True)
