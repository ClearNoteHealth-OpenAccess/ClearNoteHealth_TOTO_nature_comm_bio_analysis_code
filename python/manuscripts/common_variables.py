import pickle, os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib as mpl
import scipy as sp
import seaborn as sns
import statsmodels.api as sm  
from nheatmap import nhm
from DensityPlot import density2d
from matplotlib_venn import venn2

def load_drg(fn, drg, group='', sep=','):
    drg[group] = pd.read_csv(fn, sep=sep, index_col=0)
    drg[group]['sig'] = (drg[group]['FDR'] < 0.05) * (drg[group]['FC'] != 0)
    drg[group]['kind'] = 'not significant'
    drg[group].loc[np.logical_and(drg[group]['sig'], drg[group]['FC'] > 0),
                   'kind'] = 'up-regulated'
    drg[group].loc[np.logical_and(drg[group]['sig'], drg[group]['FC'] < 0),
                   'kind'] = 'down-regulated'
    drg[group]['color'] = drg[group]['kind'].map(colors)

translate = {
    'pancreas': 'pancreatic',
    'breast': 'breast',
    'all': 'all',
    'colon': 'colorectal',
    'ovary': 'ovarian',
    'all': 'all',
    'lung': 'lung'
}

rev_translate = {y: x for x, y in zip(translate.keys(), translate.values())}

colors = {
    'not significant': '#969696',
    'up-regulated': '#B31B21',
    'down-regulated': '#1465AC'
}
cmap_primary = {
    x: 'C{:}'.format(i)
    for i, x in enumerate(['breast', 'colon', 'lung', 'ovary', 'pancreas', 'all'])
}
cmaps = {
    'has_cancer': {
        'Cancer': 'red',
        'Healthy': 'blue'
    },
    'cohort': {
        'Tumor': '#B31B21',
        'Normal': '#1465AC'
    },
    'primary_diagnosis': cmap_primary,
    'tissue_type': {
        'NAT': 'forestgreen',
        'normal': 'limegreen',
        'tumor': 'salmon'
    }
}

venn_c = ('#B31B21', '#1465AC')

remap_stage_dt = {
    'IIIB': 'III',
    'IA': 'I',
    'IIA': 'II',
    'IIB': 'II',
    'I': 'I',
    'II': 'II',
    'III': 'III',
    'IV': 'IV',
    'NAT': 'Normal',
    'Healthy': 'Normal',
    'IVA': 'IV',
    'Unknown': 'Unknown',
    'IIIC': 'III',
    'IIIA': 'III',
    'IB': 'I'
}
remap_stage_nat_dt = remap_stage_dt.copy()
remap_stage_nat_dt['NAT'] = 'NAT'

def barplot_annotate_brackets(ax, x, y, text, dx=0.15, fs=12, barh=0.025, bar_pad=1, dh=0.05):
    """
    Annotate barplot with p-values.

    Parameters
    ----------
    ax : matplotlib axis
    x : array-like
        x-coordinates of the brackets
    y : array-like
        y-coordinates of the brackets
    text : array-like
    dx : float
        x offset of the text
    fs : float
        font size
    barh : float
        height of the bracket lines
    bar_pad : float
        distance between brackets and text
    dh : float
        height of the bracket lines
    """

    ax_y0, ax_y1 = ax.get_ylim()
    dh *= (ax_y1 - ax_y0)
    barh *= (ax_y1 - ax_y0)

    lx = x - dx
    rx = x + dx
    y = y + bar_pad

    barx = [lx, lx, rx, rx]
    bary = [y, y+barh, y+barh, y]
    mid = ((lx+rx)/2, y+barh)

    ax.plot(barx, bary, c='black')

    kwargs = dict(ha='center', va='bottom')
    if fs is not None:
        kwargs['fontsize'] = fs

    ax.text(*mid, text, **kwargs)
    
def add_identity(axes, ls='--', color='k', xlim=None, ylim=None, *line_args, **line_kwargs):
    if 'color' not in line_kwargs and 'c' not in line_kwargs:
        line_kwargs['color'] = color
    identity, = axes.plot([], [], ls=ls, *line_args, **line_kwargs)
    def callback(axes):
        if xlim is None:
            low_x, high_x = axes.get_xlim()
        else:
            low_x, high_x = xlim[0], xlim[1]
        if ylim is None:
            low_y, high_y = axes.get_ylim()
        else:
            low_y, high_y = ylim[0], ylim[1]
        low = max(low_x, low_y)
        high = min(high_x, high_y)
        identity.set_data([low, high], [low, high])
    callback(axes)
    axes.callbacks.connect('xlim_changed', callback)
    axes.callbacks.connect('ylim_changed', callback)
    return axes
    
def calculate_specificity(df):
    """
    Calculate the specificity for each class in a confusion matrix.

    Parameters
    ----------
    df : pandas.DataFrame
    confusion matrix. rows are predicted labels, columns are true labels.

    Returns
    -------
    The specificity for each class.
    """
    specificity = {}
    for label in df.index:
        total = df.sum().sum()
        tp = df.loc[label, label]
        fn = df[label].sum() - tp
        fp = df.loc[label].sum() - tp
        tn = total - tp - fn - fp
        specificity[label] = tn / (tn + fp)
    return pd.Series(specificity)

def prob_intersect_hypergeom(M, table, **args):
    """
    Perform a hypergeometric test to assess statistical significance of
    intersection between the elements of two vectors.

    Parameters
    ----------
        M: total length of the sample space (e.g. number of genes > 0 expression in the genome)
        table: 2 x 2 contingency table of DEG set size in two sets. True:True is the intersection size.
    """
    n = table.loc[True, True]
    N1 = table.sum(1)[True]
    N2 = table.sum(0)[True]
    N = N1 + N2 - n
    rv = sp.stats.hypergeom(M, n, N, **args)
    pval = rv.sf(n - 1)
    return (rv, pval)

_output_formats = [".pdf"]
_savefig_args = {
    'dpi': 500,
    'bbox_inches': 'tight',
    'pad_inches': 0.05,
    'transparent': True
}

def savefig(fig, name, rasterized=False):
    assert name[-1] != '/', 'cannot provide a directory name'
    for output_format in _output_formats:
        fig.savefig(name + output_format, **_savefig_args)
    return None

