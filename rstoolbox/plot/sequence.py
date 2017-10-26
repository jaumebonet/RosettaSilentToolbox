from rstoolbox.analysis import sequence_frequency_matrix, binary_overlap
from .color_schemes import color_scheme

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.patheffects
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, PathPatch
from matplotlib import transforms
from matplotlib.font_manager import FontProperties
from matplotlib.text import TextPath

class Scale(matplotlib.patheffects.RendererBase):
    def __init__(self, sx, sy=None):
        self._sx = sx
        self._sy = sy

    def draw_path(self, renderer, gc, tpath, affine, rgbFace):
        affine = affine.identity().scale(self._sx, self._sy)+affine
        renderer.draw_path(gc, tpath, affine, rgbFace)

def _letterAt(letter, x, y, yscale=1, ax=None, globscale=1.35, LETTERS=None, COLOR_SCHEME=None):
    text = LETTERS[letter]

    t = mpl.transforms.Affine2D().scale(1*globscale, yscale*globscale) + \
        mpl.transforms.Affine2D().translate(x,y) + ax.transData
    p = PathPatch(text, lw=0, fc=COLOR_SCHEME[letter],  transform=t)
    if ax != None:
        ax.add_artist(p)
    return p

def _dataframe2logo( data ):
    aa = list(data)
    odata = []
    for index, pos in data.iterrows():
        pdata = []
        for k in aa:
            if pos[k] > 0.0:
                #pdata.append( ( k, float(format(pos[k] * 2.23, '.2f')) ) )
                pdata.append( ( k, float(format(pos[k], '.2f')) ) )
                #pdata.append( ( k, float(pos[k]) ) )
        odata.append(sorted(pdata, key=lambda x: x[1]))
    return odata

def barcode_plot( df, column_name, axis, color="blue" ):
    result = binary_overlap( df, column_name )
    pd.Series(result).plot("bar", ax=axis, ylim=(0,1), grid=False, color=color, width=1 )
    axis.yaxis.set_ticks([])
    axis.xaxis.set_ticks(np.arange(0, len(result)+1, 10))
    axis.xaxis.set_ticklabels(np.arange(0, len(result)+1, 10) + 1, rotation=45)
    axis.set_xlabel("sequence")

def sequence_frequency_plot( df, column_name, axis, ref_seq=None, key_residues=None, colormap = "Blues", border_color="green", nobar=False, cbar_ax=None, orientation="horizontal" ):
    order = ["A","V","I","L","M","F","Y","W","S","T","N","Q","R","H","K","D","E","C","G","P"]
    data = sequence_frequency_matrix( df, column_name ).transpose().reindex(order)
    if key_residues is not None:
        data = data.ix[:, key_residues]
        if ref_seq is not None:
            tmp_seq = ""
            for k in key_residues:
                tmp_seq += ref_seq[k-1]
            ref_seq = tmp_seq
    if cbar_ax == None:
        sns.heatmap(data, ax=axis, square=True, cbar=not nobar, cbar_kws={"orientation": orientation}, linewidths=1, cmap=colormap)
    else:
        sns.heatmap(data, ax=axis, square=True, cbar_ax=cbar_ax, cbar_kws={"orientation": orientation}, linewidths=1, cmap=colormap)
    order.reverse()
    axis.yaxis.set_ticklabels(order, rotation=0)
    axis.xaxis.set_ticklabels(data.columns.values.tolist())
    axis.set_ylabel("residue type")
    if ref_seq is not None:
        for i in range(len(ref_seq)):
            axis.add_patch(Rectangle((i, order.index(ref_seq[i])), 1, 1, fill=False, edgecolor=border_color, lw=2))

def logo_plot( df, column_name, ref_seq=None, outfile=None, key_residues=None, colors="WEBLOGO" ):
    # Graphical Properties of resizable letters
    fp = FontProperties(family="DejaVu Sans Mono", weight="bold")
    globscale = 1.35
    letters_shift = -0.3

    LETTERS = {}
    for aa in color_scheme(colors):
        LETTERS[aa] = TextPath((letters_shift, 0), aa, size=1, prop=fp)

    data = sequence_frequency_matrix( df, column_name )
    if key_residues is not None:
        data = data.loc[key_residues]
        if ref_seq is not None:
            tmp_seq = ""
            for k in key_residues:
                tmp_seq += ref_seq[k-1]
            ref_seq = tmp_seq
    if ref_seq is not None:
        assert len(ref_seq) == len(data)

    fig, ax = plt.subplots(figsize=(len(data) * 2, 2.3 * 2))

    font = FontProperties()
    font.set_size(80)
    font.set_weight('bold')

    ax.set_xticks(range(1, len(data) + 2))
    ax.set_yticks( range(0, 2) )
    ax.set_xticklabels( data.index.values )
    ax.set_yticklabels( np.arange( 0, 2 , 1 ) )
    if ref_seq is not None:
        ax2 = ax.twiny()
        ax2.set_xticks(ax.get_xticks())
        ax2.set_xticklabels(list(ref_seq))
    sns.despine(ax=ax, trim=True)
    ax.grid(False)
    if ref_seq is not None:
        sns.despine(ax=ax2, top=False, right=True, left=True, trim=True)
        ax2.grid(False)
    ax.lines = []
    wdata = _dataframe2logo( data )
    all_scores = wdata
    x = 1
    maxi = 0
    for scores in all_scores:
        y = 0
        for base, score in scores:
            _letterAt(base, x,y, score, ax, globscale, LETTERS, color_scheme(colors))
            y += score
        x += 1
        maxi = max(maxi, y)
    if outfile is not None:
        fig.savefig( outfile )
    return fig, ax

def per_residue_value( df ):
    pass
