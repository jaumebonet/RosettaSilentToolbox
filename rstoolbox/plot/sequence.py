import os

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

from rstoolbox.analysis import sequence_frequency_matrix, binary_overlap
from rstoolbox.components import DesignFrame, SequenceFrame
from .color_schemes import color_scheme

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
            if pos[k] > 0.0000000:
                pdata.append( ( k, float(pos[k]) ) )
        odata.append(sorted(pdata, key=lambda x: x[1]))
    return odata

def barcode_plot( df, column_name, ax, color="blue" ):
    result = binary_overlap( df, column_name )
    pd.Series(result).plot("bar", ax=ax, ylim=(0,1), grid=False, color=color, width=1 )
    ax.yaxis.set_ticks([])
    ax.xaxis.set_ticks(np.arange(0, len(result)+1, 10))
    ax.xaxis.set_ticklabels(np.arange(0, len(result)+1, 10) + 1, rotation=45)
    ax.set_xlabel("sequence")

def sequence_frequency_plot( df, seqID, ax, refseq=True, key_residues=None, border_color="green", **kwargs ):
    """
    Makes a heatmap subplot into the provided axis showing the sequence distribution of each residue type
    for each position. A part from the function arguments, any argument that can be provided to the
    seaborn.heatmap function can also be provided here.
    As a tip:
    (1) Do you want to set the orientation of the color bar horizontal?
    Add the parameter: cbar_kws={"orientation": "horizontal"}
    (2) Do you want to put the color bar in a different axis?
    Add the parameter: cbar_ax=[second_axis]
    (3) Yo don't want a color bar?
    Add the parameter: cbar=False

    :param DataFrame df: Ideally, a :py:class:`.DesignFrame` or :py:class:`.SequenceFrame`. Data content.
        requires the existence of a "sequence_{seqID}" column with the sequence to plot.
    :param str seqID: Identifier of the query sequence.
    :param axis ax: matplotlib axis to which we will plot.
    :param bool refseq: if True (default), mark the original residues according to the reference sequence.
    :param list key_residues: List to limit the plotted positions to those of interest.
    :param str border_color: Color to use to mark the original residue types.
    :raises: ValueError if input is not a DataFrame derived object.
    """

    order = ["A","V","I","L","M","F","Y","W","S","T","N","Q","R","H","K","D","E","C","G","P"]
    data = df

    # Data type management.
    if not isinstance(data, pd.DataFrame):
        raise ValueError("Input data must be in a DataFrame, DesignFrame or SequenceFrame")
    else:
        if not isinstance(data, (DesignFrame, SequenceFrame)):
            data = DesignFrame( data )
    if isinstance(data, DesignFrame):
        data = data.sequence_frequencies( seqID )
    if isinstance(data, SequenceFrame):
        if not data.is_transposed():
            data = data.transpose().reindex(order)
        else:
            data = data.reindex(order)

    # Refseq and key_residues management.
    ref_seq = data.key_reference_sequence(key_residues, False)
    if key_residues is not None:
        data = data[key_residues]

    # heatmap parameters
    if not "cmap" in kwargs:        # define the color-range of the plot
        kwargs["cmap"] = "Blues"
    kwargs["linewidths"] = 1        # linewidths are fixed to 1, overwrite user selection
    kwargs["square"] = True         # square is True, overwrite user selection
    if not "cbar_kws" in kwargs:    # by default the color bar is horizontal
        kwargs["cbar_kws"] = {"orientation": "horizontal"}

    # plot
    sns.heatmap(data, ax=ax, **kwargs)

    # styling plot
    order.reverse()
    ax.yaxis.set_ticklabels(order, rotation=0)
    ax.xaxis.set_ticklabels(data.columns.values.tolist())
    ax.set_ylabel("residue type")

    # marking reference sequence
    if ref_seq is not "" and refseq:
        for i in range(len(ref_seq)):
            ax.add_patch(Rectangle((i, order.index(ref_seq[i])), 1, 1, fill=False, edgecolor=border_color, lw=2))

def logo_plot( df, column_name, ref_seq=None, outfile=None, key_residues=None, colors="WEBLOGO" ):
    mpl.rcParams['svg.fonttype'] = 'none'
    # Graphical Properties of resizable letters
    path = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        '../components/square.ttf'
    )

    fp = FontProperties(fname=path, weight="bold")
    globscale = 1.22
    letters_shift = -0.5

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
    x = 1
    maxi = 0
    for scores in wdata:
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
