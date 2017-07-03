# @Author: Jaume Bonet <bonet>
# @Date:   29-Jun-2017
# @Email:  jaume.bonet@gmail.com
# @Filename: logo.py
# @Last modified by:   bonet
# @Last modified time: 30-Jun-2017

from __future__ import print_function
import sys

import matplotlib as mpl
mpl.use('tkagg')
import seaborn
import matplotlib.pyplot as plt
plt.style.use('seaborn-ticks')
from matplotlib import transforms
import matplotlib.patheffects
from matplotlib.font_manager import FontProperties
import numpy as np

class Scale(matplotlib.patheffects.RendererBase):
    def __init__(self, sx, sy=None):
        self._sx = sx
        self._sy = sy

    def draw_path(self, renderer, gc, tpath, affine, rgbFace):
        affine = affine.identity().scale(self._sx, self._sy)+affine
        renderer.draw_path(gc, tpath, affine, rgbFace)

COLOR_SCHEME = {
    'A': '#CCFF00', 'C': '#FFFF00', 'D': '#FF0000', 'E': '#FF0066',
    'F': '#00FF66', 'G': '#FF9900', 'H': '#0066FF', 'I': '#66FF00',
    'K': '#6600FF', 'L': '#33FF00', 'M': '#00FF00', 'N': '#CC00FF',
    'P': '#FFCC00', 'Q': '#FF00CC', 'R': '#0000FF', 'S': '#FF3300',
    'T': '#FF6600', 'V': '#99FF00', 'W': '#00CCFF', 'Y': '#00FFCC'
}

def dataframe2logo( data ):
    aa = list(data)
    odata = []
    for index, pos in data.iterrows():
        pdata = []
        for k in aa:
            if pos[k] > 0.0:
                pdata.append( ( k, float(format(pos[k] * 2.23, '.2f')) ) )
        odata.append(sorted(pdata, key=lambda x: x[1]))
    return odata

def plot_logo( data, filename, ref_seq = None ):
    if ref_seq is not None:
        assert len(ref_seq) == len(data)

    fontfamily='Arial'
    size=80
    mpl.rcParams['font.family'] = fontfamily
    fig, ax = plt.subplots(figsize=(len(data), 2.3))
    font = FontProperties()
    font.set_size(size)
    font.set_weight('bold')

    ax.set_xticks(range(1, len(data) + 1))
    ax.set_yticks( range(0, 2) )
    ax.set_xticklabels( data.index.values )
    ax.set_yticklabels( np.arange( 0, 2 , 1 ) )
    seaborn.despine(ax=ax, trim=True)

    trans_offset = transforms.offset_copy(ax.transData,
                                          fig=fig,
                                          x=1,
                                          y=0,
                                          units='dots')
    wdata = dataframe2logo( data )
    sys.stdout.write("plotting: ")
    for index, scores in enumerate( wdata ):
        yshift = 0
        if index % 10 != 0: sys.stdout.write('.')
        else: sys.stdout.write('*')
        sys.stdout.flush()
        for aa, score in scores:
            txt = ax.text( index + 1, 0, aa,
                          transform=trans_offset,
                          fontsize=80,
                          color=COLOR_SCHEME[aa],
                          ha='center',
                          fontproperties=font,
                         )
            txt.set_path_effects([Scale(1.0, score)])
            fig.canvas.draw()
            window_ext = txt.get_window_extent(txt._renderer)
            yshift = window_ext.height * score / 1.78
            trans_offset = transforms.offset_copy(txt._transform,
                                                  fig=fig,
                                                  y=yshift,
                                                  units='points')
        trans_offset = transforms.offset_copy(ax.transData,
                                                  fig=fig,
                                                  x=1,
                                                  y=0,
                                                  units='points')

    fig.savefig( filename )
