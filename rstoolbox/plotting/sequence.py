from rstoolbox.analysis import sequence_frequency_matrix

import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib.patches import Rectangle

def barcode_plot( df, column_name, axis, color="blue" ):
    a = df[column_name].values
    x = len(a[0])
    result = [0] * x
    for seq in a:
        for _, b in enumerate(seq):
            if bool(int(b)): result[_] = 1
    pd.Series(result).plot("bar", ax=axis, ylim=(0,1), grid=False, color=color, width=1 )
    axis.yaxis.set_ticks([])
    axis.xaxis.set_ticks(np.arange(0, x+1, 10))
    axis.xaxis.set_ticklabels(np.arange(0, x+1, 10) + 1, rotation=45)
    axis.set_xlabel("sequence")

def sequence_frequency_plot( df, column_name, axis, ref_seq=None, key_residues=None, colormap = "Blues", border_color="green", nobar=False, cbar_ax=None ):
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
        sns.heatmap(data, ax=axis, square=True, cbar=not nobar, cbar_kws={"orientation": "horizontal"}, linewidths=1, cmap=colormap)
    else:
        sns.heatmap(data, ax=axis, square=True, cbar_ax=cbar_ax, cbar_kws={"orientation": "horizontal"}, linewidths=1, cmap=colormap)
    order.reverse()
    axis.yaxis.set_ticklabels(order, rotation=0)
    axis.xaxis.set_ticklabels(data.columns.values.tolist())
    axis.set_ylabel("residue type")
    if ref_seq is not None:
        for i in range(len(ref_seq)):
            axis.add_patch(Rectangle((i, order.index(ref_seq[i])), 1, 1, fill=False, edgecolor=border_color, lw=2))
