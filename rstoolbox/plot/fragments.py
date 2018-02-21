# @Author: Jaume Bonet <bonet>
# @Date:   17-Jan-2018
# @Email:  jaume.bonet@gmail.com
# @Filename: fragments.py
# @Last modified by:   bonet
# @Last modified time: 20-Feb-2018


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import rstoolbox.utils as ru
import rstoolbox.analysis as ra
from .sequence import positional_sequence_similarity_plot
from .structure import positional_structural_similarity_plot

def _sse_frequencies_match(sse_fit, frags):

    data_sse = {"position": range(1, len(sse_fit) + 1), "sse": list(sse_fit)}
    df_sse = pd.DataFrame(data_sse)

    ssea = frags[["position", "neighbors", "sse"]].groupby(["position", "sse"]).size().reset_index(name='counts')
    sseb = frags[["position", "neighbors"]].groupby("position").size().reset_index(name='counts')
    sse = pd.merge(ssea, sseb, how="left", on="position")
    sse["percs"] = sse['counts_x']/sse['counts_y']
    sse = sse.drop( columns=["counts_x", "counts_y"] )
    sse = pd.merge( df_sse, sse, how="left", on=["position", "sse"])
    sse["position"] = sse["position"].apply( lambda x: x - 1 )

    return sse.fillna(0)

def _seq_frequencies_match(seq_fit, frags):

    data_sse = {"position": range(1, len(seq_fit) + 1), "aa": list(seq_fit)}
    df_sse = pd.DataFrame(data_sse)

    ssea = frags[["position", "neighbors", "aa"]].groupby(["position", "aa"]).size().reset_index(name='counts')
    sseb = frags[["position", "neighbors"]].groupby("position").size().reset_index(name='counts')
    sse = pd.merge(ssea, sseb, how="left", on="position")
    sse["percs"] = sse['counts_x']/sse['counts_y']
    sse = sse.drop( columns=["counts_x", "counts_y"] )
    sse = pd.merge( df_sse, sse, how="left", on=["position", "aa"])
    sse["position"] = sse["position"].apply( lambda x: x - 1 )

    return sse.fillna(0)

def plot_fragment_profiles( fig, small_frags, large_frags, ref_seq, ref_sse, matrix="BLOSUM62" ):

    # make subplots
    grid = (4, 2)
    ax00 = plt.subplot2grid(grid, (0, 0), fig=fig)
    ax10 = plt.subplot2grid(grid, (1, 0), sharex=ax00, fig=fig)
    ax20 = plt.subplot2grid(grid, (2, 0), rowspan=2, sharex=ax00, fig=fig)

    ax01 = plt.subplot2grid(grid, (0, 1), sharey=ax00, fig=fig)
    ax11 = plt.subplot2grid(grid, (1, 1), sharex=ax01, fig=fig)
    ax21 = plt.subplot2grid(grid, (2, 1), rowspan=2, sharex=ax01, fig=fig)

    # fill subplots
    plot_fragments( small_frags, large_frags, ax20, ax21, titles=None )
    positional_structural_similarity_plot(
        pd.concat([ra.positional_structural_count(small_frags), ra.positional_structural_identity(small_frags, ref_sse=ref_sse)], axis=1),
        ax10)
    positional_structural_similarity_plot(
        pd.concat([ra.positional_structural_count(large_frags), ra.positional_structural_identity(large_frags, ref_sse=ref_sse)], axis=1),
        ax11)
    positional_sequence_similarity_plot( ra.positional_sequence_similarity( small_frags, "A", ref_seq, matrix ), ax00 )
    positional_sequence_similarity_plot( ra.positional_sequence_similarity( large_frags, "A", ref_seq, matrix ), ax01 )

    # fix axis
    plt.setp(ax00.get_xticklabels(), visible=False)
    ax00.set_ylabel("aa freq")
    ax00.set_xlim(xmin=-0.5)
    ax00.set_ylim(0, 1.01)
    plt.setp(ax01.get_xticklabels(), visible=False)
    ax01.set_ylabel("aa freq")
    ax01.set_xlim(xmin=-0.5)
    ax01.set_ylim(0, 1.01)
    plt.setp(ax10.get_xticklabels(), visible=False)
    ax10.set_ylabel("sse freq")
    ax10.set_ylim(0, 1.01)
    plt.setp(ax11.get_xticklabels(), visible=False)
    ax11.set_ylabel("sse freq")
    ax11.set_ylim(0, 1.01)

    # fix display
    ru.add_top_title(ax00, "{}mers".format(small_frags["size"].values[0]))
    ru.add_top_title(ax01, "{}mers".format(large_frags["size"].values[0]))
    fig.subplots_adjust(wspace=0.1, hspace=0.1)

    fig.legend(handles=[
        mpatches.Patch(color="green",     label="aa identity"),
        mpatches.Patch(color="orange",    label="aa similarity"),
        mpatches.Patch(color="royalblue", label="sse helix content"),
        mpatches.Patch(color="tomato",    label="sse beta content"),
        mpatches.Patch(color="black",     label="sse match to expected")
    ], ncol=5, loc='lower center', borderaxespad=0.)


def plot_fragments(small_frags, large_frags, small_ax, large_ax, small_color=0, large_color=0,
                   small_max=None, large_max=None, titles="top", seq_fit=None,
                   small_seq_color=1, large_seq_color=1, sse_fit=None, small_sse_color=2,
                   large_sse_color=2, **kwargs):
    """
    Plot a pair of :py:class:`.FragmentFrame`s in two provided axis. Thought to more easily print
    both small and large fragments.

    :param :py:class:`.FragmentFrame` small_frags: Data for the small fragments.
    :param :py:class:`.FragmentFrame` large_frags: Data for the large fragments.
    :param axis small_ax: Axis where to print the small fragments.
    :param axis large_ax: Axis where to print the large fragments.
    :param small_color: string or int. Color to use on the small fragments. If string,
        that is the assumed color. If integer, it will provide that position for the
        currently active color palette in seaborn.
    :param large_color: string or int. Color to use on the large fragments. If string,
        that is the assumed color. If integer, it will provide that position for the
        currently active color palette in seaborn.
    :param float small_max: Max value for the y (RMSD) axis of the small fragments. If
        not provided, the system picks it according to the given data.
    :param float large_max: Max value for the y (RMSD) axis of the large fragments. If
        not provided, the system picks it according to the given data.
    :param string titles: Title placement. Options are "top" or "right". Other options
        will result in no titles added to the plot.
    """

    # Color management
    if isinstance(small_color, int):
        small_color = sns.color_palette()[small_color]
    if isinstance(large_color, int):
        large_color = sns.color_palette()[large_color]
    if isinstance(small_seq_color, int):
        small_seq_color = sns.color_palette()[small_seq_color]
    if isinstance(large_seq_color, int):
        large_seq_color = sns.color_palette()[large_seq_color]
    if isinstance(small_sse_color, int):
        small_sse_color = sns.color_palette()[small_sse_color]
    if isinstance(large_sse_color, int):
        large_sse_color = sns.color_palette()[large_sse_color]

    # Data compactness
    small_frags_ = small_frags[small_frags["position"] == small_frags["frame"]]
    large_frags_ = large_frags[large_frags["position"] == large_frags["frame"]]

    sns.boxplot(x="frame", y="rmsd", data=small_frags_, ax=small_ax, color=small_color, **kwargs)
    sns.boxplot(x="frame", y="rmsd", data=large_frags_, ax=large_ax, color=large_color, **kwargs)

    # Sequence Variability
    if seq_fit is not None or sse_fit is not None:
        small_twin = small_ax.twinx()
        small_twin.yaxis.grid(False)
        small_twin.xaxis.grid(False)
        large_twin = large_ax.twinx()
        large_twin.yaxis.grid(False)
        large_twin.xaxis.grid(False)

    if seq_fit is not None:
        # small_seq = small_frags[["position", "aa"]].drop_duplicates(["position", "aa"]).groupby("position").count().reset_index()
        # small_seq["aa"] = small_seq["aa"].apply(lambda x: float(x)/21.0 )
        # small_seq["position"] = small_seq["position"].apply( lambda x: x - 1 )
        # small_twin.plot( small_seq["position"], small_seq["aa"], linestyle="solid", linewidth=3, color=small_seq_color)
        # large_seq = large_frags[["position", "aa"]].drop_duplicates(["position", "aa"]).groupby("position").count().reset_index()
        # large_seq["aa"] = large_seq["aa"].apply(lambda x: float(x)/21.0 )
        # large_seq["position"] = large_seq["position"].apply( lambda x: x - 1 )
        # large_twin.plot( large_seq["position"], large_seq["aa"], linestyle="solid", linewidth=3, color=large_seq_color)
        small_seq = _seq_frequencies_match(seq_fit, small_frags)
        large_seq = _seq_frequencies_match(seq_fit,large_frags)
        small_twin.plot( small_seq["position"], small_seq["percs"], linestyle="solid", linewidth=3, color=small_seq_color)
        large_twin.plot( large_seq["position"], large_seq["percs"], linestyle="solid", linewidth=3, color=large_seq_color)

    if sse_fit is not None:
        small_sse = _sse_frequencies_match(sse_fit, small_frags)
        large_sse = _sse_frequencies_match(sse_fit, large_frags)
        small_twin.plot( small_sse["position"], small_sse["percs"], linestyle="dashed", linewidth=3, color=small_sse_color)
        large_twin.plot( large_sse["position"], large_sse["percs"], linestyle="dashed", linewidth=3, color=large_sse_color)

    # Basic formating
    small_ax.set_xticks(range(0, max(small_frags["frame"]), 5))
    small_ax.set_xticklabels(range(1, max(small_frags["frame"]) + 1, 5))
    small_ax.set_xlabel("sequence")
    small_ax.set_ylabel("RMSD")
    if small_max is not None:
        small_ax.set_ylim(0, small_max)
    else:
        small_ax.set_ylim(ymin=0)
    small_ax.yaxis.grid(False)
    small_ax.xaxis.grid(True)
    if seq_fit is not None or sse_fit is not None:
        small_twin.set_ylim(0, 1.01)
        small_twin.set_ylabel("frequency")
    small_ax.set_axisbelow(True)

    large_ax.set_xticks(range(0, max(large_frags["frame"]), 5))
    large_ax.set_xticklabels(range(1, max(large_frags["frame"]) + 1, 5))
    large_ax.set_xlabel("sequence")
    large_ax.set_ylabel("RMSD")
    if large_max is not None:
        large_ax.set_ylim(0, large_max)
    else:
        large_ax.set_ylim(ymin=0)
    large_ax.yaxis.grid(False)
    large_ax.xaxis.grid(True)
    if seq_fit is not None or sse_fit is not None:
        large_twin.set_ylim(0, 1.01)
        large_twin.set_ylabel("frequency")
    large_ax.set_axisbelow(True)


    # Titles
    if titles is not None:
        if titles.lower() == "top":
            ru.add_top_title(small_ax, "{}mers".format(small_frags["size"].values[0]))
            ru.add_top_title(large_ax, "{}mers".format(large_frags["size"].values[0]))
        elif titles.lower() == "right":
            ru.add_right_title(small_ax, "{}mers".format(small_frags["size"].values[0]))
            ru.add_right_title(large_ax, "{}mers".format(large_frags["size"].values[0]))
        else:
            pass
