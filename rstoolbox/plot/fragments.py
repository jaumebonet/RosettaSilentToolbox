# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: plot_fragment_profiles
.. func:: plot_fragments
"""
# Standard Libraries

# External Libraries
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# This Library
import rstoolbox.utils as ru
import rstoolbox.analysis as ra
from .sequence import positional_sequence_similarity_plot
from .structure import positional_structural_similarity_plot


__all__ = ['plot_fragment_profiles', 'plot_fragments']


def plot_fragment_profiles( fig, small_frags, large_frags, ref_seq, ref_sse, matrix="BLOSUM62" ):
    """Plots a full summary of the a :class:`.FragmentFrame` quality with sequence and expected
    secondary structure match.

    :param fig: Figure into which the data is going to be plotted.
    :type fig: :class:`~matplotlib.figure.Figure`
    :param small_frags: Data for the small fragments.
    :type small_frags: :class:`.FragmentFrame`
    :param large_frags: Data for the large fragments.
    :type large_frags: :class:`.FragmentFrame`
    :param str ref_seq: Reference sequence over which to compare.
    :param str ref_sse: Reference secondary structure over which to compare.
    :param str matrix: Sequence similarity matrix to use for calculations.
        Defualt is ``BLOSUM62``.

    :return: :func:`list` of :class:`~matplotlib.axes.Axes`

    .. seealso::
        :func:`.plot_fragments`

    .. rubric:: Example

    .. ipython::
        :okwarning:

        In [1]: from rstoolbox.io import parse_rosetta_fragments
           ...: from rstoolbox.plot import plot_fragment_profiles
           ...: import matplotlib.pyplot as plt
           ...: df3 = parse_rosetta_fragments("../rstoolbox/tests/data/wauto.200.3mers.gz")
           ...: df9 = parse_rosetta_fragments("../rstoolbox/tests/data/wauto.200.9mers.gz")
           ...: df3 = df3.add_quality_measure(None)
           ...: df9 = df9.add_quality_measure(None)
           ...: fig = plt.figure(figsize=(25, 10))
           ...: seq = "ETPYAIALNDRVIGSSMVLPVDLEEFGAGFLFGQGYIKKAEEIREILVCPQGRISVYA"
           ...: sse = "LEEEEEEELLEEEEEEEELLLLHHHHHHHHHHHHLLLLLLLLLLLEEEELLLEEEELL"
           ...: axs = plot_fragment_profiles(fig, df3, df9, seq, sse)
           ...: plt.tight_layout()

        @savefig plot_fragment_profiles_docs.png width=5in
        In [2]: plt.show()

        In [3]: plt.close()
    """

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
    ref_sse.replace("C", "L")
    positional_structural_similarity_plot(
        pd.concat([ra.positional_structural_count(small_frags),
                   ra.positional_structural_identity(small_frags, ref_sse=ref_sse)], axis=1),
        ax10)
    positional_structural_similarity_plot(
        pd.concat([ra.positional_structural_count(large_frags),
                   ra.positional_structural_identity(large_frags, ref_sse=ref_sse)], axis=1),
        ax11)
    positional_sequence_similarity_plot(ra.positional_sequence_similarity(small_frags, "A",
                                                                          ref_seq, matrix=matrix ),
                                        ax00 )
    positional_sequence_similarity_plot(ra.positional_sequence_similarity(large_frags, "A",
                                                                          ref_seq, matrix=matrix ),
                                        ax01 )

    # fix axis
    plt.setp(ax00.get_xticklabels(), visible=False)
    ax00.set_ylabel("aa freq")
    ax00.set_xlim(-0.5, max(small_frags["position"]))
    ax00.set_ylim(0, 1.01)
    plt.setp(ax01.get_xticklabels(), visible=False)
    ax01.set_ylabel("aa freq")
    ax01.set_xlim(-0.5, max(large_frags["position"]))
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

    return [ax00, ax10, ax20, ax01, ax11, ax21]


def plot_fragments( small_frags, large_frags, small_ax, large_ax, small_color=0, large_color=0,
                    small_max=None, large_max=None, titles="top", **kwargs ):
    """
    Plot RMSD quality of a pair of :class:`.FragmentFrame` in two provided axis.
    Thought to more easily print both small and large fragments together.

    On plotting, fragment RMSD values are assigned to the first position of the fragments.
    This means that the plots will have a length of::

        :math:`len(sequence) - len(fragment set)`

    :param small_frags: Data for the small fragments.
    :type small_frags: :class:`.FragmentFrame`
    :param large_frags: Data for the large fragments.
    :type large_frags: :class:`.FragmentFrame`
    :param small_ax: Axis where to print the small fragments.
    :type small_ax: :class:`~matplotlib.axes.Axes`
    :param large_ax: Axis where to print the large fragments.
    :type large_ax: :class:`~matplotlib.axes.Axes`
    :param small_color: Color to use on the small fragments. If string,
        that is the assumed color. If integer, it will provide that position for the
        currently active color palette in seaborn.
    :type small_color: Union[:class:`str`, :class:`int`]
    :param large_color: Color to use on the large fragments. If string,
        that is the assumed color. If integer, it will provide that position for the
        currently active color palette in seaborn.
    :type large_color: Union[:class:`str`, :class:`int`]
    :param float small_max: Max value for the y (RMSD) axis of the small fragments. If
        not provided, the system picks it according to the given data.
    :param float large_max: Max value for the y (RMSD) axis of the large fragments. If
        not provided, the system picks it according to the given data.
    :param str titles: Title placement. Options are "top" or "right". Other options
        will result in no titles added to the plot.

    .. seealso::
        :func:`.plot_fragment_profiles`

    .. rubric:: Example

    .. ipython::
        :okwarning:

        In [1]: from rstoolbox.io import parse_rosetta_fragments
           ...: from rstoolbox.plot import plot_fragments
           ...: import matplotlib.pyplot as plt
           ...: df3 = parse_rosetta_fragments("../rstoolbox/tests/data/wauto.200.3mers.gz")
           ...: df9 = parse_rosetta_fragments("../rstoolbox/tests/data/wauto.200.9mers.gz")
           ...: df3 = df3.add_quality_measure(None)
           ...: df9 = df9.add_quality_measure(None)
           ...: fig = plt.figure(figsize=(35, 10))
           ...: ax00 = plt.subplot2grid((1, 2), (0, 0))
           ...: ax01 = plt.subplot2grid((1, 2), (0, 1))
           ...: plot_fragments(df3, df9, ax00, ax01)
           ...: plt.tight_layout()

        @savefig plot_fragments_docs.png width=5in
        In [2]: plt.show()

        In [3]: plt.close()
    """

    # Color management
    if isinstance(small_color, int):
        small_color = sns.color_palette()[small_color]
    if isinstance(large_color, int):
        large_color = sns.color_palette()[large_color]

    # Data compactness
    small_frags_ = small_frags[small_frags["position"] == small_frags["frame"]]
    large_frags_ = large_frags[large_frags["position"] == large_frags["frame"]]

    sns.boxplot(x="frame", y="rmsd", data=small_frags_, ax=small_ax, color=small_color, **kwargs)
    sns.boxplot(x="frame", y="rmsd", data=large_frags_, ax=large_ax, color=large_color, **kwargs)

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
    large_ax.set_axisbelow(True)

    # Titles
    if titles is not None:
        if titles.lower() == "top":
            ru.add_top_title(small_ax, "{}mers".format(small_frags["size"].values[0]))
            ru.add_top_title(large_ax, "{}mers".format(large_frags["size"].values[0]))
        elif titles.lower() == "right":
            ru.add_left_title(small_ax, "{}mers".format(small_frags["size"].values[0]))
            ru.add_left_title(large_ax, "{}mers".format(large_frags["size"].values[0]))
        else:
            pass
