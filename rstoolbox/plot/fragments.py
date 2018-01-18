# @Author: Jaume Bonet <bonet>
# @Date:   17-Jan-2018
# @Email:  jaume.bonet@gmail.com
# @Filename: fragments.py
# @Last modified by:   bonet
# @Last modified time: 18-Jan-2018


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

import rstoolbox.utils as ru

def plot_fragments(small_frags, large_frags, small_ax, large_ax, small_color=0, large_color=0,
                   small_max=None, large_max=None, titles="top", **kwargs):
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

    # Data compactness
    small_frags = small_frags[small_frags["position"] == small_frags["frame"]]
    large_frags = large_frags[large_frags["position"] == large_frags["frame"]]

    sns.boxplot(x="frame", y="rmsd", data=small_frags, ax=small_ax, color=small_color, **kwargs)
    sns.boxplot(x="frame", y="rmsd", data=large_frags, ax=large_ax, color=large_color, **kwargs)

    # Basic formating
    small_ax.set_xticks(range(0, max(small_frags["frame"]), 5))
    small_ax.set_xticklabels(range(1, max(small_frags["frame"]) + 1, 5))
    small_ax.set_xlabel("sequence")
    small_ax.set_ylabel("RMSD")
    if small_max is not None:
        small_ax.set_ylim(0, small_max)
    else:
        small_ax.set_ylim(ymin=0)
    small_ax.grid(which='both')

    large_ax.set_xticks(range(0, max(large_frags["frame"]), 5))
    large_ax.set_xticklabels(range(1, max(large_frags["frame"]) + 1, 5))
    large_ax.set_xlabel("sequence")
    large_ax.set_ylabel("RMSD")
    if large_max is not None:
        large_ax.set_ylim(0, large_max)
    else:
        large_ax.set_ylim(ymin=0)
    large_ax.grid(which='both')

    # Titles
    if titles.lower() == "top":
        ru.add_top_title(small_ax, "{}mers".format(small_frags["size"].values[0]))
        ru.add_top_title(large_ax, "{}mers".format(large_frags["size"].values[0]))
    elif titles.lower() == "right":
        ru.add_right_title(small_ax, "{}mers".format(small_frags["size"].values[0]))
        ru.add_right_title(large_ax, "{}mers".format(large_frags["size"].values[0]))
    else:
        pass
