# @Author: Jaume Bonet <bonet>
# @Date:   20-Feb-2018
# @Email:  jaume.bonet@gmail.com
# @Filename: structure.py
# @Last modified by:   bonet
# @Last modified time: 20-Feb-2018


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


def positional_structural_similarity_plot( df, ax ):
    """
    Generates a bar plot for positional prevalence of secondary structure elements.
    Input data can/should be generated with :py:func:`.positional_structural_similarity`.
    If there is a `identity_perc` column present, it will also print a line showing how
    often the secondary structure matches the expected/reference one. One can obtain
    this column with the :py:func:`.positional_structural_identity` function and mixing
    both :py:class:`~pandas.DataFrame` with::

        pd.concat([df1, df2], axis=1)

    :param df: Input data, where rows are positions and columns are `H`, `E and `L`
    :type df: :py:class:`~pandas.DataFrame`
    :param ax: matplotlib axis to which we will plot.
    :type ax: :py:class:`~matplotlib.axes.Axes`

    """

    # @todo Expand controls for positional_structural_similarity_plot
    # @body add attributes to change positive and identity colors

    h = df["H"].values
    e = df["E"].values

    ax.bar(range(len(h)), h, color="royalblue" )
    ax.bar(range(len(h)), e, bottom=h, color="tomato")
    if "identity_perc" in df:
        ax.plot(range(len(h)), df["identity_perc"].values, linestyle="solid", linewidth="3", color="black")

    ax.set_ylim(0, 1)
    ax.set_xlim(-0.5, len(h))
