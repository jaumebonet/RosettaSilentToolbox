# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: multiple_distributions
"""
# Standard Libraries
import itertools

# External Libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# This Library
from rstoolbox.utils import add_top_title, add_column


__all__ = ['multiple_distributions']


def multiple_distributions( df, fig, grid, values="*", titles=None, labels=None,
                            refdata=None, ref_equivalences=None, violins=True, **kwargs ):
    """Automatically plot boxplot distributions for multiple score types of the
    decoy population.

    A part from the fixed options, the function accepst any option of
    :func:`~seaborn.boxplot`, except for ``y``, ``data`` and ``ax``, which
    are used internally by this function.

    :param df: Data container
    :type df: :class:`~pandas.DataFrame`
    :param fig: Figure into which the data is going to be plotted.
    :type fig: :class:`~matplotlib.figure.Figure`
    :param grid: Shape of the grid to plot the values in the figure (rows x columns).
    :type grid: :class:`tuple` with two :class:`int`
    :param values: Contents from the data container that are expected to be plotted.
    :type values: :func:`list` of :class:`str`
    :param titles: Titles to assign to the value of each plot (if provided).
    :type titles: :func:`list` of :class:`str`
    :param labels: Y labels to assign to the value of each plot. By default this will be
        the name of the value.
    :type labels: :func:`list` of :class:`str`
    :param str refdata: Data content to use as reference.
    :type refdata: :class:`~pandas.DataFrame`
    :param dict ref_equivalences: When names between the query data and the provided data are the
        same, they will be directly assigned. Here a dictionary ``db_name``:``query_name`` can be
        provided otherwise.
    :param bool violins: When :data:`True`, plot refdata comparisson with violins, otherwise do it
        with kdplots.

    :return: :func:`list` of :class:`~matplotlib.axes.Axes`

    :raises:
        :ValueError: If columns are requested that do not exist in the :class:`~pandas.DataFrame`.
        :ValueError: If the given grid does not have enought positions for all the requested values.
        :ValueError: If the number of values and titles do not match.
        :ValueError: If the number of values and labels do not match.
        :ValueError: If ``refdata`` is not :class:`~pandas.DataFrame`.

    .. rubric:: Example 1: Raw design population data.

    .. ipython::
        :okwarning:

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: from rstoolbox.plot import multiple_distributions
           ...: import matplotlib.pyplot as plt
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz")
           ...: values = ["score", "hbond_sr_bb", "B_ni_rmsd", "hbond_bb_sc",
           ...:           "cav_vol", "design_score", "packstat", "rmsd_drift"]
           ...: fig = plt.figure(figsize=(25, 10))
           ...: axs = multiple_distributions(df, fig, (2, 4), values)
           ...: plt.tight_layout()

        @savefig multiple_distributions_docs.png width=5in
        In [2]: plt.show()

        In [3]: plt.close()

    .. rubric:: Example 2: Design population data vs. DB reference.

    .. ipython::
        :okwarning:

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: from rstoolbox.plot import multiple_distributions
           ...: from rstoolbox.utils import load_refdata
           ...: import matplotlib.pyplot as plt
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz",
           ...:                         {'sequence': 'A'})
           ...: slength = len(df.iloc[0]['sequence_A'])
           ...: refdf = load_refdata('scop2')
           ...: refdf = refdf[(refdf['length'] >= slength - 5) &
           ...:               (refdf['length'] <= slength + 5)]
           ...: values = ["score", "hbond_sr_bb", "B_ni_rmsd", "hbond_bb_sc",
           ...:           "cav_vol", "design_score", "packstat", "rmsd_drift"]
           ...: fig = plt.figure(figsize=(25, 10))
           ...: axs = multiple_distributions(df, fig, (2, 4), values, refdata=refdf)
           ...: plt.tight_layout()

        @savefig multiple_distributions_docs2.png width=5in
        In [2]: plt.show()

        In [3]: plt.close()
    """
    if values == "*":
        values = df.select_dtypes(include=[np.number]).columns.tolist()
    if len(set(values).difference(set(list(df.columns)))) > 0:
        raise ValueError("Some of the requested values do not exist "
                         "in the data container.")
    if grid[0] * grid[1] < len(values):
        raise ValueError("The grid does not provide enought positions for all"
                         " requested values.")
    if titles is not None and len(titles) != len(values):
        raise ValueError("Number of expected plots and titles do not match.")
    if labels is not None and len(labels) != len(values):
        raise ValueError("Number of expected labels and titles do not match.")

    if refdata is not None:
        if not isinstance(refdata, pd.DataFrame):
            raise ValueError('Unknown reference data format.')
        if ref_equivalences is not None:
            refdata = refdata.rename(columns=ref_equivalences)
        refvalues = refdata.select_dtypes(include=[np.number]).columns.tolist()
    else:
        refvalues = []

    kwargs.pop("y", None)
    kwargs.pop("data", None)
    kwargs.pop("axis", None)

    axis = []
    for _, x in enumerate(itertools.product(*[range(grid[0]), range(grid[1])])):
        if _ >= len(values):
            break
        ax = plt.subplot2grid(grid, x, fig=fig)
        if values[_] not in refvalues:
            sns.boxplot(y=values[_], data=df, ax=ax, **kwargs)
        else:
            s1 = add_column(pd.DataFrame(df[values[_]]), 'target', 'query')
            s1 = add_column(s1, 'violinx', 1)
            s2 = add_column(pd.DataFrame(refdata[values[_]]), 'target', 'reference')
            s2 = add_column(s2, 'violinx', 1)
            qd = pd.concat([s1, s2])
            if violins:
                sns.violinplot(x='violinx', y=values[_], hue='target', data=qd, ax=ax,
                               hue_order=["query", "reference"], split=True)
                ax.get_legend().remove()
                ax.set_xlabel('')
                ax.set_xticklabels('')
            else:
                sns.kdeplot(s1[values[_]], ax=ax, shade=True)
                sns.kdeplot(s2[values[_]], ax=ax, shade=True)
                ax.get_legend().remove()
                ax.set_xlabel(values[_])
        if titles is not None:
            add_top_title(ax, titles[_])
        if labels is not None:
            ax.set_ylabel(labels[_])
        axis.append(ax)

    return axis
