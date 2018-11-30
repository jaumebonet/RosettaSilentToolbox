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
from rstoolbox.utils import add_top_title, add_column, load_refdata


__all__ = ['multiple_distributions']


def multiple_distributions( df, fig, grid, values="*", titles=None, labels=None,
                            ref=None, seqID=None, ref_range=5, ref_equivalences=None, **kwargs ):
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
    :param str ref: Use reference data against the one in the :class:`.DesignFrame`. Available
        options are ``cath``, ``scop``, ``scop2`` and ``chain``.
    :param str seqID: |seqID_param|. This is necessary if ``ref`` is requested in order to assess
        the length of the query.
    :param int ref_range: When picking reference data, define a protein-length window of
        domain/chain to compare against.
    :param dict ref_equivalences: When names between the query data and the provided data are the
        same, they will be directly assigned. Here a dictionary ``db_name``:``query_name`` can be
        provided otherwise.

    :return: :func:`list` of :class:`~matplotlib.axes.Axes`

    :raises:
        :ValueError: If columns are requested that do not exist in the :class:`~pandas.DataFrame`.
        :ValueError: If the given grid does not have enought positions for all the requested values.
        :ValueError: It the number of values and titles do not match.
        :ValueError: It the number of values and labels do not match.
        :ValueError: It an unknown reference is requested.

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

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: from rstoolbox.plot import multiple_distributions
           ...: import matplotlib.pyplot as plt
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz",
           ...:                         {'sequence': 'A'})
           ...: values = ["score", "hbond_sr_bb", "B_ni_rmsd", "hbond_bb_sc",
           ...:           "cav_vol", "design_score", "packstat", "rmsd_drift"]
           ...: fig = plt.figure(figsize=(25, 10))
           ...: axs = multiple_distributions(df, fig, (2, 4), values, ref='scop2', seqID='A')
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

    if ref is not None:
        if ref.lower() not in ['cath', 'scop', 'scop2', 'chain']:
            raise ValueError('Unknown reference requested.')
        else:
            refdata = load_refdata(ref)
            if seqID not in df.get_available_sequences():
                raise ValueError('column for sequence {} is expected '
                                 'to calculate size'.format(seqID))
            slength = len(df.iloc[0].get_sequence(seqID))
            refdata = refdata[(refdata['length'] >= slength - ref_range) &
                              (refdata['length'] <= slength + ref_range)]
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
            sns.violinplot(x='violinx', y=values[_], hue='target', data=qd, ax=ax,
                           hue_order=["query", "reference"], split=True)
            ax.get_legend().remove()
            ax.set_xlabel('')
            ax.set_xticklabels('')
        if titles is not None:
            add_top_title(ax, titles[_])
        if labels is not None:
            ax.set_ylabel(labels[_])
        axis.append(ax)

    return axis
