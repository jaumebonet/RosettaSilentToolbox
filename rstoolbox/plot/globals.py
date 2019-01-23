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
from rstoolbox.utils import add_top_title, add_column, discrete_cmap_from_colors


__all__ = ['multiple_distributions', 'plot_in_context', 'distribution_quality']


def multiple_distributions( df, fig, grid, igrid=None, values="*", titles=None, labels=None,
                            refdata=None, ref_equivalences=None, violins=True, legends=False, **kwargs ):
    """Automatically plot boxplot distributions for multiple score types of the
    decoy population.

    A part from the fixed options, the function accepst any option of
    :func:`~seaborn.boxplot`, except for ``y``, ``data`` and ``ax``, which
    are used internally by this function.

    :param df: Data container.
    :type df: :class:`~pandas.DataFrame`
    :param fig: Figure into which the data is going to be plotted.
    :type fig: :class:`~matplotlib.figure.Figure`
    :param grid: Shape of the grid to plot the values in the figure (rows x columns).
    :type grid: :class:`tuple` with two :class:`int`
    :param igrid: Initial position of the grid. Defaults to (0, 0)
    :type igrid: :class:`tuple` with two :class:`int`
    :param values: Contents from the data container that are expected to be plotted.
    :type values: :func:`list` of :class:`str`
    :param titles: Titles to assign to the value of each plot (if provided).
    :type titles: :func:`list` of :class:`str`
    :param labels: Y labels to assign to the value of each plot. By default this will be
        the name of the value.
    :type labels: :func:`list` of :class:`str`
    :param refdata: Data content to use as reference.
    :type refdata: :class:`~pandas.DataFrame`
    :param dict ref_equivalences: When names between the query data and the provided data are the
        same, they will be directly assigned. Here a dictionary ``db_name``:``query_name`` can be
        provided otherwise.
    :param bool violins: When :data:`True`, plot refdata comparisson with violins, otherwise do it
        with kdplots.
    :param bool legends: When :data:`True`, show the legends of each axis.

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
           ...: axs = multiple_distributions(df, fig, (2, 4), values=values)
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
           ...: axs = multiple_distributions(df, fig, (2, 4), values=values, refdata=refdf)
           ...: plt.tight_layout()

        @savefig multiple_distributions_docs2.png width=5in
        In [2]: plt.show()

        In [3]: plt.close()
    """
    if igrid is None:
        igrid = (0, 0)
    if values == "*":
        values = df.select_dtypes(include=[np.number]).columns.tolist()
    if len(set(values).difference(set(list(df.columns)))) > 0:
        raise ValueError("Some of the requested values do not exist "
                         "in the data container.")
    if (grid[0] * grid[1]) - (igrid[0] * igrid[1]) < len(values):
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

    axis_args = {'rowspan': kwargs.pop('rowspan', 1), 'colspan': kwargs.pop('colspan', 1)}

    axis = []
    for _, pgrid in enumerate(itertools.product(*[range(grid[0]), range(grid[1])])):
        if _ >= len(values):
            break
        pgrid = list(pgrid)
        pgrid[0] += igrid[0]
        pgrid[1] += igrid[1]
        ax = plt.subplot2grid(grid, pgrid, fig=fig, rowspan=axis_args['rowspan'])
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
                if not legends:
                    ax.get_legend().remove()
                ax.set_xlabel('')
                ax.set_xticklabels('')
            else:
                sns.kdeplot(s1[values[_]], ax=ax, shade=True)
                sns.kdeplot(s2[values[_]], ax=ax, shade=True)
                if not legends:
                    ax.get_legend().remove()
                ax.set_xlabel(values[_])
        if titles is not None:
            add_top_title(ax, titles[_])
        if labels is not None:
            ax.set_ylabel(labels[_])
        axis.append(ax)

    return axis


def plot_in_context( df, fig, grid, refdata, igrid=None, values='*', ref_equivalences=None,
                     legends=False, **kwargs ):
    """Plot position of decoys in a backgroud reference dataset.

    .. note::
        Personalization is possible by providing argument keys for :func:`~seaborn.kdeplot`
        with the prefix ``kde_`` and for :func:`~matplotlib.plot` (for the points) with the
        prefix ``point_``.

    :param df: Data container.
    :type df: :class:`~pandas.DataFrame`
    :param fig: Figure into which the data is going to be plotted.
    :type fig: :class:`~matplotlib.figure.Figure`
    :param grid: Shape of the grid to plot the values in the figure (rows x columns).
    :type grid: :class:`tuple` with two :class:`int`
    :param refdata: Data content to use as reference.
    :type refdata: :class:`~pandas.DataFrame`
    :param igrid: Initial position of the grid. Defaults to (0, 0)
    :type igrid: :class:`tuple` with two :class:`int`
    :param values: Contents from the data container that are expected to be plotted.
    :type values: :func:`list` of :class:`str`
    :param dict ref_equivalences: When names between the query data and the provided data are the
        same, they will be directly assigned. Here a dictionary ``db_name``:``query_name`` can be
        provided otherwise.
    :param bool legends: When :data:`True`, show the legends of each axis.

    :return: :func:`list` of :class:`~matplotlib.axes.Axes`

    :raises:
        :ValueError: If columns are requested that do not exist in the :class:`~pandas.DataFrame` of
            data **and** reference.
        :ValueError: If the given grid does not have enought positions for all the requested values.
        :ValueError: If ``refdata`` or ``df`` are not :class:`~pandas.DataFrame`.

    .. rubric:: Example:

    .. ipython::
        :okwarning:

        In [1]: from rstoolbox.plot import plot_in_context
           ...: from rstoolbox.utils import load_refdata
           ...: import matplotlib.pyplot as plt
           ...: df = load_refdata('scop2')
           ...: slength = 100
           ...: refdf = load_refdata('scop2', 50)
           ...: refdf = refdf[(refdf['length'] >= slength - 5) &
           ...:               (refdf['length'] <= slength + 5)]
           ...: values = ["score", "hbond_sr_bb", "avdegree", "hbond_bb_sc",
           ...:           "cavity", "CYDentropy", "pack", "radius"]
           ...: fig = plt.figure(figsize=(25, 10))
           ...: axs = plot_in_context(df, fig, (2, 4), refdata=refdf, values=values)
           ...: plt.tight_layout()

        @savefig plot_in_context_docs1.png width=5in
        In [2]: plt.show()

        In [3]: plt.close()
    """
    if igrid is None:
        igrid = (0, 0)
    if not isinstance(df, pd.DataFrame):
        raise ValueError('Unknown data format.')
    if not isinstance(refdata, pd.DataFrame):
        raise ValueError('Unknown reference data format.')
    if values == "*":
        values = df.select_dtypes(include=[np.number]).columns.tolist()
    if ref_equivalences is not None:
        refdata = refdata.rename(columns=ref_equivalences)
    if len(set(values).difference(set(list(df.columns)))) > 0:
        raise ValueError("Some of the requested values do not exist "
                         "in the data container.")
    if len(set(values).difference(set(list(refdata.columns)))) > 0:
        raise ValueError("Some of the requested values do not exist "
                         "in the reference data container.")
    if (grid[0] * grid[1]) - (igrid[0] * igrid[1]) < len(values):
        raise ValueError("The grid does not provide enought positions for all"
                         " requested values.")

    pk = [k for k in kwargs if k.startswith('point_')]
    kk = [k for k in kwargs if k.startswith('kde_')]
    kwargs_point = {}
    kwargs_kde = {}

    for p in pk:
        kwargs_point.setdefault(p.replace('point_', ''), kwargs[p])
    kwargs_point.setdefault('marker', 'o')
    kwargs_point.setdefault('color', 'orange')
    for k in kk:
        kwargs_kde.setdefault(k.replace('kde_', ''), kwargs[k])
    kwargs_kde.setdefault('shade', True)

    axis = []
    for _, pgrid in enumerate(itertools.product(*[range(grid[0]), range(grid[1])])):
        pgrid = list(pgrid)
        pgrid[0] += igrid[0]
        pgrid[1] += igrid[1]
        if _ >= len(values):
            break
        ax = plt.subplot2grid(grid, pgrid, fig=fig)
        kde = sns.kdeplot(refdata[values[_]], ax=ax, **kwargs_kde)
        data_x, data_y = kde.lines[0].get_data()
        for __, row in df.iterrows():
            ref_x = row[values[_]]
            ref_y = np.interp(ref_x, data_x, data_y)
            kde.plot([ref_x], [ref_y], **kwargs_point)
        if not legends:
            ax.get_legend().remove()
        ax.set_xlabel(values[_])
        axis.append(ax)

    return axis


def distribution_quality( df, refdata, values, ascending, names, fig):
    """Locate the quantile position of each putative :class:`.DesingSerie`
    in a list of score distributions.

    :param df: Data container.
    :type df: :class:`~pandas.DataFrame`

    :param grid: Shape of the grid to plot the values in the figure (rows x columns).
    :type grid: :class:`tuple` with two :class:`int`
    :param refdata: Data content to use as reference.
    :type refdata: :class:`~pandas.DataFrame`
    :param values: Contents from the data container that are expected to be plotted.
    :type values: :func:`list` of :class:`str`
    :param ascending: Way the data should be sorted. :data:`True` if the score is better
        when lower, :data:`False` otherwise.
    :type ascending: :func:`list` of :class:`bool`
    :param names: Columns to use as identifiers for the query data.
    :type names: :func:`list` of :class:`str`
    :param fig: Figure into which the data is going to be plotted.
    :type fig: :class:`~matplotlib.figure.Figure`

    :return: :class:`~matplotlib.axes.Axes`

    :raises:
        :ValueError: If columns are requested that do not exist in the :class:`~pandas.DataFrame` of
            data **and** reference.
        :ValueError: If there isn't a ``ascending`` definition for each ``value``.
        :ValueError: If ``refdata`` or ``df`` are not :class:`~pandas.DataFrame`.
        :valueError: If the requested names do not exist in the input data.

    .. rubric:: Example:

    .. ipython::
        :okwarning:

        In [1]: from rstoolbox.plot import distribution_quality
           ...: from rstoolbox.utils import load_refdata
           ...: import matplotlib.pyplot as plt
           ...: df = load_refdata('scop')
           ...: qr = pd.DataFrame([['2F4V', 'C'], ['3BFU', 'B'], ['2APJ', 'C'],
           ...:                    ['2C37', 'V'], ['2I6E', 'H']],
           ...:                   columns=['pdb', 'chain'])
           ...: qr = qr.merge(df, on=['pdb', 'chain'])
           ...: refs = []
           ...: for i, t in qr.iterrows():
           ...:     refs.append(df[(df['length'] >= (t['length'] - 5)) &
           ...:                    (df['length'] <= (t['length'] + 5))])
           ...: fig = plt.figure(figsize=(25, 6))
           ...: ax = distribution_quality(df=qr, refdata=refs,
           ...:                           values=['score', 'pack', 'avdegree',
           ...:                                   'cavity', 'psipred'],
           ...:                           ascending=[True, False, True, True, False],
           ...:                           names=['pdb', 'chain'], fig=fig)
           ...: plt.tight_layout()

        @savefig distribution_quality_docs1.png width=5in
        In [2]: plt.show()

        In [3]: plt.close()
    """
    if not isinstance(df, pd.DataFrame):
        raise ValueError('Unknown data format.')
    if not isinstance(refdata, (pd.DataFrame, list)):
        raise ValueError('Unknown reference data format.')
    if len(set(values).difference(set(list(df.columns)))) > 0:
        raise ValueError("Some of the requested values do not exist "
                         "in the data container.")
    if len(set(names).difference(set(list(df.columns)))) > 0:
        raise ValueError("Some of the requested identifiers do not exist "
                         "in the data container.")
    if isinstance(refdata, list):
        if len(refdata) != df.shape[0]:
            raise ValueError('If multiple references are provided, '
                             'there should be the same as queries.')
        for i, x in enumerate(refdata):
            if not isinstance(x, pd.DataFrame):
                raise ValueError('Unknown reference {} data format.'.format(i))
            if len(set(values).difference(set(list(x.columns)))) > 0:
                raise ValueError("Some of the requested values do not exist "
                                 "in the {} reference container.".format(i))
    else:
        if len(set(values).difference(set(list(refdata.columns)))) > 0:
            raise ValueError("Some of the requested values do not exist "
                             "in the {} reference container.".format(i))
        refdata = [refdata, ] * len(df.shape[0])

    if len(values) != len(ascending):
        raise ValueError("Number of values and orders should match.")

    ax = plt.subplot2grid((1, 1), (0, 0), fig=fig)
    cmap = discrete_cmap_from_colors([(144.0 / 255, 238.0 / 255, 144.0 / 255),
                                      (135.0 / 255, 206.0 / 255, 250.0 / 255),
                                      (255.0 / 255, 165.0 / 255, 0.0 / 255),
                                      (205.0 / 255, 92.0 / 255, 92.0 / 255)])

    data = []
    labs = []
    identifiers = df[names[0]].map(str)
    for i in range(1, len(names)):
        identifiers += '_' + df[names[i]].map(str)
    df = df.reset_index(drop=True)
    for i, row in df.iterrows():
        data.append([])
        labs.append([])
        for isc, sc in enumerate(values):
            qt = refdata[i][sc].quantile([.25, .5, .75])
            if row[sc] <= qt[.25]:
                data[-1].append(.12 if ascending[isc] else .87)
                labs[-1].append('Q1' if ascending[isc] else 'Q4')
            elif row[sc] <= qt[.5]:
                data[-1].append(.37 if ascending[isc] else .67)
                labs[-1].append('Q2' if ascending[isc] else 'Q3')
            elif row[sc] <= qt[.75]:
                data[-1].append(.67 if ascending[isc] else .37)
                labs[-1].append('Q3' if ascending[isc] else 'Q2')
            else:
                data[-1].append(.87 if ascending[isc] else .12)
                labs[-1].append('Q4' if ascending[isc] else 'Q1')

    df = pd.DataFrame(data, columns=values, index=identifiers)
    sns.heatmap(df, square=True, cmap=cmap, cbar=False, annot=pd.DataFrame(labs), fmt='s', ax=ax)
    plt.setp( ax.yaxis.get_majorticklabels(), rotation=0 )
    return ax
