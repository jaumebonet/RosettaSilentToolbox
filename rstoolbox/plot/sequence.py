# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: plot_sequence_frequency_graph
.. func:: plot_alignment
.. func:: logo_plot
.. func:: positional_sequence_similarity_plot
.. func:: sequence_frequency_plot
.. func:: per_residue_matrix_score_plot
"""
# Standard Libraries
from distutils.version import LooseVersion
import os
import copy
import math
import operator

# External Libraries
import pandas as pd
import numpy as np
import seaborn as sns
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, PathPatch
from matplotlib.font_manager import FontProperties
from matplotlib.text import TextPath

# This Library
from rstoolbox.analysis import binary_overlap, sequence_similarity
from rstoolbox.analysis.SimilarityMatrix import SimilarityMatrix
from rstoolbox.components import DesignFrame, DesignSeries, SequenceFrame
from rstoolbox.components import get_selection
from rstoolbox.utils.getters import _check_column
from rstoolbox.utils import discrete_cmap_from_colors, add_column
from .color_schemes import color_scheme

__all__ = ["plot_sequence_frequency_graph", "plot_alignment", "logo_plot",
           "positional_sequence_similarity_plot", "sequence_frequency_plot",
           "per_residue_matrix_score_plot"]


def barcode_plot( df, column_name, ax, color="blue" ):
    result = binary_overlap( df, column_name )
    pd.Series(result).plot("bar", ax=ax, ylim=(0, 1), grid=False, color=color, width=1 )
    ax.yaxis.set_ticks([])
    ax.xaxis.set_ticks(np.arange(0, len(result) + 1, 10))
    ax.xaxis.set_ticklabels(np.arange(0, len(result) + 1, 10) + 1, rotation=45)
    ax.set_xlabel("sequence")


def plot_sequence_frequency_graph( G, ax ):
    """Given a sequence frequency graph as obtained through
    :meth:`.FragmentFrame.make_frequency_network` or
    :meth:`.FragmentFrame.make_per_position_frequency_network`,
    generate a plot representing the possible transitions between
    nodes.

    :param G: Sequence frequency graph
    :type G: :class:`~networkx.DiGraph`
    :param ax: Where to plot the network.
    :type ax: :class:`~matplotlib.axes.Axes`

    """
    alphabet = "ARNDCQEGHILKMFPSTWYV"
    all_pos = {}
    for node, data in G.nodes(data=True):
        if data['type'] in alphabet:
            all_pos.setdefault(node, (data['order'], alphabet.index(data['type'])))
        else:
            all_pos.setdefault(node, (data['order'], len(alphabet) / 2))
    nx.draw_networkx(G, ax=ax, pos=all_pos, with_labels=False, arrows=False)
    ax.set_yticks(range(0, len(alphabet)))
    ax.set_yticklabels(list(alphabet))
    ax.set_xlim(G.node["0X"]['order'] - 1, G.node["-1X"]['order'] + 1)


def sequence_frequency_plot( df, seqID, ax, aminosY=True, clean_unused=-1,
                             refseq=True, key_residues=None, border_color="green",
                             border_width=2, labelsize=None, xrotation=0, yrotation=0,
                             **kwargs ):
    """Makes a heatmap subplot into the provided axis showing the sequence distribution
    of each residue type for each position.

    A part from the function arguments, any argument that can be provided to the
    :func:`seaborn.heatmap` function can also be provided here.

    By default, the heatmap generated will have the residue types as y-axis and the
    sequence positions as x-axis.

    Some tips:

    #. **Do you want to set the orientation of the color bar vertical?** \
        Add the parameter: ``cbar_kws={"orientation": "vertical"}``
    #. **Do you want to put the color bar in a different axis?** \
        This is quite recommendable, as the color bar in the same axis does not \
        tend to look that good. Add the parameter: ``cbar_ax=[second_axis]``
    #. **You don't want a color bar?** \
        Add the parameter: ``cbar=False``

    :param df: Data container.
    :type df: Union[:class:`.DesignFrame`, :class:`.SequenceFrame`]
    :param str seqID: |seqID_param|.
    :param ax: Where to plot the heatmap.
    :type ax: :class:`~matplotlib.axes.Axes`
    :param bool aminosY: Set to :data:`False` to get invert the orientation of the heatmap.
    :param float clean_unused: Remove amino acids from the plot when they never get represented
        over the given frequency. Residues present in the reference sequence are not taken
        into account.
    :param rbool refseq: if :data:`True` (default), mark the original residues according to
        the reference sequence.
    :param key_residues: |keyres_param|.
    :type key_residue: |keyres_types|
    :param border_color: Color to use to mark the original residue types.
    :type border_color: Union[:class:`int`, :class:`str`]
    :param int border_width: Line width used to mark the original residue types.
    :param int labelsize: Change the size of the text in the axis.
    :param float xrotation: Rotation to apply in the x-axis text (degrees).
    :param float yrotation: Rotation to apply in the y-axis text (degrees).

    .. note::

        Attribute ``clean_unused``, if applied deletes the full column/row assigned to an
        unrepresented residue type, this means that if that residue type is part of the
        ``refseq``, it will not be labeled.

    :raises:
        :ValueError: if input is not a :class:`~pandas.DataFrame` derived object.
        :KeyError: |reference_error|.

    .. rubric:: Example

    .. ipython::
        :okwarning:

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: from rstoolbox.plot import sequence_frequency_plot
           ...: import matplotlib.pyplot as plt
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz",
           ...:                         {"sequence": "B"})
           ...: fig = plt.figure(figsize=(25, 10))
           ...: ax = plt.subplot2grid((1, 1), (0, 0))
           ...: sequence_frequency_plot(df, "B", ax, refseq=False, cbar=False, xrotation=90)

        @savefig sequence_frequency_plot_docs.png width=5in
        In [2]: plt.show()

        In [3]: plt.close()
    """

    order = ["A", "V", "I", "L", "M", "F", "Y", "W", "S", "T", "N",
             "Q", "R", "H", "K", "D", "E", "C", "G", "P"]
    data = copy.deepcopy(df)

    fp = FontProperties()
    fp.set_family("monospace")

    # Data type management.
    if not isinstance(data, pd.DataFrame):
        raise ValueError("Input data must be in a DataFrame, DesignFrame or SequenceFrame")
    else:
        if not isinstance(data, (DesignFrame, SequenceFrame)):
            if len(set(data.columns.values).intersction(set(order))) == len(order):
                data = SequenceFrame(data)
            else:
                data = DesignFrame(data)
    if isinstance(data, DesignFrame):
        data = data.sequence_frequencies(seqID)
    if isinstance(data, SequenceFrame):
        order = sorted(data.columns.values.tolist(), key=lambda x: order.index(x))
        if not data.is_transposed():
            data = data.transpose().reindex(order)
        else:
            data = data.reindex(order)

    # Refseq and key_residues management.
    ref_seq = data.get_reference_sequence(seqID, key_residues) if refseq else ""

    # data and key_residues management.
    data = data.get_key_residues(key_residues)

    if clean_unused >= 0:
        data.delete_empty(clean_unused)
        data = data.clean()
        order = sorted(data.index.values.tolist(), key=lambda x: order.index(x))
        data = data.reindex(order)

    # heatmap parameters and others
    kwargs.setdefault("cmap", "Blues")  # define the color-range of the plot
    kwargs.setdefault("linewidths", 1)  # linewidths are fixed to 1
    kwargs.setdefault("square", True)   # square is True if user don't say otherwise
    # by default the color bar is horizontal
    kwargs.setdefault("cbar_kws", {"orientation": "horizontal"})

    # plot
    if not aminosY:
        data = data.transpose()
    sns.heatmap(data, ax=ax, **kwargs)

    # styling plot
    # seaborn made a change in the ticks from 0.7 to 0.8,
    # this should take care that both versions work ok.
    if LooseVersion(sns.__version__) < LooseVersion("0.8"):
        order.reverse()
    if aminosY:
        ax.yaxis.set_ticks(np.arange(0.5, len(order) + 0.5))
        ax.yaxis.set_ticklabels(order, rotation=yrotation)
        for label in ax.get_yticklabels():
            label.set_fontproperties(fp)
        ax.xaxis.set_ticks(np.arange(0.5, len(data.columns.values.tolist()) + 0.5))
        ax.xaxis.set_ticklabels(data.columns.values.tolist(), rotation=xrotation)
        ax.set_ylabel("residue type")
        if labelsize is not None:
            ax.tick_params(labelsize=labelsize)
    else:
        ax.xaxis.set_ticks(np.arange(0.5, len(order) + 0.5))
        ax.xaxis.set_ticklabels(order, rotation=xrotation)
        for label in ax.get_xticklabels():
            label.set_fontproperties(fp)
        ax.yaxis.set_ticks(np.arange(0.5, len(data.index.values.tolist()) + 0.5))
        ax.yaxis.set_ticklabels(data.index.values.tolist(), rotation=yrotation)
        ax.set_xlabel("residue type")
        if labelsize is not None:
            ax.tick_params(labelsize=labelsize)

    # marking reference sequence
    if ref_seq is not "" and refseq:
        if isinstance(border_color, int):
            border_color = sns.color_palette()[border_color]
        for i, aa in enumerate(ref_seq):
            if aminosY:
                if aa in order:
                    aa_position = (i, order.index(aa))
            else:
                if aa in order:
                    aa_position = (order.index(aa), i)
            ax.add_patch(Rectangle(aa_position, 1, 1, fill=False, clip_on=False,
                                   edgecolor=border_color, lw=border_width, zorder=100))


def plot_alignment( df, seqID, ax, line_break=None, matrix=None ):
    """Make an image representing the alignment of sequences with higlights to mutant positions.

    :param df: Data container.
    :type df: Union[:class:`.DesignFrame`, :class:`.SequenceFrame`]
    :param str seqID: |seqID_param|.
    :param ax: Where to plot the heatmap. If a list of axis is provided, it assumes that one wants
        to split the alignment in that many pieces.
    :type ax: Union[:class:`~matplotlib.axes.Axes`, :func:`list` of
        :class:`~matplotlib.axes.Axes`]
    :param int line_break: Set a line break after x residues. To execute that,
        multiple axes must be provided.
    :param str matrix: Identifier of the matrix used to evaluate similarity. Default is
        :data:`None` highlight differences.

    :raises:
        :ValueError: if input is not a :class:`~pandas.DataFrame` derived object.
        :KeyError: if reference sequence is requested but the data container
            does not have one.
        :KeyError: |reference_error|.
    """
    def score(col, matrix):
        if matrix is None:
            return pd.Series([0 if _ == col.name else 1 for _ in col.values])
        else:
            data = []
            for _ in col.values:
                if _ == col.name:
                    data.append(0)
                else:
                    data.append(1 if int(matrix.get_value(col.name, _)) >= 0 else -1)
            return pd.Series(data)

    def chunks(l, n):
        """Yield successive n-sized chunks from l."""
        # https://stackoverflow.com/a/312464/2806632
        for i in range(0, len(l), n):
            yield l[i:i + n]

    def split(a, n):
        k, m = divmod(len(a), n)
        return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n))

    if not isinstance(df, pd.DataFrame):
        raise ValueError("Input data must be in a DataFrame or DesignFrame")
    if not df.has_reference_sequence(seqID):
        raise KeyError("Reference sequence is needed")

    seq = list(df.get_reference_sequence(seqID))
    pos = df.get_reference_shift(seqID)
    ids = list(df.get_id().values)
    if isinstance(pos, int):
        pos = list(range(pos, len(seq) + pos))

    cname = _check_column(df, "sequence", seqID)
    df = df[cname].apply(lambda x: pd.Series(list(x)))
    df.columns = seq
    if matrix is not None:
        matrix = SimilarityMatrix.get_matrix(matrix)
    ref = pd.DataFrame([seq])
    ref.columns = seq
    df = pd.concat([ref, df], ignore_index=True)
    ids.insert(0, "reference")
    df2 = df.apply(lambda col: score(col, matrix), axis=0)
    df2.index = ids

    if matrix is None:
        cmap = discrete_cmap_from_colors([(255.0 / 255, 255.0 / 255, 255.0 / 255),
                                          (255.0 / 255, 255.0 / 255, 0.0 / 255)])
    else:
        cmap = discrete_cmap_from_colors([(255.0 / 255, 0.0 / 255, 0.0 / 255),
                                          (255.0 / 255, 255.0 / 255, 255.0 / 255),
                                          (0.0 / 255, 153.0 / 255, 0.0 / 255)])

    if not isinstance(ax, list):
        sns.heatmap(df2, ax=ax, annot=df, cbar=False, square=True, fmt='', cmap=cmap)
    else:
        if line_break is None:
            chnk = list(split(range(df2.shape[1]), len(ax)))
        else:
            chnk = list(chunks(range(df2.shape[1]), line_break))
            if len(chnk) > len(ax):
                raise ArithmeticError("Not enough axis for the number of requested pieces")
        max_len = 0
        for i, axis in enumerate(ax):
            counts = df2.iloc[:, chnk[i][0]:chnk[i][-1]]
            names  = df.iloc[:, chnk[i][0]:chnk[i][-1]]
            if i == 0:
                max_len = counts.shape[1]
            if i == len(ax) - 1:
                this_len = counts.shape[1]
                for _ in range(max_len - this_len):
                    counts = add_column(counts, "", 0)
                    names = add_column(names, "", "")
            sns.heatmap(counts, ax=axis, cbar=False, fmt='', annot=names, square=True, cmap=cmap)
            # @TODO: Add numerical count on ticks
            # @BODY: pick the numbers from the reference_shift and add them as ticks every 10 or so
            if i == len(ax) - 1:
                axis.set_xticks([_ + 0.5 for _ in list(range(this_len))])
            axis.set_xticklabels([])
            axis.set_xticks([])


def positional_sequence_similarity_plot( df, ax, identity_color="green",
                                         similarity_color="orange" ):
    """Generates a plot covering the amount of identities and positives matches from a population
    of designs to a reference sequence according to a substitution matrix.

    Input data can/should be generated with :func:`.positional_sequence_similarity`.

    :param df: Input data, where rows are positions and columns are
        `identity_perc` and `positive_perc`
    :type df: :py:class:`~pandas.DataFrame`
    :param ax: matplotlib axis to which we will plot.
    :type ax: :py:class:`~matplotlib.axes.Axes`
    """

    # Color management
    if isinstance(identity_color, int):
        identity_color = sns.color_palette()[identity_color]
    if isinstance(similarity_color, int):
        similarity_color = sns.color_palette()[similarity_color]

    y = df["positive_perc"].values
    ax.plot(range(len(y)), y, color="orange", linestyle="solid", linewidth=2)
    ax.fill_between(range(len(y)), 0, y, color=similarity_color, alpha=1)

    y = df["identity_perc"].values
    ax.plot(range(len(y)), y, color="green", linestyle="solid", linewidth=2)
    ax.fill_between(range(len(y)), 0, y, color=identity_color, alpha=1)

    ax.set_ylim(0, 1)
    ax.set_xlim(0, len(y) - 1)


def per_residue_matrix_score_plot( df, seqID, ax, matrix="BLOSUM62",
                                   selections=None, **kwargs ):
    """Plot a linear representation of the scoring obtained by applying a
    substitution matrix.

    Applies to a single decoy against the ``reference_sequence``.

    Parameters to control the properties of the plotted line (``color``,
    ``linestyle``...) can be provided too.

    :param df: |df_param|
    :type df: :class:`.DesignSeries`
    :param str seqID: |seqID_param|
    :param ax: matplotlib axis to which we will plot.
    :type ax: :py:class:`~matplotlib.axes.Axes`
    :param str matrix: |matrix_param|
    :param selections: List of regions to highlight; each position should be
        a selector and a color.
    :type selections: :func:`list` of :class:`tuple` with |keyres_types| and
        a color (:class:`str` or :class:`int`)

    :raises:
        :ValueError: If the data container is not :class:`.DesignSeries` or it
            does not have a ``reference_sequence``.

    .. rubric:: Example

    .. ipython::
        :okwarning:

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: from rstoolbox.plot import per_residue_matrix_score_plot
           ...: import matplotlib.pyplot as plt
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz",
           ...:                         {"sequence": "B"})
           ...: df.add_reference_sequence('B', df.iloc[0]['sequence_B'])
           ...: df.add_reference_shift('B', 10)
           ...: seles = [('15-25', 'red'), ('45B-60B', 'green')]
           ...: fig = plt.figure(figsize=(25, 10))
           ...: ax0 = plt.subplot2grid((2, 1), (0, 0))
           ...: per_residue_matrix_score_plot(df.iloc[1], "B", ax0)
           ...: ax1 = plt.subplot2grid((2, 1), (1, 0))
           ...: per_residue_matrix_score_plot(df.iloc[1], "B", ax1, selections=seles)

        @savefig per_residue_matrix_score_plot_docs.png width=5in
        In [2]: plt.show()

        In [3]: plt.close()
    """
    if not isinstance(df, DesignSeries) or not df.has_reference_sequence(seqID):
        raise ValueError("Data must be a DesignSeries with reference for the requested seqID")

    shift = df.get_reference_shift(seqID)
    refsq = df.get_reference_sequence(seqID)

    column = '{0}_{1}_per_res'.format(matrix.lower(), seqID)
    if column not in df.index:
        df = sequence_similarity(df.to_frame().T, seqID, matrix=matrix).iloc[0]

    ax.plot(range(0, len(refsq)), [0, ] * len(refsq), color='grey', linestyle='dashed')
    ax.plot(range(0, len(refsq)), df[column], **kwargs)

    ax.set_xlim(0, len(refsq) - 1)
    ax.set_xticks(range(0, len(refsq), 5))
    if isinstance(shift, int):
        ax.set_xticklabels([_ + shift for _ in range(0, len(refsq), 5)])
    else:
        ax.set_xticklabels(shift[0::5])

    axb = ax.twiny()
    axb.set_xticks(range(0, len(refsq)))
    axb.set_xticklabels(list(df['{0}_{1}_ali'.format(matrix.lower(), seqID)].replace('.', ' ')))
    axb.tick_params('x', top=False, pad=0)

    axlim = ax.get_ylim()
    if selections is None:
        selections = []
    for s in selections:
        s_ = get_selection(s[0], seqID, shift, len(refsq))
        # plot starts in 0
        ax.fill([s_[0] - 1, s_[-1] - 1, s_[-1] - 1, s_[0] - 1],
                [axlim[0] - 1, axlim[0] - 1, axlim[1] + 1, axlim[1] + 1],
                color=s[1], alpha=0.2, zorder=-100)
    ax.set_ylim(axlim[0], axlim[1])

    ax.set_ylabel(matrix.upper())


def logo_plot( df, seqID, refseq=True, key_residues=None, line_break=None,
               hight_prop=4, font_size=35, refplot=False, colors="WEBLOGO" ):
    """Generates classic **LOGO** plots.

    :param df: Data container.
    :type df: Union[:class:`.DesignFrame`, :class:`.SequenceFrame`]
    :param str seqID: |seqID_param|.
    :param bool refseq: if :data:`True` (default), mark the original residues according to
        the reference sequence.
    :param key_residues: |keyres_param|.
    :type key_residue: |keyres_param|
    :param int line_break: Force a line-change in the plot after n residues are plotted.
    :param int hight_prop: Hight proportion for each row of the plot.
    :param float font_size: Expected size of the axis font.
    :param bool refplot: When :data:`True`, it will res-order the residues in each position
        so that the reference residue will be on the bottom and setting a two-color scheme
        (provide a single color name in ``colors``) that allows to quickly identify the reference
        type in each position.
    :param colors: Colors to assign; it can be the name of a available color set or
        a dictionary with a color for each type. Available color schemes are: Weblogo
        (default), Hydrophobicity, Chemistry, and Charge.
    :type colors: Union[:class:`str`, :class:`dict`]

    :return: :class:`~matplotlib.figure.Figure` and
        :func:`list` of :class:`~matplotlib.axes.Axes`

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: from rstoolbox.plot import logo_plot
           ...: import matplotlib.pyplot as plt
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz",
           ...:                         {"sequence": "B"})
           ...: df.add_reference_sequence("B", df.get_sequence("B")[0])
           ...: fig, axes = logo_plot(df, "B", refseq=True, line_break=50)
           ...: plt.tight_layout()

        @savefig sequence_logo_plot_docs.png width=5in
        In [2]: plt.show()

        In [3]: plt.close()
    """
    def _letterAt( letter, x, y, yscale=1, ax=None, globscale=1.35,
                   LETTERS=None, COLOR_SCHEME=None ):
        text = LETTERS[letter]
        t = mpl.transforms.Affine2D().scale(1 * globscale, yscale * globscale) + \
            mpl.transforms.Affine2D().translate(x, y) + ax.transData
        p = PathPatch(text, lw=0, fc=COLOR_SCHEME[letter],  transform=t)
        if ax is not None:
            ax.add_artist(p)
        return p

    def _dataframe2logo( data ):
        aa = list(data)
        odata = []
        for _, pos in data.iterrows():
            pdata = []
            for k in aa:
                if pos[k] > 0.0000000:
                    pdata.append( ( k, float(pos[k]) ) )
            odata.append(sorted(pdata, key=operator.itemgetter(1, 0)))
        return odata

    def _chunks(l, n):
        for i in range(0, len(l), n):
            yield l[i:i + n]

    order = ["A", "V", "I", "L", "M", "F", "Y", "W", "S", "T", "N",
             "Q", "R", "H", "K", "D", "E", "C", "G", "P"]
    data = copy.deepcopy(df)
    if data.empty:
        raise ValueError("Provided data container is empty. Nothing to plot.")

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
    if isinstance(colors, dict):
        for aa in colors:
            LETTERS[aa] = TextPath((letters_shift, 0), aa, size=1, prop=fp)
    elif isinstance(colors, str):
        for aa in color_scheme(colors):
            LETTERS[aa] = TextPath((letters_shift, 0), aa, size=1, prop=fp)
    else:
        raise ValueError("Colors need to either be a string representing the name of a available "
                         "color set or a dictionary with a color for each type.")

    # Data type management.
    if not isinstance(data, pd.DataFrame):
        raise ValueError("Input data must be in a DataFrame, DesignFrame or SequenceFrame")
    else:
        if not isinstance(data, (DesignFrame, SequenceFrame)):
            if len(set(data.columns.values).intersection(set(order))) == len(order):
                data = SequenceFrame(data)
            else:
                data = DesignFrame(data)
    if isinstance(data, DesignFrame):
        data = data.sequence_frequencies(seqID)

    # key_residues management.
    length = len(data.get_reference_sequence(seqID)) if refseq else None
    key_residues = get_selection(key_residues, seqID, list(data.index.values), length)

    # Plot
    if line_break is None:
        figsize = (len(data) * 2, 2.3 * hight_prop)
        grid = (1, 1)
        fig  = plt.figure(figsize=figsize)
        axs  = [plt.subplot2grid(grid, (0, 0)), ]
        krs  = [key_residues, ]
    else:
        rows = int(math.ceil(float(len(data)) / line_break))
        figsize = (float(len(data) * 2 ) / rows, 2.3 * hight_prop * rows)
        grid = (rows, 1)
        fig  = plt.figure(figsize=figsize)
        axs  = [plt.subplot2grid(grid, (_, 0)) for _ in range(rows)]
        krs  = list(_chunks(key_residues, line_break))

    font = FontProperties()
    font.set_size(font_size)
    font.set_weight('bold')

    for _, ax in enumerate(axs):
        # Refseq and key_residues management.
        ref_seq = data.get_reference_sequence(seqID, krs[_]) if refseq else ""
        # data and key_residues management.
        _data = data.get_key_residues(krs[_])

        maxv = int(math.ceil(data.max_hight()))

        ticks = len(_data)
        if line_break is not None and len(_data) < line_break:
            ticks = line_break
        ax.set_xticks(np.arange(0.5, ticks + 1))
        ax.set_yticks( range(0, maxv + 1) )
        ax.set_xticklabels( _data.index.values )
        ax.set_yticklabels( np.arange( 0, maxv + 1, 1 ) )
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
        wdata = _dataframe2logo( _data )
        x = 0.5
        maxi = 0
        for scores in wdata:
            y = 0
            for base, score in scores:
                if isinstance(colors, dict):
                    _letterAt(base, x, y, score, ax, globscale, LETTERS, colors)
                else:
                    _letterAt(base, x, y, score, ax, globscale, LETTERS, color_scheme(colors))
                y += score
            x += 1
            maxi = max(maxi, y)
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontproperties(font)
        if ref_seq is not None:
            for label in (ax2.get_xticklabels() + ax2.get_yticklabels()):
                label.set_fontproperties(font)

    return fig, axs
