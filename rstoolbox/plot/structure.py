# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: positional_structural_similarity_plot
.. func:: sequence_similarity
.. func:: positional_sequence_similarity
.. func:: binary_similarity
.. func:: binary_overlap
"""
# Standard Libraries

# External Libraries
import seaborn as sns

# This Library

__all__ = ['positional_structural_similarity_plot']


def positional_structural_similarity_plot( df, ax, alpha_color='royalblue',
                                           beta_color='tomato', identity_color='black',
                                           identity_width=3 ):
    """Generates a bar plot for positional prevalence of secondary structure elements.

    Input data can/should be generated with :func:`.positional_structural_count`.

    If there is a ``identity_perc`` column present, which can be obtained by running
    :func:`.positional_structural_identity`, it will also print a line showing how
    often the secondary structure matches the expected/reference one.

    Both :class:`~pandas.DataFrame` obtained through the two functions can be simply
    merged with::

        pd.concat([df1, df2], axis=1)

    :param df: |df_param|, where rows are sequence positions and columns
        represent the secondary structure type: ``H``, ``E`` and ``L``
    :type df: :class:`~pandas.DataFrame`
    :param ax: |axis_param|.
    :type ax: :class:`~matplotlib.axes.Axes`
    :param alpha_color: Color to represent helices.
    :type alpha_color: |color_types|
    :param beta_color: Color to represent beta strands.
    :type beta_color: |color_types|
    :param identity_color: Color to represent the secondary structure
        identity with the expected ``reference_structure``.
    :type identity_color: |color_types|
    :param int identity_width: Width of the line representing the
        identity with the expected ``reference_structure``.

    .. seealso::
        :func:`.positional_structural_count`
        :func:`.positional_structural_identity`

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: from rstoolbox.analysis import positional_structural_count
           ...: from rstoolbox.analysis import positional_structural_identity
           ...: from rstoolbox.plot import positional_structural_similarity_plot
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_ssebig.minisilent.gz",
           ...:                         {'scores': ['score'], 'structure': 'C'})
           ...: df.add_reference_structure('C', df.get_structure('C').values[0])
           ...: df1 = positional_structural_identity(df.iloc[1:], 'C')
           ...: df2 = positional_structural_count(df.iloc[1:], 'C')
           ...: fig = plt.figure(figsize=(35, 10))
           ...: ax00 = plt.subplot2grid((1, 1), (0, 0))
           ...: positional_structural_similarity_plot(pd.concat([df1, df2], axis=1), ax00)
           ...: plt.tight_layout()

        @savefig plot_positional_structural_similarity.png width=5in
        In [2]: plt.show()
    """

    # Color management
    if isinstance(alpha_color, int):
        alpha_color = sns.color_palette()[alpha_color]
    if isinstance(beta_color, int):
        beta_color = sns.color_palette()[beta_color]
    if isinstance(identity_color, int):
        identity_color = sns.color_palette()[identity_color]

    h = df["H"].values
    e = df["E"].values

    ax.bar(range(len(h)), h, color=alpha_color )
    ax.bar(range(len(h)), e, bottom=h, color=beta_color)
    if "identity_perc" in df:
        ax.plot(range(len(h)), df["identity_perc"].values,
                linestyle="solid", linewidth=str(identity_width),
                color=identity_color)

    ax.set_ylim(0, 1)
    ax.set_xlim(-0.5, len(h))
    ax.set_xticks(range(0, len(h), 5))
    ax.set_xticklabels(range(df.index.values[0],
                             df.index.values[-1], 5))
