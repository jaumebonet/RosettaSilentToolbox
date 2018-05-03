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
from rstoolbox.components import DesignSeries
from rstoolbox.utils import add_top_title

# External Libraries
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import colors

# This Library

__all__ = ['positional_structural_similarity_plot', 'plot_ramachandran']


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

def plot_ramachandran( df, seqID, fig):
    """Generates a ramachandran plot in RAMPAGE style.

    For more details and sources please refert to:

        https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/

    Background distribution data taken from:

        https://github.com/S-John-S/Ramachandran-Plot.git

    The phi - psi dihedrals should be present in the DesignSerie. If this is not the case,
    consider computing them using for example the :func:`.get_sequence_and_structure`.
    Note that this function will only plot the ramachandran for a single decoy. If one
    would like to compute it for mutiple decoys, please see the example below.

    :param df: |df_param|, where ONE column cotains the phi and a second column
        the psi angles.
    :type df: :class:`~pandas.Series`
    :param ax: |axis_param|.
    :type ax: :class:`~matplotlib.axes.Axes`

    .. seealso::
        :func:`.get_sequence_and_structure`
        :func:`.get_dihedrals`
        :func:`.get_phi`
        :func:`.get_psi`

    .. rubric:: Example

    .. ipython::

        In [1]: import rstoolbox as rb
           ...: import pandas as pd
           ...: plt.style.use('ggplot')
           ...: pylab.rcParams["figure.figsize"] = (20, 5)
           ...: definitions = {
                               "scores": ["score"],
                               "sequence" : "A",
                               "psipred" : "*",
                               "structure" : "*",
                               "dihedrals": "*"
                               }
           ...: dsf = rb.io.parse_rosetta_file(
                    "~/RosettaSilentToolbox/rstoolbox/tests/data/input_3ssepred.minisilent.gz",
                    definitions )
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_ssebig.minisilent.gz",
           ...:                         {'scores': ['score'], 'structure': 'C'})
           ...: figure = plt.figure(figsize=(15,10))
           ...: rb.plot.plot_ramachandran(dsf.iloc[0], "A", figure)
           ...: plt.tight_layout()
           ...: fig.subplots_adjust(top=1.2)

        @savefig plot_ramachandran.png width=5in
        In [2]: plt.show()
    """
    # Data type management.
    if not isinstance(df, pd.Series):
        raise ValueError("Input data must be in a Series or DesignSeries")
    if not isinstance(df, DesignSeries):
        df = DesignSeries(df)
    if not isinstance(df.get_phi(seqID), np.ndarray):
        raise ValueError(
        "Ramachandran plot function can only be applied on one decoy at once."
        )
    if not isinstance(df.get_psi(seqID), np.ndarray):
        raise ValueError(
        "Ramachandran plot function can only be applied on one decoy at once."
        )

    # General variable for the background preferences.
    cwd = os.path.dirname(__file__)
    rama_preferences = {
        "General": {
            "file": os.path.join(cwd, "rama_bgdists/pref_general.data"),
            "cmap": colors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
            "bounds": [0, 0.0005, 0.02, 1],
        },
        "GLY": {
            "file": os.path.join(cwd, "rama_bgdists/pref_glycine.data"),
            "cmap": colors.ListedColormap(['#FFFFFF', '#FFE8C5', '#FFCC7F']),
            "bounds": [0, 0.002, 0.02, 1],
        },
        "PRO": {
            "file": os.path.join(cwd, "rama_bgdists/pref_proline.data"),
            "cmap": colors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
            "bounds": [0, 0.002, 0.02, 1],
        },
        "PRE-PRO": {
            "file": os.path.join(cwd, "rama_bgdists/pref_preproline.data"),
            "cmap": colors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
            "bounds": [0, 0.002, 0.02, 1],
        }
    }

    # Read in the expected torsion angles.
    __location__ = os.path.join(cwd, "rama_bgdists/")
    rama_pref_values = {}
    for key, val in rama_preferences.items():
        rama_pref_values[key] = np.full((360, 360), 0, dtype=np.float64)
        with open(os.path.join(__location__, val["file"])) as fn:
            for line in fn:
                if not line.startswith("#"):
                    # Preference file has values for every second position only
                    rama_pref_values[key][int(float(line.split()[1])) + 180][int(float(line.split()[0])) + 180] = float(
                        line.split()[2])
                    rama_pref_values[key][int(float(line.split()[1])) + 179][int(float(line.split()[0])) + 179] = float(
                        line.split()[2])
                    rama_pref_values[key][int(float(line.split()[1])) + 179][int(float(line.split()[0])) + 180] = float(
                        line.split()[2])
                    rama_pref_values[key][int(float(line.split()[1])) + 180][int(float(line.split()[0])) + 179] = float(
                        line.split()[2])

    # Ramachandran residue classification.
    seq = df.get_sequence(seqID)
    rama_types = {
                  "GLY": [],
                  "PRO": [],
                  "PRE-PRO": [],
                  "GENERAL": []
    }
    for i,aa in enumerate(seq):
        if aa == "G":
            rama_types["GLY"].append(i)
        elif aa == "P":
            rama_types["PRO"].append(i)
        elif i+1 < len(seq) and seq[i+1] == "P":
            rama_types["PRE-PRO"].append(i)
        else:
            rama_types["GENERAL"].append(i)

    # Generate the plots
    all_phi = df.get_phi(seqID)
    all_psi = df.get_psi(seqID)
    order = ["GENERAL", "GLY", "PRE-PRO", "PRO"]
    grid = (2, 2)
    ax = [plt.subplot2grid(grid, (0, 0), fig=fig),
          plt.subplot2grid(grid, (0, 1), fig=fig),
          plt.subplot2grid(grid, (1, 0), fig=fig),
          plt.subplot2grid(grid, (1, 1), fig=fig)]
    for i, (key, val) in enumerate(sorted(rama_preferences.items(), key=lambda x: x[0].lower())):
        phi = all_phi[rama_types[order[i]]]
        psi = all_psi[rama_types[order[i]]]

        ax[i].imshow(rama_pref_values[key],
                  cmap=rama_preferences[key]["cmap"],
                  norm=colors.BoundaryNorm(rama_preferences[key]["bounds"],
                  rama_preferences[key]["cmap"].N),
                  extent=(-180, 180, 180, -180))
        ax[i].scatter(phi,
                      psi,
                      color="black")
        add_top_title( ax[i], order[i])
        ax[i].set_xlim([-180, 180])
        ax[i].set_ylim([-180, 180])
        ax[i].plot([-180, 180], [0, 0], color="black")
        ax[i].plot([0, 0], [-180, 180], color="black")
        ax[i].locator_params(axis='x', nbins=7)
        ax[i].set_xlabel(r'$\phi$ [degrees]')
        ax[i].set_ylabel(r'$\psi$ [degrees]')
        ax[i].grid()
    #plt.tight_layout()
    #fig.subplots_adjust(top=1.2)
