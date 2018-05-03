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
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import colors

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

def _correct_phi_psi( df ):
    """Corrects the phi psi angles computed (i.e. puts them into [-180, 180])"""
    phi_corr = []
    psi_corr = []
    for phi, psi in zip(df["phi"], df["psi"]):
        # phi angles
        if  phi > 180.:
            phi_corr.append(phi - 360)
        elif phi < -180.:
            phi_corr.append(phi + 360)
        else:
            phi_corr.append(phi)
        # psi angles
        if psi > 180.:
            psi_corr.append(psi - 360)
        elif psi < -180.:
            psi_corr.append(psi + 360)
        else:
            psi_corr.append(psi)
    df["phi_corr"] = phi_corr
    df["psi_corr"] = psi_corr
    return df

def _rama_grid( df ):
    """Categorizing into 4 types (general, proline, glycine, pre-proline) and adds it to the dataframe."""
    rama_types = []
    for i,aa in enumerate(df["amino acid"]):
        if aa == "G":
            rama_types.append("GLY")
        elif aa == "P":
            rama_types.append("PRO")
        elif i+1 < len(df) and df["amino acid"][i+1] == "P":
            rama_types.append("PRE-PRO")
        else:
            rama_types.append("General")
    df["rama_types"] = rama_types
    return df

def plot_ramachandran( df, title=None, save_figure=None, prefix="",
                       expected_torsion_angles="rama_bgdists/" ):
    """Plots ramachandran map in RAMPAGE style.

        :input df:      A dataframe with "amino acids" (single letter codes), "phi", "psi" columns
        :return:        A ramachandran plot in RAMPAGE style.
        """
    # General variable for the background preferences
    rama_preferences = {
        "General": {
            "file": "phi_psi_pref/pref_general.data",
            "cmap": colors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
            "bounds": [0, 0.0005, 0.02, 1],
        },
        "GLY": {
            "file": "phi_psi_pref/pref_glycine.data",
            "cmap": colors.ListedColormap(['#FFFFFF', '#FFE8C5', '#FFCC7F']),
            "bounds": [0, 0.002, 0.02, 1],
        },
        "PRO": {
            "file": "phi_psi_pref/pref_proline.data",
            "cmap": colors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
            "bounds": [0, 0.002, 0.02, 1],
        },
        "PRE-PRO": {
            "file": "phi_psi_pref/pref_preproline.data",
            "cmap": colors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
            "bounds": [0, 0.002, 0.02, 1],
        }
    }

    # Read in the expected torsion angles
    __location__ = expected_torsion_angles
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
    # Check dataframe compatiblility
    if ("phi_corr" or "psi_corr") not in df.columns:
        df = _correct_phi_psi(df)
    if "rama_types" not in df.columns:
        df = _rama_grid(df)

    # Generate the plots
    fig = plt.figure(figsize=(15,10))
    for idx, (key, val) in enumerate(sorted(rama_preferences.items(), key=lambda x: x[0].lower())):
        phi = df[df["rama_types"] == key]["phi_corr"]
        psi = df[df["rama_types"] == key]["psi_corr"]

        plt.subplot(2, 2, idx + 1)

        plt.title(key)
        plt.imshow(rama_pref_values[key], cmap=rama_preferences[key]["cmap"],
                   norm=colors.BoundaryNorm(rama_preferences[key]["bounds"], rama_preferences[key]["cmap"].N),
                   extent=(-180, 180, 180, -180))
        plt.scatter(phi, psi, color="black")
        #plt.scatter(outliers[key]["x"], outliers[key]["y"], color="red")
        plt.xlim([-180, 180])
        plt.ylim([-180, 180])
        plt.plot([-180, 180], [0, 0], color="black")
        plt.plot([0, 0], [-180, 180], color="black")
        plt.locator_params(axis='x', nbins=7)
        plt.xlabel(r'$\phi$ [degrees]')
        plt.ylabel(r'$\psi$ [degrees]')
        plt.grid()
    plt.tight_layout()
    if title:
        fig.subplots_adjust(top=0.93)
        plt.suptitle("{}".format(title))

    if save_figure:
        plt.savefig("{}_rama.png".format(prefix), format="png", dpi=300, transparent=True, bbox_inches="tight")
    #plt.show()
