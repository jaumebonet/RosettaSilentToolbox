#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

Generate a plot for fragment RMSD evaluation. Requires the output from the
``r_fraq_qual`` application in **Rosetta**.

.. note::
    Depends on :ref:`rosetta.path <options>` and :ref:`rosetta.compilation <options>`,
    if the quality file is not provided.

.. seealso::
    :meth:`.FragmentFrame.add_quality_measure`.
"""
# Standard Libraries
import argparse
import os

# External Libraries
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

# This Library
from rstoolbox.io import parse_rosetta_fragments
from rstoolbox.plot import plot_fragments

# Configuration properties
sns.set(font_scale=2)
mpl.rcParams['svg.fonttype'] = 'none'
sns.set_style("whitegrid")


def make_parser( *args, **kwds ):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-in:frag:small', dest='fsmall', action='store',
                        help='Name of the small fragments file.', default=None)
    parser.add_argument('-in:qual:small', dest='qsmall', action='store',
                        help='Name of the small fragments quality file.', default=None)
    parser.add_argument('-in:frag:large', dest='flarge', action='store',
                        help='Name of the large fragments file.', default=None)
    parser.add_argument('-in:qual:large', dest='qlarge', action='store',
                        help='Name of the large fragments quality file.', default=None)
    parser.add_argument('-in:pdb', dest='pdb', action='store',
                        help='PDB file. If quality files are not provided, this input is '
                        'needed so that ``r_fraq_qual`` can be executed.', default=None)
    parser.add_argument('-out:silent', dest='silent', action='store_true',
                        help='If True, do not plot on screen.', default=False)
    parser.add_argument('-out:format', dest='format', action='store',
                        help='Make plot vertical (v) or horizontal (h).', default="h")
    parser.add_argument('-out:file', dest='ofile',  action='store',
                        help='Output image file (.png/.svg). If None, no image is created.',
                        default=None)
    return parser


def get_options( parser ):

    options = parser.parse_args()

    if options.fsmall is None:
        raise IOError("A filename has to be provided.")
    if not os.path.isfile(options.fsmall):
        raise IOError("{} not found".format(options.fsmall))

    if options.flarge is None:
        raise IOError("A filename has to be provided.")
    if not os.path.isfile(options.flarge):
        raise IOError("{} not found".format(options.flarge))

    # loading quality files if they exist with the appropiate name
    if options.qsmall is None and os.path.isfile(options.fsmall + ".qual" ):
        options.qsmall = options.fsmall + ".qual"
    if options.qlarge is None and os.path.isfile(options.flarge + ".qual" ):
        options.qlarge = options.flarge + ".qual"

    if options.format not in ["v", "h"]:
        raise AttributeError("-out:format can only accept as value 'v' or 'h'.")

    return options


def main( options ):
    # Read Fragment Files
    small_f = parse_rosetta_fragments(options.fsmall)
    large_f = parse_rosetta_fragments(options.flarge)

    # Read or calculate Fragment Quality
    small_f = small_f.add_quality_measure(options.qsmall, options.pdb)
    large_f = large_f.add_quality_measure(options.qlarge, options.pdb)

    # Plot
    fig  = plt.figure(figsize=(40, 10) if options.format == "h" else (20, 20))
    grid = (1, 2) if options.format == "h" else (2, 1)

    ax00 = plt.subplot2grid(grid, (0, 0), fig=fig)
    ax01 = plt.subplot2grid(grid, (0, 1) if options.format == "h" else (1, 0), fig=fig)

    ax00.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax01.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    plot_fragments(small_f, large_f, small_ax=ax00, large_ax=ax01,
                   showfliers=False, titles="top" if options.format == "h" else "right")

    if options.format == "h":
        plt.tight_layout(pad=2)
    else:
        plt.tight_layout(rect=(0.037, 0, 1, 1))

    # Write to file
    if options.ofile is not None:
        plt.savefig(options.ofile)

    # Show on screen
    if not options.silent:
        plt.show()
    return fig


if __name__ == '__main__':
    _ = main( get_options( make_parser() ) )
