#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

Plot any two variables from the Rosetta Silent/Score file.
"""
# Standard Libraries
import argparse

# External Libraries
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

# This Library
from rstoolbox.io import parse_rosetta_file
from rstoolbox.utils import add_top_title

# Configuration properties
sns.set(font_scale=2)
mpl.rcParams['svg.fonttype'] = 'none'
sns.set_style("whitegrid")


def make_parser( *args, **kwds ):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-in:file', dest='ifile', action='store',
                        help='Input silent file (single file)', default=None)
    parser.add_argument('-in:files', dest='ifiles', action='store',
                        help='Prefix of input silent files (multiple files). '
                        'Incompatible with ``-in:file``.', default=None)
    parser.add_argument('-in:y', dest='y', action='store',
                        help='Name of the score for Y axis.', default=None)
    parser.add_argument('-in:x', dest='x', action='store',
                        help='Name of the score for X axis.', default=None)
    parser.add_argument('-plot:title', dest='title', action='store',
                        help='Title of the plot.', default=None)
    parser.add_argument('-plot:color', dest='color', action='store',
                        help='Color of the plot. If number, it refers to '
                        'the seaborn color palette.', default=0  )
    parser.add_argument('-plot:ylab', dest='ylab', action='store',
                        help='Label for Y axis.', default=None)
    parser.add_argument('-plot:xlab', dest='xlab', action='store',
                        help='Label for X axis.', default=None)
    parser.add_argument('-plot:ylim', dest='ylim', action='store',
                        help='Range for Y axis.', nargs=2, default=[None, None])
    parser.add_argument('-plot:xlim', dest='xlim', action='store',
                        help='Range for X axis.', nargs=2,  default=[None, None])
    parser.add_argument('-plot:fsize', dest='fsize', action='store',
                        help='Size of the figure.', nargs=2, default=[None, None])
    parser.add_argument('-out:silent', dest='silent', action='store_true',
                        help='If True, do not plot on screen.', default=False)
    parser.add_argument('-out:file', dest='ofile', action='store',
                        help='Output image file. If None, no file is created.',
                        default=None)
    return parser


def get_options( parser ):

    options = parser.parse_args()

    if options.ifile is None and options.ifiles is None:
        raise AttributeError("A filename or a prefix for multiple filename have to be provided.")
    if options.ifile is not None and options.ifiles is not None:
        raise AttributeError("Provide only ONE file or a prefix for multiple files, not both.")
    if options.y is None:
        raise AttributeError("A score name must be provided for the Y axis.")
    if options.x is None:
        raise AttributeError("A score name must be provided for the X axis.")
    return options


def main( options ):
    infile = options.ifile if options.ifile is not None else options.ifiles
    df = parse_rosetta_file(infile, multi=options.ifiles is not None)

    # Plot
    if options.fsize[0] is None:
        fig = plt.figure()
    else:
        fig = plt.figure(figsize=[float(x) for x in options.fsize])
    ax = plt.subplot2grid((1, 1), (0, 0), fig=fig)

    sns.regplot(x=options.x, y=options.y, data=df, fit_reg=False,
                ax=ax, color=sns.color_palette()[int(options.color)])
    if options.ylim[0] is not None:
        ax.set_ylim(bottom=float(options.ylim[0]), top=float(options.ylim[1]))
    if options.xlim[0] is not None:
        ax.set_xlim(left=float(options.xlim[0]), right=float(options.xlim[1]))
    if options.ylab is not None:
        ax.set_ylabel(options.ylab)
    if options.xlab is not None:
        ax.set_xlabel(options.xlab)

    add_top_title(ax, options.title)

    plt.tight_layout()

    # Write to file
    if options.ofile is not None:
        plt.savefig(options.ofile)

    # Show on screen
    if not options.silent:
        plt.show()
    return fig


if __name__ == '__main__':
    _ = main( get_options( make_parser() ) )
