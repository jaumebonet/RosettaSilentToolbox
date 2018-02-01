# @Author: Jaume Bonet <bonet>
# @Date:   17-Jan-2018
# @Email:  jaume.bonet@gmail.com
# @Filename: plot_fragments_rmsd.py
# @Last modified by:   bonet
# @Last modified time: 01-Feb-2018


import argparse
import gzip
import os

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

from rstoolbox.io import parse_rosetta_fragments
from rstoolbox.plot import plot_fragments
from rstoolbox.utils import add_right_title, add_top_title

sns.set(font_scale=2)
mpl.rcParams['svg.fonttype'] = 'none'
sns.set_style("whitegrid")

def get_options(*args, **kwds):

    parser = argparse.ArgumentParser(description="Generate a plot for fragment RMSD evaluation.")

    parser.add_argument( '-in:frag:small', dest='fsmall', action='store',      help='Name of the small fragments file.',         default=None )
    parser.add_argument( '-in:qual:small', dest='qsmall', action='store',      help='Name of the small fragments quality file.', default=None )
    parser.add_argument( '-in:frag:large', dest='flarge', action='store',      help='Name of the large fragments file.',         default=None )
    parser.add_argument( '-in:qual:large', dest='qlarge', action='store',      help='Name of the large fragments quality file.', default=None )
    parser.add_argument( '-in:pdb',        dest='pdb',    action='store',      help='PDB file (in case we need it).',            default=None )
    parser.add_argument( '-out:silent',    dest='silent', action='store_true', help='If True, do not plot on screen',            default=False  )
    parser.add_argument( '-out:format',    dest='format', action='store',      help='Make plot vertical (v) or horizontal (h).', default="h"  )
    parser.add_argument( '-out:file',      dest='ofile',  action='store',      help='Output image file (.png/.svg).',            default=None )

    options = parser.parse_args()

    if options.fsmall is None:
        raise IOError("A filename has to be provided.")
    if not os.path.isfile(options.fsmall):
        raise IOError("{} not found".format(options.fsmall))

    if options.flarge is None:
        raise IOError("A filename has to be provided.")
    if not os.path.isfile(options.flarge):
        raise IOError("{} not found".format(options.flarge))

    # Loading quality if not provided but exists
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

    ax00 = plt.subplot2grid(grid, (0, 0))
    ax01 = plt.subplot2grid(grid, (0, 1) if options.format == "h" else (1, 0))

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

if __name__ == '__main__':
    main( get_options() )
