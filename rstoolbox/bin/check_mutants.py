# @Author: Jaume Bonet <bonet>
# @Date:   12-Apr-2018
# @Email:  jaume.bonet@gmail.com
# @Filename: check_mutants.py
# @Last modified by:   bonet
# @Last modified time: 12-Apr-2018

# flake8: noqa

# Standard Libraries
import argparse

# External Libraries
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

# This Library
import rstoolbox
from rstoolbox.io import parse_rosetta_file, read_fasta
from rstoolbox.io import write_clustalw, write_mutant_alignments
from rstoolbox.plot import logo_plot

# Configuration properties
sns.set(font_scale=2)
mpl.rcParams['svg.fonttype'] = 'none'
sns.set_style("whitegrid")
rstoolbox.core.set_option("system", "overwrite", True)


def get_options( *args, **kwds ):

    parser = argparse.ArgumentParser(
        description="Check the sequences in the silent files.")

    parser.add_argument('-in:file', dest='ifile', action='store',
                        help='Input silent file', default=None)
    parser.add_argument('-in:files', dest='ifiles', action='store',
                        help='Prefix of input silent files', default=None)
    parser.add_argument('-in:fasta', dest='ffile', action='store',
                        help='Rerefence sequence (single FASTA)', default=None)
    parser.add_argument('-in:seqID', dest='seqID', action='store',
                        help='Sequence chain identifier (def: A)', default='A')

    parser.add_argument('-out:prefix', dest='ofile', action='store',
                        help='Prefix for the outputs', default=None)
    parser.add_argument('-out:format', dest='iformat', action='store', choices=['png', 'svg'],
                        help='Ouptut images format', default='png')
    parser.add_argument('-out:font', dest='ifont', action='store',
                        help='Font size for logo plot (might need manual adjustment)',
                        default=35)

    options = parser.parse_args()

    if options.ifile is None and options.ifiles is None:
        raise AttributeError("A filename or a prefix for multiple filename have to be provided.")
    if options.ifile is not None and options.ifiles is not None:
        raise AttributeError("Provide only ONE file or a prefix for multiple files, not both.")
    if options.ofile is None:
        options.ofile = "seqcompare"

    return options

def main( options ):
    infile = options.ifile if options.ifile is not None else options.ifiles
    defs = {"sequence": options.seqID}
    df = parse_rosetta_file(infile, defs, multi=options.ifiles is not None)

    if options.ffile is not None:
        refseq = read_fasta(options.ffile).get_sequence("A").values[0]
        df.add_reference_sequence(options.seqID, refseq)

    # Alignment file
    alif = options.ofile + ".clw"
    write_clustalw(df, options.seqID, alif)

    # Mutation list file
    if options.ffile is not None:
        mutf = options.ofile + "_mutants.clw"
        write_mutant_alignments(df, options.seqID, mutf)

    # Logo Plot
    logof = options.ofile + "_logo" + "." + options.iformat
    _, _ = logo_plot(df, options.seqID, refseq=options.ffile is not None,
                     line_break=50, font_size=int(options.ifont) )
    plt.savefig(logof)

    # Alignment plot
    # ...


if __name__ == '__main__':
    main( get_options() )
