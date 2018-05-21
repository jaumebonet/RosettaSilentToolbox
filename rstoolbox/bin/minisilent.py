#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

Transform fully-fledged silent file(s) into portable minisilent
by storing only the content that can be processed by the ``rstoolbox``.

It is very convenient to store mid-step data in git repositories.
"""
# Standard Libraries
import argparse
import gzip
import os

# External Libraries

# This Library
from rstoolbox.io import open_rosetta_file


def make_parser( *args, **kwds ):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-in:file', dest='ifile', action='store',
                        help='Input silent file (single entry).', default=None)
    parser.add_argument('-in:files', dest='ifiles', action='store',
                        help='Prefix of input silent files (multiple files).'
                        'Incompatible with ``-in:file``.', default=None)
    parser.add_argument('-out:file', dest='ofile', action='store',
                        help='Output name for minisilent (can be gzipped).',
                        default=None)
    parser.add_argument('-overwrite', dest='force', action='store_true',
                        help='Allows overwriting existing file.', default=False)
    return parser


def get_options( parser ):

    options = parser.parse_args()

    if options.ifile is None and options.ifiles is None:
        raise AttributeError("A filename or a prefix for multiple filename have to be provided.")
    if options.ifile is not None and options.ifiles is not None:
        raise AttributeError("Provide only ONE file or a prefix for multiple files, not both.")
    if options.ofile is None:
        raise AttributeError("Output filename must be provided")
    if os.path.isfile( options.ofile ) and not options.force:
        raise IOError("File {0} exists and will not be overwritten.".format( options.ofile ) )
    return options


def main( options ):
    ofile = options.ofile
    is_gz = ofile.endswith(".gz")
    fd = gzip.open( ofile, "wb" ) if is_gz else open( ofile, "w" )
    infile = options.ifile if options.ifile is not None else options.ifiles
    for line, _, _, _ in open_rosetta_file( infile, options.ifile is None, check_symmetry=False ):
        fd.write( line.encode('utf-8') if is_gz else line )


if __name__ == '__main__':
    main( get_options( make_parser() ) )
