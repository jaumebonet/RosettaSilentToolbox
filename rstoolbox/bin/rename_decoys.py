#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

Renames decoys with a new prefix identifier and the same numbering
schema that Rosetta uses.
"""
# Standard Libraries
import argparse
import gzip
import os

# External Libraries

# This Library
from rstoolbox.io import parse_rosetta_file


def make_parser( *args, **kwds ):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-in:file', dest='ifile', action='store',
                        help='Input silent file.', default=None )
    parser.add_argument('-prefix', dest='prefix', action='store',
                        help='Prefix for the new naming schema.', default=None )
    parser.add_argument('-out:file', dest='ofile', action='store',
                        help='Output silent file.', default=None )
    parser.add_argument('-overwrite', dest='force', action='store_true',
                        help='Allow overwrite if output file exists.', default=False)
    return parser


def get_options( parser ):

    options = parser.parse_args()

    if options.ifile is None:
        raise AttributeError("A filename has to be provided.")
    if options.ofile is None:
        raise AttributeError("Output filename must be provided")
    if options.prefix is None:
        raise AttributeError("A prefix for the naming schema needs to be provided")
    if os.path.isfile( options.ofile ) and not options.force:
        raise IOError("File {0} exists and will not be overwritten.".format( options.ofile ) )
    return options


def new_names( count, prefix ):
    return "{0}_{1:05d}".format(prefix, count)


def main( options ):
    # Get names and make new ones
    names = parse_rosetta_file( options.ifile, {"scores": ["description"]} )
    names["count"] = names.index + 1
    names["new"]   = names.apply( lambda row: new_names(row["count"], options.prefix), axis=1 )

    # Load the silentfile and change names
    is_gz = options.ifile.endswith(".gz")
    fd   = gzip.open( options.ifile ) if is_gz else open( options.ifile )
    if not is_gz:
        data = "".join(fd.readlines())
    else:
        data = "".join([_.decode('utf8') for _ in fd.readlines()])
    fd.close()

    for _, row in names.iterrows():
        data = data.replace(row["description"], row["new"])

    # Save and write
    is_gz = options.ofile.endswith(".gz")
    fd = gzip.open( options.ofile, "w" ) if is_gz else open( options.ofile, "w" )
    fd.write(data.encode('utf-8') if is_gz else data)
    fd.close()


if __name__ == '__main__':
    main( get_options( make_parser() ) )
