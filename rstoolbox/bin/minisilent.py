# @Author: Jaume Bonet <bonet>
# @Date:   12-Oct-2017
# @Email:  jaume.bonet@gmail.com
# @Filename: minisilent.py
# @Last modified by:   bonet
# @Last modified time: 23-Jan-2018


import argparse
import gzip
import os

from rstoolbox.io import open_rosetta_file

def get_options(*args, **kwds):

    parser = argparse.ArgumentParser(description="Transform fully-fledged silent file(s) into portable minisilent")

    parser.add_argument( '-in:file',   dest='ifile',  action='store',      help='Input silent file',            default=None  )
    parser.add_argument( '-in:files',  dest='ifiles', action='store',      help='Prefix of input silent files', default=None  )
    parser.add_argument( '-out:file',  dest='ofile',  action='store',      help='Output minisilent',            default=None  )
    parser.add_argument( '-overwrite', dest='force',  action='store_true', help='Allow overwrite',              default=False )

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
    fd = gzip.open( options.ofile, "w" ) if options.ofile.endswith(".gz") else open( options.ofile, "w" )
    infile = options.ifile if options.ifile is not None else options.ifiles
    for line, is_header, count in open_rosetta_file( infile, options.ifile is None ):
        fd.write( line )

if __name__ == '__main__':
    main( get_options() )
