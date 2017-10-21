import os
import sys
import glob
import gzip
import json

import pandas as pd

import rstoolbox.core as core

_headers = ["SCORE", "REMARK", "RES_NUM", "FOLD_TREE", "RT", "ANNOTATED_SEQUENCE", "NONCANONICAL_CONNECTION", "CHAIN_ENDINGS"]

def _file_vs_json( data ):
    if isinstance( data, str ):
        if not os.path.isfile( data ):
            raise IOError("{0}: file not found.".format(data))
        fd = gzip.open( data ) if data.endswith(".gz") else open( data )
        data = json.loads([x.strip() for x in fd])
    return data

def open_rosetta_file( filename, multi=False ):

    files = []
    if not multi:
        if not os.path.isfile( filename ):
            raise IOError("{0}: file not found.".format(filename))
        files.append( filename )
    else:
        files = glob.glob( filename + "*" )
    print files

    for f in files:
        fd = gzip.open( f ) if f.endswith(".gz") else open( f )
        for line in fd:
            if line.strip().split()[0].strip(":") in _headers:
                yield line
        fd.close()

def parse_rosetta_file( filename, description, multi=False ):

    description = _file_vs_json( description )


def make_structures( data, silentfiles, column="description", outdir=None, multi=False, tagsfilename="tags", keep_tagfile=True ):

    if outdir is None:
        outdir = core.get_option("system", "output")
    if not os.path.isdir( outdir ):
        os.mkdir( outdir )
    if not outdir.endswith("/"):
        outdir += "/"
    tagsfilename = os.path.join( outdir, tagsfilename )
    if os.path.isfile( tagsfilename ) and not core.get_option("system", "overwrite"):
        raise IOError("Filename {0} already exists and cannot be overwrite.".format(tagsfilename))

    exe = os.path.join( core.get_option("rosetta", "path"), "extract_pdbs.{0}".format(core.get_option("rosetta", "compilation")))
    if not os.path.isfile(exe):
        raise IOError("The expected Rosetta executable {0} is not found".format(exe))

    data[[column]].to_csv( tagsfilename, index=False, header=False)
    if not os.path.isfile(exe):
        raise IOError("Something went wrong writing the file {0}".format(tagsfilename))

    sfiles = []
    if not multi:
        if not os.path.isfile( silentfiles ):
            raise IOError("The expected silent file input {0} is not found".format(silentfiles))
        sfiles.append( silentfiles )
    else:
        sfiles = glob.glob( silentfiles + "*" )
    if len(sfiles) == 0:
        raise IOError("No files found with the pattern {0}".format(silentfiles))

    sfiles = " ".join(sfiles)
    command = "{0} -in:file:silent {1} -in:file:tagfile {2} -out:prefix {3}".format( exe, sfiles, tagsfilename, outdir )
    sys.stdout.write("Executing Rosetta's extract_pdbs app\n")
    sys.stdout.write("(depending on the total number of decoys and how many have been requested this might take a while...)\n")
    error = os.system( command )
    if not bool(error):
        sys.stdout.write("Execution has finished\n")
    else:
        sys.stdout.write("Execution has failed\n")

    if not keep_tagfile:
        os.unlink( tagsfilename )
