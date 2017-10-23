import os
import sys
import re
import glob
import gzip
import json
from collections import OrderedDict

import pandas as pd

import rstoolbox.core as core
import rstoolbox.components as cp

_headers = ["SCORE", "REMARK", "RES_NUM", "FOLD_TREE", "RT", "ANNOTATED_SEQUENCE", "NONCANONICAL_CONNECTION", "CHAIN_ENDINGS"]

def _check_type( value ):
    try:
        int(value)
    except ValueError:
        try:
            float(value)
        except ValueError:
            return value
        else:
            return float(value)
    else:
        return int(value)

def open_rosetta_file( filename, multi=False ):

    files = []
    if not multi:
        if not os.path.isfile( filename ):
            raise IOError("{0}: file not found.".format(filename))
        files.append( filename )
    else:
        files = glob.glob( filename + "*" )

    for f in files:
        fd = gzip.open( f ) if f.endswith(".gz") else open( f )
        for line in fd:
            if line.strip().split()[0].strip(":") in _headers:
                yield line, line.strip().split()[-1] == "description"
        fd.close()

def parse_rosetta_file( filename, description, multi=False ):

    desc   = cp.Description( description )
    header = []
    data   = {}
    chains = {"id": "", "seq": "", "stc": "", "done": False}
    for line, is_header in open_rosetta_file( filename, multi ):
        if is_header:
            header = line.strip().split()[1:]
            continue
        if line.startswith("SCORE"):
            chains = {"id": "", "seq": "", "stc": "", "done": False}
            for cv, value in enumerate( line.strip().split()[1:] ):
                if desc.is_requested_key( header[cv] ):
                    data.setdefault( desc.get_expected_key( header[cv]), [] ).append( _check_type( value ) )
            for namingID, namingVL in desc.get_naming_pairs( line.strip().split()[-1] ):
                data.setdefault( namingID, [] ).append( _check_type( namingVL ) )
            continue
        if line.startswith("RES_NUM"):
            chains["id"] = "".join(list(OrderedDict.fromkeys("".join([x.split(":")[0] for x in line.split()[1:-1]]))))
            continue
        if line.startswith("ANNOTATED_SEQUENCE"):
            chains["seq"] = re.sub( r'\[[^]]*\]', '', line.strip().split()[1] )
            if len(chains["id"]) == 1:
                chains["stc"] = "".join([chains["id"]] * len(chains["seq"]))
                for seqname, seq in desc.get_expected_sequences( chains ):
                    data.setdefault( seqname, [] ).append( seq )
                chains["done"] = True
            continue
        if line.startswith("CHAIN_ENDINGS") and not chains["done"]:
            endings = [int(x) for x in line.split()[1:-1]]
            txains = []
            endings.append(len(chains["seq"]))
            for x in range(len(endings)):
                if x > 0: endings[x] -= endings[x-1]
                txains.extend([chains["id"][x]] * endings[x])
            chains["stc"] = "".join(txains)
            for seqname, seq in desc.get_expected_sequences( chains ):
                data.setdefault( seqname, [] ).append( seq )
            continue
    return pd.DataFrame( data )

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
