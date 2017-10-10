# @Author: Jaume Bonet <bonet>
# @Date:   29-Jun-2017
# @Email:  jaume.bonet@gmail.com
# @Filename: input.py
# @Last modified by:   bonet
# @Last modified time: 10-Jul-2017

import os
import gzip

from ..components import *

def read_silent_file( filename, design_type="", allow_repeats=False ):

    if not os.path.isfile( filename ):
        raise IOError( "Cannot find file {0}\n".format( filename ) )

    dsgList = read_score_file( filename, design_type, allow_repeats )
    if filename.endswith('.gz'):
        fd = gzip.open( filename )
    else:
        fd = open( filename )
    for line in fd:
        if line.startswith("ANNOTATED_SEQUENCE"):
            spline = line.strip().split()[1:]
            if dsgList.has_design( spline[-1] ):
                dsgList.get_design( spline[-1] ).add_sequence( spline[0] )
        if line.startswith("RES_NUM"):
            spline = line.strip().split()[1:]
            if dsgList.has_design( spline[-1] ):
                dsgList.get_design( spline[-1] ).add_chains( spline[:-1] )
        if line.startswith("CHAIN_ENDINGS"):
            spline = line.strip().split()[1:]
            if dsgList.has_design( spline[-1] ):
                dsgList.get_design( spline[-1] ).set_endings( spline[:-1] )
    fd.close()
    return dsgList

def read_score_file( filename, design_type="", allow_repeats=False ):

    if not os.path.isfile( filename ):
        raise IOError( "Cannot find file {0}\n".format( filename ) )

    dsgList = DesignList()
    if filename.endswith('.gz'):
        fd = gzip.open( filename )
    else:
        fd = open( filename )
    headers = []
    for line in fd:
        if line.startswith("SCORE"):
            spline = line.strip().split()[1:]
            if spline[-1] == "description":
                headers = spline
            else:
                if spline[0] == "-nan" or spline[0] == "nan": continue
                dsg = Design( design_type )
                for _ in range( len( spline ) ):
                    dsg.add( headers[_], spline[_] )
                dsgList.add_design( dsg )
    dsgList.sanity_check( allow_repeats )
    fd.close()
    return dsgList
