# @Author: Jaume Bonet <bonet>
# @Date:   29-Jun-2017
# @Email:  jaume.bonet@gmail.com
# @Filename: DesignList.py
# @Last modified by:   bonet
# @Last modified time: 29-Jun-2017

import os
import re
from collections import Counter

class DesignList( list ):
    def __init__( self ):
        self._index = {}

    def add_design( self, design ):
        if design.sanity_check():
            self.append( design )
            self._index[ design.tag() ] = len( self ) - 1

    def has_design( self, design_id ):
        return design_id in self._index

    def get_design( self, design_id ):
        return self[ self._index[ design_id ] ]

    def filter_by_tags_file( self, filename ):
        if not os.path.isfile( filename ):
            raise IOError("Cannot find file {0}\n".format( filename ) )

        tags = []
        with open( filename ) as fd:
            tags = [ line.strip() for line in fd ]
        return self.filter_by_tags( tags )

    def filter_by_tags ( self, tags ):
        if len( tags ) == 0 : return DesignList()
        newList = DesignList()
        for tag in tags:
            if self.has_design( tag ):
                newList.add_design( self.get_design( tag ) )
            else:
                print( "{0} tag is not in the DesignList\n".format( tag ) )
        return newList

    def sanity_check( self ):
        ids = Counter([ _.tag() for _ in self ])
        repeats = 0
        for _ in ids:
            if ids[_] > 1:
                print("Decoy id {0} has {1} repetitions".format( _, ids[_] ) )
                repeats += 1
        if repeats > 0:
            print( "Fix your input file to manage repeated decoy ids!!\n" )
