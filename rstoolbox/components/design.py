# @Author: Jaume Bonet <bonet>
# @Date:   29-Jun-2017
# @Email:  jaume.bonet@gmail.com
# @Filename: Design.py
# @Last modified by:   bonet
# @Last modified time: 30-Jun-2017

import os
import re

class Design( object ):
    def __init__( self, design_type="decoy" ):
        self._scores   = {}
        self._scores.setdefault( "design_type", [] ).append( design_type )
        self._sequence = {}
        self._chains   = ""
        self._endings  = []

    def tag( self ):
        return self[ "description" ]

    def add( self, key, value ):
        self._scores.setdefault( key, [] ).append( self._check_type( value ) )

    def alter( self, key, value ):
        self._scores[key][0] = value

    def add_sequence( self, sequence ):
        if "[" in sequence:  # is annotated sequence
            self._sequence["annotated"] = sequence
        self._sequence["sequence"] = re.sub( r'\[[^]]*\]', '', sequence )
        assert len(self._chains) == len(self._sequence["sequence"])

    def get_sequence( self, seq_format="string", chains=None, selection=None ):
        assert seq_format in ["string", "list"]
        seq = self._sequence["sequence"]
        if chains is not None:
            nseq = []
            for rc, r in enumerate(seq):
                if self._chains[rc] in chains:
                    nseq.append(r)
            seq = "".join(nseq)
        if selection is not None:
            nseq = []
            for cx, x in enumerate(selection):
                nseq.append(seq[x - 1])
            seq = "".join(nseq)
        if seq_format == "list":
            seq = seq.split()
        return seq

    def get_annotated_sequence( self, seq_format="string", chains=None ):
        raise NotImplementedError  #TODO

    def add_chains( self, in_chains ):
        chains = []
        for _ in in_chains:
            chain  = _.split(":")[0]
            length =  _.split(":")[1]
            if not "-" in length:
                length = 1
            else:
                length = int( length.split("-")[1] ) - int( length.split("-")[0] ) + 1
            chains.extend( [chain] * length )
        self._chains = "".join( chains )

    def set_endings( self, endings ):
        self._endings = [self._check_type(_) - 1 for _ in endings]

    def get_endings( self ):
        return self._endings

    def has( self, key ):
        return key in self._scores

    def per_residue( self, energy ):
        residue_h = []
        if energy == "ddg":
            residue_h = [ k for k in self._scores if k.startswith( "residue_ddg" ) ]
        if len( residue_h ) == 0:
            return []
        residue_v = [ 0 ] * len( residue_h )
        for _ in residue_h:
            position = int( _.split( "_" )[-1] ) - 1
            residue_v[ position ] = self[_]
        return residue_v

    def sanity_check( self ):
        for _ in self:
            if len( self._scores[_] ) > 1:
                return False
        return True

    def __iter__( self ):
        for k in self._scores:
            yield k

    def __getitem__( self, key ):
        return self._scores[key][0]

    def __str__( self ):
        return self.description()

    def _check_type( self, value ):
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
