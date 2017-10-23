import os
import json

import numpy as np

class Description( object ):
    def __init__( self, ddict ):
        self.keys    = {}
        self.keysbr  = {}
        self.namings = {}
        self.chains  = []
        self._process_description( self._file_vs_json( ddict ) )

    def is_requested_key( self, key ):
        return key in self.keys

    def is_requested_per_residue_key( self, key ):
        for k in self.keysbr:
            if key.startswith(k):
                return k
        return False

    def get_expected_per_residue_key( self, key ):
        for k in self.keysbr:
            if key.startswith(k):
                return self.keysbr[k]
        raise KeyError("Per residue key {} not found".format(key))

    def get_expected_key( self, key ):
        return self.keys[key]

    def get_naming_pairs( self, description ):
        if len(self.namings) > 0:
            d = description.split("_")
            c = 0
            for cn, n in enumerate(d):
                if cn in self.namings:
                    c += 1
                    yield self.namings[cn], n
                if c == len(self.namings):
                    break

    def get_expected_sequences( self, chainsdct ):
        guide  = np.array(list(chainsdct["stc"]))
        seq    = chainsdct["seq"]
        for ch in self.chains:
            if ch in chainsdct["id"]:
                guidep = np.where(guide == ch)[0]
                yield "sequence_" + ch, seq[guidep[0]:guidep[-1] + 1]

    def _file_vs_json( self, data ):
        if isinstance( data, str ):
            if not os.path.isfile( data ):
                raise IOError("{0}: file not found.".format(data))
            fd = gzip.open( data ) if data.endswith(".gz") else open( data )
            text = "".join([x.strip() for x in fd])
            data = json.loads(text)
        return data

    def _process_description( self, ddict ):
        # scores is a list of scores to keep with that name
        if "scores" in ddict:
            for sc in ddict["scores"]:
                self.keys[sc] = sc
        # scores_rename is a list of scores to keep with a different name
        if "scores_rename" in ddict:
            for sc in ddict["scores_rename"]:
                self.keys[sc] = ddict["scores_rename"][sc]
        # scores_by_residue is the prefix of a score that has to be saved by residue
        if "scores_by_residue" in ddict:
            for sc in ddict["scores_by_residue"]:
                sk = ddict["scores_by_residue"][sc]
                if self.is_requested_key( sk ):
                    raise KeyError( "The per residue score {} already exists as a score column".format(sk))
                self.keysbr[sc] = sk
        if "naming" in ddict:
            for cname, name in enumerate(ddict["naming"]):
                if name != "":
                    if self.is_requested_key( name ):
                        raise KeyError( "The naming split {} already exists as a score column".format(name))
                    self.namings[cname] = name
        if "sequence" in ddict:
            self.chains = ddict["sequence"]
        if "split" in ddict:
            pass
