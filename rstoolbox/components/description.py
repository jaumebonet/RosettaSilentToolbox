import os
import json

import numpy as np

class Description( object ):
    def __init__( self, ddict ):
        self.keys    = {}
        self.keysbr  = {}
        self.perresi = []
        self.namings = {}
        self.ignores = []
        self.chains  = []
        self.labels  = []
        self._auto   = False
        self._process_description( self._file_vs_json( ddict ) )

    def add_per_residues_keys( self, per_residues ):
        self.perresi = per_residues

    def fill_if_empty_scores( self, header ):
        if self._auto:
            for h in header:
                if self.is_per_residue_key( h ):
                    hh = self.get_per_residue_name( h )
                    if not hh in self.keysbr and not hh in self.ignores:
                        self.keysbr[hh]=hh
                        continue
                if not h in self.keys and not h in self.ignores:
                    self.keys[h] = h

    def is_requested_key( self, key ):
        return key in self.keys

    def is_per_residue_key( self, key ):
        for k in self.perresi:
            if key.startswith(k):
                return True
        return False

    def get_per_residue_name( self, key ):
        for k in self.perresi:
            if key.startswith(k):
                return k

    def is_requested_per_residue_key( self, key ):
        for k in self.keysbr:
            if key.startswith(k):
                return True
        return False

    def is_requested_label( self, label ):
        return label.upper() in self.labels

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
        # Option 0: cover no definition: Then, all scores get picked.
        if ddict is None:
            self._auto = True
            return
        # Option 1: scores is a list of scores to keep with the same name
        # Thus, it can have a list or "*" if all scores are meant to be kept.
        # Alternatively, one can pick only scores to ignore with scores_ignore.
        # Both options have to be incompatible.
        if "scores" in ddict:
            if isinstance(ddict["scores"], list):
                for sc in ddict["scores"]:
                    self.keys[sc] = sc
            else:
                if ddict["scores"] == "*":
                    self._auto = True
        elif "scores_ignore" in ddict:
            self.ignores = ddict["scores_ignore"]
            self._auto = True
        if "scores_ignore" in ddict:
            self.ignores = ddict["scores_ignore"]
        # Option 2: scores_rename is a list of scores to keep with a different name
        if "scores_rename" in ddict:
            for sc in ddict["scores_rename"]:
                self.keys[sc] = ddict["scores_rename"][sc]
        # Option 3: scores_by_residue is the prefix of a score that has to be saved
        # by residue.
        if "scores_by_residue" in ddict:
            for sc in ddict["scores_by_residue"]:
                sk = ddict["scores_by_residue"][sc]
                if self.is_requested_key( sk ):
                    raise KeyError( "The per residue score {} already exists as a score column".format(sk))
                self.keysbr[sc] = sk
        # naming splits the definition by "_" and sets the names as new columns
        if "naming" in ddict:
            for cname, name in enumerate(ddict["naming"]):
                if name != "":
                    if self.is_requested_key( name ):
                        raise KeyError( "The naming split {} already exists as a score column".format(name))
                    self.namings[cname] = name
        # sequence define which sequence(s) the user wants to keep
        if "sequence" in ddict:
            self.chains = ddict["sequence"]
        # labels selects from the REMARK LABELS line which are of interest
        if "labels" in ddict:
            self.labels = [x.upper() for x in ddict["labels"]]
        if "split" in ddict:
            pass
