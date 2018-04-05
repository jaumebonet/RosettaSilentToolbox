# @Author: Jaume Bonet <bonet>
# @Date:   19-Feb-2018
# @Email:  jaume.bonet@gmail.com
# @Filename: description.py
# @Last modified by:   bonet
# @Last modified time: 05-Apr-2018


import re

import six
import numpy as np

__all__ = ["Description"]


class Description( object ):

    _PARAMS = [
        "scores", "scores_ignore", "scores_rename", "scores_by_residue",
        "naming", "sequence", "structure", "psipred", "labels"
    ]

    def __init__( self, **params ):
        # Avoid false attributes and misspellings
        for k in params:
            if k not in self._PARAMS:
                raise AttributeError(
                    "Wrong attribute provided. Accepted keys are: "
                    "{}".format(",".join(self._PARAMS)))

        # Load selections
        self.scores            = params.get("scores",            "*" )
        self.scores_ignore     = params.get("scores_ignore",     None)
        self.scores_rename     = params.get("scores_rename",     None)
        self.scores_by_residue = params.get("scores_by_residue", None)
        self.naming            = params.get("naming",            None)
        self.sequence          = params.get("sequence",          None)
        self.structure         = params.get("structure",         None)
        self.psipred           = params.get("psipred",           None)
        self.labels            = params.get("labels",            None)

        # Others
        self._per_residues = ["residue_ddg_"]
        self._scores       = {}

        # Check that the "per_residue" scores requested are known.
        if self.scores_by_residue is not None:
            for k in self.scores_by_residue:
                if k not in self._per_residues:
                    raise AttributeError("Unknown per-residue score: {}".format(k))

        # Fix user input in case strings instead of lists were provided
        if isinstance(self.scores, six.string_types) and self.scores not in ["*", "-"]:
            self.scores = self.scores.split(",")

        # Negate all scores
        if self.scores_ignore == "*" or self.scores == "-":
            self.scores = None

        # Labels will be set in UPPERCASE
        if self.labels is not None:
            self.labels = [x.upper() for x in self.labels]

    def wanted_score( self, score_name ):
        # skip per-residue values when nod asked for
        for k in self._per_residues:
            if score_name.startswith(k) and score_name not in self.scores:
                return False
        if self.scores_ignore is not None:
            if self.scores_ignore == "*" or score_name in self.scores_ignore:
                return False
            elif any(re.match(s, score_name) for s in self.scores_ignore if "*" in s):
                return False
        if self.scores_rename is not None:
            if score_name in self.scores_rename:
                return True
        if self.scores is None:
            return False
        if self.scores == "*" or score_name in self.scores:
            return True
        return False

    def score_name( self, score_name ):
        if self.scores_rename is not None:
            if score_name in self.scores_rename:
                return self.scores_rename[score_name]
        return score_name

    def wanted_per_residue_score( self, score_name ):
        if self.scores_by_residue is None:
            return False
        for k in self._per_residues:
            if score_name.startswith(k):
                return True

    def check_naming( self, header ):
        if self.naming is None:
            return
        for h in self.naming:
            if h in header and self.wanted_score(h):
                if h == self.score_name(h):
                    raise AttributeError(
                        "New column {} ".format(h) +
                        "cannot be created, a score has that name")

    def get_naming_pairs( self, description ):
        if self.naming is not None:
            d = description.split("_")
            for cn, n in enumerate(d):
                if cn == len(self.naming):
                    break
                if self.naming[cn] != "":
                    yield self.naming[cn], n

    def setup_labels( self, data ):
        if self.labels is not None:
            for l in self.labels:
                data.setdefault("lbl_" + l.upper(), []).append("0")
        return data

    def get_expected_sequences( self, chains ):
        if self.sequence is not None:
            sele = set(chains["id"]) if self.sequence == "*" else set(self.sequence)
            seq  = chains["seq"]
            if len(sele.difference(chains["id"])) > 0:
                raise ValueError(
                    "Requested a chain not present in the file. "
                    "Available chain are {}".format(",".join(list(set(chains["id"])))))
            guide  = np.array(list(chains["id"]))
            for ch in sele:
                guidep = np.where(guide == ch)[0]
                yield "sequence_" + ch, "".join(seq[guidep[0]:guidep[-1] + 1])

    def get_expected_structures( self, chains ):
        if self.structure is not None:
            sele = set(chains["id"]) if self.structure == "*" else set(self.structure)
            seq  = chains["dssp"]
            if len(sele.difference(chains["id"])) > 0:
                raise ValueError(
                    "Requested a chain not present in the file. "
                    "Available chain are {}".format(",".join(list(set(chains["id"])))))
            guide  = np.array(list(chains["id"]))
            for ch in sele:
                guidep = np.where(guide == ch)[0]
                yield "structure_" + ch, "".join(seq[guidep[0]:guidep[-1] + 1])

    def get_expected_psipred( self, chains ):
        if self.psipred is not None:
            sele = set(chains["id"]) if self.psipred == "*" else set(self.psipred)
            seq  = chains["psipred"]
            if len(sele.difference(chains["id"])) > 0:
                raise ValueError(
                    "Requested a chain not present in the file. "
                    "Available chain are {}".format(",".join(list(set(chains["id"])))))
            guide  = np.array(list(chains["id"]))
            for ch in sele:
                guidep = np.where(guide == ch)[0]
                yield "psipred_" + ch, "".join(seq[guidep[0]:guidep[-1] + 1])

    def to_json( self ):
        data = {}
        for a in self._PARAMS:
            if getattr(self, a) is not None:
                data[a] = getattr(self, a)
        return data
