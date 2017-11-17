# @Author: Jaume Bonet <bonet>
# @Date:   29-Jun-2017
# @Email:  jaume.bonet@gmail.com
# @Filename: design.py
# @Last modified by:   bonet
# @Last modified time: 17-Nov-2017

# Standard Libraries
import os
import re
import string

# External Libraries
import pandas as pd

# This Library
from .sequenceFrame import SequenceFrame


class DesignFrame( pd.DataFrame ):
    """
    The :py:class:`.DesignFrame` extends :py:class:`pandas.DataFrame`
    adding some functionalities in order to improve its usability in
    the analysis of sets of design decoys.
    Filled through the functions provided through this library, each
    row represents a decoy while each column represents the scores
    attached to it.

    """
    _metadata = ['_reference_sequence']
    _reference_sequence = {}

    def reference_sequence( self, seqID, sequence=None ):
        """
        Setter/Getter for a reference sequence attached to a particular
        sequence ID.

        :param str seqID: Identifier of the reference sequence
        :param str sequence: Reference sequence. By default is
            :py:data:`None`, which turns the function into a getter.
        :return: str
        """
        if sequence is not None:
            self._reference_sequence[seqID] = sequence
        else:
            if seqID not in self._reference_sequence:
                raise KeyError("There is no reference sequence with ID: {}\n".format(seqID))
        return self._reference_sequence[seqID]

    def has_reference_sequence( self, seqID ):
        """
        Checks if there is a reference sequence for the provided
        sequence ID.
        :param str seqID: Identifier of the reference sequence
        :return: bool
        """
        return seqID in self._reference_sequence

    def get_sequence_with( self, column, selection, confidence=1 ):
        """
        Selects those decoys with a particular set of residue matches.
        Basically, is meant to find, for example, all the decoys in which
        position 25 is A and position 46 is T.

        :param str column: Name of the target column containing sequences
        :param str selection: List of tuples with position and residue
            type (in 1 letter code)
        :param float confidence: Percentage of the number of the selection
            rules that we expect the matches to fulfill. Default is 1 (all).
        :return: Filtered :py:class:`.DesignFrame`
        """
        def match_residues( value, selection, confidence ):
            t = 0
            for s in selection:
                t += 1 if value[s[0]-1] == s[1] else 0
            return t/float(len(selection)) >= float(confidence)
        return self.loc[ self.apply(
                    lambda row: match_residues(row[column], selection, confidence ),
                    axis=1
                )]

    def sequence_frequencies( self, seqID, seqType="protein", cleanExtra=True, cleanUnused=False ):
        """
        Generates a :py:class:`.SequenceFrame` for the frequencies of
        the sequences in the __DesignFrame__ with seqID identifier.
        If there is a reference_sequence for this seqID, it will also
        be attached to the __SequenceFrame__.
        All letters in the sequence will be capitalized. All symbols that
        do not belong to string.ascii_uppercase will be transformed to "*"
        as this is the symbol recognized by the substitution matrices.

        :param str seqID: Identifier of the sequence sets of interest.
        :param str seqType: Type of sequence: protein, dna, rna.
        :param bool cleanExtra: Remove from the SequenceFrame the non-regular
        amino/nucleic acids if they are empty for all positions.
        :param bool cleanUnused: Remove from the SequenceFrame the regular
        amino/nucleic acids if they are empty for all positions
        :return: :py:class:`.SequenceFrame`
        """
        sserie = self["sequence_{0}".format(seqID)].values
        table, extra = self._get_sequence_table( seqType )
        for x in range(len(sserie[0])):
            for k in table:
                table[k].append(float(0))
            for y in range(len(sserie)):
                aa = sserie[y][x].upper()
                if aa in string.whitespace or aa in string.punctuation:
                    aa = "*"
                table[aa][-1] += float(1)
        for k in table:
            for x in range(len(table[k])):
                if table[k][x] != 0:
                    table[k][x] /= float(len(sserie))

        df = SequenceFrame(table)
        df.measure("frequency")
        df.extras( extra )
        if self.has_reference_sequence(seqID):
            df.reference_sequence(self.reference_sequence(seqID))
        df.delete_extra( cleanExtra )
        df.delete_empty( cleanUnused )
        df.index = df.index + 1
        return df

    def sequence_bits( self, seqID, seqType="protein", cleanExtra=True, cleanUnused=False ):
        """
        Generates a :py:class:`.SequenceFrame` for the bits of
        the sequences in the __DesignFrame__ with seqID identifier.
        If there is a reference_sequence for this seqID, it will also
        be attached to the __SequenceFrame__.
        :param str seqID: Identifier of the sequence sets of interest.
        :param str seqType: Type of sequence: protein, dna, rna.
        :param bool cleanExtra: Remove from the SequenceFrame the non-regular
        amino/nucleic acids if they are empty for all positions.
        :param bool cleanUnused: Remove from the SequenceFrame the regular
        amino/nucleic acids if they are empty for all positions
        :return: :py:class:`.SequenceFrame`
        """
        sserie = self["sequence_{0}".format(seqID)].values
        table, extra = self._get_sequence_table( seqType )
        pass

    def _get_sequence_table( self, seqType ):
        """
        Generates the table to fill sequence data in order to create
        a :py:class:`.SequenceFrame`
        :param str seqType: Type of sequence: protein, dna, rna.
        :return: dict
        :raise ValueError: If seqType is not known
        """
        table = {}
        extra = []
        if seqType.lower() == "protein":
            table = {
                'C' : [], 'D' : [], 'S' : [], 'Q' : [], 'K' : [],
                'I' : [], 'P' : [], 'T' : [], 'F' : [], 'N' : [],
                'G' : [], 'H' : [], 'L' : [], 'R' : [], 'W' : [],
                'A' : [], 'V' : [], 'E' : [], 'Y' : [], 'M' : [],
                'X' : [], '*' : [], 'B' : [], 'Z' : []
            }
            extra = ['X', '*', 'B', 'Z']
        elif seqType.lower() in ["dna", "rna"]:
            # B = C or G or T  # D = A or G or T
            # H = A or C or T  # K = G or T
            # M = A or C       # N = A or C or G or T
            # R = A or G       # S = C or G
            # V = A or C or G  # W = A or T
            # Y = C or T
            table = {
                'C' : [], 'A' : [], 'T' : [], 'G' : [], 'X' : [], '*' : [],
                'B' : [], 'D' : [], 'H' : [], 'K' : [], 'M' : [], 'N' : [],
                'R' : [], 'S' : [], 'V' : [], 'W' : [], 'Y' : []
            }
            if seqType.lower() == "rna":
                table.setdefault( 'U' , [])
                table.pop('T', None)
            extra = ['X', '*', 'B', 'D', 'H', 'K', 'M', 'N', 'R', 'S', 'V', 'W', 'Y']
        else:
            raise ValueError("sequence type {0} unknown".format(seqType))
        return table, extra

    @property
    def _constructor(self):
        return DesignFrame
