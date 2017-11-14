# @Author: Jaume Bonet <bonet>
# @Date:   29-Jun-2017
# @Email:  jaume.bonet@gmail.com
# @Filename: design.py
# @Last modified by:   bonet
# @Last modified time: 10-Nov-2017

# Standard Libraries
import os
import re

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

    def sequence_frequencies( self, seqID, seqType="protein" ):
        """
        Generates a :py:class:`.SequenceFrame` for the frequencies of
        the sequences in the __DesignFrame__ with seqID identifier.
        If there is a reference_sequence for this seqID, it will also
        be attached to the __SequenceFrame__.
        :param str seqID: Identifier of the sequence sets of interest.
        :param str seqType: Type of sequence: protein, dna, rna.
        :return: :py:class:`.SequenceFrame`
        """
        sserie = self["sequence_{0}".format(seqID)].values
        table = self._get_sequence_table( seqType )
        for x in range(len(sserie[0])):
            for k in table:
                table[k].append(float(0))
            for y in range(len(sserie)):
                aa = sserie[y][x]
                table[aa][-1] += float(1)
        for k in table:
            for x in range(len(table[k])):
                if table[k][x] != 0:
                    table[k][x] /= float(len(sserie))

        df = SequenceFrame(table)
        if self.has_reference_sequence(seqID):
            df.reference_sequence(self.reference_sequence(seqID))
        df.index = df.index + 1
        return df

    def sequence_bits( self, seqID, seqType="protein" ):
        """
        Generates a :py:class:`.SequenceFrame` for the bits of
        the sequences in the __DesignFrame__ with seqID identifier.
        If there is a reference_sequence for this seqID, it will also
        be attached to the __SequenceFrame__.
        :param str seqID: Identifier of the sequence sets of interest.
        :param str seqType: Type of sequence: protein, dna, rna.
        :return: :py:class:`.SequenceFrame`
        """
        sserie = self["sequence_{0}".format(seqID)].values
        table = self._get_sequence_table( seqType )
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
        if seqType.lower() == "protein":
            table = {
                'C' : [], 'D' : [], 'S' : [], 'Q' : [], 'K' : [],
                'I' : [], 'P' : [], 'T' : [], 'F' : [], 'N' : [],
                'G' : [], 'H' : [], 'L' : [], 'R' : [], 'W' : [],
                'A' : [], 'V' : [], 'E' : [], 'Y' : [], 'M' : []
            }
        elif seqType.lower() == "dna":
            table = {
                'C' : [], 'A' : [], 'T' : [], 'G' : []
            }
        elif seqType.lower() == "rna":
            table = {
                'C' : [], 'A' : [], 'U' : [], 'G' : []
            }
        else:
            raise ValueError("sequence type {0} unknown".format(seqType))
        return table

    @property
    def _constructor(self):
        return DesignFrame
