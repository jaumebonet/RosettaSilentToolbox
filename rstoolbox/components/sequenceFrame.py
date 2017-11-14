import os
import re
import pandas as pd

class SequenceFrame( pd.DataFrame ):
    """
    The :py:class:`.SequenceFrame` extends :py:class:`pandas.DataFrame`.
    Ideally, in this DataFrame, each column represents a ResidueType (or
    NucleotideType) and each row a sequence position.
    A reference sequence can be provided, but it has to be the exact same
    length as the sequence expressed in the DataFrame.
    As it represents sequence, unless the SequenceFrame is generated by
    slicing another, one should expect the index to start in 1.
    """
    _metadata = ['_reference_sequence', '_measure']
    _reference_sequence = ""
    _measure = ""

    def reference_sequence( self, sequence=None ):
        """
        Setter/Getter for a reference sequenceD.
        :param str sequence: Reference sequence. By default is
            :py:data:`None`, which turns the function into a getter.
        """
        if sequence is not None:
            assert(len(sequence) == self.shape[0])
            self._reference_sequence = sequence
        return self._reference_sequence

    def measure( self, measure=None ):
        known = ["frequency", "bits"]
        if measure is not None:
            if measure.lower() not in known:
                raise ValueError("Measure type {0} not in {1}".format(measure, ",".join(known)))
            self._measure = measure
        return self._measure


    def has_reference_sequence( self ):
        """
        Checks if there is a reference sequence.
        :param str seqID: Identifier of the reference sequence
        :return: bool
        """
        return self._reference_sequence != ""

    def select_key_sequence_positions( self, selector ):
        """
        By providing a list of key positions to keep, the
        SequenceFrame is filtered so that only those positions
        are kept
        :param list selector: List of integers of the positions
            of interest.
        :return: :py:class:`.SequenceFrame`
        """
        pass

    @property
    def _constructor( self ):
        return SequenceFrame
