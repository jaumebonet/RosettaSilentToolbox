# @Author: Jaume Bonet <bonet>
# @Date:   22-Feb-2018
# @Email:  jaume.bonet@gmail.com
# @Filename: reference.py
# @Last modified by:   bonet
# @Last modified time: 22-Feb-2018

import copy
import warnings

import pandas as pd


def has_reference_sequence( self, seqID ):
    """
    Checks if there is a reference sequence for the provided
    sequence ID.

    :param seqID: Identifier of the reference sequence
    :type seqID: :py:class:`str`

    :return: bool

    :raises:
        :TypeError: If applied over a data container without `_reference` attribute.
    """
    try:
        return seqID in self._reference and self._reference[seqID]["seq"] != ""
    except AttributeError:
        raise TypeError("The query object does not store reference data.")

def add_reference_sequence( self, seqID, sequence ):
    """
    Add a reference sequence attached to a particular sequence ID.

    :param seqID: Identifier of the reference sequence
    :type seqID: :py:class:`str`
    :param sequence: Reference sequence.
    :type sequence: :py:class:`str`

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame` or :py:class:`~pandas.Series`
        :KeyError: If the data container does not contain sequence data for the given seqID.
        :IndexError: If a reference structure exists and sequence length do not match
    """
    if not isinstance(self, (pd.DataFrame, pd.Series)):
        raise TypeError("Data container has to be a DataFrame/Series or a derived class.")
    if not seqID in self.get_available_sequences():
        raise KeyError("Data container does not have data for sequence {}".format(seqID))

    if seqID in self._reference:
        str = len(self._reference[seqID]["str"])
        if str > 0 and len(sequence) != str:
            raise IndexError("Structure length do not match sequence length")
        self._reference[seqID]["seq"] = sequence
    else:
        self._reference.setdefault(seqID, {"seq": sequence, "str": "", "stf": 1})

def get_reference_sequence( self, seqID ):
    """
    Get a reference sequence attached to a particular sequence ID.

    :param seqID: Identifier of the reference sequence
    :type seqID: :py:class:`str`

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame` or :py:class:`~pandas.Series`
        :KeyError: If there is no reference sequence for seqID.
    """
    if not isinstance(self, (pd.DataFrame, pd.Series)):
        raise TypeError("Data container has to be a DataFrame/Series or a derived class.")
    if not self.has_reference_sequence(seqID):
        raise KeyError("Data container does not have reference data for sequence {}".format(seqID))

    return self._reference[seqID]["seq"]

def has_reference_structure( self, seqID ):
    """
    Checks if there is a reference structure for the provided
    sequence ID.

    :param seqID: Identifier of the reference structure
    :type seqID: :py:class:`str`

    :return: bool

    :raises:
        :TypeError: If applied over a data container without `_reference` attribute.
    """
    try:
        return seqID in self._reference and self._reference[seqID]["str"] != ""
    except AttributeError:
        raise TypeError("The query object does not store reference data.")

def add_reference_structure( self, seqID, structure ):
    """
    Add a reference structure attached to a particular sequence ID.

    :param seqID: Identifier of the reference structure
    :type seqID: :py:class:`str`
    :param structure: Reference structure.
    :type structure: :py:class:`str`

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame` or :py:class:`~pandas.Series`
        :KeyError: If the data container does not contain structure data for the given seqID.
        :IndexError: If a reference sequence exists and structure length do not match
    """
    if not isinstance(self, (pd.DataFrame, pd.Series)):
        raise TypeError("Data container has to be a DataFrame/Series or a derived class.")
    if not seqID in self.get_available_structures():
        raise KeyError("Data container does not have data for structure {}".format(seqID))

    if seqID in self._reference:
        seq = len(self._reference[seqID]["seq"])
        if seq > 0 and len(structure) != seq:
            raise IndexError("Structure length do not match sequence length")
        self._reference[seqID]["str"] = structure
    else:
        self._reference.setdefault(seqID, {"str": structure, "seq": "", "stf": 1})

def get_reference_structure( self, seqID ):
    """
    Get a reference structure attached to a particular sequence ID.

    :param seqID: Identifier of the reference structure
    :type seqID: :py:class:`str`

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame` or :py:class:`~pandas.Series`
        :KeyError: If there is no reference structure for seqID.
    """
    if not isinstance(self, (pd.DataFrame, pd.Series)):
        raise TypeError("Data container has to be a DataFrame/Series or a derived class.")
    if not self.has_reference_structure(seqID):
        raise KeyError("Data container does not have reference data for structure {}".format(seqID))

    return self._reference[seqID]["str"]

def add_reference_shift( self, seqID, shift ):
    """
    Add a reference shift attached to a particular sequence ID.

    **What is shift?** In case the sequence does not start in 1,
    shift defines how to count it. It is a way to keep plotting
    and analysis showing residue number compatible with the PDB.
    There are two main ways to set the shift:

    #. Provide the number of the first residue of the chain; the \
    rest will be set up from there.
    #. If the original PDB has breaks, one will need to provide an \
    array with the numbers of each position.

    :param seqID: Identifier of the reference sequence
    :type seqID: :py:class:`str`
    :param shift: Starting residue number or per-residue number assignment.
    :type shift: Union[:py:class:`int`, :py:class:`list`]

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame` or :py:class:`~pandas.Series`
        :KeyError: If there is no reference structure or sequence for seqID.
        :KeyError: If shift is a list and the data container does not contain structure or sequence data for
        the given seqID.
        :IndexError: If shift is a list and the length is different than the reference sequence/structure
    """
    if not isinstance(self, (pd.DataFrame, pd.Series)):
        raise TypeError("Data container has to be a DataFrame/Series or a derived class.")
    if not seqID in self.get_available_structures() and not seqID in self.get_available_sequences():
        raise KeyError("Data container does not have data for structure {}".format(seqID))

    if isinstance(shift, list):
        if not self.has_reference_structure(seqID) and not self.has_reference_sequence(seqID):
            raise KeyError("Data container does not have reference data for structure {}".format(seqID))
        current_length = max(len(self.get_reference_sequence(seqID)), len(self.get_reference_structure(seqID)))
        if len(shift) != current_length:
            raise IndexError("Number of positions do not match reference sequence/structure length")
        self._reference[seqID]["sft"] = shift
    if isinstance(shift, int):
        if seqID in self._reference:
            self._reference[seqID]["sft"] = shift
        else:
            self._reference.setdefault(seqID, {"str": "", "seq": "", "stf": shift})

def get_reference_shift( self, seqID ):
    """
    Get a reference shift attached to a particular sequence ID.
    If none was provided will return 1 as default.

    :param seqID: Identifier of the reference sequence
    :type seqID: :py:class:`str`

    :type shift: Union[:py:class:`int`, :py:class:`list`]

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame` or :py:class:`~pandas.Series`
        :KeyError: If there is no reference structure or sequence for seqID.
    """
    if not isinstance(self, (pd.DataFrame, pd.Series)):
        raise TypeError("Data container has to be a DataFrame/Series or a derived class.")
    if not seqID in self.get_available_structures() and not seqID in self.get_available_sequences():
        raise KeyError("Data container does not have data for structure {}".format(seqID))

    if seqID in self._reference:
        return self._reference[seqID]["sft"]
    else:
        return 1

def add_reference( self, seqID, sequence="", structure="", shift=1):
    """
    Single access to :py:func:`.add_reference_sequence`, :py:func:`.add_reference_structure`
    and :py:func:`.add_reference_shift`.

    :param seqID: Identifier of the reference sequence
    :type seqID: :py:class:`str`
    :param sequence: Reference sequence.
    :type sequence: :py:class:`str`
    :param structure: Reference structure.
    :type structure: :py:class:`str`
    :param shift: Starting residue number or per-residue number assignment.
    :type shift: Union[:py:class:`int`, :py:class:`list`]

    :raises:
        :AttributeError: If trying to add a reference type to a class without it
    """
    if len(sequence) > 0:
        try:
            self.add_reference_sequence( seqID, sequence )
        except AttributeError:
            warnings.warn( "Trying to assign reference sequence to an object without the propery" )
    if len(structure) > 0:
        try:
            self.add_reference_structure( seqID, structure )
        except AttributeError:
            warnings.warn( "Trying to assign reference structure to an object without the propery" )
    try:
        self.add_reference_shift( seqID, shift )
    except AttributeError:
            warnings.warn( "Trying to assign reference shift to an object without the propery" )

def transfer_reference( self, df ):
    """
    Transfer reference from one data container to another. This **overwrittes** previous references completely.
    Use with care.

    :param df: Data container. Derives from :py:class:`~pandas.DataFrame` or :py:class:`~pandas.Series`.
    :type df: :py:class:`str`

    :raises:
        :AttributeError: If trying to add a reference type to/from a class without it
    """
    try:
        self._reference = copy.deepcopy(df._reference)
    except AttributeError:
        warnings.warn( "Either the target or source object cannot hold a reference." )

