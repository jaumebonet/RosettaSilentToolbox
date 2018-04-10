# @Author: Jaume Bonet <bonet>
# @Date:   22-Feb-2018
# @Email:  jaume.bonet@gmail.com
# @Filename: reference.py
# @Last modified by:   bonet
# @Last modified time: 10-Apr-2018

import copy
import warnings

import six
import pandas as pd
import numpy as np


def _has_reference( obj, ctype, seqID ):
    try:
        return seqID in obj._reference and obj._reference[seqID][ctype[:3]] != ""
    except AttributeError:
        raise TypeError("The query object does not store reference data.")


def _get_reference( obj, ctype, seqID ):
    if not isinstance(obj, (pd.DataFrame, pd.Series)):
        raise TypeError("Data container has to be a DataFrame/Series or a derived class.")
    if not _has_reference(obj, ctype, seqID):
        raise KeyError("No reference found for {0} {1}".format(ctype, seqID))
    return obj._reference[seqID][ctype[:3]]


def _get_key_reference( obj, ctype, seqID, key_residues ):
    from rstoolbox.components import get_selection

    seq = _get_reference(obj, ctype, seqID)
    sft = _get_reference(obj, "sft", seqID)
    if key_residues is None:
        return seq

    kr = get_selection(key_residues, seqID, sft, len(seq))

    return "".join(np.array(list(seq))[kr - 1])


def get_available_references( self ):
    """
    Identify for which seqIDs there is some reference data.
    """
    return self._reference.keys()


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
    return _has_reference(self, "sequence", seqID)


def add_reference_sequence( self, seqID, sequence ):
    """
    Add a reference sequence attached to a particular sequence ID.

    :param seqID: Identifier of the reference sequence
    :type seqID: :py:class:`str`
    :param sequence: Reference sequence.
    :type sequence: :py:class:`str`

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame`
            or :py:class:`~pandas.Series`
        :ValueError: If ``sequence`` is not a :class:`str`
        :KeyError: If the data container does not contain sequence data for the given seqID.
        :IndexError: If a reference structure exists and sequence length do not match
    """
    if not isinstance(self, (pd.DataFrame, pd.Series)):
        raise TypeError("Data container has to be a DataFrame/Series or a derived class.")
    if self._subtyp != "sequence_frame" and seqID not in self.get_available_sequences():
        raise KeyError("Data container does not have data for sequence {}".format(seqID))
    if not isinstance(sequence, six.string_types):
        raise ValueError("Reference sequence must be a string.")

    if seqID in self._reference:
        seq = len(self._reference[seqID]["str"])
        if seq > 0 and len(sequence) != seq:
            raise IndexError("Structure length do not match sequence length")
        self._reference[seqID]["seq"] = sequence
    else:
        self._reference.setdefault(seqID, {"seq": sequence, "str": "", "sft": 1})


def get_reference_sequence( self, seqID, key_residues=None ):
    """
    Get a reference sequence attached to a particular sequence ID.

    :param seqID: Identifier of the reference sequence
    :type seqID: :py:class:`str`
    :param key_residues: Residues of interest to pick
    :type key_residues: Union[:py:class:`int`,
        :py:class:`list`(:py:class:`int`), :py:class:`.Selection` ]

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame`
        or :py:class:`~pandas.Series`
        :KeyError: If there is no reference sequence for seqID.
    """
    if key_residues is None:
        return _get_reference(self, "sequence", seqID)
    else:
        return _get_key_reference(self, "sequence", seqID, key_residues)


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
    return _has_reference(self, "structure", seqID)


def add_reference_structure( self, seqID, structure ):
    """
    Add a reference structure attached to a particular sequence ID.

    :param seqID: Identifier of the reference structure
    :type seqID: :py:class:`str`
    :param structure: Reference structure.
    :type structure: :py:class:`str`

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame`
            or :py:class:`~pandas.Series`
        :ValueError: If ``sequence`` is not a :class:`str`
        :KeyError: If the data container does not contain structure data for the given seqID.
        :IndexError: If a reference sequence exists and structure length do not match
    """
    if not isinstance(self, (pd.DataFrame, pd.Series)):
        raise TypeError("Data container has to be a DataFrame/Series or a derived class.")
    if self._subtyp != "sequence_frame" and seqID not in self.get_available_structures():
        raise KeyError("Data container does not have data for structure {}".format(seqID))
    if not isinstance(structure, six.string_types):
        raise ValueError("Reference structure must be a string.")

    if seqID in self._reference:
        seq = len(self._reference[seqID]["seq"])
        if seq > 0 and len(structure) != seq:
            raise IndexError("Structure length do not match sequence length")
        self._reference[seqID]["str"] = structure
    else:
        self._reference.setdefault(seqID, {"str": structure, "seq": "", "sft": 1})


def get_reference_structure( self, seqID, key_residues=None ):
    """
    Get a reference structure attached to a particular sequence ID.

    :param seqID: Identifier of the reference structure
    :type seqID: :py:class:`str`
    :param key_residues: Residues of interest to pick
    :type key_residues: Union[:py:class:`int`,
        :py:class:`list`(:py:class:`int`), :py:class:`.Selection` ]

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame`
        or :py:class:`~pandas.Series`
        :KeyError: If there is no reference structure for seqID.
    """
    if key_residues is None:
        return _get_reference( self, "structure", seqID )
    else:
        return _get_key_reference(self, "structure", seqID, key_residues)


def add_reference_shift( self, seqID, shift, shift_labels=False ):
    """
    Add a reference shift attached to a particular sequence ID.

    **What is shift?** In case the sequence does not start in 1,
    shift defines how to count it. It is a way to keep plotting
    an analysis showing residue number compatible with the PDB.
    There are two main ways to set the shift:

    #. Provide the number of the first residue of the chain; the \
    rest will be set up from there.
    #. If the original PDB has breaks, one will need to provide an \
    array with the numbers of each position.

    :param seqID: Identifier of the reference sequence
    :type seqID: :py:class:`str`
    :param shift: Starting residue number or per-residue number assignment.
    :type shift: Union[:py:class:`int`, :py:class:`list`]
    :type shift: :py:class:`bool`
    :param shift_labels: When adding the shift, apply it to any label
    :type shift_labels: :py:class:`bool`

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame`
            or :py:class:`~pandas.Series`
        :KeyError: If there is no reference structure or sequence for seqID.
        :KeyError: If shift is a list and the data container does not contain structure
            or sequence data for the given seqID.
        :IndexError: If shift is a list and the length is different than the reference
            sequence/structure
    """
    if not isinstance(self, (pd.DataFrame, pd.Series)):
        raise TypeError("Data container has to be a DataFrame/Series or a derived class.")
    if self._subtyp != "sequence_frame" and (seqID not in self.get_available_structures() and
                                             seqID not in self.get_available_sequences()):
        raise KeyError("No reference found for sequece/structure {}".format(seqID))

    if isinstance(shift, list):
        if not self.has_reference_structure(seqID) and not self.has_reference_sequence(seqID):
            raise KeyError("No reference data for sequence/structure {}".format(seqID))
        current_length = 0
        if self.has_reference_sequence(seqID):
            current_length = len(self.get_reference_sequence(seqID))
        elif self.has_reference_structure(seqID):
            current_length = len(self.get_reference_structure(seqID))
        if len(shift) != current_length:
            raise IndexError("Number of positions do not match reference sequence/structure length")
        self._reference[seqID]["sft"] = shift
    if isinstance(shift, int):
        if seqID in self._reference:
            self._reference[seqID]["sft"] = shift
        else:
            self._reference.setdefault(seqID, {"str": "", "seq": "", "sft": shift})

    if shift_labels:
        labels = self.get_available_labels()
        for lbl in labels:
            clnm = "lbl_{}".format(lbl)
            if isinstance(self, pd.DataFrame):
                self.apply(lambda x: x[clnm].shift(seqID, shift), axis=1)
            else:
                self[clnm].shift(seqID, shift)


def get_reference_shift( self, seqID ):
    """
    Get a reference shift attached to a particular sequence ID.
    If none was provided will return 1 as default.

    :param seqID: Identifier of the reference sequence
    :type seqID: :py:class:`str`

    :type shift: Union[:py:class:`int`, :py:class:`list`]

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame`
        or :py:class:`~pandas.Series`
        :KeyError: If there is no reference structure or sequence for seqID.
    """
    if not isinstance(self, (pd.DataFrame, pd.Series)):
        raise TypeError("Data container has to be a DataFrame/Series or a derived class.")
    if self._subtyp != "sequence_frame" and (seqID not in self.get_available_structures() and
                                             seqID not in self.get_available_sequences()):
        raise KeyError("Data container does not have data for structure {}".format(seqID))

    if seqID in self._reference:
        return self._reference[seqID]["sft"]
    else:
        return 1


def add_reference( self, seqID, sequence="", structure="", shift=1, shift_labels=False ):
    """
    Single access to :py:meth:`.add_reference_sequence`, :py:meth:`.add_reference_structure`
    and :py:meth:`.add_reference_shift`.

    :param seqID: Identifier of the reference sequence
    :type seqID: :py:class:`str`
    :param sequence: Reference sequence.
    :type sequence: :py:class:`str`
    :param structure: Reference structure.
    :type structure: :py:class:`str`
    :param shift: Starting residue number or per-residue number assignment.
    :type shift: Union[:py:class:`int`, :py:class:`list`]
    :param shift_labels: When adding the shift, apply it to any label
    :type shift_labels: :py:class:`bool`

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
        self.add_reference_shift( seqID, shift, shift_labels )
    except AttributeError:
            warnings.warn( "Trying to assign reference shift to an object without the propery" )


def transfer_reference( self, df ):
    """
    Transfer reference from one data container to another. This **overwrittes** previous
    references completely. Use with care.

    :param df: Data container. Derives from :py:class:`~pandas.DataFrame`
    or :py:class:`~pandas.Series`.
    :type df: :py:class:`str`

    :raises:
        :AttributeError: If trying to add a reference type to/from a class without it
    """
    try:
        self._reference = copy.deepcopy(df._reference)
    except AttributeError:
        warnings.warn( "Either the target or source object cannot hold a reference." )


def delete_reference( self, seqID, shift_labels=False ):
    """
    Remove all reference data regarding a particular seqID.

    :param seqID: Identifier of the reference sequence
    :type seqID: :py:class:`str`
    :param shift_labels: When removing the shift, apply it to any label that has been shifted
    :type shift_labels: :py:class:`bool`

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame`
        or :py:class:`~pandas.Series`
    """
    if not isinstance(self, (pd.DataFrame, pd.Series)):
        raise TypeError("Data container has to be a DataFrame/Series or a derived class.")

    if seqID not in self._reference:
        return

    shift = self.get_reference_shift(seqID)
    del(self._reference[seqID])

    if shift_labels:
        labels = self.get_available_labels()
        for lbl in labels:
            clnm = "lbl_{}".format(lbl)
            if isinstance(self, pd.DataFrame):
                self.apply(lambda x: x[clnm].unshift(seqID, shift), axis=1)
            else:
                self[clnm].unshift(seqID, shift)
