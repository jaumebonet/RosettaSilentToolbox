# @Author: Jaume Bonet <bonet>
# @Date:   22-Feb-2018
# @Email:  jaume.bonet@gmail.com
# @Filename: getters.py
# @Last modified by:   bonet
# @Last modified time: 27-Mar-2018


import pandas as pd
import numpy as np


def _check_type( obj ):
    if not isinstance(obj, (pd.DataFrame, pd.Series)):
        raise TypeError("Data container has to be a DataFrame/Series or a derived class.")


def _check_column( obj, ctype, seqID ):
    if "{0}_{1}".format(ctype, seqID) not in obj:
        raise KeyError(
            "{0} for {1} not found in data set. ".format(ctype.capitalize(), seqID) +
            "Column `{0}_{1}` is missing.".format(ctype, seqID))
    return "{0}_{1}".format(ctype, seqID)


def _get_available( obj, ctype ):
    if isinstance(obj, pd.DataFrame):
        return ["_".join(x.split("_")[1:]) for x in obj.columns.values if x.startswith(ctype)]
    else:
        return ["_".join(x.split("_")[1:]) for x in obj.index.values if x.startswith(ctype)]


def _get_key_sequence( obj, ctype, seqID, key_residues ):
    from rstoolbox.components import get_selection
    from .reference import _get_reference

    seq = obj[_check_column(obj, ctype, seqID)]
    sft = _get_reference(obj, "sft", seqID)

    if isinstance(obj, pd.Series):
        length = len(seq)
    else:
        length = len(seq.iloc[0])
    kr = get_selection(key_residues, seqID, sft, length)

    if isinstance(obj, pd.Series):
        return "".join(np.array(list(seq))[kr - 1])
    else:
        return seq.apply(lambda seq: "".join(np.array(list(seq))[kr - 1]))


def get_id( self ):
    """
    Return identifier data for the design(s).

    :return: :py:class:`str` or :py:class:`~pandas.Series`

    :raises:
        :TypeError: If the data container is not `~pandas.DataFrame`
        or :py:class:`~pandas.Series`
        :KeyError: If the column `sequence_[seqID]` is cannot be found.
    """
    _check_type(self)
    if "description" not in self:
        raise KeyError(
            "Identifiers not found in data set. "
            "Column `description` is missing.")
    return self["description"]


def get_available_sequences( self ):
    """
    List which sequence identifiers are available in the data container

    :return: :py:class:`list`

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame`
        or :py:class:`~pandas.Series`
    """
    _check_type(self)
    return _get_available(self, "sequence_")


def get_sequence( self, seqID, key_residues=None ):
    """
    Return the sequence data for `seqID` available in the container.

    :param seqID: Identifier of the sequence of interest.
    :type seqID: :py:class:`str`
    :param key_residues: Residues of interest. Are affected by the reference shift (if any).
    :type key_residues: Union[:py:class:`list`[:py:class:`int`], :py:class:`str`]

    :return: :py:class:`str` or :py:class:`~pandas.Series`

    :raises:
        :TypeError: If the data container is not `~pandas.DataFrame`
        or :py:class:`~pandas.Series`
        :KeyError: If the column `sequence_[seqID]` is cannot be found.
    """
    _check_type(self)
    if key_residues is None:
        return self[_check_column(self, "sequence", seqID)]
    else:
        return _get_key_sequence(self, "sequence", seqID, key_residues)


def get_available_structures( self ):
    """
    List which structure identifiers are available in the data container

    :return: :py:class:`list`

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame`
        or :py:class:`~pandas.Series`
    """
    _check_type(self)
    return _get_available(self, "structure_")


def get_structure( self, seqID, key_residues=None ):
    """
    Return the structure(s) data.

    :param seqID: Identifier of the structure of interest.
    :type seqID: :py:class:`str`
    :param key_residues: Residues of interest. Are affected by the reference shift (if any).
    :type key_residues: Union[:py:class:`list`[:py:class:`int`], :py:class:`str`]

    :return: py:class:`str` or :py:class:`~pandas.Series`

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame`
        or :py:class:`~pandas.Series`
        :KeyError: If the column `structure_[seqID] is cannot be found.
    """
    _check_type(self)
    if key_residues is None:
        return self[_check_column(self, "structure", seqID)]
    else:
        return _get_key_sequence(self, "structure", seqID, key_residues)


def get_available_structure_predictions( self ):
    """
    List which structure prediction identifiers are available in the data container

    :return: :py:class:`list`

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame`
        or :py:class:`~pandas.Series`
    """
    _check_type(self)
    return _get_available(self, "psipred_")


def get_structure_prediction( self, seqID, key_residues=None ):
    """
    Return the structure prediction(s) data.

    :param seqID: Identifier of the structure of interest.
    :type seqID: :py:class:`str`
    :param key_residues: Residues of interest. Are affected by the reference shift (if any).
    :type key_residues: Union[:py:class:`list`[:py:class:`int`], :py:class:`str`]

    :return: py:class:`str` or :py:class:`~pandas.Series`

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame`
        or :py:class:`~pandas.Series`
        :KeyError: If the column `psipred_[seqID] is cannot be found.
    """
    _check_type(self)
    if key_residues is None:
        return self[_check_column(self, "psipred", seqID)]
    else:
        return _get_key_sequence(self, "psipred", seqID, key_residues)


def get_sequential_data( self, query, seqID ):
    """
    Provides data on the requested query.

    :param query: Query type: `sequence`, `structure`, `structure_prediction`.
    :type query: :py:class:`str`
    :param seqID: Identifier of the structure of interest.
    :type seqID: :py:class:`str`

    :return: py:class:`str` or :py:class:`~pandas.Series`

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame`
        or :py:class:`~pandas.Series`
        :KeyError: If `query` has a non-accepted value.
    """
    queries = ["sequence", "structure", "structure_prediction"]
    if query.lower() not in queries:
        raise KeyError("Available queries are: {}".format(",".join(queries)))

    if query.lower() == "sequence":
        return self.get_sequence(seqID)
    if query.lower() == "structure":
        return self.get_structure(seqID)
    if query.lower() == "structure_prediction":
        return self.get_structure_prediction(seqID)


def get_available_labels( self ):
    """
    List which slabels are available in the data container.

    :return: :py:class:`list`

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame`
        or :py:class:`~pandas.Series`
    """
    _check_type(self)
    return _get_available(self, "lbl_")


def get_label( self, label, seqID=None ):
    """
    Return the content(s) of the labels of interest. As a py:class:`.Selection`
    for a given sequece. If only one seqID is available, it will automatically
    pick labels for that, even if other data is present; otherwise, seqID must
    be provided. This takes into account availability of sequence, structure and
    psipred data.

    :param label: Label identifier. Will be uppercased.
    :type label: :py:class:`str`
    :param seqID: Identifier of the sequence of interest.
    :type seqID: :py:class:`str`

    :return: py:class:`.Selection`

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame`
        or :py:class:`~pandas.Series`
        :KeyError: If [seqID] is not specified and more than one seqID is possible.
        :KeyError: If the column `lbl_[label]` is cannot be found.
    """
    _check_type(self)
    if seqID is None:
        ss = set(get_available_sequences(self))
        sr = set(get_available_structures(self))
        sp = set(get_available_structure_predictions(self))
        sa = ss.union(sr).union(sp)
        if len(sa) > 1:
            raise KeyError("Information for multiple seqID is present. Choose one.")
        else:
            seqID = sa.pop()
    data = self[_check_column(self, "lbl", label.upper())]
    if isinstance(self, pd.DataFrame):
        return data.apply(lambda x: x[seqID])
    else:
        return data[seqID]
