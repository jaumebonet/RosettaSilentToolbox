# @Author: Jaume Bonet <bonet>
# @Date:   22-Feb-2018
# @Email:  jaume.bonet@gmail.com
# @Filename: getters.py
# @Last modified by:   bonet
# @Last modified time: 05-Mar-2018


import pandas as pd


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


def get_sequence( self, seqID ):
    """
    Return the sequence data for `seqID` available in the container.

    :param seqID: Identifier of the sequence of interest.
    :type seqID: :py:class:`str`

    :return: :py:class:`str` or :py:class:`~pandas.Series`

    :raises:
        :TypeError: If the data container is not `~pandas.DataFrame`
        or :py:class:`~pandas.Series`
        :KeyError: If the column `sequence_[seqID]` is cannot be found.
    """
    _check_type(self)
    return self[_check_column(self, "sequence", seqID)]


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


def get_structure( self, seqID ):
    """
    Return the structure(s) data.

    :param seqID: Identifier of the structure of interest.
    :type seqID: :py:class:`str`

    :return: py:class:`str` or :py:class:`~pandas.Series`

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame`
        or :py:class:`~pandas.Series`
        :KeyError: If the column `structure_[seqID] is cannot be found.
    """
    _check_type(self)
    return self[_check_column(self, "structure", seqID)]


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


def get_structure_prediction( self, seqID ):
    """
    Return the structure prediction(s) data.

    :param seqID: Identifier of the structure of interest.
    :type seqID: :py:class:`str`

    :return: py:class:`str` or :py:class:`~pandas.Series`

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame`
        or :py:class:`~pandas.Series`
        :KeyError: If the column `psipred_[seqID] is cannot be found.
    """
    _check_type(self)
    return self[_check_column(self, "psipred", seqID)]


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
