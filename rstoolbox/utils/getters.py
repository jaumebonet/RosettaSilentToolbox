# @Author: Jaume Bonet <bonet>
# @Date:   22-Feb-2018
# @Email:  jaume.bonet@gmail.com
# @Filename: getters.py
# @Last modified by:   bonet
# @Last modified time: 22-Feb-2018


import pandas as pd

def get_available_sequences( self ):
    """
    List which sequence identifiers are available in the data container

    :return: :py:class:`list`

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame` or :py:class:`~pandas.Series`
    """
    if not isinstance(self, (pd.DataFrame, pd.Series)):
        raise TypeError("Data container has to be a DataFrame/Series or a derived class.")
    return [x.split("_")[-1] for x in self.columns.values if x.startswith("sequence_")]

def get_sequence( self, seqID ):
    """
    Return the sequence(s) data.

    :param seqID: Identifier of the sequence of interest.
    :type seqID: :py:class:`str`

    :return: py:class:`str` or :py:class:`~pandas.Series`

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame` or :py:class:`~pandas.Series`
        :KeyError: If the column `sequence_[seqID] is cannot be found.
    """
    col = "sequence_{}".format(seqID)
    if not isinstance(self, (pd.DataFrame, pd.Series)):
        raise TypeError("Data container has to be a DataFrame/Series or a derived class.")
    if not col in self:
        raise KeyError("Sequence for {0} not found in data set.  Column `sequence_{0}` is missing.".format(seqID))
    return self[col]

def get_available_structures( self ):
    """
    List which structure identifiers are available in the data container

    :return: :py:class:`list`

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame` or :py:class:`~pandas.Series`
    """
    if not isinstance(self, (pd.DataFrame, pd.Series)):
        raise TypeError("Data container has to be a DataFrame/Series or a derived class.")
    return [x.split("_")[-1] for x in self.columns.values if x.startswith("structure_")]

def get_structure( self, seqID ):
    """
    Return the structure(s) data.

    :param seqID: Identifier of the structure of interest.
    :type seqID: :py:class:`str`

    :return: py:class:`str` or :py:class:`~pandas.Series`

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame` or :py:class:`~pandas.Series`
        :KeyError: If the column `structure_[seqID] is cannot be found.
    """
    col = "structure_{}".format(seqID)
    if not isinstance(self, (pd.DataFrame, pd.Series)):
        raise TypeError("Data container has to be a DataFrame/Series or a derived class.")
    if not col in self:
        raise KeyError("Structure for {0} not found in data set. Column `structure_{0}` is missing.".format(seqID))
    return self[col]

def get_available_structure_predictions( self ):
    """
    List which structure prediction identifiers are available in the data container

    :return: :py:class:`list`

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame` or :py:class:`~pandas.Series`
    """
    if not isinstance(self, (pd.DataFrame, pd.Series)):
        raise TypeError("Data container has to be a DataFrame/Series or a derived class.")
    return [x.split("_")[-1] for x in self.columns.values if x.startswith("psipred_")]

def get_structure_prediction( self, seqID ):
    """
    Return the structure prediction(s) data.

    :param seqID: Identifier of the structure of interest.
    :type seqID: :py:class:`str`

    :return: py:class:`str` or :py:class:`~pandas.Series`

    :raises:
        :TypeError: If the data container is not :py:class:`~pandas.DataFrame` or :py:class:`~pandas.Series`
        :KeyError: If the column `psipred_[seqID] is cannot be found.
    """
    col = "psipred_{}".format(seqID)
    if not isinstance(self, (pd.DataFrame, pd.Series)):
        raise TypeError("Data container has to be a DataFrame/Series or a derived class.")
    if not col in self:
        raise KeyError("Structure prediction for {0} not found in data set. Column `psipred_{0}` is missing.".format(seqID))
    return self[col]
