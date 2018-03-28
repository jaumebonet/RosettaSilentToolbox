# @Author: Jaume Bonet <bonet>
# @Date:   20-Feb-2018
# @Email:  jaume.bonet@gmail.com
# @Filename: structure.py
# @Last modified by:   bonet
# @Last modified time: 23-Mar-2018


import collections

import pandas as pd


def positional_structural_count( df, seqID=None ):
    """
    Count the number of `H`, `E` in the structural provided data

    :param df: Input data.
    :type df: Union[:py:class:`.DesignFrame`, :py:class:`.FragmentFrame`]
    :param seqID: Identifier of the sequence sets of interest. Required when input is :py:class:`.DesignFrame`
    :type seqID: :py:class:`str`

    :return: :py:class:`~pandas.DataFrame` where rows are positions and columns are `H`, `E` and `L`.

    :raises:
        :AttributeError: if the data passed is not in Union[:py:class:`.DesignFrame`, :py:class:`.FragmentFrame`].
        :AttributeError: if input is :py:class:`.DesignFrame` and ``seqID`` is not provided.
        :KeyError: if input is :py:class:`.DesignFrame` and ``seqID`` cannot be found.
    """
    from rstoolbox.components import DesignFrame, FragmentFrame
    data = {"H": [], "E": [], "L": []}

    # @todo: Code positional_structural_count for DesignFrame
    # @body: Should behave the same way it does for the FragmentFrame
    if isinstance(df, DesignFrame):
        if seqID is None:
            raise AttributeError("seqID needs to be provided")
        if not "structure_{}".format(seqID) in df:
            raise KeyError("Structure {} not found in decoys.".format(seqID))

    elif isinstance(df, FragmentFrame):
        for i in df["position"].drop_duplicates().values:
            qseq = "".join(df[df["position"] == i]["sse"].values).upper()
            sse = collections.Counter(qseq)
            data["H"].append(float(sse["H"]) / float(len(qseq)))
            data["E"].append(float(sse["E"]) / float(len(qseq)))
            data["L"].append(float(sse["L"]) / float(len(qseq)))

    else:
        raise AttributeError("Input data has to be a DesignFrame with a reference sequence or a FragmentFrame.")

    return pd.DataFrame(data)


def positional_structural_identity( df, seqID=None, ref_sse=None ):
    """
    Evaluate per position how many times the provided data matches the expected secondary structure.

    :param df: Input data.
    :type df: Union[:py:class:`.DesignFrame`, :py:class:`.FragmentFrame`]
    :param seqID: Identifier of the sequence sets of interest. Required when input is :py:class:`.DesignFrame`
    :type seqID: :py:class:`str`
    :param ref_sse: Reference sequence. Required.
    :type ref_sse: :py:class:`str`

    :return: :py:class:`~pandas.DataFrame` where rows are positions and columns are `sse`, `max_sse` and `identity_perc`.

    :raises:
        :AttributeError: if the data passed is not in Union[:py:class:`.DesignFrame`, :py:class:`.FragmentFrame`].
        :AttributeError: if input is :py:class:`.DesignFrame` and ``seqID`` is not provided.
        :KeyError: if input is :py:class:`.DesignFrame` and ``seqID`` cannot be found.
        :AttributeError if ``ref_sse`` is not provided.
    """
    from rstoolbox.components import DesignFrame, FragmentFrame
    data = {"sse": [], "max_sse": [], "identity_perc": []}

    # @todo: Code positional_structural_identity for DesignFrame
    # @body: Should behave the same way it does for the FragmentFrame
    if ref_sse is None:
        raise AttributeError("ref_sse needs to be provided")

    if isinstance(df, DesignFrame):
        if seqID is None:
            raise AttributeError("seqID needs to be provided")
        if not "structure_{}".format(seqID) in df:
            raise KeyError("Structure {} not found in decoys.".format(seqID))

    elif isinstance(df, FragmentFrame):
        for i in df["position"].drop_duplicates().values:
            qseq = "".join(df[df["position"] == i]["sse"].values).upper()
            sse = collections.Counter(qseq)
            data["sse"].append(ref_sse[i - 1])
            data["max_sse"].append(sse.most_common(1)[0][0])
            data["identity_perc"].append(float(sse[ref_sse[i - 1]]) / float(len(qseq)))

    else:
        raise AttributeError("Input data has to be a DesignFrame with a reference sequence or a FragmentFrame.")

    return pd.DataFrame(data)
