# @Author: Jaume Bonet <bonet>
# @Date:   10-Oct-2017
# @Email:  jaume.bonet@gmail.com
# @Filename: sequence.py
# @Last modified by:   bonet
# @Last modified time: 20-Feb-2018

import string
import os

import pandas as pd
import numpy as np

from rstoolbox.components import DesignFrame, FragmentFrame
from SimilarityMatrix import SimilarityMatrix as SM


def _extract_key_residue_sequence( seq, key_residues=None ):
    if key_residues is None:
        return seq
    tmp_seq = ""
    for k in key_residues:
        tmp_seq += seq[k-1]
    return tmp_seq

def _pick_key_residues( seq, key_residues=None ):
    if key_residues is None:
        return seq
    tmp_seq = ""
    for k in key_residues:
        tmp_seq += seq[k]
    return tmp_seq

def _calculate_linear_sequence_similarity( qseq, rseq, matrix, key_residues=None ):
    score = 0;
    qseq = _extract_key_residue_sequence( qseq, key_residues )
    rseq = _extract_key_residue_sequence( rseq, key_residues )
    assert len(qseq) == len(rseq)
    for i in range(len(qseq)):
        score += matrix.get_value(qseq[i], rseq[i])
    return score

def _sequence_similarity( qseq, rseq, matrix ):
    assert len(qseq) == len(rseq)
    raw, idn, pos, neg = 0, 0, 0, 0
    for i in range(len(qseq)):
        sc = matrix.get_value(qseq[i], rseq[i])
        raw += sc
        if qseq[i] == rseq[i]:
            idn += 1
        if sc > 0:
            pos += 1
        else:
            neg += 1
    return raw, idn, pos, neg

def _positional_similarity( qseq, rseq, matrix ):
    raw, idn, pos, neg = 0, 0, 0, 0
    for i in range(len(qseq)):
        sc = matrix.get_value(qseq[i], rseq)
        raw += sc
        if qseq[i] == rseq:
            idn += 1
        if sc > 0:
            pos += 1
        else:
            neg += 1
    return raw, idn, pos, neg

def _calculate_binary_sequence_similarity( qseq, rseq, matrix, key_residues=None ):
    new_seq = ""
    qseq = _extract_key_residue_sequence( qseq, key_residues )
    rseq = _extract_key_residue_sequence( rseq, key_residues )
    assert len(qseq) == len(rseq)
    for i in range(len(qseq)):
        if  matrix.get_value(qseq[i], rseq[i]) >= 1:
            new_seq += "1"
        else:
            new_seq += "0"
    return new_seq

def pick_key_residues( df, seqID, key_residues=None ):
    """
    Select the residues of interest from sequence seqID.

    :param df: designs data.
    :type df: :py:class:`.DesignFrame`
    :param seqID: Identifier of the sequence sets of interest.
    :type seqID: :py:class:`str`
    :param key_residues: Residues of interest. Are affected by the reference shift (if any).
    One can provide a list of residues to work with or a string indicating a LABEL to use.
    :type key_residues: Union[:py:class:`list`[:py:class:`int`], :py:class:`str`]

    :return: :py:class:`~pandas.DataFrame` with a single column named `sequence_[seqID]`.

    :raises:
        :AttributeError: if the data passed is not a :py:class:`.DesignFrame`.
        :KeyError: if ``seqID`` cannot be found.
        :KeyError: if a label is requested and not found .
    """

    if not isinstance(df, DesignFrame):
        raise AttributeError("Input data has to be a DesignFrame with a reference sequence.")
    if not "sequence_{}".format(seqID) in df:
        raise KeyError("Sequence {} not found in decoys.".format(seqID))

    if key_residues is None:
        return pd.DataFrame(df["sequence_{}".format(seqID)])

    selected_residues = np.array(key_residues) - df.reference_shift(seqID)

    if np.any(selected_residues < 0):
        raise IndexError("Selected residue out of specified bound")
    if np.any(selected_residues > len(df["sequence_{}".format(seqID)].values[0]) + df.reference_shift(seqID) - 1):
        raise IndexError("Selected residue out of specified bound")

    sr = df.apply(lambda x : _pick_key_residues(x["sequence_{}".format(seqID)], selected_residues), axis=1)
    return pd.DataFrame(sr).rename({0:"sequence_{}".format(seqID)}, axis="columns")


def sequence_similarity( df, seqID, key_residues=None, matrix="BLOSUM62" ):
    """
    Evaluate the sequence similarity between each decoy and the reference sequence for a given seqID.
    Will create several new columns:

    #. **<matrix>_<seqID>_score_raw:** Contains the score obtained by comparing the sequence to the reference.
    #. **<matrix>_<seqID>_score_perc:** Contains the score obtained by comparing the sequence to the reference.\
    in relation to the maximum possible score (reference sequence against itself).
    #. **<matrix>_<seqID>_identity:**  Count of the number of identity matches.
    #. **<matrix>_<seqID>_positive:**  Count of the number of positive matches; those that score positively in the \
    similarity matrix used. This will, by definition, include the identities too.
    #. **<matrix>_<seqID>_negative:** Count of the number of negative matches, representing unfavoured substitutions \
    according the to similarity matrix. Includes matches scored as 0. Logically, _positives_ + _negatives_ adds up \
    to the total number of residues in the sequence.

    Running this function multiple times (different key_residue selections, for example) adds suffix to the previously
    mentioned columns following pandas' merge naming logic (_x, _y, _z, ...).

    :param df: designs data.
    :type df: :py:class:`.DesignFrame`
    :param seqID: Identifier of the sequence sets of interest.
    :type seqID: :py:class:`str`
    :param matrix: Identifier of the matrix used to evaluate similarity. Default is BLOSUM62.
    :type matrix: :py:class:`str`

    :return: :py:class:`.DesignFrame`.

    :raises:
        :AttributeError: if the data passed is not a :py:class:`.DesignFrame`.
        :KeyError: if ``seqID`` cannot be found.
        :AttributeError: if there is no reference for ``seqID``.
    """
    if not isinstance(df, DesignFrame):
        raise AttributeError("Input data has to be a DesignFrame with a reference sequence.")
    if not df.has_reference_sequence(seqID):
        raise AttributeError("There is no reference sequence for seqID {}".format(seqID))
    if not "sequence_{}".format(seqID) in df:
        raise KeyError("Sequence {} not found in decoys.".format(seqID))

    mat = SM.get_matrix(matrix)
    ref_seq = df.key_reference_sequence(seqID, key_residues)
    ref_raw, ref_id, ref_pos, ref_neg = _sequence_similarity(ref_seq, ref_seq, mat)
    df2 = pd.DataFrame(pick_key_residues(df, seqID, key_residues).apply( lambda x: _sequence_similarity(x["sequence_{}".format(seqID)], ref_seq, mat), axis=1).values.tolist())

    df2[4] = df2[0] / ref_raw
    df2 = df2.rename(pd.Series(["{0}_{1}_raw".format(matrix.lower(), seqID),
                                "{0}_{1}_identity".format(matrix.lower(), seqID),
                                "{0}_{1}_positive".format(matrix.lower(), seqID),
                                "{0}_{1}_negative".format(matrix.lower(), seqID),
                                "{0}_{1}_perc".format(matrix.lower(), seqID)]), axis="columns").rename(pd.Series(df.index.values), axis="rows")
    return df.merge(df2, left_index=True, right_index=True)

def positional_similarity( df, seqID=None, ref_seq=None, matrix="BLOSUM62" ):
    """
    Generates per-position match data of the provided sequence over the reference sequence.

    :param df: Input data.
    :type df: Union[:py:class:`.DesignFrame`, :py:class:`.FragmentFrame`]
    :param seqID: Identifier of the sequence sets of interest. Required when input is :py:class:`.DesignFrame`
    :type seqID: :py:class:`str`
    :param ref_seq: Reference sequence. Required when input is :py:class:`.FragmentFrame`
    :type ref_seq: :py:class:`str`
    :param matrix: Identifier of the matrix used to evaluate similarity. Default is BLOSUM62.
    :type matrix: :py:class:`str`

    :return: :py:class:`.DesignFrame` where rows are positions and columns are `identity_perc` and `positive_perc`.

    :raises:
        :AttributeError: if the data passed is not in Union[:py:class:`.DesignFrame`, :py:class:`.FragmentFrame`].
        :AttributeError: if input is :py:class:`.DesignFrame` and ``seqID`` is not provided.
        :KeyError: if input is :py:class:`.DesignFrame` and ``seqID`` cannot be found.
        :AttributeError: if input is :py:class:`.DesignFrame` and there is no reference for ``seqID``.
        :AttributeError if input is :py:class:`.FragmentFrame` and ``ref_seq`` is not provided.
    """
    data = {"identity_perc": [], "positive_perc": []}
    mat = SM.get_matrix(matrix)
    if isinstance(df, DesignFrame):
        if seqID is None:
            raise AttributeError("seqID needs to be provided")
        if not df.has_reference_sequence(seqID):
            raise AttributeError("There is no reference sequence for seqID {}".format(seqID))
        if not "sequence_{}".format(seqID) in df:
            raise KeyError("Sequence {} not found in decoys.".format(seqID))

    elif isinstance(df, FragmentFrame):
        if ref_seq is None:
            raise AttributeError("ref_seq needs to be provided")

        for i in df["position"].drop_duplicates().values:
            qseq = "".join(df[df["position"] == i]["aa"].values)
            raw, idn, pos, neg = _positional_similarity( qseq, ref_seq[i-1], mat )
            data["identity_perc"].append(float(idn)/float(len(qseq)))
            data["positive_perc"].append(float(pos)/float(len(qseq)))

    else:
        raise AttributeError("Input data has to be a DesignFrame with a reference sequence or a FragmentFrame.")

    return pd.DataFrame(data)

def binary_similarity( df, ref_seq, matrix="IDENTITY", seq_column="sequence", prefix=None, key_residues=None ):
    mat       = SM.get_matrix(matrix)
    all_seqs  = df[[seq_column]]
    sims      = []
    if prefix is None:
        prefix = matrix.lower()
    else:
        prefix = str(prefix) + matrix.lower()
    for v in all_seqs.values:
        sims.append( _calculate_binary_sequence_similarity( v[0], ref_seq, mat, key_residues ) )
    wdf = df.copy()
    if prefix + "_binary" in df:
        wdf = wdf.drop([prefix + "_binary"], axis=1)
    wdf.insert( wdf.shape[1], prefix + "_binary", pd.Series( sims, index=wdf.index ) )
    return wdf

def binary_overlap( df, column_name="identity_binary" ):
    a = df[column_name].values
    x = len(a[0])
    result = [0] * x
    for seq in a:
        for _, b in enumerate(seq):
            if bool(int(b)): result[_] = 1
    return result

def sequence_frequency_matrix( series, seq_column="sequence" ):
    sserie = series[seq_column].values
    table = {
        'C' : [], 'D' : [], 'S' : [], 'Q' : [], 'K' : [],
        'I' : [], 'P' : [], 'T' : [], 'F' : [], 'N' : [],
        'G' : [], 'H' : [], 'L' : [], 'R' : [], 'W' : [],
        'A' : [], 'V' : [], 'E' : [], 'Y' : [], 'M' : []
    }

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

    df = pd.DataFrame( table )
    df.index = df.index + 1
    return df
