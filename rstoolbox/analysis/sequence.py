# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: sequential_frequencies
.. func:: sequence_similarity
.. func:: positional_sequence_similarity
.. func:: binary_similarity
.. func:: binary_overlap
.. func:: selector_percentage
.. func:: label_percentage
.. func:: positional_enrichment
"""
# Standard Libraries
import copy
import collections
import re

# External Libraries
import pandas as pd
import numpy as np

# This Library
from .SimilarityMatrix import SimilarityMatrix as SM

__all__ = ['sequential_frequencies', 'sequence_similarity',
           'positional_sequence_similarity', 'binary_similarity',
           'binary_overlap', 'selector_percentage', 'label_percentage',
           'positional_enrichment']


def _get_sequential_table( seqType ):
    """
    Generates the table to fill sequence data in order to create
    a :class:`.SequenceFrame`

    :param seqType: Type of sequence: ``protein``, ``protein_sse``,
        ``dna``, ``rna``.
    :type seqType: :class:`str`

    :return: :class:`dict`

    :raise:
        :ValueError: If ``seqType`` is not known.
    """
    table = {}
    extra = []
    if seqType.lower() == "protein":
        # X = UNKNOWN  # * = GAP
        # B = N or D   # Z = E or Q
        table = {
            'C': [], 'D': [], 'S': [], 'Q': [], 'K': [],
            'I': [], 'P': [], 'T': [], 'F': [], 'N': [],
            'G': [], 'H': [], 'L': [], 'R': [], 'W': [],
            'A': [], 'V': [], 'E': [], 'Y': [], 'M': [],
            'X': [], '*': [], 'B': [], 'Z': []
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
            'C': [], 'A': [], 'T': [], 'G': [], 'X': [], '*': [],
            'B': [], 'D': [], 'H': [], 'K': [], 'M': [], 'N': [],
            'R': [], 'S': [], 'V': [], 'W': [], 'Y': []
        }
        if seqType.lower() == "rna":
            table.setdefault( 'U', [])
            table.pop('T', None)
        extra = ['X', '*', 'B', 'D', 'H', 'K', 'M', 'N', 'R', 'S', 'V', 'W', 'Y']
    elif seqType.lower() == "protein_sse":
        table = {
            'H': [], 'E': [], 'L': [], '*': [], 'G': []
        }
        extra = ['*', 'G']
    else:
        raise ValueError("sequence type {0} unknown".format(seqType))
    return table, extra


def _sequence_similarity( qseq, rseq, matrix ):
    if len(qseq) != len(rseq):
        raise ValueError("Comparable sequences have to be the same size.")
    raw, idn, pos, neg = 0, 0, 0, 0
    ali = []
    pres = []
    for i, qseqi in enumerate(qseq):
        sc = matrix.get_value(qseqi, rseq[i])
        pres.append(sc)
        raw += sc
        if qseqi == rseq[i]:
            idn += 1
            pos += 1
            ali.append(rseq[i])
        elif sc > 0:
            pos += 1
            ali.append("+")
        else:
            neg += 1
            ali.append(".")
    return raw, idn, pos, neg, "".join(ali), pres


def _positional_similarity( qseq, rseq, matrix ):
    raw, idn, pos, neg = 0, 0, 0, 0
    for _, qseqi in enumerate(qseq):
        sc = matrix.get_value(qseqi, rseq)
        raw += sc
        if qseqi == rseq:
            idn += 1
        if sc > 0:
            pos += 1
        else:
            neg += 1
    return raw, idn, pos, neg


def sequential_frequencies( df, seqID, query="sequence", seqType="protein",
                            cleanExtra=True, cleanUnused=-1 ):
    """Generates a :class:`.SequenceFrame` for the frequencies of the sequences in the
    :class:`.DesignFrame` with ``seqID`` identifier.

    If there is a ``reference_sequence`` for this ``seqID``, it will also
    be attached to the :class:`.SequenceFrame`.

    All letters in the sequence will be capitalized. All symbols that
    do not belong to ``string.ascii_uppercase`` will be transformed to `"*"`
    as this is the symbol recognized by the substitution matrices as ``gap``.

    This function is directly accessible through some :class:`.DesignFrame` methods.

    :param df: |df_param|.
    :type df: Union[:class:`.DesignFrame`, :class:`~pandas.DataFrame`]
    :param str seqID: |seqID_param|.
    :param str query: |query_param|.
    :param str seqType: |seqType_param| and ``protein_sse``.
    :param bool cleanExtra: |cleanExtra_param|.
    :param float cleanUnused: |cleanUnused_param|.

    :return: :class:`.SequenceFrame`

    .. seealso::
        :meth:`.DesignFrame.sequence_frequencies`
        :meth:`.DesignFrame.sequence_bits`
        :meth:`.DesignFrame.structure_frequencies`
        :meth:`.DesignFrame.structure_bits`

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: from rstoolbox.analysis import sequential_frequencies
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz",
           ...:                         {'scores': ['score'], 'sequence': 'AB'})
           ...: df = sequential_frequencies(df, 'B')
           ...: df.head()
    """
    from rstoolbox.components import SequenceFrame

    def count_instances( seq, table ):
        t = copy.deepcopy(table)
        c = collections.Counter(seq)
        for aa in table:
            _ = c[aa]
            if _ > 0:
                t[aa] = float(_) / len(seq)
            else:
                t[aa] = 0
        return t

    # Cast if possible, so that we can access the different methods of DesignFrame
    if df._subtyp != 'design_frame' and isinstance(df, pd.DataFrame):
        from rstoolbox.components import DesignFrame
        df = DesignFrame(df)

    # Get all sequences; exclude empty ones (might happen) and uppercase all residues.
    sserie = df.get_sequential_data(query, seqID).replace('', np.nan).dropna().str.upper()
    # Get the table to fill
    table, extra = _get_sequential_table( seqType )
    # Fill the table with the frequencies
    sserie = sserie.apply(lambda x: pd.Series(list(x)))
    sserie = sserie.apply(lambda x: pd.Series(count_instances(x.str.cat(), table))).T

    # Create the SequenceFrame
    dfo = SequenceFrame(sserie)
    dfo.measure("frequency")
    dfo.extras( extra )
    # Attach the reference sequence if there is any
    if df.has_reference_sequence(seqID):
        dfo.add_reference(seqID, sequence=df.get_reference_sequence(seqID),
                          shift=df.get_reference_shift(seqID))
    dfo.delete_extra( cleanExtra )
    dfo.delete_empty( cleanUnused )
    dfo.clean()
    shft = df.get_reference_shift(seqID)
    # Shift the index so that the index of the SequenceFrame == PDB count
    if isinstance(shft, int):
        dfo.index = dfo.index + shft
    else:
        dfo.index = shft
    return dfo


def sequence_similarity( df, seqID, key_residues=None, matrix="BLOSUM62" ):
    """Evaluate the sequence similarity between each decoy and the ``reference_sequence``
    for a given ``seqID``.

    Sequence similarity is understood in the context of substitution matrices. Thus, a part from
    identities, also similarities can be evaluated.

    It will return the input data container with several new columns:

    ===============================  ===================================================
    New Column                       Data Content
    ===============================  ===================================================
    **<matrix>_<seqID>_raw**         Score obtained by applying ``<matrix>``
    **<matrix>_<seqID>_perc**        Score obtained by applying ``<matrix>`` over score \
                                     of reference_sequence against itself
    **<matrix>_<seqID>_identity**    Total identity matches
    **<matrix>_<seqID>_positive**    Total positive matches according to ``<matrix>``
    **<matrix>_<seqID>_negative**    Notal negative matches according to ``<matrix>``
    **<matrix>_<seqID>_ali**         Representation of aligned residues
    **<matrix>_<seqID>_per_res**     Per position score of applying ``<matrix>``
    ===============================  ===================================================

    Matrix name in each new column is setup in lowercase.

    .. tip::
        If ``key_residues`` are applied, the scoring is only used on those, but nothing in the
        naming of the columns will indicate a partial evaluation. It is important to keep that in
        mind moving forward on whatever analysis you are performing.

    Running this function multiple times (different key_residue selections, for example)
    adds suffix to the previously mentioned columns following pandas' merge naming
    logic (_x, _y, _z, ...).

    :param df: |df_param|.
    :type df: Union[:class:`.DesignFrame`, :class:`~pandas.DataFrame`]
    :param str seqID: |seqID_param|.
    :param key_residues: |keyres_param|.
    :type key_residues: |keyres_types|
    :param str matrix: |matrix_param|. Default is ``BLOSUM62``.

    :return: :class:`.DesignFrame`.

    :raises:
        :AttributeError: |designframe_cast_error|.
        :KeyError: |seqID_error|.
        :AttributeError: |reference_error|.

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: from rstoolbox.analysis import sequence_similarity
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz",
           ...:                         {'scores': ['score'], 'sequence': 'B'})
           ...: df.add_reference_sequence('B', df.get_sequence('B').values[0])
           ...: df = sequence_similarity(df.iloc[1:], 'B')
           ...: df.head()

    """
    from rstoolbox.components import DesignFrame

    # We don't need to try to cast, as reference_sequence is needed anyway
    if not isinstance(df, DesignFrame):
        raise AttributeError("Input data has to be a DesignFrame with a reference sequence.")
    if not df.has_reference_sequence(seqID):
        raise AttributeError("There is no reference sequence for seqID {}".format(seqID))
    if not "sequence_{}".format(seqID) in df:
        raise KeyError("Sequence {} not found in decoys.".format(seqID))

    # Get matrix data
    mat = SM.get_matrix(matrix)
    # Get total score of the reference (depending on the matrix, identities != 1)
    ref_seq = df.get_reference_sequence(seqID, key_residues)
    ref_raw, _, _, _, _, _ = _sequence_similarity(ref_seq, ref_seq, mat)
    # Get only the key residues and apply similarity analysis
    df2 = df._constructor(df.get_sequence(seqID, key_residues))
    df2 = df2.apply(
        lambda x: pd.Series(_sequence_similarity(x.get_sequence(seqID), ref_seq, mat)), axis=1)
    df2[6] = df2[0] / ref_raw
    df2 = df2.rename(pd.Series(["{0}_{1}_raw".format(matrix.lower(), seqID),
                                "{0}_{1}_identity".format(matrix.lower(), seqID),
                                "{0}_{1}_positive".format(matrix.lower(), seqID),
                                "{0}_{1}_negative".format(matrix.lower(), seqID),
                                "{0}_{1}_ali".format(matrix.lower(), seqID),
                                "{0}_{1}_per_res".format(matrix.lower(), seqID),
                                "{0}_{1}_perc".format(matrix.lower(), seqID)]),
                     axis="columns").rename(pd.Series(df.index.values), axis="rows")
    return pd.concat([df.reset_index(drop=True),
                      df2.reset_index(drop=True)], axis=1)


def positional_sequence_similarity( df, seqID=None, ref_seq=None,
                                    key_residues=None, matrix="BLOSUM62" ):
    """Per position identity and similarity against a ``reference_sequence``.

    Provided a data container with a set of sequences, it will evaluate the percentage of
    identities and similarities that the whole set has against a ``reference_sequence``.
    It would do so by sequence position instead that by each individual sequence.

    In a way, this generates an extreme simplification from a :class:`.SequenceFrame`.

    :param df: |df_param|.
    :type df: Union[:class:`.DesignFrame`, :class:`.FragmentFrame`]
    :param str seqID: |seqID_param|. Required when input is :class:`.DesignFrame`.
    :param str ref_seq: Reference sequence. Required when input is :class:`.FragmentFrame`.
        Will overwrite the reference sequence of :class:`.DesignFrame` if provided.
    :param key_residues: |keyres_param|.
    :type key_residues: |keyres_types|
    :param str matrix: |matrix_param|. Default is ``BLOSUM62``.


    :return: :class:`~pandas.DataFrame` - where rows are sequence positions and
        columns are ``identity_perc`` and ``positive_perc``.

    :raises:
        :AttributeError: if the data passed is not in Union[:class:`.DesignFrame`,
            :class:`.FragmentFrame`]. It will *not* try to cast a provided
            :class:`~pandas.DataFrame`, as it would not be possible to know into which of
            the two possible inputs it needs to be casted.
        :AttributeError: if input is :class:`.DesignFrame` and ``seqID`` is not provided.
        :KeyError: |seqID_error| when input is :class:`.DesignFrame`.
        :AttributeError: |reference_error| when input is :class:`.DesignFrame`.
        :AttributeError:  if input is :class:`.FragmentFrame` and ``ref_seq`` is not provided.

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: from rstoolbox.analysis import positional_sequence_similarity
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz",
           ...:                         {'scores': ['score'], 'sequence': 'B'})
           ...: df.add_reference_sequence('B', df.get_sequence('B').values[0])
           ...: df = positional_sequence_similarity(df.iloc[1:], 'B')
           ...: df.head()
    """
    from rstoolbox.components import DesignFrame, FragmentFrame
    from rstoolbox.components import get_selection

    data = {"identity_perc": [], "positive_perc": []}
    # Get matrix data
    mat = SM.get_matrix(matrix)

    if isinstance(df, DesignFrame):
        if seqID is None:
            raise AttributeError("seqID needs to be provided")
        if not df.has_reference_sequence(seqID):
            raise AttributeError("There is no reference sequence for seqID {}".format(seqID))
        if not "sequence_{}".format(seqID) in df:
            raise KeyError("Sequence {} not found in decoys.".format(seqID))

        ref_seq = ref_seq if ref_seq is not None else df.get_reference_sequence(seqID)
        seqdata = df.get_sequence(seqID)
        seqdata = seqdata.apply(lambda x: pd.Series(list(x)))
        for _, i in enumerate(seqdata.columns.values):
            qseq = "".join(seqdata[i].tolist())
            _, idn, pos, _ = _positional_similarity( qseq, ref_seq[_], mat )
            data["identity_perc"].append(float(idn) / float(len(qseq)))
            data["positive_perc"].append(float(pos) / float(len(qseq)))

    elif isinstance(df, FragmentFrame):
        if ref_seq is None:
            raise AttributeError("ref_seq needs to be provided")

        for i in df["position"].drop_duplicates().values:
            qseq = "".join(df[df["position"] == i]["aa"].values)
            _, idn, pos, _ = _positional_similarity( qseq, ref_seq[i - 1], mat )
            data["identity_perc"].append(float(idn) / float(len(qseq)))
            data["positive_perc"].append(float(pos) / float(len(qseq)))

    else:
        raise AttributeError("Input data has to be a DesignFrame with a "
                             "reference sequence or a FragmentFrame.")

    dfo = pd.DataFrame(data)
    # Get shift only from DesignFrame; FragmentFrame does not have one
    shft = df.get_reference_shift(seqID) if isinstance(df, DesignFrame) else 1
    # Shift the index so that index == PDB count
    if isinstance(shft, int):
        dfo.index = dfo.index + shft
    else:
        dfo.index = shft
    return dfo.loc[list(get_selection(key_residues, seqID, list(dfo.index)))]


def binary_similarity( df, seqID, key_residues=None, matrix="IDENTITY"):
    """Binary profile for each design sequence against the ``reference_sequence``.

    Makes a :class:`DesignFrame` with a new column to map binary identity (0/1) with
    the ``reference_sequence``. If a different matrix than ``IDENTITY`` is provides,
    the binary sequence sets to ``1`` all the positive values.

    ===============================  ===================================================
    New Column                       Data Content
    ===============================  ===================================================
    **<matrix>_<seqID>_binary**      Binary representation of the match with the
                                     ``reference_sequence``.
    ===============================  ===================================================

    :param df: |df_param|.
    :type df: Union[:class:`.DesignFrame`, :class:`~pandas.DataFrame`]
    :param str seqID: |seqID_param|.
    :param key_residues: |keyres_param|.
    :type key_residues: |keyres_types|
    :param str matrix: |matrix_param|. Default is ``IDENTITY``.

    :return: :class:`.DesignFrame`.

    :raises:
        :AttributeError: |designframe_cast_error|.
        :KeyError: |seqID_error|.
        :AttributeError: |reference_error|.

    .. seealso::
        :func:`.sequence_similarity`

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: from rstoolbox.analysis import binary_similarity
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz",
           ...:                         {'scores': ['score'], 'sequence': 'B'})
           ...: df.add_reference_sequence('B', df.get_sequence('B').values[0])
           ...: df = binary_similarity(df.iloc[1:], 'B')
           ...: df.head()

    """
    dfss = sequence_similarity( df, seqID, key_residues, matrix=matrix )
    alicolumn = "{0}_{1}_ali".format(matrix.lower(), seqID)
    bincolumn = "{0}_{1}_binary".format(matrix.lower(), seqID)

    dfss[bincolumn] = dfss.apply(lambda row: re.sub('\D', '1', re.sub('\.', '0', row[alicolumn])),
                                 axis=1)

    return pd.concat([df.reset_index(drop=True),
                      dfss[bincolumn].reset_index(drop=True)],
                     axis=1)


def binary_overlap( df, seqID, key_residues=None, matrix="IDENTITY" ):
    """Overlap the binary similarity representation of all decoys in a
    :class:`.DesignFrame`.

    :param df: |df_param|.
    :type df: Union[:class:`.DesignFrame`, :class:`~pandas.DataFrame`]
    :param str seqID: |seqID_param|.
    :param key_residues: |keyres_param|.
    :type key_residues: |keyres_types|
    :param str matrix: |matrix_param|. Default is ``IDENTITY``.

    :return: :func:`list` of :class:`int` - ones and zeros for each
        position of the length of the sequence

    .. seealso::
        :func:`.binary_similarity`

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: from rstoolbox.analysis import binary_overlap
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz",
           ...:                         {'scores': ['score'], 'sequence': 'B'})
           ...: df.add_reference_sequence('B', df.get_sequence('B').values[0])
           ...: binoverlap = binary_overlap(df.iloc[1:], 'B')
           ...: "".join([str(_) for _ in binoverlap])
    """
    bincolumn = "{0}_{1}_binary".format(matrix.lower(), seqID)
    if bincolumn not in df.columns.values:
        df = binary_similarity(df, seqID, key_residues, matrix)

    a = df[bincolumn].values
    x = len(a[0])
    result = [0] * x
    for seq in a:
        for _, b in enumerate(seq):
            if bool(int(b)):
                result[_] = 1
    return result


def selector_percentage( df, seqID, key_residues, selection_name='selection' ):
    """Calculate the percentage coverage of a :class:`.Selection`
    over the sequence.

    Depends on sequence information for the ``seqID``.

    Adds a new column to the data container:

    ====================================  =======================================================
    New Column                             Data Content
    ====================================  =======================================================
    **<selection_name>_<seqID>_perc**     Percentage of the sequence covered by the key_residues.
    ====================================  =======================================================

    :param df: |df_param|.
    :type df: Union[:class:`.DesignFrame`, :class:`.DesignSeries`]
    :param str seqID: |seqID_param|.
    :param key_residues: |keyres_param|.
    :type key_residues: |keyres_types|
    :param str selection_name: Prefix to add to the selection. Default is ``selection``.

    :return: Union[:class:`.DesignFrame`, :class:`.DesignSeries`]

    :raises:
        :NotImplementedError: if the data passed is not in Union[:class:`.DesignFrame`,
            :class:`.DesignSeries`].
        :KeyError: |seqID_error|.

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: from rstoolbox.analysis import selector_percentage
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_ssebig.minisilent.gz",
           ...:                         {'scores': ['score'], 'sequence': 'C'})
           ...: df = selector_percentage(df, 'C', '1-15')
           ...: df.head()
    """
    from rstoolbox.components import DesignFrame, DesignSeries

    colname = '{0}_{1}_perc'.format(selection_name, seqID)

    if isinstance(df, DesignFrame):
        df2 = df.apply(lambda row: selector_percentage(row, seqID, key_residues, selection_name),
                       axis=1, result_type='expand')
        return df2
    elif isinstance(df, DesignSeries):
        seq1 = list(df.get_sequence(seqID))
        seq2 = list(df.get_sequence(seqID, key_residues))
        return df.append(pd.Series([float(len(seq2)) / len(seq1)], [colname]))
    else:
        raise NotImplementedError


def label_percentage( df, seqID, label ):
    """Calculate the percentage coverage of a ``label`` over the sequence.

    Depends on sequence information and label data for the ``seqID``.

    Adds a new column to the data container:

    ===========================  ====================================================
    New Column                   Data Content
    ===========================  ====================================================
    **<label>_<seqID>_perc**     Percentage of the sequence covered by the ``label``.
    ===========================  ====================================================

    :param df: |df_param|.
    :type df: Union[:class:`.DesignFrame`, :class:`.DesignSeries`]
    :param str seqID: |seqID_param|.
    :param str lable: Label identifier.
    :param key_residues: |keyres_param|.
    :type key_residues: |keyres_types|

    :return: Union[:class:`.DesignFrame`, :class:`.DesignSeries`]

    :raises:
        :NotImplementedError: if the data passed is not in Union[:class:`.DesignFrame`,
            :class:`.DesignSeries`].
        :KeyError: |lblID_error|.

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: from rstoolbox.analysis import label_percentage
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz",
           ...:                         {'scores': ['score'], 'sequence': '*',
           ...:                          'labels': ['MOTIF']})
           ...: df = label_percentage(df, 'B', 'MOTIF')
           ...: df.head()
    """
    from rstoolbox.components import DesignFrame, DesignSeries
    colname = '{0}_{1}_perc'.format(label.upper(), seqID)

    if isinstance(df, DesignFrame):
        df2 = df.apply(lambda row: label_percentage(row, seqID, label),
                       axis=1, result_type='expand')
        return df2
    elif isinstance(df, DesignSeries):
        seq1 = list(df.get_sequence(seqID))
        seq2 = list(df.get_sequence(seqID, df.get_label(label, seqID)))
        return df.append(pd.Series([float(len(seq2)) / len(seq1)], [colname]))
    else:
        raise NotImplementedError


def positional_enrichment(df, other, seqID):
    """Calculates per-residue enrichment from sequences in the first :class:`.DesignFrame`
    with respect to the second.

    .. note::
        Position / AA type pairs present in ``df`` but not ``other`` will have a value of
        :data:`~np.inf`.

    :param df: |df_param|.
    :type df: Union[:class:`.DesignFrame`, :class:`~pandas.DataFrame`]
    :param other: |df_param|.
    :type other: Union[:class:`.DesignFrame`, :class:`~pandas.DataFrame`]
    :param str seqID: |seqID_param|.

    :return: :class:`.FragmentFrame` - with enrichment percentages.

    :raises:
        :NotImplementedError: if the data passed is not in Union[:class:`.DesignFrame`,
            :class:`~pandas.DataFrame`].
        :KeyError: |seqID_error|.
    """
    from rstoolbox.components import DesignFrame

    for i, x in enumerate([df, other]):
        if not isinstance(x, DesignFrame):
            if not isinstance(x, pd.DataFrame):
                raise NotImplementedError('Unknow input format')
            else:
                if i == 0:
                    df = DesignFrame(df)
                else:
                    other = DesignFrame(other)
    result = df.sequence_frequencies(seqID) / other.sequence_frequencies(seqID)
    if df._reference == other._reference:
        result.transfer_reference(df)
    return result.replace(np.nan, 0)
