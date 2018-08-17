# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: positional_structural_count
.. func:: positional_structural_identity
.. func:: secondary_structure_percentage
"""

# Standard Libraries
import collections

# External Libraries
import pandas as pd

# This Library

__all__ = ['positional_structural_count', 'positional_structural_identity',
           'secondary_structure_percentage']


def positional_structural_count( df, seqID=None, key_residues=None ):
    """Percentage of secondary structure types for each sequence position of all
    decoys.

    The secondary structure dictionary is a minimized one: ``H``, ``E`` and ``L``.

    :param df: |df_param|.
    :type df: Union[:py:class:`.DesignFrame`, :py:class:`.FragmentFrame`]
    :param str seqID: |seqID_param|. Required when input is :class:`.DesignFrame`.
    :param key_residues: |keyres_param|.
    :type key_residues: |keyres_types|

    :return: :class:`~pandas.DataFrame` - where rows are sequence positions and
        columns are the secondary structure identifiers ``H``, ``E``, ``L``.

    :raises:
        :AttributeError: if the data passed is not in Union[:class:`.DesignFrame`,
            :class:`.FragmentFrame`]. It will *not* try to cast a provided
            :class:`~pandas.DataFrame`, as it would not be possible to know into which of
            the two possible inputs it needs to be casted.
        :AttributeError: if input is :class:`.DesignFrame` and ``seqID`` is not provided.
        :KeyError: |sseID_error| when input is :class:`.DesignFrame`.

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: from rstoolbox.analysis import positional_structural_count
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_ssebig.minisilent.gz",
           ...:                         {'scores': ['score'], 'structure': 'C'})
           ...: df = positional_structural_count(df.iloc[1:], 'C')
           ...: df.head()
    """
    from rstoolbox.components import DesignFrame, FragmentFrame
    from rstoolbox.components import get_selection
    data = {"H": [], "E": [], "L": []}

    if isinstance(df, DesignFrame):
        if seqID is None:
            raise AttributeError("seqID needs to be provided")
        if not "structure_{}".format(seqID) in df:
            raise KeyError("Structure {} not found in decoys.".format(seqID))
        seqdata = df.get_sequential_data('structure', seqID)
        seqdata = seqdata.apply(lambda x: pd.Series(list(x)))
        for _, i in enumerate(seqdata.columns.values):
            qseq = "".join(seqdata[i].tolist())
            sse = collections.Counter(qseq)
            data["H"].append(float(sse["H"]) / float(len(qseq)))
            data["E"].append(float(sse["E"]) / float(len(qseq)))
            data["L"].append(float(sse["L"]) / float(len(qseq)))

    elif isinstance(df, FragmentFrame):
        for i in df["position"].drop_duplicates().values:
            qseq = "".join(df[df["position"] == i]["sse"].values).upper()
            sse = collections.Counter(qseq)
            data["H"].append(float(sse["H"]) / float(len(qseq)))
            data["E"].append(float(sse["E"]) / float(len(qseq)))
            data["L"].append(float(sse["L"]) / float(len(qseq)))

    else:
        raise AttributeError("Input data has to be a DesignFrame or a FragmentFrame.")

    dfo = pd.DataFrame(data)
    # Get shift only from DesignFrame; FragmentFrame does not have one
    shft = df.get_reference_shift(seqID) if isinstance(df, DesignFrame) else 1
    # Shift the index so that index == PDB count
    if isinstance(shft, int):
        dfo.index = dfo.index + shft
    else:
        dfo.index = shft
    return dfo.loc[list(get_selection(key_residues, seqID, list(dfo.index)))]


def positional_structural_identity( df, seqID=None, ref_sse=None, key_residues=None ):
    """Per position evaluation of how many times the provided data matches the expected
    ``reference_structure``.

    :param df: |df_param|.
    :type df: Union[:class:`.DesignFrame`, :class:`.FragmentFrame`]
    :param str seqID: |seqID_param|. Required when input is :class:`.DesignFrame`
    :param str ref_sse: Reference sequence. Required when input is :class:`.FragmentFrame`.
        Will overwrite the reference sequence of :class:`.DesignFrame` if provided.
    :param key_residues: |keyres_param|.
    :type key_residues: |keyres_types|

    :return: :class:`~pandas.DataFrame` - where rows are sequence positions and
        columns are ``sse`` (expected secondary structure),
        ``max_sse`` (most represented secondary structure) and
        ``identity_perc`` (percentage of matched secondary structure).

    :raises:
        :AttributeError: if the data passed is not in Union[:class:`.DesignFrame`,
            :class:`.FragmentFrame`]. It will *not* try to cast a provided
            :class:`~pandas.DataFrame`, as it would not be possible to know into which of
            the two possible inputs it needs to be casted.
        :AttributeError: if input is :class:`.DesignFrame` and ``seqID`` is not provided.
        :KeyError: |sseID_error| when input is :class:`.DesignFrame`.
        :AttributeError: if input is :class:`.FragmentFrame` and ``ref_sse`` is not provided.

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: from rstoolbox.analysis import positional_structural_identity
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_ssebig.minisilent.gz",
           ...:                         {'scores': ['score'], 'structure': 'C'})
           ...: df.add_reference_structure('C', df.get_structure('C').values[0])
           ...: df = positional_structural_identity(df.iloc[1:], 'C')
           ...: df.head()
    """
    from rstoolbox.components import DesignFrame, FragmentFrame
    from rstoolbox.components import get_selection
    data = {"sse": [], "max_sse": [], "identity_perc": []}

    if isinstance(df, DesignFrame):
        if seqID is None:
            raise AttributeError("seqID needs to be provided")
        if not df.has_reference_structure(seqID):
            raise AttributeError("There is no reference structure for seqID {}".format(seqID))
        if not "structure_{}".format(seqID) in df:
            raise KeyError("Structure {} not found in decoys.".format(seqID))
        ref_sse = ref_sse if ref_sse is not None else df.get_reference_structure(seqID)
        seqdata = df.get_structure(seqID)
        seqdata = seqdata.apply(lambda x: pd.Series(list(x)))
        for _, i in enumerate(seqdata.columns.values):
            qseq = "".join(seqdata[i].tolist())
            sse = collections.Counter(qseq)
            data["sse"].append(ref_sse[i])
            data["max_sse"].append(sse.most_common(1)[0][0])
            data["identity_perc"].append(float(sse[ref_sse[i - 1]]) / float(len(qseq)))

    elif isinstance(df, FragmentFrame):
        if ref_sse is None:
            raise AttributeError("ref_sse needs to be provided")

        for i in df["position"].drop_duplicates().values:
            qseq = "".join(df[df["position"] == i]["sse"].values).upper()
            sse = collections.Counter(qseq)
            data["sse"].append(ref_sse[i - 1])
            data["max_sse"].append(sse.most_common(1)[0][0])
            data["identity_perc"].append(float(sse[ref_sse[i - 1]]) / float(len(qseq)))

    else:
        raise AttributeError("Input data has to be a DesignFrame with a reference sequence "
                             "or a FragmentFrame.")

    dfo = pd.DataFrame(data)
    # Get shift only from DesignFrame; FragmentFrame does not have one
    shft = df.get_reference_shift(seqID) if isinstance(df, DesignFrame) else 1
    # Shift the index so that index == PDB count
    if isinstance(shft, int):
        dfo.index = dfo.index + shft
    else:
        dfo.index = shft
    return dfo.loc[list(get_selection(key_residues, seqID, list(dfo.index)))]


def secondary_structure_percentage( df, seqID, key_residues=None ):
    """Calculate the percentage of the different secondary structure types.

    Requires secondary structure data.

    Adds 3 new columns to the data container:

    ===============================  ===================================================
    New Column                       Data Content
    ===============================  ===================================================
    **structure_<seqID>_H**          Percentage of **alpha helices** in the structure.
    **structure_<seqID>_E**          Percentage of **beta sheets** in the structure.
    **structure_<seqID>_L**          Percentage of **loops** in the structure.
    ===============================  ===================================================

    :param df: |df_param|.
    :type df: Union[:py:class:`.DesignFrame`, :py:class:`.DesignSeries`]
    :param str seqID: |seqID_param|.
    :param key_residues: |keyres_param|.
    :type key_residues: |keyres_types|

    :return: Union[:py:class:`.DesignFrame`, :py:class:`.DesignSeries`]

    :raises:
        :NotImplementedError: if the data passed is not in Union[:class:`.DesignFrame`,
            :class:`.DesignSeries`].
        :KeyError: |sseID_error|.

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: from rstoolbox.analysis import secondary_structure_percentage
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_ssebig.minisilent.gz",
           ...:                         {'scores': ['score'], 'structure': 'C'})
           ...: df = secondary_structure_percentage(df, 'C')
           ...: df.head()
    """
    from rstoolbox.components import DesignFrame, DesignSeries
    H = 'structure_{}_H'.format(seqID)
    E = 'structure_{}_E'.format(seqID)
    L = 'structure_{}_L'.format(seqID)

    if isinstance(df, DesignFrame):
        df2 = df.apply(lambda row: secondary_structure_percentage(row,
                                                                  seqID,
                                                                  key_residues),
                       axis=1, result_type='expand')
        df2.transfer_reference(df)
        return df2
    elif isinstance(df, DesignSeries):
        sse = list(df.get_structure(seqID, key_residues))
        csse = collections.Counter(sse)
        dfp = df.append(pd.Series([float(csse['H']) / len(sse), float(csse['E']) / len(sse),
                        float(csse['L']) / len(sse)], [H, E, L]))
        dfp.transfer_reference(df)
        return dfp
    else:
        raise NotImplementedError
