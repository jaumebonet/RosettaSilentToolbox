# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: get_available_references
.. func:: has_reference_sequence
.. func:: add_reference_sequence
.. func:: get_reference_sequence
.. func:: has_reference_structure
.. func:: add_reference_structure
.. func:: get_reference_structure
.. func:: get_reference_shift
.. func:: add_reference_shift
.. func:: add_reference
.. func:: transfer_reference
.. func:: delete_reference
"""
# Standard Libraries
import copy
import warnings

# External Libraries
import six
import pandas as pd
import numpy as np

# This Library


__all__ = ['get_available_references', 'has_reference_sequence',
           'add_reference_sequence', 'get_reference_sequence',
           'has_reference_structure', 'add_reference_structure',
           'get_reference_structure', 'get_reference_shift',
           'add_reference_shift', 'add_reference', 'transfer_reference',
           'delete_reference']


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

    # -1 as we are accessing string count
    return "".join(np.array(list(seq))[kr - 1])


def get_available_references( self ):
    """List which **decoy chain** identifiers have some kind of reference data.

    :return: :func:`list` of :class:`str`

    :TypeError: If applied over a data container without ``_reference`` attribute.

    .. seealso::
        :meth:`.DesignFrame.has_reference_sequence`
        :meth:`.DesignFrame.has_reference_structure`
        :meth:`.DesignSeries.has_reference_sequence`
        :meth:`.DesignSeries.has_reference_structure`
    """
    return list(self._reference.keys())


def has_reference_sequence( self, seqID ):
    """Checks if there is a ``reference_sequence`` for ``sequID``.

    :param str seqID: |seqID_param|.

    :return: :class:`bool`

    :raises:
        :TypeError: If applied over a data container without ``_reference`` attribute.

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz",
           ...:                         {'sequence': 'AB'})
           ...: df.add_reference_sequence('A', df.iloc[0].get_sequence('A'))
           ...: df.has_reference_sequence('A')
           ...: df.has_reference_sequence('B')
    """
    return _has_reference(self, "sequence", seqID)


def add_reference_sequence( self, seqID, sequence ):
    """Add a ``reference_sequence`` attached to chain ``seqID``.

    :param str seqID: |seqID_param|.
    :param str sequence: Reference sequence.

    :raises:
        :TypeError: |indf_error|.
        :ValueError: if ``sequence`` is not a :class:`str`
        :KeyError: |seqID_error|.
        :IndexError: If a ``reference_structure`` for ``seqID`` exists
            but length does not match.

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz",
           ...:                         {'sequence': 'AB'})
           ...: df.add_reference_sequence('A', df.iloc[0].get_sequence('A'))
    """
    if not isinstance(self, (pd.DataFrame, pd.Series)):
        raise TypeError("Data container has to be a DataFrame/Series or a derived class.")
    if self._subtyp != "sequence_frame" and seqID not in self.get_available_sequences():
        raise KeyError("Data container does not have data for sequence {}".format(seqID))
    if not isinstance(sequence, six.string_types):
        if isinstance(sequence, pd.Series) and sequence.size == 1:
            sequence = sequence.values[0]
        else:
            raise ValueError("Reference sequence must be a string.")

    if seqID in self._reference:
        seq = len(self._reference[seqID]["str"])
        if seq > 0 and len(sequence) != seq:
            raise IndexError("Structure length do not match sequence length")
        self._reference[seqID]["seq"] = sequence
    else:
        self._reference.setdefault(seqID, {"seq": sequence, "str": "", "sft": 1})


def get_reference_sequence( self, seqID, key_residues=None ):
    """Get the ``reference_sequence`` attached to chain ``seqID``.

    :param str seqID: |seqID_param|.
    :param key_residues: |keyres_param|.
    :type key_residues: |keyres_types|

    :raises:
        :TypeError: |indf_error|.
        :KeyError: |reference_error|.

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz",
           ...:                         {'sequence': 'AB'})
           ...: df.add_reference_sequence('A', df.iloc[0].get_sequence('A'))
           ...: df.get_reference_sequence('A')
    """
    if key_residues is None:
        return _get_reference(self, "sequence", seqID)
    else:
        return _get_key_reference(self, "sequence", seqID, key_residues)


def has_reference_structure( self, seqID ):
    """Checks if there is a ``reference_structure`` for ``sequID``.

    :param str seqID: |seqID_param|.

    :return: :class:`bool`

    :raises:
        :TypeError: If applied over a data container without ``_reference`` attribute.

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_ssebig.minisilent.gz",
           ...:                         {'sequence': 'C', 'structure': 'C'})
           ...: df.add_reference_structure('C', df.iloc[0].get_structure('C'))
           ...: df.has_reference_structure('C')
    """
    return _has_reference(self, "structure", seqID)


def add_reference_structure( self, seqID, structure ):
    """Add a ``reference_structure`` attached to chain ``seqID``.

    :param str seqID: |seqID_param|.
    :param str structure: Reference structure.

    :raises:
        :TypeError: |indf_error|.
        :ValueError: if ``structure`` is not a :class:`str`
        :KeyError: |seqID_error|.
        :IndexError: If a ``reference_sequence`` for ``seqID`` exists
            but length does not match.

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_ssebig.minisilent.gz",
           ...:                         {'sequence': 'C', 'structure': 'C'})
           ...: df.add_reference_structure('C', df.iloc[0].get_structure('C'))
    """
    if not isinstance(self, (pd.DataFrame, pd.Series)):
        raise TypeError("Data container has to be a DataFrame/Series or a derived class.")
    if self._subtyp != "sequence_frame" and seqID not in self.get_available_structures():
        raise KeyError("Data container does not have data for structure {}".format(seqID))
    if not isinstance(structure, six.string_types):
        if isinstance(structure, pd.Series) and structure.size == 1:
            structure = structure.values[0]
        else:
            raise ValueError("Reference structure must be a string.")

    if seqID in self._reference:
        seq = len(self._reference[seqID]["seq"])
        if seq > 0 and len(structure) != seq:
            raise IndexError("Structure length do not match sequence length")
        self._reference[seqID]["str"] = structure
    else:
        self._reference.setdefault(seqID, {"str": structure, "seq": "", "sft": 1})


def get_reference_structure( self, seqID, key_residues=None ):
    """Get the ``reference_structure`` attached to chain ``seqID``.

    :param str seqID: |seqID_param|.
    :param key_residues: |keyres_param|.
    :type key_residues: |keyres_types|

    :raises:
        :TypeError: |indf_error|.
        :KeyError: |referencestr_error|.

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_ssebig.minisilent.gz",
           ...:                         {'sequence': 'C', 'structure': 'C'})
           ...: df.add_reference_structure('C', df.iloc[0].get_structure('C'))
           ...: df.get_reference_structure('C')
    """
    if key_residues is None:
        return _get_reference( self, "structure", seqID )
    else:
        return _get_key_reference(self, "structure", seqID, key_residues)


def add_reference_shift( self, seqID, shift, shift_labels=False ):
    """Add a ``reference_shift`` attached to a chain ``seqID``.

    **What is shift?** In case the sequence does not start in 1,
    shift defines how to count it. It is a way to keep plotting
    an analysis showing residue number compatible with the PDB.
    There are two main ways to set the shift:

    #. Provide the number of the first residue of the chain; the \
    rest will be set up from there. This is the simplest option, \
    when its just a matter of the structure not actually starting \
    at the begining of the real protein sequence.
    #. If the original PDB has breaks, one will need to provide an \
    array with the numbers of each position. This is the only solution \
    to consistently track those positions.

    :param str seqID: |seqID_param|.
    :param shift: Starting residue number or per-residue number assignment.
    :type shift: Union[:class:`int`, :func:`list` of :class:`int`]
    :param bool shift_labels: When adding the shift, should it be automatically
        applied to any label present in the data container? (Default is :data:`False`).

    :raises:
        :TypeError: |indf_error|.
        :KeyError: If there is no reference structure or sequence for seqID.
        :KeyError: If shift is a list and the data container does not contain structure
            or sequence data for the given seqID.
        :IndexError: If shift is a list and the length is different than the reference
            sequence/structure

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_ssebig.minisilent.gz",
           ...:                         {'sequence': 'C', 'structure': 'C'})
           ...: df.add_reference_structure('C', df.iloc[0].get_structure('C'))
           ...: df.add_reference_shift('C', 3)
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
    """Get a ``reference_shift`` attached to a particular ``seqID``.

    If none was provided, it will return **1** as default.

    :param str seqID: |seqID_param|.

    :type shift: Union[:class:`int`, :class:`list`]

    :raises:
        :TypeError: |indf_error|.

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_ssebig.minisilent.gz",
           ...:                         {'sequence': 'C', 'structure': 'C'})
           ...: df.add_reference_structure('C', df.iloc[0].get_structure('C'))
           ...: df.add_reference_shift('C', 3)
           ...: df.get_reference_shift('C')
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
    """Single access to :meth:`.add_reference_sequence`, :meth:`.add_reference_structure`
    and :meth:`.add_reference_shift`.

    :param str seqID: |seqID_param|.
    :param str sequence: Reference sequence.
    :param str structure: Reference structure.
    :param shift: Starting residue number or per-residue number assignment.
    :type shift: Union[:class:`int`, :class:`list`]
    :param bool shift_labels: When adding the shift, should it be automatically
        applied to any label present in the data container? (Default is :data:`False`).

    :raises:
        :AttributeError: if trying to add a reference type to/from a class without it.

    .. seealso::
        :meth:`DesignFrame.add_reference_sequence`
        :meth:`DesignFrame.add_reference_structure`
        :meth:`DesignFrame.add_reference_shift`
        :meth:`DesignSeries.add_reference_sequence`
        :meth:`DesignSeries.add_reference_structure`
        :meth:`DesignSeries.add_reference_shift`
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
    """Transfer reference data from one container to another.

    This **overwrittes** previous references completely. **Use with care.**

    :param df: |df_param|.
    :type df: Union[:class:`~pandas.DataFrame`, :class:`~pandas.Series`]

    :raises:
        :AttributeError: if trying to add a reference type to/from a class without it.
    """
    try:
        self._reference = copy.deepcopy(df._reference)
    except AttributeError:
        warnings.warn( "Either the target or source object cannot hold a reference." )


def delete_reference( self, seqID, shift_labels=False ):
    """Remove all reference data regarding a particular seqID.

    :param str seqID: |seqID_param|.
    :param bool shift_labels: When removing the shift, should automatically
        revert it to any label present in the data container? (Default is :data:`False`).

    :raises:
        :TypeError: |indf_error|
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
