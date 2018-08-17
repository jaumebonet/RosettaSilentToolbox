# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. class:: DesignSeries
.. class:: DesignFrame
"""
# Standard Libraries
from distutils.version import LooseVersion
import itertools

# External Libraries
import pandas as pd
import numpy as np

# This Library
from .rsbase import RSBaseDesign
import rstoolbox.analysis as ra


__all__ = ["DesignSeries", "DesignFrame"]


if LooseVersion(pd.__version__) < LooseVersion("0.23"):
    raise ImportError(
        'pandas>=0.23 is required!\nSeems that the library '
        'was not force-updated by setup. Execute: \n'
        '   pip install -U pandas\n'
        'to obtain the latest version.'
    )


def _metadata_defaults(name):
    if name == "_source_files":
        return set()
    if name == "_reference":
        return {}
    return None


class DesignSeries( pd.Series, RSBaseDesign ):
    """
    The :class:`.DesignSeries` extends the :class:`~pandas.Series`
    adding some functionalities in order to improve its usability in
    the analysis of a single design decoys.

    It is generated as the **reduced dimensionality version** of the
    :class:`.DesignFrame`.

    .. seealso::
        :class:`.DesignFrame`
    """

    _metadata = ['_reference', '_source_files']
    _subtyp = 'design_series'

    def __init__( self, *args, **kwargs ):
        reference = kwargs.pop('reference', {})
        super(DesignSeries, self).__init__(*args, **kwargs)
        self._reference = reference

    @property
    def _constructor( self ):
        return DesignSeries

    @property
    def _constructor_expanddim( self ):
        def f(*args, **kwargs):
            return DesignFrame(*args, **kwargs).__finalize__(self, method='inherit')
        return f

    def __finalize__(self, other, method=None, **kwargs):
        if method == "inherit":
            # Avoid columns from DesignFrame to become DesignSeries
            if not isinstance(self.name, (int, np.int64)):
                return pd.Series(self)

        for name in self._metadata:
            setattr(self, name, getattr(other, name, _metadata_defaults(name)))
        return self


class DesignFrame( pd.DataFrame, RSBaseDesign ):
    """
    The :class:`.DesignFrame` extends the :class:`~pandas.DataFrame`
    adding some functionalities in order to improve its usability in
    the analysis of sets of design decoys.

    Filled through the functions provided through this library, each
    row represents a decoy while each column represents the scores
    attached to it.

    As a rule, it is assumed that the object:

    #. has a column named **description** \
    that stores the identifier of the corresponding decoy.
    #. holds sequences (a design decoy might be composed of \
    multiple chains) in columns named **sequence_<seqID>**.

    This two assumptions are easily adapted if casting a
    :class:`~pandas.DataFrame` into the class, and several functions
    of the library depend on them.

    .. note::
        This assumptions are automatically fulfilled when the data container is
        loaded through :func:`.parse_rosetta_file`. To obtain sequence information
        is is necessary to request for that particular data, as described in
        :ref:`tutorial: reading Rosetta <readrosetta>`.

    The :class:`.DesignFrame` basically contains four extra attributes
    (accessible through the appropiate functions):

    #. **reference_sequence:** A reference sequence can be added for each ``seqID`` \
    present in the :class:`.DesignFrame`. By adding this sequence, other functions \
    of the library can add that information to its calculations.
    #. **reference_structure:** A reference secondary structure can be added for each \
    ``seqID`` present in the :class:`.DesignFrame`. By adding this sequence, other \
    functions of the library can add that information to its calculations.
    #. **reference_shift:** A reference shift can be added for each ``seqID`` \
    present in the :class:`.DesignFrame`. In short, this would be the initial number \
    of the protein in the source PDB. This allows working with the right numbering. This \
    value is, by default, 1 in all ``seqID``. A more complex alternative allows for a \
    list of numbers to also be assigned as ``reference_shift``. This is usefull when \
    the original structure does not have a continuous numbering schema.
    #. **source_files:** The object stores the source files from which it has been loaded \
    (as long as it is loaded with :func:`.parse_rosetta_file`). This information can \
    be used to extract the structures from the silent files.

    .. seealso::
        :func:`.parse_rosetta_file`

    """
    _metadata = ['_reference', '_source_files']
    _subtyp = 'design_frame'

    def __init__(self, *args, **kwargs):

        if len(args) > 0 or 'data' in kwargs:
            d = args[0] if 'data' not in kwargs else kwargs['data']
            if isinstance(d, (DesignSeries, DesignFrame)):
                kwargs.setdefault('reference', d._reference)
                kwargs.setdefault('source', d._source_files)

        reference = kwargs.pop('reference', {})
        source    = kwargs.pop('source', set())
        super(DesignFrame, self).__init__(*args, **kwargs)
        self._reference = reference
        self._source_files = source

    def get_sequence_with( self, seqID, selection, confidence=1 ):
        """Selects those decoys with a particular set of residue matches.

        Basically, is meant to find, for example, all the decoys in which
        position 25 is A and position 46 is T.

        :param str seqID: |seqID_param|.
        :param selection: List of tuples with position and residue
            type (in 1 letter code).
        :type selection: :func:`list` of :class:`tuple`
        :param float confidence: Percentage of the number of the selection
            rules that we expect the matches to fulfill. Default is 1 (all).

        :return: :class:`.DesignFrame` - filtered by the requested sequence

        .. rubric:: Example

        .. ipython::

            In [1]: from rstoolbox.io import parse_rosetta_file
               ...: import pandas as pd
               ...: pd.set_option('display.width', 1000)
               ...: pd.set_option('display.max_columns', 500)
               ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz",
               ...:                         {'scores': ['score'], 'sequence': 'B'})
               ...: df.get_sequence_with('B', [(1, 'T')])
        """
        from .selection import get_selection

        def match_residues( value, seqID, confidence ):
            t = 0
            for s in selection:
                # -1 as we access string position directly
                t += 1 if value[get_selection(s[0], seqID)[0] - 1] == s[1] else 0
            return t / float(len(selection)) >= float(confidence)

        return self.loc[ self.apply(
            lambda row: match_residues(row["sequence_{0}".format(seqID)], selection, confidence ),
            axis=1
        )]

    def sequence_distance( self, seqID, other=None ):
        """Make identity sequence distance between the selected decoys.

        Generate a matrix counting the distance between each pair of sequences in the
        :class:`.DesignFrame`. This is a time-consuming operation; better to execute it over a
        specific set of selected decoys that over all your designs.

        If ``other`` is provided as a second :class:`.DesignFrame`, distances are calculated between
        the sequences of the current :class:`.DesignFrame` against the sequence of the other.

        :param str seqID: |seqID_param|.
        :param other: Secondary data container. Optional.
        :type other: :class:`.DesignFrame`

        return: :class:`~pandas.DataFrame` - table with the sequence distances.

        :raises:
            :KeyError: |seqID_error|.
            :KeyError: if ``description`` column cannot be found.
            :ValueError: if sequence of ``self`` and ``other`` are of different length.
            :ValueError: if data container only has one sequence and no ``other`` is provided.

        .. rubric:: Example

        .. ipython::

            In [1]: from rstoolbox.io import parse_rosetta_file
               ...: import pandas as pd
               ...: pd.set_option('display.width', 1000)
               ...: pd.set_option('display.max_columns', 500)
               ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz",
               ...:                         {'scores': ['score', 'description'], 'sequence': 'B'})
               ...: df.sequence_distance('B')
        """
        # https://stackoverflow.com/questions/49235538/pandas-series-compare-values-all-vs-all
        def own_distance( self, seqID ):
            if self.shape[0] == 1:
                raise ValueError("More than one sequence is needed to compare.")
            a = np.array(list(map(list, self.get_sequence(seqID).values)))
            ids = self.get_id().values
            n = len(a)
            d = {(i, j): np.sum(a[i] != a[j]) for i in range(n) for j in range(n) if j > i}
            res = np.zeros((n, n))
            keys = list(zip(*d.keys()))
            res[keys[0], keys[1]] = list(d.values())
            res += res.T
            return pd.DataFrame(res, columns=ids, index=ids, dtype=int)

        def vs_distance( self, seqID, other ):
            a = np.array(list(map(list, self.get_sequence(seqID).values)))
            b = np.array(list(map(list, other.get_sequence(seqID).values)))
            aids = self.get_id().values
            bids = other.get_id().values
            na = len(a)
            nb = len(b)
            if len(a[0]) != len(b[0]):
                raise ValueError('Comparable sequence have to be of the same size')
            d = {(i, j): np.sum(a[i] != b[j]) for i in range(na) for j in range(nb)}
            res = np.zeros((na, nb))
            keys = list(zip(*d.keys()))
            res[keys[0], keys[1]] = list(d.values())
            return pd.DataFrame(res, columns=bids, index=aids, dtype=int)

        if other is None:
            return own_distance(self, seqID)
        else:
            return vs_distance(self, seqID, other)

    def sequence_frequencies( self, seqID, seqType="protein", cleanExtra=True, cleanUnused=-1 ):
        """Create a frequency-based :class:`.SequenceFrame`.

        Generates a :class:`.SequenceFrame` for the frequencies of
        the sequences in the :class:`.designFrame` with ``seqID`` identifier.
        If there is a reference_sequence for this ``seqID``, it will also
        be attached to the :class:`.SequenceFrame`.

        All letters in the sequence will be capitalized. All symbols that
        do not belong to ``string.ascii_uppercase`` will be transformed to "*"
        as this is the symbol recognized by the substitution matrices.

        :param str seqID: |seqID_param|.
        :param str seqType: Type of sequence: ``protein``, ``dna``, ``rna``.
        :param bool cleanExtra: Remove from the class:`.SequenceFrame` the non-regular
            amino/nucleic acids if they are empty for all positions.
        :param int cleanUnused: Remove from the class:`.SequenceFrame` the regular
            amino/nucleic acids if they frequency is equal or under the value . Default is -1,
            so nothing is deleted.

        :return: :class:`.SequenceFrame`

        .. seealso::
            :meth:`.DesignFrame.sequence_bits`
            :meth:`.DesignFrame.structure_frequencies`
            :meth:`.DesignFrame.structure_bits`
            :func:`.sequential_frequencies`
        """
        return ra.sequential_frequencies(self, seqID, "sequence", seqType, cleanExtra, cleanUnused)

    def sequence_bits( self, seqID, seqType="protein", cleanExtra=True, cleanUnused=False ):
        """Create a bit-based :class:`.SequenceFrame`.

        Generates a :class:`.SequenceFrame` for the bits of
        the sequences in the :class:`.designFrame` with ``seqID`` identifier.
        If there is a reference_sequence for this ``seqID``, it will also
        be attached to the :class:`.SequenceFrame`.

        Bit calculation is performed as explained in http://www.genome.org/cgi/doi/10.1101/gr.849004
        such as::

            Rseq = Smax - Sobs = log2 N - (-sum(n=1,N):pn * log2 pn)

        Where:
            - N is the total number of options (4: DNA/RNA; 20: PROTEIN).
            - pn is the observed frequency of the symbol n.

        :param str seqID: |seqID_param|.
        :param str seqType: Type of sequence: ``protein``, ``dna``, ``rna``.
        :param bool cleanExtra: Remove from the class:`.SequenceFrame` the non-regular
            amino/nucleic acids if they are empty for all positions.
        :param int cleanUnused: Remove from the class:`.SequenceFrame` the regular
            amino/nucleic acids if they frequency is equal or under the value . Default is -1,
            so nothing is deleted.

        :return: :class:`.SequenceFrame`

        .. seealso::
            :meth:`.DesignFrame.sequence_frequencies`
            :meth:`.DesignFrame.structure_frequencies`
            :meth:`.DesignFrame.structure_bits`
            :func:`.sequential_frequencies`
        """
        df = self.sequence_frequencies(seqID, seqType, cleanExtra, cleanUnused)
        return df.to_bits()

    def structure_frequencies( self, seqID, seqType="protein", cleanExtra=True, cleanUnused=-1 ):
        """Create a frequency-based :class:`.SequenceFrame` for secondary structure assignation.

        Generates a :class:`.SequenceFrame` for the frequencies of
        the secondary structure in the :class:`.designFrame` with ``seqID`` identifier.
        If there is a reference_structure for this ``seqID``, it will also
        be attached to the :class:`.SequenceFrame`.

        All letters in the secondary structure will be capitalized. All symbols that
        do not belong to string.ascii_uppercase will be transformed to "*"
        as this is the symbol recognized by the substitution matrices.

        :param str seqID: |seqID_param|.
        :param str seqType: Type of sequence: ``protein``, ``dna``, ``rna``.
        :param bool cleanExtra: Remove from the class:`.SequenceFrame` the non-regular
            amino/nucleic acids if they are empty for all positions.
        :param int cleanUnused: Remove from the class:`.SequenceFrame` the regular
            amino/nucleic acids if they frequency is equal or under the value . Default is -1,
            so nothing is deleted.

        :return: :class:`.SequenceFrame`

        .. seealso::
            :meth:`.DesignFrame.sequence_frequencies`
            :meth:`.DesignFrame.sequence_bits`
            :meth:`.DesignFrame.structure_bits`
            :func:`.sequential_frequencies`
        """
        seqType = seqType + "_sse"
        return ra.sequential_frequencies(self, seqID, "structure", seqType, cleanExtra, cleanUnused)

    def structure_bits( self, seqID, seqType="protein", cleanExtra=True, cleanUnused=False ):
        """Create a bit-based :class:`.SequenceFrame` for secondary structure assignation.

        Generates a :class:`.SequenceFrame` for the bits of
        the secondary structure in the :class:`.designFrame` with ``seqID`` identifier.
        If there is a reference_structure for this ``seqID``, it will also
        be attached to the :class:`.SequenceFrame`.

        Bit calculation is performed as explained in http://www.genome.org/cgi/doi/10.1101/gr.849004
        such as::

            Rseq = Smax - Sobs = log2 N - (-sum(n=1,N):pn * log2 pn)

        Where:
            - N is the total number of options (4: DNA/RNA; 20: PROTEIN).
            - pn is the observed frequency of the symbol n.

        :param str seqID: |seqID_param|.
        :param str seqType: Type of sequence: ``protein``, ``dna``, ``rna``.
        :param bool cleanExtra: Remove from the class:`.SequenceFrame` the non-regular
            amino/nucleic acids if they are empty for all positions.
        :param int cleanUnused: Remove from the class:`.SequenceFrame` the regular
            amino/nucleic acids if they frequency is equal or under the value . Default is -1,
            so nothing is deleted.

        :return: :class:`.SequenceFrame`

        .. seealso::
            :meth:`.DesignFrame.sequence_frequencies`
            :meth:`.DesignFrame.sequence_bits`
            :meth:`.DesignFrame.structure_frequencies`
            :func:`.sequential_frequencies`
        """
        df = self.structure_frequencies(seqID, seqType, cleanExtra, cleanUnused)
        return df.to_bits()

    #
    # Implement pandas methods
    #

    @property
    def _constructor(self):
        return DesignFrame

    @property
    def _constructor_sliced(self):
        def f(*args, **kwargs):
            return DesignSeries(*args, **kwargs).__finalize__(self, method='inherit')
        return f

    def __finalize__(self, other, method=None, **kwargs):
        """propagate metadata from other to self """
        # strict  = kwargs["strict"] if "strict" in kwargs else True
        # concat operation:
        #   (1) accumulate _source_files from all the concatenated DesignFrames (if there is any)
        #   (2) check reference_sequence and reference_shift, merge and keep only
        #       if they are the same.
        if method == 'concat':
            source_files = set()
            refseqs = []
            reference_sequence = {}
            for i, o in enumerate(other.objs):
                source_files.update(getattr(o, "_source_files", set()))
                refseqs.append(getattr(o, "_reference", {}))
            # _source_files
            setattr(self, "_source_files", source_files)
            # _reference
            ids = list(set(itertools.chain.from_iterable([x.keys() for x in refseqs])))
            for r in refseqs:
                for i in ids:
                    if i in r:
                        if i not in reference_sequence:
                            reference_sequence.setdefault(i, r[i])
                        else:
                            if r[i] != reference_sequence[i]:
                                raise ValueError("Concatenating designFrames with "
                                                 "different ref sequence for the same seqID.")
            setattr(self, "_reference", reference_sequence)
        # merge operation:
        # Keep metadata of the left object.
        elif method == 'merge':
            for name in self._metadata:
                setattr(self, name, getattr(other.left, name, _metadata_defaults(name)))
        # inherit operation:
        # Keep metadata of the other object.
        elif method == 'inherit':
            for name in self._metadata:
                setattr(self, name, getattr(other, name, _metadata_defaults(name)))
        else:
            for name in self._metadata:
                setattr(self, name, getattr(other, name, _metadata_defaults(name)))
        return self
