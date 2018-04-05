# @Author: Jaume Bonet <bonet>
# @Date:   29-Jun-2017
# @Email:  jaume.bonet@gmail.com
# @Filename: designFrame.py
# @Last modified by:   bonet
# @Last modified time: 28-Mar-2018

# Standard Libraries
import itertools

# External Libraries
import pandas as pd
from pandas.core.index import Index
from pandas.core.dtypes.common import is_scalar
import numpy as np

# This Library
from .rsbase import RSBaseDesign
from .apply import frame_apply
import rstoolbox.analysis as ra


__all__ = ["DesignSeries", "DesignFrame"]


class DesignSeries( pd.Series, RSBaseDesign ):
    """
    The :py:class:`.DesignSeries` extends the :py:class:`~pandas.Series`
    adding some functionalities in order to improve its usability in
    the analysis of a single design decoys.

    It is generated as the *reduced dimensionality version* of the
    :py:class:`.DesignFrame`.
    """

    _metadata = ['_reference']
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
        return DesignFrame

    def __finalize__(self, other, method=None, **kwargs):
        if method == "inherit":
            # Avoid columns from DesignFrame to become DesignSeries
            if not isinstance(self.name, int):
                return pd.Series(self)

        for name in self._metadata:
            setattr(self, name, getattr(other, name, {}))
        return self


class DesignFrame( pd.DataFrame, RSBaseDesign ):
    """
    The :py:class:`.DesignFrame` extends the :py:class:`~pandas.DataFrame`
    adding some functionalities in order to improve its usability in
    the analysis of sets of design decoys.

    Filled through the functions provided through this library, each
    row represents a decoy while each column represents the scores
    attached to it.

    As a rule, it is assumed that the object:

    #. has a column named **description** \
    that stores the identifier of the corresponding decoy.
    #. holds sequences (a design decoy might be composed of \
    multiple chains) in columns named **sequnece_[seqID]**.

    This two assumptions are easily adapted if casting a
    :py:class:`~pandas.DataFrame` into the class, and several functions
    of the library depend on them.

    The :py:class:`.DesignFrame` basically contains three extra attributes
    (accessible through the appropiate functions):

    #. **reference_sequence:** A reference sequence can be added for each ``seqID`` \
    present in the :py:class:`.DesignFrame`. By adding this sequence, other functions \
    of the library can add that information to its calculations.
    #. **reference_shift:** A reference shift can be added for each ``seqID`` \
    present in the :py:class:`.DesignFrame`. In short, this would be the initial number \
    of the protein in the source PDB. This allows working with the right numbering. This \
    value is, by default, 1 in all ``seqID``.
    #. **source_files:** The object stores the source files from which it has been loaded \
    (as long as it is loaded with :py:func:`.parse_rosetta_file`). This information can \
    be used to extract the structures from the silent files.
    """
    _metadata = ['_reference', '_source_files']
    _subtyp = 'design_frame'

    def __init__(self, *args, **kwargs):
        reference = kwargs.pop('reference', {})
        source    = kwargs.pop('source', set())
        super(DesignFrame, self).__init__(*args, **kwargs)
        self._reference = reference
        self._source_files = source

    def add_source_file( self, file ):
        """
        Adds a source file to the :py:class:`.DesignFrame`. This can be used to know where to
        extract the structure from if needed.

        :param str file: Name of the file to add.
        """
        self._source_files.add( file )

    def replace_source_files( self, files ):
        """
        Replaces source files of the :py:class:`.DesignFrame`. These can be used to know where to
        extract the structure from if needed.

        :param files: List of names of the files to add.
        """
        self._source_files = set( files )

    def add_source_files( self, files ):
        """
        Adds source files to the :py:class:`.DesignFrame`. These can be used to know where to
        extract the structure from if needed.

        :param files: List of names of the files to add.
        """
        self._source_files = self._source_files.union( files )

    def get_source_files( self ):
        """
        Returns the set of source files linked to this :py:class:`.DesignFrame`.

        :return: Set of files.
        """
        return self._source_files

    def has_source_files( self ):
        """
        Checks if there are source files added.

        :return: bool
        """
        return bool(self._source_files)

    def get_sequence_with( self, seqID, selection, confidence=1 ):
        """
        Selects those decoys with a particular set of residue matches.
        Basically, is meant to find, for example, all the decoys in which
        position 25 is A and position 46 is T.

        :param str seqID: Identifier of the sequence of interest.
        :param str selection: List of tuples with position and residue
            type (in 1 letter code)
        :param float confidence: Percentage of the number of the selection
            rules that we expect the matches to fulfill. Default is 1 (all).
        :return: Filtered :py:class:`.DesignFrame`
        """
        def match_residues( value, seqID, confidence ):
            t = 0
            for s in selection:
                t += 1 if value[s[0] - 1] == s[1] else 0
            return t / float(len(selection)) >= float(confidence)

        return self.loc[ self.apply(
            lambda row: match_residues(row["sequence_{0}".format(seqID)], selection, confidence ),
            axis=1
        )]

    def sequence_distance( self, seqID ):
        """
        Generate a matrix counting the distance between each pair of sequences in the
        :py:class:`.designFrame`. This is a time-consuming operation; better to execute it over a
        specific set of selected decoys that over all your designs.

        :param seqID: Identifier of the sequence of interest.
        :type seqID: :py:class:`str`

        return: :py:class:`~pandas.DataFrame`.

        :raises:
            :KeyError: if ``seqID`` cannot be found.
            :KeyError: if ``description`` column cannot be found.
        """
        # https://stackoverflow.com/questions/49235538/pandas-series-compare-values-all-vs-all
        a = np.array(list(map(list, self.get_sequence(seqID).values)))
        ids = self.get_id().values

        n = len(a)
        d = {(i, j): np.sum(a[i] != a[j]) for i in range(n) for j in range(n) if j > i}

        res = np.zeros((n, n))
        keys = list(zip(*d.keys()))

        res[keys[0], keys[1]] = list(d.values())
        res += res.T

        return pd.DataFrame(res, columns=ids, index=ids, dtype=int)

    def sequence_frequencies( self, seqID, seqType="protein", cleanExtra=True, cleanUnused=-1 ):
        """
        Generates a :py:class:`.SequenceFrame` for the frequencies of
        the sequences in the __DesignFrame__ with seqID identifier.
        If there is a reference_sequence for this seqID, it will also
        be attached to the __SequenceFrame__.
        All letters in the sequence will be capitalized. All symbols that
        do not belong to string.ascii_uppercase will be transformed to "*"
        as this is the symbol recognized by the substitution matrices.

        :param str seqID: Identifier of the sequence sets of interest.
        :param str seqType: Type of sequence: protein, dna, rna.
        :param bool cleanExtra: Remove from the SequenceFrame the non-regular
            amino/nucleic acids if they are empty for all positions.
        :param int cleanUnused: Remove from the SequenceFrame the regular
            amino/nucleic acids if they frequency is equal or under the value . Default is -1,
            so nothing is deleted.
        :return: :py:class:`.SequenceFrame`
        """
        return ra.sequential_frequencies(self, seqID, "sequence", seqType, cleanExtra, cleanUnused)

    def sequence_bits( self, seqID, seqType="protein", cleanExtra=True, cleanUnused=False ):
        """
        Generates a :py:class:`.SequenceFrame` for the bits of
        the sequences in the __DesignFrame__ with seqID identifier.
        If there is a reference_sequence for this seqID, it will also
        be attached to the __SequenceFrame__.
        Bit calculation is performed as explained in http://www.genome.org/cgi/doi/10.1101/gr.849004
        such as:

        Rseq = Smax - Sobs = log2 N - (-sum(n=1,N):pn * log2 pn)

        Where:
            - N is the total number of options (4: DNA/RNA; 20: PROTEIN).
            - pn is the observed frequency of the symbol n.

        :param str seqID: Identifier of the sequence sets of interest.
        :param str seqType: Type of sequence: protein, dna, rna.
        :param int shift: Numbering assign to the first residue of the chain. Default is None, pick
            shift from the reference sequence
        :param bool cleanExtra: Remove from the SequenceFrame the non-regular
        amino/nucleic acids if they are empty for all positions.
        :param bool cleanUnused: Remove from the SequenceFrame the regular
        amino/nucleic acids if they are empty for all positions

        :return: :py:class:`.SequenceFrame`
        """
        df = self.sequence_frequencies(seqID, seqType, cleanExtra, cleanUnused)
        return df.to_bits()

    def structure_frequencies( self, seqID, seqType="protein", cleanExtra=True, cleanUnused=-1 ):
        """
        Generates a :py:class:`.SequenceFrame` for the frequencies of
        the secondary structure in the __DesignFrame__ with seqID identifier.
        If there is a reference_structure for this seqID, it will also
        be attached to the __SequenceFrame__.
        All letters in the secondary structure will be capitalized. All symbols that
        do not belong to string.ascii_uppercase will be transformed to "*"
        as this is the symbol recognized by the substitution matrices.

        :param str seqID: Identifier of the secondary structure sets of interest.
        :param str seqType: Type of sequence: protein.
        :param bool cleanExtra: Remove from the SequenceFrame the non-regular
            amino/nucleic acids if they are empty for all positions.
        :param int cleanUnused: Remove from the SequenceFrame the regular
            amino/nucleic acids if they frequency is equal or under the value . Default is -1,
            so nothing is deleted.
        :return: :py:class:`.SequenceFrame`
        """
        seqType = seqType + "_sse"
        return ra.sequential_frequencies(self, seqID, "structure", seqType, cleanExtra, cleanUnused)

    def structure_bits( self, seqID, seqType="protein", cleanExtra=True, cleanUnused=False ):
        """
        Generates a :py:class:`.SequenceFrame` for the bits of
        the secondary structure in the __DesignFrame__ with seqID identifier.
        If there is a reference_structure for this seqID, it will also
        be attached to the __SequenceFrame__.
        Bit calculation is performed as explained in http://www.genome.org/cgi/doi/10.1101/gr.849004
        such as:

        Rseq = Smax - Sobs = log2 N - (-sum(n=1,N):pn * log2 pn)

        Where:
            - N is the total number of options (4: DNA/RNA; 20: PROTEIN).
            - pn is the observed frequency of the symbol n.

        :param str seqID: Identifier of the sequence sets of interest.
        :param str seqType: Type of sequence: protein, dna, rna.
        :param int shift: Numbering assign to the first residue of the chain. Default is None, pick
            shift from the reference sequence
        :param bool cleanExtra: Remove from the SequenceFrame the non-regular
        amino/nucleic acids if they are empty for all positions.
        :param bool cleanUnused: Remove from the SequenceFrame the regular
        amino/nucleic acids if they are empty for all positions

        :return: :py:class:`.SequenceFrame`
        """
        df = self.structure_frequencies(seqID, seqType, cleanExtra, cleanUnused)
        return df.to_bits()

    def _metadata_defaults(self, name):
        if name == "_source_files":
            return set()
        if name == "_reference":
            return {}
        return None

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

    def _box_col_values(self, values, items):
        """ provide boxed values for a column """
        return self._constructor_sliced(values, index=self.index,
                                        name=items, fastpath=True)

    def _ixs(self, i, axis=0):
        """
        i : int, slice, or sequence of integers
        axis : int
        """

        # irow
        if axis == 0:
            """
            Notes
            -----
            If slice passed, the resulting data will be a view
            """

            if isinstance(i, slice):
                return self[i]
            else:
                label = self.index[i]
                if isinstance(label, Index):
                    # a location index by definition
                    result = self.take(i, axis=axis)
                    copy = True
                else:
                    new_values = self._data.fast_xs(i)
                    if is_scalar(new_values):
                        return new_values

                    # if we are a copy, mark as such
                    copy = (isinstance(new_values, np.ndarray) and
                            new_values.base is None)
                    result = self._constructor_sliced(new_values,
                                                      index=self.columns,
                                                      name=self.index[i],
                                                      dtype=new_values.dtype)
                result._set_is_copy(self, copy=copy)
                return result

        # icol
        else:
            """
            Notes
            -----
            If slice passed, the resulting data will be a view
            """

            label = self.columns[i]
            if isinstance(i, slice):
                # need to return view
                lab_slice = slice(label[0], label[-1])
                return self.loc[:, lab_slice]
            else:
                if isinstance(label, Index):
                    return self._take(i, axis=1, convert=True)

                index_len = len(self.index)

                # if the values returned are not the same length
                # as the index (iow a not found value), iget returns
                # a 0-len ndarray. This is effectively catching
                # a numpy error (as numpy should really raise)
                values = self._data.iget(i)

                if index_len and not len(values):
                    values = np.array([np.nan] * index_len, dtype=object)
                result = self._constructor_sliced(
                    values, index=self.index, name=label, fastpath=True)

                # this is a cached value, mark it so
                result._set_as_cached(label, self)

                return result

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
                setattr(self, name, getattr(other.left, name, self._metadata_defaults(name)))
        else:
            for name in self._metadata:
                setattr(self, name, getattr(other, name, self._metadata_defaults(name)))
        return self

    # this is a fix until the issue is solved in pandas
    def apply(self, func, axis=0, broadcast=None, raw=False, reduce=None,
              result_type=None, args=(), **kwds):
            op = frame_apply(self,
                             func=func,
                             axis=axis,
                             broadcast=broadcast,
                             raw=raw,
                             reduce=reduce,
                             result_type=result_type,
                             args=args,
                             kwds=kwds)
            return op.get_result()
