# @Author: Jaume Bonet <bonet>
# @Date:   29-Jun-2017
# @Email:  jaume.bonet@gmail.com
# @Filename: designFrame.py
# @Last modified by:   bonet
# @Last modified time: 12-Feb-2018

# Standard Libraries
import itertools
import os
import re
import string

# External Libraries
import pandas as pd
import numpy as np

# This Library
from .sequenceFrame import SequenceFrame
from rstoolbox.utils import add_column


class DesignFrame( pd.DataFrame ):
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
    _metadata = ['_reference_sequence', '_source_files']

    def __init__(self, *args, **kw):
        super(DesignFrame, self).__init__(*args, **kw)
        self._reference_sequence = self._metadata_defaults("_reference_sequence")
        self._source_files = self._metadata_defaults("_source_files")

    def reference_sequence( self, seqID, sequence=None, shift=1 ):
        """
        Setter/Getter for a reference sequence attached to a particular
        sequence ID. It also allows to provide the shift of the sequence
        count with respect to a linear numbering (Rosetta numbering).

        :param str seqID: Identifier of the reference sequence
        :param str sequence: Reference sequence. By default is
            :py:data:`None`, which turns the function into a getter.
        :param int shift: In case the sequence does not start in 1, how much
            do we need to shift? Basically provide the number of the first
            residue of the chain. Default is 1.
        :return: str
        :raise KeyError: If seqID does not exist.
        """
        if sequence is not None:
            self._reference_sequence.setdefault(seqID, {"seq": sequence, "sft": shift})
        else:
            if seqID not in self._reference_sequence:
                raise KeyError("There is no reference sequence with ID: {}\n".format(seqID))
        return self._reference_sequence[seqID]["seq"]

    def reference_shift( self, seqID, shift=None ):
        """
        Setter/Getter for a reference shift attached to a particular
        sequence ID. If there is no reference sequence attached, it returns
        1, as "no reference shift".

        :param str seqID: Identifier of the reference sequence
        :param int shift: In case the sequence does not start in 1, how much
            do we need to shift? Basically provide the number of the first
            residue of the chain. By default is :py:data:`None`, which turns
            the function into a getter.
        :return: int
        """
        if seqID in self._reference_sequence:
            if shift is not None:
                self._reference_sequence[seqID]["sft"] = shift
            return self._reference_sequence[seqID]["sft"]
        else:
            return 1

    def has_reference_sequence( self, seqID ):
        """
        Checks if there is a reference sequence for the provided
        sequence ID.

        :param str seqID: Identifier of the reference sequence
        :return: bool
        """
        return seqID in self._reference_sequence

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
                t += 1 if value[s[0]-1] == s[1] else 0
            return t/float(len(selection)) >= float(confidence)

        return self.loc[ self.apply(
                    lambda row: match_residues(row["sequence_{0}".format(seqID)], selection, confidence ),
                    axis=1
                )]

    def identify_mutants( self, seqID, shift=None ):
        """
        Checks the sequence in column sequence_<seqID> againts the reference_sequence.
        Adds to the :py:class:`.designFrame` two new columns: mutants_<seqID>, which lists
        the mutations of the particular decoy vs. the reference_sequence, mutant_positions_<seqID>
        just with those same positions and mutant_count_<seqID> with the count of the number of
        mutations. Reference and design sequence must be of the same length.

        :param str seqID: Identifier of the sequence of interest.
        :param int shift: Numbering assign to the first residue of the chain. Default is None, pick
            shift from the reference sequence
        :return: Changes :py:class:`.designFrame` and returns it
        """
        def mutations( reference, sequence, shift=1 ):
            data = []
            assert len(reference) == len(sequence)
            for i in range(len(reference)):
                if reference[i].upper() != sequence[i].upper():
                    data.append(reference[i].upper() + str(i + shift) + sequence[i].upper())
            return ",".join(data)
        def count_muts( mutations ):
            if mutations == "": return 0
            return len(mutations.split(","))

        this_shift = shift if shift is not None else self.reference_shift(seqID)
        self["mutants_{0}".format(seqID)] = self.apply(
            lambda row: mutations(self.reference_sequence(seqID), row["sequence_{0}".format(seqID)], this_shift),
            axis=1 )
        self["mutant_positions_{0}".format(seqID)] = self["mutants_{0}".format(seqID)].str.replace(r"[a-zA-Z]","")
        self["mutant_count_{0}".format(seqID)]     = self.apply(lambda row: count_muts(row["mutants_{0}".format(seqID)]), axis=1)

        return self

    def sequence_distance( self, seqID ):
        """
        Generate a matrix counting the distance between each pair of sequences in the :py:class:`.designFrame`.

        :param str seqID: Identifier of the sequence of interest.
        :type seqID: :py:class:`str`

        return: :py:class:`~pandas.DataFrame`

        :raises:
            :KeyError: if ``seqID`` cannot be found.
        """
        def count_differences( sequence, df ):
            return df.apply(lambda x : sum(1 for i, j in zip(x["sequence_{}".format(seqID)], sequence) if i != j), axis=1)

        try:
            df = self.apply(lambda x: count_differences( x["sequence_{}".format(seqID)], self), axis=1)
            return df.rename(self["description"], axis="columns").rename(self["description"], axis="rows")
        except KeyError:
            raise KeyError("There is no seqID {} to be found".format(seqID))

    def generate_mutant_variants( self, seqID, mutations ):
        """
        Expands the selected sequences by ensuring that all the provided mutant combinations.
        Thus, something such as::

            df.generate_mutant_variants("A", [(20, "AIV"), (31, "EDQR")])

        Will generate all the variant sequences of "A" that combine the expected mutations in
        position 20 and 31 (according to the reference shift). That would be a total of 3*4
        combinations plus the original sequence.

        Alters the names of the designs in **description**.

        :param seqID: Identifier of the sequence sets of interest.
        :type seqID: :py:class:`str`
        :param mutations: List of mutations to generate in a format (position, variants)
        :type mutations: :py:class:`list`[(:py:class:`int`, :py:class:`str`),..]

        :return: :py:class:`.DesignFrame`
        """
        def multiplex( row, seqID, muts, refshift):
            seqNM = "sequence_{}".format(seqID)
            seq   = list(row[seqNM])
            for p in muts:
                seq[p[0]-1] = p[1]
            data = {seqNM: ["".join(x) for x in itertools.product(*seq)]}
            data[seqNM].insert(0, row[seqNM])
            data["description"] = [row["description"] + "_{0:04d}".format(x) for x in range(len(data[seqNM]))]
            for col in row.index:
                if col not in [seqNM, "description"]:
                    data[col] = [row[col]] * len(data["description"])
            df = DesignFrame(data)
            df.reference_sequence(seqID, row[seqNM], refshift)
            df = df.identify_mutants(seqID)
            del(df._reference_sequence[seqID])
            return df

        designs = []
        refshift = self.reference_shift(seqID)
        for i, row in self.iterrows():
            designs.append(multiplex(row, seqID, mutations, refshift))
        df = pd.concat(designs)
        if self.has_reference_sequence(seqID):
            df.reference_sequence(seqID, self.reference_sequence(seqID), self.reference_shift(seqID))
        df.drop_duplicates(["sequence_{}".format(seqID)], inplace=True)
        return df.reset_index(drop=True)

    def sequence_frequencies( self, seqID, seqType="protein", shift=None, cleanExtra=True, cleanUnused=-1 ):
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
        :param int shift: Numbering assign to the first residue of the chain. Default is None, pick
            shift from the reference sequence
        :param bool cleanExtra: Remove from the SequenceFrame the non-regular
            amino/nucleic acids if they are empty for all positions.
        :param int cleanUnused: Remove from the SequenceFrame the regular
            amino/nucleic acids if they frequency is equal or under the value . Default is -1,
            so nothing is deleted.
        :return: :py:class:`.SequenceFrame`
        """
        sserie = self["sequence_{0}".format(seqID)].values
        table, extra = self._get_sequence_table( seqType )
        for x in range(len(sserie[0])):
            for k in table:
                table[k].append(float(0))
            for y in range(len(sserie)):
                aa = sserie[y][x].upper()
                if aa in string.whitespace or aa in string.punctuation:
                    aa = "*"
                table[aa][-1] += float(1)
        for k in table:
            for x in range(len(table[k])):
                if table[k][x] != 0:
                    table[k][x] /= float(len(sserie))

        df = SequenceFrame(table)
        df.measure("frequency")
        df.extras( extra )
        if self.has_reference_sequence(seqID):
            df.reference_sequence(self.reference_sequence(seqID), self.reference_shift(seqID))
        df.delete_extra( cleanExtra )
        df.delete_empty( cleanUnused )
        df.clean()
        df.index = df.index + (shift if shift is not None else self.reference_shift(seqID))
        return df

    def sequence_bits( self, seqID, seqType="protein", shift=None, cleanExtra=True, cleanUnused=False ):
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
        df = self.sequence_frequencies(seqID, seqType, shift, cleanExtra, cleanUnused)
        return df.to_bits()

    def _get_sequence_table( self, seqType ):
        """
        Generates the table to fill sequence data in order to create
        a :py:class:`.SequenceFrame`

        :param str seqType: Type of sequence: protein, dna, rna.
        :return: dict
        :raise ValueError: If seqType is not known
        """
        table = {}
        extra = []
        if seqType.lower() == "protein":
            # X = UNKNOWN  # * = GAP
            # B = N or D   # Z = E or Q
            table = {
                'C' : [], 'D' : [], 'S' : [], 'Q' : [], 'K' : [],
                'I' : [], 'P' : [], 'T' : [], 'F' : [], 'N' : [],
                'G' : [], 'H' : [], 'L' : [], 'R' : [], 'W' : [],
                'A' : [], 'V' : [], 'E' : [], 'Y' : [], 'M' : [],
                'X' : [], '*' : [], 'B' : [], 'Z' : []
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
                'C' : [], 'A' : [], 'T' : [], 'G' : [], 'X' : [], '*' : [],
                'B' : [], 'D' : [], 'H' : [], 'K' : [], 'M' : [], 'N' : [],
                'R' : [], 'S' : [], 'V' : [], 'W' : [], 'Y' : []
            }
            if seqType.lower() == "rna":
                table.setdefault( 'U' , [])
                table.pop('T', None)
            extra = ['X', '*', 'B', 'D', 'H', 'K', 'M', 'N', 'R', 'S', 'V', 'W', 'Y']
        else:
            raise ValueError("sequence type {0} unknown".format(seqType))
        return table, extra

    def _metadata_defaults(self, name):
        if name == "_source_files":
            return set()
        if name == "_reference_sequence":
            return {}
        return None
    #
    # Implement pandas methods
    #

    @property
    def _constructor(self):
        return DesignFrame

    def __finalize__(self, other, method=None, **kwargs):
        """propagate metadata from other to self """
        strict  = kwargs["strict"] if "strict" in kwargs else True
        # concat operation:
        #   (1) accumulate _source_files from all the concatenated DesignFrames (if there is any)
        #   (2) check reference_sequence and reference_shift, merge and keep only if they are the same.
        if method == 'concat':
            source_files = set()
            refseqs = []
            reference_sequence = {}
            for i, o in enumerate(other.objs):
                source_files.update(getattr(o, "_source_files", set()))
                refseqs.append(getattr(o, "_reference_sequence", {}))
            # _source_files
            setattr(self, "_source_files", source_files)
            # _reference_sequence
            ids = list(set(itertools.chain.from_iterable([x.keys() for x in refseqs])))
            for r in refseqs:
                for i in ids:
                    if i in r:
                        if not i in reference_sequence:
                            reference_sequence.setdefault(i, r[i])
                        else:
                            if r[i] != reference_sequence[i]:
                                raise ValueError("Concatenating designFrames with different ref sequence for the same seqID.")
            setattr(self, "_reference_sequence", reference_sequence)
        # merge operation:
        # Keep metadata of the left object.
        elif method == 'merge':
            for name in self._metadata:
                setattr(self, name, getattr(other.left, name, self._metadata_defaults(name)))
        else:
            for name in self._metadata:
                setattr(self, name, getattr(other, name, self._metadata_defaults(name)))
        return self
