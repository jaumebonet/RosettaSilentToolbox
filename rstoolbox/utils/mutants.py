# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: get_identified_mutants
.. func:: get_mutations
.. func:: get_mutation_positions
.. func:: get_mutation_count
.. func:: identify_mutants
.. func:: generate_mutant_variants
.. func:: generate_mutants_from_matrix
.. func:: generate_wt_reversions
.. func:: score_by_pssm
.. func:: make_resfile
.. func:: view_mutants_alignment
"""
# Standard Libraries
import os
import sys
import itertools
import re
import tempfile

# External Libraries
import pandas as pd
import numpy as np

# This Library
from .getters import _check_type, _get_available, _check_column
import rstoolbox.core as core
from rstoolbox.utils import make_rosetta_app_path, execute_process


def get_identified_mutants( self ):
    """List for which sequence identifiers mutants have been calculated.

    This will generate a list with all the available identifiers in no particular
    order.

    :return: :func:`list` of :class:`str`

    :raises:
        :TypeError: If the data container is not :class:`~pandas.DataFrame`
            or :class:`~pandas.Series`

    .. seealso::
        :meth:`.DesignFrame.identify_mutants`
        :meth:`.DesignSeries.identify_mutants`
    """
    _check_type(self)
    return _get_available(self, "mutants_")


def get_mutations( self, seqID ):
    """Return the mutantions data for ``seqID`` available in the container.

    For each available decoy, this will be a coma-separated string of mutant
    positions with the code ``<original><position><mutation>``.

    :param str seqID: |seqID_param|.

    :return: :class:`str` or :class:`~pandas.Series` - depending on the input

    :raises:
        :TypeError: |indf_error|.
        :KeyError: |mutID_error|.

    .. seealso::
        :meth:`.DesignFrame.identify_mutants`
        :meth:`.DesignSeries.identify_mutants`
    """
    _check_type(self)
    return self[_check_column(self, "mutants", seqID)]


def get_mutation_positions( self, seqID ):
    """Return the mutantion positions data for `seqID` available in the container.

    For each available decoy, this will be a coma-separated string with the sequence
    position of the mutated residues.

    :param str seqID: |seqID_param|.

    :return: :class:`str` or :class:`~pandas.Series` - depending on the input

    :raises:
        :TypeError: |indf_error|.
        :KeyError: |mutID_error|.

    .. seealso::
        :meth:`.DesignFrame.identify_mutants`
        :meth:`.DesignSeries.identify_mutants`
    """
    _check_type(self)
    return self[_check_column(self, "mutant_positions", seqID)]


def get_mutation_count( self, seqID ):
    """ Return the number of mutantion positions data for `seqID` available in the container.

    Basically this contains the total number of mutations in a particular decoy with respect
    to its ``reference_sequence``.

    :param str seqID: |seqID_param|.

    :return: :class:`int` or :class:`~pandas.Series` - depending on the input

    :raises:
        :TypeError: |indf_error|.
        :KeyError: |mutID_error|.

    .. seealso::
        :meth:`.DesignFrame.identify_mutants`
        :meth:`.DesignSeries.identify_mutants`
    """
    _check_type(self)
    return self[_check_column(self, "mutant_count", seqID)]


def identify_mutants( self, seqID ):
    """Assess mutations of each decoy for sequence ``seqID`` againt the ``reference_sequence``.

    Adds to the container two new columns:

    ============================  ===========================================================
    Column                                                Data Content
    ============================  ===========================================================
    **mutants_<seqID>**           Lists the **mutations** of the particular decoy
    **mutant_positions_<seqID>**  Lists the **positions** of mutations in the particular decoy
    **mutant_count_<seqID>**      **Count** of the number of mutations
    ============================  ===========================================================

    .. tip::
        ``reference_sequence`` and design sequence must be of the same length. If that is **not**
        the case, it could be solved with the use of a non ``string.ascii_uppercase`` character like
        `"*"`.

    :param str seqID: |seqID_param|.

    :return: Union[:class:`.DesignSeries`, :class:`.DesignFrame`] -
        a copy of the data container with the new columns.

    :raise:
        :ValueError: If length of ``reference_sequence`` and decoy are not the same.

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz",
           ...:                         {'scores': ['score'], 'sequence': 'B'})
           ...: df.add_reference_sequence('B', df.get_sequence('B').values[0])
           ...: df.iloc[1:].identify_mutants('B')
    """
    from rstoolbox.components import get_selection

    def mutations( reference, row, seqID ):
        data = []
        datn = []
        sequence = row.get_sequence(seqID)
        if len(reference) != len(sequence):
            raise ValueError("Sequence lengths do not match")
        for i, refi in enumerate(reference):
            if refi.upper() != sequence[i].upper():
                shift = get_selection(i + 1, seqID, row.get_reference_shift(seqID))[0]
                data.append(refi.upper() + str(shift) + sequence[i].upper())
                datn.append(str(shift))
        return ",".join(data), ",".join(datn), len(data)

    refseq  = self.get_reference_sequence(seqID)
    mutants = "mutants_{0}".format(seqID)
    mposits = "mutant_positions_{0}".format(seqID)
    mcounts = "mutant_count_{0}".format(seqID)
    df = self.copy()
    if isinstance(self, pd.DataFrame):
        df[[mutants, mposits, mcounts]] = df.apply(
            lambda row: mutations(refseq, row, seqID),
            axis=1, result_type="expand" )
    elif isinstance(df, pd.Series):
        a, b, c = mutations(refseq, df, seqID)
        df[mutants], df[mposits], df[mcounts] = a, b, c

    return df


def generate_mutant_variants( self, seqID, mutations, keep_scores=False ):
    """Expands selected decoy sequences generating all the provided mutant combinations.

    For all the new mutations provided, it will generate all the possible combinations with
    those mutations and annotate them with respect to the ``reference_sequence``.

    A mutation will be specified as a :class:`tuple` of ``length=2``. The first position will
    be the sequence position to target (``reference_shift`` aware) and the second will be a
    string with all the desired residue types. Multiple positions can be provided in a
    :func:`list`::

        mutants = [(20, "AIV"), (31, "EDQR")]

    Lastly, when multiple changes are provided for a position, this will translate into an
    **insertion**.

    .. tip::
        The number of positions and mutations for position produce an exponential increment
        of the generated sequences. Thus, the previous example will generate ``3 * 4`` new
        sequences. Depending on the input this can explode pretty fast, be aware.

    .. tip::
        ``*`` will call all 20 regular amino acids for a given position.

    Alters the names of the designs in **description** by adding a ``_v<number>`` suffix.

    By providing multiple input decoys, sequence can be repeated. Thus, repeated sequences
    will be filtered; the provided copy will be the first instance of the sequence.

    :param str seqID: |seqID_param|.
    :param mutations: List of mutations to generate in a format (position, variants)
    :type mutations: :func:`list` of :class:`tuple` (:class:`int`, :class:`str`)
    :param bool keep_scores: New variants inherit scores from their source sequence.
        This is **not recommended**, as it can get confusing (Default: :data:`False`).

    :return: :class:`.DesignFrame`

    .. seealso::
        :meth:`.DesignFrame.generate_mutants_from_matrix`
        :meth:`.DesignFrame.make_resfile`
        :meth:`.DesignSeries.generate_mutants_from_matrix`
        :meth:`.DesignSeries.make_resfile`

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz",
           ...:                         {'scores': ['score', 'description'], 'sequence': 'B'})
           ...: df.add_reference_sequence('B', df.get_sequence('B').values[0])
           ...: mutants = [(20, "AIV"), (31, "EDQR")]
           ...: df.iloc[1].generate_mutant_variants('B', mutants)
    """
    from rstoolbox.components import get_selection

    def multiplex( row, seqID, muts ):
        seqNM = _check_column(row, "sequence", seqID)
        idNM = "description"
        seq   = list(row[seqNM])
        for p in reversed(muts):
            # -1 because we are going to access string positions.
            shift = get_selection(p[0], seqID, row.get_reference_shift(seqID))[0] - 1
            seq[shift] = p[1] if p[1] != "*" else "ARNDCQEGHILKMFPSTWYV"
        data = {seqNM: ["".join(x) for x in itertools.product(*seq)]}
        data[seqNM].insert(0, row[seqNM])
        name = row.get_id()
        if not bool(re.search("_v\d+$", row.get_id())):
            data[idNM] = [name + "_v{0:04d}".format(x) for x in range(len(data[seqNM]))]
            data[idNM][0] = name
        else:
            data[idNM] = [name + "_v{0:04d}".format(x) for x in range(1, len(data[seqNM]) + 1)]
        if keep_scores:
            for col in row.index:
                if col not in [seqNM, idNM]:
                    data[col] = [row[col]] * len(data[idNM])
        else:
            for seq in row.get_available_sequences():
                if seq != seqID:
                    data[_check_column(row, "sequence", seq)] = row.get_sequence(seq)
        df = row._constructor_expanddim(data)
        return df

    mutations.sort(key=lambda tup: tup[0])
    if isinstance(self, pd.DataFrame):
        designs = []
        for _, row in self.iterrows():
            designs.append(multiplex(row, seqID, mutations))
        df = pd.concat(designs)
    elif isinstance(self, pd.Series):
        df = multiplex(self, seqID, mutations)
    else:
        raise NotImplementedError

    df.transfer_reference(self)

    avail_seqs = df.get_available_sequences()
    avail_refs = []
    seqs = [_check_column(df, "sequence", seq) for seq in avail_seqs]
    df.drop_duplicates(seqs, inplace=True)
    for seq in avail_seqs:
        if df.has_reference_sequence(seq):
            df = df.identify_mutants(seq)
            avail_refs.append(seq)
    if len(avail_refs) > 1:
        muts = ["mutant_count_{}".format(seq) for seq in avail_refs]
        df["mutant_count_all"] = df[muts].sum(axis=1)

    return df.reset_index(drop=True)


def generate_mutants_from_matrix( self, seqID, matrix, count,
                                  key_residues=None, limit_refseq=False ):
    """From a provided positional frequency matrix, generates ``count`` random variants.

    It takes into account the individual frequency assigned to each residue type and
    position. It does **not** generate the highest possible scored sequence according to
    the matrix, but picks randomly at each position according to the frequencies in for
    that position.

    For each :class:`.DesignSeries`, it will generate a :class:`.DesignFrame` in which the
    original sequence becomes the ``reference_sequence``, inheriting the ``reference_shift``.

    .. warning::
        This is a **computationaly expensive** function. Take this in consideration when trying
        to run it.

    Each :class:`.DesignFrame` will have the following structure:

    ======================  ============================================
    Column                                                Data Content
    ======================  ============================================
    **description**         Identifier fo the mutant
    **sequence_<seqID>**    Sequence content
    **pssm_score_<seqID>**  Score obtained by applying ``matrix``
    ======================  ============================================

    :param str seqID: |seqID_param|
    :param matrix: Positional frequency matrix. **column:** residue type; **index:**
        sequence position.
    :type matrix: :class:`~pandas.DataFrame`
    :param int count: Expected number of **unique** generated combinations. If the number is
        bigger than the possible options, it will default to the total amount of options.
    :param key_residues: |keyres_param|.
    :type key_residues: |keyres_types|
    :param bool limit_refseq: When :data:`True`, pick only residue types with probabilities
        equal or higher to the source sequence.

    :return: :func:`list` of :class:`.DesignFrame` - New set of design sequences.

    :raises:
        :ValueError: if matrix rows do not match sequence length.

    .. seealso::
        :meth:`.DesignFrame.generate_mutant_variants`
        :meth:`.DesignFrame.score_by_pssm`
        :meth:`.DesignSeries.generate_mutant_variants`
        :meth:`.DesignSeries.score_by_pssm`

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: from rstoolbox.tests.helper import random_frequency_matrix
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz",
           ...:                         {'scores': ['score', 'description'], 'sequence': 'B'})
           ...: df.add_reference_sequence('B', df.get_sequence('B').values[0])
           ...: matrix = random_frequency_matrix(len(df.get_reference_sequence('B')), 0)
           ...: key_res = [3,5,8,12,15,19,25,27]
           ...: mutants = df.iloc[1].generate_mutants_from_matrix('B', matrix, 5, key_res)
           ...: mutants[0].identify_mutants('B')

    """
    from rstoolbox.components import get_selection
    from rstoolbox.components import DesignSeries, DesignFrame

    def max_options( matrix, seq, key_residues, limit_refseq):
        if limit_refseq is False:
            return np.power(20, len(key_residues))
        else:
            ori_index = matrix.index
            matrix = matrix.copy()
            matrix.index = range(0, matrix.shape[0])
            options = (matrix.apply(lambda row: np.sum(row >= row[seq[row.name]]), axis=1))
            options.index = ori_index
            return np.prod(options[key_residues])

    data = []
    if isinstance(self, pd.DataFrame):
        for _, row in self.iterrows():
            data.extend(row.generate_mutants_from_matrix(seqID, matrix, count,
                                                         key_residues, limit_refseq))
        return data

    if matrix.shape[0] != len(self.get_sequence(seqID)):
        raise ValueError("Matrix rows and sequence length should match.")
    # Make sure index and sequence shift match
    matrix = matrix.copy()
    shift = self.get_reference_shift(seqID)
    matrix.index = get_selection(None, seqID, shift, length=matrix.shape[0])

    if key_residues is not None:
        key_residues = get_selection(key_residues, seqID, shift, matrix.shape[0])
    else:
        key_residues = list(matrix.index.values)

    seqnm = "sequence_{}".format(seqID)
    data.append(DesignFrame([], columns=["description", seqnm]))
    name  = self.get_id()

    options = max_options(matrix, self.get_sequence(seqID), key_residues, limit_refseq)
    # some numbers are just too big for python...
    if options <= 0:
        options = count + 1

    while data[-1].shape[0] < min(count, options):
        seqaa = list(self.get_sequence(seqID))
        thisname = name + "_v{0:04d}".format(data[-1].shape[0] + 1)
        for aap in key_residues:
            matI = matrix.loc[aap].copy()
            if limit_refseq:
                matI[matI < matI[seqaa[aap - 1]]] = 0
                matI = matI / matI.sum()
            seqaa[aap - 1] = np.random.choice(matI.index.values, 1, p=list(matI))[0]
        if "".join(seqaa) == self.get_sequence(seqID):
            continue
        data[-1] = data[-1].append(DesignSeries([thisname, "".join(seqaa)],
                                                ["description", seqnm]),
                                   ignore_index=True)
        data[-1].drop_duplicates([seqnm])
    data[-1].add_reference(seqID, self.get_sequence(seqID), shift=self.get_reference_shift(seqID))
    data[-1] = data[-1].score_by_pssm(seqID, matrix)
    return data


def generate_wt_reversions( self, seqID, key_residues=None ):
    """Generate all variant that revert decoy sequences to the ``reference_sequence``.

    Expand the selected sequences by generating all the combinatorial options to revert to
    the reference (WT) sequence.

    Alters the names of the designs in **description** column.

    .. warning::
        This is a **computationaly expensive** function. Take this in consideration when trying
        to run it. When trying to revert too many mutations, it can run out of memory.

    :param str seqID: |seqID_param|.
    :param key_residues: |keyres_param|
    :type key_residues: |keyres_types|

    :return: :class:`.DesignFrame`

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz",
           ...:                         {'scores': ['score', 'description'], 'sequence': 'B'})
           ...: df.add_reference_sequence('B', df.get_sequence('B').values[0])
           ...: key_res = [3, 5, 8, 12, 15, 19, 25, 27]
           ...: df.iloc[1].generate_wt_reversions('B', key_res).identify_mutants('B')
    """
    from rstoolbox.components import get_selection

    def format_mutations( row, seqID, key_residues ):
        shift = row.get_reference_shift(seqID)
        seq = row.get_sequence(seqID)
        kr = get_selection(key_residues, seqID, shift, len(seq))
        mutations = row.get_mutations(seqID).split(",")
        muts = []
        if mutations != ['']:
            for m in mutations:
                m = m.strip()
                pos = int(re.search("(\d+)", m).group(1))
                if pos in kr:
                    muts.append((pos, "".join([m[0], m[-1]])))
        return muts

    adf = []
    if "mutants_{}".format(seqID) not in self:
        idf = self.identify_mutants(seqID)
    else:
        idf = self.copy()

    if isinstance(idf, pd.DataFrame):
        designs = []
        for _, row in idf.iterrows():
            mutations = format_mutations(row, seqID, key_residues)
            designs.append(row.generate_mutant_variants(seqID, mutations))
        df = pd.concat(designs)
    elif isinstance(idf, pd.Series):
        mutations = format_mutations(idf, seqID, key_residues)
        df = idf.generate_mutant_variants(seqID, mutations)
    else:
        raise NotImplementedError
    adf.append(df)
    adf = pd.concat(adf)

    avail_seqs = adf.get_available_sequences()
    seqs = [_check_column(adf, "sequence", seq) for seq in avail_seqs]
    adf.drop_duplicates(seqs, inplace=True)
    return adf.reset_index(drop=True)


def score_by_pssm( self, seqID, matrix ):
    """Score sequences according to a provided PSSM matrix.

    Generates new column by applying the PSSM score to each position
    of the requested sequences:

    ======================  =====================================
    New Column                                       Data Content
    ======================  =====================================
    **pssm_score_<seqID>**  Score obtained by applying ``matrix``
    ======================  =====================================

    :param str seqID: |seqID_param|.
    :param matrix: Positional frequency matrix. **column:** residue type; **index:**
        sequence position.
    :type matrix: :class:`~pandas.DataFrame`

    :return: Union[:class:`.DesignSeries`, :class:`.DesignFrame`]
        - Itself with the new column.

    :raises:
        :NotImplementedError: |indf_error|.
        :ValueError: if matrix rows do not match sequence length.

    .. seealso::
        :meth:`.DesignFrame.generate_mutants_from_matrix`
        :meth:`.DesignSeries.generate_mutants_from_matrix`

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: from rstoolbox.tests.helper import random_frequency_matrix
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz",
           ...:                         {'scores': ['score', 'description'], 'sequence': 'B'})
           ...: matrix = random_frequency_matrix(len(df.get_sequence('B')[0]), 0)
           ...: df.score_by_pssm('B', matrix)
    """
    def evaluate_sequence(seq, pssm):
        score = 0
        if len(seq) != pssm.shape[0]:
            raise ValueError("Lenght of sequence and matrix do not match")
        for i, aa in enumerate(seq):
            if aa in list(pssm.columns):
                score += pssm.iloc[i][aa]
        return score

    outcol = "pssm_score_{}".format(seqID)
    if isinstance(self, pd.Series):
        self[outcol] = evaluate_sequence(self.get_sequence(seqID), matrix)
    elif isinstance(self, pd.DataFrame):
        self[outcol] = self.apply(lambda row: evaluate_sequence(row.get_sequence(seqID), matrix),
                                  axis=1)
    else:
        raise NotImplementedError

    return self


def make_resfile( self, seqID, header, filename, write=True ):
    """Generate a Rosetta `resfile
    <https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/resfiles>`_
    to match the design's sequence assuming the ``reference_sequence`` as the starting
    point.

    The function picks the mutant positions between the design sequence
    and the ``reference_sequence`` and generates a resfile that will
    transform one into the other.

    If more than one design is provided, the resfile name (``filename``) is going
    to be modified with a counter. A new column will be added to the data container
    in order to keep track of the filename assignation:

    ========================  =========================
    New Column                             Data Content
    ========================  =========================
    **resfile_<seqID>**             Name of the resfile
    ========================  =========================

    :param str seqID: |seqID_param|.
    :param str header: Header content for the resfile; defines default behaviour.
    :param str filename: Identifier of the resfile. Will be altered with a numerical
        suffix if the data container holds more thant one sequence.
    :param bool write: **Testing attribute**. When :data:`False`, resfiles are not
        actually created

    :return: Union[:class:`.DesignSeries`, :class:`.DesignFrame`]
        - Itself with the new column.

    :raise:
        :NotImplementedError: |indf_error|.
        :KeyError: |reference_error|
        :IOError: |overwrite_error|

    .. note::
        Depends on :ref:`system.overwrite <options>` and
        :ref:`system.output <options>`.

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz",
           ...:                         {'scores': ['score', 'description'], 'sequence': 'B'})
           ...: df.add_reference_sequence('B', df.get_sequence('B').values[0])
           ...: dfwt = df.iloc[0].generate_mutant_variants('B', [(1, "TGP"), (6, "ERG"),
           ...:                                                  (14, "MAT")])
           ...: # Call in test-mode
           ...: dfwt = dfwt.make_resfile("B", "NATAA", "mutants.resfile", write=False )
           ...: dfwt.head()
    """
    if not self.has_reference_sequence(seqID):
        raise KeyError("A reference sequence for {} is needed.".format(seqID))

    def resfile( row, seqID, header, filename, suffix):
        if not isinstance(row, pd.Series):
            raise NotImplementedError

        if suffix is not None:
            filename = list(os.path.splitext(filename))
            filename[0] += "_{:>04d}".format(suffix)
            filename = "".join(filename)

        if core.get_option('system', 'output') != "./":
            if os.path.basename(filename) == filename:
                filename = os.path.join(core.get_option('system', 'output'), filename)
        if not core.get_option('system', 'overwrite'):
            if os.path.isfile(filename):
                raise IOError('{} already exists and cannot be overwriten'.format(filename))

        data = [header, "START\n"]
        df = row.copy()
        if seqID not in df.get_identified_mutants():
            df = df.identify_mutants(seqID)
        shift = df.get_reference_shift(seqID)
        if len(df.get_mutations(seqID)) > 0:
            for mutation in df.get_mutations(seqID).split(","):
                if isinstance(shift, int):
                    position = str(int(mutation[1:-1]) + shift - 1)
                else:
                    position = str(shift[int(mutation[1:-1]) - 1])
                data.append(str(" ".join([position, seqID, "PIKAA", mutation[-1]])))

        if write:
            with open(filename, 'w') as fd:
                fd.write("\n".join(data))
        return filename

    outcol = "resfile_{}".format(seqID)
    if isinstance(self, pd.Series):
        self[outcol] = resfile(self, seqID, header, filename, None)
    elif isinstance(self, pd.DataFrame):
        self[outcol] = self.apply(lambda row: resfile(row, seqID, header, filename, row.name),
                                  axis=1)
    else:
        raise NotImplementedError
    return self


def apply_resfile( self, seqID, filename, rscript=None, keep_input_scores=False ):
    """Apply a generated Rosetta `resfile
    <https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/resfiles>`_
    to the decoy.

    This function needs to be created after the appropiate mutant variants have been created
    and their corresponding **resfiles** have been written.

    .. note::
        Depends on :ref:`rosetta.path <options>` and :ref:`rosetta.compilation <options>`,
        if the ``filename`` does not exist.

    .. attention::
        This function **REQUIRES** a local installation of **Rosetta**.

    To execute this function it is important that the ``source_file`` assigned to the
    :class:`.DesignFrame` is an original silent file and **not a minisilent**, as the
    original structure of the decoy needs to be used in order to generate the variants.
    If that is not the case, use :class:`.DesignFrame.replace_source_files`.

    :param str seqID: |seqID_param|
    :param str filename: Name of the final silent file that will contain all the variant's data.
        If the file exists, it is assumed that the data was already created and data will be
        directly loaded from that file.
    :param str rscript: By default, the script executed will be the one generated by
        :func:`.mutations`. One can provide its own script (either as the file name of the
        script or as a string of the content itself) **as long as it fulfills two conditions**:
        (1) It must contain the **AddJobPairData Mover** and (2) it should accept the script
        variable ``resfile``. An example on how to use these two conditions can be extrapolated
        from :func:`.mutations`.
    :param bool keep_input_scores: When :data:`True` (default :data:`False`), it will keep the
        score terms present in the source decoy (as they appear in the original silent file)
        for the variants.

    :return: :class:`.DesignFrame` with the scores for the mutants.

    :raise:
        :SystemError: If all variants faile to be generated or if they cannot be merged.
        :IOError: If Rosetta path cannot be found.
        :AttributeError: If the resfiles for the variants were not previously created.

    .. seealso:
        :meth:`.DesignFrame.generate_mutant_variants`
        :meth:`.DesignFrame.generate_mutants_from_matrix`
        :meth:`.DesignFrame.generate_wt_reversions`
        :meth:`.DesignFrame.make_resfile`
        :meth:`.DesignSeries.generate_mutant_variants`
        :meth:`.DesignSeries.generate_mutants_from_matrix`
        :meth:`.DesignSeries.generate_wt_reversions`
        :meth:`.DesignSeries.make_resfile`

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz",
           ...:                         {'scores': ['score', 'description'], 'sequence': 'B'})
           ...: df.add_reference_sequence('B', df.get_sequence('B').values[0])
           ...: dfwt = df.iloc[0].generate_mutant_variants('B', [(1, "TGP"), (6, "ERG"),
           ...:                                                  (14, "MAT")])
           ...: # Call in test-mode
           ...: dfwt = dfwt.make_resfile("B", "NATAA", "mutants.resfile", write=False )
           ...: dfwt2 = dfwt.iloc[:3].apply_resfile("B",
           ...:                                     "../rstoolbox/tests/data/variants.silent.gz")
           ...: dfwt2
    """
    from rstoolbox.components import DesignSeries, DesignFrame
    from rstoolbox.io import parse_rosetta_file
    from rstoolbox.utils import mutations

    if isinstance(self, DesignSeries):
        self = DesignFrame(self).T

    resfile = 'resfile_{}'.format(seqID)
    if not os.path.isfile(filename):
        wdir = tempfile.mkdtemp()
        exe = make_rosetta_app_path('rosetta_scripts')
        if resfile not in self.columns:
            raise AttributeError("Resfiles are needed to execute this function.")
        if rscript is None:
            rscript = mutations(seqID)
        if not os.path.isfile(rscript):
            fd = open(os.path.join(wdir, 'script.xml'), 'w')
            fd.write(rscript)
            fd.close()
            rscript = os.path.join(wdir, 'script.xml')

        command = ['{0}', '-parser:protocol {1}', '-in:file:silent {2}', '-in:file:tags {3}',
                   '-out:file:silent {4}', '-parser:script_vars resfile={5}']
        if not keep_input_scores:
            command.append('-keep_input_scores false')
        command = ' '.join(command)
        outfiles = []
        errors = 0
        sys.stdout.write("Running Rosetta\n")
        for _, row in self.iterrows():
            if re.search('_v\d{4}$', row['description']):
                origin = "_".join(row['description'].split('_')[:-1])
            else:
                origin = row['description']
            outfiles.append(os.path.join(wdir, row['description'] + '.silent'))
            cmd = command.format(exe, rscript, " ".join(self.get_source_files()),
                                 origin, outfiles[-1], row[resfile])
            sys.stdout.write(cmd + "\n")
            error = execute_process( command )
            if bool(error):
                errors += 1
                sys.stdout.write("Execution for variant {} has failed\n".format(row['description']))

        if errors < self.shape[0]:
            exe = make_rosetta_app_path('combine_silent')
            command = ['{0}', '-in:file:silent {1}', '-out:file:silent {2}']
            command = ' '.join(command)
            cmd = command.format(exe, " ".join(outfiles), filename)
            sys.stdout.write("Merging all silent files\n")
            sys.stdout.write(cmd + "\n")
            error = execute_process( command )
            if bool(error):
                raise SystemError("A file with the new variants could not be created.")
        else:
            raise SystemError("All variants failed to be generated.")

    df = parse_rosetta_file(filename)
    df = df.drop(columns=['description'])
    return self.merge(df, on=resfile, how='left')


def view_mutants_alignment( self, seqID, mutants_bg_color="IndianRed", mutants_text_color="white",
                            identities_bg_color="YellowGreen", identities_text_color="black"):
    """Generates a pretty representation alignment of the mutations in **Jupyter Notebooks**.

    :param str seqID: |seqID_param|.
    :param str mutants_bg_color: Color to apply to the background of mutants.
    :param str mutants_text_color: Color to apply to the text of mutants.
    :param str identities_bg_color: Color to apply to the background of identities.
    :param str identities_text_color: Color to apply to the text of identities.

    :raise:
        :KeyError: |reference_error|.
    """
    if not self.has_reference_sequence(seqID):
        raise KeyError("A reference sequence for {} is needed.".format(seqID))

    def highlight_mutants( row, refseq, background="IndianRed", text="white" ):
        attr = 'background-color: {0}; color: {1}'.format(background, text)
        is_mut = row != refseq
        return [attr if v else '' for v in is_mut]

    def highlight_identities( row, refseq, background="YellowGreen", text="black" ):
        attr = 'background-color: {0}; color: {1}'.format(background, text)
        is_mut = row == refseq
        return [attr if v else '' for v in is_mut]

    seq = list(self.get_reference_sequence(seqID))
    pos = self.get_reference_shift(seqID)
    if isinstance(pos, int):
        pos = list(range(pos, len(seq) + pos))
    cols = ["\n".join(list("{:>03d}".format(x))) + "\n" + y for x, y in zip(pos, seq)]
    if seqID not in self.get_identified_mutants():
        df = self.identify_mutants(seqID)
    else:
        df = self.copy()

    df = df[_check_column(df, "sequence", seqID)].apply(lambda x: pd.Series(list(x)))
    df.columns = cols
    df.index = self["description"].values

    return df.style.set_caption('sequence alignment') \
             .set_table_attributes('class="dataframe"') \
             .apply(highlight_mutants, refseq=seq, background=mutants_bg_color,
                    text=mutants_text_color, axis=1) \
             .apply(highlight_identities, refseq=seq, background=identities_bg_color,
                    text=identities_text_color, axis=1)
