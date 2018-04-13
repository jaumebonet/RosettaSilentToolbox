# @Author: Jaume Bonet <bonet>
# @Date:   19-Feb-2018
# @Email:  jaume.bonet@gmail.com
# @Filename: sequence.py
# @Last modified by:   bonet
# @Last modified time: 13-Apr-2018


import os
import re
import gzip
import bisect

import rstoolbox.core as core
import rstoolbox.components as cp
from rstoolbox.io.rosetta import _gather_file_list
from rstoolbox.utils.getters import _check_column


__all__ = ["read_fasta", "write_fasta", "write_clustalw", "write_mutant_alignments"]


def read_fasta( filename, expand=False, multi=False ):
    """
    Reads one or more **FASTA** files and returns the appropiate object
    containing the requested data: the :class:`.DesignFrame`.

    The default generated :class:`.DesignFrame` will contain two columns:

    ===============  ===================================================
    Column Name      Data Content
    ===============  ===================================================
    **description**  Sequence identifier.
    **sequence_A**   Sequence content.
    ===============  ===================================================

    The sequence column assigned as ``sequence_A`` is an arbitrary decision that
    has to do compatibility issues with the rest of functions and methods of
    :class:`.DesignFrame`.

    .. ipython::

        In [1]: from rstoolbox.io import read_fasta
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: df = read_fasta("../rstoolbox/tests/data/*fa", multi=True)
           ...: df

    If the **FASTA** comes or is formated as **PDB FASTA** (as in the example avobe),
    it is possible to better assign the column names to the actual sequence ID. To force
    that behaviour, activate the ``expand`` option.

    .. ipython::

        In [1]: from rstoolbox.io import read_fasta
           ...: df = read_fasta("../rstoolbox/tests/data/*fa", expand=True, multi=True)
           ...: df

    .. note::
        Notice everything from the original ``description`` after the ``|`` symbol
        is lost after that process.

    :param filename: file name or file pattern to search.
    :type filename: :class:`str`
    :param expand: Try to better associate sequence ID if format is **PDB FASTA**.
    :type expand: :class:`bool`
    :param multi: When :data:`True`, indicates that data is readed from
        multiple files.
    :type multi: :class:`bool`

    :return: :py:class:`.DesignFrame`.

    :raises:
        :IOError: if ``filename`` cannot be found.

    .. note::
        Depends on :ref:`system.output <options>`.

    .. seealso::
        :func:`~.write_fasta`
    """
    files = _gather_file_list( filename, multi )
    data = {"description": [], "sequence_A": []}
    for file_count, f in enumerate( files ):
        fd = gzip.open( f ) if f.endswith(".gz") else open( f )
        for line in fd:
            line = line.decode('utf8') if f.endswith(".gz") else line
            line = line.strip()
            if line.startswith(">"):
                line = line.strip(">")
                data["description"].append(line)
                data["sequence_A"].append("")
            elif len(line) > 0:
                data["sequence_A"][-1] += line

    df = cp.DesignFrame( data )
    if expand and bool(re.search("^\S{4}\:\S{1}", df.iloc[0]["description"])):
        df["description"] = df["description"].apply(lambda col: col.split("|")[0])
        df[['description', 'seq']] = df['description'].str.split(':', expand=True)
        df = df.pivot('description', 'seq',
                      'sequence_A').add_prefix("sequence_").rename_axis(None,
                                                                        axis=1).reset_index()
        df = cp.DesignFrame(df)
    df.add_source_files( files )
    return df


def write_fasta( df, seqID, separator=None, filename=None, split=False ):
    """
    Writes fasta files of the selected decoys.

    It assumes that the provided data is contained in a :class:`.DesignFrame`
    or a :class:`~pandas.DataFrame`.

    Mandatory columns are:

    ====================  ===================================================
    Column Name           Data Content
    ====================  ===================================================
    **description**       Sequence identifier.
    **sequence_<seqID>**  Sequence content.
    ====================  ===================================================

    .. ipython::

        In [1]: from rstoolbox.io import read_fasta, write_fasta
           ...: df = read_fasta("../rstoolbox/tests/data/*fa", multi=True)
           ...: print write_fasta(df, "A")

    When working with multiple ``seqID``, one can select which ones to be printed;
    empty sequences will be skipped.

    .. ipython::

        In [1]: from rstoolbox.io import read_fasta, write_fasta
           ...: df = read_fasta("../rstoolbox/tests/data/*fa", expand=True, multi=True)
           ...: print write_fasta(df, "AC")

    :param df: Data content.
    :type df: Union[:class:`.DesignFrame`, :class:`~pandas.DataFrame`]
    :param seqID: Identifier(s) of the sequences expected to be printed.
    :type seqID: :class:`str`
    :param separator: Add ``seqID`` to sequence identifier through a particular
        string separator. If multiple ``seqID`` are provided, it defaults to ``:``.
    :type separator: :class:`str`
    :param filename: Output file name.
    :type filename: :class:`str`
    :param split: Split each fasta in a different file. ``filename`` first part of the filename
        is used as `prefix`, with a following enumeration.
    :type split: :class:`bool`

    :return: :class:`str` - **FASTA** formated string.

    :raises:
        :IOError: If ``filename`` exists and global option :ref:`system.overwrite <options>`
            is not :data:`True`.
        :AttributeError: If any of the requested seqID cannot be found.

    .. note::
        Depends on :ref:`system.overwrite <options>` and
        :ref:`system.output <options>`.

    .. seealso::
        :func:`~.write_fasta`
    """
    def nomenclator(row, seqID, separator):
        sequence = row.get_sequence(seqID)
        if sequence is None or len(sequence) == 0:
            return ""
        name = ">" + row.get_id()
        if separator is not None:
            name = name + separator + seqID
        return name + "\n" + row.get_sequence(seqID)

    if filename is not None:
        if os.path.isfile(filename) and not core.get_option("system", "overwrite"):
            raise IOError("File {} already exists".format(filename))
    if not isinstance(df, cp.DesignFrame):
        df = cp.DesignFrame(df)
    if len(seqID) > 0 and separator is None:
        separator = ":"

    data = []
    for chain in seqID:
        eachfa = df.apply(lambda row: nomenclator(row, chain, separator), axis=1)
        data.extend(eachfa.values)

    if filename is not None:
        if not split:
            fd = open(filename, "w") if not filename.endswith(".gz") else gzip.open(filename, "wb")
            fd.write("\n".join(data).strip() + "\n")
            fd.close()
        else:
            suffix = "_f{0:04d}"
            cplxname = os.path.splitext(filename)
            for i, sequence in enumerate(data):
                fname = cplxname[0] + suffix.format(i + 1) + cplxname[1]
                fd = open(fname, "w") if not fname.endswith(".gz") else gzip.open(fname, "wb")
                fd.write(sequence + "\n")
                fd.close()

    return "\n".join(data).strip() + "\n"


def write_clustalw( df, seqID, filename=None ):
    """
    Write sequences of selected designs as a CLUSTALW alignment.

    If a ``reference_sequence`` exists, it is set up as the first sequence
    of the alignment. The name assigned to it will be the multipl longest common
    subsequence of all the decoys ``description``. If none is found, or if it
    actually matches one of the already existing identifiers, then it will default
    to *reference*.

    :param df: Data content.
    :type df: Union[:class:`.DesignFrame`, :class:`~pandas.DataFrame`]
    :param seqID: Identifier(s) of the sequences expected to be printed.
    :type seqID: :class:`str`
    :param filename: Output file name.
    :type filename: :class:`str`

    :return: :class:`str` - **CLUSTALW** formated string.

    :raises:
        :IOError: If ``filename`` exists and global option :ref:`system.overwrite <options>`
            is not :data:`True`.
        :AttributeError: If any of the requested seqID cannot be found.

    .. note::
        Depends on :ref:`system.overwrite <options>` and
        :ref:`system.output <options>`.
    """
    def chunkstring(string, length):
        return [string[0 + i:length + i] for i in range(0, len(string), length)]

    if filename is not None:
        if os.path.isfile(filename) and not core.get_option("system", "overwrite"):
            raise IOError("File {} already exists".format(filename))

    data = ["CLUSTAL W(1.83) multiple sequence alignment\n"]
    names = list(df.get_id().values)
    seqs = list(df.get_sequence(seqID).values)
    if df.has_reference_sequence(seqID):
        refname = mlcs(names)
        if len(refname) > 0 and refname not in names:
            if re.match('[\w\d\_\.]+\_[0]*$', str(refname)):
                refname = refname.rstrip("0").rstrip("_")
        else:
            refname = "reference"
        names.insert(0, refname)
        seqs.insert(0, df.get_reference_sequence(seqID))
    seqs = [chunkstring(_, 50) for _ in seqs]
    chunks = len(seqs[0])
    nm_len = max([len(_) for _ in names])
    for j in range(chunks):
        for i, nm in enumerate(names):
            line = ("{:<" + str(nm_len) + "}").format(nm) + " " + seqs[i][j]
            data.append(line)
        data.append("\n")

    if filename is not None:
        fd = open(filename, "w") if not filename.endswith(".gz") else gzip.open(filename, "wb")
        fd.write("\n".join(data))
        fd.close()

    return "\n".join(data)


def write_mutant_alignments( df, seqID, filename=None ):
    """
    Writes a text file containing only the positions changed with respect to
    the ``reference_sequence``.

    The format will be CLUSTALW but residues that repeat the ``reference_sequence``
    are substituted by ".".

    :param df: Data content.
    :type df: Union[:class:`.DesignFrame`, :class:`~pandas.DataFrame`]
    :param seqID: Identifier(s) of the sequences expected to be printed.
    :type seqID: :class:`str`
    :param filename: Output file name.
    :type filename: :class:`str`

    :return: :class:`str` - **CLUSTALW** formated string.

    :raises:
        :IOError: If ``filename`` exists and global option :ref:`system.overwrite <options>`
            is not :data:`True`.
        :AttributeError: If any of the requested seqID cannot be found.

    .. note::
        Depends on :ref:`system.overwrite <options>` and
        :ref:`system.output <options>`.
    """
    def mask_row(row, seqID):
        pos = row.get_mutation_positions(seqID)
        if len(pos) == 0:
            return row.get_sequence(seqID)
        pos = [int(_) for _ in row.get_mutation_positions(seqID).split(",")]
        seq = list(row.get_sequence(seqID))
        sft = row.get_reference_shift(seqID)
        return "".join([aa if (i + sft) in pos else "." for i, aa in enumerate(seq)])

    df = df.identify_mutants(seqID)
    df[_check_column(df, "sequence", seqID)] = df.apply(lambda row: mask_row(row, seqID), axis=1)

    return write_clustalw(df, seqID, filename)


def mlcs(strings):
    """Return a long common subsequence of the strings.

    Uses a greedy algorithm, so the result is not necessarily the
    longest common subsequence.

    """
    # https://codereview.stackexchange.com/a/90381/163119
    if not strings:
        raise ValueError("mlcs() argument is an empty sequence")
    strings = list(set(strings))  # deduplicate
    alphabet = set.intersection(*(set(s) for s in strings))

    # indexes[letter][i] is list of indexes of letter in strings[i].
    indexes = {letter: [[] for _ in strings] for letter in alphabet}
    for i, s in enumerate(strings):
        for j, letter in enumerate(s):
            if letter in alphabet:
                indexes[letter][i].append(j)

    # pos[i] is current position of search in strings[i].
    pos = [len(s) for s in strings]

    # Generate candidate positions for next step in search.
    def candidates():
        for letter, letter_indexes in indexes.items():
            distance, candidate = 0, []
            for ind, p in zip(letter_indexes, pos):
                i = bisect.bisect_right(ind, p - 1) - 1
                q = ind[i]
                if i < 0 or q > p - 1:
                    break
                candidate.append(q)
                distance += (p - q)**2
            else:
                yield distance, letter, candidate

    result = []
    while True:
        try:
            # Choose the closest candidate position, if any.
            _, letter, pos = min(candidates())
        except ValueError:
            return ''.join(reversed(result))
        result.append(letter)
