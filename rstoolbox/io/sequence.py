# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: read_fasta
.. func:: write_fasta
.. func:: write_clustalw
.. func:: write_mutant_alignments
.. func:: read_hmmsearch
"""
# Standard Libraries
import os
import re
import gzip
import bisect

# External Libraries
import pandas as pd

# This Library
import rstoolbox.core as core
import rstoolbox.components as cp
from rstoolbox.io.rosetta import _gather_file_list
from rstoolbox.utils.getters import _check_column


__all__ = ['read_fasta', 'write_fasta', 'write_clustalw', 'write_mutant_alignments',
           'read_hmmsearch']


def read_fasta( filename, expand=False, multi=False, defchain='A' ):
    """Reads one or more **FASTA** files and returns the appropiate object
    containing the requested data: the :class:`.DesignFrame`.

    The default generated :class:`.DesignFrame` will contain two columns:

    ====================  ===================================================
    Column Name            Data Content
    ====================  ===================================================
    **description**        Sequence identifier.
    **sequence_<chain>**   Sequence content.
    ====================  ===================================================

    The sequence column assigned as ``sequence_A`` is an arbitrary decision that
    has to do compatibility issues with the rest of functions and methods of
    :class:`.DesignFrame`.

    .. ipython::

        In [1]: from rstoolbox.io import read_fasta
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = read_fasta("../rstoolbox/tests/data/*fa$", multi=True)
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

    :param str filename: file name or file pattern to search.
    :param bool expand: Try to better associate sequence ID if format is **PDB FASTA**.
    :param bool multi: When :data:`True`, indicates that data is readed from
        multiple files.
    :param str defchain: Default chain to use. If not provided that is 'A'.

    :return: :class:`.DesignFrame`.

    :raises:
        :IOError: if ``filename`` cannot be found.

    .. seealso::
        :func:`~.write_fasta`
    """
    seqcol = "sequence_{}".format(defchain)
    files = _gather_file_list( filename, multi )
    data = {"description": [], seqcol: []}
    for _, f in enumerate( files ):
        fd = gzip.open(f, 'rt') if f.endswith('.gz') else open( f )
        for line in fd:
            line = line.strip('\n')
            if line.startswith(">"):
                line = line.strip(">")
                data["description"].append(line)
                data[seqcol].append("")
            elif len(line) > 0:
                data[seqcol][-1] += line

    df = cp.DesignFrame( data )
    if expand and bool(re.search("^\S{4}\:\S{1}", df.iloc[0]["description"])):
        df["description"] = df["description"].apply(lambda col: col.split("|")[0])
        df[['description', 'seq']] = df['description'].str.split(':', expand=True)
        df = df.pivot('description', 'seq',
                      seqcol).add_prefix("sequence_").rename_axis(None,
                                                                  axis=1).reset_index()
        df = cp.DesignFrame(df)
    df.add_source_files( files )
    return df


def write_fasta( df, seqID, separator=None, filename=None, split=False ):
    """Writes fasta files of the selected decoys.

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
           ...: print(write_fasta(df, "A"))

    When working with multiple ``seqID``, one can select which ones to be printed;
    empty sequences will be skipped.

    .. ipython::

        In [1]: from rstoolbox.io import read_fasta, write_fasta
           ...: df = read_fasta("../rstoolbox/tests/data/*fa", expand=True, multi=True)
           ...: print(write_fasta(df, "AC"))

    :param df: Data content.
    :type df: Union[:class:`.DesignFrame`, :class:`~pandas.DataFrame`]
    :param str seqID: |seqID_param|.
    :param str separator: Add ``seqID`` to sequence identifier through a particular
        string separator. If multiple ``seqID`` are provided, it defaults to ``:``.
    :param str filename: Output file name.
    :param bool split: Split each fasta in a different file. ``filename`` first part of the filename
        is used as `prefix`, with a following enumeration.

    :return: :class:`str` - **FASTA** formated string.

    :raises:
        :IOError: If ``filename`` exists and global option :ref:`system.overwrite <options>`
            is not :data:`True`.
        :AttributeError: |seqID_error|.

    .. note::
        Depends on :ref:`system.overwrite <options>` and :ref:`system.output <options>`.

    .. seealso::
        :func:`~.read_fasta`
    """
    def nomenclator(row, seqID, separator):
        sequence = row.get_sequence(seqID)
        if sequence is None or isinstance(sequence, float) or len(sequence) == 0:
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
    """Write sequences of selected designs as a CLUSTALW alignment.

    If a ``reference_sequence`` exists, it is set up as the first sequence
    of the alignment. The name assigned to it will be the multipl longest common
    subsequence of all the decoys ``description``. If none is found, or if it
    actually matches one of the already existing identifiers, then it will default
    to *reference*.

    :param df: Data content.
    :type df: Union[:class:`.DesignFrame`, :class:`~pandas.DataFrame`]
    :param str seqID: |seqID_param|.
    :type seqID: :class:`str`
    :param filename: Output file name.
    :type filename: :class:`str`

    :return: :class:`str` - **CLUSTALW** formated string.

    :raises:
        :IOError: If ``filename`` exists and global option :ref:`system.overwrite <options>`
            is not :data:`True`.
        :AttributeError: |seqID_error|.

    .. note::
        Depends on :ref:`system.overwrite <options>` and :ref:`system.output <options>`.
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
    """Writes a text file containing only the positions changed with respect to
    the ``reference_sequence``.

    The format will be **CLUSTALW** but residues that repeat the ``reference_sequence``
    are substituted by ".".

    :param df: Data content.
    :type df: Union[:class:`.DesignFrame`, :class:`~pandas.DataFrame`]
    :param str seqID: |seqID_param|.
    :param str filename: Output file name.

    :return: :class:`str` - **CLUSTALW** formated string.

    :raises:
        :IOError: If ``filename`` exists and global option :ref:`system.overwrite <options>`
            is not :data:`True`.
        :AttributeError: |seqID_error|.

    .. note::
        Depends on :ref:`system.overwrite <options>` and :ref:`system.output <options>`.
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


def read_hmmsearch( filename ):
    """Read output from ``hmmsearch`` or ``hmmscan``.

    Processess the output of Hidden Markov Models search over a set
    of sequences with `hmmsearch <http://hmmer.org/>`_.

    Will return a :class:`~pandas.DataFrame` with the following columns:

    ====================  ===================================================
    Column Name           Data Content
    ====================  ===================================================
    **query**             Query identifier.
    **description**       In ``hmmseach`` only. Sequence identifier.
    **domain**            In ``hmmscan`` only.  Domain identifier.
    **definition**        In ``hmmscan`` only. Definition of the domain name.
    **full-e-value**      E-value for full sequence match.
    **full-score**        Score for full sequence match.
    **full-bias**         Bias for full sequence.
    **dom-e-value**       E-values for best scored domain.
    **dom-score**         Score for best scored domain.
    **dom-bias**          Bias for best scored domain.
    **dom-exp**           Expected number of domains.
    **dom-N**             Actual number of domains.
    **sequence**          Part of the sequence aligned on its own.
    ====================  ===================================================

    It will include more columns related with the alignment data itself.

    :param str filename: Name of the ``hmmsearch`` output file.

    :return: :class:`~pandas.DataFrame`

    :raises:
        :IOError: if ``filename`` does not exist.

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import read_hmmsearch
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = read_hmmsearch("../rstoolbox/tests/data/search.hmm.gz")
           ...: df.head()
    """
    if not os.path.isfile(filename):
        raise IOError('{} cannot be found'.format(filename))

    ali = re.compile('\d+\s[\!\?][\s\S]*')
    fd = open(filename)if not filename.endswith("gz") else gzip.open(filename, 'rt')
    searches = fd.read().split("\n//")
    fd.close()

    dfs = []
    mode = ''
    for _, search in enumerate(searches):
        data = {'full-e-value': [], 'full-score': [], 'full-bias': [],
                'dom-e-value': [], 'dom-score': [], 'dom-bias': [],
                'dom-exp': [], 'dom-N': []}
        dat2 = {'score': [], 'bias': [], 'c-Evalue': [],
                'i-Evalue': [], 'hmmfrom': [], 'hmmto': [], 'alifrom': [],
                'alito': [], 'envfrom': [], 'envto': [], 'acc': [], 'sequence': []}
        nam2 = ''
        seq = ''
        read = 0
        nohits = None
        for line in search.split('\n'):
            if line.startswith('#'):
                if 'hmmsearch' in line:
                    mode = 'hmmsearch'
                if 'hmmscan' in line:
                    mode = 'hmmscan'
                continue
            if len(line.strip()) == 0:
                continue
            if line.startswith('Query:'):
                query = [s for s in line.split(" ") if s != ''][1]
                continue
            if line.strip().startswith('------') and read != 2:
                read = 1
                continue
            if line.startswith('Domain annotation for each sequence'):
                # hmmsearch
                read = 2
                continue
            if line.startswith('Domain annotation for each model (and alignments):'):
                # hmmscan
                read = 2
                continue
            if "[No hits detected that satisfy reporting thresholds]" in line:
                # hmmscan - No hits, empty DataFrame
                nohits = True
                data.setdefault('query', []).append(query)
                cols = data
                cols.update(dat2)
                cols.setdefault('domain', [])
                cols.setdefault('definition', [])
                for col in cols:
                    if col == "query":
                        continue
                    else:
                        data[col].append('-')
                continue
            if read == 1:
                line = line.strip().split()
                data['full-e-value'].append(float(line[0]))
                data['full-score'].append(float(line[1]))
                data['full-bias'].append(float(line[2]))
                data['dom-e-value'].append(float(line[3]))
                data['dom-score'].append(float(line[4]))
                data['dom-bias'].append(float(line[5]))
                data['dom-exp'].append(float(line[6]))
                data['dom-N'].append(float(line[7]))
                data.setdefault('query', []).append(query)

                if mode == 'hmmsearch':
                    data.setdefault('description', []).append(line[8])
                if mode == 'hmmscan':
                    data.setdefault('domain', []).append(line[8])
                    data.setdefault('definition', []).append(" ".join(line[9:]))
            if read == 2:
                if line.startswith('>>'):
                    nam2 = line.split()[1].strip()
                    continue
                if re.match(ali, line.strip()):
                    if mode == 'hmmsearch':
                        dat2.setdefault('description', []).append(nam2)
                    if mode == 'hmmscan':
                        dat2.setdefault('domain', []).append(nam2)
                    lnp = line.strip().split()
                    dat2['score'].append(float(lnp[2]))
                    dat2['bias'].append(float(lnp[3]))
                    dat2['c-Evalue'].append(float(lnp[4]))
                    dat2['i-Evalue'].append(float(lnp[5]))
                    dat2['hmmfrom'].append(int(lnp[6]))
                    dat2['hmmto'].append(int(lnp[7]))
                    dat2['alifrom'].append(int(lnp[9]))
                    dat2['alito'].append(int(lnp[10]))
                    dat2['envfrom'].append(int(lnp[12]))
                    dat2['envto'].append(int(lnp[13]))
                    dat2['acc'].append(float(lnp[15]))
                    continue
                if line.startswith('  == domain '):
                    if seq != '':
                        dat2['sequence'].append(seq)
                        seq = ''
                    continue
                if mode == 'hmmsearch':
                    if nam2 in line.strip() and nam2 != '':
                        seq += line.strip().split()[2]
                        continue
                if mode == 'hmmscan':
                    if query in line.strip() and query != '':
                        seq += line.strip().split()[2]
                        continue

        if seq != '' and nohits is None:
            dat2['sequence'].append(seq)
        if nohits == True:
            df = pd.DataFrame(cols)
        else:
            if len(dat2['sequence']) == 0:
                dat2['sequence'] = ['', ] * len(dat2['score'])
            onid = 'description'
            if mode == 'hmmscan':
                onid = 'domain'
            if pd.DataFrame(dat2).empty or pd.DataFrame(data).empty:
                continue
            df = pd.merge(pd.DataFrame(data), pd.DataFrame(dat2),
                          on=onid, how='outer').fillna(0)
        dfs.append(df)
    return pd.concat(dfs, sort=True)


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
