# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. class:: FragmentFrame
"""
# Standard Libraries
import os
import math
import sys
from collections import Counter
import itertools

# External Libraries
import pandas as pd
import numpy as np
import networkx as nx

# This Library
from rstoolbox.utils import make_rosetta_app_path, execute_process

__all__ = ["FragmentFrame"]


def _metadata_defaults(name):
    if name == "_source_file":
        return None
    return None


class FragmentFrame( pd.DataFrame ):
    """Data container for Fragment data.

    Extends :py:class:`pandas.DataFrame` adding some functionalities in order
    to facilitate analysis, plotting and comparisson.

    Filled through the functions provided through this library, each
    row represents a position of a fragment, while columns describe the
    properties.

    .. seealso::
        :func:`.parse_rosetta_fragments`
        :func:`.plot_fragments`
        :func:`.plot_fragment_profiles`

    . rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_fragments
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_fragments("../rstoolbox/tests/data/wauto.200.3mers.gz")
           ...: df.head()
    """
    _metadata = ['_source_file']
    _internal_names = pd.DataFrame._internal_names + ['_crunched']
    _internal_names_set = set(_internal_names)

    def __init__(self, *args, **kw):
        source_file = kw["file"] if "file" in kw else None
        if "file" in kw:
            del(kw["file"])
        super(FragmentFrame, self).__init__(*args, **kw)
        self._crunched = {}
        self._source_file = source_file

    def add_source_file( self, file ):
        """Adds a source file to the :class:`.FragmentFrame`.

        If loaded through :func:`.parse_rosetta_fragments`, the source file is
        automatically added.

        This can be used to automatically generate the fragment RMSD quality.

        :param str file: Name of the file.

        .. seealso:
            :meth:`FragmentFrame.add_quality_measure`
        """
        self._source_file = os.path.abspath(file)

    def get_source_file( self ):
        """Returns the file name linked to this :class:`.FragmentFrame`.

        :return: :class:`str` - Source file path.
        """
        return self._source_file

    def has_source_file( self ):
        """Checks if there is a source file added.

        :return: bool
        """
        return self._source_file is not None

    def is_comparable( self, df ):
        """Evaluate if the current :class:`.FragmentFrame` is comparable to
        the provided one.

        This is checked on terms of (a) covered sequence length -range- and
        (2) fragment size.

        :return: :class:`bool` - True if the two :class:`.FragmentFrame` are comparable.
        """
        if max(self["position"]) != max(df["position"]):
            return False
        if min(self["position"]) != min(df["position"]):
            return False
        if self["size"].values[0] != df["size"].values[0]:
            return False
        return True

    def add_quality_measure( self, filename, pdbfile=None ):
        """Add RMSD quality measure to the fragment data.

        The RMSD quality measurement is performed by the ``r_fraq_qual`` application
        from **Rosetta**. It can be called as:

        .. code-block:: bash

           r_fraq_qual.linuxgccrelease -in:file:native <pdb> -f <fragfile> -out:qual <output>

        :param str filename: Name containing the quality measure. If ``filename`` is None,
            it assumes that RMSD quality has not been calculated yet, so it'll run the
            ``r_fraq_qual`` application as long as the :class:`.FragmentFrame`
            has a ``source_file``. Standart output will be the name of the source
            file with the extension ".qual". If a file with this naming schema exists,
            it will be automatically picked. To be able to run *Rosetta* it will need
            a ``pdbfile``.
        :param str pdbfile: In case the quality has to be calculated. Provide the
            PDB over which to calculate it. Default is :data:`None`.

        :raises:
            :IOError: if ``filename`` does not exist.
            :IOError: if ``pdbfile`` is provided and does not exist.
            :IOError: if the *Rosetta* executable is not found.
            :AttributeError: if ``filename`` is :data:`None` and there is no attached
                ``source_file`` to the object.
            :ValueError: if no rmsd data is assigned. Might indicate that wrong data is
                trying to be assigned.

        .. note::
            Depends on :ref:`rosetta.path <options>` and :ref:`rosetta.compilation <options>`,
            if the quality file is not provided.

        .. attention::
            Some configurations of this function require a local installation of **Rosetta**.
        """
        if filename is None and not self.has_source_file():
            raise AttributeError("No quality file is provided and no source file can be found.")

        # Make the quality fragmet eval if needed.
        if filename is None:
            sofi = self._source_file
            if sofi.endswith(".gz"):
                sofi = ".".join(sofi.split(".")[:-1])
            filename = sofi + ".qual"
            if not os.path.isfile(filename) and not os.path.isfile(filename + ".gz"):
                # Check rosetta executable
                exe = make_rosetta_app_path('r_frag_quality')
                if not os.path.isfile(pdbfile):
                    raise IOError("{0} not found".format(pdbfile))
                command = "{0} -in:file:native {1} -f {2} -out:qual {3}".format(
                          exe, pdbfile, self._source_file, filename )
                error = execute_process( command )
                if not bool(error):
                    sys.stdout.write("Execution has finished\n")
                else:
                    sys.stdout.write("Execution has failed\n")
            elif os.path.isfile(filename + ".gz"):
                filename = filename + ".gz"

        # Load the data
        df = pd.read_csv(filename, header=None, sep="\s+",
                         names=["size", "frame", "neighbor", "rmsd", "_null1", "_null2"],
                         usecols=["size", "frame", "neighbor", "rmsd"])

        df = self.merge(df, how='left', on=["size", "frame", "neighbor"])
        if (df['rmsd'].isnull()).all():
            raise ValueError('No rmsd was assigned to the fragment data. '
                             'Check that the correct quality data is being assigned.')
        return df

    def select_quantile( self, quantile=0.25 ):
        """Returns fragments under the rmsd threshold of the specified uantile.

        :param float quantile: Quantile maximum limit.

        :return: :class:`.FragmentFrame` - The filtered data.

        :raises:
            :KeyError: if the ``rmsd`` column cannot be found.

        .. seealso::
            :meth:`~.FragmentFrame.add_quality_measure`
        """
        def _select_quantile(group, quantile):
            qtl = group["rmsd"].quantile(.25)
            return group[group["rmsd"] <= qtl]

        df = self.groupby("frame").apply(lambda g: _select_quantile(g, quantile))
        df.index = df.index.get_level_values(1)
        df._source_file = self._source_file
        return df

    def make_sequence_matrix( self, frequency=False, round_data=False ):
        """Generate a PSSM-like matrix from the fragments.

        The matrix will contain, for each position the relative enrichment of each
        residue type. By default, this will be:

        :math:`logodd(f(ni)/f(bi))`

        unless ``frequency`` is requested.

        :param bool frequency: Return the matrix with frequency values..
        :param bool round_data: Round-floor the values.

        :return: :class:`~pandas.DataFrame`
        """
        alphabet = "ARNDCQEGHILKMFPSTWYV"
        baseline = dict.fromkeys(alphabet, 1.0)
        total = sum(baseline.values())
        for k in baseline:
            baseline[k] = float(baseline[k]) / total

        matrix = {}
        for i in range(1, max(self["position"].values) + 1):
            qseq = Counter(self[self["position"] == i]["aa"].values)
            qttl = sum(qseq.values(), 0.0)
            for k in qseq:
                qseq[k] /= qttl
            for aa in baseline:
                q = qseq[aa]
                if not frequency:
                    if q > 0:
                        logodds = math.log(q / baseline[aa], 2)
                    else:
                        logodds = -9
                    matrix.setdefault(aa, []).append(logodds)
                else:
                    matrix.setdefault(aa, []).append(q)
        df = pd.DataFrame(matrix)
        if round_data:
            df = df.applymap(np.around).astype(np.int64).reindex(columns=list(alphabet))
        else:
            df = df.reindex(columns=list(alphabet))
        df.index = range(1, df.shape[0] + 1)
        return df

    def make_per_position_frequency_network( self ):
        """Generate a graph representation of the per-residue frequency.

        Generate a directed graph of class :class:`~networkx.DiGraph` in which each node is
        a residue type and each edge the frequency expected for that residue type in that
        position according to the frequency matrix. As a matter of fact, the edge is the
        *inverted frequency*, as the idea is to use the graph to calculate the shortest paths
        (i.e. the more probable sequences).

        :return: :class:`~networkx.DiGraph` - with as many nodes as the 20 possible amino acids
            by the length of the sequence.

        .. seealso::
            :meth:`.FragmentFrame.make_frequency_network`
        """
        matrix = self.make_sequence_matrix(frequency=True)

        g = nx.DiGraph()

        for i, row in matrix.iterrows():
            if i == matrix.iloc[0].name:
                nterm = ["0X", ]
            invrow = (1 - row[row > 0])
            cterm = [str(i) + _ for _ in list(invrow.index)]
            for n, c in itertools.product(nterm, cterm):
                g.add_edge(n, c, weight=invrow[c[-1]])
            nterm = cterm
        for n in nterm:
            g.add_edge(n, "-1X", weight=0)
        return g

    def make_frequency_network( self, use_rmsd=False ):
        """Generate a graph with per-position frequency between residue types at each position.

        Generate a :class:`~networkx.DiGraph` in which each node is a residue type and each
        edge the frequency expected for the transition between residue-type in position i and
        residue-type in position i+1.

        :param bool use_rmsd: When :data:`True`, correct the pair counts by the RMSD. Basically,
            the smaller the RMSD, the more the count weights.

        :return: :class:`~networkx.DiGraph` - with as many nodes as the 20 possible amino acids
            by the length of the sequence.

        .. seealso::
            :meth:`.FragmentFrame.make_per_position_frequency_network`
        """
        G = nx.DiGraph()
        min_nodes = []
        max_nodes = []

        def window(seq, n=2):
            it = iter(seq)
            result = tuple(itertools.islice(it, n))
            if len(result) == n:
                yield result
            for elem in it:
                result = result[1:] + (elem,)
                yield result

        data = {}
        for tp, df in self.groupby(["frame", "neighbor"]):
            transitions = list(window(list(df["position"].unique())))
            rmsd = self[self["frame"] == tp[0]]["rmsd"].max() - df["rmsd"].values[0]
            rmsd = int(rmsd * 100)
            for tr in transitions:
                n = df[df["position"] == tr[0]]["aa"].values[0]
                c = df[df["position"] == tr[1]]["aa"].values[0]
                if use_rmsd:
                    data.setdefault(tr[0], []).extend([(n, c), ] * rmsd)
                else:
                    data.setdefault(tr[0], []).append((n, c))

        for k in data:
            options = len(data[k])
            data[k] = Counter(data[k])
            options = data[k].most_common(1)[0][1]
            for p in data[k]:
                n = str(k) + p[0]
                G.add_node(n, order=k, type=p[0])
                c = str(k + 1) + p[1]
                G.add_node(c, order=k + 1, type=p[1])
                G.add_edge(n, c, weight=(options - data[k][p]) / options)
                if k == min(data):
                    min_nodes.append(n)
                if k == max(data):
                    max_nodes.append(c)

        for node in min_nodes:
            G.add_node("0X", order=0, type="X")
            G.add_edge("0X", node, weight=0)
        for node in max_nodes:
            G.add_node("-1X", order=max(data) + 2, type="X")
            G.add_edge(node, "-1X", weight=0)

        return G

    def quick_consensus_sequence( self ):
        """Consensus sequence with the highest representative per position.

        :return: :class:`str` - consensus sequence
        """
        consensus = []
        for i in range(1, max(self["position"].values) + 1):
            values = self[self["position"] == i]["aa"].values
            qseq = sorted(Counter(values).most_common(), key=lambda x: (-x[1], x[0]))[0]
            consensus.append(qseq[0])
        return "".join(consensus)

    def quick_consensus_secondary_structure( self ):
        """Consensus secondary structure with the highest representative per position.

        :return: :class:`str` - consensus secondary structure
        """
        consensus = []
        for i in range(1, max(self["position"].values) + 1):
            values = self[self["position"] == i]["sse"].values
            qseq = sorted(Counter(values).most_common(), key=lambda x: (-x[1], x[0]))[0]
            consensus.append(qseq[0])
        return "".join(consensus)

    #
    # Implement pandas methods
    #
    @property
    def _constructor(self):
        return FragmentFrame

    def __finalize__(self, other, method=None, **kwargs):

        if method == 'merge':
            for name in self._metadata:
                setattr(self, name, getattr(other.left, name, _metadata_defaults(name)))
        else:
            for name in self._metadata:
                setattr(self, name, getattr(other, name, _metadata_defaults(name)))

        return self
