# @Author: Jaume Bonet <bonet>
# @Date:   09-Jan-2018
# @Email:  jaume.bonet@gmail.com
# @Filename: fragmentFrame.py
# @Last modified by:   bonet
# @Last modified time: 09-Apr-2018

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
import rstoolbox.core as core

__all__ = ["FragmentFrame"]


class FragmentFrame( pd.DataFrame ):
    """
    The :py:class:`.FragmentFrame` extends :py:class:`pandas.DataFrame`
    adding some functionalities in order to facilitate comparissons.
    Filled through the functions provided through this library, each
    row represents a position of a fragment, while columns describe the
    properties.

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
        """
        Adds a source file to the :py:class:`.FragmentFrame`. This can be used to automatically
        generate the fragment RMSD quality.

        :param str file: Name of the file.
        """
        self._source_file = file

    def get_source_file( self ):
        """
        Returns the file name linked to this :py:class:`.FragmentFrame`.

        :return: File name.
        """
        return self._source_file

    def has_source_file( self ):
        """
        Checks if there is a source file added.

        :return: bool
        """
        return self._source_file is not None

    def is_comparable( self, df ):
        """
        Evaluate if the current :py:class:`.FragmentFrame` is comparable to
        the provided one. This is checked on terms of (a) covered sequence
        length -range- and (2) fragment size.

        :return: True if the two :class:`.FragmentFrame` are comparable.
        """
        if max(self["position"]) != max(df["position"]):
            return False
        if min(self["position"]) != min(df["position"]):
            return False
        if self["size"].values[0] != df["size"].values[0]:
            return False
        return True

    def add_quality_measure( self, filename, pdbfile=None ):
        """
        Add to the data the RMSD quality measure provided
        by the r_fraq_qual app from Rosetta.

        :param str filename: Name containing the quality measure.
            If filename is None, it assumes that there is RMSD quality
            yet, so it'll run it as long as the :py:class:`.FragmentFrame`
            has a source file. Standart output will be the name of the source
            file with the extension ".qual". If this file exists, it will be
            automatically picked.
        :param str pdbfile: In case the quality has to be calculated. Provide the
            PDB over which to calculate it. Default is None.
        :raise: IOError if filename does not exist.
        :raise: IOError if pdbfile is provided and does not exist.
        :raise: IOError if the rosetta executable is not found.
            Depends on rosetta.path and rosetta.compilation
        :raise: AttributeError if filename is None and there is no attached
            source file to the object.
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
                exe = os.path.join(core.get_option("rosetta", "path"),
                                   "r_frag_quality.{0}".format(core.get_option("rosetta",
                                                                               "compilation")))
                if not os.path.isfile(exe):
                    raise IOError("The expected Rosetta executable {0} is not found".format(exe))
                if not os.path.isfile(pdbfile):
                    raise IOError("{0} not found".format(pdbfile))
                command = "{0} -in:file:native {1} -f {2} -out:qual {3}".format(
                          exe, pdbfile, self._source_file, filename )
                error = os.system( command )
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

        return self.merge(df, how='left', on=["size", "frame", "neighbor"])

    def select_quantile( self, quantile=0.25 ):
        """
        Returns only the fragments under the rmsd threshold of the specified
        quantile.

        :param quantile: Quantile maximum limit.
        :type quantile: :class:`float`

        :return: :class:`.FragmentFrame` - The filtered data.

        :raises:
            :KeyError: if the `rmsd` column cannot be found.

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

    def make_sequence_matrix( self, frequency=False, round=False ):
        """
        Generate a PSSM-like matrix from the fragments.

        :param frequency: Return the matrix with frequency values? Default will
            return the values as :math:`logodd(f(ni)/f(bi))`.
        :type frequency: :class:`bool`
        :param round: Round-floor the values.
        :type round: :class:`bool`

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
        if round:
            df = df.applymap(np.around).astype(np.int64).reindex(columns=list(alphabet))
        else:
            df = df.reindex(columns=list(alphabet))
        df.index = range(1, df.shape[0] + 1)
        return df

    def make_per_position_frequency_network( self ):
        """
        Generate a :class:`~networkx.DiGraph` in which each node is a residue type and each
        edge the frequency expected for that residue type in that position according to the
        frequency matrix. As a matter of fact, the edge is the *inverted frequency*, as the
        idea is to use the graph to calculate the shortest paths (i.e. the more probable
        sequences).

        :return: :class:`~networkx.DiGraph` - with as many nodes as the 20 possible amino acids
            by the length of the sequence.
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
        """
        Generate a :class:`~networkx.DiGraph` n which each node is a residue type and each
        edge the frequency expected for the transition between residue-type in position i and
        residue-type in position i+1.

        :param use_rmsd: When :data:`True`, correct the pair counts by the RMSD. Basically,
            the smaller the RMSD, the more the count weights.
        :type use_rmsd: :class:`bool`

        :return: :class:`~networkx.DiGraph` - with as many nodes as the 20 possible amino acids
            by the length of the sequence.
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
        """
        Generate a consensus sequence as the most common representative for each position.

        :return: :class:`str` - consensus sequence
        """
        consensus = []
        for i in range(1, max(self["position"].values) + 1):
            values = self[self["position"] == i]["aa"].values
            qseq = sorted(Counter(values).most_common(), key=lambda x: (-x[1], x[0]))[0]
            consensus.append(qseq[0])
        return "".join(consensus)

    def quick_consensus_secondary_structure( self ):
        """
        Generate a consensus secondary structure as the most common representative for each
        position.

        :return: :class:`str` - consensus secondary structure
        """
        consensus = []
        for i in range(1, max(self["position"].values) + 1):
            values = self[self["position"] == i]["sse"].values
            qseq = sorted(Counter(values).most_common(), key=lambda x: (-x[1], x[0]))[0]
            consensus.append(qseq[0])
        return "".join(consensus)

    def angle_coverage( self, df, threshold=5 ):
        """
        Check, for frame (window) how many shared fragments do to fragment set
        have. Two fragments are considered equal if all their angles are under a
        threshold.
        :param FragmentFrame df: :py:class:`.FragmentFrame` to compare to.
        :param float threshold: Degrees of difference allowed.
        :return: List. For each frame the number of shared fragments, the number of
            fragments in that frame for the current :py:class:`.FragmentFrame` and
            for  the provided one.
        :raises: AssertionError if the two sets are not comparable.
        :raises: ValueError if the other DataFrame is not a :py:class:`.FragmentFrame`
        """
        # We assume none of the two are crunched
        if "_anglesAll" not in self:
            assert self.is_comparable(df), "The two sets cannot be compared."
            assert isinstance(df, FragmentFrame), "Compare must be done between FragmentFrames."
            return self._crunch("angles").angle_coverage(df._crunch("angles"), threshold)

        def angle(a1, a2):
            return 180 - abs((abs(a1 - a2) % 360) - 180)

        def compare_angles(a1, a2, thr):
            return not (pd.Series(zip(a1, a2)).apply(lambda row: (angle(row[0], row[1]) <= thr).all()) == False).any()

        def compare_frame( line, df, thr ):
            return df.apply(lambda row: compare_angles(row["_anglesAll"], line, thr), axis=1)

        if self["frame"].nunique() == 1:
            counter = 0
            counter = self.apply(lambda row: compare_frame(row["_anglesAll"], df, threshold), axis=1).apply(lambda row: (row == True).any() ).sum()
            return counter, self.shape[0], df.shape[0]
        else:
            sys.stdout.write("Angle coverage: This process is computationaly expensive. Will take a bit.")
            data = []
            for i in self["frame"].unique():
                data.append(self[self["frame"] == i].angle_coverage(df[df["frame"] == i]))
            return data

    def _crunch( self, what ):
        """
        Crunching is a process by which the :py:class:`.FragmentFrame` is compressed
        so flatten the requested values
        :param str what: what to crunch: angles, seq, sse
        :return: The crunched :py:class:`.FragmentFrame`.
        """
        df = self.copy()
        if what in self._crunched:
            return self._crunched[what]
        if what.lower() == "angles":
            df["_anglesAll"] = df.apply(lambda row: list(row[["phi", "psi", "omega"]]), axis=1)
            df = df.drop(["phi", "psi", "omega"], axis=1)
            df = FragmentFrame(df.groupby(["frame", "neighbor"])["_anglesAll"].apply(np.hstack).reset_index())

        self._crunched[what] = df
        return self._crunched[what]

    def _metadata_defaults(self, name):
        if name == "_source_file":
            return None
        return None

    #
    # Implement pandas methods
    #
    @property
    def _constructor(self):
        return FragmentFrame

    def __finalize__(self, other, method=None, **kwargs):

        if method == 'merge':
            for name in self._metadata:
                setattr(self, name, getattr(other.left, name, self._metadata_defaults(name)))
        else:
            for name in self._metadata:
                setattr(self, name, getattr(other, name, self._metadata_defaults(name)))

        return self
