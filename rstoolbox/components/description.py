# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

One of the key advantadges of ``rstoolbox`` is the ability to control the amount
and type of data that is loaded from a silent/score file. This control is managed
through a ``definition``, a dictionary that describes the type of data that can be
loaded.

.. note::

    ``definition`` is meant to be applied to :func:`.parse_rosetta_file`.

As of now, there are 10 different options that can be convined into a ``definition``:

========================= ===================================================
definition term           description
========================= ===================================================
:ref:`scores`             Basic selection of the scores to store. Default is
                          all scores.
:ref:`scores_ignore`      Selection of specific scores to ignore.
:ref:`scores_rename`      Rename some score names to others.
:ref:`scores_by_residue`  Pick score by residue types into a single array
                          value.
:ref:`naming`             Use the decoy identifier's name to create extra
                          score terms.
:ref:`sequence`           Pick sequence data from the silent file.
:ref:`structure`          Pick structural data from the silent file.
:ref:`psipred`            Pick PSIPRED data from the silent file.
:ref:`dihedrals`          Retrieve dihedral data from the silent file.
:ref:`labels`             Retrieve residue labels from the silent file.
========================= ===================================================

.. tip::

    ``definition`` can be passed directly as a dictionary or can be saved as a
    **JSON** or **YAML** file and loaded from there.

.. _scores:

scores
------

This is the most basic parameter, and refer to regular scores in the silent/score
file. It allows to select just the scores that are wanted for the analysis.
There are three main ways to define ``scores``, provide a list naming the
scores of interest::

    {'scores': ['score', 'packstat', 'description']}

add a string asterisc if all scores all wanted (this is the default value for this
parameter)::

    {'scores': '*'}

or add a minus sign, which will ignore all scores::

    {'scores': '-'}

.. rubric:: Example

.. ipython::

    In [1]: from rstoolbox.io import parse_rosetta_file
       ...: import pandas as pd
       ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
       ...: ifile = '../rstoolbox/tests/data/input_2seq.minisilent.gz'
       ...: definition1 = {'scores': ['score', 'packstat', 'description']}
       ...: df = parse_rosetta_file(ifile, definition1)
       ...: df.head()

    In [2]: definition2 = {'scores': '*'}
       ...: df1 = parse_rosetta_file(ifile, definition2)
       ...: df2 = parse_rosetta_file(ifile)
       ...: df1.head()

    In [3]: (df1.columns == df2.columns).all()

    In [4]: definition3 = {'scores': '-'}
       ...: df3 = parse_rosetta_file(ifile, definition3)
       ...: df3.head()

.. _scores_ignore:

scores_ignore
-------------

This is basically the oposite from the previous one, meant in case particular scores are
to be ignored. This is usefull when expecting to mix data from different experiments but
one has generated more score data than the other (for example, one needed loop closure and
includes some ``loop_closure`` scores). Extra scores can affect the concatenation of the
different data containers as per :mod:`pandas` constraints. There are two ways to define
``scores_ignore``, either list the scores to skip::

    {'scores_ignore': ['loop_closure', 'packstat']}

or a string asterisc if aiming to ignore all scores::

    {'scores_ignore': '*'}

.. _scores_rename:

scores_rename
-------------

Allows to retrieve a particular score term with a different name. Again, usefull to
merge data from multiple runs when some naming is not matching. It is defined as a
dictionary in which keys are the original score names and values the new naming
schema::

    {'scores_rename': {'grmsd': 'globalRMSD', 'lrmsd': 'localRMSD'}}

.. _scores_by_residue:

scores_by_residue
-----------------

Scores that target per-residue values are normally ignored by the library. This parameter
allows for them to be captured as a single vector score. As of now, the only available
score is ``residue_ddg_``::

    {'scores_by_residue': ['residue_ddg_']}

as provided by **Rosetta**'s' ``ddG`` mover::

    <ddG name="(&string)" per_residue_ddg="1" "/>

.. _naming:

naming
------

Naming conventions are important. As long as one keeps that in mind, one could target
particular identifiers inside a decoy description as new score columns, thus allowing
to cluster different decoys for analysis. As an example, let's assume that the naming
of a set of decoys is such as::

    nubinitio_auto_2015_binder_2pw9C_0001

In this case, three of the forut first elements of the identifier relate to the
conditions of the experiment. To capture them, an array needs to be defined. As long
as an identifier is set for an element, that element will be captured as a new score.
Elements can be skiped with an empty string. The array has to be long enough to capture
all the elements of interest, but it does not need to have as many fields as elements
in the definition. Thus::

    {'naming': ['experiment', 'fragments', '', 'binder']}

will create the three new scores with value ``nubinitio``, ``auto``, ``binder``.

.. _sequence:

sequence
--------

Sequence data is integrated into the silent file as default. It can be captured
and used to evaluate mutants and sequence drift amongst others. To retrieve that
data, a string with the identifiers of all the chains of interest need to be
provided, thus::

    {'sequence': 'AB'}

will allow to capture the sequence for chains ``A`` and ``B``. Sequence data can
then be accessed through the appropiate getter functions of :class:`.DesignFrame`
and :class:`.DesignSeries`, and is expected for any of the sequence analysis
functions and plots.

.. _structure:

structure
---------

Secondary structure data can be loaded into the silent file by means of **Rosetta**'s
``WriteSSEMover``::

    <WriteSSEMover name="(&string;)" dssp="1" />

Similarly to ``sequence``, it can be loaded as extra data by calling the chains of
interest::

    {'structure': 'AB'}

.. _psipred:

psipred
-------

Secondary structure prediction data can be loaded into the silent file by means
of **Rosetta**'s ``WriteSSEMover``::

    <WriteSSEMover name="(&string;)" cmd="/path/to/psipred" />

Similarly to ``sequence`` and ``structure``, it can be loaded as extra data by
calling the chains of interest::

    {'psipred': 'AB'}

.. _dihedrals:

dihedrals
---------

Phi and psi dihedral angle data can be loaded into the silent file by means
of **Rosetta**'s ``WriteSSEMover``::

    <WriteSSEMover name="(&string;)" write_phipsi="1" />

Angle data from specific chains can be loaded as::

    {'dihedrals': 'AB'}

The data will be loaded as a single array of floats.

.. _labels:

labels
------

Residue labels allow to target residues through a simulation by a particular tag.
They can also be used to highlight residues with particular properties during the
simulation. They can be saved into the silent file with **Rosetta**'s
``DisplayPoseLabelsMover``::

    <DisplayPoseLabelsMover name="(&string;)" write="1" />

To retrieve that data, a list with the names of the labels of interest can be
provided::

    {'labels': ['MOTIF', 'CONTEXT', 'HOTSPOTS']}

"""
# Standard Libraries
import re

# External Libraries
import six
import numpy as np

# This Library

__all__ = ["Description"]


class Description( object ):

    _PARAMS = [
        "scores", "scores_ignore", "scores_rename", "scores_by_residue",
        "naming", "sequence", "structure", "psipred", "dihedrals", "labels"
    ]

    def __init__( self, **params ):
        # Avoid false attributes and misspellings
        for k in params:
            if k not in self._PARAMS:
                raise AttributeError(
                    "Wrong attribute provided. Accepted keys are: "
                    "{}".format(",".join(self._PARAMS)))

        # Load selections
        self.scores            = params.get("scores",            "*" )
        self.scores_ignore     = params.get("scores_ignore",     None)
        self.scores_rename     = params.get("scores_rename",     None)
        self.scores_by_residue = params.get("scores_by_residue", None)
        self.naming            = params.get("naming",            None)
        self.sequence          = params.get("sequence",          None)
        self.structure         = params.get("structure",         None)
        self.psipred           = params.get("psipred",           None)
        self.dihedrals         = params.get("dihedrals",         None)
        self.labels            = params.get("labels",            None)

        # Others
        self._per_residues = ["residue_ddg_"]
        self._scores       = {}

        # Check that the "per_residue" scores requested are known.
        if self.scores_by_residue is not None:
            for k in self.scores_by_residue:
                if k not in self._per_residues:
                    raise AttributeError("Unknown per-residue score: {}".format(k))

        # Fix user input in case strings instead of lists were provided
        if isinstance(self.scores, six.string_types) and self.scores not in ["*", "-"]:
            self.scores = self.scores.split(",")

        # Negate all scores
        if self.scores_ignore == "*" or self.scores == "-":
            self.scores = None

        # Labels will be set in UPPERCASE
        if self.labels is not None:
            self.labels = [x.upper() for x in self.labels]

    def wanted_score( self, score_name ):
        # skip per-residue values when nod asked for
        for k in self._per_residues:
            if score_name.startswith(k) and score_name not in self.scores:
                return False
        if self.scores_ignore is not None:
            if self.scores_ignore == "*" or score_name in self.scores_ignore:
                return False
            elif any(re.match(s, score_name) for s in self.scores_ignore if "*" in s):
                return False
        if self.scores_rename is not None:
            if score_name in self.scores_rename:
                return True
        if self.scores is None:
            return False
        if self.scores == "*" or score_name in self.scores:
            return True
        return False

    def score_name( self, score_name ):
        if self.scores_rename is not None:
            if score_name in self.scores_rename:
                return self.scores_rename[score_name]
        return score_name

    def wanted_per_residue_score( self, score_name ):
        if self.scores_by_residue is None:
            return False
        for k in self._per_residues:
            if score_name.startswith(k):
                return True

    def check_naming( self, header ):
        if self.naming is None:
            return
        for h in self.naming:
            if h in header and self.wanted_score(h):
                if h == self.score_name(h):
                    raise AttributeError(
                        "New column {} ".format(h) +
                        "cannot be created, a score has that name")

    def get_naming_pairs( self, description ):
        if self.naming is not None:
            d = description.split("_")
            for cn, n in enumerate(d):
                if cn == len(self.naming):
                    break
                if self.naming[cn] != "":
                    yield self.naming[cn], n

    def setup_labels( self, data ):
        if self.labels is not None:
            for l in self.labels:
                data.setdefault("lbl_" + l.upper(), []).append("")
        return data

    def get_expected_sequences( self, chains ):
        if self.sequence is not None:
            sele = set(chains["id"]) if self.sequence == "*" else set(self.sequence)
            seq  = chains["seq"]
            if len(sele.difference(chains["id"])) > 0:
                raise ValueError(
                    "Requested a chain not present in the file. "
                    "Available chain are {}".format(",".join(list(set(chains["id"])))))
            guide  = np.array(list(chains["id"]))
            for ch in sele:
                guidep = np.where(guide == ch)[0]
                yield "sequence_" + ch, "".join(seq[guidep[0]:guidep[-1] + 1])

    def get_expected_structures( self, chains ):
        if self.structure is not None:
            sele = set(chains["id"]) if self.structure == "*" else set(self.structure)
            seq  = chains["dssp"]
            if len(sele.difference(chains["id"])) > 0:
                raise ValueError(
                    "Requested a chain not present in the file. "
                    "Available chain are {}".format(",".join(list(set(chains["id"])))))
            guide  = np.array(list(chains["id"]))
            for ch in sele:
                guidep = np.where(guide == ch)[0]
                yield "structure_" + ch, "".join(seq[guidep[0]:guidep[-1] + 1])

    def get_expected_psipred( self, chains ):
        if self.psipred is not None:
            sele = set(chains["id"]) if self.psipred == "*" else set(self.psipred)
            seq  = chains["psipred"]
            if len(sele.difference(chains["id"])) > 0:
                raise ValueError(
                    "Requested a chain not present in the file. "
                    "Available chain are {}".format(",".join(list(set(chains["id"])))))
            guide  = np.array(list(chains["id"]))
            for ch in sele:
                guidep = np.where(guide == ch)[0]
                yield "psipred_" + ch, "".join(seq[guidep[0]:guidep[-1] + 1])

    def get_expected_dihedrals( self, chains, angle ):
        if self.dihedrals is not None:
            sele = set(chains["id"]) if self.dihedrals == "*" else set(self.dihedrals)
            seq  = chains[angle]
            if len(sele.difference(chains["id"])) > 0:
                raise ValueError(
                    "Requested a chain not present in the file. "
                    "Available chain are {}".format(",".join(list(set(chains["id"])))))
            guide  = np.array(list(chains["id"]))
            for ch in sele:
                guidep = np.where(guide == ch)[0]
                yield angle + "_" + ch, np.array(seq[guidep[0]:guidep[-1] + 1])

    def to_json( self ):
        data = {}
        for a in self._PARAMS:
            if getattr(self, a) is not None:
                data[a] = getattr(self, a)
        return data
