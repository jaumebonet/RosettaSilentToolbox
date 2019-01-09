# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: baseline
"""
# Standard Libraries
import textwrap

# External Libraries

# This Library

__all__ = ['baseline', 'mutations']


def baseline( minimize=False ):  # pragma: no cover
    """RosettaScript to calculate DSSP secondary structure and
    phi-psi angles.

    :param bool minimize: If :data:`True`, apply minimization before
        evaluating the structure.

    :return: :class:`str`

    .. seealso::
        :func:`.get_sequence_and_structure`

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.utils import baseline
           ...: print(baseline())
    """
    mover = '\n\t<Add mover="m" />' if minimize else ''
    return textwrap.dedent("""\
    <ROSETTASCRIPTS>
        <MOVERS>
            <WriteSSEMover dssp="1" name="w" write_phipsi="true" />
            <MinMover name="m" bb="true" chi="true" tolerance=".1" />
        </MOVERS>
        <FILTERS>
            <CavityVolume name="cavity" confidence="0" />
            <PackStat name="pack" confidence="0" />
            <AverageDegree name="avdegree" confidence="0" />
        </FILTERS>
        <PROTOCOLS>{}
            <Add mover="w" />
            <Add filter="cavity" />
            <Add filter="pack" />
            <Add filter="avdegree" />
        </PROTOCOLS>
    </ROSETTASCRIPTS>
    """).format(mover)


def mutations( seqID='A' ):  # pragma: no cover
    """RosettaScript to execute a
    `RESFILE <https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/resfiles>`_.

    :param str seqID: |seqID_param|

    :return: :class:`str`

    .. seealso::
        :func:`.DesignFrame.apply_resfile`

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.utils import mutations
           ...: print(mutations())
    """
    return textwrap.dedent("""\
    <ROSETTASCRIPTS>
        <TASKOPERATIONS>
            <ReadResfile name="targets" filename="%%resfile%%"/>
        </TASKOPERATIONS>
        <MOVERS>
            <PackRotamersMover name="packrot" task_operations="targets" />
            <AddJobPairData name="annotate" value_type="string"
                            key="resfile_{}" value="%%resfile%%" />
        </MOVERS>
        <PROTOCOLS>
            <Add mover="packrot" />
            <Add mover="annotate" />
        </PROTOCOLS>
    </ROSETTASCRIPTS>
    """).format(seqID)
