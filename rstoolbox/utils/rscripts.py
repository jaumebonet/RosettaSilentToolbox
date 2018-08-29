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


def baseline():
    """RosettaScript to calculate DSSP secondary structure and
    phi-psi angles.

    :return: :class:`str`

    .. seealso::
        :func:`.get_sequence_and_structure`

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.utils import baseline
           ...: print(baseline())
    """
    return textwrap.dedent("""\
    <ROSETTASCRIPTS>
        <MOVERS>
            <WriteSSEMover dssp="1" name="w" write_phipsi="true" />
        </MOVERS>
        <PROTOCOLS>
            <Add mover="w" />
        </PROTOCOLS>
    </ROSETTASCRIPTS>
    """)


def mutations( seqID='A' ):
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
