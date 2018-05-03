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

__all__ = ['baseline']


def baseline():
    """RosettaScript to calculate DSSP secondary structure and
    phi-psi angles.

    :return: :class:`str`

    .. seealso::
        :func:`.get_sequence_and_structure`

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.utils import baseline
           ...: print baseline()
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
