# @Author: Jaume Bonet <bonet>
# @Date:   26-Mar-2018
# @Email:  jaume.bonet@gmail.com
# @Filename: rscripts.py
# @Last modified by:   bonet
# @Last modified time: 06-Apr-2018


import textwrap


def baseline():
    """
    RosettaScript to calculate DSSP secondary structure.

    .. ipython::

        In [1]: from rstoolbox.utils import baseline
           ...: print baseline()

    :return: :class:`str`

    .. seealso::
        :func:`.get_sequence_and_structure`
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
