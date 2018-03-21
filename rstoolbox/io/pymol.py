import pandas as pd

import rstoolbox.core as core
import rstoolbox.components as cp


def pymol_mutant_selector( df ):
    """
    Given a :class:`.DesignFrame` with columns mutant_positions_<seqID> and a
    description column with the name of the decoy, it generates a list of commands
    to select the residues considered as mutants in pymol.

    :param df:
    :type df: :class:`.DesignFrame`

    :return: :func:`list` of :class:`strp` - one command of selection for provided decoy
    """

    def make_selector( name, row, chains ):
        comm = "sele {0}_mut, {0} ".format(name)
        sels = []
        for i, h in enumerate(list(row)):
            sele = "(c. {} and (".format(chains[i])
            sele += " or ".join(["i. {}".format(n) for n in h.split(",")])
            sele += "))"
            sels.append(sele)
        comm += "and (" + " or ".join(sels) + ")"
        return comm

    headers = [x for x in list(df) if x.startswith("mutant_positions_")]
    chains  = [h[-1] for h in headers]
    coms = df.apply( lambda row: make_selector(row["description"], row[headers], chains), axis=1)
    return coms
