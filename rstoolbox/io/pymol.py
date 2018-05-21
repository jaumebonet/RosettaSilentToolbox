# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: pymol_mutant_selector
"""
# Standard Libraries

# External Libraries
import pandas as pd

# This Library
import rstoolbox.core as core
import rstoolbox.components as cp


__all__ = ['pymol_mutant_selector']


def pymol_mutant_selector( df ):
    """Generate selectors for the mutations in target decoys.

    Given a :class:`.DesignFrame` with columns ``mutant_positions_<seqID>`` and a
    description column with the name of the decoy, it generates a list of commands
    to select the residues considered as mutants in pymol.

    .. note::
        This function requires that :meth:`.DesignFrame.identify_mutants` has been
        previously run on the data container.

    :param df: |df_param|
    :type df: Union[:class:`.DesignFrame`, :class:`.DesignSeries`]

    :return: :func:`list` of :class:`str` - one command of selection for provided decoy
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
    coms = df.apply(lambda row: make_selector(row["description"], row[headers], chains), axis=1)
    return coms
