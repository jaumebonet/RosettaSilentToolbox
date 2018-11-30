# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: load_refdata
"""
# Standard Libraries
import os
from ast import literal_eval

# External Libraries
import pandas as pd

# This Library


__all__ = ['load_refdata']


def load_refdata( ref ):
    """Load the reference data from ``cath``, ``scop``, ``scop2`` or
    ``chain``.

    :param str ref: Reference data to load.

    :return: :class:`~pandas.DataFrame`

    :raises:
        :ValueError: It an unknown reference is requested.
    """
    if ref.lower() not in ['cath', 'scop', 'scop2', 'chain']:
        raise ValueError('Unknown reference requested.')

    cwd = os.path.dirname(os.path.abspath(__file__))
    df = pd.read_csv(os.path.join(cwd, 'baseline', '{}.gz'.format(ref.lower())))

    df['selectors'] = df['selectors'].apply(literal_eval)
    df['phi_$'] = df['phi_$'].apply(literal_eval)
    df['psi_$'] = df['psi_$'].apply(literal_eval)
    df['node_id'] = df['node_id'].apply(literal_eval)
    df['node_name'] = df['node_name'].apply(literal_eval)

    return df
