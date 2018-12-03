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


def load_refdata( ref, homology=None ):
    """Load the predefined reference data from ``cath``, ``scop``, ``scop2`` or
    ``chain``.

    Learn in the tutorials more about how to :ref:`contextualize your data
    <contextualize>`.

    :param str ref: Reference data to load.
    :param float homology: Allowed redundancy threshold. Derived from the PDB's
        pre-calculated homologies. Thus, available similarities are:
        [``30``, ``40``, ``50``, ``70``, ``90``, ``95``, ``100``].

    :return: :class:`~pandas.DataFrame`

    :raises:
        :ValueError: If an unknown ``reference`` is requested.
        :ValueError: If an unknown ``homology`` is requested.
    """
    if ref.lower() not in ['cath', 'scop', 'scop2', 'chain']:
        raise ValueError('Unknown reference requested.')
    if homology is not None and homology not in [30, 40, 50, 70, 90, 95, 100]:
        raise ValueError('Unknown homology threshold.')

    cwd = os.path.dirname(os.path.abspath(__file__))
    df = pd.read_csv(os.path.join(cwd, 'baselines', '{}.gz'.format(ref.lower())))

    df['selectors'] = df['selectors'].apply(literal_eval)
    if 'phi_$' in df:
        df['phi_$'] = df['phi_$'].apply(literal_eval)
    if 'psi_$' in df:
        df['psi_$'] = df['psi_$'].apply(literal_eval)
    if 'node_id' in df:
        df['node_id'] = df['node_id'].apply(literal_eval)
    if 'node_name' in df:
        df['node_name'] = df['node_name'].apply(literal_eval)

    if homology is not None:
        if ref == 'cath':
            df['pdb'] = df['pdb'].str.upper()
        hmdf = pd.read_csv(os.path.join(cwd, 'baselines', 'redundancies.gz'))
        df = df.merge(hmdf, on=['pdb', 'chain'])
        df = (df.sort_values('score')
                .groupby('c{}'.format(homology))
                .first().reset_index()
                .drop(columns=['c30', 'c40', 'c50', 'c70', 'c90', 'c95', 'c100']))
    return df
