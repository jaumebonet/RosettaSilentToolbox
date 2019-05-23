# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: parse_master_file
"""
# Standard Libraries
import os
from ast import literal_eval

# External Libraries
import pandas as pd
import numpy as np

# This Library


def parse_master_file( filename, max_rmsd=None, piece_count=18, shift_0=False ):
    """Load data obtained from a `MASTER <https://grigoryanlab.org/master/>`_ search.

    :param str filename: Output file.
    :param float max_rmsd: Maximum RMSD value to recover.
    :param int piece_count: If known, specify the number of structural pieces
        in the MASTER search. Otherwise, it is assumed.
    :param bool shift_0: MASTER matches start counting in 0. If the flag
        is :data:`.True`, change it to start with 1.

    :return: :class:`~pandas.DataFrame`

    Columns of the returned :class:`~pandas.DataFrame` are:

    =========== ===================================================
    column name description
    =========== ===================================================
    rmsd        RMSD value between query and match
    pds_path    Path to the PDS file containing the match
    pdb         PDB identifier (conditional)
    chain       Chain identifier (conditional)
    match       List of ranges for the match (Rosetta count)
    =========== ===================================================

    The generation of the ``pdb`` and ``chain`` columns assumes that the PDS files
    basename has the standard nomenclature ``<pdbid>_<chain>.pds``. If that is not
    the case, these columns will not be processed.
    """
    def shift(x):
        return np.asarray(np.asarray(x) + 1).tolist()

    df = pd.read_csv(filename,
                     names=list(range(piece_count + 2)), engine='python',
                     sep=r'\s+', header=None).dropna(axis=1, how='all')
    df['match'] = df[df.columns[2:]].astype(str).sum(axis=1).apply(literal_eval)
    if shift_0:
        df['match'] = df['match'].apply(shift)
    df = df.rename(columns={0: 'rmsd', 1: 'pds_path'})
    df = df.drop(columns=[i for i in df.columns if isinstance(i, int)])

    try:
        df[['pdb', 'chain']] = (pd.DataFrame(list(df['pds_path'].str.replace('.pds', '')
                                .apply(lambda x: os.path.basename(x).split('_')).values)))
        rcols = ['rmsd', 'pds_path', 'pdb', 'chain', 'match']
    except Exception:
        rcols = ['rmsd', 'pds_path', 'match']

    if max_rmsd is not None:
        df = df[(df['rmsd'] <= max_rmsd)]
    return df[rcols]
