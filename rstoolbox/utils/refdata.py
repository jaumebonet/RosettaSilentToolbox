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
from ftplib import FTP  # nosec
from ast import literal_eval
from io import BytesIO

# External Libraries
import pandas as pd

# This Library

__all__ = ['load_refdata', 'make_redundancy_table']


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
        hmdf = make_redundancy_table(precalculated=True)
        df = df.merge(hmdf, on=['pdb', 'chain'])
        df = (df.sort_values('score')
                .groupby('c{}'.format(homology))
                .first().reset_index()
                .drop(columns=['c30', 'c40', 'c50', 'c70', 'c90', 'c95', 'c100']))
    return df


def make_redundancy_table( precalculated=False, select=None ):
    """Query into the PDB to retrieve the pre-calculated homology tables.

    .. warning::
        This function requires an active Internet connection.

    :param bool precalculated: When :data:`True`, provide the homology table provided
        by the library; otherwise, query on PDB.
    :param select: List of homology thresholds to include. By default is :data:`None`.
        This option is ignored if ``precalculated`` is :data:`True`.
    :type select: :func:`list` of :class:`int`

    :return: :class:`~pandas.DataFrame`

    :raises:
        :ValueError: When trying to select an homology threshold not provided by the PDB
            database.
    """
    if precalculated:
        cwd = os.path.dirname(os.path.abspath(__file__))
        return pd.read_csv(os.path.join(cwd, 'baselines', 'redundancies.gz'))

    known_hom = [30, 40, 50, 70, 90, 95, 100]
    if select is not None:
        for s in select:
            if s not in known_hom:
                wrn1 = 'Homology threshold {} not provided by the PDB. '.format(s)
                wrn2 = 'Available options are: {}'.format(','.join(known_hom))
                raise ValueError(wrn1 + wrn2)

    all_ = []
    ftp = FTP('resources.rcsb.org')  # nosec
    ftp.login()
    for ident in known_hom:
        if select is not None and ident not in select:
            continue
        r = BytesIO()
        ftp.retrbinary('RETR /sequence/clusters/bc-{}.out'.format(ident), r.write)
        idata = [x.strip() for x in r.getvalue().decode('UTF-8').split('\n')]
        hmdf = []
        for i, line in enumerate(idata):
            ldata = []
            ld = line.split()
            for cnt in ld:
                pdb = cnt.split('_')
                ldata.append([pdb[0], pdb[1], i])
            hmdf.append(pd.DataFrame(ldata, columns=['pdb', 'chain', 'c{}'.format(ident)]))
        all_.append(pd.concat(hmdf))

    if len(all_) == 0:
        return pd.DataFrame()
    if len(all_) == 1:
        return all_[0]

    data = all_[0]
    for i in range(1, len(all_)):
        data = data.merge(all_[i], on=['pdb', 'chain'])
    return data
