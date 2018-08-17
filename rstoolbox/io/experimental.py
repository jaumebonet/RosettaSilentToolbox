# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. codeauthor:: Fabian Sesterhenn <sesterhenn.fabian@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: read_CD
.. func:: read_SPR
.. func:: read_MALS
.. func:: read_fastq
"""
# Standard Libraries
import gzip
import os

# External Libraries
import pandas as pd
import numpy as np
from six import StringIO

# This Library
import rstoolbox.components as rc
import rstoolbox.utils as ru

__all__ = ['read_CD', 'read_SPR', 'read_MALS', 'read_fastq']


class CDFrame( pd.DataFrame ):
    _subtyp = 'cd_frame'

    @property
    def _constructor(self):
        return CDFrame


class SPRFrame( pd.DataFrame ):
    _subtyp = 'spr_frame'

    @property
    def _constructor(self):
        return SPRFrame


def read_CD( dirname, prefix=None, invert_temp=False, model='J-815'):
    """Read `Circular Dichroism <https://www.wikiwand.com/en/Circular_dichroism>`_ data
    for multiple temperatures.

    Assumes that all files of a experiment are in the same folder.

    .. tip::
        The provided :class:`~pandas.DataFrame` will actually be casted to **CDFrame**.
        This class has no other purpose that identification to help other library functions
        and, thus, works as a normal :class:`~pandas.DataFrame`.

    :param str dirname: Folder containing all the CD files.
    :param str prefix: Prefix of the files inside the folder.
    :param bool invert_temp: Temperature assignation might be inverted. Switch it with this
        option.
    :param str model: Format of the data. Available models are: ``['J-815']``

    :return:  :class:`~pandas.DataFrame`.

    :raise:
        :IOError: If ``dirname`` is not a directory.
        :ValueError: If there are inconsistencies between the files (title or temperatures).
        :AttributeError: Unknown ``model`` provided.

    .. seealso::
        :func:`.plot_CD`

    """
    if model == 'J-815':
        return CDFrame(_read_CD_J815(dirname, prefix, invert_temp))
    else:
        raise AttributeError('Unknown CD format')


def _read_CD_J815( dirname, prefix, invert_temp):
    """CD read method for the J-815 output format.

    .. seealso::
        :func:`.read_CD`
    """
    df = []
    if not os.path.isdir(dirname):
        raise IOError('Directory {} not found'.format(dirname))
    for root, _, files in os.walk(dirname):
        temp = [-100, -100]
        title = None
        counter = 0
        for f in files:
            if prefix is not None and not f.startswith(prefix):
                continue
            if not f.startswith('.'):
                with open(os.path.join(root, f)) as fd:
                    data = []
                    inforead = True
                    for ln in fd:
                        if ln.startswith(('XYDATA', '##### Extended Information')):
                            inforead = not inforead
                            continue
                        if not inforead:
                            data.append(ln)
                        else:
                            if ln.startswith(('End Temperature.', 'Start Temperature.')):
                                i = 0 if ln.startswith('End Temperature.') else 1
                                tmp = float(ln.split()[2])
                                if temp[i] == -100:
                                    temp[i] = tmp
                                elif temp[i] != tmp:
                                    raise ValueError('Temp range between files does not match')
                            if ln.startswith('TITLE'):
                                tl = ln.split()[-1].strip()
                                if title is None:
                                    title = tl
                                else:
                                    if tl != title:
                                        raise ValueError('Experiment titles do not match')
                df.append(pd.read_csv(StringIO("".join(data)), header=None,
                                      names=['Wavelength', 'MRE', 'voltage'], sep='\t'))

                df[-1] = ru.add_column(df=df[-1], name='title', value=title)
                counter += 1
                df[-1] = ru.add_column(df=df[-1], name='bin', value=counter)
                minw = df[-1][(df[-1]['MRE'] == df[-1]['MRE'].min())].iloc[0]['Wavelength']
                df[-1] = ru.add_column(df=df[-1], name='minw', value=minw)
    df = pd.concat(df, ignore_index=True)

    bins = len(df['bin'].unique())
    if invert_temp:
        temp.reverse()
    splits = round(float(temp[1] - temp[0]) / bins)
    if splits != 0:
        tmp = pd.DataFrame([(x + 1, y) for x, y in enumerate(np.arange(temp[0],
                                                                       temp[1] + splits,
                                                                       (splits)))],
                           columns=['bin', 'Temp'])
        df = df.merge(tmp, on=['bin'])
    return df


def read_SPR( filename ):
    """Reads **Surface Plasmon Resonance** data.

    The input data should be a comma-separated file with two types of header, one for raw data::

        Run 1; Ch 1; Cy 6; RefSub; name; conc=0.0585_X ,...

    And one for fitted data::

        Run 1; Ch 1; Cy 16; RefSub; 1kx8_d02; conc=60 fitted curve_X ,...

    each ``conc`` condition will have a ``_X`` and ``_Y`` value.

    This seems to be a pretty standard format as output from the machine.

    .. tip::
        The provided :class:`~pandas.DataFrame` will actually be casted to **SPRFrame**.
        This class has no other purpose that identification to help other library functions
        and, thus, works as a normal :class:`~pandas.DataFrame`.

    :param str filename: Input file.

    :return: :class:`~pandas.DataFrame` with ``MultiIndex``.

    .. seealso::
        :func:`.plot_SPR`

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import read_SPR
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = read_SPR("../rstoolbox/tests/data/spr_data.csv.gz")
           ...: df.head(2)
    """
    df = pd.read_csv(filename, sep='\t', low_memory=False).dropna(how='all')
    if len(df.columns) == 1:
        df = pd.read_csv(filename, low_memory=False).dropna(how='all')
    df.columns = [_.split(';')[-1].replace('conc=', '') for _ in df.columns]
    condition = ['raw' if 'fitted' not in _ else 'fitted' for _ in df.columns]
    concentration = [_.split('_')[0].replace('fitted curve', '').strip() for _ in df.columns]
    axis = [_.split('_')[-1].strip() for _ in df.columns]
    newidx = [(x, y, z) for x, y, z in zip(condition, concentration, axis)]
    df.columns = pd.MultiIndex.from_tuples(newidx)
    df.columns.names = ['data', 'concentration', 'axis']
    df.index = range(0, df.shape[0])
    return SPRFrame( df.astype('float64') )


def read_MALS( filename, mmfile=None ):
    """Read data from
    `Multi-Angle Light Scattering <https://www.wikiwand.com/en/Multiangle_light_scattering>`_
    data.

    It is assumed that the format of this file would be 4 columns alternating time and
    LV and UV. MW can also be present.

    :param str filename: Input file.
    :param str mmfile: Molecular Mass file.

    :return: :class:`~pandas.DataFrame` with ``MultiIndex``.

    .. seealso::
        :func:`.plot_MALS`
    """
    df = pd.read_csv(filename)
    if df.shape[1] == 1:
        df = pd.read_csv(filename, sep='\t')
    df1 = df.T.iloc[:2].T
    df1.columns = [x.replace('.1', '') for x in df1.columns]
    df2 = df.T.iloc[2:].T
    df2.columns = [x.replace('.1', '') for x in df2.columns]
    df = pd.concat([df1, df2], sort=False)
    df = df.rename(columns={'time (min)': 'Time'})
    if mmfile is not None:
        mm = pd.read_csv(mmfile)
        mm = mm.T.iloc[:2].T
        mm = mm.rename(columns={'time (min)': 'Time',
                                'Molar Mass 1': 'MW',
                                'Molar Mass 2': 'MW'})
        df = pd.concat([df, mm], sort=False)
    df.index = list(range(df.shape[0]))
    return df


def read_fastq( filename, seqID='A'):
    """Reads a FASTQ file and stores the ID together with the sequence.

    The default generated :class:`.DesignFrame` will contain two columns:

    ====================  ===================================================
    Column Name            Data Content
    ====================  ===================================================
    **description**        Sequence identifier.
    **sequence_<chain>**   Sequence content.
    ====================  ===================================================

    :param str filename: FASTQ filename.
    :param str seqID: |seqID_param|

    :return: :class:`.DesignFrame`

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import read_fastq
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = read_fastq("../rstoolbox/tests/data/cdk2_rand_001.fasq.gz")
           ...: df.head(8)
    """
    # Empty array to store tuples of ID & sequence information
    fastq = []
    idq = []

    # Create a file handle for parsing
    is_gz = filename.endswith('gz')
    fastq_file = gzip.open(filename) if is_gz else open(filename)
    for line in fastq_file:
        line = line.decode('utf8') if is_gz else line
        if line.startswith('@'):
            idq.append(str(line.split(':')[0].split(';')[0][1:]))
        if '@' in line or '+' in line or any(c.islower() for c in line):
            continue
        if len(line) == 0:
            continue
        fastq.append(str(line.strip()))
    return rc.DesignFrame({'description': idq, 'sequence_{}'.format(seqID): fastq})
