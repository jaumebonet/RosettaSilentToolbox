# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: format_Ipython
.. func:: add_column
.. func:: split_values
"""
# Standard Libraries
import copy
import textwrap

# External Libraries
import pandas as pd

# This Library


__all__ = ['format_Ipython', 'add_column', 'split_values']


def format_Ipython():
    """Ensure ``monospace`` representation of :class:`~pandas.DataFrame`
    in **Jupyter Notebooks**.

    Just need to call it after importing the library.

    .. note::
        In order for this function to work, it is important that is the last
        one in the Jupyter cell to be called.
    """
    pd.set_option("display.max_columns", None)
    pd.set_option("display.max_rows", None)
    pd.set_option("display.max_seq_items", 3)
    pd.set_option("display.max_colwidth", -1)
    from IPython.core.display import HTML
    CSS = textwrap.dedent("""
        table.dataframe {
            font-family: monospace;
        }
    """)
    return HTML('<style>{}</style>'.format(CSS))


def add_column( df, name, value ):
    """Adds a new column to the DataFrame with the given value.

    :param df: Data container.
    :type df: :class:`~pandas.DataFrame`
    :param str name: Name of the new column
    :param value: Value that will be given to all rows of the new column (any type)

    :return: :class:`~pandas.DataFrame` - The data container with the new column
    """
    data = pd.Series([value] * df.shape[0])
    return df.assign(_placeholder=data).rename(columns={"_placeholder": name})


def split_values( df, keys ):
    """Reshape the data to aide plotting of multiple comparable scores.

    .. note::
        This might change the data in a way that a decoy would be repeated
        multiple times.

    The dictionary that needs to be provided to split the data container has three
    main keys:

    #. ``keep``: Identity the columns to keep (they cannot be the ones that split). \
        If not provided, all columns are kept.
    #. ``split``: List with columns to split. Each position is a tuple. The first position \
        is the name of the column to split and the rest will be the value names that will be \
        used to identify it.
    #. ``names``: Names of the columns. The first one will be the name of the column where the \
        values will be assigned, the rest will be the names of the columns for the rest of the \
        identifiers.

    :param df: Data container.
    :type df: :class:`~pandas.DataFrame`
    :param dict keys: Selection of the columns to keep and split.

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: from rstoolbox.utils import split_values
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: ifile = '../rstoolbox/tests/data/input_2seq.minisilent.gz'
           ...: scorel = ['score', 'GRMSD2Target', 'GRMSD2Template', 'LRMSD2Target',
           ...:           'LRMSDH2Target', 'LRMSDLH2Target', 'description']
           ...: df = parse_rosetta_file(ifile, {'scores': scorel})
           ...: df

        In [2]: split1 = {'split': [('GRMSD2Target', 'grmsdTr'), ('GRMSD2Template', 'grmsdTp'),
           ...:                     ('LRMSD2Target', 'lrmsdTp'), ('LRMSDH2Target', 'lrmsdh2'),
           ...:                     ('LRMSDLH2Target', 'lrmsdlh2')],
           ...: 'names': ['rmsd', 'rmsd_type']}
           ...: split_values(df, split1)

        In [3]: split2 = {'split': [('GRMSD2Target', 'global', 'target'),
           ...:                     ('GRMSD2Template', 'global', 'template'),
           ...:                     ('LRMSD2Target', 'local', 'target'),
           ...:                     ('LRMSDH2Target', 'local', 'helix2'),
           ...:                     ('LRMSDLH2Target', 'local', 'lhelix2')],
           ...: 'names': ['rmsd', 'rmsd_type', 'rmsd_target']}
           ...: split_values(df, split2)
    """
    split_columns = [_[0] for _ in keys['split']]
    if 'keep' not in keys:
        keys.setdefault('keep', list(set(df.columns).difference(set(split_columns))))
        keys['keep'].sort(key=lambda x: list(df.columns.values).index(x))
    dataframes = []
    for k in keys["split"]:
        colIDs = copy.copy(keys["keep"])
        colIDs.append(k[0])
        wdf = df[colIDs]
        wdf = wdf.assign(tmpkey1=pd.Series([k[1]] * len(wdf[colIDs[0]])).values).copy(True)
        wdf = wdf.rename(index=str, columns={
            k[0]: keys["names"][0],
            "tmpkey1": keys["names"][1]
        })
        if ( len(k) > 2 ):
            wdf = wdf.assign(tmpkey2=pd.Series([k[2]] * len(wdf[colIDs[0]])).values).copy(True)
            wdf = wdf.rename(index=str, columns={
                "tmpkey2": keys["names"][2]
            })
        dataframes.append(wdf)
    return pd.concat(dataframes)
