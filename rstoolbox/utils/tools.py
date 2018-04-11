# @Author: Jaume Bonet <bonet>
# @Date:   19-Feb-2018
# @Email:  jaume.bonet@gmail.com
# @Filename: tools.py
# @Last modified by:   bonet
# @Last modified time: 11-Apr-2018


import copy
import textwrap

import pandas as pd


def format_Ipython():
    """Ensure ``monospace`` representation of :class:`~pandas.DataFrame`
    in **Jupyter Notebooks**.

    Just need to call it after importing the library.
    """
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
    :param name: Name of the new column
    :type name: :class:`str`
    :param value: Value that will be given to all rows of the new column (any type)

    :return: :class:`~pandas.DataFrame` - The data container with the new column
    """
    data = pd.Series([value] * df.shape[0])
    return df.assign(_placeholder=data).rename(columns={"_placeholder": name})


def split_values( df, keys ):
    dataframes = []
    for k in keys["split"]:
        colIDs = copy.copy(keys["keep"])
        colIDs.append(k[0])
        wdf = df[colIDs]
        wdf = wdf.assign(temporarykey1=pd.Series([k[1]]*len(wdf[colIDs[0]])).values).copy(True)
        wdf = wdf.rename(index=str, columns={
            k[0]: keys["names"][0],
            "temporarykey1": keys["names"][1]
        })
        if ( len(k) > 2 ):
            wdf = wdf.assign(temporarykey2=pd.Series([k[2]]*len(wdf[colIDs[0]])).values).copy(True)
            wdf = wdf.rename(index=str, columns={
                "temporarykey2": keys["names"][2]
            })
        dataframes.append(wdf)
    return pd.concat(dataframes)
