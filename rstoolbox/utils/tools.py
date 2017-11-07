import copy
import pandas as pd

def add_column( df, name, value ):
    """Adds a new column to the DataFrame with the given value.

    :param str name: Name of the new column
    :param value: Value that will be given to all rows of the new column (any type)
    :return: The new DataFrame with the added column
    :rtype: DataFrame
    """
    return df.assign(_placeholder=pd.Series([value]*df.shape[0])).rename(columns={"_placeholder": name})

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
