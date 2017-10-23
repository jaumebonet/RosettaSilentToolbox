import pandas as pd

def add_column( df, name, value ):
    return df.assign(_placeholder=pd.Series([value]*df.shape[0])).rename(columns={"_placeholder": name})
