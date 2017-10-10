# @Author: Jaume Bonet <bonet>
# @Date:   29-Jun-2017
# @Email:  jaume.bonet@gmail.com
# @Filename: analyse.py
# @Last modified by:   bonet
# @Last modified time: 29-Sep-2017

import pandas as pd

def sequence_frequency_matrix( series ):
    sserie = series["sequence"].values
    table = {
        'C' : [], 'D' : [], 'S' : [], 'Q' : [], 'K' : [],
        'I' : [], 'P' : [], 'T' : [], 'F' : [], 'N' : [],
        'G' : [], 'H' : [], 'L' : [], 'R' : [], 'W' : [],
        'A' : [], 'V' : [], 'E' : [], 'Y' : [], 'M' : []
    }

    for x in range(len(sserie[0])):
        for k in table:
            table[k].append(float(0))
        for y in range(len(sserie)):
            aa = sserie[y][x]
            table[aa][-1] += float(1)
    for k in table:
        for x in range(len(table[k])):
            if table[k][x] != 0:
                table[k][x] /= float(len(sserie))

    keys = series["ranges"].iloc[0]
    if keys is not None:
        df = pd.DataFrame( table, index=keys )
    else:
        df = pd.DataFrame( table )

    df.index = df.index + 1
    return df
