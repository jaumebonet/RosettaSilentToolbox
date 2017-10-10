from SimilarityMatrix import SimilarityMatrix as SM
import pandas as pd
import numpy as np

def _calculate_linear_sequence_similarity( qseq, rseq, matrix ):
    score = 0;
    assert len(qseq) == len(rseq)
    for i in range(len(qseq)):
        score += matrix[qseq[i]][rseq[i]]
    return score

def _calculate_binary_sequence_similarity( qseq, rseq, matrix ):
    new_seq = ""
    assert len(qseq) == len(rseq)
    for i in range(len(qseq)):
        if  matrix[qseq[i]][rseq[i]] > 0:
            new_seq += "1"
        else:
            new_seq += "0"
    return new_seq

def linear_sequence_similarity( df, ref_seq, matrix="BLOSUM62", new_columns_prefix="similarity", seq_column="sequence" ):
    mat       = SM.get_matrix(matrix)
    all_seqs  = df[[seq_column]]
    sims      = []
    max_value = _calculate_linear_sequence_similarity( ref_seq, ref_seq, mat )
    for v in all_seqs.values:
        sims.append( _calculate_linear_sequence_similarity( v[0], ref_seq, mat ) )
    simsperc = np.array(sims) / float(max_value)
    df.insert( df.shape[1], new_columns_prefix + "_raw", pd.Series( sims, index=df.index ) )
    df.insert( df.shape[1], new_columns_prefix + "_perc", pd.Series( simsperc, index=df.index ) )
    return df

def binary_similarity( df, ref_seq, matrix="IDENTITY", seq_column="sequence" ):
    mat       = SM.get_matrix(matrix)
    all_seqs  = df[[seq_column]]
    sims      = []
    for v in all_seqs.values:
        sims.append( _calculate_binary_sequence_similarity( v[0], ref_seq, mat ) )
    df.insert( df.shape[1], matrix.lower() + "_binary", pd.Series( sims, index=df.index ) )
    return df

def sequence_frequency_matrix( series, seq_column="sequence" ):
    sserie = series[seq_column].values
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
