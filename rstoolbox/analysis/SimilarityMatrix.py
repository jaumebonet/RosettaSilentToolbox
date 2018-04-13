# @Author: Jaume Bonet <bonet>
# @Date:   11-May-2015
# @Email:  jaume.bonet@gmail.com
# @Filename: SimilarityMatrix.py
# @Last modified by:   bonet
# @Last modified time: 13-Apr-2018


from collections import deque
import os
import string


class SimilarityMatrix(object):
    '''
    Convert a biomatrix as downloaded from
    ftp://ftp.ncbi.nih.gov/blast/matrices
    to a python dictionary
    '''
    def __init__(self, data):
        self._data = data

    @staticmethod
    def get_matrix(matrixID):
        matdir = os.path.join(os.path.normpath(os.path.dirname(__file__)),
                              'matrices')
        matfile = os.path.join(matdir, matrixID.upper())
        if not os.path.isfile(matfile):
            raise ValueError("The provided SimilarityMatrix name does not exist in the database.")

        return SimilarityMatrix._parse_matrix(matfile)

    @staticmethod
    def _parse_matrix(matrix_file):
        data   = {}
        aalist = []
        fd = open(matrix_file)
        for line in fd:
            if line.startswith(' '):
                ldata = line.strip().split()
                for aa in ldata:
                    data.setdefault(aa, {})
                    aalist.append(aa)
            else:
                ldata  = deque(line.strip().split())
                aa = ldata.popleft()
                for i in range(len(aalist)):
                    data[aalist[i]][aa] = int(ldata[i])

        fd.close()
        return SimilarityMatrix( data )

    def get_value( self, k1, k2 ):
        if k1 in string.whitespace or k1 in string.punctuation:
            k1 = "*"
        if k2 in string.whitespace or k2 in string.punctuation:
            k2 = "*"
        return self._data[k1][k2]
