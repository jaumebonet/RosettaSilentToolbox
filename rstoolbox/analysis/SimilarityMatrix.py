from collections import deque
import os


class SimilarityMatrix(object):
    '''
    Convert a biomatrix as downloaded from
    ftp://ftp.ncbi.nih.gov/blast/matrices
    to a python dictionary
    '''

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
                l = line.strip().split()
                for aa in l:
                    data.setdefault(aa, {})
                    aalist.append(aa)
            else:
                l  = deque(line.strip().split())
                aa = l.popleft()
                for i in range(len(aalist)):
                    data[aalist[i]][aa] = int(l[i])

        fd.close()
        return data
