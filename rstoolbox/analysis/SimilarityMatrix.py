# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. class:: SimilarityMatrix
"""
# Standard Libraries
from collections import deque
import os
import string

# External Libraries

# This Library


__all__ = ['SimilarityMatrix']


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
                for i, i_ in enumerate(aalist):
                    data[i_][aa] = int(ldata[i])

        fd.close()
        return SimilarityMatrix( data )

    def get_value( self, k1, k2 ):
        if k1 in string.whitespace or k1 in string.punctuation:
            k1 = "*"
        if k2 in string.whitespace or k2 in string.punctuation:
            k2 = "*"
        return self._data[k1][k2]
