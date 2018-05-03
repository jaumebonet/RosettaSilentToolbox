# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: random_frequency_matrix
"""
# Standard Libraries

# External Libraries
import pandas as pd
import numpy as np

# This Library

__all__ = ['random_frequency_matrix']


def random_frequency_matrix(size, seed=None):
    """Generate a random frequency matrix.

    :param int size: How many rows should the matrix have.
    :param in seed: If provided, set a fixed seed.
    """
    if seed is not None:
        np.random.seed(seed)
    alphabet = list("ARNDCQEGHILKMFPSTWYV")
    data = []
    for _ in range(size):
        data.append(np.random.dirichlet(np.ones(20), size=1)[0])
    return pd.DataFrame(data, columns=alphabet)
