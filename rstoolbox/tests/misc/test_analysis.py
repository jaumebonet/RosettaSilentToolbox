# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import os

# External Libraries
import pytest
import numpy as np
import matplotlib as mpl
if os.environ.get('DISPLAY', '') == '':
    mpl.use('Agg')
import matplotlib.pyplot as plt

# This Library
import rstoolbox.analysis as ra
from rstoolbox.tests.helper import baseline_test_dir


class TestAnalysis( object ):
    """
    Test utilities in analysis.
    """
    @pytest.mark.mpl_image_compare(baseline_dir=baseline_test_dir(),
                                   filename='plot_cumulative.png')
    def test_cumulative( self ):
        np.random.seed(0)
        data = np.random.rand(1000)
        fig = plt.figure(figsize=(25, 25))
        ax00 = plt.subplot2grid((2, 2), (0, 0), fig=fig)
        ax01 = plt.subplot2grid((2, 2), (0, 1), fig=fig)
        ax10 = plt.subplot2grid((2, 2), (1, 0), fig=fig)
        ax11 = plt.subplot2grid((2, 2), (1, 1), fig=fig)
        raw, y, x = ra.cumulative(data)
        ax00.plot(x, y)
        ax00.set_title('cumulative')
        raw, y, x = ra.cumulative(data, cumulative=0)
        ax01.plot(x, y)
        ax01.set_title('non-cumulative')
        raw, y, x = ra.cumulative(data, cumulative=-1)
        ax10.plot(x, y)
        ax10.set_title('reverse-cumulative')
        raw, y, x = ra.cumulative(data)
        ax11.plot(x, raw)
        ax11.set_title('raw data')
        plt.tight_layout()
        return fig
