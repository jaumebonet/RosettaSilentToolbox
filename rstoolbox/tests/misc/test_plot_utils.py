# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries

# External Libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pytest

# This Library
import rstoolbox.utils as ru
import rstoolbox.plot as rp


class TestPlotUtils( object ):
    """
    Test utilities in plots.
    """
    @pytest.mark.mpl_image_compare(baseline_dir='../baseline_images',
                                   filename='plot_titles.png')
    def test_plot_titles( self ):
        fig  = plt.figure(figsize=(10, 10))
        grid = (1, 1)
        ax00 = plt.subplot2grid(grid, (0, 0), fig=fig)
        ax00.plot([1, 2, 3], [1, 2, 3])
        ru.add_right_title(ax00, 'right title text', rotation=-90)
        ru.add_top_title(ax00, 'top title text')
        ru.add_left_title(ax00, 'left title text', rotation=90)
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir='../baseline_images',
                                   filename='plot_96wells_blanc.png')
    def test_plot_96wells_blanc( self ):
        fig, ax = rp.plot_96wells()
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir='../baseline_images',
                                   filename='plot_96wells_color.png')
    def test_plot_96wells_color( self ):
        np.random.seed(0)
        df = pd.DataFrame(np.random.randn(8, 12))
        fig, ax = rp.plot_96wells(cdata=df)
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir='../baseline_images',
                                   filename='plot_96wells_size.png')
    def test_plot_96wells_size( self ):
        np.random.seed(0)
        df = pd.DataFrame(np.random.randn(8, 12))
        fig, ax = rp.plot_96wells(sdata=-df)
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir='../baseline_images',
                                   filename='plot_96wells_bool.png')
    def test_plot_96wells_bool( self ):
        np.random.seed(0)
        df = pd.DataFrame(np.random.randn(8, 12))
        fig, ax = rp.plot_96wells(bdata=df < 0)
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir='../baseline_images',
                                   filename='plot_96wells_all.png')
    def test_plot_96wells_all( self ):
        np.random.seed(0)
        df = pd.DataFrame(np.random.randn(8, 12))
        fig, ax = rp.plot_96wells(cdata=df, sdata=-df, bdata=df < 0)
        return fig
