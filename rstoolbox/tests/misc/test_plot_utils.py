# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries

# External Libraries
import matplotlib.pyplot as plt
import pytest

# This Library
import rstoolbox.utils as ru


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
