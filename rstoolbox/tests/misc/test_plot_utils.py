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
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pytest

# This Library
import rstoolbox.utils as ru
import rstoolbox.plot as rp
import rstoolbox.io as ri


class TestPlotUtils( object ):
    """
    Test utilities in plots.
    """
    def setup_method( self, method ):
        self.dirpath = os.path.join(os.path.dirname(__file__), '..', 'data')

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

    @pytest.mark.mpl_image_compare(baseline_dir='../baseline_images',
                                   filename='plot_spr.png')
    def test_spr( self ):
        df = ri.read_SPR(os.path.join(self.dirpath, 'spr_data.csv.gz'))
        fig = plt.figure(figsize=(10, 6.7))
        ax = plt.subplot2grid((1, 1), (0, 0))
        rp.plot_SPR(df, ax, datacolor='black', fitcolor='red')
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir='../baseline_images',
                                   filename='plot_cd.png')
    def test_cd( self ):
        df = pd.read_csv(os.path.join(self.dirpath, 'cd.csv'))
        fig = plt.figure(figsize=(10, 6.7))
        ax = plt.subplot2grid((1, 1), (0, 0))
        rp.plot_CD(df, ax)
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir='../baseline_images',
                                   filename='plot_mals.png')
    def test_mals( self ):
        df = pd.read_csv(os.path.join(self.dirpath, 'mals.csv'))
        fig = plt.figure(figsize=(10, 6.7))
        ax = plt.subplot2grid((1, 1), (0, 0))
        rp.plot_MALS(df, ax)
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir='../baseline_images',
                                   filename='plot_tm.png')
    def test_thermal_melt( self ):
        df = pd.read_csv(os.path.join(self.dirpath, 'thermal_melt.csv'))
        fig = plt.figure(figsize=(10, 6.7))
        ax = plt.subplot2grid((1, 1), (0, 0))
        rp.plot_thermal_melt(df, ax)
        return fig

    def test_multi_fastq( self ):
        indat = {'binder1': {'conc1': os.path.join(self.dirpath, 'cdk2_rand_001.fasq.gz'),
                             'conc2': os.path.join(self.dirpath, 'cdk2_rand_002.fasq.gz'),
                             'conc3': os.path.join(self.dirpath, 'cdk2_rand_003.fasq.gz')},
                 'binder2': {'conc1': os.path.join(self.dirpath, 'cdk2_rand_004.fasq.gz'),
                             'conc2': os.path.join(self.dirpath, 'cdk2_rand_005.fasq.gz'),
                             'conc3': os.path.join(self.dirpath, 'cdk2_rand_006.fasq.gz')}}
        enrich = {'binder1': ['conc1', 'conc3'],
                  'binder2': ['conc1', 'conc3']}
        bounds = ['GAS', 'PGT']
        matches = ['ALKKI']
        df = ru.sequencing_enrichment(indat, enrich, bounds, matches)
        assert 'binder2_conc1' in df.columns
        assert 'binder1_conc3' in df.columns
        assert 'enrichment_binder1' in df.columns
        assert df.shape == (20, 11)
        assert df['enrichment_binder2'].mean() == pytest.approx(1.13, rel=1e-3)
        assert df['enrichment_binder1'].max() == pytest.approx(5, rel=1e-3)
