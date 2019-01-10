# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import os
import random

# External Libraries
import numpy as np
import pandas as pd
import matplotlib as mpl
if os.environ.get('DISPLAY', '') == '':
    mpl.use('Agg')
import matplotlib.pyplot as plt
import pytest

# This Library
import rstoolbox.utils as ru
import rstoolbox.plot as rp
import rstoolbox.io as ri
import rstoolbox.components as rc
from rstoolbox.tests.helper import baseline_test_dir


class TestPlotUtils( object ):
    """
    Test utilities in plots.
    """
    def setup_method( self, method ):
        self.dirpath = os.path.join(os.path.dirname(__file__), '..', 'data')

    @pytest.mark.mpl_image_compare(baseline_dir=baseline_test_dir(),
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

    @pytest.mark.mpl_image_compare(baseline_dir=baseline_test_dir(),
                                   filename='plot_96wells_blanc.png')
    def test_plot_96wells_blanc( self ):
        fig, ax = rp.plot_96wells()
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir=baseline_test_dir(),
                                   filename='plot_96wells_color.png')
    def test_plot_96wells_color( self ):
        np.random.seed(0)
        df = pd.DataFrame(np.random.randn(8, 12))
        fig, ax = rp.plot_96wells(cdata=df)
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir=baseline_test_dir(),
                                   filename='plot_96wells_size.png')
    def test_plot_96wells_size( self ):
        np.random.seed(0)
        df = pd.DataFrame(np.random.randn(8, 12))
        fig, ax = rp.plot_96wells(sdata=-df)
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir=baseline_test_dir(),
                                   filename='plot_96wells_bool.png')
    def test_plot_96wells_bool( self ):
        np.random.seed(0)
        df = pd.DataFrame(np.random.randn(8, 12))
        fig, ax = rp.plot_96wells(bdata=df < 0)
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir=baseline_test_dir(),
                                   filename='plot_96wells_all.png')
    def test_plot_96wells_all( self ):
        np.random.seed(0)
        df = pd.DataFrame(np.random.randn(8, 12))
        fig, ax = rp.plot_96wells(cdata=df, sdata=-df, bdata=df < 0)
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir=baseline_test_dir(),
                                   filename='plot_spr.png')
    def test_spr( self ):
        df = ri.read_SPR(os.path.join(self.dirpath, 'spr_data.csv.gz'))
        fig = plt.figure(figsize=(10, 6.7))
        ax = plt.subplot2grid((1, 1), (0, 0))
        rp.plot_SPR(df, ax, datacolor='black', fitcolor='red')
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir=baseline_test_dir(),
                                   filename='plot_cd.png')
    def test_cd( self ):
        df = pd.read_csv(os.path.join(self.dirpath, 'cd.csv'))
        fig = plt.figure(figsize=(10, 6.7))
        ax = plt.subplot2grid((1, 1), (0, 0))
        rp.plot_CD(df, ax)
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir=baseline_test_dir(),
                                   filename='plot_cd2.png')
    def test_cd_read( self ):
        def sampling( m, n ):
            return [i * n // m + n // (2 * m) for i in range(m)]

        df = ri.read_CD(os.path.join(self.dirpath, 'CD'), prefix='kx8', model='J-815')
        assert len(df['bin'].unique()) == 36
        assert sampling(5, 35) == [3, 10, 17, 24, 31]
        fig = plt.figure(figsize=(10, 6.7))
        ax = plt.subplot2grid((1, 1), (0, 0))
        rp.plot_CD(df, ax, sample=5)
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir=baseline_test_dir(),
                                   filename='plot_cd_chirascan.png')
    def test_cd_read_chirascan( self ):
        df = ri.read_CD(os.path.join(self.dirpath, 'chirascan_cd.csv'),  model='chirascan')
        fig  = plt.figure(figsize=(15, 15))
        grid = (3, 2)
        for i, sample in enumerate(sorted(df.keys())):
            ax = plt.subplot2grid(grid, (int(i / 2), i % 2), fig=fig)
            rp.plot_CD(df[sample], ax, sample=5)
            ru.add_top_title(ax, sample)
        plt.tight_layout()
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir=baseline_test_dir(),
                                   filename='plot_mals.png')
    def test_mals( self ):
        df = pd.read_csv(os.path.join(self.dirpath, 'mals.csv'))
        fig = plt.figure(figsize=(10, 6.7))
        ax = plt.subplot2grid((1, 1), (0, 0))
        rp.plot_MALS(df, ax)
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir=baseline_test_dir(),
                                   filename='plot_mals2.png')
    def test_mals_read( self ):
        df = ri.read_MALS(filename=os.path.join(self.dirpath, 'mota_1kx8_d2.csv'),
                          mmfile=os.path.join(self.dirpath, 'mota_1kx8_d2_mm.csv'))
        fig = plt.figure(figsize=(10, 6.7))
        ax = plt.subplot2grid((1, 1), (0, 0))
        rp.plot_MALS(df, ax)
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir=baseline_test_dir(),
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

    @pytest.mark.mpl_image_compare(baseline_dir=baseline_test_dir(),
                                   filename='plot_color_hydrophobicity.png')
    def test_color_scheme_hydrophobicity( self ):
        df = rc.DesignFrame(pd.read_csv(os.path.join(self.dirpath, 'logo_plot_sequence.csv'),
                                        header=None).rename(columns={0: 'sequence_A'}))
        fig, axs = rp.logo_plot(df, "A", refseq=False, font_size=10, hight_prop=2,
                                colors='HYDROPHOBICITY')
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir=baseline_test_dir(),
                                   filename='plot_color_chemistry.png')
    def test_color_scheme_chemistry( self ):
        df = rc.DesignFrame(pd.read_csv(os.path.join(self.dirpath, 'logo_plot_sequence.csv'),
                                        header=None).rename(columns={0: 'sequence_A'}))
        fig, axs = rp.logo_plot(df, "A", refseq=False, line_break=50, font_size=10, hight_prop=2,
                                colors="CHEMISTRY")
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir=baseline_test_dir(),
                                   filename='plot_color_charge.png')
    def test_color_scheme_charge( self ):
        df = rc.DesignFrame(pd.read_csv(os.path.join(self.dirpath, 'logo_plot_sequence.csv'),
                                        header=None).rename(columns={0: 'sequence_A'}))
        fig, axs = rp.logo_plot(df, "A", refseq=False, line_break=50, font_size=10, hight_prop=2,
                                colors="CHARGE")
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir=baseline_test_dir(),
                                   filename='plot_color_custom.png')
    def test_color_scheme_custom( self ):
        custom = {
            'A': '#e6194b', 'C': '#3cb44b', 'D': '#ffe119', 'E': '#ffe119',
            'F': '#f58231', 'G': '#911eb4', 'H': '#46f0f0', 'I': '#f032e6',
            'K': '#d2f53c', 'L': '#d2f53c', 'M': '#008080', 'N': '#e6beff',
            'P': '#aa6e28', 'Q': '#fffac8', 'R': '#800000', 'S': '#aaffc3',
            'T': '#808000', 'V': '#ffd8b1', 'W': '#000080', 'Y': '#808080'
        }
        df = rc.DesignFrame(pd.read_csv(os.path.join(self.dirpath, 'logo_plot_sequence.csv'),
                                        header=None).rename(columns={0: 'sequence_A'}))
        fig, axs = rp.logo_plot(df, "A", refseq=False, line_break=50, font_size=10, hight_prop=2,
                                colors=custom)
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir=baseline_test_dir(),
                                   filename='plot_logo_sse.png')
    def test_sse_logo(self):
        custom = {
            'E': '#0000FF', 'H': '#00FF00', 'L': '#FF0000'
        }
        ff = os.path.join(self.dirpath, 'input_3ssepred.minisilent.gz')
        df = ri.parse_rosetta_file(ff, {'structure': 'A'})
        fs = df.structure_bits('A')
        fig, axs = rp.logo_plot(fs, "A", refseq=False, line_break=50, font_size=10, hight_prop=2,
                                colors=custom)
        return fig

    def test_plot_labels( self ):
        plt.plot([random.randint(0, 100) for i in range(100)], label='text1')
        plt.plot([random.randint(0, 100) for i in range(100)], label='text2')
        ax = plt.gca()
        ax.legend()
        inilabs = [x.get_text() for x in ax.get_legend().texts]
        newlabs = ['text3', 'text4']
        ru.edit_legend_text(ax, newlabs, 'lines')
        endlabs = [x.get_text() for x in ax.get_legend().texts]
        with pytest.raises(IndexError):
            ru.edit_legend_text(ax, ['text1', 'text2', 'text3'], 'lines')
        plt.close()

        assert newlabs == endlabs
        assert endlabs != inilabs

    def test_colors( self ):
        red = [255, 0, 0]
        newred = ru.color_variant(red, brightness_offset=1)
        assert newred == '#ff0101'

        cmap = ru.add_white_to_cmap(color='blue')
        assert cmap.name == 'FromWhite'
        assert cmap.N == 256
