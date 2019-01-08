# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import os
from argparse import Namespace

# External Libraries
import pytest

# This Library
from rstoolbox.bin.minisilent import main as minisilent_main
from rstoolbox.bin.rename_decoys import main as rename_main
from rstoolbox.bin.check_mutants import main as check_mutants_main
from rstoolbox.bin.plot_fragments_rmsd import main as fragment_main
from rstoolbox.bin.regplot_rosetta import main as regplot_main
from rstoolbox.tests.helper import baseline_test_dir


class TestExecutables( object ):
    """
    Test usage for the stand alone executables.
    """
    def setup_method( self, method ):
        self.dirpath = os.path.join(os.path.dirname(__file__), '..', 'data')
        self.silent1 = os.path.join(self.dirpath, 'input_2seq.minisilent.gz')
        self.silent2 = os.path.join(self.dirpath, 'input_sse.minsilent.gz')
        self.silent3 = os.path.join(self.dirpath, 'input_ssebig.minisilent.gz')
        self.silent4 = os.path.join(self.dirpath, 'input_3ssepred.minisilent.gz')
        self.fastawt = os.path.join(self.dirpath, 'input_2seq.wt.seq')
        self.frag3   = os.path.join(self.dirpath, 'wauto.200.3mers.gz')
        self.frag3q  = os.path.join(self.dirpath, 'wauto.200.3mers.qual.gz')
        self.frag9   = os.path.join(self.dirpath, 'wauto.200.9mers.gz')
        self.frag9q  = os.path.join(self.dirpath, 'wauto.200.9mers.qual.gz')

    @pytest.fixture(autouse=True)
    def setup( self, tmpdir ):
        self.tmpdir = tmpdir.strpath

    def test_exe_minisilent_gz(self):
        options = Namespace(ifile=self.silent1, ifiles=None, force=False,
                            ofile=os.path.join(self.tmpdir, "minisilent.gz"))
        minisilent_main(options)

    def test_exe_minisilent(self):
        options = Namespace(ifile=self.silent1, ifiles=None, force=False,
                            ofile=os.path.join(self.tmpdir, "minisilent.sc"))
        minisilent_main(options)

    def test_exe_rename_gz(self):
        options = Namespace(ifile=self.silent1, prefix='test', force=False,
                            ofile=os.path.join(self.tmpdir, "renamed.gz"))
        rename_main(options)

    def test_exe_rename(self):
        options = Namespace(ifile=self.silent1, prefix='test', force=False,
                            ofile=os.path.join(self.tmpdir, "renamed.sc"))
        rename_main(options)

    @pytest.mark.mpl_image_compare(baseline_dir=baseline_test_dir(),
                                   filename='plot_exe_check_mutants_logo.png')
    def test_exe_check_mutants_logo(self):
        options = Namespace(ifile=self.silent1, ifiles=None, ifasta=None, seqID='B',
                            ffile=self.fastawt, ofile=os.path.join(self.tmpdir, 'mutants_'),
                            iformat='png', ifont=35)
        lfig, afig = check_mutants_main(options)
        return lfig

    @pytest.mark.mpl_image_compare(baseline_dir=baseline_test_dir(),
                                   filename='plot_exe_check_mutants_ali.png')
    def test_exe_check_mutants_ali(self):
        options = Namespace(ifile=self.silent1, ifiles=None, ifasta=None, seqID='B',
                            ffile=self.fastawt, ofile=os.path.join(self.tmpdir, 'mutants_'),
                            iformat='png', ifont=35)
        lfig, afig = check_mutants_main(options)
        return afig

    @pytest.mark.mpl_image_compare(baseline_dir=baseline_test_dir(),
                                   filename='plot_exe_fragments.png')
    def test_exe_plot_fragments(self):
        options = Namespace(fsmall=self.frag3, qsmall=self.frag3q, flarge=self.frag9,
                            qlarge=self.frag9q, pdb=None, silent=True, format='h', ofile=None)
        return fragment_main(options)

    @pytest.mark.mpl_image_compare(baseline_dir=baseline_test_dir(),
                                   filename='plot_exe_regplot.png')
    def test_exe_regplot(self):
        options = Namespace(ifile=self.silent3, ifiles=None, x='finalRMSD', y='score',
                            title='test plot', color=0, xlab='rmsd', ylab='score',
                            ylim=[-80, -20], xlim=[0, 6], fsize=(20, 20), silent=True,
                            ofile=None)
        return regplot_main(options)
