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
import matplotlib.pyplot as plt
import pytest

# This Library
from rstoolbox.bin.minisilent import main as minisilent_main


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

    @pytest.fixture(autouse=True)
    def setup( self, tmpdir ):
        self.tmpdir = tmpdir.strpath

    def test_minisilent_gz(self):
        options = Namespace(ifile=self.silent1, ifiles=None, force=False,
                            ofile=os.path.join(self.tmpdir, "minisilent.gz"))
        minisilent_main(options)

    def test_minisilent(self):
        options = Namespace(ifile=self.silent1, ifiles=None, force=False,
                            ofile=os.path.join(self.tmpdir, "minisilent.sc"))
        minisilent_main(options)
