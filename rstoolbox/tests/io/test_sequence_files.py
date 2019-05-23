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
import matplotlib as mpl
if os.environ.get('DISPLAY', '') == '':
    mpl.use('Agg')
import matplotlib.pyplot as plt

# This Library
from rstoolbox.io import read_fasta, write_fasta, read_hmmsearch
from rstoolbox.io import parse_rosetta_file, pymol_mutant_selector
from rstoolbox.io import parse_master_file
from rstoolbox.plot import sequence_frequency_plot
from rstoolbox.tests.helper import baseline_test_dir


class TestReadSilentFiles( object ):
    """
    Test reading silent files.
    Checks: apply different description and data retrival logic.
    """

    @pytest.fixture(autouse=True)
    def setup( self, tmpdir ):
        self.dirpath = os.path.join(os.path.dirname(__file__), '..', 'data')
        self.tmpdir = tmpdir.strpath

    @pytest.mark.mpl_image_compare(baseline_dir=baseline_test_dir(),
                                   filename='seq_freq_plot_fasta.png')
    def test_fasta(self):
        # Test simple read
        plain_id_string = "{}|PDBID|CHAIN|SEQUENCE"
        plain_ids = ["2TEP:A", "2TEP:B", "2TEP:C", "2TEP:D", "3TP2:A", "3TP2:B"]
        plain_ids = [plain_id_string.format(_) for _ in plain_ids]
        df1 = read_fasta(os.path.join(self.dirpath, "*.fa$"), multi=True)

        assert sorted(plain_ids) == sorted(list(df1['description'].values))
        assert len(df1['sequence_A'].values[0]) == 236
        assert len(df1['sequence_A'].values[-1]) == 229
        assert df1.shape == (6, 2)
        assert len(df1['sequence_A'].unique()) == 2

        # Test expanded read
        expand_ids = ["2TEP", "3TP2"]
        df2 = read_fasta(os.path.join(self.dirpath, "*.fa$"), expand=True, multi=True)
        assert sorted(expand_ids) == sorted(list(df2['description'].values))
        assert df2.shape == (2, 5)
        assert 'sequence_A' in df2
        assert 'sequence_B' in df2
        assert 'sequence_C' in df2
        assert 'sequence_D' in df2

        # Test write
        all_text = [
            ">2TEP:A|PDBID|CHAIN|SEQUENCE:A",
            "AETVSFNFNSFSEGNPAINFQGDVTVLSNGNIQLTNLNKVNSVGRVLYAMPVRIWSSATGNVASFLTSFSFEMKDIKD"
            "YDPADGIIFFIAPEDTQIPAGSIGGGTLGVSDTKGAGHFVGVEFDTYSNSEYNDPPTDHVGIDVNSVDSVKTVPWNSV"
            "SGAVVKVTVIYDSSTKTLSVAVTNDNGDITTIAQVVDLKAKLPERVKFGFSASGSLGGRQIHLIRSWSFTSTLITTTRRS",
            ">2TEP:B|PDBID|CHAIN|SEQUENCE:A",
            "AETVSFNFNSFSEGNPAINFQGDVTVLSNGNIQLTNLNKVNSVGRVLYAMPVRIWSSATGNVASFLTSFSFEMKDIKDY"
            "DPADGIIFFIAPEDTQIPAGSIGGGTLGVSDTKGAGHFVGVEFDTYSNSEYNDPPTDHVGIDVNSVDSVKTVPWNSVSG"
            "AVVKVTVIYDSSTKTLSVAVTNDNGDITTIAQVVDLKAKLPERVKFGFSASGSLGGRQIHLIRSWSFTSTLITTTRRS",
            ">2TEP:C|PDBID|CHAIN|SEQUENCE:A",
            "AETVSFNFNSFSEGNPAINFQGDVTVLSNGNIQLTNLNKVNSVGRVLYAMPVRIWSSATGNVASFLTSFSFEMKDIKDY"
            "DPADGIIFFIAPEDTQIPAGSIGGGTLGVSDTKGAGHFVGVEFDTYSNSEYNDPPTDHVGIDVNSVDSVKTVPWNSVSG"
            "AVVKVTVIYDSSTKTLSVAVTNDNGDITTIAQVVDLKAKLPERVKFGFSASGSLGGRQIHLIRSWSFTSTLITTTRRS",
            ">2TEP:D|PDBID|CHAIN|SEQUENCE:A",
            "AETVSFNFNSFSEGNPAINFQGDVTVLSNGNIQLTNLNKVNSVGRVLYAMPVRIWSSATGNVASFLTSFSFEMKDIKDY"
            "DPADGIIFFIAPEDTQIPAGSIGGGTLGVSDTKGAGHFVGVEFDTYSNSEYNDPPTDHVGIDVNSVDSVKTVPWNSVSG"
            "AVVKVTVIYDSSTKTLSVAVTNDNGDITTIAQVVDLKAKLPERVKFGFSASGSLGGRQIHLIRSWSFTSTLITTTRRS",
            ">3TP2:A|PDBID|CHAIN|SEQUENCE:A",
            "GAMTSWRDKSAKVQVKESELPSSIPAQTGLTFNIWYNKWSQGFAGNTRFVSPFALQPQLHSGKTRGDNDGQLFFCLFFA"
            "KGMCCLGPKCEYLHHIPDEEDIGKLALRTEVLDCFGREKFADYREDMGGIGSFRKKNKTLYVGGIDGALNSKHLKPAQI"
            "ESRIRFVFSRLGDIDRIRYVESKNCGFVKFKYQANAEFAKEAMSNQTLLLPSDKEWDDRREGTGLLVKWAN",
            ">3TP2:B|PDBID|CHAIN|SEQUENCE:A",
            "GAMTSWRDKSAKVQVKESELPSSIPAQTGLTFNIWYNKWSQGFAGNTRFVSPFALQPQLHSGKTRGDNDGQLFFCLFFA"
            "KGMCCLGPKCEYLHHIPDEEDIGKLALRTEVLDCFGREKFADYREDMGGIGSFRKKNKTLYVGGIDGALNSKHLKPAQI"
            "ESRIRFVFSRLGDIDRIRYVESKNCGFVKFKYQANAEFAKEAMSNQTLLLPSDKEWDDRREGTGLLVKWAN"
        ]
        assert write_fasta(df1, "A") == "\n".join(all_text) + "\n"

        # Pick chains from expanded read
        picked_text = [
            ">2TEP:A",
            "AETVSFNFNSFSEGNPAINFQGDVTVLSNGNIQLTNLNKVNSVGRVLYAMPVRIWSSATGNVASFLTSFSFEMKDIKDYDPA"
            "DGIIFFIAPEDTQIPAGSIGGGTLGVSDTKGAGHFVGVEFDTYSNSEYNDPPTDHVGIDVNSVDSVKTVPWNSVSGAVVKVT"
            "VIYDSSTKTLSVAVTNDNGDITTIAQVVDLKAKLPERVKFGFSASGSLGGRQIHLIRSWSFTSTLITTTRRS",
            ">3TP2:A",
            "GAMTSWRDKSAKVQVKESELPSSIPAQTGLTFNIWYNKWSQGFAGNTRFVSPFALQPQLHSGKTRGDNDGQLFFCLFFAKGM"
            "CCLGPKCEYLHHIPDEEDIGKLALRTEVLDCFGREKFADYREDMGGIGSFRKKNKTLYVGGIDGALNSKHLKPAQIESRIRF"
            "VFSRLGDIDRIRYVESKNCGFVKFKYQANAEFAKEAMSNQTLLLPSDKEWDDRREGTGLLVKWAN",
            ">2TEP:C",
            "AETVSFNFNSFSEGNPAINFQGDVTVLSNGNIQLTNLNKVNSVGRVLYAMPVRIWSSATGNVASFLTSFSFEMKDIKDYDPA"
            "DGIIFFIAPEDTQIPAGSIGGGTLGVSDTKGAGHFVGVEFDTYSNSEYNDPPTDHVGIDVNSVDSVKTVPWNSVSGAVVKVT"
            "VIYDSSTKTLSVAVTNDNGDITTIAQVVDLKAKLPERVKFGFSASGSLGGRQIHLIRSWSFTSTLITTTRRS"
        ]
        assert write_fasta(df2, "AC") == "\n".join(picked_text) + "\n"

        # Individual prints
        write_fasta(df2, "A", filename=os.path.join(self.tmpdir, "singles.fa"), split=True)

        for _ in ["singles_f0001.fa", "singles_f0002.fa"]:
            df = read_fasta(os.path.join(self.tmpdir, _))
            assert df.shape == (1, 2)
            newseq = df["sequence_A"].values[0]
            newid = df["description"].values[0]
            expected = df2[df2["description"] == newid.split(":")[0]]["sequence_A"].values[0]
            assert expected == newseq

        # plot
        fig = plt.figure(figsize=(25, 10))
        ax = plt.subplot2grid((1, 1), (0, 0))
        sequence_frequency_plot(df1, 'A', ax, refseq=False, key_residues='12-35', clean_unused=0)
        return fig

    def test_hmm( self ):
        df = read_hmmsearch(os.path.join(self.dirpath, 'search.hmm.gz'))
        assert df.shape[0] == 4932
        assert len(df['description'].unique()) == 4927
        assert df[df['full-e-value'] < 10].shape[0] == 2650

        df = read_hmmsearch(os.path.join(self.dirpath, 'search2.hmm.gz'))
        assert df.shape[0] == 11
        assert len(df.iloc[0]['sequence']) == 87

        df = read_hmmsearch(os.path.join(self.dirpath, 'scan.hmm.gz'))
        assert df.shape[0] == 9
        assert len(df.iloc[0]['sequence']) == 99

    def test_pymol( self ):
        df = parse_rosetta_file(os.path.join(self.dirpath, 'input_2seq.minisilent.gz'),
                                {'sequence': 'B'})
        df.add_reference_sequence('B', df.iloc[0].get_sequence('B'))
        df = df.identify_mutants('B').head()

        pick1 = ""
        pick2 = "sele test_3lhp_binder_labeled_00002_mut, test_3lhp_binder_labeled_00002 and " \
                "((c. B and (i. 1-2 or i. 7-9 or i. 11-12 or i. 14-17 or i. 19 or i. 21-23 or " \
                "i. 25-27 or i. 31-33 or i. 35-39 or i. 42 or i. 45 or i. 48 or i. 52 or " \
                "i. 64-68 or i. 70-75 or i. 77 or i. 79-82 or i. 84-86 or i. 88-89 or " \
                "i. 91-102 or i. 104-111 or i. 113-116)))"
        sel = pymol_mutant_selector(df)
        assert len(sel[0]) == 0
        assert sel[0] == pick1
        assert len(sel[1]) != 0
        assert sel[1] == pick2

    def test_master( self ):
        df = parse_master_file(os.path.join(self.dirpath, 'master.search'),
                               max_rmsd=1.4, piece_count=2, shift_0=True)
        assert df.rmsd.max() == 1.3967
        assert df.shape == (42, 5)
        assert df.iloc[-1].match == [[34, 40], [42, 48]]
