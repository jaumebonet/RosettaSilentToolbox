# @Author: Jaume Bonet <bonet>
# @Date:   07-Apr-2018
# @Email:  jaume.bonet@gmail.com
# @Filename: test_read_sequence_files.py
# @Last modified by:   bonet
# @Last modified time: 07-Apr-2018


# Standard Libraries
import os

# External Libraries
import pytest

# This Library
from rstoolbox.io import read_fasta, write_fasta


class TestReadSilentFiles( object ):
    """
    Test reading silent files.
    Checks: apply different description and data retrival logic.
    """

    @pytest.fixture(autouse=True)
    def setup( self, tmpdir ):
        self.dirpath = os.path.join(os.path.dirname(__file__), '..', 'data')
        self.tmpdir = tmpdir.strpath

    def test_fasta(self):
        # Test simple read
        plain_id_string = "{}|PDBID|CHAIN|SEQUENCE"
        plain_ids = ["2TEP:A", "2TEP:B", "2TEP:C", "2TEP:D", "3TP2:A", "3TP2:B"]
        plain_ids = [plain_id_string.format(_) for _ in plain_ids]
        df1 = read_fasta(os.path.join(self.dirpath, "*.fa"), multi=True)

        assert sorted(plain_ids) == sorted(list(df1['description'].values))
        assert len(df1['sequence_A'].values[0]) == 236
        assert len(df1['sequence_A'].values[-1]) == 229
        assert df1.shape == (6, 2)
        assert len(df1['sequence_A'].unique()) == 2

        # Test expanded read
        expand_ids = ["2TEP", "3TP2"]
        df2 = read_fasta(os.path.join(self.dirpath, "*.fa"), expand=True, multi=True)
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
