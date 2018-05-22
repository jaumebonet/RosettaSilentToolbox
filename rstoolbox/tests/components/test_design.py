# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import os
import copy

# External Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pytest

# This Library
import rstoolbox.io as ri
import rstoolbox.components as rc
import rstoolbox.plot as rp
import rstoolbox.analysis as ra
from rstoolbox.tests.helper import random_frequency_matrix


class TestDesign( object ):
    """
    Test usage of the DesignSeries/DesignFrame components.
    Checks: apply the attached functions of the objects. This includes
    setters/getters, reference data management and other utilities.
    In the reference data test, it includes a test for the SelectionContainer.
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

    def test_getters( self ):
        """
        Test usage of the getter functions.
        """

        # Assert types. Rows are DesignSeries, columns are not
        sc_des  = {"labels": ["MOTIF", "CONTACT", "CONTEXT"], "sequence": "AB"}
        df = ri.parse_rosetta_file(self.silent1, sc_des)
        assert isinstance(df, rc.DesignFrame)
        sr = df.iloc[0]
        assert isinstance(sr, rc.DesignSeries)
        assert not isinstance(df["description"], rc.DesignSeries)
        assert isinstance(df["description"], pd.Series)

        # Check working with sequence getters
        # We check everything both for DesignSeries and DesignFrame
        # DesignFrame returns Series, while DesignSeries returns the
        # actual data.
        assert sorted(df.get_available_sequences()) == ["A", "B"]
        assert sorted(sr.get_available_sequences()) == ["A", "B"]
        assert len(df.get_sequence("A")) == 6
        assert len(sr.get_sequence("A")) == 157
        assert df.get_sequence("B")[0] == sr.get_sequence("B")

        # Check working with label getters
        # We check everything both for DesignSeries and DesignFrame
        # DesignFrame returns Series, while DesignSeries returns the
        # actual data.
        assert sorted(df.get_available_labels()) == sorted(sc_des["labels"])
        assert sorted(sr.get_available_labels()) == sorted(sc_des["labels"])
        with pytest.raises(KeyError):
            sr.get_label("MOTIF")
        assert isinstance(df.get_label("MOTIF", "A")[0], rc.Selection)
        assert isinstance(sr.get_label("MOTIF", "A"), rc.Selection)
        assert str(df.get_label("CONTEXT", "B")[0]) == ""
        assert str(sr.get_label("CONTEXT", "A")) == "1-157"
        assert str(sr.get_label("CONTEXT", "B")) != str(sr.get_label("CONTEXT", "A"))

        # Check working with structure getters
        # We check everything both for DesignSeries and DesignFrame
        # DesignFrame returns Series, while DesignSeries returns the
        # actual data.
        sc_des  = {"sequence": "C", "structure": "C"}
        df = ri.parse_rosetta_file(self.silent2, sc_des)
        sr = df.iloc[0]
        assert df.get_available_structures() == ["C"]
        assert sr.get_available_structures() == ["C"]
        with pytest.raises(KeyError):
            assert len(df.get_structure("B")) == 6
        with pytest.raises(KeyError):
            assert len(sr.get_structure("B")) == 157
        assert df.get_structure("C")[0] == sr.get_structure("C")

        # Check working with structure prediction getters
        # We check everything both for DesignSeries and DesignFrame
        # DesignFrame returns Series, while DesignSeries returns the
        # actual data.
        assert df.get_available_structure_predictions() == []
        with pytest.raises(KeyError):
            assert len(df.get_structure_prediction("C")) == 6

        sc_des  = {'sequence': 'A', 'structure': 'A', 'psipred': 'A', 'dihedrals': 'A'}
        df = ri.parse_rosetta_file(self.silent4, sc_des)
        sr = df.iloc[0]
        assert df.get_available_structure_predictions() == ['A']
        assert df.get_structure_prediction('A')[0] == sr.get_structure_prediction('A')
        assert len(df.get_structure_prediction('A')[0]) == 88

        assert isinstance(df.get_dihedrals("A"), pd.DataFrame)
        assert isinstance(sr.get_dihedrals("A"), list)
        for e in sr.get_dihedrals("A"):
            assert isinstance(e, np.ndarray)
        assert np.array_equal(df.get_dihedrals("A").iloc[0][0], sr.get_dihedrals("A")[0])

        # these are the ranges of the rosetta angles.
        assert sr.get_phi("A").max() <= 180
        assert sr.get_phi("A").min() >= -180
        assert sr.get_psi("A").max() <= 180
        assert sr.get_psi("A").min() >= -180

    def test_reference( self ):
        """
        Test reference data usage.
        """
        # Without sequence/structure data, there are no references.
        df = ri.parse_rosetta_file(self.silent1)
        sr = df.iloc[0]
        refseq = "AYSTREILLALCIRDSRVHGNGTLHPVLELAARETPLRLSPEDTVVLRYHVLLEEIIERN" \
                 "SETFTETWNRFITHTEHVDLDFNSVFLEIFHRGDPSLGRALAWMAWCMHACRTLCCNQST" \
                 "PYYVVDLSVRGMLEASEGLDGWIHQQGGWSTLIEDNI"
        with pytest.raises(KeyError):
            df.add_reference_sequence("A", refseq)
        with pytest.raises(KeyError):
            sr.add_reference_sequence("A", refseq)
        with pytest.raises(KeyError):
            df.add_reference_shift("A", 2)

        # Get label and sequence/structure data to play with integer shift.
        sc_des  = {"labels": ["MOTIF", "CONTACT", "CONTEXT"], "sequence": "A"}
        _a = "9-26,28-29,31-32,35,37-40,67-68,70-71,89,91-116"
        _b = "AYSTREILLALCIRDSRVH"
        df = ri.parse_rosetta_file(self.silent1, sc_des)
        df.add_reference_sequence("A", refseq)
        sr = df.iloc[0]
        assert df.get_reference_sequence("A") == sr.get_reference_sequence("A")
        # Shift tests
        assert str(df.get_label("CONTACT", "A")[0]) == "1-19"
        ctcopy = copy.deepcopy(df.get_label("CONTACT", "A")[0])
        assert str(df.get_label("CONTACT", "B")[0]) == _a
        assert df.get_reference_sequence("A", df.get_label("CONTACT", "A")[0]) == _b
        df.add_reference_shift("A", 5, shift_labels=True)
        # Expected behaviour: all DesignSeries from a DesignFrame share reference data
        # and SelectionContainer
        assert df.get_reference_shift("A") == 5
        assert df.get_reference_shift("A") == sr.get_reference_shift("A")
        assert str(df.get_label("CONTACT", "A")[0]) == "5A-23A"
        assert str(df.get_label("CONTACT", "B")[0]) == _a
        assert str(df.get_label("CONTACT", "A")[0]) == str(sr.get_label("CONTACT", "A"))
        assert df.get_reference_sequence("A", df.get_label("CONTACT", "A")[0]) == _b
        assert str(ctcopy) == "1-19"
        assert df.get_reference_sequence("A", ctcopy) == _b
        df.delete_reference("A", shift_labels=True)
        assert str(df.get_label("CONTACT", "A")[0]) == "1-19"
        assert str(df.get_label("CONTACT", "B")[0]) == _a
        assert df.get_reference_shift("A") == 1
        with pytest.raises(KeyError):
            df.get_reference_sequence("A")

        # Let's work with an array-type shift
        ashift = list(range(1, len(refseq) + 1))
        ashift[30:] = list(np.array(ashift[30:]) + 5)
        with pytest.raises(ValueError):
            ashift.index(32)

        df = ri.parse_rosetta_file(self.silent1, sc_des)
        _c = "LHPVLELAARETPLRLSPEDTVVLRYHVLLEEI"
        df.add_reference_sequence("A", refseq)
        sr = df.iloc[1]
        assert str(sr.get_label("CONTACT", "A")) == "24-56"
        assert sr.get_reference_sequence("A", sr.get_label("CONTACT", "A")) == _c
        df.add_reference_shift("A", ashift, shift_labels=True)
        assert str(sr.get_label("CONTACT", "A")) == "24A-30A,36A-61A"
        assert sr.get_reference_sequence("A", sr.get_label("CONTACT", "A")) == _c
        df.delete_reference("A", shift_labels=True)
        assert str(sr.get_label("CONTACT", "A")) == "24-56"
        assert df.get_reference_shift("A") == 1
        with pytest.raises(KeyError):
            df.get_reference_sequence("A")

    def test_labels(self):
        sc_des  = {"scores": ["score"], "labels": ["MOTIF", "CONTACT", "CONTEXT"],
                   "sequence": "AB"}
        df = ri.parse_rosetta_file(self.silent1, sc_des)
        df = ra.selector_percentage(df, "A", "10-25", "test")
        df = ra.selector_percentage(df, "B", "12-20", "test")
        assert set(df.columns) == set(['score', 'lbl_MOTIF', 'lbl_CONTACT',
                                       'lbl_CONTEXT', 'sequence_A', 'sequence_B',
                                       'test_A_perc', 'test_B_perc'])
        assert len(df['test_A_perc'].unique()) == 1
        assert len(df['test_B_perc'].unique()) == 1
        assert df['test_A_perc'].values[0] == pytest.approx(0.1019, rel=1e-3)
        assert df['test_B_perc'].values[0] == pytest.approx(0.07758, rel=1e-3)

        df = ra.label_percentage(df, "A", "CONTEXT")
        df = ra.label_percentage(df, "A", "CONTACT")
        df = ra.label_percentage(df, "A", "MOTIF")
        df = ra.label_percentage(df, "B", "CONTACT")
        df = ra.label_percentage(df, "B", "MOTIF")
        df = ra.label_percentage(df, "B", "CONTEXT")
        assert len(df['CONTEXT_A_perc'].unique()) == 1
        assert df['CONTEXT_A_perc'].values[0] == 1
        assert len(df['MOTIF_A_perc'].unique()) == 1
        assert df['MOTIF_A_perc'].values[0] == 0
        assert len(df['CONTACT_A_perc'].unique()) > 1
        assert df['CONTACT_A_perc'].mean() == pytest.approx(0.0552, rel=1e-3)
        assert len(df['CONTEXT_B_perc'].unique()) == 1
        assert df['CONTEXT_B_perc'].values[0] == 0
        assert len(df['MOTIF_B_perc'].unique()) == 1
        assert df['MOTIF_B_perc'].values[0] == pytest.approx(0.1896, rel=1e-3)
        assert len(df['CONTACT_B_perc'].unique()) > 1
        assert df['CONTACT_B_perc'].mean() == pytest.approx(0.4669, rel=1e-3)

    def test_getseqs(self):
        sc_des  = {"sequence": "B"}

        # Start test
        df = ri.parse_rosetta_file(self.silent1, sc_des)
        assert df.shape[0] == 6
        df.get_sequence_with('B', [(1, 'T')]).shape[0] == 3

    @pytest.mark.mpl_image_compare(baseline_dir='../baseline_images',
                                   filename='plot_mutants_alignment.png')
    def test_mutants(self):
        # Static data
        refseq = "GSISDIRKDAEVRMDKAVEAFKNKLDKFKAAVRKVFPTEERIDMRPEIWIAQELRRIGDE" \
                 "FNAYRDANDKAAALGKDKEINWFDISQSLWDVQKLTDAAIKKIEAALADMEAWLTQ"
        columns = ["mutants_B", "mutant_count_B", "mutant_positions_B"]
        mut_number = [97, 91, 88, 90, 92, 92]
        mut_type = [
            "G1T,S2R,I3P,S4E,D5E,I6A,K8E,D9R,E11W,V12R,R13L,M14A,D15E,K16I,V18M,E19R,A20K,F21G,"
            "K22W,N23E,K24E,L25H,D26E,K27R,F28E,K29W,A30E,A31W,V32W,R33K,K34R,V35A,F36S,P37K,"
            "T38G,E39R,R41E,I42R,R45L,I48R,W49M,Q52A,E53A,R56A,D59E,E60I,Y64E,R65W,D66Q,A67M,N68R"
            ",D69L,K70E,A71M,A72E,A73K,L74E,G75R,D77N,K78P,E79N,I80A,N81G,W82E,F83E,D84K,I85M,S86K,"
            "Q87E,S88Q,L89K,W90K,D91E,V92A,Q93W,L95I,T96A,D97Y,A98Y,A99W,I100G,K101L,K102M,I103A,"
            "E104A,A105Y,A106W,L107I,A108K,D109Q,M110H,E111R,A112E,W113K,L114E,T115R,Q116K",
            "G1P,S2K,I3P,S4E,D5E,I6A,R7M,K8R,D9E,E11Y,V12K,R13L,M14I,D15K,A17Y,V18M,E19L,A20K,F21A,"
            "K22Q,N23K,K24E,L25A,D26Q,K27E,F28E,K29W,A30E,A31R,V32M,K34R,V35T,F36D,P37G,E39K,R41E,"
            "I42K,R45F,I48K,W49M,E53A,R56A,D59E,E60I,R65Y,D66W,N68F,D69L,A71L,A72Q,A73E,L74F,G75K,"
            "D77Y,K78P,E79S,I80V,N81R,F83E,D84E,I85Q,S86E,Q87E,S88A,L89R,W90K,D91R,V92L,Q93K,K94I,"
            "L95M,T96M,D97K,A98I,A99G,I100A,K101E,K102W,I103A,E104R,A105E,A106I,L107A,A108R,D109E,"
            "E111K,A112E,W113R,L114I,T115K,Q116R",
            "G1T,S2K,I3P,S4E,D5E,I6M,R7A,K8R,D9E,E11Y,V12K,D15L,V18L,E19K,A20Q,F21G,K22E,N23E,K24E,"
            "L25M,D26K,K27R,F28M,K29Y,A30E,A31Q,V32M,R33K,V35G,F36V,P37D,T38S,E39K,R41E,I42R,R45E,"
            "I48K,W49M,Q52I,E53A,R56A,D59E,E60L,Y64W,R65M,D66K,N68L,D69R,K70H,A71M,A72K,A73E,G75R,"
            "D77L,K78G,E79T,I80S,N81G,W82P,F83K,D84E,I85E,S86E,Q87K,S88H,L89W,W90R,D91W,V92I,Q93F,"
            "K94E,T96H,D97R,A98W,I100G,K101E,K102E,E104Q,A105R,L107A,A108E,D109I,M110Q,A112R,W113K,"
            "L114A,T115R,Q116W",
            "G1T,S2K,I3P,S4E,D5E,I6W,R7A,K8R,D9W,E11Y,V12K,R13E,M14H,D15L,A17M,V18A,A20K,F21H,K22R,"
            "N23K,K24E,L25M,D26E,K27I,F28E,K29W,A30E,A31E,V32L,R33K,K34R,V35R,F36D,P37G,T38K,R41E,"
            "I42K,R45W,I48R,W49M,Q52M,E53A,R56A,D59E,E60L,A63H,Y64H,R65M,D66Y,N68E,D69M,K70R,A72K,"
            "A73E,L74E,G75K,D77K,K78P,I80A,N81K,W82T,F83E,D84E,I85A,S86R,Q87R,S88A,L89R,W90R,D91E,"
            "V92I,Q93M,L95Y,T96H,D97H,A98E,I100G,K101R,K102L,A105E,L107M,A108R,D109R,M110L,E111M,"
            "A112E,W113R,L114H,T115K,Q116K",
            "G1K,S2K,I3W,S4E,D5E,I6M,R7M,K8R,D9E,V12R,R13Q,M14G,D15K,K16E,A17Y,V18A,E19Q,A20K,F21A,"
            "K22W,N23K,K24E,L25A,D26L,K27L,F28E,K29W,A30K,A31W,V32M,V35R,F36P,P37V,R41M,I42K,R45A,"
            "I48W,W49M,Q52A,E53A,R56A,D59E,E60H,A63I,R65W,D66Q,A67Q,N68K,D69L,K70E,A71H,A72E,A73K,"
            "G75R,D77I,K78P,E79N,I80V,N81P,W82E,F83E,D84E,I85L,S86E,Q87K,S88G,L89K,W90E,D91E,V92L,"
            "Q93K,K94R,L95I,T96E,D97E,A98E,I100A,K101R,K102M,I103A,A105K,A106Y,L107M,A108Q,D109E,"
            "M110L,E111R,A112K,W113K,L114M,T115E,Q116S",
            "G1P,S2R,I3P,S4E,D5E,I6M,R7A,K8R,D9F,E11K,V12E,R13E,D15H,A17H,V18E,A20K,F21A,K22Y,N23R"
            ",K24E,L25F,D26L,K27L,F28E,K29Y,A30E,A31L,V32A,R33I,K34R,V35K,F36N,R41P,I42K,R45Q,I48W"
            ",W49A,Q52A,E53A,R56A,D59E,E60I,A63Q,Y64W,R65M,D66Y,A67H,N68L,D69L,K70E,A71I,A72R,A73K"
            ",L74E,G75N,K76G,D77S,K78S,E79H,I80T,N81R,W82Y,F83E,D84E,I85R,S86E,Q87K,S88Y,L89R,W90K"
            ",D91L,V92A,Q93K,K94R,T96H,D97E,A98E,I100A,K102E,E104W,A105K,A106F,L107M,A108H,D109E,"
            "M110A,E111M,A112R,W113R,L114F,T115E,Q116S"
        ]
        mut_pos = [",".join([_[1:-1] for _ in m.split(",")]) for m in mut_type]
        sc_des  = {"labels": ["MOTIF", "CONTACT"], "sequence": "B"}

        # Start test
        df = ri.parse_rosetta_file(self.silent1, sc_des)
        df.add_reference_sequence("B", refseq)

        df = df.identify_mutants("B")
        for col in columns:
            assert col in df

        sr = df.iloc[0]
        assert df.get_reference_sequence("B") == sr.get_reference_sequence("B")
        assert df.get_identified_mutants() == ["B", ]

        for i, row in df.iterrows():
            # Check number of mutations
            assert row.get_mutation_count("B") == mut_number[i]
            # Check type of mutations
            assert row.get_mutations("B") == mut_type[i]
            # Check position of mutations
            assert row.get_mutation_positions("B") == mut_pos[i]

        # Make new variants
        dfm2 = df.iloc[0].generate_mutant_variants('B', [(1, "TGAP"), (14, "MAPT")])
        assert dfm2.shape[0] == 16
        assert 0 in dfm2.get_mutation_count('B')

        # Revert to WT
        dfwt = df.iloc[0:2].generate_wt_reversions('B', [1, 14])
        assert dfwt.shape[0] == 8

        dfwt = rc.DesignFrame({"description": ["reference"], "sequence_B": [refseq]})
        dfwt.add_reference_sequence('B', refseq)
        dfwt = dfwt.generate_mutant_variants('B', [(1, "TGP"), (6, "ERG"), (14, "MAT")])
        assert dfwt.shape[0] == 28
        dfwt = dfwt.generate_wt_reversions('B').identify_mutants('B')
        assert dfwt.shape[0] == 36
        assert 0 in dfwt.get_mutation_count('B').values
        assert refseq in dfwt.get_sequence('B').values

        # Make mutants from Matrix
        dfwt = rc.DesignFrame({"description": ["reference"], "sequence_B": [refseq]})
        dfwt.add_reference_sequence('B', refseq)
        matrix = random_frequency_matrix(len(df.get_reference_sequence('B')), 0)
        key_res = [3, 5, 8, 12, 15, 19, 25, 27]
        mutants = dfwt.generate_mutants_from_matrix('B', matrix, 5, key_res)
        assert isinstance(mutants, list)
        assert len(mutants) == 1
        mutants = mutants[0].identify_mutants('B')
        assert mutants.shape[0] == 5
        assert mutants.pssm_score_B.mean() != 0

        # write to resfiles
        df.make_resfile("B", "NATAA", os.path.join(self.tmpdir, "mutanttest.resfile"))
        for i, row in df.iterrows():
            newfile = os.path.join(self.tmpdir, "mutanttest" + "_{:>04d}".format(i) + ".resfile")
            assert row["resfile_B"] == newfile
            assert os.path.isfile(newfile)

        # write alignment
        ri.write_mutant_alignments(df, "B", os.path.join(self.tmpdir, "mutanttest.clw"))
        assert os.path.isfile(os.path.join(self.tmpdir, "mutanttest.clw"))

        # plot mutant
        fig = plt.figure(figsize=(30, 10))
        ax = plt.subplot2grid((1, 1), (0, 0), fig=fig)
        rp.plot_alignment(df, "B", ax, matrix="BLOSUM62")
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir='../baseline_images',
                                   filename='plot_summary.png')
    def test_summary_plot(self):
        # Start test
        df = ri.parse_rosetta_file(self.silent1)
        fig = plt.figure(figsize=(30, 30))

        rp.multiple_distributions(df, fig, (3, 3),
                                  ['score', 'GRMSD2Target', 'GRMSD2Template',
                                   'LRMSD2Target', 'LRMSDH2Target', 'LRMSDLH2Target',
                                   'design_score', 'packstat', 'rmsd_drift'])
        plt.tight_layout()
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir='../baseline_images',
                                   filename='plot_logo.png')
    def test_logo_plot(self):
        refseq = "GSISDIRKDAEVRMDKAVEAFKNKLDKFKAAVRKVFPTEERIDMRPEIWIAQELRRIGDE" \
                 "FNAYRDANDKAAALGKDKEINWFDISQSLWDVQKLTDAAIKKIEAALADMEAWLTQ"

        sc_des  = {"sequence": "B"}

        # Start test
        df = ri.parse_rosetta_file(self.silent1, sc_des)
        df.add_reference_sequence("B", refseq)

        fig, _ = rp.logo_plot( df, "B", refseq=True, line_break=50 )
        plt.tight_layout()
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir='../baseline_images',
                                   filename='plot_logo_noref.png')
    def test_logo_plot_noref(self):
        refseq = "GSISDIRKDAEVRMDKAVEAFKNKLDKFKAAVRKVFPTEERIDMRPEIWIAQELRRIGDE" \
                 "FNAYRDANDKAAALGKDKEINWFDISQSLWDVQKLTDAAIKKIEAALADMEAWLTQ"

        sc_des  = {"sequence": "B"}

        # Start test
        df = ri.parse_rosetta_file(self.silent1, sc_des)

        fig, _ = rp.logo_plot( df, "B", refseq=False, line_break=50 )
        plt.tight_layout()
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir='../baseline_images',
                                   filename='plot_logo_bits.png')
    def test_logo_plot_bits(self):
        refseq = "GSISDIRKDAEVRMDKAVEAFKNKLDKFKAAVRKVFPTEERIDMRPEIWIAQELRRIGDE" \
                 "FNAYRDANDKAAALGKDKEINWFDISQSLWDVQKLTDAAIKKIEAALADMEAWLTQ"

        sc_des  = {"sequence": "B"}

        # Start test
        df = ri.parse_rosetta_file(self.silent1, sc_des)
        df.add_reference_sequence("B", refseq)
        df = df.sequence_bits('B')

        fig, _ = rp.logo_plot( df, "B", refseq=True, line_break=50 )
        plt.tight_layout()
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir='../baseline_images',
                                   filename='plot_logo_bits_noref.png')
    def test_logo_plot_bits_noref(self):
        refseq = "GSISDIRKDAEVRMDKAVEAFKNKLDKFKAAVRKVFPTEERIDMRPEIWIAQELRRIGDE" \
                 "FNAYRDANDKAAALGKDKEINWFDISQSLWDVQKLTDAAIKKIEAALADMEAWLTQ"

        sc_des  = {"sequence": "B"}

        # Start test
        df = ri.parse_rosetta_file(self.silent1, sc_des)
        df = df.sequence_bits('B')

        fig, _ = rp.logo_plot( df, "B", refseq=False, line_break=50 )
        plt.tight_layout()
        return fig

    def test_sequence_similarities(self):
        refseq = "GSISDIRKDAEVRMDKAVEAFKNKLDKFKAAVRKVFPTEERIDMRPEIWIAQELRRIGDE" \
                 "FNAYRDANDKAAALGKDKEINWFDISQSLWDVQKLTDAAIKKIEAALADMEAWLTQ"
        diff1  = "....+.R+.A....+.A+.....+.++.....++.....E..DM.PE..IA..LR.IG+." \
                 "FNA......+.....K+.......+.+...+..K+...........+........+"
        diff2  = "000000100100000010000000000000000000000100110110011001101100" \
                 "11100000000000010000000000000000010000000000000000000000"
        diff3  = "000000100110110110100000000000001100111100110110011101101100" \
                 "11110010011001010010010000000000011000101011010001100000"

        sc_des  = {"scores": ["score"], "sequence": "B"}

        new_cols = ["blosum62_B_raw", "blosum62_B_perc", "blosum62_B_identity",
                    "blosum62_B_positive", "blosum62_B_negative", "blosum62_B_ali"]

        # Start test
        df = ri.parse_rosetta_file(self.silent1, sc_des)
        df.add_reference_sequence("B", refseq)

        # global sequence similarity
        dfss = ra.sequence_similarity( df, "B" )
        assert len(dfss.columns) == len(df.columns) + 6
        assert len(set(dfss.columns).difference(set(df.columns))) == len(new_cols)
        assert df.shape[0] == dfss.shape[0]
        assert dfss.blosum62_B_raw.mean() == 41.0
        assert dfss.blosum62_B_perc.mean() == pytest.approx(0.0692, rel=1e-3)
        assert dfss.blosum62_B_identity.mean() == pytest.approx(24.333, rel=1e-3)
        assert dfss.blosum62_B_positive.mean() == pytest.approx(46.166, rel=1e-3)
        assert dfss.blosum62_B_negative.mean() == pytest.approx(69.833, rel=1e-3)
        assert dfss.blosum62_B_ali.values[0] == diff1

        # local sequence similarity
        dfps = ra.positional_sequence_similarity(df, "B")
        assert dfps.shape == (len(refseq), 2)
        assert list(dfps.index.values) == list(range(1, len(refseq) + 1))
        assert dfps.identity_perc.mean() < dfps.positive_perc.mean()
        assert dfps.identity_perc.mean() == pytest.approx(0.2097, rel=1e-3)
        assert dfps.positive_perc.mean() == pytest.approx(0.3979, rel=1e-3)

        # binary similarity
        df01 = ra.binary_similarity(df, "B")
        assert len(df01.columns) == len(df.columns) + 1
        assert df01.identity_B_binary.values[0] == diff2

        # binary overlap
        assert "".join([str(_) for _ in ra.binary_overlap(df01, "B")]) == diff3

    def test_structure_similarities(self):
        sse_ref = "LEEEEEEELLLEEEEEEELLLLHHHHHHHHHHHHLLLLLLLLLLLEEEELLLEEEELL"
        diff1   = "LEEEEEEELLEEEEEEEELLLLHHHHHHHHHHHHLLLLLLLLLLEEEEELLLEEEEEL"
        sc_des  = {"scores": ["score"], "structure": "C"}

        # Start test
        df = ri.parse_rosetta_file(self.silent3, sc_des)
        df.add_reference_structure("C", sse_ref)

        # secondary structure distribution
        dfsse = ra.positional_structural_count(df, 'C')
        assert set(dfsse.columns.values) == set(['H', 'E', 'L'])
        assert dfsse.shape[0] == len(sse_ref)
        assert dfsse.H.mean() == pytest.approx(0.2033, rel=1e-3)
        assert dfsse.E.mean() == pytest.approx(0.4038, rel=1e-3)
        assert dfsse.L.mean() == pytest.approx(0.3927, rel=1e-3)

        # secondary structure match
        dfsm = ra.positional_structural_identity(df, 'C')
        assert set(dfsm.columns.values) == set(['identity_perc', 'sse', 'max_sse'])
        assert dfsm.shape[0] == len(sse_ref)
        assert "".join(list(dfsm.sse.values)) == sse_ref
        assert "".join(list(dfsm.max_sse.values)) == diff1
        assert dfsm.identity_perc.mean() == pytest.approx(0.8121, rel=1e-3)

        # percentages
        dfpc = ra.secondary_structure_percentage(df, 'C')
        assert 'structure_C_H' in dfpc.columns
        assert 'structure_C_E' in dfpc.columns
        assert 'structure_C_L' in dfpc.columns
        assert dfpc['structure_C_H'].max() == pytest.approx(0.2413, rel=1e-3)
        assert dfpc['structure_C_E'].mean() == pytest.approx(0.4038, rel=1e-3)
        assert dfpc['structure_C_L'].min() == pytest.approx(0.3275, rel=1e-3)

    @pytest.mark.mpl_image_compare(baseline_dir='../baseline_images',
                                   filename='plot_sse_profile.png')
    def test_sse_profile_plot(self):
        sse_ref = "LEEEEEEELLLEEEEEEELLLLHHHHHHHHHHHHLLLLLLLLLLLEEEELLLEEEELL"
        sc_des  = {"scores": ["score"], "structure": "C"}

        # Start test
        df = ri.parse_rosetta_file(self.silent3, sc_des)
        df.add_reference_structure("C", sse_ref)

        df1 = ra.positional_structural_count(df, 'C')
        df2 = ra.positional_structural_identity(df, 'C')

        fig = plt.figure(figsize=(35, 10))
        ax00 = plt.subplot2grid((1, 1), (0, 0))
        rp.positional_structural_similarity_plot(pd.concat([df1, df2], axis=1), ax00)
        plt.tight_layout()
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir='../baseline_images',
                                   filename='plot_ramachandran.png')
    def test_ramachandran_plot(self):
        # Start test
        sa_des = {"scores": ["score"], "sequence": "*", "dihedrals": "*"}
        df = ri.parse_rosetta_file(self.silent4, sa_des)

        fig = plt.figure(figsize=(15, 10))
        fig2 = plt.figure(figsize=(15, 10))
        with pytest.raises(ValueError):
            rp.plot_ramachandran(df, "A", fig2)
        rp.plot_ramachandran(df.iloc[0], "A", fig)
        plt.tight_layout()
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir='../baseline_images',
                                   filename='plot_dssp_vs_psipred.png')
    def test_plot_dssp_vs_psipred(self):
        # Start test
        sa_des = {"scores": ["score"], "psipred": "*", "structure": "*"}
        df = ri.parse_rosetta_file(self.silent4, sa_des)

        fig = plt.figure(figsize=(15, 10))
        ax = plt.gca()
        rp.plot_dssp_vs_psipred( df.iloc[0], "A", ax )
        plt.tight_layout()
        return fig
