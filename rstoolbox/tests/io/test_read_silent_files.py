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
import six
import pytest
import pandas as pd

# This Library
import rstoolbox.io as ri
import rstoolbox.components as rc


if six.PY3:
    def cmp(a, b):
        if a == b:
            return 0
        else:
            return 1


class TestReadSilentFiles( object ):
    """
    Test reading silent files.
    Checks: apply different description and data retrival logic.
    """

    def setup_method( self, method ):
        self.dirpath = os.path.join(os.path.dirname(__file__), '..', 'data')
        self.silent1 = os.path.join(self.dirpath, 'input_2seq.minisilent.gz')
        self.silent2 = os.path.join(self.dirpath, 'input_sse.minsilent.gz')
        self.silent3 = os.path.join(self.dirpath, 'input_symmetry.minisilent.gz')
        self.silent4 = os.path.join(self.dirpath, 'motifgraft.scores.gz')
        self.contact = os.path.join(self.dirpath, 'contacts.csv.gz')
        self.jsonscr = os.path.join(self.dirpath, 'score.json.gz')
        self.defhead = ["score", "fa_atr", "fa_rep", "fa_sol", "fa_intra_rep", "fa_elec",
                        "pro_close", "hbond_sr_bb", "hbond_lr_bb", "hbond_bb_sc", "hbond_sc",
                        "dslf_fa13", "rama", "omega", "fa_dun", "p_aa_pp", "yhh_planarity",
                        "ref", "BUNS", "B_ni_mtcontacts", "B_ni_rmsd", "B_ni_rmsd_threshold",
                        "B_ni_trials", "GRMSD2Target", "GRMSD2Template", "LRMSD2Target",
                        "LRMSDH2Target", "LRMSDLH2Target", "cav_vol", "design_score", "packstat",
                        "rmsd_drift", "time", "description"]
        self.symhead = ["score", "fa_atr", "fa_rep", "fa_sol", "fa_intra_rep", "fa_elec",
                        "pro_close", "hbond_sr_bb", "hbond_lr_bb", "hbond_bb_sc", "hbond_sc",
                        "dslf_fa13", "rama", "omega", "fa_dun", "p_aa_pp", "yhh_planarity", "ref",
                        "bsa", "ddg", "ddg_filter", "dsasa", "sasa", "sasa_h", "sasa_p", "shape",
                        "complex_normalized", "dG_cross", "dG_cross/dSASAx100", "dG_separated",
                        "dG_separated/dSASAx100", "dSASA_hphobic", "dSASA_int", "dSASA_polar",
                        "delta_unsatHbonds", "hbond_E_fraction", "hbonds_int", "nres_all",
                        "nres_int", "packstat", "per_residue_energy_int", "sc_value",
                        "shape_int_area", "side1_normalized", "side1_score", "side2_normalized",
                        "side2_score", "time", "description"]

    def test_read_default( self ):
        """
        What do we pick when nothing is defined.
        """
        df = ri.parse_rosetta_file(self.silent1)

        assert list(df.columns.values) == self.defhead
        assert list(df.shape) == [6, len(self.defhead)]
        assert df["score"].mean() == pytest.approx(-207.9, 0.2)
        assert df["packstat"].mean() == pytest.approx(0.59, 0.02)

    def test_read_json( self ):
        df = ri.parse_rosetta_json(self.jsonscr)
        assert df.shape == (88, 39)

    def test_read_pdb( self ):
        """
        Read a PDB file.
        """
        df = []
        for x in range(1, 4):
            df.append(ri.parse_rosetta_pdb(os.path.join(self.dirpath, 'INPUT_000{}.pdb'.format(x))))
        df = pd.concat(df).reset_index(drop=True)
        assert df.shape == (3, 23)
        assert 'label' not in df
        assert 'sequence_A' in df

        df2 = ri.parse_rosetta_pdb(os.path.join(self.dirpath, 'INPUT_0001.pdb'), keep_weights=True)
        assert df2.shape == (2, 24)
        assert 'label' in df2

    def test_read_score_and_complete( self ):
        """
        Read a score file and add sequences from the PDB.
        """
        df = ri.parse_rosetta_file(os.path.join(self.dirpath, 'INPUT_score.sc'))
        assert df.shape == (3, 53)
        assert 'sequence_A' not in df

        df = df.retrieve_sequences_from_pdbs(prefix=self.dirpath)
        assert df.shape == (3, 55)
        assert 'sequence_A' in df
        assert 'sequence_B' in df

    def test_select_scores( self ):
        """
        Select only some scores of interest.
        """
        # pick some
        sc_des = {"scores": ["score", "rama", "omega", "packstat", "rmsd_drift"]}
        df = ri.parse_rosetta_file(self.silent1, sc_des)
        assert list(df.columns.values) == sc_des["scores"]
        assert df["packstat"].mean() == pytest.approx(0.59, 0.02)

        # pick all
        sc_des = {"scores": "*"}
        df = ri.parse_rosetta_file(self.silent1, sc_des)
        assert list(df.columns.values) == self.defhead

    def test_ignore_scores( self ):
        """
        Use the ignore_scores parameter.
        """
        # ignore point columns
        defhead = self.defhead[:]
        sc_des = {"scores_ignore": ["dslf_fa13", "rama", "omega", "fa_dun"]}
        fl_des = os.path.join(self.dirpath, 'description_ignore.yaml')
        for x in sc_des["scores_ignore"]:
            defhead.remove(x)

        df = ri.parse_rosetta_file(self.silent1, fl_des)
        assert list(df.columns.values) == defhead
        assert df["packstat"].mean() == pytest.approx(0.59, 0.02)
        with pytest.raises(KeyError):
            df["dslf_fa13"]

        # ignore by widlcard
        defhead = [x for x in self.defhead if not x.startswith("fa_")]
        sc_des = {"scores_ignore": ["fa_*"]}
        df = ri.parse_rosetta_file(self.silent1, sc_des)
        assert list(df.columns.values) == defhead
        with pytest.raises(KeyError):
            df["fa_dun"]

    def test_rename_scores( self ):
        """
        Rename scores into something different.
        """
        sc_des = os.path.join(self.dirpath, 'description_rename.json')
        defhead = self.defhead[:]
        defhead[defhead.index("packstat")] = "inscore"
        defhead[defhead.index("rama")] = "dingong"
        df = ri.parse_rosetta_file(self.silent1, sc_des)
        assert len(df.columns.values) == len(self.defhead)
        assert list(df.columns.values) == defhead
        with pytest.raises(KeyError):
            assert df["packstat"].mean() == pytest.approx(0.59, 0.02)
        assert df["inscore"].mean() == pytest.approx(0.59, 0.02)

    def test_scores_by_residue( self ):
        """
        Pick score_by_residue data.
        """
        # when non-requested, do not pick per-residue data
        df = ri.parse_rosetta_file(self.silent3)
        assert set(df.columns.values) == set(self.symhead)

        # one can request individual positions, though
        sc_des = {"scores": ["residue_ddg_66"]}
        df = ri.parse_rosetta_file(self.silent3, sc_des)
        assert len(df.columns.values) == 1
        assert df.iloc[0]["residue_ddg_66"] == pytest.approx(-3.47, 0.02)

        # get all per-residue values
        sc_des = {"scores": "-", "scores_by_residue": ["residue_ddg_"]}
        df = ri.parse_rosetta_file(self.silent3, sc_des)
        assert len(df.columns.values) == 1

        # request wrong per-residue value
        with pytest.raises(AttributeError):
            sc_des = {"scores_by_residue": ["residue_score_"]}
            df = ri.parse_rosetta_file(self.silent3, sc_des)

    def test_naming( self ):
        """
        Generate new data columsn from the design's description.
        """
        sc_des = {"naming": ["", "source", "", "status", "dcount"]}
        df = ri.parse_rosetta_file(self.silent1, sc_des)
        assert len(df.columns.values) == len(self.defhead) + 3
        assert df.iloc[0]["status"] == "labeled"
        assert df["dcount"].mean() == 3.5

    def test_sequence_data( self ):
        """
        Load data from sequence, structure or psipred
        """
        with pytest.raises(ValueError):
            sc_des = {"sequence": "A", "structure": "A"}
            df = ri.parse_rosetta_file(self.silent2, sc_des)

        sc_des = {"sequence": "C", "structure": "C"}
        df = ri.parse_rosetta_file(self.silent2, sc_des)
        assert len(df.columns.values) == 24
        assert len(df["sequence_C"]) == len(df["structure_C"])

        sc_des = {"scores_ignore": "*", "sequence": "C", "structure": "C"}
        df = ri.parse_rosetta_file(self.silent2, sc_des)
        assert len(df.columns.values) == 2
        assert len(df["sequence_C"]) == len(df["structure_C"])

    def test_read_labels( self ):
        """
        Check how labels are read and loaded.
        """
        sc_des = {"labels": ["MOTIF", "CONTACT", "CONTEXT"]}
        motif   = ("43-64", "A:#(0),B:#(22)")
        contact = ("9-26,28-29,31-32,35,37-40,67-68,70-71,89,91-135", "A:#(19),B:#(58)")
        context = ("117-273", "A:#(157),B:#(0)")
        df = ri.parse_rosetta_file(self.silent1, sc_des)

        defhead = self.defhead + ["lbl_MOTIF", "lbl_CONTACT", "lbl_CONTEXT"]
        assert list(df.columns.values) == defhead

        # Check the label values and types
        s = df.iloc[0]
        assert str(s["lbl_MOTIF"]) == motif[1]
        assert str(s["lbl_CONTACT"]) == contact[1]
        assert str(s["lbl_CONTEXT"]) == context[1]
        assert isinstance(s["lbl_MOTIF"], rc.SelectionContainer)

        # Check the internal values
        assert s["lbl_CONTACT"]["B"] == "9-26,28-29,31-32,35,37-40,67-68,70-71,89,91-116"
        assert s["lbl_CONTACT"]["A"] == "1-19"

        # Check that data is not repeated when it should not be
        f = df.iloc[1]
        assert cmp(s["lbl_MOTIF"], f["lbl_MOTIF"]) == 0
        assert cmp(s["lbl_CONTACT"], f["lbl_CONTACT"]) != 0
        assert cmp(s["lbl_CONTEXT"], f["lbl_CONTEXT"]) == 0

    def test_symmetry( self ):
        """
        Check on sequence capture with symmetry silent files.
        """
        sc_des = {"sequence": "AB"}
        df = ri.parse_rosetta_file(self.silent3, sc_des)

        assert set(df.columns.values) == set(self.symhead + ["sequence_A", "sequence_B"])
        assert df.iloc[0]["sequence_A"] == df.iloc[0]["sequence_B"]

    def test_motifgraft( self ):
        """
        Check reading from a 2-motif MotifGraftMover output with missing values.
        """
        sc_des = {'graft_ranges': 2,
                  'scores_missing': ['rama_per_res_filter'],
                  'scores_ignore': ['graft_out_scaffold_ranges']}
        df = ri.parse_rosetta_file(self.silent4, sc_des)

        assert df.columns[df.isna().any()].tolist() == ['rama_per_res_filter']
        assert 'graft_out_scaffold_ranges01' not in df.columns
        assert 'graft_scaffold_size_change01' in df.columns
        assert len(df.columns) == 45

    def test_contacts( self ):
        """
        Check reading contact files.
        """
        df, irow, icol = ri.parse_rosetta_contacts(self.contact)
        assert df.shape[0] == df.shape[1]
        assert irow == icol
