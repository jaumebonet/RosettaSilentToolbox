# @Author: Jaume Bonet <bonet>
# @Date:   05-Mar-2018
# @Email:  jaume.bonet@gmail.com
# @Filename: test_design.py
# @Last modified by:   bonet
# @Last modified time: 12-Mar-2018


import os
import copy

import pandas as pd
import numpy as np
import pytest

import rstoolbox.io as ri
import rstoolbox.components as rc


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

        # Let's work with an array-type shift
        ashift = range(1, len(refseq) + 1)
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
