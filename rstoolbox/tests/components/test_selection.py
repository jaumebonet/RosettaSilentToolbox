# @Author: Jaume Bonet <bonet>
# @Date:   27-Feb-2018
# @Email:  jaume.bonet@gmail.com
# @Filename: test_selection.py
# @Last modified by:   bonet
# @Last modified time: 06-Mar-2018


import pytest
import numpy as np

import rstoolbox.components as rc


class TestSelection( object ):

    def test_input_array( self ):
        """
        Check how Selection works when an array is given
        """
        arr = [2, 3, 4, 5, 15, 21, 22, 23, 24, 68, 72]
        shw = "2-5,15,21-24,68,72"

        s = rc.Selection(arr)
        assert s.to_list() == arr
        assert s.to_string() == shw
        assert s.is_empty() == False
        assert s.is_shifted() == False

    def test_input_string( self ):
        """
        Check how Selection works when a string is given
        """
        arr = [2, 3, 4, 5, 15, 21, 22, 23, 24, 68, 72]
        shw = "2-5,15,21-24,68,72"

        s = rc.Selection(shw)
        assert s.to_list() == arr
        assert s.to_string() == shw
        assert s.seqID() is None
        assert s.is_empty() == False
        assert s.is_shifted() == False

        s = rc.Selection("5")
        assert s.to_list() == [5]
        assert s.to_string() == "5"
        assert s.is_empty() == False
        assert s.is_shifted() == False

    def test_input_empty( self ):
        """
        Check how Selection works when nothing is given
        """
        arr = []
        shw = ""

        s = rc.Selection(arr)
        assert s.to_list() == arr
        assert s.to_string() == shw
        assert s.is_empty() == True
        assert s.is_shifted() == False

        s = rc.Selection(shw)
        assert s.to_list() == arr
        assert s.to_string() == shw
        assert s.is_empty() == True
        assert s.is_shifted() == False

    def test_input_PDBcount( self ):
        """
        Check how selection works when a PDB count string
        is provided
        """
        arr = [2, 3, 4, 5, 15, 21, 22, 23, 24, 68, 72]
        shw = "2A-5A,15A,21A-24A,68A,72A"

        s = rc.Selection(shw)
        assert s.to_list() == arr
        assert s.to_string() == shw
        assert s.seqID() == "A"
        assert s.is_empty() == False
        assert s.is_shifted() == True

        s = rc.Selection("5A")
        assert s.to_list() == [5]
        assert s.to_string() == "5A"
        assert s.is_empty() == False
        assert s.is_shifted() == True

        with pytest.raises(AttributeError):
            shw = "2A-5A,15B,21A-24A,68A,72A"
            s = rc.Selection(shw)

            shw = "2A-5A,15,21A-24A,68A,72A"
            s = rc.Selection(shw)

    def test_operations( self ):

        a1 = [2, 3, 4, 5, 15, 21, 22, 23, 24, 68, 72]
        _12 = [2, 3, 4, 5, 6, 15, 21, 22, 23, 24, 68, 72]
        _13 = [2, 3, 4, 5, 21, 22, 23, 24, 68, 72]
        _14 = [1, 2, 3, 4, 14, 20, 21, 22, 23, 67, 71]
        _15 = [3, 4, 5, 6, 16, 22, 23, 24, 25, 69, 73]
        a2 = [2, 15, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27]
        _21 = [2, 3, 4, 5, 15, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 68, 72]
        _22 = [3, 4, 5, 68, 72]
        _23 = [2, 15, 21, 22, 23, 24]
        t1 = "2A-5A,15A,21A-24A,68A,72A"
        t2 = "2A,15A-18A,21A-27A"
        s1 = rc.Selection(a1)
        s2 = rc.Selection(a2)
        s3 = rc.Selection(t1)
        s4 = rc.Selection(t2)

        # equal
        assert s1 == a1
        assert s1 != a2

        # invert
        sn1 = ~s1
        assert len(set(sn1.to_list()).intersection(s1.to_list())) == len(s1)
        assert sn1 != s1

        # add, sub
        s1a = s1 + 6   # operate w/ integer
        assert s1a == _12
        s1a = s1 - 15  # operate w/ integer
        assert s1a == _13
        s2a = s1 + s2  # operate w/ Selection
        assert s2a == _21
        s2a = s1 - s2  # operate w/ Selection
        assert s2a == _22
        s2a = s1a + a2  # operate w/ list
        assert s2a == _21
        s2a = s1a - a2  # operate w/ list
        assert s2a == _22
        s3a = s3 + s4   # operate w/ Selection w/ seqID
        assert s3a.to_list() == _21
        s3a = s3 - s4   # operate w/ Selection w/ seqID
        assert s3a.to_list() == _22
        with pytest.raises(KeyError):
            s3a = s3 + 6  # operate w/ in to Select w/ seqID
        with pytest.raises(KeyError):
            s3a = s3 + a2  # operate w/ in to Select w/ seqID
        s3a = s3 + "6A"   # opearte w/ string
        assert s3a.to_list() == _12
        s3a = s3 - "15A"  # opearte w/ string
        assert s3a.to_list() == _13
        # operate with invert
        assert (s1 + s2) == (s1 - ~s2)
        assert (s1 - s2) == (s1 + ~s2)

        # shift
        s4a = s1 << 1
        assert s4a == _14
        s4a = s1 >> 1
        assert s4a == _15

        # logical operations
        s5a = s1 & s2
        assert s5a == _23
        s5a = s1 | s2
        assert s5a == _21

    def test_map_to_sequences( self ):
        a1 = [2, 3, 4, 5, 15, 21, 22, 23, 24, 50, 51, 68, 72, 110, 111, 112, 113]
        a2 = [2, 3, 210]
        _1 = [2, 3, 4, 5, 15, 21, 22, 23, 24, 50]
        _2 = list(np.array([51, 68, 72]) - 50)
        _3 = list(np.array([110, 111, 112, 113]) - 100)
        t1 = "2A-5A,15A,21A-24A,68A,72A, 110A-113A"

        s1 = rc.Selection(a1)
        s2 = rc.Selection(t1)
        s3 = rc.Selection(a2)

        smap = ["A", ] * 50
        smap.extend(["B", ] * 50)
        smap.extend(["C", ] * 50)
        smap.extend(["D", ] * 50)
        cntr = "A:#(10),B:#(3),C:#(4),D:#(0)"

        s1a = s1.map_to_sequences(smap)
        assert s1a["A"] == _1
        assert s1a["B"] == _2
        assert s1a["C"] == _3
        assert str(s1a) == cntr

        # shift assumes first position is 1, not 0
        s1a.shift("A", 5)
        s1a.shift("B", 2)
        assert s1a["A"].to_list() == list(np.array(_1) + (5 - 1))
        assert s1a["B"].to_list() == list(np.array(_2) + (2 - 1))
        assert s1a["C"] == _3
        assert str(s1a) == cntr

        # Selection w/ seqID should raise an error
        with pytest.raises(KeyError):
            s2.map_to_sequences(smap)

        # Going out of range raises error
        with pytest.raises(IndexError):
            s3.map_to_sequences(smap)
