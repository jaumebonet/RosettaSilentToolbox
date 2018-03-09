# @Author: Jaume Bonet <bonet>
# @Date:   13-Dec-2017
# @Email:  jaume.bonet@gmail.com
# @Filename: Selection.py
# @Last modified by:   bonet
# @Last modified time: 09-Mar-2018


import copy
from itertools import groupby, count
import re

from pandas import Series
import numpy as np


class Selection( object ):
    """
    :py:class:`.Selection` allows for the complex management of residue selection
    from a sequence. It can be used in any function that accepts the `key_residue`
    parameter.

    It accepts both numerical and string data. Thus, a :py:class:`.Selection`
    can be declared::

        # From an array of numbers
        sn = Selection([3, 4, 5, 13, 14, 15, 21, 25])

        # From a string representation of numbers
        ss = Selection("3-5,13-15,21,25")

        # From a string representation of PDB numbering.
        sp = Selection("4A-6A,14A-16A,22A,26A")

    If a :py:class:`~pandas.Series` is provided, :py:class:`.Selection` will try to
    extract appropiate content.
    If a PDB numbering schema is provided, :py:class:`.Selection` will consider that
    the shift of the original PDB is already taken into account, and will ignore the
    reference shift set into other components of the library. Instead it will make sure
    that it is called for the appropiate `seqID`.

    Current limitations of :py:class:`.Selection` include:

    #. Cannot combine PDB and regular numbering.
    #. Cannot add PDB numbering from multiple chains.

    Multiple operations are available for Selection::

        # Negate selection: will call 'select all except'
        not_sele = ~sele

        # Addition and substract will join or find the difference between
        # two selections.
        new_sele = sele1 - sele2  # Res in sele1 not in sele2
        new_sele = sele1 + sele2  # Res in both sele1 and sele2

        # Logical operations
        new_sele = sele1 & sele2  # Res in sele1 that are also in sele2
        new_sele = sele1 | sele2  # Res in both sele1 and sele2

        # Shift
        new_sele = sele << 2  # Shift all residue selection by -2
        new_sele = sele >> 2  # Shift all residue selection by +2

    :param selection: This can be a list of integers with the residues
        of interest or a string. If it is the last one, residue numbers
        have to be separated by space and ranges by '-'. Both things can
        be mixed.
    :type selecton: Union[:py:class:`str`, :py:class:`list`, :py:class:`~pandas.Series`]

    :raises:
        :AttributeError: Whenever the input does not meet the previously mentioned
        limitations.
        :AttributeError: If the provided selection object cannot be processed.
        :ValueError: If values provided cannot be properly converted to an integer list.

    """
    def __init__( self, selection=None ):

        self._asarr = []     # Selected Residues.
        self._seqID = None   # Sequence ID; if present do NOT apply shift.
        self._revrs = False  # Select all but Selection.

        if isinstance(selection, Series):
            if selection.shape[0] == 1:
                selection = selection.values[0]
            else:
                selection = list(selection.values)
        if isinstance(selection, basestring):
            if len(selection) > 0:
                self._asarr = self._string_to_list(selection)
        elif isinstance(selection, (list, np.ndarray)):
            self._asarr = sorted(list(set(selection)))
        elif selection is None:
            pass
        else:
            raise AttributeError("Unable to processs the provided selection")

    def to_list( self ):
        """
        Provide the values of the :py:class:`.Selection` as a
        list of integers.

        :return: :py:class:`list`(:py:class:`int`)
        """
        return self._asarr

    def to_string( self ):
        """
        Provide the values of the :py:class:`.Selection` as a
        string.

        :return: :py:class:`str`
        """
        return self._list_to_string()

    def seqID( self ):
        """
        Identifier of the sequence to which :py:class:`.Selection`
        is assigned.

        :return: :py:class:`str`
        """
        return self._seqID

    def is_empty( self ):
        """
        Evaluate if :py:class:`.Selection` is empty.

        :return: :py:class:`bool`
        """
        return len(self._asarr) == 0

    def is_shifted( self ):
        """
        Evaluate if :py:class:`.Selection` is shifted.

        :return: :py:class:`bool`
        """
        return bool(self._seqID)

    def map_to_sequences( self, sequence_map ):
        """
        Mainly called through the :py:func:`.parse_rosetta_file`.

        For a :py:class:`.Selection` without `seqID`, a `sequence_map`
        can be provided. A `sequence_map` is an array as long as all the
        sequences available with the `seqID` per position. These are
        created when parsing Rosetta files. It will return a dict with
        :py:class:`.Selection`s for each available `seqID`. The
        :py:class:`.Selection`s won't have `seqID`, as the shift is not
        added.
        If a `seqID` is present in the `sequence_map` but has no residue
        selected, an empty :py:class:`.Selection` will be generated.

        :param sequence_map: List with `seqID` per position
        :type sequence_map: :py:class:`list`(:py:class:`str`)

        :return: :py:class:`dict`{:py:class:`str`: :py:class:`.Selection`}

        :raises:
            :KeyError: If the :py:class:`.Selection` has a `seqID`
            :IndexError: If position exceeds the sequencee map
        """
        if self._seqID is not None:
            raise KeyError("The Selection has already an assigned sequence id")

        s = list(Series(self.to_list()).apply(lambda x: sequence_map[x - 1]))
        x = SelectionContainer()
        for n, c in zip(self, s):
            x.setdefault(c, []).append(n - sequence_map.index(c))
        for c in x:
            x[c] = Selection(x[c])
        for _ in set(sequence_map).difference(x):
            x[_] = Selection()
        return x

    def shift( self, seqID, shift ):
        """
        Shifts the :py:class:`.Selection` according to a shift and sets
        up the belonging chain.

        :param seqID: Identifier of the reference sequence
        :type seqID: :py:class:`str`
        :param value: Identifier of the reference sequence
        :type value: Union[:py:class:`int`, :py:class:`list`(:py:class:`int`)]

        :return: New shifted :py:class:`.Selection`.
        """
        newsele = Selection()
        if isinstance(shift, int):
            newsele = self >> (shift - 1)
        if isinstance(shift, list) and len(self) > 0:
            newsele._asarr = [shift[x - 1] for x in self._asarr]
        newsele._seqID = seqID
        return newsele

    def unshift( self, seqID, shift ):
        """
        Unshifts the :py:class:`.Selection` according to a shift and sets
        up the belonging chain.

        :param seqID: Identifier of the reference sequence
        :type seqID: :py:class:`str`
        :param value: Identifier of the reference sequence
        :type value: Union[:py:class:`int`, :py:class:`list`(:py:class:`int`)]

        :return: New shifted :py:class:`.Selection`.
        """
        newsele = Selection()
        if isinstance(shift, int):
            newsele = self << (shift - 1)
        if isinstance(shift, list) and len(self) > 0:
            newsele._asarr = [shift.index(x) + 1 for x in self._asarr]
        newsele._seqID = seqID
        return newsele

    #
    # PRIVATE METHODS
    #
    def _string_to_list( self, selection ):
        """
        Will transform the string definition inside the object to an array.

        :param selection: Representation of the selection positions
        :type selection: :py:class:`str`

        :return: :py:class:`list`(:py:class:`int`)

        :raises:
            :AttributeError: If more than one seqID is provided.
        """
        if re.search('[a-zA-Z]', selection):
            seqID = set(re.sub(r'[0-9\-\,\s]+', '', selection))
            if len(seqID) > 1:
                raise AttributeError("More than one chain ID is provided.")
            self._seqID = seqID.pop()

        o = []
        for x in selection.split(","):
            if "-" not in x:
                o.append(self._evaluate_number(x))
            else:
                xx = x.split("-")
                for i in range(self._evaluate_number(xx[0]), self._evaluate_number(xx[1]) + 1):
                    o.append(i)
        return sorted(list(set(o)))

    def _evaluate_number( self, number ):
        """
        Return integer value from a string taking into account possible
        seqID assignation.

        :return: :py:class:`int`

        :raises:
            :AttributeError: If it finds a seqID that does not match with
            the expected ones.
        """
        number = number.strip()
        if self._seqID is None:
            return int(number)

        if number[-1] != self._seqID:
            raise AttributeError("Mixed PDB and regular count.")
        return int(number[:-1])

    def _list_to_string( self ):
        """
        Will transform the list of selected residues into a string definition.

        :return: :py:class:`str`
        """
        def as_range(iterable, seqID):
            nums = list(iterable)
            if len(nums) > 1:
                return '{0}{2}-{1}{2}'.format(nums[0], nums[-1], seqID)
            else:
                return '{0}{1}'.format(nums[0], seqID)

        if len(self) == 0:
            return ""
        value = groupby(self._asarr, key=lambda n, c=count(): n - next(c))
        seqID = self._seqID if self._seqID is not None else ""
        return ','.join(as_range(g, seqID) for _, g in value)

    def _compressed_str( self ):
        if self.is_shifted():
            return "@({})".format(len(self))
        else:
            return "#({})".format(len(self))

    #
    # MAGIC METHODS
    #
    def __str__( self ):
        return self._list_to_string()

    def __repr__( self ):
        return self._list_to_string()

    def __iter__( self ):
        return iter(self._asarr)

    def __len__( self ):
        return len(self._asarr)

    def __invert__( self ):
        s = copy.deepcopy(self)
        s._revrs = True
        return s

    def __hash__( self ):
        mark = "~" if self._revrs else ""
        return mark + self._list_to_string()

    def __eq__( self, other ):
        if isinstance(other, Selection):
            if self._seqID != other._seqID:
                raise KeyError("Cannot compare Selections with different seqID")
            ilen = len(set(self._asarr).intersection(other._asarr))
            return (ilen == len(self)) and (self._revrs == other._revrs)
        if isinstance(other, (Series, basestring, list)):
            return self == Selection(other)
        raise NotImplementedError

    def __ne__( self, other ):
        return not self.__eq__(other)

    def __lshift__( self, other ):
        if isinstance(other, int):
            s = Selection(list(np.array(self._asarr) - other))
            s._seqID = self._seqID
            return s
        raise NotImplementedError

    def __rshift__( self, other ):
        if isinstance(other, int):
            s = Selection(list(np.array(self._asarr) + other))
            s._seqID = self._seqID
            return s
        raise NotImplementedError

    def __or__( self, other ):
        return self.__add__(other)

    def __and__( self, other ):
        if isinstance(other, Selection):
            if self._seqID != other._seqID:
                raise KeyError("Cannot add Selections with different seqID")
            if self._revrs == other._revrs:
                s = Selection(list(set(self.to_list()).intersection(other.to_list())))
            else:
                s = Selection(list(set(self.to_list()).difference(other.to_list())))
            s._seqID = self._seqID
            s._revrs = self._revrs
            return s
        if isinstance(other, int):
            return self + Selection([other, ])
        if isinstance(other, (Series, basestring, list)):
            return self + Selection(other)
        raise NotImplementedError

    def __add__( self, other ):
        if isinstance(other, Selection):
            if self._seqID != other._seqID:
                raise KeyError("Cannot add Selections with different seqID")
            if self._revrs == other._revrs:
                s = Selection(list(set(self.to_list()).union(other.to_list())))
            else:
                s = Selection(list(set(self.to_list()).difference(other.to_list())))
            s._seqID = self._seqID
            s._revrs = self._revrs
            return s
        if isinstance(other, int):
            return self + Selection([other, ])
        if isinstance(other, (Series, basestring, list)):
            return self + Selection(other)
        raise NotImplementedError

    def __sub__( self, other ):
        if isinstance(other, Selection):
            if self._seqID != other._seqID:
                raise KeyError("Cannot add Selections with different seqID")
            if self._revrs == other._revrs:
                s = Selection(list(set(self.to_list()).difference(other.to_list())))
            else:
                s = Selection(list(set(self.to_list()).union(other.to_list())))
            s._seqID = self._seqID
            s._revrs = self._revrs

            return s
        if isinstance(other, int):
            return self - Selection([other, ])
        if isinstance(other, (Series, basestring, list)):
            return self - Selection(other)
        raise NotImplementedError


class SelectionContainer( object ):
    """
    Helper class to manage representation of selectors in pandas.
    """
    def __init__( self, *args ):
        self._content = dict(args)

    def shift( self, seqID, value ):
        """
        Helper to ease the apply function. Shifts by value the labels
        assigned to a given seqID.

        :param seqID: Identifier of the reference sequence
        :type seqID: :py:class:`str`
        :param value: Identifier of the reference sequence
        :type value: Union[:py:class:`int`, :py:class:`list`(:py:class:`int`)]

        :raises:
            :ValueError: If the :py:class:`.Selection` is already shifted
        """
        if seqID in self:
            if self[seqID].is_shifted():
                raise ValueError("Selection is alreay shifted.")
            self[seqID] = self[seqID].shift(seqID, value)

    def setdefault( self, k, d ):
        self._content.setdefault(k, d)
        return self._content[k]

    def __getitem__( self, key ):
        return self._content[key]

    def __setitem__(self, key, value):
            self._content[key] = value

    def __iter__( self ):
        return self._content.__iter__()

    def __cmp__( self, other ):
        if isinstance(other, SelectionContainer):
            return cmp(self._content, other._content)
        else:
            raise NotImplementedError

    def __str__( self ):
        return ",".join(["{0}:{1}".format(x, self[x]._compressed_str()) for x in sorted(self)])

    def __repr__( self ):
        return str(self)
