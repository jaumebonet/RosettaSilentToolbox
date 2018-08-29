# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. class:: FragmentFrame
"""
# Standard Libraries
import copy
from itertools import groupby, count
import re

# External Libraries
import six
from pandas import Series
import numpy as np

# This Library


__all__ = ["Selection", "SelectionContainer", "get_selection"]


def get_selection( key_residues, seqID, shift=1, length=None ):
    """*Internal function*; global management and casting of :class:`.Selection`.

    Call this function in any function with the ``key_residues`` parameter. It will
    manage the type of input and return the appropiate list of selected residues.

    :param key_residues: |keyres_param|.
    :type key_residues: |keyres_types|
    :param str seqID: |seqID_param|.
    :param shift: Starting residue number or per-residue number assignment.
    :type shift: Union[:class:`int`, :func:`list` of :class:`int`]
    :param int length: Length of the sequence over which it will be applied.
        Needed in case of inverted :class:`.Selection`.

    :return: :class:`~numpy.ndarray`

    :raises:
        :NotImplementedError: If ``key_residues`` is of a non-expected type.
    """
    listtypes = (list, np.ndarray, six.string_types)
    if six.PY3:
        listtypes = (list, np.ndarray, six.string_types, range)

    if key_residues is None:
        if isinstance(shift, list):
            key_residues = Selection(shift)
            key_residues._seqID = seqID
        elif length is not None:
            key_residues = range(1, length + 1)
    if isinstance(key_residues, int):
        key_residues = [key_residues, ]
    if isinstance(key_residues, listtypes):
        key_residues = Selection(key_residues)
    if isinstance(key_residues, SelectionContainer):
        key_residues = key_residues[seqID]
    if isinstance(key_residues, Selection):
        if key_residues.is_shifted():
            kr = key_residues.unshift(seqID, shift)
        else:
            kr = key_residues.unshift(None, 1)
        kr = np.array(kr.to_list(length))
    else:
        raise NotImplementedError
    return kr


class Selection( object ):
    """Complex management of residue selection from a sequence.

    It can be used in any function that accepts the ``key_residue`` parameter.

    It accepts both numerical and string data.
    Thus, a :class:`.Selection` can be declared in multiple ways.

    .. ipython::

        In [1]: from rstoolbox.components import Selection
           ...: # From an array of numbers
           ...: sn = Selection([3, 4, 5, 13, 14, 15, 21, 25])
           ...: # From a string representation of numbers
           ...: ss = Selection("3-5,13-15,21,25")
           ...: # From a string representation of PDB numbering.
           ...: sp = Selection("4A-6A,14A-16A,22A,26A")

    If a :class:`~pandas.Series` is provided, :class:`.Selection` will try to
    extract the appropiate content.

    If **regular numbering** is provided, it will assume that it refers to direct
    sequence positions.

    .. note::
        :class:`.Selection` works with **sequence positioning**; that means that it
        expects the first position to be *1*, not *0*.

    If a **PDB numbering** schema is provided, :class:`.Selection` will consider that
    the *shift* of the original PDB is already taken into account and will correct accordingly
    when applied to the different functons and data containers.

    .. note::
        **PDB numbering** cannot combine selections from multiple chains.

    :param selection: Residue positions that define the sequence selection.
    :type selection: Union[:class:`str`, :func:`list`, :class:`~pandas.Series`]

    :raises:
        :AttributeError: if the provided selection object cannot be processed.
        :ValueError: if values provided cannot be properly converted to an integer list.

    Multiple operations are available for Selection.

    .. ipython::

        In [1]: sele = Selection([3, 4, 5, 13, 14, 15, 21, 25])
           ...: # Negate selection: will call 'select all except'
           ...: not_sele = ~sele
           ...: sele.to_list()

        In [1]: not_sele.to_list(25)


        In [1]: # Addition and substract will join or find the difference between
           ...: # two selections.
           ...: sele1 = Selection([3, 4, 5, 13, 14])
           ...: sele2 = Selection([14, 15, 21, 25])
           ...: sele1 - sele2  # Res in sele1 not in sele2

        In [1]: sele1 + sele2  # Res in both sele1 and sele2

        In [1]: # Logical operations
           ...: sele1 & sele2  # Res in sele1 that are also in sele2

        In [1]: new_sele = sele1 | sele2  # Res in both sele1 and sele2

        In [1]: # Shift
           ...: sele << 2  # Shift all residue selection by -2

        In [1]: sele >> 2  # Shift all residue selection by +2


    .. warning::

        *Shifted* and *Unshifted* :class:`.Selection` cannot operate between them.

    """
    def __init__( self, selection=None ):

        self._asarr = []     # Selected Residues.
        self._seqID = None   # Sequence ID; if present do NOT apply shift.
        self._revrs = False  # Select all but Selection.
        self._isarr = None   # Reverse selection, if needed
        self._ialen = None   # Length used on reversed

        listtypes = (list, np.ndarray)
        if six.PY3:
            listtypes = (list, np.ndarray, range)

        if isinstance(selection, Series):
            if selection.shape[0] == 1:
                selection = selection.values[0]
            else:
                selection = list(selection.values)
        if isinstance(selection, six.string_types):
            if len(selection) > 0:
                self._asarr = self._string_to_list(selection)
        elif isinstance(selection, listtypes):
            self._asarr = sorted(list(set(selection)))
        elif selection is None:
            pass
        elif isinstance(selection, Selection):
            self._asarr = selection._asarr
            self._seqID = selection._seqID
            self._revrs = selection._revrs
            self._isarr = selection._isarr
            self._ialen = selection._ialen
        else:
            raise AttributeError("Unable to processs the provided selection")

    def to_list( self, length=None ):
        """Provide the values of the :class:`.Selection` as a list of integers.

        :param length: Expected total length of the sequence to which the
            :class:`.Selection` will be applied.
        :type length: :class:`int`

        :return: (:func:`list` of :class:`int`) -- Numerical list representation
            of the selection positions.

        :raises:
            :AttributeError: If the :class:`.Selection` is reversed and no
                ``length`` is provided.

        If the :class:`.Selection` **is reversed**, the ``length`` over which it will
        be applied needs to be provided, so that the actual selected positions can
        be determined.

        .. ipython::

            In [2]: sr = ~ss
               ...: ss.to_list(25)

        .. rubric:: Example

        .. ipython::

            In [1]: from rstoolbox.components import Selection
               ...: ss = Selection("3-5,13-15,21,25")
               ...: ss.to_list()
        """
        if not self._revrs:
            return self._asarr
        else:
            if self._isarr is None and length is None:
                raise AttributeError("Reversed Selections need to know the "
                                     "sequence length")
            if self._isarr is None or self._ialen != length:
                self._isarr = [x for x in range(1, length + 1) if x not in self._asarr]
                self._ialen = length
            return sorted(self._isarr)

    def to_string( self ):
        """Provide the values of the :class:`.Selection` as a string.

        :return: :class:`str`

        .. rubric:: Example

        .. ipython::

            In [1]: from rstoolbox.components import Selection
               ...: ss = Selection([3, 4, 5, 13, 14, 15, 21, 25])
               ...: ss.to_string()

        If the :class:`.Selection` is :meth:`~Selection.is_shifted`, it
        will write in **PDB annotation**.

        .. ipython::

            In [1]: ss.shift("A", 3).to_string()


        """
        return self._list_to_string()

    def seqID( self ):
        """Identifier of the sequence to which :class:`.Selection` is assigned.

        A :class:`.Selection` only has ``seqID`` if :meth:`~Selection.is_shifted`,
        otherwise ``seqID`` is :data:`None`.

        :return: :class:`str`

        .. rubric:: Example

        .. ipython::

            In [1]: from rstoolbox.components import Selection
               ...: ss = Selection([3, 4, 5, 13, 14, 15, 21, 25])
               ...: ss.seqID()

            In [1]: ss.shift("A", 3).seqID()
        """
        return self._seqID

    def is_empty( self ):
        """Evaluate if :class:`.Selection` is empty.

        :return: :class:`bool`

        .. rubric:: Example

        .. ipython::

            In [1]: from rstoolbox.components import Selection
               ...: ss = Selection([3, 4, 5, 13, 14, 15, 21, 25])
               ...: ss.is_empty()

            In [1]: Selection().is_empty()
        """
        return len(self._asarr) == 0

    def is_shifted( self ):
        """Evaluate if :class:`.Selection` is shifted.

        Is the selection is shifted, it will need to be assigned to a seqID.

        :return: :class:`bool`

        .. rubric:: Example

        .. ipython::

            In [1]: from rstoolbox.components import Selection
               ...: ss = Selection([3, 4, 5, 13, 14, 15, 21, 25])
               ...: ss.is_shifted()

            In [1]: ss.shift("A", 3).is_shifted()
        """
        return bool(self._seqID)

    def map_to_sequences( self, sequence_map ):
        """Generator for :func:`.parse_rosetta_file`.

        This function is not really expected to be directly accessed by the user.

        For a :class:`.Selection` without ``seqID``, a ``sequence_map``
        can be provided. A ``sequence_map`` is an array as long as all the
        sequences available with the ``seqID`` per position. These are
        created when parsing Rosetta files. It will return a dict with
        :class:`.Selection` for each available ``seqID``. The
        :class:`.Selection` won't have ``seqID``, as the shift is not
        added.

        .. note::
            If a ``seqID`` is present in the ``sequence_map`` but has no residue
            selected, an empty :class:`.Selection` will be generated for that key.

        :param sequence_map: List with ``seqID`` per position.
        :type sequence_map: :func:`list` of :class:`str`

        :return: :class:`.SelectionContainer`

        :raises:
            :KeyError: If the :class:`.Selection` has a ``seqID``.
            :IndexError: If position exceeds the sequencee map.

        .. seealso::
            :func:`.parse_rosetta_file`

        .. warning::
            It cannot be applied to an already shifted :class:`.Selection`.

        .. rubric:: Example

        .. ipython::

            In [1]: from rstoolbox.components import Selection
               ...: seq = ["A",] * 14 + ["B",] * 11
               ...: "".join(seq)

            In [1]: ss = Selection([3, 4, 5, 13, 14, 15, 21, 25])
               ...: sd = ss.map_to_sequences(seq)
               ...: for seqID in sd:
               ...:     print(seqID, sd[seqID])
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
        """Shifts the :class:`.Selection` according to a value and sets
        up the chain to which it is associated.

        There are *two* ways in which a shift can be provided:

        * An :class:`int` specifying the numerical identity of the first position \
            of the sequence. As :class:`.Selection` works with **sequence positioning**, \
            shifting by 1 is the same as not shifting at all (although then the ``seqID`` \
            will be added and the :class:`.Selection` will be considered as shifted).

        .. ipython::

            In [1]: from rstoolbox.components import Selection
               ...: ss = Selection([3, 4, 5, 13, 14, 15, 21, 25])
               ...: ss.shift("A", 3)

            In [1]: ss.shift("A", 1)

        * A :func:`list` of :class:`int` specifying the identity of each position of the \
            sequence. This is usefull if your reference sequence has gaps in its numbering.

        .. ipython::

            In [1]: from rstoolbox.components import Selection
               ...: ss = Selection([3, 4, 5, 13, 14, 15, 21, 25])
               ...: seq = range(1, 31)
               ...: seq = list(seq[:14]) + list(seq[19:])
               ...: ",".join([str(_) for _ in seq])

            In [1]: ss.shift("A", seq)

        :param str seqID: |seqID_param|.
        :param shift: Expected displacement
        :type shift: Union[:class:`int`, :func:`list` of :class:`int`]

        :return: :class:`.Selection` - new shifted selection

        .. seealso::
            :meth:`~Selection.unshift`
        """
        newsele = Selection()
        if isinstance(shift, int):
            newsele = self >> (shift - 1)
        if isinstance(shift, list) and len(self) > 0:
            newsele._asarr = [shift[x - 1] for x in self._asarr]
        newsele._seqID = seqID
        newsele._revrs = self._revrs
        return newsele

    def unshift( self, seqID, shift ):
        """Unhifts the :class:`.Selection` according to a value.

        Inverst the shift in the :class:`.Selection` to allow for
        direct sequence position targeting.

        There are *two* ways in which the unshift can be provided:

        * An :class:`int` specifying the numerical identity of the first position \
            of the sequence. As :class:`.Selection` works with **sequence positioning**, \
            unshifting by 1 is the same as not shifting at all (although then the ``seqID`` \
            will be removed and the :class:`.Selection` will be considered not shifted).

        .. ipython::

            In [1]: from rstoolbox.components import Selection
               ...: ss = Selection([3, 4, 5, 13, 14, 15, 21, 25])
               ...: sf = ss.shift("A", 3)
               ...: sf

            In [1]: sf.unshift("A", 3)

            In [1]: sf.unshift("A", 1)

        * A :func:`list` of :class:`int` specifying the identity of each position of the \
            sequence. This is usefull if your reference sequence has gaps in its numbering.

        .. ipython::

            In [1]: from rstoolbox.components import Selection
               ...: ss = Selection([3, 4, 5, 13, 14, 15, 21, 25])
               ...: seq = range(1, 31)
               ...: seq = list(seq[:14]) + list(seq[19:])
               ...: ",".join([str(_) for _ in seq])

            In [1]: sf = ss.shift("A", seq)
               ...: sf

            In [1]: sf.unshift("A", seq)

        :param str seqID: |seqID_param|.
        :param shift: Expected displacement.
        :type shift: Union[:class:`int`, :func:`list` of :class:`int`]

        :return: New unshifted :class:`.Selection`.

        .. seealso::
            :meth:`~Selection.shift`
        """
        newsele = Selection()
        if isinstance(shift, int):
            newsele = self << (shift - 1)
        if isinstance(shift, list) and len(self) > 0:
            newsele._asarr = [shift.index(x) + 1 for x in self._asarr]
        newsele._seqID = None
        newsele._revrs = self._revrs
        return newsele

    #
    # PRIVATE METHODS
    #
    def _string_to_list( self, selection ):
        """Will transform the string definition inside the object to an array.

        :param selection: Representation of the selection positions
        :type selection: :class:`str`

        :return: :class:`list` of class:`int`

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
        """Return integer value from a string taking into account possible
        ``seqID`` assignation.

        :return: :class:`int`

        :raises:
            :AttributeError: |seqID_error|.
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

        :return: :class:`str`
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
            show = "@({})".format(len(self))
        else:
            show = "#({})".format(len(self))
        return "~" + show if self._revrs else show

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
        if isinstance(other, (Series, six.string_types, list)):
            return self == Selection(other)
        raise NotImplementedError

    def __ne__( self, other ):
        return not self.__eq__(other)

    def __lshift__( self, other ):
        if isinstance(other, int):
            s = Selection(list(np.array(self._asarr) - other))
            s._seqID = self._seqID
            s._revrs = self._revrs
            return s
        raise NotImplementedError

    def __rshift__( self, other ):
        if isinstance(other, int):
            s = Selection(list(np.array(self._asarr) + other))
            s._seqID = self._seqID
            s._revrs = self._revrs
            return s
        raise NotImplementedError

    def __or__( self, other ):
        return self.__add__(other)

    def __and__( self, other ):
        return _operate(self, other, "intersection", "difference", "__add__")

    def __add__( self, other ):
        return _operate(self, other, "union", "difference", "__add__")

    def __sub__( self, other ):
        return _operate(self, other, "difference", "union", "__sub__")


class SelectionContainer( object ):
    """Helper class to manage representation of selectors in :mod:`pandas`.

    A :class:`.SelectionContainer` is generated when labels are read through
    the :func:`.parse_rosetta_file`, as **ResidueLabels** are saved in
    **RosettaNumbering**.

    Basically it is just a :class:`dict` mimic that allows to quickly access
    the :class:`.Selection` data avoiding raising errors if :class:`.Selection`
    for particular ``seqID`` are not present and wraps the shifting functions.

    The other main function is to minimize the representation of :class:`.Selection`
    in the printed :class:`~pandas.DataFrame`. Thus, its way of being represented by
    the length of each contained :class:`.Selection`:

    .. ipython::

        In [1]: from rstoolbox.components import Selection, SelectionContainer
           ...: sc = SelectionContainer()
           ...: sc.setdefault("A", Selection([3, 4, 5, 13, 14, 15, 21, 25]))
           ...: sc.setdefault("B", Selection("15-19,21-25"))
           ...: sc

    Representation of the individual :class:`.Selection` changes if it
    :meth:`~.Selection.is_shifted`:

    .. ipython::

        In [1]: sc.shift("A", 3)
           ...: sc

    And when it is reversed:

    .. ipython::

        In [1]: sc["A"] = ~sc["A"]
           ...: sc

    .. seealso::
        :func:`.parse_rosetta_file`
    """
    def __init__( self, *args ):
        self._content = dict(args)

    def shift( self, seqID, value ):
        """Shifts by value the labels assigned to a given ``seqID``.

        Helper to ease the apply function.

        :param str seqID: |seqID_param|.
        :param value: Identifier of the reference sequence
        :type value: Union[:class:`int`, :class:`list` of :class:`int`]

        :raises:
            :ValueError: if the :class:`.Selection` is already shifted

        .. seealso::
            :meth:`Selection.shift`
        """
        if seqID in self:
            if self[seqID].is_shifted():
                raise ValueError("Selection is alreay shifted.")
            self[seqID] = self[seqID].shift(seqID, value)

    def unshift( self, seqID, value ):
        """Unshifts by value the labels assigned to a given seqID.

        Helper to ease the apply function.

        :param str seqID: |seqID_param|.
        :param value: Identifier of the reference sequence
        :type value: Union[:class:`int`, :class:`list` of :class:`int`]

        :raises:
            :ValueError: if the :class:`.Selection` is already shifted

        .. seealso::
            :meth:`Selection.unshift`
        """
        if seqID in self:
            if self[seqID].is_shifted():
                self[seqID] = self[seqID].unshift(None, value)

    def setdefault( self, k, d ):
        self._content.setdefault(k, d)
        return self._content[k]

    def __getitem__( self, key ):
        return self._content[key]

    def __setitem__(self, key, value):
            self._content[key] = value

    def __iter__( self ):
        return self._content.__iter__()

    def __len__( self ):
        return len(self._content)

    def __cmp__( self, other ):
        if isinstance(other, SelectionContainer):
            return cmp(self._content, other._content)
        else:
            raise NotImplementedError

    def __eq__( self, other ):
        if isinstance(other, SelectionContainer):
            return self._content == other._content
        else:
            raise NotImplementedError

    def __str__( self ):
        return ",".join(["{0}:{1}".format(x, self[x]._compressed_str()) for x in sorted(self)])

    def __repr__( self ):
        return str(self)


def _operate( self, other, func1, func2, final ):
    if isinstance(other, Selection):
        if self._seqID != other._seqID:
            raise KeyError("Cannot operate Selections with different seqID")
        if self._revrs == other._revrs:
            s = Selection(list(getattr(set(self.to_list()), func1)(other.to_list())))
        else:
            s = Selection(list(getattr(set(self.to_list()), func2)(other.to_list())))
        s._seqID = self._seqID
        s._revrs = self._revrs
        return s

    if isinstance(other, int):
        if final == "__sub__":
            return self - Selection([other, ])
        if final == "__add__":
            return self + Selection([other, ])
    if isinstance(other, (Series, six.string_types, list)):
        if final == "__sub__":
            return self - Selection(other)
        if final == "__add__":
            return self + Selection(other)
    raise NotImplementedError
