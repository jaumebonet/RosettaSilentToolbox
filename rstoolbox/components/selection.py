# @Author: Jaume Bonet <bonet>
# @Date:   13-Dec-2017
# @Email:  jaume.bonet@gmail.com
# @Filename: Selection.py
# @Last modified by:   bonet
# @Last modified time: 13-Dec-2017


from itertools import groupby, count

class Selection(object):
    """
    The selection object allows to provide selected key residues
    to analyze.
    Transform it to a string will show the residues in a manner similar
    to Rosetta count. The representation of it (__repr__), will show it with
    the sequence identifier, like a PDB count.

    :param str seqID: Identifier of the sequence sets of interest.
    :param selection: This can be a list of integers with the residues
        of interest or a string. If it is the last one, residue numbers
        have to be separated by space and ranges by '-'. Both things can
        be mixed.
    """
    def __init__(self, seqID, selection):
        super(Selection, self).__init__()
        self._seqID    = seqID
        self._asstring = ""
        self._aslist   = []
        if isinstance(selection, basestring):
            self._asstring = selection
            self._aslist   = self._string_to_list()
        elif isinstance(selection, list):
            self._aslist   = sorted(list(set(selection)))
            self._asstring = self._list_to_string()

    @property
    def seqID(self):
        """
        Identifier of the sequence of interest
        :return: str
        """
        return self._seqID

    @seqID.setter
    def seqID(self, value):
        """
        Identifier of the sequence of interest
        :param str value: seqID
        """
        self._seqID = value

    def _string_to_list(self):
        """
        Will transform the string definition inside the object to an array.
        :return: array list of integers.
        """
        r = self._asstring.split(",")
        o = []
        for x in r:
            if "-" not in x:
                o.append(int(x))
            else:
                xx = x.split("-")
                for i in range(int(xx[0]),  int(xx[1]) + 1):
                    o.append(i)
        return o

    def _list_to_string(self, seqID=""):
        """
        Will transform the list of selected residues into a string definition.
        :return: str
        """
        def as_range(iterable, seqID):
            l = list(iterable)
            if len(l) > 1:
                return '{0}{2}-{1}{2}'.format(l[0], l[-1], seqID)
            else:
                return '{0}{1}'.format(l[0], seqID)

        if len(self) == 0:
            return ""
        if len(self) == 1:
            return str(self._aslist[0])
        return ','.join(as_range(g, seqID) for _, g in groupby(self._aslist, key=lambda n, c=count(): n-next(c)))

    def __str__(self):
        return self._asstring

    def __repr__(self):
        return self._list_to_string(self._seqID)

    def __iter__(self):
        return iter(self._aslist)

    def __len__(self):
        return len(self._aslist)
