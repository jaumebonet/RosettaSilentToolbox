# @Author: Jaume Bonet <bonet>
# @Date:   19-Feb-2018
# @Email:  jaume.bonet@gmail.com
# @Filename: sequence.py
# @Last modified by:   bonet
# @Last modified time: 06-Mar-2018


import os
import gzip

import rstoolbox.core as core
import rstoolbox.components as cp
from rstoolbox.io.rosetta import _gather_file_list


def read_fasta( filename, split_char=None, split_position=1, multi=False ):
    """
    Reads one or more fasta files and returns the appropiate object
    containing the requested data: the :py:class:`.DesignFrame`.

    For the purposes of the :py:class:`.DesignFrame`, the ``seqID`` assigned
    to the read sequences will be ``A``. The identifier of the sequences will
    be stored in the ``description`` column to be coherent with Rosetta reads.

    :param filename: file name or file pattern to search.
    :type filename: :py:class:`str`
    :param split_char: Split the design ID according to the given char.
    :type split_char: :py:class:`str`
    :param split_position: If split is not py:data:`None`, which position of the
        split is the name? The rest will go to a column named **info**.
    :param multi: When :py:data:`True`, indicates that data is readed from
        multiple files.
    :type multi: :py:class:`bool`

    :return: :py:class:`.DesignFrame`.

    :raises:
        :IOError: if ``filename`` cannot be found.
    """
    files = _gather_file_list( filename, multi )
    data = {"description": [], "sequence_A": [], "info": []}
    for file_count, f in enumerate( files ):
        fd = gzip.open( f ) if f.endswith(".gz") else open( f )
        for line in fd:
            line = line.strip()
            if line.startswith(">"):
                line = line.strip(">")
                if split_char is None:
                    data["description"].append(line)
                else:
                    line = line.split(split_char)
                    data["description"].append(line.pop(split_position - 1).replace(".pdb", ""))
                    data["info"].append(split_char.join(line))
                data["sequence_A"].append("")
            elif len(line) > 0:
                data["sequence_A"][-1] += line

    if len(data["info"]) == 0:
        del(data["info"])
    df = cp.DesignFrame( data )
    df.add_source_files( files )
    return df


def write_fasta( df, seqID, filename ):
    """
    Writes fasta files of the selected decoys.

    It assumes that the provided data is contained in a :py:class:`.DesignFrame`
    or a :py:class:`~pandas.DataFrame` with a **description** column for the
    decoy's names and a **sequence_[seqID]** column for the proper sequence.

    :param df: Data content.
    :type df: Union[:py:class:`.DesignFrame`, :py:class:`~pandas.DataFrame`]
    :param seqID: Identifier(s) of the sequences expected to be printed.
    :type seqID: :py:class:`str`
    :param filename: output file name.
    :type filename: :py:class:`str`

    :raises:
        :IOError: if ``filename`` exists and global option *system.overwrite* \
        is not :py:data:`True`.
        :AttributeError: if requested seqID cannot be found.
    """
    if os.path.isfile(filename) and not core.get_option("system", "overwrite"):
        raise IOError("File {} already exists".format(filename))

    data = []
    for chain in seqID:
        query = "sequence_{}".format(chain)
        name  = "description"
        if query not in df:
            raise AttributeError("seqID {} not found in data".format(chain))
        eachfa = df.apply(lambda row: ">" + row[name] + "_" + chain + "\n" + row[query], axis=1)
        data.append("\n".join(list(eachfa)))

    fd = open(filename, "w") if not filename.endswith(".gz") else gzip.open(filename, "wb")
    fd.write("\n".join(data))
    fd.close()
