# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: open_rosetta_file
.. func:: parse_rosetta_file
.. func:: parse_rosetta_json
.. func:: parse_rosetta_contacts
.. func:: parse_rosetta_fragments
.. func:: write_rosetta_fragments
.. func:: write_fragment_sequence_profiles
.. func:: get_sequence_and_structure
.. func:: make_structures
"""
# Standard Libraries
import os
import sys
import re
import glob
import gzip
import json
import string
import shutil
from collections import OrderedDict

# External Libraries
import six
import pandas as pd
import yaml

# This Library
import rstoolbox.core as core
import rstoolbox.components as rc
from rstoolbox.utils import baseline, make_rosetta_app_path, execute_process

__all__ = ['open_rosetta_file', 'parse_rosetta_file', 'parse_rosetta_contacts',
           'parse_rosetta_fragments', 'write_rosetta_fragments',
           'write_fragment_sequence_profiles', 'get_sequence_and_structure',
           'make_structures', 'parse_rosetta_json']

_headers = ["SCORE", "REMARK", "RES_NUM", "FOLD_TREE", "RT",
            "ANNOTATED_SEQUENCE", "NONCANONICAL_CONNECTION",
            "SYMMETRY_INFO", "CHAIN_ENDINGS"]


def _file_vs_json( data ):
    """
    Transform file into json if needed.
    """
    if data is None:
        return {}
    if isinstance( data, str ):
        if not os.path.isfile( data ):
            raise IOError("{0}: file not found.".format(data))
        try:
            fd = gzip.open( data ) if data.endswith(".gz") else open( data )
            text = "".join([x.strip() for x in fd])
            data = json.loads(text)
        except ValueError:
            fd = gzip.open( data ) if data.endswith(".gz") else open( data )
            data = yaml.safe_load(fd.read())
        fd.close()
    return data


def _check_type( value ):
    """
    Makes sure that, upon reading a value, it gets assigned
    the correct type.
    """
    try:
        int(value)
    except ValueError:
        try:
            float(value)
        except ValueError:
            return value
        else:
            return float(value)
    else:
        return int(value)


def _gather_file_list( filename, multi=False ):
    """
    Provided a file name or pattern, generates a list
    with all the files that are expected to be read.

    :param str filename: file name or pattern (without "*"). If filename
        ends with ``$`` it is assumed that that is the end of the name.
    :param bool multi: Tell if a file name or pattern is provided.
        Default is 'False' (single file name)

    :return: list of str names of files
    """
    files = []
    if isinstance(filename, list):
        multi = True
    if not multi:
        if not os.path.isfile( filename ):
            raise IOError("{0}: file not found.".format(filename))
        files.append( filename )
    else:
        if isinstance(filename, six.string_types):
            if filename.endswith('$'):
                filename = filename.rstrip('$')
            elif not filename.endswith("*"):
                filename = filename + "*"
            files = glob.glob( filename )
        else:
            files = filename
    if len(files) == 0:
        raise IOError("Input pattern did not find any file.")
    for _ in files:
        if not os.path.isfile( _ ):
            raise IOError("{0}: file not found.".format(_))
    return files


def _add_sequences( manager, data, chains ):
    """
    Fix and load requested sequence/structure into the data.
    Also takes labels into account.

    :return: data
    """
    # Correct by non polymer residues
    nonPoly = [x for x, v in enumerate(chains["seq"]) if v == 'Z']
    ochaini = chains["id"]
    for index in sorted(nonPoly, reverse=True):
        del chains["id"][index]
        del chains["seq"][index]
        if len(chains["dssp"]) > 0:
            del chains["dssp"][index]

    for seqname, seq in manager.get_expected_sequences( chains ):
        data.setdefault( seqname, [] ).append( seq )
    if len(chains["dssp"]) > 0:
        for ssename, str3d in manager.get_expected_structures( chains ):
            data.setdefault( ssename, [] ).append( str3d )
    if len(chains["psipred"]) > 0:
        for ssename, str3d in manager.get_expected_psipred( chains ):
            data.setdefault( ssename, [] ).append( str3d )
    if len(chains["phi"]) > 0:
        for ssename, str3d in manager.get_expected_dihedrals( chains, "phi" ):
            data.setdefault( ssename, [] ).append( str3d )
    if len(chains["psi"]) > 0:
        for ssename, str3d in manager.get_expected_dihedrals( chains, "psi" ):
            data.setdefault( ssename, [] ).append( str3d )

    for x in [i for i in data if i.startswith("lbl_")]:
        data[x][-1] = rc.Selection(data[x][-1]).map_to_sequences(ochaini)

    return data


def open_rosetta_file( filename, multi=False, check_symmetry=True ):
    """
    *Internal function*; reads through a Rosetta silent file and yields only
    the lines that the library knows how to parse.

    For each *"parsable"* line, it yields 4 values:

    ===== ============= ====================================
    order     type       content
    ===== ============= ====================================
        1 :class:`str`  data of the line
        2 :class:`bool` is line header?
        3 :class:`int`  name of the readed file
        4 :class:`bool` does the file contain symmetry info?
    ===== ============= ====================================

    :param filename: file name, file pattern to search or list of files.
    :type filename: Union[:class:`str`, :func:`list`]
    :param multi: Tell if a file name (single file) or pattern (multifile) is provided.
    :type multi: :class:`bool`
    :param check_symmetry: Check if the silent file contains symmetry info.
    :type check_symmetry: :class:`bool`

    :yields: Union[:class:`str`, :class:`bool`, :class:`int`, :class:`bool`]

    :raises:
        :IOError: if ``filename`` cannot be found.
        :IOError: if ``filename`` pattern (``multi=True``) generates no files.

    .. seealso:
        :func:`parse_rosetta_file`
    """
    symm = False
    files = _gather_file_list( filename, multi )
    for file_count, f in enumerate( files ):
        if check_symmetry:
            cmd = "zgrep SYMMETRY_INFO {}" if f.endswith(".gz") else "grep SYMMETRY_INFO {}"
            process = execute_process(cmd.format(f))
            symm = (process == 0)
        fd = gzip.open( f ) if f.endswith(".gz") else open( f )
        for line in fd:
            line = line.decode('utf8') if f.endswith(".gz") else line
            if line.strip().split()[0].strip(":") in _headers:
                yield line, line.strip().split()[-1] == "description", file_count, symm
        fd.close()


def parse_rosetta_file( filename, description=None, multi=False ):
    """Read a Rosetta score or silent file and returns the design population
    in a :class:`.DesignFrame`.

    By default, it will pick the data contained in **all the score columns**
    with the exception of positional scores (such as *per-residue ddg*). The
    user can specify scores to be ignored.

    When working with *silent files*, extra information can be picked, such as
    *sequence* and *secondary structure* data, *residue labels* or positional
    scores. The fine control of these options is explained in detail in
    :ref:`tutorial: reading Rosetta <readrosetta>`.

    Some basic usage cases::

        # (1) The default scenario, just read scores from a single file.
        df = rstoolbox.io.parse_rosetta_file("silentfile")

        # (2) Reading from multiple files. Assumes all files start with
        # the particular prefix.
        df = rstoolbox.io.parse_rosetta_file("silentfile", multi=True)

        # (3) Getting all scores and the sequence of each design.
        description = {'sequence': 'A'}
        df = rstoolbox.io.parse_rosetta_file("silentfile", description)

        # (4) Get only total_score and RMSD, and rename total_score to score.
        description = {'scores': ['RMSD'], 'scores_rename': {'total_score': 'score'}}
        df = rstoolbox.io.parse_rosetta_file("silentfile", description)

    :param filename: file name, file pattern to search or list of files.
    :type filename: Union[:class:`str`, :func:`list`]
    :param description: Parsing rules. It can be a dictionary describing
        the rules or the name of a file containing such dictionary. The
        dictionary definition is explained in :ref:`tutorial: reading Rosetta <readrosetta>`.
    :type description: Union[:class:`str`, :class:`dict`]
    :param bool multi: When :data:`True`, indicates that data is readed from multiple files.

    :return: :class:`.DesignFrame`.

    :raises:
        :IOError: if ``filename`` cannot be found.
        :IOError: if ``filename`` pattern (``multi=True``) generates no files.

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_file("../rstoolbox/tests/data/input_2seq.minisilent.gz")
           ...: df.head(2)
    """

    manager = rc.Description( **_file_vs_json( description ) )
    header  = []
    data    = OrderedDict()

    for line, is_header, _, symm in open_rosetta_file( filename, multi ):
        if is_header:
            header = line.strip().split()[1:]
            continue

        if line.startswith("SCORE"):
            per_res = {}
            chains  = {"id": [], "seq": "", "dssp": "", "psipred": "", "phi": [], "psi": []}

            # General scores
            for cv, value in enumerate( line.strip().split()[1:] ):
                hcv = header[cv]
                if manager.wanted_per_residue_score( hcv ):
                    hcvn = re.sub("\d+$", "", hcv)
                    per_res.setdefault( hcvn, {} )
                    per_res[hcvn][int(re.findall('\d+$', hcv)[0])] = _check_type( value )
                    continue
                if manager.wanted_score( hcv ):
                    data.setdefault( manager.score_name( hcv), [] ).append( _check_type( value ) )

            # Namings from the description
            manager.check_naming( header )
            for namingID, namingVL in manager.get_naming_pairs( line.strip().split()[-1] ):
                data.setdefault( namingID, [] ).append( _check_type( namingVL ) )

            # Fix per-residue
            for k in per_res:
                data.setdefault( k, [] ).append( OrderedDict(sorted(per_res[k].items())).values() )

            # Setup labels
            data = manager.setup_labels( data )
            continue

        if line.startswith("RES_NUM"):  # In multichains and not starting in A1.
            for x in line.split()[1:-1]:
                chain, numbers = x.split(":")
                nums = numbers.split("-")
                if len(nums) == 1 or nums[0] == "":
                    nums = 1
                else:
                    nums = (int(nums[1]) - int(nums[0])) + 1
                chains["id"].extend([chain, ] * nums)
            continue

        if line.startswith("SYMMETRY_INFO"):  # When working with symmetry, RES_NUM is not there...
            chain = "".join(string.ascii_uppercase[:int(line.split()[2])])
            for c in chain:
                chains["id"].extend([c, ] * int(line.split()[4]))

            data = _add_sequences( manager, data, chains )
            continue

        if line.startswith("ANNOTATED_SEQUENCE"):
            chains["seq"] = list(re.sub( r'\[[^]]*\]', '', line.strip().split()[1] ))
            if not symm:
                # When info is chain A starting in 1, it is not printed in the silent file
                if len(chains["id"]) == 0:
                    chains["id"].extend(["A", ] * len(chains["seq"]))

                data = _add_sequences( manager, data, chains )
            else:
                chains["seq"] = list("".join(chains["seq"]).rstrip("X"))

            continue

        if line.startswith("REMARK DSSP"):
            chains["dssp"] = list(line.split()[2].strip())
            continue
        if line.startswith("REMARK PSIPRED"):
            chains["psipred"] = list(line.split()[2].strip())
            continue
        if line.startswith("REMARK LABELS"):
            for label in line.split()[2].split(";"):
                labinfo = label.split(":")
                if "lbl_" + labinfo[0].upper() in data:
                    data["lbl_" + labinfo[0].upper()][-1] = labinfo[1]
            continue
        if line.startswith("REMARK PHI"):
            chains["phi"] = [float(x) for x in line.split()[2].strip().split(",")]
            continue
        if line.startswith("REMARK PSI"):
            chains["psi"] = [float(x) for x in line.split()[2].strip().split(",")]
            continue

    df = rc.DesignFrame( data )
    df.add_source_files( _gather_file_list( filename, multi ) )
    return df


def parse_rosetta_json( filename ):
    """Read a json formated rosetta score file.

    Only reads back scores, as those are the only content present in a ``JSON`` file.

    :param str filename: File containing the Rosetta score file.

    :return: :class:`.DesignFrame`.

    .. note::
        To be coherent with the silent files, the decoy id column name ``decoy`` is
        changed to ``description``.

    :raises:
        :IOError: if ``filename`` cannot be found.

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_json
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = parse_rosetta_json("../rstoolbox/tests/data/score.json.gz")
           ...: df.head(2)
    """
    is_gz = filename.endswith(".gz")
    fd = gzip.open( filename ) if is_gz else open( filename )
    data = {}
    for line in fd:
        if is_gz:
            dt = json.loads(line.decode('utf8').strip())
        else:
            dt = json.loads(line.strip())
        for k in dt:
            data.setdefault(k, []).append(dt[k])
    df = rc.DesignFrame( data )
    df.rename(columns={'decoy': 'description'})
    df.add_source_file( filename )
    return df


def parse_rosetta_contacts( filename ):
    """Read a residue contact file as generated by **ContactMapMover**.

    Returns three objects:

    ===== ============================ =================================================
    order                         type                                           content
    ===== ============================ =================================================
        1 :class:`~pandas.DataFrame`   matrix with boolean :data:`True` in
                                       contacts; **Rosetta numbering** (no ``seqID``)
        2 :func:`list` of :class:`str` list of 3-letter code amino acids for row axis
        3 :func:`list` of :class:`str` list of 3-letter code amino acids for column axis
    ===== ============================ =================================================

    In a regular run for the **ContactMapMover** without selectors,
    list 2 and 3 will be identical.

    :param str filename: File containing the Rosetta fragments.

    :return: Union[:class:`~pandas.DataFrame`,
        :func:`list` of :class:`str`,
        :func:`list` of :class:`str`]

    :raises:
        :IOError: if ``filename`` cannot be found.
    """
    if not os.path.isfile(filename):
        raise IOError("{} not found!".format(filename))

    df    = pd.read_csv(filename, comment="#", delim_whitespace=True)
    rows  = [x[:3] for x in df.index.values]
    irows = [int(x[3:]) for x in df.index.values]
    cols  = [x[:3] for x in df.columns.values]
    icols = [int(x[3:]) for x in df.columns.values]

    df.columns = icols
    df.index = irows

    return df, rows, cols


def parse_rosetta_fragments( filename ):
    """Read a Rosetta fragment-file and return the appropiate :class:`.FragmentFrame`.

    It supports both old and new fragment formats.
    Does not support varying size fragment sets. When working with both *small* and
    *large* fragments sets, one would need to call the function twice and save the
    data of each in a different variable.

    :param str filename: File containing the Rosetta fragments.

    :return: :class:`.FragmentFrame`.

    :raises:
        :IOError: if ``filename`` cannot be found.

    .. seealso::
        :func:`.plot_fragments`
        :func:`.plot_fragment_profiles`
    """
    fformat = 1  # formats are identified as 1 and 0
    data = OrderedDict({"pdb": [], "frame": [], "neighbors": [], "neighbor": [], "position": [],
                        "size": [], "aa": [], "sse": [], "phi": [], "psi": [], "omega": []})

    if not os.path.isfile(filename):
        raise IOError("{} not found!".format(filename))

    fframe, fne, fpos, fsize = None, None, None, None
    fsaved_size = 0
    nei = 0
    was_space = False
    fd = gzip.open( filename ) if filename.endswith(".gz") else open( filename )
    for line in fd:
        line = line.decode('utf8') if filename.endswith(".gz") else line
        line = line.strip()
        if line == "":
            fpos = int(fframe)
            if fsize != 0:
                fsaved_size = fsize
            fsize = 0
            if not was_space:
                nei += 1
            was_space = True
            continue
        line = line.split()
        if line[0] == "FRAME":
            fformat = 1
            fframe = line[1]
            fpos = int(fframe)
            fsize = 0
            nei = 0
        elif line[0] == "position:":
            fformat = 0
            fframe = line[1]
            fne = line[-1]
            nei = -1
        else:
            if fformat == 1:
                data["pdb"].append(str(line[2]))
            else:
                data["pdb"].append(str(line[0]))
            data["frame"].append(int(fframe))
            data["neighbors"].append(fne)
            data["neighbor"].append(nei + 1)
            data["position"].append(int(line[0] if bool(fformat) else fpos))
            data["aa"].append(line[3 + fformat])
            data["sse"].append(line[4 + fformat])
            data["phi"].append(float(line[5 + fformat]))
            data["psi"].append(float(line[6 + fformat]))
            data["omega"].append(float(line[7 + fformat]))

            fpos += 1
            fsize += 1
        was_space = False
    fd.close()
    data["size"] = [fsaved_size, ] * len(data["frame"])

    df = rc.FragmentFrame(data, file=filename)
    if bool(fformat):
        df2 = df.groupby(["frame", "size"]).size().reset_index(name="neighbors")
        df2["neighbors"] = (df2["neighbors"] / df2["size"]).astype(int)
        df = df.merge(df2, on=["frame", "size"])
        df = df.drop(["neighbors_x"], axis=1)
        df = df.rename({"neighbors_y": "neighbors"}, axis=1)
    return df.reindex(["pdb", "frame", "neighbors", "neighbor", "position", "size",
                       "aa", "sse", "phi", "psi", "omega"], axis=1)


def write_rosetta_fragments( df, frag_size, n_frags=200 ):
    """Writes a Rosetta fragment-file (new format) from an appropiate :class:`.FragmentFrame`.

    Supports varying size fragment sets.

    :param df: Selected set of fragments that have to be written.
    :type df: :class:`.DesignFrame`
    :param int frag_size: Size of the fragments.
    :param int n_frags: Number of fragments per frame.
        Default are 200 fragments per frame.
    """
    _STRING = " {:4s} {:1s} {:5d} {:1s} {:1s} {:8.3f} {:8.3f} {:8.3f}\n"
    _HEADER = "position:            {} neighbors:          {}\n\n"
    with open("rosetta_frags.{}mers".format(frag_size), "w") as f:
        frame_count = 0
        for i in range(len(df)):
            if i % ((frag_size * n_frags)) == 0:
                frame_count += 1
                f.write(_HEADER.format(frame_count, n_frags))
            f.write(_STRING.format( df.loc[i]["pdb"], "X", int(0), df.loc[i]["aa"],
                    df.loc[i]["sse"], df.loc[i]["phi"], df.loc[i]["psi"], df.loc[i]["omega"]) )
            if i != 0 and (i + 1) % frag_size == 0:
                f.write("\n")


def write_fragment_sequence_profiles( df, filename=None, consensus=None ):
    """Write a sequence profile from :class:`.FragmentFrame` to load into
    Rosetta's **SeqprofConsensus**.

    Format mimicks as much as possible **BLAST PSSM**.

    :param df: Fragments from which to create the matrix or the matrix itself.
    :type df: Union[:class:`.FragmentFrame`, :class:`~pandas.DataFrame`]
    :param str filename: Output file name.
    :param str consensus: Consensus sequence to show.

    :return: :class:`str` - the expected file content if ``filename`` is :data:`None`.

    :raises:
        :IOError: if ``filename`` exists and :ref:`system.overwrite <options>` is
            :data:`False`.
        :ValueError: if ``consensus`` lenth differs from the expected by the
            given fragments.
    """
    def format_row(row, aa, row0):
        val1 = "".join(row.apply("{0:>3d}".format))
        val0 = "".join(row0.apply("{0:>4d}".format))
        data = "{0:>5d} {1}  ".format(row.name + 1, aa)
        return data + val1 + " " + val0 + \
            "{0:>6.2f}{0:>8.2f}".format(0)

    if isinstance(df, rc.FragmentFrame):
        matrix = df.make_sequence_matrix(round=True)
    else:
        matrix = df.copy()
    if consensus is None:
        consensus = df.quick_consensus_sequence()
    if len(consensus) != matrix.shape[0]:
        raise ValueError('Sequence need to be the same length.')
    matrix2 = matrix.copy()
    matrix2[:] = 0
    data = list(matrix.apply(lambda row: format_row(row, consensus[row.name],
                                                    matrix2.iloc[row.name]),
                             axis=1))
    head = "{0:>11}".format(" ") + \
           "  ".join(list(matrix.columns)) + "   " + "   ".join(list(matrix.columns))
    data.insert(0, head)
    prefix = "Last position-specific scoring matrix computed, weighted observed percentages " + \
             "rounded down, information per position, and relative weight of gapless real " + \
             "matches to pseudocounts"
    data.insert(0, prefix)
    data.insert(0, "")
    data = "\n".join(data)

    if filename is not None:
        if os.path.isfile(filename) and not core.get_option("system", "overwrite"):
            raise IOError("File {} already exists".format(filename))
        with open(filename, "w") as fd:
            fd.write(data)
    else:
        return data


def get_sequence_and_structure( pdbfile, mk_minisilent=True, ignore_unrecognized_res=True ):
    """Provided a PDB file, it will run a small **RosettaScript** to capture its sequence and
    structure, i.e. dssp and phi-psi dihedrals.

    .. note::
        Depends on :ref:`rosetta.path <options>` and :ref:`rosetta.compilation <options>`,
        if the corresponding silent file does not exist.

    .. attention::
        This function **REQUIRES** a local installation of **Rosetta**.

    It will generate an output file called: ``<pdbfile>.dssp.minisilent``. If this file
    already exists, it will be directly read. You can choose to compress it;
    ``<pdbfile>.dssp.minisilent.gz`` will also work.

    :param str pdbfile: Name of the input structure.
    :param bool mk_minisilent: If :data:`True`, transeform output into ``minisilent`` format.
    :param bool ignore_unrecognized_res: If :data:`True`, **Rosetta** ignores non-recognizable
        residues.

    :return: :class:`.DesignFrame`.

    :raises:
        :IOError: if ``pdbfile`` cannot be found.
        :IOError: if Rosetta executable cannot be found.
        :ValueError: if Rosetta execution fails

    The script that is runned (and that can be runned outside and then brough back to
    your computer) is:

    .. ipython::

        In [1]: from rstoolbox.utils import baseline
           ...: print(baseline())

    .. seealso::
        :func:`.baseline`
    """
    if not os.path.isfile( pdbfile ):
        raise IOError("Structure {} cannot be found".format(pdbfile))
    if mk_minisilent:
        minisilent = re.sub("\.pdb|\.cif$", "", re.sub("\.gz$", "", pdbfile)) + ".dssp.minisilent"
    else:
        minisilent = re.sub("\.pdb|\.cif$", "", re.sub("\.gz$", "", pdbfile)) + ".dssp.silent"

    if os.path.isfile(minisilent):
        return parse_rosetta_file(minisilent,
                                  {"sequence": "*", "structure": "*", "dihedrals": "*"})
    elif os.path.isfile(minisilent + ".gz"):
        return parse_rosetta_file(minisilent + ".gz",
                                  {"sequence": "*", "structure": "*", "dihedrals": "*"})

    sys.stdout.write("Generating file {}\n".format(minisilent))

    with open("dssp.xml", "w") as fd:
        fd.write(baseline())

    # Check rosetta executable & run
    exe = make_rosetta_app_path('rosetta_scripts')
    command = ['{0}', '-parser:protocol {1}', '-s {2}', '-out:file:silent {3}',
               '-ignore_zero_occupancy off']
    if ignore_unrecognized_res:
        command.append('-ignore_unrecognized_res')
    command = ' '.join(command)
    command = command.format( exe, "dssp.xml", pdbfile, str(os.getpid()) + "_" )
    sys.stdout.write("Running Rosetta\n")
    sys.stdout.write(command + "\n")
    error = execute_process( command )
    os.unlink("dssp.xml")
    if not bool(error):
        if os.path.isfile(str(os.getpid()) + "_"):
            sys.stdout.write("Execution has finished\n")
            if mk_minisilent:
                sys.stdout.write("Making minisilent\n")
                fd = open( minisilent, "w" )
                for line, _, _, _ in open_rosetta_file( str(os.getpid()) + "_" ):
                    fd.write( line )
                fd.close()
            else:
                sys.stdout.write("Keeping original silent output\n")
                shutil.copy( str(os.getpid()) + "_", minisilent)
            os.unlink(str(os.getpid()) + "_")
            return get_sequence_and_structure( pdbfile, mk_minisilent )
        else:
            raise ValueError("Execution has failed\n")
    else:
        raise ValueError("Execution has failed\n")


def make_structures( df, outdir=None, tagsfilename="tags", prefix=None, keep_tagfile=True ):
    """Extract the selected decoys (if any).

    .. note::
        Depends on :ref:`rosetta.path <options>` and :ref:`rosetta.compilation <options>`.
        Depends on :ref:`system.overwrite <options>` and :ref:`system.output <options>`.

    .. attention::
        This function **REQUIRES** a local installation of **Rosetta**.

    It basicall runs the ``extract_pdbs`` **Rosetta** application over the selected decoys.

    .. code-block:: bash

       extract_pdbs.linuxgccrelease -in:file:silent <pdb> -tags <selected>

    It requires the :class:`.DesignFrame` to have ``source_file`` attached identifying the
    silent files from which the data can be extracted. **minisilent files will not work here**.
    This should happen by default with the library, if one reads from actual silent files,
    but can be set up with::

        # (1) Read from a minisilent file that does not contain structural data: substitute
        ``df.replace_source_files(["file1", "file2", ])``
        # (2) Add files to a recently casted DesignFrame
        ``df.add_source_files(["file1", "file2", ])``

    :param df: Selected set of decoy that have to be extracted.
    :type df: :class:`.DesignFrame`
    :param str outdir: Directory in which to save the PDB files. If none is provided,
        it will be loaded from the :ref:`system.output <options>` global option.
    :param str tagsfilename: Name of the file containing the ids of the decoys of interest.
        It will be created in the ``outdir``. An previously existing file will not be
        overwritten if the global option :ref:`system.overwrite <options>` is :data:`False`.
    :param str prefix: If provided, a prefix is added to the PDB files.
    :param bool keep_tagfile: If :data:`True`, do not delete the tag file after using it.

    :raises:
        :ValueError: if the provided data does not have a **description** column.
        :ValueError: if the **description** column has repeated identifiers.
        :AttributeError: if silent files from where to extract structures are not found.
        :IOError: if the attached silent files do not exist.
        :IOError: when trying to overwrite the ``tagsfilename`` if *system.overwrite* is
            :py:data:`False`.
        :IOError: if the rosetta executable is not found. Depends on *rosetta.path* and
            *rosetta.compilation*.
    """
    # Check that the selection has at least one decoy
    if df.shape[0] == 0:
        sys.stdout.write("There are no decoys that fullfill the selection criteria.")
        return

    # Check that a column named "description", from which the IDs are (must be unique)
    column = "description"
    if column not in df.columns:
        raise ValueError("Identifiers of the decoys must be assigned to the column 'description'.")
    if True in df.duplicated(column).value_counts().index:
        raise ValueError("There are repeated identifiers. This might indicate the merging of files "
                         "with identical prefixes and will be an issue with extracting "
                         " the structures.")

    # Check that we have associated silent files to extract the data from
    if not isinstance(df, rc.DesignFrame) or len(df.get_source_files()) == 0:
        raise AttributeError("There are not source files from where to extract the structures.")
    sfiles = list(df.get_source_files())

    # Manage output directory
    outdir = outdir if outdir is not None else core.get_option("system", "output")
    if not os.path.isdir( outdir ):
        os.makedirs( outdir )
    if not outdir.endswith("/"):
        outdir += "/"

    # Manage tag file name
    tagsfilename = os.path.join( outdir, tagsfilename )
    if os.path.isfile( tagsfilename ) and not core.get_option("system", "overwrite"):
        raise IOError("Filename {0} already exists and cannot be overwrite.".format(tagsfilename))

    # Manage prefix
    if prefix is not None:
        outdir = os.path.join(outdir, prefix)

    # Check rosetta executable
    exe = make_rosetta_app_path('extract_pdbs')

    # Print the tag file
    df[[column]].to_csv( tagsfilename, index=False, header=False)
    if not os.path.isfile(tagsfilename):
        raise IOError("Something went wrong writing the file {0}".format(tagsfilename))

    # Run process
    sfiles = " ".join(sfiles)
    command = "{0} -in:file:silent {1} -in:file:tagfile {2} -out:prefix {3}"
    command = command.format( exe, sfiles, tagsfilename, outdir )
    sys.stdout.write("Executing Rosetta's extract_pdbs app\n")
    sys.stdout.write("(depending on the total number of decoys and how many have "
                     "been requested this might take a while...)\n")
    error = execute_process( command )
    if not bool(error):
        sys.stdout.write("Execution has finished\n")
    else:
        sys.stdout.write("Execution has failed\n")

    # Remove extra files if requested
    if not keep_tagfile:
        os.unlink( tagsfilename )
