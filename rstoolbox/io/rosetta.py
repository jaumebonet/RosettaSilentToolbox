import os
import sys
import re
import glob
import gzip
import json
import string
import subprocess
from collections import OrderedDict

import pandas as pd

import rstoolbox.core as core
import rstoolbox.components as cp

_headers = ["SCORE", "REMARK", "RES_NUM", "FOLD_TREE", "RT",
            "ANNOTATED_SEQUENCE", "NONCANONICAL_CONNECTION",
            "SYMMETRY_INFO", "CHAIN_ENDINGS"]
_per_residues = ["residue_ddg_"]

def _file_vs_json( data ):
    """
    Transform file into json if needed.
    """
    if data is None:
        return {}
    if isinstance( data, str ):
        if not os.path.isfile( data ):
            raise IOError("{0}: file not found.".format(data))
        fd = gzip.open( data ) if data.endswith(".gz") else open( data )
        text = "".join([x.strip() for x in fd])
        data = json.loads(text)
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

def _gather_file_list( filename,  multi=False ):
    """
    Provided a file name or pattern, generates a list
    with all the files that are expected to be read.
    :param str filename: file name or pattern (without "*")
    :param bool multi: Tell if a file name or pattern is provided.
        Default is 'False' (single file name)
    :return: list of str names of files
    """
    files = []
    if not multi:
        if not os.path.isfile( filename ):
            raise IOError("{0}: file not found.".format(filename))
        files.append( filename )
    else:
        if isinstance(filename, basestring):
            if not filename.endswith("*"):
                filename = filename + "*"
            files = glob.glob( filename )
        else:
            files = filename
    return files

def open_rosetta_file( filename, multi=False ):
    """
    Reads through a Rosetta score or silent file yielding only the lines
    that can be parsed by the rstoolbox.

    For each "parsable" line, 4 different values are provided:

    #. line content: :py:class:`str` with the actual data of the file.
    #. is header: a :py:class:`bool` identifying the provided line as header or not.
    #. file counter: a :py:class:`int` indicating which file is being read (for multi file input)
    #. is symmetry: a :py:class:`bool` indicating if the silent file contains symmetry info.

    :param str filename: file name or file pattern to search.
    :type filename: :py:class:`str`
    :param multi: Tell if a file name (single file) or pattern (multifile) is provided.
    :type multi: :py:class:`bool`

    :yields: Union[:py:class:`str`, :py:class:`bool`, :py:class:`int`, :py:class:`bool`]

    :raises:
        :IOError: if ``filename`` cannot be found.

    """
    files = _gather_file_list( filename, multi )
    for file_count, f in enumerate( files ):
        cmd = "zgrep SYMMETRY_INFO {} |wc" if f.endswith(".gz") else "grep SYMMETRY_INFO {} |wc"
        process = subprocess.Popen(cmd.format(f), stdout=subprocess.PIPE, shell=True)
        symm = int(process.communicate()[0].strip().split()[0]) > 0
        fd = gzip.open( f ) if f.endswith(".gz") else open( f )
        for line in fd:
            if line.strip().split()[0].strip(":") in _headers:
                yield line, line.strip().split()[-1] == "description", file_count, symm
        fd.close()

def parse_rosetta_file( filename, description=None, multi=False ):
    """
    Reads a Rosetta score or silent file and returns the appropiate object
    containing the requested data: the :py:class:`.DesignFrame`.

    The user can specify the columns of interest, change names on columns to
    facilitate merging with other data sets and request extra information such
    as sequence of the designs or residue labels.

    :param filename: file name or file pattern to search.
    :type filename: :py:class:`str`
    :param description: Parsing rules. It can be a dictionary describing
        the rules or the name of a file containing such dictionary.
    :type description: Union[:py:class:`str`, :py:class:`float`]
    :param multi: When :py:data:`True`, indicates that data is readed from
        multiple files.
    :type multi: :py:class:`bool`

    :return: :py:class:`.DesignFrame`.

    :raises:
        :IOError: if ``filename`` cannot be found.

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
    """

    manager = cp.Description( **_file_vs_json( description ) )
    header  = []
    data    = OrderedDict()

    for line, is_header, count, symm in open_rosetta_file( filename, multi ):
        if is_header:
            header = line.strip().split()[1:]
            continue

        if line.startswith("SCORE"):
            per_res = {}
            chains  = {"id": [], "seq": "", "dssp": "", "psipred": ""}

            # General scores
            for cv, value in enumerate( line.strip().split()[1:] ):
                if manager.wanted_per_residue_score( header[cv] ):
                    per_res.setdefault( re.sub("\d+$", "", header[cv]), {} )
                    per_res[re.sub("\d+$", "", header[cv])][int(re.findall('\d+$', header[cv])[0])] = _check_type( value )
                    continue
                if manager.wanted_score( header[cv] ):
                    data.setdefault( manager.score_name( header[cv]), [] ).append( _check_type( value ) )

            # Namings from the dscription
            manager.check_naming( header )
            for namingID, namingVL in manager.get_naming_pairs( line.strip().split()[-1] ):
                data.setdefault( namingID, [] ).append( _check_type( namingVL ) )

            # Fix per-residue
            for k in per_res:
                data.setdefault( k, [] ).append( OrderedDict(sorted(per_res[k].items())).values() )

            # Setup labels
            data = manager.setup_labels( data )
            continue

        if line.startswith("RES_NUM"): # In multichains and not starting in A1.
            for x in line.split()[1:-1]:
                chain, numbers = x.split(":")
                nums = numbers.split("-")
                if len(nums) == 1 or nums[0] == "": nums = 1
                else: nums = (int(nums[1]) - int(nums[0])) + 1
                chains["id"].extend([chain,] * nums)
            continue

        if line.startswith("SYMMETRY_INFO"): # When working with symmetry, RES_NUM is not there...
            chain = "".join(string.uppercase[:int(line.split()[2])])
            for c in chain:
                chains["id"].extend([c, ] * int(line.split()[4]))

            # Correct by non polymer residues
            nonPoly = [x for x, v in enumerate(chains["seq"]) if v == 'Z']
            for index in sorted(nonPoly, reverse=True):
                del chains["id"][index]
                del chains["seq"][index]
                if len(chains["dssp"]) > 0:
                    del chains["dssp"][index]

            for seqname, seq in manager.get_expected_sequences( chains ):
                data.setdefault( seqname, [] ).append( seq )
            if len(chains["dssp"]) > 0:
                for ssename, str in manager.get_expected_structures( chains ):
                    data.setdefault( ssename, [] ).append( str )
            if len(chains["psipred"]) > 0:
                for ssename, str in manager.get_expected_psipred( chains ):
                    data.setdefault( ssename, [] ).append( str )
            continue

        if line.startswith("ANNOTATED_SEQUENCE"):
            chains["seq"] = list(re.sub( r'\[[^]]*\]', '', line.strip().split()[1] ))
            if not symm:
                if len(chains["id"]) == 0: # When info is chain A starting in 1, it is not printed in the silent file
                    chains["id"].extend(["A",] * len(chains["seq"]))

                # Correct by non polymer residues
                nonPoly = [x for x, v in enumerate(chains["seq"]) if v == 'Z']
                for index in sorted(nonPoly, reverse=True):
                    del chains["id"][index]
                    del chains["seq"][index]
                    if len(chains["dssp"]) > 0:
                        del chains["dssp"][index]

                for seqname, seq in manager.get_expected_sequences( chains ):
                    data.setdefault( seqname, [] ).append( seq )
                if len(chains["dssp"]) > 0:
                    for ssename, str in manager.get_expected_structures( chains):
                        data.setdefault( ssename, [] ).append( str )
                if len(chains["psipred"]) > 0:
                    for ssename, str in manager.get_expected_psipred( chains ):
                        data.setdefault( ssename, [] ).append( str )
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
                if labinfo[0].upper() in data:
                    data[labinfo[0].upper()][-1] = labinfo[1]
            continue

    df = cp.DesignFrame( data )
    df.add_source_files( _gather_file_list( filename, multi ) )
    return df

def parse_rosetta_fragments( filename ):
    """
    Read a Rosetta fragment-file and return the appropiate :py:class:`.FragmentFrame`.
    It supports both old and new fragment formats. Does not support varying size fragment sets.

    :param filename: File containing the Rosetta fragments.
    :type filename: :py:class:`str`
    :return: :py:class:`.FragmentFrame`.

    :raises:
        :IOError: if ``filename`` cannot be found.
    """
    fformat = 1 # formats are identified as 1 and 0
    data = OrderedDict({"frame":[], "neighbors":[], "neighbor":[], "position":[], "size":[],
                        "aa":[], "sse":[], "phi":[], "psi":[], "omega":[]})

    if not os.path.isfile(filename):
        raise IOError("{} not found!".format(filename))

    fframe, fne, fpos, fsize, faa, fsse, fphi, fpsi, fomega = None, None, None, None, None, None, None, None, None
    fsaved_size = 0
    nei = 0
    was_space = False
    fd = gzip.open( filename ) if filename.endswith(".gz") else open( filename )
    for line in fd:
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
    data["size"] = [fsaved_size,] * len(data["frame"])

    df = cp.FragmentFrame(data, file=filename)
    if bool(fformat):
        df2 = df.groupby(["frame", "size"]).size().reset_index(name="neighbors")
        df2["neighbors"] = (df2["neighbors"]/df2["size"]).astype(int)
        df = df.merge(df2, on=["frame", "size"])
        df = df.drop(["neighbors_x"], axis=1)
        df = df.rename({"neighbors_y": "neighbors"}, axis=1)
    return df.reindex(["frame", "neighbors", "neighbor","position", "size", "aa", "sse", "phi", "psi", "omega"], axis=1)

def get_sequence_and_structure( pdbfile ):
    """
    Provided a PDB file, it will run a small RosettaScript to capture the sequence and structure that PDB.

    It will generate an output file called: ``[pdbfile].dssp.minisilent``. If this file already exists, it will
    be directly read. You can choose to compress it; ``[pdbfile].dssp.minisilent.gz`` will also work.

    :param pdbfile: Name of the input structure.
    :type pdbfile: :py:class:`str`

    :return: :py:class:`.DesignFrame`.

    :raises:
        :IOError: if ``pdbfile`` cannot be found.
        :IOError: if Rosetta executable cannot be found.
        :ValueError: if Rosetta execution fails
    """
    if not os.path.isfile( pdbfile ):
        raise IOError("Structure {} cannot be found".format(pdbfile))
    minisilent = re.sub("\.pdb|\.cif$", "", re.sub("\.gz$", "", pdbfile)) + ".dssp.minisilent"
    if os.path.isfile(minisilent):
        return parse_rosetta_file(minisilent, {"sequence": "*", "structure": "*"})
    elif os.path.isfile(minisilent + ".gz"):
        return parse_rosetta_file(minisilent + ".gz", {"sequence": "*", "structure": "*"})

    script = "<ROSETTASCRIPTS><MOVERS><WriteSSEMover dssp=\"1\" name=\"w\" /></MOVERS><PROTOCOLS><Add mover=\"w\" /></PROTOCOLS></ROSETTASCRIPTS>"
    with open("dssp.xml", "w") as fd:
        fd.write(script)

    # Check rosetta executable & run
    exe = os.path.join( core.get_option("rosetta", "path"), "rosetta_scripts.{0}".format(core.get_option("rosetta", "compilation")))
    if not os.path.isfile(exe):
        raise IOError("The expected Rosetta executable {0} is not found".format(exe))
    command = "{0} -parser:protocol {1} -s {2} -out:file:silent {3} -ignore_unrecognized_res".format( exe, "dssp.xml", pdbfile, str(os.getpid()) + "_" )
    sys.stdout.write("Running Rosetta\n")
    sys.stdout.write(command + "\n")
    error = os.system( command )
    os.unlink("dssp.xml")
    if not bool(error):
        if os.path.isfile(str(os.getpid()) + "_"):
            sys.stdout.write("Execution has finished\n")
            fd = open( minisilent, "w" )
            for line, is_header, count, symm in open_rosetta_file( str(os.getpid()) + "_" ):
                fd.write( line )
            fd.close()
            os.unlink(str(os.getpid()) + "_")
            return get_sequence_and_structure( pdbfile )
        else:
            raise ValueError("Execution has failed\n")
    else:
        raise ValueError("Execution has failed\n")


def make_structures( df, outdir=None, tagsfilename="tags", prefix=None, keep_tagfile=True ):
    """
    Extract the selected decoys (if any).
    It takes several assumptions:

    #. There is a local instalation of Rosetta.
    #. The global options *rosetta.path* and *rosetta.compilation* are correctly set up.
    #. The provided data is contained in a :py:class:`.DesignFrame`. This should be direct if using this \
    library, otherwise it is easy to cast from a :py:class:`~pandas.DataFrame`::

        df = rstoolbox.components.DesignFrame(df).rename(columns={"identifier_col_name": "description"})

    #. The :py:class:`.DesignFrame` has *silent files* attached from which to extract the data. Again, this \
    should happen by default with the library, but can be set up with::

        # (1) Read from a minisilent file that does not contain structural data: substitute
        df.replace_source_files(["file1", "file2", ])
        # (2) Add files to a recently casted DesignFrame
        df.add_source_files(["file1", "file2", ])

    :param df: Selected set of decoy that have to be extracted.
    :type df: :py:class:`.DesignFrame`
    :param outdir: Directory in which to save the PDB files. If none is provided, it will be loaded
        from the *system.ouput* global option.
    :type outdir: :py:class:`str`
    :param tagsfilename: Name of the file containing the ids of the decoys of interest. It will be created
        in the ``outdir``. An previously existing file will not be overwritten if the global option
        *system.overwrite* is :py:data:`False`.
    :type tagsfilename: :py:class:`str`
    :param prefix: If provided, a prefix is added to the PDB files.
    :type prefix: :py:class:`str`
    :param keep_tagfile: If :py:data:`True`, do not delete the tag file after using it.
    :type keep_tagfile: :py:class:`bool`

    :raises:
        :ValueError: if the provided data does not have a **description** column.
        :ValueError: if the **description** column has repeated identifiers.
        :AttributeError: if silent files from where to extract structures are not found.
        :IOError: if the attached silent files do not exist.
        :IOError: when trying to overwrite the ``tagsfilename`` if *system.overwrite* is :py:data:`False`.
        :IOError: if the rosetta executable is not found. Depends on *rosetta.path* and *rosetta.compilation*.

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
        raise ValueError("There are repeated identifiers. This might indicate the merging of files with identical prefixes "
        "and will be an issue with extracting the structures.")

    # Check that we have associated silent files to extract the data from
    if not isinstance(df, cp.DesignFrame) or len(df.get_source_files()) == 0:
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
    exe = os.path.join( core.get_option("rosetta", "path"), "extract_pdbs.{0}".format(core.get_option("rosetta", "compilation")))
    if not os.path.isfile(exe):
        raise IOError("The expected Rosetta executable {0} is not found".format(exe))

    # Print the tag file
    df[[column]].to_csv( tagsfilename, index=False, header=False)
    if not os.path.isfile(tagsfilename):
        raise IOError("Something went wrong writing the file {0}".format(tagsfilename))

    # Run process
    sfiles = " ".join(sfiles)
    command = "{0} -in:file:silent {1} -in:file:tagfile {2} -out:prefix {3}".format( exe, sfiles, tagsfilename, outdir )
    sys.stdout.write("Executing Rosetta's extract_pdbs app\n")
    sys.stdout.write("(depending on the total number of decoys and how many have been requested this might take a while...)\n")
    error = os.system( command )
    if not bool(error):
        sys.stdout.write("Execution has finished\n")
    else:
        sys.stdout.write("Execution has failed\n")

    # Remove extra files if requested
    if not keep_tagfile:
        os.unlink( tagsfilename )
