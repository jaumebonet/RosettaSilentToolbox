# @Author: Jaume Bonet <bonet>
# @Date:   29-Jun-2017
# @Email:  jaume.bonet@gmail.com
# @Filename: api.py
# @Last modified by:   bonet
# @Last modified time: 21-Sep-2017

import os
import copy

from .io import read_silent_file, read_score_file
from .components import *
from .utils import *

import yaml
import pandas as pd

def read_rosetta_silent(filename, design_type="decoy", filterfile=None, allow_repeats=False):
    dlist = read_silent_file( filename, design_type, allow_repeats )
    if filterfile is not None:
        return dlist.filter_by_tags_file(filterfile)
    return dlist

def process_from_yaml_file( dlist, yamlfile ):
    assert os.path.isfile( yamlfile )
    defsdic = yaml.load("".join(open(yamlfile).readlines()))
    if "sequence" in defsdic:
        defsdic["sequence"]["range"] = selection2list(defsdic["sequence"]["range"]) if "range" in defsdic["sequence"] else None
    return process_from_definitions( dlist, defsdic )

def process_from_definitions( dlist, definition_dict ):
    assert len(dlist) > 0
    assert "scores" in definition_dict
    assert len(definition_dict["scores"]) > 0
    defsdic = copy.deepcopy(definition_dict)
    does_not_want_description = "description" not in defsdic["scores"]
    if does_not_want_description:
        defsdic["scores"]["description"] = "description"
    data = {}
    for d in dlist:
        for score in defsdic["scores"]:
            name = defsdic["scores"][score]
            data.setdefault(name, []).append(d[score])

        if "sequence" in defsdic:
            seqinf = defsdic["sequence"]
            chains = seqinf["chains"] if "chains" in seqinf else None
            selection  = seqinf["range"] if "range" in seqinf else None
            data.setdefault("sequence", []).append(d.get_sequence(chains=chains, selection=selection))
            data.setdefault("ranges", []).append(selection)

        if "naming" in defsdic:
            name = data[defsdic["scores"]["description"]][-1].split("_")
            assert len(defsdic["naming"]) == len(name)
            for cn, n in enumerate(defsdic["naming"]):
                if n != "":
                    data.setdefault(n, []).append(name[cn])

    df = pd.DataFrame( data )
    if does_not_want_description:
        df = df.drop("description", axis=1)

    if "select" in defsdic:
        sele = defsdic["select"]
        assert "value" in sele
        assert "number" in sele
        assert sele["value"] in data
        if "sort" in sele:
            assert sele["sort"] in ["rise", "drop"]
        ascending = False if "sort" in sele and sele["sort"] == "drop" else True
        df = df.sort_values( sele["value"], ascending=ascending ).head(n=sele["number"])


    return df

def split_columns( df, keys ):
    dataframes = []
    for k in keys["split"]:
        colIDs = copy.copy(keys["keep"])
        colIDs.append(k[0])
        wdf = df[colIDs]
        wdf = wdf.assign(temporarykey1=pd.Series([k[1]]*len(wdf[colIDs[0]])).values).copy(True)
        wdf = wdf.rename(index=str, columns={
            k[0]: keys["names"][0],
            "temporarykey1": keys["names"][1]
        })
        if ( len(k) > 2 ):
            wdf = wdf.assign(temporarykey2=pd.Series([k[2]]*len(wdf[colIDs[0]])).values).copy(True)
            wdf = wdf.rename(index=str, columns={
                "temporarykey2": keys["names"][2]
            })
        dataframes.append(wdf)
    return pd.concat(dataframes)
