# @Author: Jaume Bonet <bonet>
# @Date:   29-Jun-2017
# @Email:  jaume.bonet@gmail.com
# @Filename: api.py
# @Last modified by:   bonet
# @Last modified time: 29-Jun-2017

import os

from .io import read_silent_file, read_score_file
from .components import *

import yaml
import pandas as pd

def read_rosetta_silent(filename, design_type="decoy", filterfile=None):
    dlist = read_silent_file( filename, design_type )
    if filterfile is not None:
        return dlist.filter_by_tags_file(filterfile)
    return dlist

def process_from_yaml_file( dlist, yamlfile ):
    assert os.path.isfile( yamlfile )
    defsdic = yaml.load("".join(open(yamlfile).readlines()))
    return process_from_definitions( dlist, defsdic )

def process_from_definitions( dlist, defsdic ):
    assert len(dlist) > 0
    assert "scores" in defsdic
    assert len(defsdic["scores"]) > 0
    if "description" not in defsdic["scores"]:
        defsdic["scores"]["description"] = "description"
    data = {}
    for d in dlist:
        for score in defsdic["scores"]:
            name = defsdic["scores"][score]
            data.setdefault(name, []).append(d[score])

        if "sequence" in defsdic:
            seqinf = defsdic["sequence"]
            chains = seqinf["chains"] if "chains" in seqinf else None
            data.setdefault("sequence", []).append(d.get_sequence(chains=chains))
    df = pd.DataFrame( data )

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
