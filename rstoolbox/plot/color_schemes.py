# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: color_scheme
"""
# Standard Libraries

# External Libraries

# This Library


def color_scheme( name ):
    if name.upper() == "WEBLOGO":
        return {
            'A': '#CCFF00', 'C': '#FFFF00', 'D': '#FF0000', 'E': '#FF0066',
            'F': '#00FF66', 'G': '#FF9900', 'H': '#0066FF', 'I': '#66FF00',
            'K': '#6600FF', 'L': '#33FF00', 'M': '#00FF00', 'N': '#CC00FF',
            'P': '#FFCC00', 'Q': '#FF00CC', 'R': '#0000FF', 'S': '#FF3300',
            'T': '#FF6600', 'V': '#99FF00', 'W': '#00CCFF', 'Y': '#00FFCC'
        }
    # Color schemes resembling the default schemes available through WebLogo
    elif name.upper() == "HYDROPHOBICITY":
        return {
            'A': '#058100', 'C': '#000000', 'D': '#3603FF', 'E': '#3603FF',
            'F': '#000000', 'G': '#058100', 'H': '#058100', 'I': '#000000',
            'K': '#3603FF', 'L': '#000000', 'M': '#000000', 'N': '#3603FF',
            'P': '#058100', 'Q': '#3603FF', 'R': '#3603FF', 'S': '#058100',
            'T': '#058100', 'V': '#000000', 'W': '#000000', 'Y': '#000000'
        }
    elif name.upper() == "CHEMISTRY":
        return {
            'A': '#000000', 'C': '#058100', 'D': '#FA2202', 'E': '#FA2202',
            'F': '#000000', 'G': '#058100', 'H': '#3603FF', 'I': '#000000',
            'K': '#3603FF', 'L': '#000000', 'M': '#000000', 'N': '#970C83',
            'P': '#000000', 'Q': '#970C83', 'R': '#3603FF', 'S': '#058100',
            'T': '#058100', 'V': '#000000', 'W': '#000000', 'Y': '#058100'
        }
    elif name.upper() == "CHARGE":
        return {
            'A': '#000000', 'C': '#000000', 'D': '#FA2202', 'E': '#FA2202',
            'F': '#000000', 'G': '#000000', 'H': '#3603FF', 'I': '#000000',
            'K': '#3603FF', 'L': '#000000', 'M': '#000000', 'N': '#000000',
            'P': '#000000', 'Q': '#000000', 'R': '#3603FF', 'S': '#000000',
            'T': '#000000', 'V': '#000000', 'W': '#000000', 'Y': '#000000'
        }
    else:
        raise ValueError("Color scheme {} not found".format(name))
