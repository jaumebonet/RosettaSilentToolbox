# @Author: Jaume Bonet <bonet>
# @Date:   22-Feb-2018
# @Email:  jaume.bonet@gmail.com
# @Filename: rsbase.py
# @Last modified by:   bonet
# @Last modified time: 05-Mar-2018


import sys

import rstoolbox.utils as ru


class RSBaseDesign(object):
    pass


if (sys.version_info > (3, 0)):
    RSBaseDesign.get_sequence                        = ru.get_sequence
    RSBaseDesign.get_available_sequences             = ru.get_available_sequences
    RSBaseDesign.get_structure                       = ru.get_structure
    RSBaseDesign.get_available_structures            = ru.get_available_structures
    RSBaseDesign.get_structure_prediction            = ru.get_structure_prediction
    RSBaseDesign.get_available_structure_predictions = ru.get_available_structure_predictions
    RSBaseDesign.get_label                           = ru.get_label
    RSBaseDesign.get_available_labels                = ru.get_available_labels

    RSBaseDesign.has_reference_sequence              = ru.has_reference_sequence
    RSBaseDesign.add_reference_sequence              = ru.add_reference_sequence
    RSBaseDesign.get_reference_sequence              = ru.get_reference_sequence
    RSBaseDesign.has_reference_structure             = ru.has_reference_structure
    RSBaseDesign.add_reference_structure             = ru.add_reference_structure
    RSBaseDesign.get_reference_structure             = ru.get_reference_structure
    RSBaseDesign.add_reference_shift                 = ru.add_reference_shift
    RSBaseDesign.get_reference_shift                 = ru.get_reference_shift
    RSBaseDesign.add_reference                       = ru.add_reference
    RSBaseDesign.transfer_reference                  = ru.transfer_reference
else:
    from types import MethodType
    RSBaseDesign.get_sequence = MethodType(
        ru.get_sequence, None, RSBaseDesign)
    RSBaseDesign.get_available_sequences = MethodType(
        ru.get_available_sequences, None, RSBaseDesign)
    RSBaseDesign.get_structure = MethodType(
        ru.get_structure, None, RSBaseDesign)
    RSBaseDesign.get_available_structures = MethodType(
        ru.get_available_structures, None, RSBaseDesign)
    RSBaseDesign.get_structure_prediction = MethodType(
        ru.get_structure_prediction, None, RSBaseDesign)
    RSBaseDesign.get_available_structure_predictions = MethodType(
        ru.get_available_structure_predictions, None, RSBaseDesign)
    RSBaseDesign.get_label = MethodType(
        ru.get_label, None, RSBaseDesign)
    RSBaseDesign.get_available_labels = MethodType(
        ru.get_available_labels, None, RSBaseDesign)

    RSBaseDesign.has_reference_sequence = MethodType(
        ru.has_reference_sequence, None, RSBaseDesign)
    RSBaseDesign.add_reference_sequence = MethodType(
        ru.add_reference_sequence, None, RSBaseDesign)
    RSBaseDesign.get_reference_sequence = MethodType(
        ru.get_reference_sequence, None, RSBaseDesign)
    RSBaseDesign.has_reference_structure = MethodType(
        ru.has_reference_structure, None, RSBaseDesign)
    RSBaseDesign.add_reference_structure = MethodType(
        ru.add_reference_structure, None, RSBaseDesign)
    RSBaseDesign.get_reference_structure = MethodType(
        ru.get_reference_structure, None, RSBaseDesign)
    RSBaseDesign.add_reference_shift = MethodType(
        ru.add_reference_shift, None, RSBaseDesign)
    RSBaseDesign.get_reference_shift = MethodType(
        ru.get_reference_shift, None, RSBaseDesign)
    RSBaseDesign.add_reference = MethodType(
        ru.add_reference, None, RSBaseDesign)
    RSBaseDesign.transfer_reference = MethodType(
        ru.transfer_reference, None, RSBaseDesign)
