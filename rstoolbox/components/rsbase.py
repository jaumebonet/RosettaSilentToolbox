# @Author: Jaume Bonet <bonet>
# @Date:   22-Feb-2018
# @Email:  jaume.bonet@gmail.com
# @Filename: rsbase.py
# @Last modified by:   bonet
# @Last modified time: 14-Mar-2018


import sys

import rstoolbox.utils as ru


class RSBaseDesign( object ):
    pass


class RSBaseFrequency( object ):
    pass


if (sys.version_info > (3, 0)):
    RSBaseDesign.get_id                              = ru.get_id
    RSBaseDesign.get_sequence                        = ru.get_sequence
    RSBaseDesign.get_available_sequences             = ru.get_available_sequences
    RSBaseDesign.get_structure                       = ru.get_structure
    RSBaseDesign.get_available_structures            = ru.get_available_structures
    RSBaseDesign.get_structure_prediction            = ru.get_structure_prediction
    RSBaseDesign.get_available_structure_predictions = ru.get_available_structure_predictions
    RSBaseDesign.get_sequential_data                 = ru.get_sequential_data
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
    RSBaseDesign.delete_reference                    = ru.delete_reference

    RSBaseDesign.get_identified_mutants              = ru.get_identified_mutants
    RSBaseDesign.get_mutations                       = ru.get_mutations
    RSBaseDesign.get_mutation_positions              = ru.get_mutation_positions
    RSBaseDesign.get_mutation_count                  = ru.get_mutation_count
    RSBaseDesign.identify_mutants                    = ru.identify_mutants
    RSBaseDesign.generate_mutant_variants            = ru.generate_mutant_variants

    RSBaseFrequency.has_reference_sequence           = ru.has_reference_sequence
    RSBaseFrequency.add_reference_sequence           = ru.add_reference_sequence
    RSBaseFrequency.get_reference_sequence           = ru.get_reference_sequence
    RSBaseFrequency.has_reference_structure          = ru.has_reference_structure
    RSBaseFrequency.add_reference_structure          = ru.add_reference_structure
    RSBaseFrequency.get_reference_structure          = ru.get_reference_structure
    RSBaseFrequency.add_reference_shift              = ru.add_reference_shift
    RSBaseFrequency.get_reference_shift              = ru.get_reference_shift
    RSBaseFrequency.add_reference                    = ru.add_reference
    RSBaseFrequency.transfer_reference               = ru.transfer_reference
    RSBaseFrequency.delete_reference                    = ru.delete_reference

else:
    from types import MethodType
    RSBaseDesign.get_id = MethodType(
        ru.get_id, None, RSBaseDesign)
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
    RSBaseDesign.get_sequential_data = MethodType(
        ru.get_sequential_data, None, RSBaseDesign)
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
    RSBaseDesign.delete_reference = MethodType(
        ru.delete_reference, None, RSBaseDesign)

    RSBaseDesign.get_identified_mutants = MethodType(
        ru.get_identified_mutants, None, RSBaseDesign)
    RSBaseDesign.get_mutations = MethodType(
        ru.get_mutations, None, RSBaseDesign)
    RSBaseDesign.get_mutation_positions = MethodType(
        ru.get_mutation_positions, None, RSBaseDesign)
    RSBaseDesign.get_mutation_count = MethodType(
        ru.get_mutation_count, None, RSBaseDesign)
    RSBaseDesign.identify_mutants = MethodType(
        ru.identify_mutants, None, RSBaseDesign)
    RSBaseDesign.generate_mutant_variants = MethodType(
        ru.generate_mutant_variants, None, RSBaseDesign)

    RSBaseFrequency.has_reference_sequence = MethodType(
        ru.has_reference_sequence, None, RSBaseFrequency)
    RSBaseFrequency.add_reference_sequence = MethodType(
        ru.add_reference_sequence, None, RSBaseFrequency)
    RSBaseFrequency.get_reference_sequence = MethodType(
        ru.get_reference_sequence, None, RSBaseFrequency)
    RSBaseFrequency.has_reference_structure = MethodType(
        ru.has_reference_structure, None, RSBaseFrequency)
    RSBaseFrequency.add_reference_structure = MethodType(
        ru.add_reference_structure, None, RSBaseFrequency)
    RSBaseFrequency.get_reference_structure = MethodType(
        ru.get_reference_structure, None, RSBaseFrequency)
    RSBaseFrequency.add_reference_shift = MethodType(
        ru.add_reference_shift, None, RSBaseFrequency)
    RSBaseFrequency.get_reference_shift = MethodType(
        ru.get_reference_shift, None, RSBaseFrequency)
    RSBaseFrequency.add_reference = MethodType(
        ru.add_reference, None, RSBaseFrequency)
    RSBaseFrequency.transfer_reference = MethodType(
        ru.transfer_reference, None, RSBaseFrequency)
    RSBaseFrequency.delete_reference = MethodType(
        ru.delete_reference, None, RSBaseFrequency)
