rstoolbox.components.DesignFrame
================================

.. currentmodule:: rstoolbox.components

.. autoclass:: DesignFrame

  .. rubric:: Getters
  .. autosummary::
     :toctree: DesignFrame

      ~DesignFrame.get_id
      ~DesignFrame.get_available_sequences
      ~DesignFrame.get_sequence
      ~DesignFrame.get_available_structures
      ~DesignFrame.get_structure
      ~DesignFrame.get_available_structure_predictions
      ~DesignFrame.get_structure_prediction
      ~DesignFrame.get_sequential_data
      ~DesignFrame.get_dihedrals
      ~DesignFrame.get_phi
      ~DesignFrame.get_psi
      ~DesignFrame.get_available_labels
      ~DesignFrame.get_label


  .. rubric:: Reference Data
  .. autosummary::
     :toctree: DesignFrame

      ~DesignFrame.has_reference_sequence
      ~DesignFrame.add_reference_sequence
      ~DesignFrame.get_reference_sequence
      ~DesignFrame.has_reference_structure
      ~DesignFrame.add_reference_structure
      ~DesignFrame.get_reference_structure
      ~DesignFrame.add_reference_shift
      ~DesignFrame.get_reference_shift
      ~DesignFrame.get_available_references
      ~DesignFrame.add_reference
      ~DesignFrame.transfer_reference
      ~DesignFrame.delete_reference


  .. rubric:: Source Files
  .. autosummary::
     :toctree: DesignFrame

      ~DesignFrame.add_source_file
      ~DesignFrame.add_source_files
      ~DesignFrame.get_source_files
      ~DesignFrame.has_source_files
      ~DesignFrame.replace_source_files

  .. rubric:: Frequencies
  .. autosummary::
     :toctree: DesignFrame

     ~DesignFrame.sequence_bits
     ~DesignFrame.sequence_distance
     ~DesignFrame.sequence_frequencies
     ~DesignFrame.structure_bits
     ~DesignFrame.structure_frequencies


  .. rubric:: Mutation Methods
  .. autosummary::
     :toctree: DesignFrame

      ~DesignFrame.identify_mutants
      ~DesignFrame.get_identified_mutants
      ~DesignFrame.get_mutation_count
      ~DesignFrame.get_mutation_positions
      ~DesignFrame.get_mutations
      ~DesignFrame.get_sequence_with
      ~DesignFrame.generate_mutant_variants
      ~DesignFrame.generate_mutants_from_matrix
      ~DesignFrame.generate_wt_reversions
      ~DesignFrame.make_resfile
      ~DesignFrame.apply_resfile
      ~DesignFrame.score_by_pssm
      ~DesignFrame.view_mutants_alignment

  .. rubric:: Miscellaneous
  .. autosummary::
     :toctree: DesignFrame

      ~DesignFrame.clean_rosetta_suffix
      ~DesignFrame.retrieve_sequences_from_pdbs