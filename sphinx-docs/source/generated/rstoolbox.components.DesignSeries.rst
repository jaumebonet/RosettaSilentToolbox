rstoolbox.components.DesignSeries
=================================

.. currentmodule:: rstoolbox.components

.. autoclass:: DesignSeries

  .. rubric:: Getters
  .. autosummary::
     :toctree: DesignSeries

      ~DesignSeries.get_id
      ~DesignSeries.get_available_sequences
      ~DesignSeries.get_sequence
      ~DesignSeries.get_available_structures
      ~DesignSeries.get_structure
      ~DesignSeries.get_available_structure_predictions
      ~DesignSeries.get_structure_prediction
      ~DesignSeries.get_sequential_data
      ~DesignSeries.get_dihedrals
      ~DesignSeries.get_phi
      ~DesignSeries.get_psi
      ~DesignSeries.get_available_labels
      ~DesignSeries.get_label


  .. rubric:: Reference Data
  .. autosummary::
     :toctree: DesignSeries

      ~DesignSeries.has_reference_sequence
      ~DesignSeries.add_reference_sequence
      ~DesignSeries.get_reference_sequence
      ~DesignSeries.has_reference_structure
      ~DesignSeries.add_reference_structure
      ~DesignSeries.get_reference_structure
      ~DesignSeries.add_reference_shift
      ~DesignSeries.get_reference_shift
      ~DesignSeries.get_available_references
      ~DesignSeries.add_reference
      ~DesignSeries.transfer_reference
      ~DesignSeries.delete_reference


  .. rubric:: Mutation Methods
  .. autosummary::
     :toctree: DesignSeries

      ~DesignSeries.identify_mutants
      ~DesignSeries.get_identified_mutants
      ~DesignSeries.get_mutation_count
      ~DesignSeries.get_mutation_positions
      ~DesignSeries.get_mutations
      ~DesignSeries.generate_mutant_variants
      ~DesignSeries.generate_mutants_from_matrix
      ~DesignSeries.generate_wt_reversions
      ~DesignSeries.make_resfile
      ~DesignSeries.score_by_pssm
      ~DesignSeries.view_mutants_alignment
