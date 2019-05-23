.. _api_ref:

.. currentmodule:: rstoolbox

API reference
=============

.. _grid_api:

Components
----------

These are the list of dedicated objects provided to manage design data. They can be called through ``rstoolbox.components``.

.. autosummary::
   :toctree: generated/
   :template: class_summary.rst

   ~components.Selection
   ~components.SelectionContainer
   ~components.DesignSeries
   ~components.DesignFrame
   ~components.SequenceFrame
   ~components.FragmentFrame


IO: Sequence
------------

Helper functions to read/write direct sequence information. They can be called through ``rstoolbox.io``.

.. autosummary::
   :toctree: generated/

   ~io.read_fasta
   ~io.write_fasta
   ~io.write_clustalw
   ~io.write_mutant_alignments
   ~io.read_hmmsearch
   ~io.pymol_mutant_selector


IO: Structure
-------------

Helper functions to read/write outputs of programs based on protein structure. They can be called through ``rstoolbox.io``.

.. autosummary::
   :toctree: generated/

   ~io.parse_master_file


IO: Rosetta
-----------

Helper functions to read/write data generated with `Rosetta <https://www.rosettacommons.org/>`_. They can be called through ``rstoolbox.io``.

.. autosummary::
   :toctree: generated/

   ~io.parse_rosetta_file
   ~io.parse_rosetta_json
   ~io.parse_rosetta_pdb
   ~io.parse_rosetta_contacts
   ~io.parse_rosetta_fragments
   ~io.write_rosetta_fragments
   ~io.write_fragment_sequence_profiles
   ~io.get_sequence_and_structure
   ~io.make_structures

IO: Experiments
---------------

Helper functions to read/write data generated through wedlab experiments. They can be called through ``rstoolbox.io``.

.. autosummary::
   :toctree: generated/

   ~io.read_SPR
   ~io.read_CD
   ~io.read_MALS
   ~io.read_fastq

Analysis
--------

Helper functions for sequence analysis. They can be called through ``rstoolbox.analysis``.

.. autosummary::
   :toctree: generated/

   ~analysis.sequential_frequencies
   ~analysis.sequence_similarity
   ~analysis.positional_sequence_similarity
   ~analysis.binary_similarity
   ~analysis.binary_overlap
   ~analysis.positional_enrichment
   ~analysis.positional_structural_count
   ~analysis.positional_structural_identity
   ~analysis.secondary_structure_percentage
   ~analysis.selector_percentage
   ~analysis.label_percentage
   ~analysis.label_sequence
   ~analysis.cumulative

Plot
----

Once the data is loaded in the different **components**, it is ready to use into any
plotting library, but some special plotting alternatives are offered through ``rstoolbox.plot``.

.. autosummary::
   :toctree: generated/

   ~plot.multiple_distributions
   ~plot.sequence_frequency_plot
   ~plot.logo_plot
   ~plot.logo_plot_in_axis
   ~plot.positional_sequence_similarity_plot
   ~plot.per_residue_matrix_score_plot
   ~plot.positional_structural_similarity_plot
   ~plot.plot_fragments
   ~plot.plot_fragment_profiles
   ~plot.plot_alignment
   ~plot.plot_ramachandran
   ~plot.plot_ramachandran_single
   ~plot.plot_dssp_vs_psipred

Plot: Experiments
-----------------

Plot data obtained from experimental procedures. Accessible through ``rstoolbox.plot``.

.. autosummary::
   :toctree: generated/

   ~plot.plot_96wells
   ~plot.plot_thermal_melt
   ~plot.plot_MALS
   ~plot.plot_CD
   ~plot.plot_SPR

Utils: Plot
-----------

Special functions to help personalise your plot easily can be loaded through ``rstoolbox.utils``.

.. autosummary::
   :toctree: generated/

   ~utils.format_Ipython
   ~utils.highlight
   ~utils.use_qgrid
   ~utils.add_left_title
   ~utils.add_right_title
   ~utils.add_top_title
   ~utils.edit_legend_text
   ~utils.add_white_to_cmap
   ~utils.color_variant

Utils: Contextualize
--------------------

Functions aimed to help assess a design population in the context of known protein structures.

.. autosummary::
   :toctree: generated/

   ~utils.load_refdata
   ~utils.make_redundancy_table
   ~plot.plot_in_context
   ~plot.distribution_quality

Utils: Transforms
-----------------

Special functions to help transform your data can be loaded through ``rstoolbox.utils``.

.. autosummary::
   :toctree: generated/

   ~utils.add_column
   ~utils.split_values
   ~utils.split_dataframe_rows
   ~utils.report
   ~utils.concat_fragments

Utils: RosettaScript
--------------------

Get the **RosettaScripts** that are called by different functions of the library with ``rstoolbox.utils``.

.. autosummary::
   :toctree: generated/

   ~utils.baseline
   ~utils.mutations

Utils: Experiments
------------------

Special functions to help obtain data from multiple **Next Generation Sequencing** data.Can be loaded through ``rstoolbox.utils``.

.. autosummary::
   :toctree: generated/

   ~utils.translate_dna_sequence
   ~utils.translate_3frames
   ~utils.adapt_length
   ~utils.sequencing_enrichment

Internals: Functions for Developers
-----------------------------------

This functions are only of interest if you plan on writing new functionalities in ``rstoolbox``.

.. autosummary::
   :toctree: generated/

   io.open_rosetta_file
   components.get_selection
   utils.make_rosetta_app_path
   tests.helper.random_frequency_matrix
   tests.helper.random_proteins
   tests.helper.random_fastq
