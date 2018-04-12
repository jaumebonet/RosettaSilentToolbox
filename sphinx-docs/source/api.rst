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

IO: Rosetta
-----------

Helper functions to read/write data generated with `Rosetta <https://www.rosettacommons.org/>`_. They can be called through ``rstoolbox.io``.

.. autosummary::
   :toctree: generated/

   ~io.parse_rosetta_file
   ~io.parse_rosetta_contacts
   ~io.parse_rosetta_fragments
   ~io.write_rosetta_fragments
   ~io.write_fragment_sequence_profiles
   ~io.get_sequence_and_structure
   ~io.make_structures

Analysis
--------

Helper functions for sequence analysis. They can be called through ``rstoolbox.analysis``.

.. autosummary::
   :toctree: generated/

   ~analysis.sequential_frequencies
   ~analysis.sequence_similarity
   ~analysis.positional_sequence_similarity
   ~analysis.cumulative
   ~analysis.positional_structural_count
   ~analysis.positional_structural_identity

Plot
----

Once the data is loaded in the different **components**, it is ready to use into any
plotting library, but some special plotting alternatives are offered through ``rstoolbox.plot``.

.. autosummary::
   :toctree: generated/

   ~plot.multiple_distributions
   ~plot.sequence_frequency_plot
   ~plot.logo_plot
   ~plot.positional_sequence_similarity_plot
   ~plot.plot_fragments
   ~plot.plot_fragment_profiles

Utils: Plot
-----------

Special functions to help personalise your plot easily can be loaded through ``rstoolbox.utils``.

.. autosummary::
   :toctree: generated/

   ~utils.format_Ipython
   ~utils.add_right_title
   ~utils.add_top_title
   ~utils.add_white_to_cmap
   ~utils.color_variant

Utils: Transforms
-----------------

Special functions to help transform your data can be loaded through ``rstoolbox.utils``.

.. autosummary::
   :toctree: generated/

   ~utils.add_column
   ~utils.split_values

Utils: RosettaScript
--------------------

Get the **RosettaScripts** that are called by different functions of the library with ``rstoolbox.utils``.

.. autosummary::
   :toctree: generated/

   ~utils.baseline

Internals: Functions for Developers
-----------------------------------

This functions are only of interest if you plan on writing new functionalities in ``rstoolbox``.

.. autosummary::
   :toctree: generated/

   io.open_rosetta_file
   components.get_selection
