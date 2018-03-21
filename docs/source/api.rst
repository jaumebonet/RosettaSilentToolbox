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
   ~components.DesignFrame


IO: Sequence
------------

Helper functions to read/write direct sequence information. They can be called through ``rstoolbox.io``.

.. autosummary::
   :toctree: generated/

   ~io.read_fasta
   ~io.write_fasta

IO: Rosetta
-----------

Helper functions to read/write data generated with `Rosetta <https://www.rosettacommons.org/>`_. They can be called through ``rstoolbox.io``.

.. autosummary::
   :toctree: generated/

   ~io.open_rosetta_file
   ~io.parse_rosetta_file
   ~io.parse_rosetta_contacts
   ~io.parse_rosetta_fragments
   ~io.write_rosetta_fragments
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

   ~plot.sequence_frequency_plot
   ~plot.positional_sequence_similarity_plot

Utils: Plot
-----------

Special functions to help personalise your plot easily can be loaded through ``rstoolbox.utils``.

.. autosummary::
   :toctree: generated/

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
