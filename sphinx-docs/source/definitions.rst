.. _definitions:

.. |df_param| replace:: Data container

.. |seqID_param| replace:: Identifier of the sequence of interest

.. |seqType_param| replace:: Type of sequence: ``protein``, ``dna``, ``rna``

.. |query_param| replace:: Content type to load from the input data ``sequence``, ``structure``, ``structure_prediction``

.. |cleanExtra_param| replace:: Remove from the :class:`.SequenceFrame` the non-regular amino/nucleic acids if they are empty for all positions;
  basically remove `ambiguous` and `gap` identifiers

.. |cleanUnused_param| replace:: Remove from the :class:`.SequenceFrame` the regular amino/nucleic acids if they frequency is equal or under the value;
  basically this targets fully empty positions to minimise the size of the matrix. The value itself represents the threshold to consider a position empty.
  Thus, ``-1`` triggers no filter while ``0.3`` would consider all the frequencies equal or lower than that value as empty

.. |matrix_param| replace:: Identifier of the matrix used to evaluate similarity

.. |keyres_param| replace:: Residues of interest
.. |keyres_types| replace:: Union[:class:`int`, :func:`list` of :class:`int`, :class:`str`, :class:`.Selection`]

.. |axis_param| replace:: ``matplotlib`` axis to which we will plot

.. |color_types| replace:: Union[:class:`int`, :class:`str`]

.. |designframe_cast_error| replace:: if the data passed is not a :class:`.DesignFrame` or cannot be casted to one

.. |seqID_error| replace:: if there is no sequence information for chain ``seqID`` of the decoys
.. |sseID_error| replace:: if there is no structure information for chain ``seqID`` of the decoys

.. |reference_error| replace:: if there is no ``reference_sequence`` for chain ``seqID`` of the decoys