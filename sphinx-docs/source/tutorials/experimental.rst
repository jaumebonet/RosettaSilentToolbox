.. _experimental:

Experimental Data
=================

Experimental data derived from the analysis of the design decoys can be added and processed through the library. For some experiments, like
**Circular Dichroism** or **Surface Plasmon Resonance**, very specific formats exist. In those cases, the proper parser (:func:`.read_CD` or
:func:`.read_SPR`) exists. Other cases provide different +*CSV* files. This can be imported to python with :func:`~pandas.read_csv, but might
need to be post-processed to work with them.

Each plotting function shows in its documentation the exact naming of the column that it expects.


Circular Dichroism
------------------



Thermal Melt
------------



Multi-Angle Light Scattering
----------------------------


Surface Plasmon Resonance
-------------------------

