.. _experimental_analysis:

.. currentmodule:: rstoolbox

Experimental Data
=================

Experimental data derived from the analysis of the design decoys can be added and processed through the library. For some experiments, like
**Circular Dichroism** or **Surface Plasmon Resonance**, very specific formats exist. In those cases, the proper parser (:func:`.read_CD` or
:func:`.read_SPR`) exists. Other cases provide different +*CSV* files. This can be imported to python with :func:`~pandas.read_csv`, but might
need to be post-processed to work with them.

Each plotting function shows in its documentation the exact naming of the column that it expects.


Circular Dichroism
------------------

Circular dichroism data can be easily worked with.

Currently the library can automatically import data compatible with the `Jasco <https://jascoinc.com/products/spectroscopy/circular-dichroism/>`_ J-815
series by means of the :func:`.read_CD` function, but any data can be read as long as it finally produces a :class:`~pandas.DataFrame` with, at least, the following
two columns:

===============  ===================================================
Column Name       Data Content
===============  ===================================================
**Wavelength**   Wavelength (nm).
**MRE**          Value at each wavelength (10 deg^2 cm^2 dmol^-1).
===============  ===================================================

:ref:`Adding a new available format to the library <contributing>` would be as easy as adding the appropriate private function to ``rstoolbox/io/experimental`` and provide access to
it through the :func:`.read_CD` function.

.. note::
  When loading data with :func:`.read_CD` notice that the :class:`~pandas.DataFrame` contains extra columns with extra information. From those, ``Temp`` is the more relevant for
  plotting and analysis, as is the one that allows to separate the different data series.

Currently, CD data can be obtained and plotted such as:

.. ipython::
  :okwarning:

  In [1]: import rstoolbox as rs
     ...: dfCD = rs.io.read_CD("../rstoolbox/tests/data/CD", model='J-815')
     ...: dfCD.columns

  In [2]: import matplotlib.pyplot as plt
     ...: fig = plt.figure(figsize=(10, 6.7))
     ...: ax = plt.subplot2grid((1, 1), (0, 0), fig=fig)
     ...: rs.plot.plot_CD(dfCD, ax, sample=6)

  @savefig tutorial_exp_plt1.png width=5in
  In [3]: plt.show()

.. tip::
  Depending on how many temperatures one might have, trying to plot all of them will create a very large legend. One can either delete it or use the ``sample`` option, which
  uses the `Bresenham's line algorithm <http://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm>`_ to sample temperatures across all possible

Thermal Melt
------------

As before, thermal melt data can be loaded by properly parsing any input file as long as the generated :class:`~pandas.DataFrame` contains at least the following columns:

===============  ===================================================
Column Name       Data Content
===============  ===================================================
**Temp**         Temperatures (celsius).
**MRE**          Value at each temperature (10 deg^2 cm^2 dmol^-1).
===============  ===================================================

If loaded from CD data containing a temp column (as is the case when loaded with :func:`.read_CD`), one can directly pick the data from there:

.. ipython::
  :okwarning:

  In [4]: dfTM = dfCD[(dfCD['Wavelength'] == 220)]
     ...: fig = plt.figure(figsize=(10, 6.7))
     ...: ax = plt.subplot2grid((1, 1), (0, 0), fig=fig)
     ...: rs.plot.plot_thermal_melt(dfTM, ax)

  @savefig tutorial_exp_plt2.png width=5in
  In [5]: plt.show()

For well-behaved proteins, the library can approximate the melting point.

.. ipython::
  :okwarning:

  In [6]: fig = plt.figure(figsize=(10, 6.7))
     ...: ax = plt.subplot2grid((1, 1), (0, 0), fig=fig)
     ...: rs.plot.plot_thermal_melt(dfTM, ax, fusion_temperature=True, temp_marker=True)

  @savefig tutorial_exp_plt3.png width=5in
  In [7]: plt.show()


Multi-Angle Light Scattering
----------------------------

MALS data can be manually loaded or read through the :func:`.read_MALS` function.
In order to properly work with the plotting functions, it needs to have at least two of the following columns:

===============  ===================================================
Column Name       Data Content
===============  ===================================================
**Time**         Time (min).
**UV**           UV data (V).
**LS**           Light Scattering data (V).
**MW**           Molecular Weight (Daltons).
===============  ===================================================

being **Time** always mandatory. Depending on the amount of data, the plot function can pick which information to print.

.. ipython::
  :okwarning:

  In [8]: import pandas as pd
     ...: dfMALS = pd.read_csv("../rstoolbox/tests/data/mals.csv")
     ...: fig = plt.figure(figsize=(10, 6.7))
     ...: ax = plt.subplot2grid((1, 1), (0, 0))
     ...: rs.plot.plot_MALS(dfMALS, ax)

  @savefig tutorial_exp_plt4.png width=5in
  In [9]: plt.show()


Surface Plasmon Resonance
-------------------------

Optimal SPR data has to contain the fitted curves on input, but can be read directly from the default output of the machine.
More details on the format are listed in the corresponding function :func:`.read_SPR`.

The fitted curves will be directly represented with :func:`.plot_SPR`.


.. ipython::
  :okwarning:

  In [1]: dfSPR = rs.io.read_SPR("../rstoolbox/tests/data/spr_data.csv.gz")
     ...: fig = plt.figure(figsize=(10, 6.7))
     ...: ax = plt.subplot2grid((1, 1), (0, 0))
     ...: rs.plot.plot_SPR(dfSPR, ax, datacolor='black', fitcolor='red')

  @savefig tutorial_exp_plt5.png width=5in
  In [2]: plt.show()

  In [3]: plt.close('all')

Deep Sequence Analysis and Population Enrichment
------------------------------------------------

Individual raw **FASTQ** data obtained from deep sequencing can be directly read with the :func:`.read_fastq` function.
In order to evaluate the enrichment of different outputs submitted to a different experimental conditions in different
concentrations, and provide the ``min`` and ``max`` concentrations to use for the calculus to :func:`.sequence_enrichment`.

.. ipython::
  :okwarning:

  In [4]: indat = {'binder1': {'conc1': '../rstoolbox/tests/data/cdk2_rand_001.fasq.gz',
     ...:                      'conc2': '../rstoolbox/tests/data/cdk2_rand_002.fasq.gz',
     ...:                      'conc3': '../rstoolbox/tests/data/cdk2_rand_003.fasq.gz'},
     ...:          'binder2': {'conc1': '../rstoolbox/tests/data/cdk2_rand_004.fasq.gz',
     ...:                      'conc2': '../rstoolbox/tests/data/cdk2_rand_005.fasq.gz',
     ...:                      'conc3': '../rstoolbox/tests/data/cdk2_rand_006.fasq.gz'}}
     ...: enrich = {'binder1': ['conc1', 'conc3'],
     ...:           'binder2': ['conc1', 'conc3']}
     ...: dfSeq = rs.utils.sequencing_enrichment(indat, enrich)
     ...: dfSeq[[_ for _ in dfSeq.columns if _ != 'sequence_A']].head()