.. _sequence_analysis:

.. currentmodule:: rstoolbox

Sequence Analysis
=================

Loading a Reference
-------------------

When analysing the outcome of a protein design set, it can be useful to retrieve the data from the source structure (template)
in order to assess, on a sequence level, which changes have happened.

For the purpose of this example, we will use a domain from a `Putative formate dehydrogenase accessory protein <https://www.rcsb.org/structure/2pw9>`_.
To load the data, we will use :func:`.get_sequence_and_structure`. To this function we are going to provide the PDB file ``2pw9C.pdb``, which contains only
the chain C from that crystal structure and it will generate a file named ``2pw9C.dssp.minisilent``. The generated object will contain all the data generated
by Rosetta over the particular structure, *including the sequence*.

.. note::
  Through all the process several times the ``chainID`` of the decoy of interest will be called. This is due to the fact that the library can manipulate
  decoys with multiple chains (designed or not), and, thus, analysis must be called upon the sequences of interest.

.. warning::
  To generate this file, the function **will run Rosetta**. If Rosetta is not locally installed, the documentation of :func:`.get_sequence_and_structure`
  provides the RosettaScript that it will run. One can run it in a different computer/cluster and place the obtained silent file output in the same directory.
  As long as the naming schema for the file is maintained, if the function finds the file it will **skip the Rosetta execution**.

.. ipython::
  :okwarning:

  In [1]: import rstoolbox as rs
     ...: import pandas as pd
     ...: import matplotlib.pyplot as plt
     ...: import seaborn as sns
     ...: plt.rcParams['svg.fonttype'] = 'none' # When plt.savefig to 'svg', text is still text object
     ...: sns.set_style('whitegrid')
     ...: pd.set_option('display.width', 1000)
     ...: pd.set_option('display.max_columns', 500)
     ...: pd.set_option("display.max_seq_items", 3)
     ...: baseline = rs.io.get_sequence_and_structure('../rstoolbox/tests/data/2pw9C.pdb')
     ...: baseline

  In [2]: baseline.get_sequence('C')

Loading the Design Data
-----------------------

The next step would be to load the data that we have obtained from designing. In this case, it is a single file with 567 decoys. This decoys come from a protein
grafting experiment as performed by the structural heterologous grafting protocol `FunFolDes <https://doi.org/10.1101/378976>`_. This basically ensures that, the
motif (grafted) region will always maintain the same sequence while the rest of the protein will be completely re-designed.

A part from the usual scores that normally can be found in a Rosetta output, some residues of interest were labeled through the use of the
`DisplayPoseLabelsMover <https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/xsd/mover_DisplayPoseLabelsMover_type>`_, as described
in :ref:`how to read Rosetta outputs <readrosetta>`. We will pick some of those residues for further analysis.

As we aim to acquire some extra information, we will need to **expand** the behaviour of :func:`.parse_rosetta_file` through a :ref:`description <readrosetta>`.
Additionally, to facilitate reading, we are going to add to the description a list of score terms that we will *ignore* for the purposes of this demonstration.

.. note::
  After loading the data, we will add the previously loaded sequence as a ``reference_sequence`` for the :class:`.DesignFrame`. As we are working with a previously
  trimmed structure, we will also add a ``shift``, this being the number of the first residue of the structure (i.e. by default ``shift==1`` as that would be
  the first position). Adding this ``shift`` helps down the line in keeping the analysis and plotting numbered to the starting structure.

  Both ``reference_sequence`` and ``reference_shift`` need to be attached to the appropriate chain identifier.

.. ipython::
  :okwarning:

  In [3]: rules = {'scores_ignore': ['fa_*', 'niccd_*', 'hbond_*', 'lk_ball_wtd', 'pro_close', 'dslf_fa13', 'C_ni_rmsd_threshold',
     ...:                            'omega', 'p_aa_pp', 'yhh_planarity', 'ref', 'rama_prepro', 'time'],
     ...:          'sequence': 'C',
     ...:          'labels': ['MOTIF', 'SSE03', 'SSE05']}
     ...: df = rs.io.parse_rosetta_file('../rstoolbox/tests/data/input_ssebig.minisilent.gz', rules)
     ...: df.add_reference_sequence('C', baseline.get_sequence('C'))
     ...: df.add_reference_shift('C', 32)
     ...: df.head(3)

.. tip::
  Notice that ``scores_ignore`` accepts wildcard definitions.

Population Overview
-------------------

Before starting analyzing the sequences, one can and should take a look at the score terms of interest in order to properly sort/select the decoys that look more
interesting.

By default, ``rstoolbox`` does not bring functions to generate basic plots. As a derived class from :class:`~pandas.DataFrame`, :class:`.DesignFrame` is already compatible
with most of the popular plotting solutions in python, such as :mod:`matplotlib` or :mod:`seaborn`. Still, it adds some quick solutions to simple plots.

.. note::
  There are two main type of plotting functions in ``rstoolbox`` either they require a :class:`~matplotlib.figure.Figure` as they will plot multiple axes on it, or they
  require an :class:`~matplotlib.axes.Axes` as they only plot one figure. In the first case, the function will always return the :class:`~matplotlib.axes.Axes` objects that
  it has generated inside the :class:`~matplotlib.figure.Figure`. This way, one can always personalise all aspects of the plots.

.. ipython::
  :okwarning:

  In [4]: fig  = plt.figure(figsize=(30, 10))
     ...: grid = [2, 6]
     ...: axes = rs.plot.multiple_distributions( df, fig, grid )
     ...: plt.tight_layout()

  @savefig tutorial_seq_plt1.png width=10in
  In [5]: plt.show()

.. tip::
  For :func:`.multiple_distributions`, the grid has to define enough positions for all the scores expected to be plotted. In this case, there are 12 numerical columns and
  thus, that is the minimum number of cells necessary to plot.

The structure of the :class:`.DesignFrame` is such that all one vs. one scores can directly be seen just be providing it to :mod:`seaborn` (let's just select
some scores to minimize the final image):

.. ipython::
  :okwarning:

  In [6]: grid = sns.pairplot(df[['score', 'packstat', 'cav_vol', 'finalRMSD']])

  @savefig tutorial_seq_plt2.png width=10in
  In [7]: plt.show()

.. note::
  As a general rule, anything that a :class:`~pandas.DataFrame` can do, a :class:`.DesignFrame` also can.

Sorting and Selecting
---------------------

One key part of the analysis of protein design datasets is the correct preselection of those designs potentially interesting. A :class:`.DesignFrame` can operate with all the
same functions for sorting/selection that a :class:`~pandas.DataFrame` can.

Let's make two populations: ``df_cr`` representing the ``top 10 scored`` decoys with a ``cav_vol <= 10`` and a ``finalRMSD <= 1.5`` and
``df_sr`` with the ``top 10 scored`` decoys with a ``packstat >= 0.6`` and a ``finalRMSD <= 1.5``.

.. ipython::
  :okwarning:

  In [8]: df_cr = df[(df['cav_vol'] <= 10) & (df['finalRMSD'] <= 1.5)].sort_values('score').head(10)
     ...: df_sr = df[(df['packstat'] >= 0.6) & (df['finalRMSD'] <= 1.5)].sort_values('score').head(10)

As can be seen in this example, selections can be concatenated with either ``&`` (**AND** - meaning that all conditions must be fulfilled) or ``|`` (**OR** - meaning that
any of the conditions need to be met). After that, one can check the overlap between the two selected groups, either by checking the common identifiers as sets:

.. ipython::
  :okwarning:

  In [9]: print(set(df_cr['description'].values).intersection(df_sr['description'].values))

Or by merging only the common decoys, which will bring the full data:

.. ipython::
  :okwarning:

  In [10]: df_cr_sr = df_cr.merge(df_sr, how='inner', on='description', suffixes=['', '_2x'])
      ...: df_cr_sr.drop(list(df_cr_sr.filter(regex = '_2x')), axis = 1, inplace = True)
      ...: df_cr_sr

There are multiple options on how this two populations can be compared, if that might be needed to decide which is the one with a higher chance of success. Some of them are
showcased in their :ref:`specific tutorial <population_comparison>`.

Identification of Mutants
-------------------------

Identifying the mutations present in a population can be key to understand the sequence signatures behind the designs. It might be also a good starting point for new sets of designs
as is explained in :ref:`new_mutants`. Data referring to the identity and similarity with a ``reference_sequence`` is can easily be obtained by simple request (we will hide columns
not currently informative for this example):

.. ipython::
  :okwarning:

  In [11]: df_sr = df_sr.identify_mutants('C')
      ...: cols = ['description']
      ...: cols.extend(list(df_sr.filter(regex = '_C')))
      ...: df_sr[cols].head(2)

.. tip::

  Notice that the mutations are displayed on a **basic sequence count**; that is, assuming the protein sequence starts in 1. This is informative to keep on working on it,
  as numerical selection is performed that way, but to make a report the actual position might be more useful to get the count in the same range as the reference structure.
  One could get that format with :func:`.report()`.

.. ipython::
  :okwarning:

  In [12]: rs.utils.report(df_sr[['description']].head(2))

Thanks to this quick identification, it is easy to select decoys with a particular sequence signature of interest through :meth:`.DesignFrame.get_sequence_with`.

.. ipython::
  :okwarning:

  In [13]: df_sr.get_sequence_with('C', ((3, 'W'), (9, 'R')))

One important issue when selecting which decoys will move to the next round or be selected for experimental validation is the sequence distance between them,
to avoid (unless wanted) trying too similar designs:

.. ipython::
  :okwarning:

  In [14]: df_inner = df_sr.sequence_distance('C')
      ...: df_inner

Or even to compare the sequences between two different selection criteria:

.. ipython::
  :okwarning:

  In [15]: df_outer = df_sr.sequence_distance('C', df_cr)
      ...: df_outer

  In [16]: fig = plt.figure(figsize=(30, 10))
      ...: ax00 = plt.subplot2grid((1, 2), (0, 0))
      ...: _ = sns.heatmap(df_inner, ax=ax00)
      ...: rs.utils.add_top_title(ax00, 'Self Comparison')
      ...: ax01 = plt.subplot2grid((1, 2), (0, 1))
      ...: _ = sns.heatmap(df_outer, ax=ax01)
      ...: rs.utils.add_top_title(ax01, 'Populations Comparison')
      ...: plt.tight_layout()

  @savefig tutorial_seq_plt2b.png width=10in
  In [17]: plt.show()

Visual Analysis of Mutations
----------------------------

Amino acid variation over all the decoys can be visualised in a variety of ways, always with the ability to select only regions of interest:

* Through a :func:`.sequence_frequency_plot`:

.. ipython::
  :okwarning:

  In [18]: fig = plt.figure(figsize=(30, 10))
      ...: ax = plt.subplot2grid((1, 1), (0, 0))
      ...: rs.plot.sequence_frequency_plot(df_cr, 'C', ax, cbar=False)

  @savefig tutorial_seq_plt3.png width=10in
  In [19]: plt.show()

* As a :func:`.logo_plot` (let's see here only for the labeled region ``SSE03``):

.. ipython::
  :okwarning:

  In [20]: sse03 = df_cr.get_label('SSE03', 'C').values[0]
      ...: fig, axes = rs.plot.logo_plot(df_cr, 'C', key_residues=sse03)

  @savefig tutorial_seq_plt4.png width=10in
  In [21]: plt.show()

* As a :func:`.positional_sequence_similarity` (defaults to ``BLOSUM62`` similarity matrix):

.. ipython::
  :okwarning:

  In [22]: fig  = plt.figure(figsize=(30, 10))
      ...: ax = plt.subplot2grid((1, 1), (0, 0))
      ...: seqsim = rs.analysis.positional_sequence_similarity(df_cr, 'C')
      ...: rs.plot.positional_sequence_similarity_plot(seqsim, ax)

  @savefig tutorial_seq_plt5.png width=10in
  In [23]: plt.show()

* And, finally, as an alignment with :func:`.plot_alignment`:

.. ipython::
  :okwarning:

  In [24]: fig  = plt.figure(figsize=(30, 10))
      ...: ax = plt.subplot2grid((1, 1), (0, 0))
      ...: rs.plot.plot_alignment(df_cr, 'C', ax)

  @savefig tutorial_seq_plt6.png width=10in
  In [25]: plt.show()

Finally, a single selected decoy can be plotted individually with :func:`.per_residue_matrix_score_plot`;
in this case we will highlight the MOTIF ``region`` of interest and ``SSE05``:

.. ipython::
  :okwarning:

  In [26]: fig  = plt.figure(figsize=(30, 10))
      ...: ax = plt.subplot2grid((1, 1), (0, 0))
      ...: motif = df_cr.get_label('MOTIF', 'C').values[0]
      ...: sse05 = df_cr.get_label('SSE05', 'C').values[0]
      ...: seles = [(motif, 'red'), (sse05, 'blue')]
      ...: rs.plot.per_residue_matrix_score_plot(df_cr.iloc[0], 'C', ax, selections=seles)

  @savefig tutorial_seq_plt7.png width=10in
  In [27]: plt.show()

  In [28]: plt.close('all')
