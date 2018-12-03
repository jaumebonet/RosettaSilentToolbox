.. _contextualize:

.. currentmodule:: rstoolbox

Contextualize your Data
=======================

The generation of a new design population can provide a lot of data depending in the number of filters and metrics used to evaluate
them. Still, although these scores allow for the sorting and selection of decoys of interest, it does not always provide a
comprehensive analysis of that population with regard to the general population of known protein structures.

The **rstoolbox** provides some helper functions and datasets in order to properly assess this issue.

Pre-calculated datasets
-----------------------

There is a total of 4 different pre-calculated datasets that can be loaded directly with the **rstoolbox**; those are ``cath``, ``scop``,
``scop2`` and ``chain`` (for PDB chains). For this datasets, homology (as provided by the PDB database) can be used to avoid over-represented
structures. The access to these datasets is provided through the :func:`.load_refdata` function:

.. ipython::
  :okwarning:

  In [1]: import rstoolbox as rs
     ...: import pandas as pd
     ...: import matplotlib.pyplot as plt
     ...: import seaborn as sns
     ...: plt.rcParams['svg.fonttype'] = 'none' # When plt.savefig to 'svg', text is still text object
     ...: sns.set_style('whitegrid')
     ...: pd.set_option('display.width', 1000)
     ...: scop2 = rs.utils.load_refdata('scop2')
     ...: scop2_70 = rs.utils.load_refdata('scop2', 70)  # 70% homology filter
     ...: scop2_70.drop(columns=['sequence_$', 'structure_$', 'psi_$', 'phi_$', 'node_name']).head(5)

These seta are only provided as help, but the user can create its own sets, as long as they are loaded as a :class:`~pandas.DataFrame`.

Distribution in context
-----------------------

Let's assume the following population dataset:

.. ipython::
  :okwarning:

  In [2]: rules = {'scores_ignore': ['fa_*', 'niccd_*', 'hbond_*', 'lk_ball_wtd', 'pro_close', 'dslf_fa13', 'C_ni_rmsd_threshold',
     ...:                            'omega', 'p_aa_pp', 'yhh_planarity', 'ref', 'rama_prepro', 'time'],
     ...:          'sequence': 'C',
     ...:          'labels': ['MOTIF', 'SSE03', 'SSE05']}
     ...: df = rs.io.parse_rosetta_file('../rstoolbox/tests/data/input_ssebig.minisilent.gz', rules)
     ...: df.head(3)

As previously demonstrated in the :ref:`sequence analysis <sequence_analysis>` tutorial, we can preview the distribution of different
scores of the design population:

.. ipython::
  :okwarning:

  In [3]: fig  = plt.figure(figsize=(30, 10))
     ...: grid = [2, 6]
     ...: axes = rs.plot.multiple_distributions( df, fig, grid )
     ...: plt.tight_layout()

  @savefig tutorial_ctx_plt1.png width=10in
  In [4]: plt.show()

But we can also see the population behave against other protein of a similar size.

.. note::
  The function will pick up that columns with the same name are the same, for others that should be the same but are called
  differently, a translation dictionary can be provided, as is the case in this example with ``cavity volume`` and ``packstat``
  data.

.. ipython::
  :okwarning:

  In [5]: slength = len(df.iloc[0].get_sequence('C'))
     ...: refdf = scop2[(scop2['length'] >= slength - 5) &
     ...:               (scop2['length'] <= slength + 5)]
     ...: fig = plt.figure(figsize=(25, 10))
     ...: axs = rs.plot.multiple_distributions(df, fig, grid, refdata=refdf, violins=False,
     ...:                                      ref_equivalences={'cavity': 'cav_vol', 'pack': 'packstat'})
     ...: plt.tight_layout()

  @savefig tutorial_ctx_plt2.png width=10in
  In [6]: plt.show()

