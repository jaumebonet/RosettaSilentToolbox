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

These sets are only provided as help, but the user can create its own sets, as long as they are loaded as a :class:`~pandas.DataFrame`.

Cleaning your own background sets
---------------------------------

In case one has its own background reference dataset, **rstoolbox** offers the option to retrieve homology clusters precalculated
from PDB with the :func:`.make_redundancy_table` function. By default, this function will call on the PDB's ftp to download those
homology cluster files, but with ``precalculated==True`` it will directly provide a quicker access to the table available in the
library (as the first option is relatively time consuming).

The generated table can be used on your reference set as long as a ``pdb`` and ``chain`` column exists.

The following example shows now this table can be applied in what would be the equivalent of ``rs.utils.load_refdata('scop2', 30)``
and creates a visual representation of the clusters.

.. ipython::
  :okwarning:

  In [1]: hmdf = rs.utils.make_redundancy_table(precalculated=True)
     ...: scop2_30 = scop2.merge(hmdf, on=['pdb', 'chain'])
     ...: scop2_30 = scop2_30.sort_values('score').groupby('c30').first().reset_index()
     ...: fig = plt.figure(figsize=(20, 20))
     ...: ax = plt.subplot2grid((1, 1), (0, 0), fig=fig)
     ...: _ = sns.heatmap(hmdf.drop(columns=['pdb', 'chain']).sort_values(['c30', 'c40', 'c50', 'c70', 'c90', 'c95', 'c100']), ax=ax)

  @savefig tutorial_ctx_plt1.png width=10in
  In [2]: plt.show()

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

  @savefig tutorial_ctx_plt2.png width=10in
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
     ...:               (scop2['length'] <= slength + 5) &
     ...:               (scop2['score'] <= 100)]
     ...: fig = plt.figure(figsize=(30, 10))
     ...: axs = rs.plot.multiple_distributions(df, fig, grid, refdata=refdf, violins=False,
     ...:                                      ref_equivalences={'cavity': 'cav_vol', 'pack': 'packstat'})
     ...: plt.tight_layout()

  @savefig tutorial_ctx_plt3.png width=10in
  In [6]: plt.show()

Selections in context
---------------------

Contrary to the previous example, in which full distributions where compared, another useful tool is the ability to place selected
structures in a given context. This can have two main uses: **(1)** see the position of selected decoys in a given distribution
or **(2)** see the quality of putative template scaffolds before working with them.

.. ipython::
  :okwarning:

  In [1]: pickdf = df.sample(10)
     ...: fig = plt.figure(figsize=(20, 5))
     ...: axs = rs.plot.plot_in_context(pickdf, fig, (1, 4), refdata=refdf, values= ['score', 'packstat', 'cav_vol', 'BUNS'],
     ...:                               ref_equivalences={'cavity': 'cav_vol', 'pack': 'packstat'})
     ...: plt.tight_layout()

  @savefig tutorial_ctx_plt4.png width=10in
  In [6]: plt.show()

Each selection in its context
-----------------------------

Multiple selected decoys or scaffolds originated from different sources might not be all comparable under the same set. The most
simple case would be when those structures do not have a similar length, which does affect **Rosetta** scores.
In those scenarios, it is possible to see, for each decoy of interest, how well it compares to its reference dataset, and rank
them according to which quantile of the distribution they belong to:

.. ipython::
  :okwarning:

  In [1]: df = rs.utils.load_refdata('scop')
     ...: qr = pd.DataFrame([['2F4V', 'C'], ['3BFU', 'B'], ['2APJ', 'C'],
     ...:                    ['2C37', 'V'], ['2I6E', 'H']], columns=['pdb', 'chain'])
     ...: qr = qr.merge(df, on=['pdb', 'chain'])
     ...: refs = []
     ...: for i, t in qr.iterrows():
     ...:   refs.append(df[(df['length'] >= (t['length'] - 5)) &
     ...:                  (df['length'] <= (t['length'] + 5))])
     ...: fig  = plt.figure(figsize=(22, 7))
     ...: rs.plot.distribution_quality(df=qr, refdata=refs,
     ...:                              values=['score', 'pack', 'avdegree', 'cavity', 'psipred'],
     ...:                              ascending=[True, False, True, True, False],
     ...:                              names=['domain_id'], fig=fig)
     ...: plt.tight_layout()

  @savefig tutorial_ctx_plt5.png width=10in
  In [6]: plt.show()
