.. _population_comparison:

.. currentmodule:: rstoolbox

Compare Decoy Populations
=========================

Comparing features between populations is one of the best ways to benchmark different conditions in an experiment.

This example works over two populations designed in the same manner, with the only difference that one was performed in the presence of a binder while the other
was design alone.

  We are ignoring some scores to facilitate general reading, but this particularity is not mentioned any more during the tutorial. Similarly, we are using only a sample of the 20K
  available sequences on each decoy population.

.. ipython::
  :okwarning:

  In [1]: import rstoolbox as rs
     ...: import pandas as pd
     ...: import numpy as np
     ...: pd.set_option('display.width', 1000)
     ...: pd.set_option('display.max_columns', 500)
     ...: pd.set_option("display.max_seq_items", 3)
     ...: binder = rs.io.parse_rosetta_file('../rstoolbox/tests/data/compare/binder.mini.gz', {'sequence': 'B', 'scores_ignore': ['time', 'fa_*', 'dslf_fa13', 'yhh_planarity', 'pro_close']}).sample(frac=0.05)
     ...: binder.columns.values

  In [2]: nobinder = rs.io.parse_rosetta_file('../rstoolbox/tests/data/compare/nobinder.mini.gz', {'sequence': 'B', 'scores_ignore': ['time', 'fa_*', 'dslf_fa13', 'yhh_planarity', 'pro_close']}).sample(frac=0.05)
     ...: nobinder.columns.values

Merging Populations
-------------------

Merging populations is a very useful way to a take advantage of the plotting and analysis tools compatible with :class:`~pandas.DataFrame` derived objects. In order
for it to be useful and possible, two conditions must meet:

The two populations should share column names
`````````````````````````````````````````````

Not all, but sure for those that are expected to be analysed comparatively. When they all come from a same experiment, that is relatively easy to obtain in the setup. otherwise, one can use the rename option to fix that.

.. ipython::
  :okwarning:

  In [3]: binder.rename(columns={'score': 'renamed_score'}).columns.values

There must be a column to identify to which population each decoy belongs
`````````````````````````````````````````````````````````````````````````

In this example, none exist. But if we take a look at the naming of the decoys:

.. ipython::
  :okwarning:

  In [4]: print("{0} <-> {1}".format(binder.iloc[0]['description'], nobinder.iloc[0]['description']))

We could see that, by loading the data not with the :ref:`definition <readrosetta>` ``{'sequence': 'B'}`` but with ``{'sequence': 'B', 'naming': ['', '', 'condition']}`` we would
have automatically obtained a column ``condition`` defining the binder presence.
As we did not do that, or in case that would not be possible, we will just add an extra column.

.. ipython::
  :okwarning:

  In [5]: binder = rs.utils.add_column(binder, 'condition', 'binder')
     ...: binder.columns.values

  In [6]: nobinder = rs.utils.add_column(nobinder, 'condition', 'nobinder')
     ...: nobinder.columns.values

Once all conditions are met, the populations can be merged:

.. ipython::
  :okwarning:

  In [7]: df = pd.concat([binder, nobinder], ignore_index=True)

Comparing Basic Scores
----------------------

The first and most straight forward analysis between populations is the comparison between their different scores distributions. For that,
:func:`.multiple_distributions` accepts any attribute of :func:`~seaborn.boxplot`, so that data can be properly separated.

.. ipython::
  :okwarning:

  In [7]: import matplotlib.pyplot as plt
     ...: import seaborn as sns
     ...: fig  = plt.figure(figsize=(35, 15))
     ...: grid = [4, 7]
     ...: axes = rs.plot.multiple_distributions( df, fig, grid, x='condition', hue='condition', showfliers=False )
     ...: nones = [x.legend_.remove() for x in axes]
     ...: plt.tight_layout()

  @savefig tutorial_duo_plt1.png width=10in
  In [8]: plt.show()

The data can even be directly used on :mod:`seaborn` functions:

.. ipython::
  :okwarning:

  In [9]: sns.set(font_scale=1)
     ...: sns.set_style("whitegrid")
     ...: _ = sns.pairplot(df[['design_score', 'GRMSD2Target', 'GRMSD2Template', 'LRMSD2Target', 'LRMSDH2Target', 'LRMSDLH2Target', 'condition']],
     ...:                  hue='condition')

  @savefig tutorial_duo_plt2.png width=10in
  In [1]: plt.show()

More specific comparisons can be performed depending on the actual aim of the project.

Sequence Distances
------------------

Let's assume that one wants to evaluate the spread of the two populations and how well they recover a given target reference sequence. For that,
we will need to provide a reference sequence.

.. ipython::
  :okwarning:

  In [2]: reference = rs.io.get_sequence_and_structure('../rstoolbox/tests/data/compare/4oyd.pdb', {'sequence': 'B'})
     ...: _ = df.add_reference_sequence('B', reference.iloc[0].get_sequence('B')[:-1])

A tutorial for :ref:`population sequence analysis <sequence_analysis>` is also available, but there are options that can be used when comparing two populations.

For instance, one can evaluate the internal sequence inside each decoy population and then compare it between the two:

.. ipython::
  :okwarning:

  In [3]: disbinder = df[df['condition'] == 'binder'].sequence_distance('B')
     ...: disbinder = disbinder.mask(np.triu(np.ones(disbinder.shape, dtype=np.bool_))) # mask diagonal and upper diagonal values
     ...: disbinder = disbinder.values.flatten() # to array
     ...: disbinder = disbinder[[~np.isnan(disbinder)]] # clean NaN
     ...:
     ...: disnobinder = df[df['condition'] == 'nobinder'].sequence_distance('B')
     ...: disnobinder = disnobinder.mask(np.triu(np.ones(disnobinder.shape, dtype=np.bool_))) # mask diagonal and upper diagonal values
     ...: disnobinder = disnobinder.values.flatten() # to array
     ...: disnobinder = disnobinder[[~np.isnan(disnobinder)]] # clean NaN

Or evaluate the distance between the two sequence populations:

.. ipython::
  :okwarning:

  In [4]: discomp = df[df['condition'] == 'binder'].sequence_distance('B', df[df['condition'] == 'nobinder'])
     ...: discomp = discomp.values.flatten() # to array
     ...: discomp = discomp[[~np.isnan(discomp)]] # clean NaN

Or the distance of each group to the reference sequence:

.. ipython::
  :okwarning:

  In [5]: disbinderref = rs.analysis.sequence_similarity(df[df['condition'] == 'binder'], 'B', matrix='IDENTITY')['identity_B_negative']
     ...: disnobinderref = rs.analysis.sequence_similarity(df[df['condition'] == 'nobinder'], 'B', matrix='IDENTITY')['identity_B_negative']

We can show all this together as:

.. ipython::
  :okwarning:

  In [6]: fig  = plt.figure(figsize=(17, 5))
     ...: grid = [1, 3]
     ...: ax00 = plt.subplot2grid(grid, (0, 0))
     ...: _ = sns.distplot( disbinder , color="skyblue", label="binder", ax=ax00)
     ...: _ = sns.distplot( disnobinder , color="red", label="nobinder", ax=ax00)
     ...: rs.utils.add_top_title(ax00, 'internal distance')
     ...: ax01 = plt.subplot2grid(grid, (0, 1))
     ...: _ = sns.distplot( discomp , color="green", label="binder vs. nobinder", ax=ax01)
     ...: rs.utils.add_top_title(ax01, 'comparative distance')
     ...: ax02 = plt.subplot2grid(grid, (0, 2))
     ...: _ = sns.distplot( disbinderref , color="skyblue", label="binder", ax=ax02)
     ...: _ = sns.distplot( disnobinderref , color="red", label="nobinder", ax=ax02)
     ...: rs.utils.add_top_title(ax02, 'distance to reference')
     ...: plt.tight_layout()

  @savefig tutorial_duo_plt3.png width=10in
  In [7]: plt.show()

To fully asses the difference between the two populations, we can evaluate the enrichment of different residue types in different positions. In
this scenario, we will measure the enrichment of variants in the binder-designed population against the other. This entails that high values
represent residue types more present in the binder-design population per position.

.. ipython::
  :okwarning:

  In [8]: result = rs.analysis.positional_enrichment(df[df['condition'] == 'binder'], df[df['condition'] == 'nobinder'], 'B')

When a position is represented in the first population but not in the second, the enrichment is :data:`~numpy.inf`. To visualize it,
one would need to transform it into a numerical value. In this case, we transform it into a value such as ``max + (max/2)``.

.. ipython::
  :okwarning:

  In [9]: maxval = result.replace(np.inf, -1).max().max()
     ...: maxval += (maxval / float(2))
     ...: result = result.replace(np.inf, maxval)
     ...: fig = plt.figure(figsize=(25, 10))
     ...: ax = plt.subplot2grid((1, 1), (0, 0))
     ...: rs.plot.sequence_frequency_plot(df, "B", ax, cbar=False, xrotation=90)

  @savefig tutorial_duo_plt4.png width=10in
  In [1]: plt.show()

Comparing enrichment between the populations (i.e. full sequence enrichment) is demonstrated in the :ref:`experimental analysis tutorial <experimental_analysis>`.