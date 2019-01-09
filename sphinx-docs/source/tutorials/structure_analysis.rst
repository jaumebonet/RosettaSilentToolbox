.. _structure_analysis:

.. currentmodule:: rstoolbox

Structural Analysis
===================

Structural data can be obtained in multiple ways, and it should be relatively easy to implement a  `DSSP <https://swift.cmbi.umcn.nl/gv/dssp/index.html>`_
or `PSIPRED <http://bioinf.cs.ucl.ac.uk/psipred/>`_ parser to include that data into the decoy population.

In line of a Rosetta pipeline, the easiest way is to include the
`WriteSSEMover <https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/xsd/mover_WriteSSEMover_type>`_ into a RosettaScript. This could
even be called from inside the library as explained in :ref:`sequence_analysis`.

Loading a Reference
-------------------

As in :ref:`sequence_analysis`, we will need to load a reference with :func:`.get_sequence_and_structure`.

.. note::
  Through all the process several times the ``chainID`` of the decoy of interest will be called. This is due to the fact that the library can manipulate
  decoys with multiple chains (designed or not), and, thus, analysis must be called upon the sequences of interest.

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
     ...: baseline.get_structure('C')

Loading the Design Data
-----------------------

Again, we are mimicking :ref:`sequence_analysis`, but in this case we load the structural data of the designs along with their dihedrals.

.. ipython::
  :okwarning:

  In [2]: rules = {'scores_ignore': ['fa_*', 'niccd_*', 'hbond_*', 'lk_ball_wtd', 'pro_close', 'dslf_fa13', 'C_ni_rmsd_threshold',
     ...:                            'omega', 'p_aa_pp', 'yhh_planarity', 'ref', 'rama_prepro', 'time'],
     ...:          'sequence': 'C', 'structure': 'C', 'dihedrals': 'C',
     ...:          'labels': ['MOTIF', 'SSE03', 'SSE05']}
     ...: df = rs.io.parse_rosetta_file('../rstoolbox/tests/data/input_ssebig.minisilent.gz', rules)
     ...: df.add_reference_structure('C', baseline.get_structure('C'))
     ...: df.add_reference_shift('C', 32)
     ...: df.head(3)


Select by Structural Properties
-------------------------------

One of the most obvious assessments is the calculation of percentages of secondary structure. We can even use that to pre-select decoys. For example,
let's assume that we want to make sure that strand 5 (``SSE05``) from the selected labels is mostly still a strand after our design protocol:

.. ipython::
  :okwarning:

  In [2]: sse05 = df.get_label('SSE05', 'C').values[0]
     ...: df = rs.analysis.secondary_structure_percentage(df, 'C', sse05)
     ...: dfsse = df[df['structure_C_E'] > 0.8]
     ...: dfsse.head(3)

Or that we want to check if there is a correlation between the final score of our decoys and the level of conservation of strand 5
(which does not seem to be the case):

.. ipython::
  :okwarning:

  In [2]: grid = sns.pairplot(df[['score', 'structure_C_E']])

  @savefig tutorial_str_plt1.png width=7in
  In [3]: plt.show()


Visualise Structural Data
-------------------------

The most global view, would be to see for all the population how often the expected secondary structure is retrieved with
:func:`.positional_structural_similarity_plot`.

.. ipython::
  :okwarning:

  In [4]: fig  = plt.figure(figsize=(30, 10))
     ...: ax = plt.subplot2grid((1, 1), (0, 0))
     ...: sseid1 = rs.analysis.positional_structural_count(df, 'C')
     ...: sseid2 = rs.analysis.positional_structural_identity(df, 'C')
     ...: rs.plot.positional_structural_similarity_plot(pd.concat([sseid1, sseid2], axis=1), ax)

  @savefig tutorial_str_plt2.png width=10in
  In [5]: plt.show()

In case of a single decoy, and thanks to our request for ``dihedrals`` in the :ref:`parsing rules <readrosetta>`, we can check whether or not
a the decoy's `Ramachandran Plot <https://www.wikiwand.com/en/Ramachandran_plot>`_ makes sense. We can observe it, for example, for the best
scored decoy under the :func:`.plot_ramachandran` function.

.. note::
  Dihedral angles are also obtained through `WriteSSEMover <https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/xsd/mover_WriteSSEMover_type>`_,
  but they can easily be added by parsing them into a list and adding them to the :class:`.DesignFrame` with the column name ``psi_<seqID>`` and ``phi_<seqID>``.

.. ipython::
  :okwarning:

  In [6]: fig  = plt.figure(figsize=(30, 30))
     ...: _ = rs.plot.plot_ramachandran(df.sort_values('score').iloc[0], 'C', fig)

  @savefig tutorial_str_plt3.png width=10in
  In [7]: plt.show()

Finally, one can compare the actual DSSP assignation against PSIPRED's sequence-based secondary structure prediction. This can give an idea of the feasibility of the sequence to actually
fold in the expected conformation.

For this, we are going to load a different dataset in which this data is present.

.. ipython::
  :okwarning:

  In [8]: fig = plt.figure(figsize=(20, 5))
     ...: rules = {"scores": ["score"], "psipred": "*", "structure": "*", "dihedrals": "*" }
     ...: df2 = rs.io.parse_rosetta_file("../rstoolbox/tests/data/input_3ssepred.minisilent.gz", rules)
     ...: ax = plt.subplot2grid((1, 1), (0, 0))
     ...: rs.plot.plot_dssp_vs_psipred(df2.iloc[0], "A", ax)
     ...: plt.tight_layout()
     ...: fig.subplots_adjust(top=1.2)

  @savefig tutorial_str_plt4.png width=5in
  In [9]: plt.show()

  In [1]: plt.close('all')
