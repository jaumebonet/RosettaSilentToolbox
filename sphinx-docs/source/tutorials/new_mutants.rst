.. _new_mutants:

.. currentmodule:: rstoolbox

Generation of New Variants
==========================

A first round of designs might just be the stepping stone towards a second generation. It can be used to learn and better stir the next generation according to whatever
is our final aim.

.. note::
  All the examples here will generate new sequences. Once those new sequences are generated, we can generate with a call to :meth:`.DesignFrame.make_resfile`
  the residue files that can be provided to Rosetta through the
  `ReadResfile <https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/ReadResfileOperation>`_ to guide
  the design process.
  We are not calling the method in this tutorial as it generates files which cannot be shown.


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
     ...: pd.set_option('display.width', 1000)
     ...: pd.set_option('display.max_columns', 500)
     ...: pd.set_option("display.max_seq_items", 3)
     ...: baseline = rs.io.get_sequence_and_structure('../rstoolbox/tests/data/2pw9C.pdb')
     ...: baseline.get_sequence('C')
     ...: baseline.add_reference_sequence('C', baseline.get_sequence('C'))
     ...: baseline.add_reference_shift('C', 32)

Loading the Design Data
-----------------------

Again, we are mimicking :ref:`sequence_analysis`.

.. ipython::
  :okwarning:

  In [2]: rules = {'scores_ignore': ['fa_*', 'niccd_*', 'hbond_*', 'lk_ball_wtd', 'pro_close', 'dslf_fa13', 'C_ni_rmsd_threshold',
     ...:                            'omega', 'p_aa_pp', 'yhh_planarity', 'ref', 'rama_prepro', 'time'],
     ...:          'sequence': 'C',
     ...:          'labels': ['MOTIF', 'SSE03', 'SSE05']}
     ...: df = rs.io.parse_rosetta_file('../rstoolbox/tests/data/input_ssebig.minisilent.gz', rules)
     ...: df.add_reference_sequence('C', baseline.get_sequence('C'))
     ...: df.add_reference_shift('C', 32)
     ...: df.head(3)

Minimal Mutant
--------------

We've seen multiple ways to identify and view mutations in :ref:`sequence_analysis`. Let's imagine that we have identified a good decoy candidate but
we want to try all the putative back mutations available. Basically, we ask if a less mutated decoy will perform as well as the one found.

For this, we will take the best scored decoy and we will try to :func:`.generate_wt_reversions` on the residues belonging to strand 3 (``SSE03``) and 5 (``SSE05``):

.. ipython::
  :okwarning:

  In [3]: kres = df.get_label('SSE03', 'C').values[0] + df.get_label('SSE05', 'C').values[0]
     ...: ds = df.sort_values('score').iloc[0].generate_wt_reversions('C', kres)
     ...: ds.shape[0]

  In [4]: ds.head(4)

.. note::
  To avoid confusion, all scores are unattached from the sequences, as they will need to be recalculated.

Mind that this will create all the combinatorial options from the selected mutant to the original wild type for the selected residues.

.. ipython::
  :okwarning:

  In [5]: _ = sns.distplot(ds['mutant_count_C'].values)

  @savefig tutorial_muts_plt1.png width=5in
  In [6]: plt.show()

Defining New Mutations of Interest
----------------------------------

Let's say that, after inspecting scores and visualising structures, there are some key positions for which seems relevant to try some particular residue types, that can also be controlled (also,
remember form :ref:`sequence_analysis` that we can select already decoys with some particular mutations with :meth:`.DesignFrame.get_sequence_with`). This functionality comes by the
hand of :meth:`.DesignFrame.generate_mutant_variants`. Let's now try it for the two best scored decoys.

.. ipython::
  :okwarning:

  In [5]: mutants = [(20, "AIV"), (31, "EDQR")]
     ...: ds = df.sort_values('score').iloc[:2].generate_mutant_variants('C', mutants)
     ...: ds.shape[0]

  In [6]: ds.head(4)

Exhaustive Point Mutant Sampling
--------------------------------

An interesting approach could also be to create the appropriate data to exhaustively explore all possible point mutant in all possible
positions. This is extremely easy to do with the ``rstoolbox``, generating a total of ``(len(sequence) * 19) + 1`` variants:

.. ipython::
  :okwarning:

  In [7]: point_mutants = []
     ...: for i in range(1, len(baseline.get_reference_sequence('C')) + 1):
     ...:     mutation = [(i, "*")]
     ...:     point_mutants.append(baseline.generate_mutant_variants('C', mutation))
     ...: point_mutants = pd.concat(point_mutants).drop_duplicates('sequence_C')
     ...: print(','.join(point_mutants['mutants_C'].sample(frac=0.035)))  # show only some...

With small variations, this same loop can be used to generate this exhaustive scanning over particular residues or to only scan
specific residue types (for example to perform an Alanine Scanning).

Learning from Amino Acids Frequency
-----------------------------------

Instead of manually checking the positions and residue types of interest, one could learn from different sources to improve the type of sequences that one could get.
This new decoy sequences will be scored according to the frequency matrix they were generated from, but any decoy can be mapped to a sequence matrix with
:meth:`.DesignFrame.score_by_pssm`.

For instance, one could get the information content of the ``15 top scored`` decoys and try to obtain a certain amount of mutants (10) from the ``two best packed`` decoy that would
actually follow the statistical rules of that set with :meth:`.DesignFrame.generate_mutants_from_matrix`:

.. ipython::
  :okwarning:

  In [7]: df15 = df.sort_values('score').iloc[:15].sequence_frequencies('C')
     ...: df.sort_values('packstat', ascending=False).iloc[:2].generate_mutants_from_matrix('C', df15, 10)

.. note::
  While the generation of **the best scored sequence** according to a matrix is quite straight forward, the generation of the N best is not. Thus,
  the way a frequency matrix is applied in this scenario is that the frequency corrects a randomised selection per position. Statistically, this will make
  most sequences obtained on the *upper side* of the matrix score.

Creating the New Mutants
------------------------

After generating new variations, those can be run in **Rosetta** to obtain their corresponding scores with :class:`.DesignFrame.apply_resfile`. By default, this function will run
*fixed backbone* design with a script such as:

.. ipython::
  :okwarning:

  In [8]: print(rs.utils.mutations())

But the user can provide more complex **RosettaScripts** to execute, as long as they follow the same restrictions as this one:

#. They contain the `AddJobPairData Mover <https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/AddJobPairDataMover>`_.
#. They target the **resfile** with the ``script_var`` ``%%resfile%%``.

