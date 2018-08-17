# -*- coding: utf-8 -*-
"""
.. codeauthor:: Fabian Sesterhenn <sesterhenn.fabian@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: translate_dna_sequence
.. func:: translate_3frames
.. func:: adapt_length
.. func:: sequencing_enrichment
"""
# Standard Libraries
import re

# External Libraries
import pandas as pd
import numpy as np
from six.moves import reduce

# This Library

__all__ = ['translate_dna_sequence', 'translate_3frames',
           'adapt_length', 'sequencing_enrichment']


def translate_dna_sequence( sequence ):
    """Translates **DNA** to **protein**.

    Assumes always that the codon starts in the first position
    of the sequence.

    :param str sequence: DNA sequence
    :return: :class:`str` - protein sequence

    .. seealso::
        :func:`.translate_3frames`

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import read_fastq
           ...: from rstoolbox.utils import translate_dna_sequence
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = read_fastq("../rstoolbox/tests/data/cdk2_rand_001.fasq.gz")
           ...: df.iloc[0]['sequence_A']

        In [1]: translate_dna_sequence(df.iloc[0]['sequence_A'])
    """
    codontable = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'}
    protein = ""
    last_codon_start = len(sequence) - 2
    for start in range(0, last_codon_start, 3):
        codon = sequence[start:start + 3]
        aa = codontable.get(codon, 'X')
        protein = protein + aa
    return protein


def translate_3frames( sequence, matches=None ):
    """Translates **DNA** to **protein** trying all possible frames.

    Tests the three possible reading frames. To decide which one to return,
    it can follow a double logic: when ``matches`` is :data:`None` it will
    return the longest sequence until a stop codon, otherwise it will return
    the longest sequence that contains the required match protein sequence. If
    none match, it will return an empty :class:`str`.

    All provided matches need to be found

    :param str sequence: DNA sequence
    :param matches: sequence pattern to match
    :type matches: :func:`list` of :class:`str`

    :return: :class:`str`

    .. rubric:: Example 1: No matches - Prepend nucleotides to change reading frame

    .. ipython::

        In [1]: from rstoolbox.io import read_fastq
           ...: from rstoolbox.utils import translate_3frames
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = read_fastq("../rstoolbox/tests/data/cdk2_rand_001.fasq.gz")
           ...: df.iloc[0]['sequence_A']

        In [1]: translate_3frames('AT' + df.iloc[0]['sequence_A'])

    .. rubric:: Example 2: With matches - Prepend nucleotides to change reading frame

    .. ipython::

        In [1]: from rstoolbox.io import read_fastq
           ...: from rstoolbox.utils import translate_3frames
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = read_fastq("../rstoolbox/tests/data/cdk2_rand_001.fasq.gz")
           ...: matches = ['GAS', 'FFG']
           ...: translate_3frames('AT' + df.iloc[0]['sequence_A'], matches)

        In [1]: translate_3frames('AT' + df.iloc[1]['sequence_A'], matches)
    """
    protein_frame = []

    for i in range(0, 3):
        protein_frame.append(translate_dna_sequence(sequence[i:]))

    protein_frame.sort(key=lambda s: len(s.split('_')[0]), reverse=True)

    if matches is None:
        return protein_frame[0]

    match_max = len(matches)
    match_count = [0, ] * len(protein_frame)
    for i, p in enumerate(protein_frame):
        for m in matches:
            if re.search(m, p):
                match_count[i] += 1
    try:
        i = match_count.index(match_max)
        return protein_frame[i]
    except ValueError:
        return ""


def adapt_length( seqlist, start, stop, inclusive=False ):
    """Pick only the sequence between the provided pattern tags.

    When ``inclusive`` is :data:`False` and the boundary tags are
    not found, the original sequence is returned, as it is assumed
    that the tags were out of the boundary of the retrieved sequence.

    When ``inclusive`` is :data:`True` and the boundary tags are
    not found, an empty sequence is returned for that position, as
    we understand that the interest was of getting them too and we could
    not.

    :param str seqlist: list of protein sequence
    :param str start: start pattern (not included in final sequence)
    :param str stop: stop pattern (not included in final sequence)
    :param bool inclusive: If :data:`False`, retrieve sequence **between**
        the protein tags, otherwise include the protein tags in the
        returned sequence.
    :return: :func:`list` of :class:`str`

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import read_fastq
           ...: from rstoolbox.utils import translate_dna_sequence
           ...: from rstoolbox.utils import adapt_length
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: df = read_fastq("../rstoolbox/tests/data/cdk2_rand_001.fasq.gz")
           ...: df['sequence_A'] = df.apply(lambda row: translate_dna_sequence(row['sequence_A']),
           ...:                             axis=1)
           ...: bounds = ['GAS', 'FFG']
           ...: df['sequence_A'].values[:5]

        In [1]: adapt_length(df['sequence_A'].values[:5], bounds[0], bounds[1])

        In [1]: adapt_length(df['sequence_A'].values[:5], bounds[0], bounds[1], True)

    """
    expression = '{0}(.*){1}'.format(start, stop)
    if inclusive:
        expression = '({0}.*{1})'.format(start, stop)

    outlist = []
    for seq in seqlist:
        m = re.search(expression, seq)
        if m:
            outlist.append(m.group(1))
        elif inclusive:
            outlist.append("")
        else:
            outlist.append(seq)
    return outlist


def sequencing_enrichment( indata, enrichment=None, bounds=None, matches=None, seqID='A' ):
    """Retrieve data from multiple
    `NGS <https://www.wikiwand.com/en/DNA_sequencing#/Next-generation_methods>`_ files.

    Allows to obtain data from multiple files while ataching them to two conditions, a primary one
    (key1) and a secondary one (key2).

    For instance, let's assume that one has data obtained through selection of sequences by two
    different binders and three different concentration of binder each; we would define a
    ``indata`` dictionary such as::

        {'binder1': {'conc1': 'file1.fastq', 'conc2': 'file2.fastq', 'conc3': 'file3.fastq'},
         'binder2': {'conc1': 'file4.fastq', 'conc2': 'file5.fastq', 'conc3': 'file6.fastq'}}

    Also, for each binder we could decide to calculate the enrichment between any two
    concentrations; we can do that by defining a ``enrichment`` dictionary such as::

        {'binder1': ['conc1', 'conc3'],
         'binder2': ['conc1', 'conc3']}

    :param dict indata: First key is binder, second key is concentration, value is fastq file.
    :param dict enrichment: Key is binder, value is list of two concentrations (min,max)
        to calculate enrichment.
    :param bounds: N and C limit of the sequences. Follow the logic of :func:`adapt_length`
        with ``inclusive`` as :data:`False`.
    :type bounds: :func:`list` of :class:`str`
    :param matches: Sequence pattern to match. Follows the same logic as in
        :func:`.translate_3frames`.
    :type matches: :func:`list` of :class:`str`
    :return: :class:`.DesignFrame` with the sequences, counts (sequence) per fastq file and
        enrichment per binder (if requested).

    .. rubric:: Example

    (We skip printing the sequence column to ease visibility of the differences)

    .. ipython::

        In [1]: from rstoolbox.io import read_fastq
           ...: from rstoolbox.utils import sequencing_enrichment
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 20)
           ...: indat = {'binder1': {'conc1': '../rstoolbox/tests/data/cdk2_rand_001.fasq.gz',
           ...:                      'conc2': '../rstoolbox/tests/data/cdk2_rand_002.fasq.gz',
           ...:                      'conc3': '../rstoolbox/tests/data/cdk2_rand_003.fasq.gz'},
           ...:          'binder2': {'conc1': '../rstoolbox/tests/data/cdk2_rand_004.fasq.gz',
           ...:                      'conc2': '../rstoolbox/tests/data/cdk2_rand_005.fasq.gz',
           ...:                      'conc3': '../rstoolbox/tests/data/cdk2_rand_006.fasq.gz'}}
           ...: df = sequencing_enrichment(indat)
           ...: df[[_ for _ in df.columns if _ != 'sequence_A']].head()

        In [1]: enrich = {'binder1': ['conc1', 'conc3'],
           ...:           'binder2': ['conc1', 'conc3']}
           ...: df = sequencing_enrichment(indat, enrich)
           ...: df[[_ for _ in df.columns if _ != 'sequence_A']].head()

    """
    from rstoolbox.components import DesignFrame

    def condition_reader(jobid, filename, bounds, matches):
        from rstoolbox.io import read_fastq
        df = read_fastq(filename)
        df['sequence_A'] = df.apply(lambda row: translate_3frames(row['sequence_A'], matches),
                                    axis=1)
        if bounds is not None:
            df['sequence_A'] = adapt_length(df['sequence_A'].values, bounds[0], bounds[1])

        df = df.merge(df.groupby('sequence_A').agg('count').reset_index(),
                      on='sequence_A',
                      how='left').drop_duplicates('sequence_A').reset_index(drop=True)

        df.rename(columns={'description_x': 'description', 'description_y': jobid}, inplace=True)
        return df.sort_values(jobid, ascending=False)

    def binder_reader(jobid, inputb, bounds, matches):
        data = []
        for cond in inputb:
            data.append(condition_reader(jobid + '_' + cond, inputb[cond], bounds, matches))
        df = reduce(lambda left, right: pd.merge(left, right, on='sequence_A', how='outer'),
                    data).fillna(0)
        return df

    data = []
    for binder in indata:
        data.append(binder_reader(binder, indata[binder], bounds, matches))
    df = reduce(lambda left, right: pd.merge(left, right, on='sequence_A', how='outer'),
                data).fillna(0)
    df['len'] = df.apply(lambda row: len(row['sequence_A']), axis=1)
    df = df.drop([_ for _ in df.columns if _.startswith('description')], axis=1)

    if enrichment is not None:
        for binder in enrichment:
            eb = enrichment[binder]
            id1 = '{0}_{1}'.format(binder, eb[0])
            id2 = '{0}_{1}'.format(binder, eb[1])
            df['enrichment_{}'.format(binder)] = df[id1] / df[id2]
    df = df.replace({np.inf: -1, -np.inf: -1}, regex=True).fillna(0)
    designf = DesignFrame(df.rename(columns={'sequence_A': 'sequence_{}'.format(seqID)}))
    designf = designf.reset_index().rename(columns={'index': 'description'})
    return designf
