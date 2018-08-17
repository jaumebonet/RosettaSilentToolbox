# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: random_frequency_matrix
.. func:: random_fastq
"""
# Standard Libraries
from random import randint, choice

# External Libraries
import pandas as pd
import numpy as np

# This Library
from rstoolbox.components import DesignFrame, Selection

__all__ = ['random_frequency_matrix', 'random_fastq', 'random_proteins']


def random_frequency_matrix(size, seed=None):
    """Generate a random frequency matrix.

    :param int size: How many rows should the matrix have.
    :param int seed: If provided, set a fixed seed.

    :return: :class:`~pandas.DataFrame`
    """
    if seed is not None:
        np.random.seed(seed)
    alphabet = list("ARNDCQEGHILKMFPSTWYV")
    data = []
    for _ in range(size):
        data.append(np.random.dirichlet(np.ones(20), size=1)[0])
    return pd.DataFrame(data, columns=alphabet)


def random_proteins(size, count):
    """Generate random protein sequences.

    :param int size: Length of the sequences.
    :param int count: Number of sequences.

    :return: :class:`.DesignFrame`
    """
    from rstoolbox.components import DesignFrame

    def make_sequence(size):
        alphabet = list("ARNDCQEGHILKMFPSTWYV")
        return ''.join(choice(alphabet) for _ in range(size))

    df = DesignFrame({'description': ['decoy_{:04d}'.format(x + 1) for x in range(count)]})
    df['sequence_A'] = df.apply(lambda row: make_sequence(size), axis=1)
    return df


def random_fastq(sequence, description, selection, btags, num_seqs, min_repeat,
                 max_repeat, num_files, prefix):
    """Generate a requested number of fastq files.

    :param str sequence: Starting protein sequence.
    :param str description: Name of the sequence.
    :param str selection: Region of the sequence to keep untouched.
    :param btags: Protein sequence border tags.
    :type btags: :func:`list` of :class:`str`
    :param int num_seqs: Number of individual sequences to generate.
    :param int min_repeat: Minimum number of repetitions per sequence.
    :param int max_repeat: Maximum number of repetitions per sequence.
    :param int num_files: Number of files to generate.
    :param str prefix: Prefixes for the files.

    :return: :func:`list` of :class:`str` - filenames generated
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
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'}
    revtable = {}
    for c in codontable:
        revtable.setdefault(codontable[c], []).append(c)

    def rtranslate(seq, tab):
        return "".join([tab[_][randint(0, len(tab[_]) - 1)] for _ in seq])

    def add_tag(row, start, stop):
        return start + row['dna_A'] + stop if row.name % 2 == 0 else row['dna_A']

    # Create DesignFrame
    df = DesignFrame({'description': [description], 'sequence_A': sequence})
    df.add_reference_sequence('A', df.get_sequence('A').values[0])

    # Create Random Sequences; keep a segment as MOTIF
    mtx = random_frequency_matrix(len(df.get_reference_sequence('A')), 0)
    sel = Selection(selection)
    mutants = df.iloc[0].generate_mutants_from_matrix('A', mtx, num_seqs, ~sel)[0]

    # Generate DNA
    mutants['dna_A'] = mutants.apply(lambda row: rtranslate(row['sequence_A'], revtable), axis=1)

    # Add a start and stop tags only to pair sequences
    start = "".join([revtable[_][randint(0, len(revtable[_]) - 1)] for _ in btags[0]])
    stop  = "".join([revtable[_][randint(0, len(revtable[_]) - 1)] for _ in btags[1]])

    filename = prefix + '_{:03d}.fasq'
    all_files = []
    for i in range(num_files):
        fname = filename.format(i + 1)
        all_files.append(fname)
        # Add weights and repeat
        mut = mutants.copy()
        mut['weight'] = [randint(1, 5) for _ in range(mut.shape[0])]
        mut = mut.loc[mut.index.repeat(mut.weight)].sample(frac=1).reset_index(drop=True)
        # Add border tags
        mut['dna_A'] = mut.apply(lambda row: add_tag(row, start, stop), axis=1)
        # Write file
        with open(fname, 'w') as fd:
            for index, row in mut.iterrows():
                fd.write('@{};TEST_MAKEUP\n'.format(row['description']))
                fd.write(row['dna_A'] + '\n')
                fd.write('+\n')
                fd.write('?A@EC?C@AC=B>A@??DEC?EEC@C@DDD:\n')
    return all_files
