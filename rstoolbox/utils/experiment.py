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

# External Libraries

# This Library

__all__ = ['translate_dna_sequence','translate_3frames','adapt_length','sequencing_enrichment']

def translate_dna_sequence(sequence):
	"""
	:param str sequence: DNA sequence
	:return: str - protein sequence
	"""
    codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',}
    protein=""
    last_codon_start=len(sequence)-2
    for start in range(0,last_codon_start,3):
        codon=sequence[start:start+3]
        aa=codontable.get(codon, 'X')
        protein = protein + aa
    return protein


def translate_3frames(sequence, matches):
	"""
	:param str sequence: DNA sequence
	:param matches: sequence pattern to match
	:type matches: list of str
	
	:return: str
	"""
    protein_frame = []
    
    for i in range(0, 3):
        protein_frame.append(translate_dna_sequence(sequence[i:]))
        
    match_max = len(matches)
    match_count = [0,] * len(protein_frame)
    for i, p in enumerate(protein_frame):
        for m in matches:
            if re.search(m, p):
                match_count[i] += 1
    try:
        i = match_count.index(match_max)
        return protein_frame[i]
    except ValueError:
        return ""

def adapt_length(seqlist,start,stop):
	"""
	:param str seqlist: list of protein sequence
	:param str start: start pattern (not included in final sequence)
	:param str stop: stop pattern (not included in final sequence)
	:return: list of str
	"""
    for i, seq in enumerate(seqlist):
        m = re.search(start + '(.*)' + stop, seq)
        if m:
            seqlist[i] = m.group(1)
    return seqlist
    
def sequencing_enrichment(input, enrichment, matches, seqID='A'):
	"""
	:param dict input: first key is binder, second key is concentration, value is fastq file 
	:param dict enrichment: first key is binder, value is list of two concentrations (min,max) to calculate enrichment
	:param matches: sequence pattern to match
	:return: dataframe with sequence, counts (sequence) per fastq file, enrichment per binder 
	"""

	def translate_all(seqlist, matches):
    	for i, seq in enumerate(seqlist):
        	seqlist[i] = translate_3frames(seq, matches)
    	return seqlist

	def condition_reader(jobid, filename, matches):
    	data = read_fastq(filename)
    	data = translate_all(data, matches)
    	data = adapt_length(data)
    	counts = Counter(data)
    	df = pd.DataFrame.from_dict(counts, orient='index')
    	df = df.sort_values(0, ascending=False)
    	df.reset_index(level=0, inplace=True)

    	return df.rename(columns={'index': 'seq', 0:jobid})
	
	def binder_reader(jobid, inputb, matches):
    	data = []
    	for cond in inputb:
        	data.append(condition_reader(jobid + '_' + cond, inputb[cond], matches))
    	df = reduce(lambda left, right: pd.merge(left, right, on='seq', how='outer'), data).fillna(0)
    	return df

    data = []
    for binder in input:
        data.append(binder_reader(binder, input[binder], matches))
    df = reduce(lambda left, right: pd.merge(left, right, on='seq', how='outer'), data).fillna(0)
    df['len'] = df.apply(lambda row: len(row['seq']), axis=1)
    for binder in enrichment:
    	# {'5C4': ['10nMFab', '1uM'], 'D25': ['100pM', '10nM']}
    	df['enrichment_{}'.format(binder)] = df['{0}_{1}'.format(binder, enrichment[binder][0])/df['{0}_{1}'.format(binder, enrichment[binder][1])]
    df.replace({'inf': '-1'}, regex=True)
    designf = rstoolbox.components.DesignFrame(df.rename(columns={'seq': 'sequence_{}'.format(seqID)}))
	designf = designf.reset_index().rename(columns={'index':'description'})
    return designf

    
    
    