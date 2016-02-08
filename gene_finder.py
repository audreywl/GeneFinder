# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Audrey Lewis

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
	"""Shuffles the characters in the input string
		NOTE: this is a helper function, you do not
		have to modify this in any way """
	return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
	""" Returns the complementary nucleotide

		nucleotide: a nucleotide (A, C, G, or T) represented as a string
		returns: the complementary nucleotide
	>>> get_complement('A')
	'T'
	>>> get_complement('C')
	'G'
	"""
	if nucleotide=='T':
		return 'A'
	elif nucleotide=='A':
		return 'T'
	elif nucleotide=='G':
		return 'C'
	elif nucleotide=='C':
		return 'G'
	else:
		return 'Not a valid nucleotide'


def get_reverse_complement(dna):
	""" Computes the reverse complementary sequence of DNA for the specfied DNA
		sequence

		dna: a DNA sequence represented as a string
		returns: the reverse complementary DNA sequence represented as a string
	>>> get_reverse_complement("ATGCCCGCTTT")
	'AAAGCGGGCAT'
	>>> get_reverse_complement("CCGCGTTCA")
	'TGAACGCGG'
	"""
	index=-1
	length=-len(dna)
	result=''
	while index>=length:
		letter=dna[index]
		result+=get_complement(letter)
		index=index-1
	return result


def rest_of_ORF(dna):
	""" Takes a DNA sequence that is assumed to begin with a start
		codon and returns the sequence up to but not including the
		first in frame stop codon.  If there is no in frame stop codon,
		returns the whole string.

		dna: a DNA sequence
		returns: the open reading frame represented as a string
	>>> rest_of_ORF("ATGTGAA")
	'ATG'
	>>> rest_of_ORF("ATGAGATAGG")
	'ATGAGA'
	"""
	i=0
	while i<len(dna):
		if dna[i:i+3]=='TAG' or dna[i:i+3]== 'TGA' or dna[i:i+3]== 'TAA':
			break
		else:
			i=i+3
	return dna[:i]




def find_all_ORFs_oneframe(dna):
	""" Finds all non-nested open reading frames in the given DNA
		sequence and returns them as a list.  This function should
		only find ORFs that are in the default frame of the sequence
		(i.e. they start on indices that are multiples of 3).
		By non-nested we mean that if an ORF occurs entirely within
		another ORF, it should not be included in the returned list of ORFs.

		dna: a DNA sequence
		returns: a list of non-nested ORFs
	>>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
	['ATGCATGAATGTAGA', 'ATGTGCCC']
	"""
	i=0
	ORFs=[]
	while i<len(dna):
		if dna[i:i+3]=='ATG':
			newORF=rest_of_ORF(dna[i:])
			ORFs.append(newORF)
			i=i+len(newORF)	
		else:
			i=i+3
	return ORFs


def find_all_ORFs(dna):
	""" Finds all non-nested open reading frames in the given DNA sequence in
		all 3 possible frames and returns them as a list.  By non-nested we
		mean that if an ORF occurs entirely within another ORF and they are
		both in the same frame, it should not be included in the returned list
		of ORFs.

		dna: a DNA sequence
		returns: a list of non-nested ORFs

	>>> find_all_ORFs("ATGCATGAATGTAG")
	['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
	"""
	
	ORFs1=find_all_ORFs_oneframe(dna)
	ORFs2=find_all_ORFs_oneframe(dna[1:])
	ORFs3=find_all_ORFs_oneframe(dna[2:])
	ORFs=ORFs1+ORFs2+ORFs3	
	return ORFs


def find_all_ORFs_both_strands(dna):
	""" Finds all non-nested open reading frames in the given DNA sequence on both
		strands.

		dna: a DNA sequence
		returns: a list of non-nested ORFs
	>>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
	['ATGCGAATG', 'ATGCTACATTCGCAT']
	"""
	dnareverse=get_reverse_complement(dna)
	ORFs1=find_all_ORFs(dna)
	ORFs2=find_all_ORFs(dnareverse)
	ORFs=ORFs1+ORFs2
	return ORFs



def longest_ORF(dna):
	""" Finds the longest ORF on both strands of the specified DNA and returns it
		as a string
	>>> longest_ORF("ATGCGAATGTAGCATCAAA")
	'ATGCTACATTCGCAT'
	"""
	ORFs=find_all_ORFs_both_strands(dna)
	ORFs_max=0
	long_ORF=''
	for i in range(len(ORFs)):
		current_length=len(ORFs[i])
		if current_length>ORFs_max:
			ORFs_max=current_length
			long_ORF=ORFs[i]
	return long_ORF


def longest_ORF_noncoding(dna, num_trials):
	""" Computes the maximum length of the longest ORF over num_trials shuffles
		of the specfied DNA sequence

		dna: a DNA sequence
		num_trials: the number of random shuffles
		returns: the maximum length longest ORF """
	longest_length=0
	for i in range(1,num_trials):
		shuffled=shuffle_string(dna)
		newORF=longest_ORF(shuffled)
		if len(newORF)>longest_length:
			longest_length=len(newORF)
	return longest_length
			



def coding_strand_to_AA(dna):
	""" Computes the Protein encoded by a sequence of DNA.  This function
		does not check for start and stop codons (it assumes that the input
		DNA sequence represents an protein coding region).

		dna: a DNA sequence represented as a string
		returns: a string containing the sequence of amino acids encoded by the
				 the input DNA fragment

		>>> coding_strand_to_AA("ATGCGA")
		'MR'
		>>> coding_strand_to_AA("ATGCCCGCTTT")
		'MPA'
	"""
	aa_list=[]
	remainder=len(dna)%3
	dna=dna[:(len(dna)-remainder)]
	for i in range(0,len(dna),3):
		current_codon=dna[i:i+3]
		current_aa=aa_table[current_codon]
		aa_list.append(current_aa)
	aa_str=''.join(aa_list)
	return aa_str

def gene_finder(dna):
	""" Returns the amino acid sequences that are likely coded by the specified dna

		dna: a DNA sequence
		returns: a list of all amino acid sequences coded by the sequence dna.
	"""
	found_genes=[]
	threshold = longest_ORF_noncoding(dna, 1500)
	ORF_list=find_all_ORFs_both_strands(dna)
	for ORF in ORF_list:
		if len(ORF)>threshold:
			current_gene=coding_strand_to_AA(ORF)
			found_genes.append(current_gene)
	return found_genes

dna = load_seq("./data/X73525.fa")

print gene_finder(dna)

#longest_ORF_noncoding(dna, 1500)

if __name__ == "__main__":
	import doctest
	doctest.testmod()

