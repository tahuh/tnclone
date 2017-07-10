#!/usr/bin/python

bases_table = {'A': 'T', 'G' : 'C' , 'C' : 'G' , 'T' : 'A' , 'N' : 'N'}
codon_table = {
    'A': ('GCT', 'GCC', 'GCA', 'GCG'),
    'C': ('TGT', 'TGC'),
    'D': ('GAT', 'GAC'),
    'E': ('GAA', 'GAG'),
    'F': ('TTT', 'TTC'),
    'G': ('GGT', 'GGC', 'GGA', 'GGG'),
    'I': ('ATT', 'ATC', 'ATA'),
    'H': ('CAT', 'CAC'),
    'K': ('AAA', 'AAG'),
    'L': ('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'),
    'M': ('ATG',),
    'N': ('AAT', 'AAC'),
    'P': ('CCT', 'CCC', 'CCA', 'CCG'),
    'Q': ('CAA', 'CAG'),
    'R': ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
    'S': ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
    'T': ('ACT', 'ACC', 'ACA', 'ACG'),
    'V': ('GTT', 'GTC', 'GTA', 'GTG'),
    'W': ('TGG',),
    'Y': ('TAT', 'TAC'),
    '*': ('TAA', 'TAG', 'TGA'),
}
REVERSED_CODON = {
	'CTT': 'L', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 
	'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 
	'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'ACA': 'T', 
	'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 
	'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 
	'ACG': 'T', 'CAA': 'Q', 'AGT': 'S', 'CAG': 'Q', 
	'CCG': 'P', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 
	'TGT': 'C', 'CGA': 'R', 'CCA': 'P', 'TCT': 'S', 
	'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 
	'GGG': 'G', 'TAG': '*', 'GGA': 'G', 'TAA': '*', 
	'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S', 
	'TTA': 'L', 'GAC': 'D', 'TCC': 'S', 'GAA': 'E', 
	'TCA': 'S', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 
	'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'TTC': 'F', 
	'GTT': 'V', 'GCT': 'A', 'ACC': 'T', 'TGA': '*', 
	'TTG': 'L', 'CGT': 'R', 'TGG': 'W', 'CGC': 'R'
}

def rev_comp(seq) :
	c = ''
	for b in seq :
		c = bases_table[b] + c
	return c

def comp(seq) :
	c = ''
	for b in seq :
		c += bases_table[c]
	return c

def reverse(seq) :
	return seq[::-1]

def translate(seq):
	seq = seq
	if len(seq) % 3 != 0 :
		good_length = len(seq) - len(seq) % 3
		seq = seq[:good_length]
	acid = ""
	for i in xrange(0,len(seq) ,3) :
		codon = seq[i:i+3]
		a = REVERSED_CODON[codon]
		acid += a

	return acid
	
