#!/usr/bin/python

"""
variant_call.py
This variant call is based on reference position

If insertion, the base report must be exported from contig sequence
If deletion, the base report must be exported from reference
Modification Note

2016-08-25
First written
"""

import argparse
import re
from SequenceParser import FAParser

parser = argparse.ArgumentParser(description="Variant Call based on SW alignment")

parser.add_argument("--sam", type=str, required=True, help="SAM file input")
parser.add_argument("--vcf", type=str, required=True, help="VCF file output")
parser.add_argument("--ref" , type=str, required = True, help="Reference")
args = parser.parse_args()
samname = args.sam
vcfname = args.vcf
refname = args.ref
sam= open(samname)
vcf = open(vcfname,"w")
reference = refname

vcf.write("#VCF(simple version.)\n")
vcf.write("#INFO:The 'short' mark in INFO field means the query is shorter than reference. INDEL or large soft clipping\n")
vcf.write("#INFO:The 'homopoly' mark in INFO field means there is homopolymer region in refrence. These variant might be caused by sequencing errors.\n")
vcf.write("#REF\tQUERY\tPOS\tREF\tALT\tINFO\n")

#print "Processing sam file"

reflength = 0
refname = ''
vcf_dict = {}

ref_dict = {}


def cigar_alpha(c) :
	return re.split("[0-9]+",c)[1:]

def cigar_num(c) :
	return [int(x) for x in re.split("[A-Z=]+",c)[:-1]]

def homopoly(seq) :
	### more than 4-mers of homopolymer assumed to cause sequencing error
	return  [(x.group(), x.start()) for x in re.finditer(r'([ATGC])\1{2,}' , seq) if len(x.group()) >= 4]

def inside_homopoly(ref_pos, homopoly_region):
	for region in homopoly_region :
		if ref_pos >= region[0] and ref_pos <= region[1]:
			return True
	return False

parser = FAParser()
parser.open(reference)
for id , desc,s in parser.parse():
	# id = record.id
	# s = str(record.seq)
	ref_dict[id] = s

for line in sam :
	if line.startswith("@") :
		#if "SN:" in line :
			#refname = line.rstrip().split("SN:")[1].split()[0]
		if "LN:" in line:
			reflength = int(line.rstrip().split("LN:")[1])
			continue
		else:
			continue

	data = line.rstrip().split("\t")
	query = data[0]
	rname = data[2]
	pos = int(data[3])
	cigar = data[5]
	seq = data[9]

	alpha = cigar_alpha(cigar)
	nums = cigar_num(cigar)
	
	refvar_pos = pos - 1
	contig_pos = 0 ### Must be initialzied zero because the first starting point of contig is mapped
	cigar_zip = zip(alpha, nums)
	refseq = ref_dict[rname]
	short = False
	if len(seq) < len(refseq):
		short = True
	homopolymer_set = homopoly(refseq)
	homopoly_region = []
	for h in homopolymer_set :
		start = h[0]
		end = h[1] + len(h[0])
		homopoly_region.append((start,end))

	for idx , z in enumerate(cigar_zip) :
		char = z[0]
		n = z[1]
		if char == 'M' or char == '=' :
			# Match
			refvar_pos += n
			contig_pos += n
		if char == 'S' :
			# Soft-clipping
			# These are not part of alignment
			# Must be skipped
			if idx == 0 :
				# case 1 : if soft clips are in the first cigar alphabet
				# Simply trim sequence

				seq = seq[n:]
			elif idx == len(cigar_zip) - 1:
				# Case 2 : if soft clip are at the end of cigar alphabet
				# Simply trim sequence
				seq = seq[:len(seq)-n]

			else:
				# Case 3 : if soft clips are in the middle of cigar alphabet
				# This is complicated task
				nums_before_soft = nums[:idx]
				nums_after_soft = nums[idx+1:]
				sum_num_before = sum(nums_before_soft)
				sum_num_after = sum(nums_after_soft)

				front_seq = seq[:sum_num_before]
				rear_seq = seq[sum_num_before + n:]
				seq = front_seq + rear_seq
		if char == 'I' :
			# Insertion
			# sequence must be retrieved from contig
			refbase = '-'
			ins_length = n
			alt_base = seq[contig_pos : contig_pos + n] # say 17I at pos 3 (1-based) of sequence ATTTTAGATGGGATAGATAGATAGATGA is TTTAGATGGGATAGATA
			judge_homopoly_in = inside_homopoly(refvar_pos , homopoly_region)
			if judge_homopoly_in :
				vcf.write(rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + alt_base + "\tHOMOPOLYMER_REGION\n")
			else:
				if short :
					vcf.write(rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + alt_base + "\tSHORT\n")
				else:
					vcf.write(rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + alt_base + "\t*\n")
			
			contig_pos += n
		if char == 'D' :
			# Deletion
			# Sequence must be retrieved from reference
			alt_base = '-'
			del_length = n
			refseq = ref_dict[rname]
			refbase = refseq[refvar_pos : refvar_pos + n]
			judge_homopoly_in = inside_homopoly(refvar_pos , homopoly_region)
			if judge_homopoly_in :
				vcf.write(rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + alt_base + "\tHOMOPOLYMER_REGION\n")
			else:
				if short :
					vcf.write(rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + alt_base + "\tSHORT\n")
				else:
					vcf.write(rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + alt_base + "\t*\n")
			refvar_pos += n
		if char == 'X' :
			# Mismatch
			refseq = ref_dict[rname]
			i = refvar_pos
			tmp_contig_pos = contig_pos
			for i in range(refvar_pos, refvar_pos + n):
				refbase = refseq[i]
				alt_base = seq[tmp_contig_pos]
				judge_homopoly_in = inside_homopoly(i , homopoly_region)
				if judge_homopoly_in :
					vcf.write(rname + "\t" + query + "\t" + str(i) + "\t" + refbase + "\t" + alt_base + "\tHOMOPOLYMER_REGION\n")
				else:
					if short :
						vcf.write(rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + alt_base + "\tSHORT\n")
					else:
						vcf.write(rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + alt_base + "\t*\n")
				tmp_contig_pos += 1
				#refbase = refseq[refvar_pos]
				#altbase = seq[contig_pos]
				#vcf.write(rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + altbase + "\n")
			refvar_pos += n
			contig_pos += n
vcf.close()
