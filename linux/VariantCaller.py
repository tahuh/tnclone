#!/usr/bin/python

"""
class VaraintCaller

Author : Sunghoon Heo

For Windows OS class
"""

import re
from SequenceParser import FAParser

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
	
	
class VariantCaller:
	def __init__(self, samname, vcfname , refname, message):
		self.samname = samname
		self.vcfname = vcfname
		self.refname = refname
		self.message = message
	def call(self):
		sam = open(self.samname)
		vcf = open(self.vcfname , "w")
		reference = self.refname
		
		vcf.write("#VCF(simple version.)\n")
		
		
		#vcf.write("#INFO:The 'short' mark in INFO field means the query is shorter than reference. INDEL or large soft clipping. See the LN field in the header\n")
		#vcf.write("#INFO:The 'homopoly' mark in INFO field means there is homopolymer region in refrence. These variant might be caused by sequencing errors.\n")
		#vcf.write("#REF\tQUERY\tPOS\tREF\tALT\tINFO\n")
		
		parser = FAParser()
		parser.open(reference)
		
		reflength = 0
		refname = ''
		vcf_dict = {}

		ref_dict = {}
		ref_lens = []
		var_info_list = []
		query_len_info_list = []
		self.message.append("[TnClone::VariantCaller] Parsing reference...")
		for id , desc,s in parser.parse():
			# id = record.id
			# s = str(record.seq)
			ref_dict[id] = s
			st = "#Ref:%s LN:%d"%(id , len(s))
			ref_lens.append(st)
			
		self.message.append("[TnClone::VariantCaller] Parsing alignment result && calling...")
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
			short = False
			query_len_info_list.append("#Query:%s LN:%d"%(query,len(seq)))
			
			alpha = cigar_alpha(cigar)
			nums = cigar_num(cigar)
			
			refvar_pos = pos - 1
			contig_pos = 0 ### Must be initialzied zero because the first starting point of contig is mapped
			cigar_zip = zip(alpha, nums)
			refseq = ref_dict[rname]
			ref_length = len(refseq)
			if len(seq) < ref_length :
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
						#vcf.write(rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + alt_base + "\tHOMOPOLYMER_REGION\n")
						st = rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + alt_base + "\tHOMOPOLYMER_REGION"
						var_info_list.append(st)
					else:
						if short :
							#vcf.write(rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + alt_base + "\tSHORT\n")
							st = rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + alt_base + "\tSHORT"
							var_info_list.append(st)
						else:
							#vcf.write(rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + alt_base + "\t*\n")
							st=rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + alt_base + "\t*"
							var_info_list.append(st)
					
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
						#vcf.write(rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + alt_base + "\tHOMOPOLYMER_REGION\n")
						st = rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + alt_base + "\tHOMOPOLYMER_REGION"
						var_info_list.append(st)
					else:
						if short :
							#vcf.write(rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + alt_base + "\tSHORT\n")
							st = rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + alt_base + "\tSHORT"
							var_info_list.append(st)
						else:
							#vcf.write(rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + alt_base + "\t*\n")
							st=rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + alt_base + "\t*"
							var_info_list.append(st)
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
							# vcf.write(rname + "\t" + query + "\t" + str(i) + "\t" + refbase + "\t" + alt_base + "\tHOMOPOLYMER_REGION\n")
							st = rname + "\t" + query + "\t" + str(i) + "\t" + refbase + "\t" + alt_base + "\tHOMOPOLYMER_REGION"
							var_info_list.append(st)
						else:
							if short :
								# vcf.write(rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + alt_base + "\tSHORT\n")
								st = rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + alt_base + "\tSHORT"
								var_info_list.append(st)
							else:
								# vcf.write(rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + alt_base + "\t*\n")
								st = rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + alt_base + "\t*\n"
								var_info_list.append(st)
						tmp_contig_pos += 1
						#refbase = refseq[refvar_pos]
						#altbase = seq[contig_pos]
						#vcf.write(rname + "\t" + query + "\t" + str(refvar_pos) + "\t" + refbase + "\t" + altbase + "\n")
					refvar_pos += n
					contig_pos += n
		L = []
		item = "#INFO:The 'short' mark in INFO field means the query is shorter than reference. INDEL or large soft clipping. See the LN field in the header"
		L.append(item)
		item = "#INFO:The 'homopoly' mark in INFO field means there is homopolymer region in refrence. These variants might be caused by sequencing errors."
		L.append(item)
		item = "#REF\tQUERY\tPOS\tREF\tALT\tINFO"
		L.append(item)
		query_len_info = "\n".join(query_len_info_list)
		vcf_header = "\n".join(L)
		var_information = "\n".join(var_info_list)
		if len(ref_lens) == 1:
			ref_len_info = ref_lens[0]
		else:
			ref_len_info = "\n".join(ref_lens)
		vcf.write(ref_len_info + "\n")
		vcf.write(query_len_info)
		vcf.write("\n")
		vcf.write(vcf_header)
		vcf.write("\n")
		vcf.write(var_information)
		vcf.write("\n")
		vcf.close()