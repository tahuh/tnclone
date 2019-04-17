#!/usr/bin/python

import re
import os
import numpy as np
from SequenceParser import FAParser, FQParser

BASE_COMP_RC = { 'A':'T', 'G':'C' , 'T' : 'A' , 'C' : 'G' , 'N' : 'N' }
def FAST_RC(seq):
	ret = []
	for b in seq:
		ret.append(BASE_COMP_RC[b])
	ret = ret[::-1]
	return ''.join(ret)
"""
usage
selector = ConfidentSeqSelector(fasta, vcf, fq1, fq2, k)
selector.Open()
selector.BuildDepthTable()
selector.ReadFASTA()
selector.ReadVCFInfo()
selector.Extract()
selector.FlushSelected(output)
selector.Close()

### Use this for MSA applied contigs
"""
class ConfidentSeqSelector:
	def __init__(self, fasta, vcf, fastq1, fastq2, k):
		self.table = {}
		self.fasta = FAParser(fasta)
		self.fastq1 = FQParser(fastq1)
		self.fastq2 = FQParser(fastq2)
		self.vcf_info = {}
		self.fasta_info = {}
		self.vcf = vcf
		self.found_leap = False
		self.k = k
		self.confident = {} ### id to sequence mapping
	def Open(self):
		self.fasta.open()
			
		self.fastq1.open()
		self.fastq2.open()
		self.vcf = open(self.vcf)

	def Close(self):
		self.fasta.close()
		self.fastq1.close()
		self.fastq2.close()

	def BuildDepthTable(self):
		k = self.k
		for id , seq, qual in self.fastq1.parse():
			if len(seq) <= self.k + 1:
				continue
			rc = FAST_RC(seq)
			for i in range(0, len(seq) - self.k):
				km = seq[i:i+k]
				km2 = rc[i:i+k]
				try:
					self.table[km] += 1
				except KeyError:
					self.table[km] = 1
				try:
					self.table[km2] += 1
				except KeyError:
					self.table[km2] = 1

		for id, seq, qual in self.fastq2.parse():
			if len(seq) <= self.k + 1:
				continue
			rc = FAST_RC(seq)
			for i in range(0, len(seq) - self.k):
				km = seq[i:i+k]
				km2 = rc[i:i+k]
				try:
					self.table[km] += 1
				except KeyError:
					self.table[km] = 1
				try:
					self.table[km2] += 1
				except KeyError:
					self.table[km2] = 1
	def ReadFASTA(self):
		for id , desc, seq in self.fasta.parse():
			if seq == '':
				continue
			self.fasta_info[id] = seq
			
	def ReadVCFInfo(self):
		self.vcf_info = []
		
		for line in self.vcf:
			if line[0] == '#':
				continue
			data = line[:-1].split()
			chrom = data[0]
			pos = int(data[1])
			self.vcf_info.append(pos-1)
		
		
				
	def Extract(self):
		if len(self.fasta_info) == 1:
			return
		result = []
		for chrom in self.fasta_info:
			leap_factor = 1.7
			depths = []
			seq = self.fasta_info[chrom]
			for pos in self.vcf_info:
				kmer = seq[pos - self.k : pos]
				if len(kmer) < self.k:
					continue
				depth = self.table[kmer]
				depths.append(float(depth))
			S = sum(depths)
			for i in range(len(depths)):
				depths[i] = depths[i] / S
			std = np.std(depths)
			mean = np.mean(depths)
			ccv = std / mean
			result.append((chrom, seq, mean, ccv))
		### Sort ascending order
		ccv_sorted = sorted(result, key=lambda x : x[-1], reverse=False)
		ip = 0
		for index in range(len(ccv_sorted)-1):
			ip = index
			cur = ccv_sorted[index][-1]
			nxt = ccv_sorted[index+1][-1]
			rate = nxt / cur
			if rate >= leap_factor :
				self.found_leap = True
				break
		if self.found_leap == False:
			ip = len(ccv_sorted)
			#### Sort by depth. decsnding
			if ip == 0 :
				self.depth_sorted = []
			else:
				self.depth_sorted = sorted(result, key=lambda x : x[-2], reverse=True)
		else:
			if ip > 2 :
				ip = 2 # force select 2 contigs
			for x in range(ip):
				self.confident[ccv_sorted[x][0]] = [ ccv_sorted[x][1] , ccv_sorted[x][2] , ccv_sorted[x][3] ]
	def FlushSelected(self, ofname):
		cnt = 0
		if len(self.vcf_info) == 0 :
			with open(ofname, "w") as F:
				for k in self.fasta_info:
					F.write(">" + k + "_length_" + str(len(self.fasta_info[k])) + "_depth_1.0_ccv_1.0\n" + self.fasta_info[k] + "\n")
					cnt += 1
			return cnt
		if len(self.fasta_info) == 1:
			with open(ofname, "w") as F:
				for k in self.fasta_info:
					F.write(">" + k + "_length_" + str(len(self.fasta_info[k])) + "_depth_1.0_ccv_1.0\n" + self.fasta_info[k] + "\n")
					cnt += 1
			return cnt
		if self.found_leap:
			with open(ofname , "w") as F:
				for chrom , info in self.confident.items():
					cnt += 1
					seq = info[0]
					depth = info[1]
					ccv = info[2]
					F.write(">" + chrom + "_length_%d_depth_%f_ccv_%f\n"%(len(seq), depth, ccv) + seq + "\n")
		else:
			with open(ofname , "w") as F:
				if len(self.depth_sorted) == 0:
					pass
				elif (len(self.depth_sorted) > 0) and (len(self.depth_sorted) <= 2):
					chrom = self.depth_sorted[0][0]
					seq = self.depth_sorted[0][1]
					depth = self.depth_sorted[0][2]
					ccv = self.depth_sorted[0][3]
					F.write(">" + chrom + "_length_%d_depth_%f_ccv_%f\n"%(len(seq), depth, ccv) + seq + "\n")
				else:
					chrom = self.depth_sorted[0][0]
					seq = self.depth_sorted[0][1]
					depth = self.depth_sorted[0][2]
					ccv = self.depth_sorted[0][3]
					F.write(">" + chrom + "_length_%d_depth_%f_ccv_%f\n"%(len(seq), depth, ccv) + seq + "\n")
					
					### Second contig. the max
					chrom = self.depth_sorted[1][0]
					seq = self.depth_sorted[1][1]
					depth = self.depth_sorted[1][2]
					ccv = self.depth_sorted[1][3]
					F.write(">" + chrom + "_length_%d_depth_%f_ccv_%f\n"%(len(seq), depth, ccv) + seq + "\n")
		return cnt

class AnalysisIndelDetector:
	def __init__(self, assem_path, pre_aln_path, selected_out_path, sample_file, ksize):
		self.assem_path = assem_path
		self.pre_aln_path = pre_aln_path
		self.selected_out_path = selected_out_path
		self.sample_file = sample_file
		self.k = ksize
		self.seqs = {}
		self.indel_query = {}
		self.normal_query = {}
	def Run(self):
		self.parse_data()
		self.flush()
	def flush(self):
		with open(self.sample_file) as IN:
			for line in IN:
				data = line.rstrip().split("\t")
				sample = data[0]
				nums = int(data[1])
				for i in range(1,nums+1):
					sample_name = sample + "_" + str(i)
					indel_file = open(self.selected_out_path + os.sep + sample_name + "_k%d"%(self.k) + ".indel.fa", "w")
					normal_file = open(self.selected_out_path + os.sep + sample_name + "_k%d"%(self.k) + ".normal.fa" , "w")
					try:
						qindels = self.indel_query[sample_name]
						for q in qindels:
							seq = self.seqs[sample_name][q]
							indel_file.write(">" + q + "\n" + seq + "\n")
					except KeyError:
						pass
					try:
						qnormals = self.normal_query[sample_name]
						for q in qnormals :
							seq = self.seqs[sample_name][q]
							normal_file.write(">" + q + "\n" + seq + "\n")
					except KeyError:
						pass
					indel_file.close()
					normal_file.close()
	def parse_data(self):
		with open(self.sample_file) as IN:
			for line in IN:
				data = line.rstrip().split("\t")
				sample = data[0]
				nums = int(data[1])
				for i in range(1, nums + 1):
					fasta = self.assem_path + os.sep + sample + "_" + str(i) + "_k%d"%(self.k) + ".contig"
					if not os.path.isfile(fasta):
						continue
					parser = FAParser(fasta)
					parser.open()
					self.seqs[sample + "_" + str(i)] = {}
					for id, desc , seq in parser.parse():
						self.seqs[sample + "_" + str(i)][id] = seq
					file_name = self.pre_aln_path + os.sep + sample + "_" + str(i) + "_k" + str(self.k) + ".pre.dna.sam"
					if not os.path.isfile(file_name):
						continue
					sample_name = sample + "_" + str(i)
					with open(file_name) as SAM:
						for line in SAM:
							if line[0] == '@' : continue
							data = line.rstrip().split("\t")
							query = data[0]
							cigar = data[5]
							alphas = re.split("[0-9]+" , cigar)[1:]
							flag = False
							for a in alphas:
								if (a == 'I') or (a == 'D'):
									flag = True
									# Query with indel
									try:
										self.indel_query[sample_name].add(query)
									except KeyError:
										self.indel_query[sample_name] = set([query])
								else:
									continue
							if not flag :
								# No indels
								try:
									self.normal_query[sample_name].add(query)
								except KeyError:
									self.normal_query[sample_name] = set([query])
class ReportAnalysis:
	pass
