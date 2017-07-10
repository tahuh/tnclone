#!/usr/bin/python

"""
Reporter generates analysis result
Each reporter will report all in one
"""

import glob
import re

BIOPY_AVAIL=False
try:
	from Bio import SeqIO
	BIOPY_AVAIL=True
except ImportError:
	from CustomParser import FastaParser
	BIOPY_AVAIL=False

def get_cigar_num(cigar) :
	return re.split("[A-Z=]+",cigar)[:-1]

def get_cigar_alpha(cigar) :
	return re.split("[0-9]+",cigar)[1:]

class ContigNumberReporter:
	def __init__(self, assem_path,report_fname,pass_dict, prefix , p_set) :
		self.assem_path = assem_path
		#self.n_contig = 0
		self.info = {}
		self.report = open(report_fname ,"aw")
		self.pass_dict = pass_dict # Gene num sample info
		self.prefix = prefix
		self.keylist = []
		self.p_values = p_set
	def run(self,ksize) :
		keyset = self.pass_dict.keys()
		keyset.sort()
		self.keylist = keyset
		for k , n in self.pass_dict.items() :
			gene = k
			num_samples = n
			self.info[gene] = {}
			ncon = 0
			for i in xrange(1,num_samples+1) :
				self.info[gene][i] = 0
				sample = gene + "_" + str(i)
				fname = self.assem_path + "/" + self.prefix +"_" +  sample + "_k" + str(ksize) + ".contig.final.fa"
				try:
					F = open(fname)
				except IOError:
					primary = self.assem_path + "/" + self.prefix +"_" +  sample + "_k" + str(ksize) + ".contig"
					primary_file = open(primary)
					header = primary_file.readline()
					if header.startswith('>'):
						self.info[gene][i] = 'ASSEMBLY FAIL(BROKEN PATH)'
					else:
						self.info[gene][i] = 'ASSEMBLY FAIL(NO SEED KMER PRESENT)'
					continue
				if BIOPY_AVAIL :
					for reocord in SeqIO.parse(F,"fasta") :
						self.info[gene][i] += 1

				else:
					parser = FastaParser(fname)
					parser.open()
					while 1:
						try:
							parser.read()
							self.info[gene][i] += 1
						except StopIteration:
							break
					parser.close()
		self.report.write("[NUMBER OF CONTIGS REPORT]\n")
		for gene in self.info :
			value = self.info[gene]
			value_key = value.keys()
			value_key.sort()
			
			self.report.write(gene)
			for i in value_key :
				if type(self.info[gene][i]) == str:
					self.report.write("\t" + str(self.info[gene][i]))
				else:
					if self.p_values == {}:
						### Empty set
						self.report.write("\t" + str(self.info[gene][i]))
					else:
						if self.p_values[gene][str(i)][0] < 0.05 :
							self.report.write("\t>=3")
						else:
							self.report.write("\t" + str(self.info[gene][i]))
			self.report.write("\n")
		self.report.close()
class AminoAcidReporter:
	def __init__(self, sampath, report_fname,pass_dict,prefix, no_seed_list, broken_list) :
		self.sam = None
		self.report = None
		self.sampath = sampath
		self.report_fname = report_fname
		self.reflen = None
		self.err_free_contigs = []
		self.err_contigs = []
		self.header = None
		self.info = None
		self.prefix = prefix
		self.pass_dict = pass_dict
		self.keylist = []
		self.info = {}
		for e in no_seed_list :
			g = e.split("_")[0]
				
			index = int(e.split("_")[1])
			if not self.info.has_key(g) :
				self.info[g] = {}
				self.info[g][index] = {e:"NOSEED"}
			else:
				self.info[g][index] = {e:"NOSEED"}
		for e in broken_list:
			g = e.split("_")[0]
			index = int(e.split("_")[1])
			if not self.info.has_key(g):
				self.info[g] = {}
				self.info[g][index] = {e:"BROKEN"}
			else:
				self.info[g][index] = {e:"BROKEN"}
	def init(self) :
		#self.sam = open(self.sampath)
		self.report = open(self.report_fname,"aw")

	def clear(self) :
		self.err_free_contigs = []
	def length_match(self, seq) :
		if len(seq) == self.reflen :
			return True
		else:
			return False

	def is_cigar_match_only(self,cigar) :
		rst = get_cigar_alpha(cigar)
		if len(rst) == 0 and rst[0] == 'M' :
			return True
		else:
			return False

	def analyse(self,ksize) :
		keyset = self.pass_dict.keys()
		keyset.sort()
		self.keylist = keyset
		for k , v in self.pass_dict.items() :
			self.info[k] = {}
			for i in xrange(1,v+1) :
				if not self.info.has_key(k):
					if not self.info[k].has_key(i):
						self.info[k][i] = {}
					else:
						continue
				else:
					if not self.info[k].has_key(i):
						self.info[k][i] = {}
					else:
						continue
				try:
					self.sam = open(self.sampath + "/" + self.prefix + "_" + k+"_" + str(i) + "_k" + str(ksize) + ".dna.sam")
				except IOError:
					continue			
				for line in self.sam :
					if line.startswith("@") :
						header = line.rstrip()
						if "LN:" in header :
							self.reflen = int(header.split("LN:")[1])
					else:
						data = line.rstrip().split("\t")
						query = data[0]
						cigar = data[5]
						seq = data[9]
						if self.length_match(seq) :
							self.info[k][i][query] = True
						else:
							self.info[k][i][query] = False
							continue

						if is_cigar_match_only(cigar):
							self.info[k][i][query] = True
						else:
							self.info[k][i][query] = False
				self.sam.close()
		self.report.write("\n[PROTEIN ERROR FREE ANALYSIS]\n")
		for key in self.keylist :
			
			self.report.write(key)
			values = self.info[key]
			values_k = values.keys()
			values_k.sort()
			for index in values_k :
				err_free = ["FREE("]
				err_free.append(")")
				noseed = []
				broken = []
				for query in self.info[key][index] :
					if self.info[key][index][query] == True: err_free.append(query)
					if self.info[key][index][query] == "NOSEED" : noseed.append("NOSEED")
					if self.info[key][index][query] == "BROKEN" : broken.append("BROKEN")
				if len(err_free) == 2 :
					self.report.write("\tERROR")
					if len(noseed) != 0 :
						self.report.write("(NOSEED)")
					if len(broken) != 0 :
						self.report.write("(BROKEN PATH)")
				else:
					self.report.write("\t%s"%(",".join(err_free)))
		self.report.write("\n")
		self.report.close()
class DNAReporter:
	def __init__(self, sampath, report_fname,pass_dict,prefix, no_seed_list, broken_list) :
		self.sam = None
		self.report = None
		self.sampath = sampath
		self.report_fname = report_fname
		self.reflen = None
		self.err_free_contigs = []
		self.err_contigs = []
		self.header = None
		self.info = None
		self.prefix = prefix
		self.pass_dict = pass_dict
		self.keylist = []
		self.info = {}
		for e in no_seed_list :
			g = e.split("_")[0]
				
			index = int(e.split("_")[1])
			if not self.info.has_key(g) :
				self.info[g] = {}
				self.info[g][index] = {e:"NOSEED"}
			else:
				self.info[g][index] = {e:"NOSEED"}
		for e in broken_list:
			g = e.split("_")[0]
			index = int(e.split("_")[1])
			if not self.info.has_key(g):
				self.info[g] = {}
				self.info[g][index] = {e:"BROKEN"}
			else:
				self.info[g][index] = {e:"BROKEN"}

	def init(self) :
		#self.sam = open(self.sampath)
		self.report = open(self.report_fname,"aw")

	def clear(self) :
		self.err_free_contigs = []
	def length_match(self, seq) :
		if len(seq) == self.reflen :
			return True
		else:
			return False

	def is_cigar_match_only(self,cigar) :
		rst = get_cigar_alpha(cigar)
		if len(rst) == 0 and rst[0] == '=' :
			return True
		else:
			return False

	def analyse(self,ksize) :
		keyset = self.pass_dict.keys()
		keyset.sort()
		self.keylist = keyset
		for k , v in self.pass_dict.items() :
			self.info[k] = {}
			for i in xrange(1,v+1) :
				if not self.info.has_key(k):
					if not self.info[k].has_key(i):
						self.info[k][i] = {}
					else:
						continue
				else:
					if not self.info[k].has_key(i):
						self.info[k][i] = {}
					else:
						continue
				try:
					self.sam = open(self.sampath + "/" + self.prefix + "_" + k+"_" + str(i) + "_k" + str(ksize) + ".dna.sam")
				except IOError:
					continue
				for line in self.sam :
					if line.startswith("@") :
						header = line.rstrip()
						if "LN:" in header :
							self.reflen = int(header.split("LN:")[1])
					else:
						data = line.rstrip().split("\t")
						query = data[0]
						cigar = data[5]
						seq = data[9]
						if self.length_match(seq) :
							self.info[k][i][query] = True
						else:
							self.info[k][i][query] = False
							continue

						if is_cigar_match_only(cigar):
							self.info[k][i][query] = True
						else:
							self.info[k][i][query] = False
				self.sam.close()

		for key in self.keylist :
			self.report.write("\n[DNA ERROR FREE ANALYSIS]\n")
			self.report.write(key)
			values = self.info[key]
			values_k = values.keys()
			values_k.sort()
			for index in values_k :
				err_free = ["FREE("]
				err_free.append(")")
				noseed = []
				broken = []
				for query in self.info[key][index] :
					if self.info[key][index][query] == True: err_free.append(query)
					if self.info[key][index][query] == "NOSEED" : noseed.append("NOSEED")
					if self.info[key][index][query] == "BROKEN" : broken.append("BROKEN")
				if len(err_free) == 2 :
					self.report.write("\tERROR")
					if len(noseed) != 0 :
						self.report.write("(NOSEED)")
					if len(broken) != 0 :
						self.report.write("(BROKEN PATH)")
				else:
					self.report.write("\t%s"%(",".join(err_free)))
		self.report.write("\n")
		self.report.close()
